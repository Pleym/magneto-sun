!> Ajustement rapide et exact quand les pixels forment des anneaux complets de
!> longitude régulière (grille synoptique, rangées entières).
!>
!> Sur un anneau de nlon points équidistants, cos(m phi) et sin(m phi) sont orthogonaux
!> pour 2m < nlon : A^T A est diagonale par blocs (m ; cos ou sin), comme le terme de
!> Tikhonov lambda^2 |c|^2. Le problème se découpe donc exactement en 2 lmax + 1
!> moindres carrés en latitude, de taille nanneaux x (lmax - m + 1). Pour l'anneau j :
!>   F_m(j) = somme_i b(i, j) cos(m phi_i)   (idem avec sin),
!>   N_m = nlon si m = 0, nlon / 2 sinon,
!> et le bloc (m, cos) minimise somme_j N_m [F_m(j) / N_m - somme_l P(l, m)(x_j) g(l, m)]^2.
!> Coût : O(nlon nanneaux lmax) + O(nanneaux lmax^3) + O(lmax^4), contre
!> O(npix lmax^4) + O(lmax^6) pour le solveur dense.
module sh_rings
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use legendre, only: schmidt_diagonal, schmidt_column
  use sh_fit, only: svd_decompose_multi, tikhonov_solve, lambda_grid, lcurve_scan, wall_seconds
  implicit none
  private
  public :: ring_fit

  !> Décomposition du bloc m, commune au cosinus et au sinus (même matrice) :
  !> beta(:, 1) pour cos(m phi), beta(:, 2) pour sin(m phi) si m > 0.
  type :: block_svd
    real(dp), allocatable :: sigma(:), vt(:, :), beta(:, :), rperp2(:)
    integer :: info = 0
  end type block_svd

contains

  !> Même rôle que ms_fit pour des anneaux : br(nlon, nanneaux), x(j) = cos(theta) de
  !> l'anneau j, phi0 = longitude (radians) du premier point de chaque anneau.
  !> timings : Fourier, blocs, courbe en L, solutions (s).
  !> info : 0 si succès, -1 si nanneaux <= lmax, -2 si 2 lmax >= nlon, > 0 LAPACK.
  subroutine ring_fit(lmax, phi0, x, br, nlambda, lambda_in, g, h, lambdas, res_norms, &
                      sol_norms, icorner, lambda_used, condition, info, timings)
    integer, intent(in) :: lmax, nlambda
    real(dp), intent(in) :: phi0, x(:), br(:, :), lambda_in
    real(dp), intent(out) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    real(dp), intent(out) :: lambdas(nlambda), res_norms(nlambda), sol_norms(nlambda)
    integer, intent(out) :: icorner, info
    real(dp), intent(out) :: lambda_used, condition, timings(4)
    real(dp), allocatable :: fc(:, :), fs(:, :), diag(:, :), sigma_all(:), beta_all(:)
    type(block_svd), allocatable :: blocks(:)
    real(dp) :: remainder, t
    integer :: m, nlon

    nlon = size(br, 1)
    info = merge(-2, merge(-1, 0, size(x) <= lmax), 2 * lmax >= nlon)
    if (info /= 0) return

    t = wall_seconds()
    call ring_fourier(lmax, phi0, br, fc, fs, remainder)
    timings(1) = wall_seconds() - t

    t = wall_seconds()
    allocate (diag(size(x), 0:lmax), blocks(0:lmax))
    call schmidt_diagonal(lmax, x, diag)
    ! Blocs indépendants ; les plus gros (m petit) d'abord pour équilibrer la charge
    !$omp parallel do schedule(dynamic)
    do m = 0, lmax
      if (m == 0) then
        call solve_block(m, lmax, nlon, x, diag(:, m), reshape(fc(m, :), [size(x), 1]), blocks(m))
      else
        call solve_block(m, lmax, nlon, x, diag(:, m), reshape([fc(m, :), fs(m, :)], [size(x), 2]), &
                         blocks(m))
      end if
    end do
    !$omp end parallel do
    info = maxval(blocks%info)
    timings(2) = wall_seconds() - t
    if (info /= 0) return

    t = wall_seconds()
    ! Spectre global : chaque bloc m > 0 compte deux fois (cos et sin), avec son beta
    sigma_all = [(spread(blocks(m)%sigma, 2, size(blocks(m)%beta, 2)), m=0, lmax)]
    beta_all = [(blocks(m)%beta, m=0, lmax)]
    lambdas = lambda_grid(maxval(sigma_all), nlambda)
    call lcurve_scan(sigma_all, beta_all, remainder + sum([(sum(blocks(m)%rperp2), m=0, lmax)]), &
                     lambdas, res_norms, sol_norms, icorner)
    lambda_used = merge(lambda_in, lambdas(icorner), lambda_in >= 0)
    icorner = icorner - 1
    condition = maxval(sigma_all) / minval(sigma_all)
    timings(3) = wall_seconds() - t

    t = wall_seconds()
    g = 0
    h = 0
    do m = 0, lmax
      call tikhonov_solve(blocks(m)%sigma, blocks(m)%vt, blocks(m)%beta(:, 1), lambda_used, &
                          g(m:lmax, m))
      if (m > 0) then
        call tikhonov_solve(blocks(m)%sigma, blocks(m)%vt, blocks(m)%beta(:, 2), lambda_used, &
                            h(m:lmax, m))
      end if
    end do
    timings(4) = wall_seconds() - t
  end subroutine ring_fit

  !> fc(m, j), fs(m, j) : sommes de b(i, j) cos(m phi_i) et sin(m phi_i) sur l'anneau j,
  !> par produit matriciel BLAS (dgemm) ; remainder = |b|^2 - somme F^2 / N_m, la part
  !> de la carte hors des fréquences m <= lmax.
  subroutine ring_fourier(lmax, phi0, br, fc, fs, remainder)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: phi0, br(:, :)
    real(dp), allocatable, intent(out) :: fc(:, :), fs(:, :)
    real(dp), intent(out) :: remainder
    real(dp), parameter :: PI = acos(-1.0_dp)
    external :: dgemm
    real(dp), allocatable :: trig_c(:, :), trig_s(:, :)
    real(dp) :: phi
    integer :: i, m, nlon, nrings

    nlon = size(br, 1)
    nrings = size(br, 2)
    allocate (trig_c(nlon, 0:lmax), trig_s(nlon, 0:lmax))
    ! ponytail: somme directe en O(nlon lmax) par anneau ; une FFT si nlon devient grand
    !$omp parallel do private(i, phi)
    do m = 0, lmax
      do i = 1, nlon
        phi = phi0 + (i - 1) * 2 * PI / nlon
        trig_c(i, m) = cos(m * phi)
        trig_s(i, m) = sin(m * phi)
      end do
    end do
    !$omp end parallel do
    allocate (fc(0:lmax, nrings), fs(0:lmax, nrings))
    ! fc = trig_c^T br : le gfortran intrinsèque (matmul) plafonnait à ~4 GFLOPS
    call dgemm('T', 'N', lmax + 1, nrings, nlon, 1.0_dp, trig_c, nlon, br, nlon, 0.0_dp, fc, lmax + 1)
    call dgemm('T', 'N', lmax + 1, nrings, nlon, 1.0_dp, trig_s, nlon, br, nlon, 0.0_dp, fs, lmax + 1)
    remainder = sum(br**2) - sum(fc(0, :)**2) / nlon &
                - sum(fc(1:lmax, :)**2 + fs(1:lmax, :)**2) / (nlon / 2.0_dp)
  end subroutine ring_fourier

  !> Bloc m : matrice sqrt(N_m) P(l, m)(x_j) ; seconds membres F(j) / sqrt(N_m), rangés
  !> f(j, k) : k = 1 pour cos(m phi), k = 2 pour sin(m phi).
  subroutine solve_block(m, lmax, nlon, x, pmm, f, block)
    integer, intent(in) :: m, lmax, nlon
    real(dp), intent(in) :: x(:), pmm(:), f(:, :)
    type(block_svd), intent(out) :: block
    real(dp), allocatable :: b(:, :)
    real(dp) :: norm
    integer :: n, nrhs

    n = lmax - m + 1
    nrhs = size(f, 2)
    norm = merge(real(nlon, dp), nlon / 2.0_dp, m == 0)
    allocate (b(size(x), n), block%sigma(n), block%vt(n, n), block%beta(n, nrhs), &
              block%rperp2(nrhs))
    call schmidt_column(m, lmax, x, pmm, b)
    b = b * sqrt(norm)
    call svd_decompose_multi(b, f / sqrt(norm), block%sigma, block%vt, block%beta, &
                             block%rperp2, block%info)
  end subroutine solve_block

end module sh_rings
