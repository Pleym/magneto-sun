!> Le solveur par anneaux doit donner exactement la même solution que le solveur
!> dense (mêmes coefficients pour un lambda donné, même courbe en L, même coin),
!> sur une grille dont on a retiré des rangées entières (pôles masqués).
program test_sh_rings
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use sh_fit
  use sh_rings, only: ring_fit
  implicit none
  real(dp), parameter :: PI = acos(-1.0_dp)
  integer, parameter :: NLON = 48, NLAT = 24, LMAX = 10, NLAMBDA = 60
  real(dp), parameter :: MAX_ABS_SINLAT = 0.8_dp
  integer :: nfail = 0

  call seed_random()
  call test_same_solution_as_dense_solver()
  call test_rejects_unresolvable_degrees()
  if (nfail > 0) error stop 'test_sh_rings : échec'
  print '(a)', 'sh_rings OK'

contains

  subroutine check(is_ok, what)
    logical, intent(in) :: is_ok
    character(*), intent(in) :: what
    if (.not. is_ok) then
      print '(2a)', 'ÉCHEC : ', what
      nfail = nfail + 1
    end if
  end subroutine check

  subroutine seed_random()
    integer :: n, i
    call random_seed(size=n)
    call random_seed(put=[(4242 + 3 * i, i=1, n)])
  end subroutine seed_random

  !> Anneaux de sin(latitude) <= MAX_ABS_SINLAT, carte = harmoniques aléatoires + bruit.
  subroutine synthetic_rings(x, br)
    real(dp), allocatable, intent(out) :: x(:), br(:, :)
    real(dp) :: g(0:LMAX, 0:LMAX), h(0:LMAX, 0:LMAX), all_x(NLAT), phi(NLON), noise(NLON)
    integer :: i, j
    call random_number(g)
    call random_number(h)
    g = g - 0.5_dp
    h = h - 0.5_dp
    h(:, 0) = 0
    all_x = [(-1 + (j - 0.5_dp) * 2 / NLAT, j=1, NLAT)]
    x = pack(all_x, abs(all_x) <= MAX_ABS_SINLAT)
    phi = [((i - 0.5_dp) * 2 * PI / NLON, i=1, NLON)]
    allocate (br(NLON, size(x)))
    do j = 1, size(x)
      call synthesis(LMAX, g, h, [(x(j), i=1, NLON)], phi, br(:, j))
      call random_number(noise)
      br(:, j) = br(:, j) + 0.05_dp * (noise - 0.5_dp)
    end do
  end subroutine synthetic_rings

  subroutine test_same_solution_as_dense_solver()
    real(dp), allocatable :: x(:), br(:, :), a(:, :), sigma(:), vt(:, :), beta(:), c(:)
    real(dp) :: g(0:LMAX, 0:LMAX), h(0:LMAX, 0:LMAX), gd(0:LMAX, 0:LMAX), hd(0:LMAX, 0:LMAX)
    real(dp) :: lambdas(NLAMBDA), res(NLAMBDA), sol(NLAMBDA)
    real(dp) :: lambdas_d(NLAMBDA), res_d(NLAMBDA), sol_d(NLAMBDA)
    real(dp) :: rperp2, lambda_used, condition, timings(4), lambda
    integer :: n, info, icorner, icorner_d, i, j

    ! Arrange : solveur dense sur les mêmes pixels
    call synthetic_rings(x, br)
    n = coef_count(LMAX)
    allocate (a(NLON * size(x), n), sigma(n), vt(n, n), beta(n), c(n))
    call design_matrix(LMAX, [((x(j), i=1, NLON), j=1, size(x))], &
                       [(((i - 0.5_dp) * 2 * PI / NLON, i=1, NLON), j=1, size(x))], a)
    call svd_decompose(a, reshape(br, [NLON * size(x)]), sigma, vt, beta, rperp2, info)
    lambdas_d = lambda_grid(sigma(1), NLAMBDA)
    call lcurve_scan(sigma, beta, rperp2, lambdas_d, res_d, sol_d, icorner_d)
    lambda = lambdas_d(NLAMBDA / 3)
    call tikhonov_solve(sigma, vt, beta, lambda, c)
    call unpack_coefs(LMAX, c, gd, hd)

    ! Act
    call ring_fit(LMAX, PI / NLON, x, br, NLAMBDA, lambda, g, h, lambdas, res, sol, icorner, &
                  lambda_used, condition, info, timings)

    ! Assert
    call check(info == 0, 'ring_fit sans erreur')
    call check(maxval(abs(lambdas / lambdas_d - 1)) < 1e-12_dp, 'même grille de lambda')
    call check(maxval(abs(res / res_d - 1)) < 1e-9_dp .and. maxval(abs(sol / sol_d - 1)) < 1e-9_dp, &
               'même courbe en L')
    call check(icorner == icorner_d - 1, 'même coin de la courbe en L')
    call check(maxval(abs(g - gd)) + maxval(abs(h - hd)) < 1e-10_dp, &
               'mêmes coefficients pour un lambda donné')
    call check(abs(condition - sigma(1) / sigma(n)) < 1e-8_dp * condition, 'même conditionnement')
  end subroutine test_same_solution_as_dense_solver

  subroutine test_rejects_unresolvable_degrees()
    real(dp), allocatable :: x(:), br(:, :)
    real(dp) :: g(0:30, 0:30), h(0:30, 0:30), lambdas(3), res(3), sol(3), lu, cond, t(4)
    integer :: icorner, info
    call synthetic_rings(x, br)
    ! 2 lmax >= nlon : cos(m phi) n'est plus résolu par 48 points par anneau
    call ring_fit(30, PI / NLON, x, br, 3, -1.0_dp, g, h, lambdas, res, sol, icorner, lu, cond, info, t)
    call check(info == -2, 'lmax >= nlon/2 refusé')
    ! moins d'anneaux (20) que de degrés (lmax + 1 = 21)
    call ring_fit(20, PI / NLON, x, br, 3, -1.0_dp, g(0:20, 0:20), h(0:20, 0:20), lambdas, res, &
                  sol, icorner, lu, cond, info, t)
    call check(info == -1, 'moins d''anneaux que de degrés refusé')
  end subroutine test_rejects_unresolvable_degrees

end program test_sh_rings
