!> Décomposition de B_r en harmoniques sphériques par moindres carrés, avec
!> régularisation de Tikhonov et choix du paramètre par la courbe en L.
!>
!> Modèle : B_r(theta, phi) = somme pour 0 <= m <= l <= lmax de
!>          P(l, m)(cos theta) [g(l, m) cos(m phi) + h(l, m) sin(m phi)]
!> (P de Schmidt, module legendre ; h(l, 0) = 0). Dans les systèmes linéaires, les
!> inconnues forment un vecteur c de taille (lmax+1)^2 rangé ainsi : pour chaque l,
!> g(l, 0), puis g(l, m), h(l, m) pour m = 1..l.
module sh_fit
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use legendre, only: schmidt_legendre
  implicit none
  private
  public :: coef_count, design_matrix, svd_decompose, tikhonov_solve, tikhonov_norms, &
            lcurve_corner, pack_coefs, unpack_coefs, synthesis

contains

  pure integer function coef_count(lmax)
    integer, intent(in) :: lmax
    coef_count = (lmax + 1)**2
  end function coef_count

  !> a(ipix, k) = valeur de l'harmonique k au pixel ipix.
  subroutine design_matrix(lmax, cos_theta, phi, a)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: cos_theta(:), phi(:)
    real(dp), intent(out) :: a(:, :)
    real(dp) :: p(0:lmax, 0:lmax)
    integer :: ipix, l, m, k

    ! ponytail: remplissage ligne par ligne (accès mémoire en pas de npix), à revoir en Phase 5
    do ipix = 1, size(cos_theta)
      call schmidt_legendre(lmax, cos_theta(ipix), p)
      k = 0
      do l = 0, lmax
        a(ipix, k + 1) = p(l, 0)
        do m = 1, l
          a(ipix, k + 2 * m) = p(l, m) * cos(m * phi(ipix))
          a(ipix, k + 2 * m + 1) = p(l, m) * sin(m * phi(ipix))
        end do
        k = k + 2 * l + 1
      end do
    end do
  end subroutine design_matrix

  !> Décompose A = (Q U) S V^T : QR de Householder A = Q R, puis SVD de la petite
  !> matrice carrée R = U S V^T. Renvoie S (sigma, décroissants), V^T, beta = (Q U)^T b
  !> et rperp2 = |b - projection de b sur l'image de A|^2. La matrice a est écrasée.
  !> info : 0 si succès, -1 si moins de pixels que d'inconnues, > 0 erreur LAPACK.
  subroutine svd_decompose(a, b, sigma, vt, beta, rperp2, info)
    real(dp), contiguous, intent(inout) :: a(:, :)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: sigma(:), vt(:, :), beta(:), rperp2
    integer, intent(out) :: info
    external :: dgeqrf, dormqr, dgesdd
    real(dp), allocatable :: tau(:), work(:), qtb(:), r(:, :), u(:, :)
    integer, allocatable :: iwork(:)
    real(dp) :: query(1)
    integer :: m, n, j

    m = size(a, 1)
    n = size(a, 2)
    info = -1
    if (m < n) return
    allocate (tau(n), qtb(m), r(n, n), u(n, n), iwork(8 * n))
    qtb = b

    call dgeqrf(m, n, a, m, tau, query, -1, info)
    allocate (work(int(query(1))))
    call dgeqrf(m, n, a, m, tau, work, size(work), info)
    if (info /= 0) return
    call dormqr('L', 'T', m, 1, n, a, m, tau, qtb, m, query, -1, info)
    call ensure_size(work, int(query(1)))
    call dormqr('L', 'T', m, 1, n, a, m, tau, qtb, m, work, size(work), info)
    if (info /= 0) return
    rperp2 = sum(qtb(n + 1:m)**2)

    r = 0
    do j = 1, n
      r(1:j, j) = a(1:j, j)
    end do
    call dgesdd('A', n, n, r, n, sigma, u, n, vt, n, query, -1, iwork, info)
    call ensure_size(work, int(query(1)))
    call dgesdd('A', n, n, r, n, sigma, u, n, vt, n, work, size(work), iwork, info)
    if (info /= 0) return
    beta = matmul(qtb(1:n), u)
  end subroutine svd_decompose

  subroutine ensure_size(work, n)
    real(dp), allocatable, intent(inout) :: work(:)
    integer, intent(in) :: n
    if (size(work) >= n) return
    deallocate (work)
    allocate (work(n))
  end subroutine ensure_size

  !> Solution de Tikhonov : minimise |A c - b|^2 + lambda^2 |c|^2.
  !> c = somme sur i de sigma_i beta_i / (sigma_i^2 + lambda^2) v_i ;
  !> lambda = 0 donne la solution des moindres carrés ordinaires.
  pure subroutine tikhonov_solve(sigma, vt, beta, lambda, c)
    real(dp), intent(in) :: sigma(:), vt(:, :), beta(:), lambda
    real(dp), intent(out) :: c(:)
    c = matmul(filtered(sigma, beta, lambda), vt)
  end subroutine tikhonov_solve

  pure function filtered(sigma, beta, lambda) result(w)
    real(dp), intent(in) :: sigma(:), beta(:), lambda
    real(dp) :: w(size(sigma))
    where (sigma**2 + lambda**2 > 0)
      w = sigma * beta / (sigma**2 + lambda**2)
    elsewhere
      w = 0
    end where
  end function filtered

  !> Norme du résidu |A c - b| et de la solution |c|, sans calculer c.
  pure subroutine tikhonov_norms(sigma, beta, rperp2, lambda, res_norm, sol_norm)
    real(dp), intent(in) :: sigma(:), beta(:), rperp2, lambda
    real(dp), intent(out) :: res_norm, sol_norm
    sol_norm = norm2(filtered(sigma, beta, lambda))
    res_norm = sqrt(sum((lambda**2 * beta / (sigma**2 + lambda**2))**2) + rperp2)
  end subroutine tikhonov_norms

  !> Coin de la courbe en L (Hansen 1992) : point de courbure maximale de
  !> (log |r|, log |c|), les lambda étant régulièrement espacés en log.
  pure integer function lcurve_corner(res_norms, sol_norms) result(icorner)
    real(dp), intent(in) :: res_norms(:), sol_norms(:)
    real(dp) :: x(size(res_norms)), y(size(res_norms)), dx, dy, ddx, ddy, kappa, best
    integer :: i

    x = log(res_norms)
    y = log(sol_norms)
    icorner = 1
    best = -huge(1.0_dp)
    do i = 2, size(x) - 1
      dx = (x(i + 1) - x(i - 1)) / 2
      dy = (y(i + 1) - y(i - 1)) / 2
      ddx = x(i + 1) - 2 * x(i) + x(i - 1)
      ddy = y(i + 1) - 2 * y(i) + y(i - 1)
      if (dx**2 + dy**2 <= 0) cycle
      kappa = (dx * ddy - ddx * dy) / (dx**2 + dy**2)**1.5_dp
      if (kappa > best) then
        best = kappa
        icorner = i
      end if
    end do
  end function lcurve_corner

  pure subroutine pack_coefs(lmax, g, h, c)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    real(dp), intent(out) :: c(:)
    integer :: l, m
    do l = 0, lmax
      c(l**2 + 1) = g(l, 0)
      do m = 1, l
        c(l**2 + 2 * m) = g(l, m)
        c(l**2 + 2 * m + 1) = h(l, m)
      end do
    end do
  end subroutine pack_coefs

  pure subroutine unpack_coefs(lmax, c, g, h)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: c(:)
    real(dp), intent(out) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    integer :: l, m
    g = 0
    h = 0
    do l = 0, lmax
      g(l, 0) = c(l**2 + 1)
      do m = 1, l
        g(l, m) = c(l**2 + 2 * m)
        h(l, m) = c(l**2 + 2 * m + 1)
      end do
    end do
  end subroutine unpack_coefs

  !> B_r aux points (cos_theta, phi) à partir des coefficients g, h.
  subroutine synthesis(lmax, g, h, cos_theta, phi, values)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax), cos_theta(:), phi(:)
    real(dp), intent(out) :: values(:)
    real(dp) :: p(0:lmax, 0:lmax), cos_m(0:lmax), sin_m(0:lmax)
    integer :: i, m

    do i = 1, size(cos_theta)
      call schmidt_legendre(lmax, cos_theta(i), p)
      cos_m = [(cos(m * phi(i)), m=0, lmax)]
      sin_m = [(sin(m * phi(i)), m=0, lmax)]
      values(i) = sum(p * (g * spread(cos_m, 1, lmax + 1) + h * spread(sin_m, 1, lmax + 1)))
    end do
  end subroutine synthesis

end module sh_fit
