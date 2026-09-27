!> Interface C des noyaux Fortran, appelée depuis C++ (voir magnetosun_fortran.hpp).
!> Coefficients : g(0:lmax, 0:lmax) rangés par colonnes, soit g(l, m) en
!> g[m * (lmax + 1) + l] côté C++ (idem pour h).
module capi
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use sh_fit
  use pfss, only: pfss_coefs
  implicit none
  private
  public :: ms_fit, ms_pfss_coefs, ms_pfss_br

  ! Plus petit lambda exploré pour la courbe en L, en fraction de sigma_max.
  real(c_double), parameter :: LAMBDA_MIN_RATIO = 1e-8_c_double

contains

  !> Ajuste g, h sur npix pixels (nlambda >= 3). lambda_in >= 0 : lambda imposé ;
  !> lambda_in < 0 : lambda au coin de la courbe en L. Renvoie aussi la courbe en L,
  !> l'indice de son coin (à partir de 0), le lambda utilisé et le conditionnement
  !> sigma_max / sigma_min. info : 0 si succès, voir svd_decompose sinon.
  subroutine ms_fit(lmax, npix, cos_theta, phi, br, nlambda, lambda_in, g, h, lambdas, &
                    res_norms, sol_norms, icorner, lambda_used, condition, info) &
                    bind(C, name="ms_fit")
    integer(c_int), value, intent(in) :: lmax, npix, nlambda
    real(c_double), intent(in) :: cos_theta(npix), phi(npix), br(npix)
    real(c_double), value, intent(in) :: lambda_in
    real(c_double), intent(out) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    real(c_double), intent(out) :: lambdas(nlambda), res_norms(nlambda), sol_norms(nlambda)
    integer(c_int), intent(out) :: icorner, info
    real(c_double), intent(out) :: lambda_used, condition
    real(c_double), allocatable :: a(:, :), sigma(:), vt(:, :), beta(:), c(:)
    real(c_double) :: rperp2
    integer :: n, i

    n = coef_count(lmax)
    allocate (a(npix, n), sigma(n), vt(n, n), beta(n), c(n))
    call design_matrix(lmax, cos_theta, phi, a)
    call svd_decompose(a, br, sigma, vt, beta, rperp2, info)
    if (info /= 0) return

    lambdas = [(sigma(1) * LAMBDA_MIN_RATIO**(1 - (i - 1.0_c_double) / (nlambda - 1)), i=1, nlambda)]
    do i = 1, nlambda
      call tikhonov_norms(sigma, beta, rperp2, lambdas(i), res_norms(i), sol_norms(i))
    end do
    icorner = lcurve_corner(res_norms, sol_norms) - 1
    lambda_used = merge(lambda_in, lambdas(icorner + 1), lambda_in >= 0)
    condition = sigma(1) / sigma(n)
    call tikhonov_solve(sigma, vt, beta, lambda_used, c)
    call unpack_coefs(lmax, c, g, h)
  end subroutine ms_fit

  !> Coefficients de B_r du modèle PFSS au rayon r (monopôle retiré).
  subroutine ms_pfss_coefs(lmax, g, h, rss, r, gr, hr) bind(C, name="ms_pfss_coefs")
    integer(c_int), value, intent(in) :: lmax
    real(c_double), intent(in) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    real(c_double), value, intent(in) :: rss, r
    real(c_double), intent(out) :: gr(0:lmax, 0:lmax), hr(0:lmax, 0:lmax)
    call pfss_coefs(lmax, g, h, r, rss, gr, hr)
  end subroutine ms_pfss_coefs

  !> B_r du modèle PFSS au rayon r (en rayons solaires), surface source en rss.
  subroutine ms_pfss_br(lmax, g, h, rss, r, n, cos_theta, phi, values) bind(C, name="ms_pfss_br")
    integer(c_int), value, intent(in) :: lmax, n
    real(c_double), intent(in) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax), cos_theta(n), phi(n)
    real(c_double), value, intent(in) :: rss, r
    real(c_double), intent(out) :: values(n)
    real(c_double) :: gr(0:lmax, 0:lmax), hr(0:lmax, 0:lmax)
    call pfss_coefs(lmax, g, h, r, rss, gr, hr)
    call synthesis(lmax, gr, hr, cos_theta, phi, values)
  end subroutine ms_pfss_br

end module capi
