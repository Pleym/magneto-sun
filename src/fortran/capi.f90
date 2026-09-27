!> Interface C des noyaux Fortran, appelée depuis C++ (voir magnetosun_fortran.hpp).
!> Coefficients : g(0:lmax, 0:lmax) rangés par colonnes, soit g(l, m) en
!> g[m * (lmax + 1) + l] côté C++ (idem pour h).
module capi
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use sh_fit
  use sh_rings, only: ring_fit
  use pfss, only: pfss_coefs
  !$ use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: ms_fit, ms_fit_rings, ms_pfss_coefs, ms_pfss_br, ms_max_threads

contains

  !> Solveur dense, pixels quelconques (nlambda >= 3). lambda_in >= 0 : lambda imposé ;
  !> lambda_in < 0 : lambda au coin de la courbe en L. Renvoie aussi la courbe en L,
  !> l'indice de son coin (à partir de 0), le lambda utilisé, le conditionnement
  !> sigma_max / sigma_min et le temps des étapes : matrice, QR + SVD, courbe en L,
  !> solution (s). info : 0 si succès, voir svd_decompose sinon.
  subroutine ms_fit(lmax, npix, cos_theta, phi, br, nlambda, lambda_in, g, h, lambdas, &
                    res_norms, sol_norms, icorner, lambda_used, condition, info, timings) &
                    bind(C, name="ms_fit")
    integer(c_int), value, intent(in) :: lmax, npix, nlambda
    real(c_double), intent(in) :: cos_theta(npix), phi(npix), br(npix)
    real(c_double), value, intent(in) :: lambda_in
    real(c_double), intent(out) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    real(c_double), intent(out) :: lambdas(nlambda), res_norms(nlambda), sol_norms(nlambda)
    integer(c_int), intent(out) :: icorner, info
    real(c_double), intent(out) :: lambda_used, condition, timings(4)
    real(c_double), allocatable :: a(:, :), sigma(:), vt(:, :), beta(:), c(:)
    real(c_double) :: rperp2, t
    integer :: n, i_corner

    n = coef_count(lmax)
    allocate (a(npix, n), sigma(n), vt(n, n), beta(n), c(n))
    t = wall_seconds()
    call design_matrix(lmax, cos_theta, phi, a)
    timings(1) = wall_seconds() - t
    t = wall_seconds()
    call svd_decompose(a, br, sigma, vt, beta, rperp2, info)
    timings(2) = wall_seconds() - t
    if (info /= 0) return

    t = wall_seconds()
    lambdas = lambda_grid(sigma(1), nlambda)
    call lcurve_scan(sigma, beta, rperp2, lambdas, res_norms, sol_norms, i_corner)
    icorner = i_corner - 1
    lambda_used = merge(lambda_in, lambdas(i_corner), lambda_in >= 0)
    condition = sigma(1) / sigma(n)
    timings(3) = wall_seconds() - t
    t = wall_seconds()
    call tikhonov_solve(sigma, vt, beta, lambda_used, c)
    call unpack_coefs(lmax, c, g, h)
    timings(4) = wall_seconds() - t
  end subroutine ms_fit

  !> Solveur par anneaux (même résultat que ms_fit, voir sh_rings) : br(nlon, nrings)
  !> rangée par anneau, x(j) = cos(theta) de l'anneau j, phi0 = longitude du premier
  !> point (radians). Étapes chronométrées : Fourier, blocs, courbe en L, solutions.
  !> info : -1 si nrings <= lmax, -2 si 2 lmax >= nlon.
  subroutine ms_fit_rings(lmax, nlon, nrings, phi0, x, br, nlambda, lambda_in, g, h, lambdas, &
                          res_norms, sol_norms, icorner, lambda_used, condition, info, timings) &
                          bind(C, name="ms_fit_rings")
    integer(c_int), value, intent(in) :: lmax, nlon, nrings, nlambda
    real(c_double), value, intent(in) :: phi0, lambda_in
    real(c_double), intent(in) :: x(nrings), br(nlon, nrings)
    real(c_double), intent(out) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    real(c_double), intent(out) :: lambdas(nlambda), res_norms(nlambda), sol_norms(nlambda)
    integer(c_int), intent(out) :: icorner, info
    real(c_double), intent(out) :: lambda_used, condition, timings(4)
    integer :: i_corner, i_info

    call ring_fit(lmax, phi0, x, br, nlambda, lambda_in, g, h, lambdas, res_norms, sol_norms, &
                  i_corner, lambda_used, condition, i_info, timings)
    icorner = i_corner
    info = i_info
  end subroutine ms_fit_rings

  !> Nombre de threads OpenMP disponibles (1 sans OpenMP).
  integer(c_int) function ms_max_threads() bind(C, name="ms_max_threads")
    ms_max_threads = 1
    !$ ms_max_threads = omp_get_max_threads()
  end function ms_max_threads

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
