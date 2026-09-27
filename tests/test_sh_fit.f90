!> Tests de l'ajustement en harmoniques sphériques sur des cartes synthétiques :
!> 1. carte complète sans bruit : les coefficients sont retrouvés exactement ;
!> 2. calottes polaires masquées + bruit : le problème est mal posé et la
!>    régularisation de Tikhonov choisie par la courbe en L le stabilise.
program test_sh_fit
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use sh_fit
  implicit none
  real(dp), parameter :: PI = acos(-1.0_dp)
  integer :: nfail = 0

  call seed_random()
  call test_exact_recovery_on_full_sphere()
  call test_tikhonov_stabilizes_masked_poles()
  if (nfail > 0) error stop 'test_sh_fit : échec'
  print '(a)', 'sh_fit OK'

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
    call random_seed(put=[(12345 + 7 * i, i=1, n)])
  end subroutine seed_random

  real(dp) function gaussian()
    real(dp) :: u(2)
    call random_number(u)
    gaussian = sqrt(-2 * log(1 - u(1))) * cos(2 * PI * u(2))
  end function gaussian

  !> Coefficients aléatoires d'amplitude 1/(l+1), h(l,0) = 0.
  subroutine random_coefs(lmax, g, h)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax)
    integer :: l, m
    g = 0
    h = 0
    do l = 0, lmax
      do m = 0, l
        g(l, m) = gaussian() / (l + 1)
        if (m > 0) h(l, m) = gaussian() / (l + 1)
      end do
    end do
  end subroutine random_coefs

  !> Centres des pixels d'une grille longitude x sinus de latitude,
  !> en ne gardant que |sin(latitude)| <= max_abs_sinlat.
  subroutine grid_pixels(nlon, nlat, max_abs_sinlat, cos_theta, phi)
    integer, intent(in) :: nlon, nlat
    real(dp), intent(in) :: max_abs_sinlat
    real(dp), allocatable, intent(out) :: cos_theta(:), phi(:)
    real(dp) :: all_x(nlon * nlat), all_phi(nlon * nlat)
    integer :: i, j, k
    k = 0
    do j = 0, nlat - 1
      do i = 0, nlon - 1
        k = k + 1
        all_x(k) = -1 + (j + 0.5_dp) * 2 / nlat
        all_phi(k) = (i + 0.5_dp) * 2 * PI / nlon
      end do
    end do
    cos_theta = pack(all_x, abs(all_x) <= max_abs_sinlat)
    phi = pack(all_phi, abs(all_x) <= max_abs_sinlat)
  end subroutine grid_pixels

  subroutine test_exact_recovery_on_full_sphere()
    integer, parameter :: LMAX = 8
    real(dp) :: g(0:LMAX, 0:LMAX), h(0:LMAX, 0:LMAX), gfit(0:LMAX, 0:LMAX), hfit(0:LMAX, 0:LMAX)
    real(dp), allocatable :: x(:), phi(:), b(:), a(:, :), sigma(:), vt(:, :), beta(:), c(:)
    real(dp) :: rperp2
    integer :: n, info

    ! Arrange
    call random_coefs(LMAX, g, h)
    call grid_pixels(72, 36, 1.0_dp, x, phi)
    allocate (b(size(x)))
    call synthesis(LMAX, g, h, x, phi, b)
    n = coef_count(LMAX)
    allocate (a(size(x), n), sigma(n), vt(n, n), beta(n), c(n))

    ! Act
    call design_matrix(LMAX, x, phi, a)
    call svd_decompose(a, b, sigma, vt, beta, rperp2, info)
    call tikhonov_solve(sigma, vt, beta, 0.0_dp, c)
    call unpack_coefs(LMAX, c, gfit, hfit)

    ! Assert
    call check(info == 0, 'LAPACK sans erreur')
    call check(maxval(abs(gfit - g)) < 1e-10_dp .and. maxval(abs(hfit - h)) < 1e-10_dp, &
               'coefficients retrouvés exactement sur la sphère complète')
    call check(rperp2 < 1e-20_dp, 'résidu nul pour une carte sans bruit')
  end subroutine test_exact_recovery_on_full_sphere

  subroutine test_tikhonov_stabilizes_masked_poles()
    ! À ce degré, masquer |latitude| > 60° rend le système mal conditionné (> 1e4)
    integer, parameter :: LMAX = 24, NLAMBDA = 200
    real(dp), parameter :: MAX_ABS_SINLAT = sqrt(3.0_dp) / 2   ! |latitude| <= 60°
    real(dp) :: g(0:LMAX, 0:LMAX), h(0:LMAX, 0:LMAX)
    real(dp) :: truth((LMAX + 1)**2), lambdas(NLAMBDA), errors(NLAMBDA)
    real(dp) :: res_norms(NLAMBDA), sol_norms(NLAMBDA), rperp2, noise, error_unregularized
    real(dp), allocatable :: x(:), phi(:), b(:), a(:, :), sigma(:), vt(:, :), beta(:), c(:)
    integer :: n, info, i, icorner
    character(96) :: label

    ! Arrange
    call random_coefs(LMAX, g, h)
    call pack_coefs(LMAX, g, h, truth)
    call grid_pixels(180, 90, MAX_ABS_SINLAT, x, phi)
    allocate (b(size(x)))
    call synthesis(LMAX, g, h, x, phi, b)
    noise = 0.01_dp * sqrt(sum(b**2) / size(b))
    b = b + [(noise * gaussian(), i=1, size(b))]
    n = coef_count(LMAX)
    allocate (a(size(x), n), sigma(n), vt(n, n), beta(n), c(n))

    ! Act
    call design_matrix(LMAX, x, phi, a)
    call svd_decompose(a, b, sigma, vt, beta, rperp2, info)
    lambdas = [(sigma(1) * 10.0_dp**(-8 + 8 * (i - 1.0_dp) / (NLAMBDA - 1)), i=1, NLAMBDA)]
    do i = 1, NLAMBDA
      call tikhonov_norms(sigma, beta, rperp2, lambdas(i), res_norms(i), sol_norms(i))
      call tikhonov_solve(sigma, vt, beta, lambdas(i), c)
      errors(i) = norm2(c - truth)
    end do
    icorner = lcurve_corner(res_norms, sol_norms)
    call tikhonov_solve(sigma, vt, beta, 0.0_dp, c)
    error_unregularized = norm2(c - truth)

    ! Assert
    write (label, '(a,es9.2,a,es9.2,a,es9.2)') 'conditionnement ', sigma(1) / sigma(n), &
      ', erreur sans régularisation ', error_unregularized, ', au coin ', errors(icorner)
    print '(a)', trim(label)
    call check(info == 0, 'LAPACK sans erreur')
    call check(sigma(1) / sigma(n) > 1e4_dp, 'calottes masquées : problème mal conditionné')
    call check(errors(icorner) < 0.1_dp * error_unregularized, &
               'la régularisation au coin de la courbe en L divise l''erreur par plus de 10')
    call check(errors(icorner) < 3 * minval(errors), &
               'le coin de la courbe en L est proche du meilleur lambda possible')
  end subroutine test_tikhonov_stabilizes_masked_poles

end program test_sh_fit
