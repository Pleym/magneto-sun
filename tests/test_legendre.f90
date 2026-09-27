!> Tests des fonctions de Legendre de Schmidt : formes explicites, puis
!> normalisation et orthogonalité par quadrature.
program test_legendre
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use legendre, only: schmidt_legendre
  implicit none
  integer, parameter :: LMAX = 40
  integer :: nfail = 0

  call test_closed_forms()
  call test_normalization(0)
  call test_normalization(1)
  call test_normalization(7)
  call test_normalization(30)
  if (nfail > 0) error stop 'test_legendre : échec'
  print '(a)', 'legendre OK'

contains

  subroutine check(is_ok, what)
    logical, intent(in) :: is_ok
    character(*), intent(in) :: what
    if (.not. is_ok) then
      print '(2a)', 'ÉCHEC : ', what
      nfail = nfail + 1
    end if
  end subroutine check

  subroutine test_closed_forms()
    real(dp) :: p(0:LMAX, 0:LMAX), x, s
    x = 0.3_dp
    s = sqrt(1 - x * x)
    call schmidt_legendre(LMAX, x, p)
    call check(abs(p(0, 0) - 1) < 1e-15_dp, 'P(0,0) = 1')
    call check(abs(p(1, 0) - x) < 1e-15_dp, 'P(1,0) = cos theta')
    call check(abs(p(1, 1) - s) < 1e-15_dp, 'P(1,1) = sin theta')
    call check(abs(p(2, 0) - (3 * x * x - 1) / 2) < 1e-15_dp, 'P(2,0)')
    call check(abs(p(2, 1) - sqrt(3.0_dp) * x * s) < 1e-15_dp, 'P(2,1)')
    call check(abs(p(2, 2) - sqrt(3.0_dp) / 2 * s**2) < 1e-15_dp, 'P(2,2)')
    call check(abs(p(3, 3) - sqrt(10.0_dp) / 4 * s**3) < 1e-15_dp, 'P(3,3)')
    call check(p(3, 5) == 0, 'P(l,m) = 0 pour m > l')
  end subroutine test_closed_forms

  !> Schmidt : intégrale de P(l,m) P(l',m) sur [-1, 1] = delta_ll' * 2(2 - delta_m0)/(2l+1).
  !> Quadrature du point milieu en theta, d'erreur relative ~ (l pi / 2N)^2 / 6.
  subroutine test_normalization(m)
    integer, intent(in) :: m
    integer, parameter :: NQUAD = 40000
    real(dp), parameter :: PI = acos(-1.0_dp)
    real(dp) :: p(0:LMAX, 0:LMAX), gram(0:LMAX, 0:LMAX), theta, expected
    integer :: i, l, l2
    character(64) :: label

    gram = 0
    do i = 0, NQUAD - 1
      theta = (i + 0.5_dp) * PI / NQUAD
      call schmidt_legendre(LMAX, cos(theta), p)
      do l2 = m, LMAX
        gram(m:LMAX, l2) = gram(m:LMAX, l2) + p(m:LMAX, m) * p(l2, m) * sin(theta) * PI / NQUAD
      end do
    end do
    do l = m, LMAX
      expected = merge(2.0_dp, 4.0_dp, m == 0) / (2 * l + 1)
      write (label, '(a,i0,a,i0)') 'normalisation l = ', l, ', m = ', m
      call check(abs(gram(l, l) / expected - 1) < 1e-5_dp, trim(label))
      do l2 = l + 1, LMAX
        write (label, '(a,i0,a,i0,a,i0)') 'orthogonalité l = ', l, ', ', l2, ', m = ', m
        call check(abs(gram(l, l2)) < 1e-5_dp, trim(label))
      end do
    end do
  end subroutine test_normalization

end program test_legendre
