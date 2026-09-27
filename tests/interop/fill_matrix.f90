!> Remplit A(i, j) = 10*i + j. Sert uniquement à tester le contrat
!> d'interface Fortran <-> C++ (voir test_interop.cpp).
subroutine fill_matrix(m, n, a) bind(C, name="fill_matrix")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none
  integer(c_int), value, intent(in) :: m, n
  real(c_double), intent(out) :: a(m, n)
  integer :: i, j

  do j = 1, n
    do i = 1, m
      a(i, j) = 10.0_c_double * i + j
    end do
  end do
end subroutine fill_matrix
