!> Fonctions de Legendre associées en semi-normalisation de Schmidt, sans phase
!> de Condon-Shortley (convention du géomagnétisme) :
!>   P(l, m) = sqrt((2 - delta_m0) (l - m)! / (l + m)!) P_l^m(x),
!> si bien que l'intégrale sur la sphère de [P(l, m) cos(m phi)]^2 vaut 4 pi / (2l + 1).
!>
!> Récurrences stables : la diagonale P(m, m) de proche en proche en m, puis, à m
!> fixé, P(m+1, m) et la récurrence en l. Les deux étapes travaillent sur un
!> tableau de points x(:) pour servir aussi le solveur par anneaux (sh_rings).
module legendre
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: schmidt_legendre, schmidt_diagonal, schmidt_column

contains

  !> d(i, m) = P(m, m)(x(i)) pour m = 0..lmax.
  pure subroutine schmidt_diagonal(lmax, x, d)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: d(size(x), 0:lmax)
    real(dp) :: s(size(x))
    integer :: m

    s = sqrt(max(0.0_dp, 1 - x * x))
    d(:, 0) = 1
    if (lmax >= 1) d(:, 1) = s
    do m = 2, lmax
      d(:, m) = sqrt((2 * m - 1.0_dp) / (2 * m)) * s * d(:, m - 1)
    end do
  end subroutine schmidt_diagonal

  !> p(i, l) = P(l, m)(x(i)) pour l = m..lmax, à partir de pmm(i) = P(m, m)(x(i)).
  pure subroutine schmidt_column(m, lmax, x, pmm, p)
    integer, intent(in) :: m, lmax
    real(dp), intent(in) :: x(:), pmm(:)
    real(dp), intent(out) :: p(size(x), m:lmax)
    integer :: l

    p(:, m) = pmm
    if (m + 1 <= lmax) p(:, m + 1) = sqrt(2 * m + 1.0_dp) * x * pmm
    do l = m + 2, lmax
      p(:, l) = ((2 * l - 1) * x * p(:, l - 1) - sqrt(real((l - 1)**2 - m**2, dp)) * p(:, l - 2)) &
                / sqrt(real(l**2 - m**2, dp))
    end do
  end subroutine schmidt_column

  !> Remplit p(l, m) pour 0 <= m <= l <= lmax, en x = cos(theta) ; p(l, m) = 0 si m > l.
  pure subroutine schmidt_legendre(lmax, x, p)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: x
    real(dp), intent(out) :: p(0:lmax, 0:lmax)
    real(dp) :: d(1, 0:lmax), column(1, 0:lmax)
    integer :: m

    p = 0
    call schmidt_diagonal(lmax, [x], d)
    do m = 0, lmax
      call schmidt_column(m, lmax, [x], d(:, m), column(:, m:lmax))
      p(m:lmax, m) = column(1, m:lmax)
    end do
  end subroutine schmidt_legendre

end module legendre
