!> Fonctions de Legendre associées en semi-normalisation de Schmidt, sans phase
!> de Condon-Shortley (convention du géomagnétisme) :
!>   P(l, m) = sqrt((2 - delta_m0) (l - m)! / (l + m)!) P_l^m(x),
!> si bien que l'intégrale sur la sphère de [P(l, m) cos(m phi)]^2 vaut 4 pi / (2l + 1).
module legendre
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: schmidt_legendre

contains

  !> Remplit p(l, m) pour 0 <= m <= l <= lmax, en x = cos(theta) ; p(l, m) = 0 si m > l.
  !> Récurrences stables : diagonale P(m, m), puis P(m+1, m), puis récurrence en l.
  pure subroutine schmidt_legendre(lmax, x, p)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: x
    real(dp), intent(out) :: p(0:lmax, 0:lmax)
    real(dp) :: s
    integer :: l, m

    p = 0
    s = sqrt(max(0.0_dp, 1 - x * x))
    p(0, 0) = 1
    if (lmax >= 1) p(1, 1) = s
    do m = 2, lmax
      p(m, m) = sqrt((2 * m - 1.0_dp) / (2 * m)) * s * p(m - 1, m - 1)
    end do
    do m = 0, lmax - 1
      p(m + 1, m) = sqrt(2 * m + 1.0_dp) * x * p(m, m)
      do l = m + 2, lmax
        p(l, m) = ((2 * l - 1) * x * p(l - 1, m) - sqrt(real((l - 1)**2 - m**2, dp)) * p(l - 2, m)) &
                  / sqrt(real(l**2 - m**2, dp))
      end do
    end do
  end subroutine schmidt_legendre

end module legendre
