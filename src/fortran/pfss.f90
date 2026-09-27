!> Modèle PFSS (Potential Field Source Surface), solution analytique mode par mode.
!>
!> Entre la surface (r = 1, en rayons solaires) et la surface source (r = rss), le
!> champ dérive d'un potentiel : B = -grad(Phi), laplacien(Phi) = 0. Pour le degré l,
!> Phi_l(r) = A [r^-(l+1) - rss^-(2l+1) r^l] s'annule en rss (champ radial en rss),
!> et A est fixé par le B_r observé en r = 1. D'où, mode par mode :
!>   B_r(r) / B_r(1) = [(l+1) r^-(l+2) + l rss^-(2l+1) r^(l-1)] / [(l+1) + l rss^-(2l+1)].
module pfss
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: pfss_factor, pfss_coefs

contains

  pure real(dp) function pfss_factor(l, r, rss)
    integer, intent(in) :: l
    real(dp), intent(in) :: r, rss
    real(dp) :: q
    q = rss**(-(2 * l + 1))
    pfss_factor = ((l + 1) * r**(-(l + 2)) + l * q * r**(l - 1)) / ((l + 1) + l * q)
  end function pfss_factor

  !> Coefficients de B_r sur la sphère de rayon r. Le monopôle (l = 0) est retiré :
  !> div B = 0 impose un flux net nul, le terme mesuré n'est qu'un biais de la carte.
  pure subroutine pfss_coefs(lmax, g, h, r, rss, gr, hr)
    integer, intent(in) :: lmax
    real(dp), intent(in) :: g(0:lmax, 0:lmax), h(0:lmax, 0:lmax), r, rss
    real(dp), intent(out) :: gr(0:lmax, 0:lmax), hr(0:lmax, 0:lmax)
    integer :: l
    gr(0, :) = 0
    hr(0, :) = 0
    do l = 1, lmax
      gr(l, :) = g(l, :) * pfss_factor(l, r, rss)
      hr(l, :) = h(l, :) * pfss_factor(l, r, rss)
    end do
  end subroutine pfss_coefs

end module pfss
