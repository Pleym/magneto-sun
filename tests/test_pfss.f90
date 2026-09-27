!> Tests du modèle PFSS contre des solutions analytiques.
program test_pfss
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use pfss, only: pfss_factor, pfss_coefs
  use sh_fit, only: synthesis
  implicit none
  real(dp), parameter :: RSS = 2.5_dp
  integer :: nfail = 0

  call test_factor_limits()
  call test_dipole_neutral_line_at_equator()
  call test_monopole_is_removed()
  if (nfail > 0) error stop 'test_pfss : échec'
  print '(a)', 'pfss OK'

contains

  subroutine check(is_ok, what)
    logical, intent(in) :: is_ok
    character(*), intent(in) :: what
    if (.not. is_ok) then
      print '(2a)', 'ÉCHEC : ', what
      nfail = nfail + 1
    end if
  end subroutine check

  subroutine test_factor_limits()
    integer :: l
    real(dp) :: expected
    logical :: is_one_at_surface, is_ok_at_rss, is_vacuum_far_from_ss
    is_one_at_surface = .true.
    is_ok_at_rss = .true.
    is_vacuum_far_from_ss = .true.
    do l = 1, 30
      is_one_at_surface = is_one_at_surface .and. abs(pfss_factor(l, 1.0_dp, RSS) - 1) < 1e-14_dp
      expected = (2 * l + 1) * RSS**(-(l + 2)) / ((l + 1) + l * RSS**(-(2 * l + 1)))
      is_ok_at_rss = is_ok_at_rss .and. abs(pfss_factor(l, RSS, RSS) / expected - 1) < 1e-12_dp
      ! Surface source rejetée à l'infini : champ potentiel du vide, en r^-(l+2)
      is_vacuum_far_from_ss = is_vacuum_far_from_ss .and. &
        abs(pfss_factor(l, 2.0_dp, 1e8_dp) / 2.0_dp**(-(l + 2)) - 1) < 1e-10_dp
    end do
    call check(is_one_at_surface, 'facteur = 1 à la surface')
    call check(is_ok_at_rss, 'facteur à Rss = (2l+1) Rss^-(l+2) / ((l+1) + l Rss^-(2l+1))')
    call check(is_vacuum_far_from_ss, 'Rss -> infini : B_r en r^-(l+2)')
  end subroutine test_factor_limits

  subroutine test_dipole_neutral_line_at_equator()
    integer, parameter :: LMAX = 3
    real(dp) :: g(0:LMAX, 0:LMAX), h(0:LMAX, 0:LMAX), gss(0:LMAX, 0:LMAX), hss(0:LMAX, 0:LMAX)
    real(dp) :: br(3)
    g = 0
    h = 0
    g(1, 0) = 1   ! dipole axial : B_r = cos(theta) à la surface
    call pfss_coefs(LMAX, g, h, RSS, RSS, gss, hss)
    call synthesis(LMAX, gss, hss, [0.02_dp, 0.0_dp, -0.02_dp], [1.0_dp, 2.0_dp, 3.0_dp], br)
    call check(br(1) > 0 .and. abs(br(2)) < 1e-15_dp .and. br(3) < 0, &
               'dipôle : ligne neutre à l''équateur sur la surface source')
    call check(abs(br(1) - 0.02_dp * pfss_factor(1, RSS, RSS)) < 1e-15_dp, &
               'dipôle : B_r(Rss) = facteur * cos(theta)')
  end subroutine test_dipole_neutral_line_at_equator

  subroutine test_monopole_is_removed()
    integer, parameter :: LMAX = 2
    real(dp) :: g(0:LMAX, 0:LMAX), h(0:LMAX, 0:LMAX), gr(0:LMAX, 0:LMAX), hr(0:LMAX, 0:LMAX)
    g = 0
    h = 0
    g(0, 0) = 5
    g(2, 1) = 0.7_dp
    h(2, 1) = -0.3_dp
    call pfss_coefs(LMAX, g, h, 1.0_dp, RSS, gr, hr)
    call check(gr(0, 0) == 0, 'monopôle retiré')
    call check(abs(gr(2, 1) - 0.7_dp) < 1e-15_dp .and. abs(hr(2, 1) + 0.3_dp) < 1e-15_dp, &
               'autres coefficients inchangés à la surface')
  end subroutine test_monopole_is_removed

end program test_pfss
