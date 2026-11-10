! This code is part of data-driven quadrature
!
! Copyright 2026-, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Gas optics, cloud optics and the RTE solvers fused into one kernel per call.
!   Each column loops over the spectral points one after the other and adds up broadband fluxes,
!   so no field has a spectral dimension. The steps are those of
!     tau_absorption_from_fits and add_tau_rayleigh                  (mo_gas_optics_ddq_kernels)
!     compute_Planck_source and get_layer_number                     (mo_gas_optics_utils)
!     compute_all_from_table, delta_scale_2str_k and increment_*     (cloud optics, mo_optical_props_kernels)
!     lw_solver_noscat_oneangle, sw_dif_and_source and adding        (mo_rte_solver_kernels)
!   with the same variable names. Arrays that the original routines keep per g-point are kept
!   here per column and layer only, and overwritten for every spectral point.
!
! -------------------------------------------------------------------------------------------------
module mo_gas_optics_ddq_fused_kernels
  use mo_rte_kind,             only: wp, wl
  use mo_gas_optics_constants, only: boltzmann_k, lightspeed, planck_h, m_h2o, m_dry, avogad, grav
  implicit none
  private
  public :: get_layer_number, tau_absorption_from_fits_invariants
  public :: lw_solver_noscat_ddq, sw_solver_2stream_ddq

  real(wp), parameter :: pi  = acos(-1._wp)
  real(wp), parameter :: eps = 3.0_wp*tiny(1.0_wp)   ! as in mo_optical_props_kernels

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Number density of dry air molecules, as get_layer_number in mo_gas_optics_utils, on the device
  !
  subroutine get_layer_number(ncol, nlay, vmr_h2o, plev, col_dry)
    integer,                           intent(in ) :: ncol, nlay
    real(wp), dimension(ncol, nlay  ), intent(in ) :: vmr_h2o
    real(wp), dimension(ncol, nlay+1), intent(in ) :: plev
    real(wp), dimension(ncol, nlay  ), intent(out) :: col_dry
    real(wp):: delta_plev, m_air, fact
    integer :: icol, ilev

    !$acc                         parallel loop    gang vector collapse(2) default(present) private(delta_plev, m_air, fact)
    !$omp target teams distribute parallel do simd collapse(2) private(delta_plev, m_air, fact)
    do ilev = 1, nlay
      do icol = 1, ncol
        delta_plev = abs(plev(icol,ilev) - plev(icol,ilev+1))
        fact = 1._wp / (1.+vmr_h2o(icol,ilev))
        m_air = (m_dry + m_h2o * vmr_h2o(icol,ilev)) * fact
        col_dry(icol,ilev) = delta_plev/grav * avogad * fact/m_air
      end do
    end do
  end subroutine get_layer_number
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! The nu-independent quantities hoisted out of the spectral loop in tau_absorption_from_fits,
  !   for the fits that use gas igas_used (0 for gases users haven't supplied, with vmr = 0).
  !   Called once per gas, so that only one volume mixing ratio is needed at a time.
  !
  subroutine tau_absorption_from_fits_invariants(ncol, nlay, nnu, igas_used, &
                  play, tlay, dry_num, vmr,                                  &
                  fax_ngas, fax_num_index, fax_c, fax_p0, fax_S, fax_vmr0,   &
                  xsec_ngas, xsec_num_index,                                 &
                  mtckd_ngas, mtckd_num_index, mtckd_T0, mtckd_p0,           &
                  fax_x, fax_w, fax_rs, xsec_w, mt_logr, mt_sw, mt_fw, inv_2kT)
    integer,  intent(in) :: ncol, nlay, nnu, igas_used
    real(wp), dimension(ncol, nlay), intent(in) :: play, tlay, dry_num, vmr
    integer,  intent(in) :: fax_ngas, fax_num_index(fax_ngas)
    real(wp), intent(in) :: fax_c(0:2, fax_ngas, nnu)
    real(wp), dimension(fax_ngas), intent(in) :: fax_p0, fax_S, fax_vmr0
    integer,  intent(in) :: xsec_ngas, xsec_num_index(xsec_ngas)
    integer,  intent(in) :: mtckd_ngas, mtckd_num_index(mtckd_ngas)
    real(wp), intent(in) :: mtckd_T0, mtckd_p0
    real(wp), dimension(ncol, nlay, fax_ngas),   intent(inout) :: fax_x, fax_w
    real(wp), dimension(fax_ngas, nnu),          intent(inout) :: fax_rs
    real(wp), dimension(ncol, nlay, xsec_ngas),  intent(inout) :: xsec_w
    real(wp), dimension(ncol, nlay, mtckd_ngas), intent(inout) :: mt_logr, mt_sw, mt_fw
    real(wp), dimension(ncol, nlay),             intent(inout) :: inv_2kT
    ! -----------------
    integer  :: igas, icol, ilay, inu
    real(wp) :: pr, tr, wn

    !$acc        data copyin(fax_num_index, fax_c, fax_p0, fax_S, fax_vmr0, xsec_num_index, mtckd_num_index)
    !$omp target data map(to:fax_num_index, fax_c, fax_p0, fax_S, fax_vmr0, xsec_num_index, mtckd_num_index)
    do igas = 1, fax_ngas
      if (fax_num_index(igas) /= igas_used) cycle
      !$acc                         parallel loop    collapse(2) default(present)
      !$omp target teams distribute parallel do simd collapse(2)
      do ilay = 1, nlay
        do icol = 1, ncol
          ! Increase pressure to account for self-broadening
          fax_x (icol, ilay, igas) = play(icol, ilay) / fax_p0(igas)                  &
                                   * (1._wp + vmr(icol, ilay) * fax_S(igas)) &
                                   / (1._wp + fax_vmr0(igas)  * fax_S(igas))
          fax_w (icol, ilay, igas) = vmr(icol, ilay) * dry_num(icol, ilay)   ! Integrated number density [mol/m**2]
        end do
      end do
    end do
    do igas = 1, xsec_ngas
      if (xsec_num_index(igas) /= igas_used) cycle
      !$acc                         parallel loop    collapse(2) default(present)
      !$omp target teams distribute parallel do simd collapse(2)
      do ilay = 1, nlay
        do icol = 1, ncol
          xsec_w(icol, ilay, igas) = vmr(icol, ilay) * dry_num(icol, ilay)
        end do
      end do
    end do
    do igas = 1, mtckd_ngas
      if (mtckd_num_index(igas) /= igas_used) cycle
      !$acc                         parallel loop    collapse(2) default(present) private(pr, tr, wn)
      !$omp target teams distribute parallel do simd collapse(2) private(pr, tr, wn)
      do ilay = 1, nlay
        do icol = 1, ncol
          pr  = play(icol, ilay)/mtckd_p0
          tr  = mtckd_T0/tlay(icol, ilay)
          wn  = vmr(icol, ilay) * dry_num(icol, ilay)
          mt_logr(icol, ilay, igas) = log(tr)
          ! Partial pressures, number density, and the (T0/T) of the foreign term,
          ! folded into the weights so the spectral loop reads two arrays instead of five
          mt_sw  (icol, ilay, igas) = pr *          vmr(icol, ilay)  * wn
          mt_fw  (icol, ilay, igas) = pr * (1._wp - vmr(icol, ilay)) * wn * tr
        end do
      end do
    end do
    if (igas_used == 0) then
      !$acc                         parallel loop    collapse(2) default(present)
      !$omp target teams distribute parallel do simd collapse(2)
      do inu = 1, nnu
        do igas = 1, fax_ngas
          fax_rs(igas, inu) = 1._wp / (1._wp + fax_c(1, igas, inu))
        end do
      end do
      !$acc                         parallel loop    collapse(2) default(present)
      !$omp target teams distribute parallel do simd collapse(2)
      do ilay = 1, nlay
        do icol = 1, ncol
          inv_2kT(icol, ilay) = (planck_h * lightspeed * 100._wp) &
                              / (2._wp * boltzmann_k * tlay(icol, ilay))
        end do
      end do
    end if
    !$acc end data
    !$omp end target data
  end subroutine tau_absorption_from_fits_invariants
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Longwave fluxes without scattering, as lw_solver_noscat_oneangle with broadband output,
  !   for clear sky and, if do_allsky, with absorption by clouds. Gas and cloud optical depths
  !   and Planck sources are computed per spectral point inside the kernel.
  !
  subroutine lw_solver_noscat_ddq(ncol, nlay, nnu, top_at_1, D, weight, nus, dnus,           &
                                  play, tlay, tlev, tsfc, sfc_emis,                          &
                                  fax_ngas, fax_a, fax_b, fax_T0, fax_c, fax_sigma0,         &
                                  xsec_ngas, xsec_p, mtckd_ngas, mtckd_cself, mtckd_cfrgn, mtckd_n, &
                                  fax_x, fax_w, fax_rs, xsec_w, mt_logr, mt_sw, mt_fw, inv_2kT, &
                                  do_clrsky, do_allsky, clwp, ciwp, reliq, reice,           &
                                  liq_nsteps, liq_step_size, radliq_lwr, lut_liq,            &
                                  ice_nsteps, ice_step_size, radice_lwr, lut_ice,            &
                                  trans, source_up,                                          &
                                  flux_up_clrsky, flux_dn_clrsky, flux_up, flux_dn)
    integer,     intent(in) :: ncol, nlay, nnu
    logical(wl), intent(in) :: top_at_1
    real(wp),    intent(in) :: D, weight              ! secant of propagation angle, quadrature weight
    real(wp),    intent(in) :: nus(nnu), dnus(nnu)
    real(wp),    intent(in) :: play(ncol, nlay), tlay(ncol, nlay), tlev(ncol, nlay+1), tsfc(ncol)
    real(wp),    intent(in) :: sfc_emis(nnu, ncol)
    integer,     intent(in) :: fax_ngas, xsec_ngas, mtckd_ngas
    real(wp), dimension(0:2, fax_ngas, nnu),     intent(in) :: fax_a, fax_b, fax_c
    real(wp), dimension(fax_ngas),               intent(in) :: fax_T0
    real(wp), dimension(fax_ngas, nnu),          intent(in) :: fax_sigma0, fax_rs
    real(wp), dimension(0:3, xsec_ngas, nnu),    intent(in) :: xsec_p
    real(wp), dimension(mtckd_ngas, nnu),        intent(in) :: mtckd_cself, mtckd_cfrgn, mtckd_n
    real(wp), dimension(ncol, nlay, fax_ngas),   intent(in) :: fax_x, fax_w
    real(wp), dimension(ncol, nlay, xsec_ngas),  intent(in) :: xsec_w
    real(wp), dimension(ncol, nlay, mtckd_ngas), intent(in) :: mt_logr, mt_sw, mt_fw
    real(wp), dimension(ncol, nlay),             intent(in) :: inv_2kT
    logical(wl), intent(in) :: do_clrsky, do_allsky
    real(wp), dimension(ncol, nlay), intent(in) :: clwp, ciwp, reliq, reice
    integer,     intent(in) :: liq_nsteps, ice_nsteps
    real(wp),    intent(in) :: liq_step_size, radliq_lwr, ice_step_size, radice_lwr
    real(wp),    intent(in) :: lut_liq(liq_nsteps, nnu, 3), lut_ice(ice_nsteps, nnu, 3) ! extinction, ssa, asymmetry
    real(wp), dimension(ncol, nlay, 2), intent(inout) :: trans, source_up      ! clear sky, all sky
    real(wp), dimension(ncol, nlay+1),  intent(inout) :: flux_up_clrsky, flux_dn_clrsky, flux_up, flux_dn
    ! -----------------
    integer  :: icol, ilay, ilev, inu, k, top_level, sfc_level, lay_index, inc_index, trans_index
    real(wp) :: tau, tau_loc, trans_s, source_dn, source_up_s, lay_source, lev_source, levp1_source, lev_inc, sfc_src
    real(wp) :: ltau, ltaussa, ltaussag, itau, itaussa, itaussag
    real(wp) :: radn_dn_clrsky, radn_dn, radn_up_clrsky, radn_up

    if(top_at_1) then
      top_level = 1
      sfc_level = nlay+1
    else
      top_level = nlay+1
      sfc_level = 1
    end if

    !$acc        data copyin(nus, dnus, fax_a, fax_b, fax_T0, fax_c, fax_sigma0, xsec_p) &
    !$acc             copyin(mtckd_cself, mtckd_cfrgn, mtckd_n, lut_liq, lut_ice)
    !$omp target data map(to:nus, dnus, fax_a, fax_b, fax_T0, fax_c, fax_sigma0, xsec_p) &
    !$omp             map(to:mtckd_cself, mtckd_cfrgn, mtckd_n, lut_liq, lut_ice)
    !$acc                         parallel loop    gang vector default(present) &
    !$acc   private(tau, tau_loc, trans_s, source_dn, source_up_s, lay_source, lev_source, levp1_source, lev_inc, sfc_src) &
    !$acc   private(ltau, ltaussa, ltaussag, itau, itaussa, itaussag, radn_dn_clrsky, radn_dn, radn_up_clrsky, radn_up) &
    !$acc   private(lay_index, inc_index, trans_index)
    !$omp target teams distribute parallel do simd &
    !$omp   private(tau, tau_loc, trans_s, source_dn, source_up_s, lay_source, lev_source, levp1_source, lev_inc, sfc_src) &
    !$omp   private(ltau, ltaussa, ltaussag, itau, itaussa, itaussag, radn_dn_clrsky, radn_dn, radn_up_clrsky, radn_up) &
    !$omp   private(lay_index, inc_index, trans_index)
    do icol = 1, ncol
      do ilev = 1, nlay+1
        if(do_clrsky) flux_up_clrsky(icol,ilev) = 0._wp
        if(do_clrsky) flux_dn_clrsky(icol,ilev) = 0._wp
        if(do_allsky) flux_up       (icol,ilev) = 0._wp
        if(do_allsky) flux_dn       (icol,ilev) = 0._wp
      end do
      do inu = 1, nnu
        !
        ! Transport down from the top without incident flux, keeping transmissivity and upward source.
        !   Broadband fluxes are summed over spectral points as they are computed.
        !
        radn_dn_clrsky = 0._wp
        radn_dn        = 0._wp
        lev_inc = B_nu(tlev(icol,top_level), nus(inu)) * dnus(inu)
        do k = 1, nlay
          if(top_at_1) then
            lay_index   = k
            inc_index   = lay_index
            trans_index = lay_index+1
          else
            lay_index   = nlay-k+1
            inc_index   = lay_index+1
            trans_index = lay_index
          end if
          lay_source = B_nu(tlay(icol,lay_index),   nus(inu)) * dnus(inu)
          if(top_at_1) then
            lev_source   = lev_inc
            levp1_source = B_nu(tlev(icol,trans_index), nus(inu)) * dnus(inu)
            lev_inc      = levp1_source
          else
            levp1_source = lev_inc
            lev_source   = B_nu(tlev(icol,trans_index), nus(inu)) * dnus(inu)
            lev_inc      = lev_source
          end if
          tau = tau_absorption_from_fits(icol, lay_index, inu, ncol, nlay, nnu, nus(inu), play, tlay, &
                  fax_ngas, fax_a, fax_b, fax_T0, fax_c, fax_sigma0, fax_rs, fax_x, fax_w,           &
                  xsec_ngas, xsec_p, xsec_w, mtckd_ngas, mtckd_cself, mtckd_cfrgn, mtckd_n,          &
                  mt_logr, mt_sw, mt_fw, inv_2kT)
          if(do_clrsky) then
            tau_loc = tau*D
            trans_s = exp(-tau_loc)
            call lw_source_noscat(top_at_1, lay_source, lev_source, levp1_source, &
                                  tau_loc, trans_s, source_dn, source_up_s)
            trans    (icol,lay_index,1) = trans_s
            source_up(icol,lay_index,1) = source_up_s
            radn_dn_clrsky = trans_s*radn_dn_clrsky + source_dn
            flux_dn_clrsky(icol,trans_index) = flux_dn_clrsky(icol,trans_index) + radn_dn_clrsky
          end if
          if(do_allsky) then
            !
            ! Absorption optical depth of clouds = (1-ssa) * tau = tau - taussa, as in the cloud optics
            !
            call compute_all_from_table(liq_nsteps, nnu, inu, liq_step_size, radliq_lwr, lut_liq, &
                                        clwp(icol,lay_index), reliq(icol,lay_index), ltau, ltaussa, ltaussag)
            call compute_all_from_table(ice_nsteps, nnu, inu, ice_step_size, radice_lwr, lut_ice, &
                                        ciwp(icol,lay_index), reice(icol,lay_index), itau, itaussa, itaussag)
            tau_loc = (tau + ((ltau - ltaussa) + (itau - itaussa)))*D
            trans_s = exp(-tau_loc)
            call lw_source_noscat(top_at_1, lay_source, lev_source, levp1_source, &
                                  tau_loc, trans_s, source_dn, source_up_s)
            trans    (icol,lay_index,2) = trans_s
            source_up(icol,lay_index,2) = source_up_s
            radn_dn = trans_s*radn_dn + source_dn
            flux_dn(icol,trans_index) = flux_dn(icol,trans_index) + radn_dn
          end if
        end do
        !
        ! Surface reflection and emission, then transport up
        !
        sfc_src = B_nu(tsfc(icol), nus(inu)) * dnus(inu)
        radn_up_clrsky = radn_dn_clrsky*(1._wp - sfc_emis(inu,icol)) + sfc_src*sfc_emis(inu,icol)
        radn_up        = radn_dn       *(1._wp - sfc_emis(inu,icol)) + sfc_src*sfc_emis(inu,icol)
        if(do_clrsky) flux_up_clrsky(icol,sfc_level) = flux_up_clrsky(icol,sfc_level) + radn_up_clrsky
        if(do_allsky) flux_up       (icol,sfc_level) = flux_up       (icol,sfc_level) + radn_up
        do k = nlay, 1, -1
          if(top_at_1) then
            lay_index = k
            inc_index = lay_index
          else
            lay_index = nlay-k+1
            inc_index = lay_index+1
          end if
          if(do_clrsky) then
            radn_up_clrsky = trans(icol,lay_index,1)*radn_up_clrsky + source_up(icol,lay_index,1)
            flux_up_clrsky(icol,inc_index) = flux_up_clrsky(icol,inc_index) + radn_up_clrsky
          end if
          if(do_allsky) then
            radn_up = trans(icol,lay_index,2)*radn_up + source_up(icol,lay_index,2)
            flux_up(icol,inc_index) = flux_up(icol,inc_index) + radn_up
          end if
        end do
      end do
      !
      ! Conversion from intensity to flux assuming azimuthal isotropy and quadrature weight
      !
      do ilev = 1, nlay+1
        if(do_clrsky) flux_up_clrsky(icol,ilev) = pi * weight * flux_up_clrsky(icol,ilev)
        if(do_clrsky) flux_dn_clrsky(icol,ilev) = pi * weight * flux_dn_clrsky(icol,ilev)
        if(do_allsky) flux_up       (icol,ilev) = pi * weight * flux_up       (icol,ilev)
        if(do_allsky) flux_dn       (icol,ilev) = pi * weight * flux_dn       (icol,ilev)
      end do
    end do
    !$acc end data
    !$omp end target data
  end subroutine lw_solver_noscat_ddq
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Shortwave fluxes, as sw_solver_2stream with broadband output, for clear sky and, if do_allsky,
  !   with clouds. flux_dn is total (direct and diffuse); flux_dir, if do_dir, is the direct beam. Clear and all sky are solved one after the
  !   other for each spectral point, sharing the gas optical properties kept in tau and ssa.
  !   Optionally, surface fluxes weighted by sfc_frac per spectral point, for the last sky solved.
  !
  subroutine sw_solver_2stream_ddq(ncol, nlay, nnu, top_at_1, nus, dnus, solar_source, mu0, toa_factor, &
                                   sfc_alb_dir, sfc_alb_dif, play, tlay, dry_num,                      &
                                   fax_ngas, fax_a, fax_b, fax_T0, fax_c, fax_sigma0,                  &
                                   xsec_ngas, xsec_p, mtckd_ngas, mtckd_cself, mtckd_cfrgn, mtckd_n,   &
                                   fax_x, fax_w, fax_rs, xsec_w, mt_logr, mt_sw, mt_fw, inv_2kT,       &
                                   has_rayleigh, rayleigh_xsec,                                        &
                                   do_clrsky, do_allsky, clwp, ciwp, reliq, reice,                    &
                                   liq_nsteps, liq_step_size, radliq_lwr, lut_liq,                     &
                                   ice_nsteps, ice_step_size, radice_lwr, lut_ice,                     &
                                   nfrac, sfc_frac, sfc_dn_dir, sfc_dn_dif, sfc_up,                    &
                                   tau, ssa, Rdif, Tdif, source_dn, source_up, albedo, src,           &
                                   flux_up_clrsky, flux_dn_clrsky, flux_up, flux_dn, do_dir, flux_dir)
    integer,     intent(in) :: ncol, nlay, nnu
    logical(wl), intent(in) :: top_at_1
    real(wp),    intent(in) :: nus(nnu), dnus(nnu), solar_source(nnu)
    real(wp),    intent(in) :: mu0(ncol), toa_factor(ncol)
    real(wp),    intent(in) :: sfc_alb_dir(nnu, ncol), sfc_alb_dif(nnu, ncol)
    real(wp),    intent(in) :: play(ncol, nlay), tlay(ncol, nlay), dry_num(ncol, nlay)
    integer,     intent(in) :: fax_ngas, xsec_ngas, mtckd_ngas
    real(wp), dimension(0:2, fax_ngas, nnu),     intent(in) :: fax_a, fax_b, fax_c
    real(wp), dimension(fax_ngas),               intent(in) :: fax_T0
    real(wp), dimension(fax_ngas, nnu),          intent(in) :: fax_sigma0, fax_rs
    real(wp), dimension(0:3, xsec_ngas, nnu),    intent(in) :: xsec_p
    real(wp), dimension(mtckd_ngas, nnu),        intent(in) :: mtckd_cself, mtckd_cfrgn, mtckd_n
    real(wp), dimension(ncol, nlay, fax_ngas),   intent(in) :: fax_x, fax_w
    real(wp), dimension(ncol, nlay, xsec_ngas),  intent(in) :: xsec_w
    real(wp), dimension(ncol, nlay, mtckd_ngas), intent(in) :: mt_logr, mt_sw, mt_fw
    real(wp), dimension(ncol, nlay),             intent(in) :: inv_2kT
    logical(wl), intent(in) :: has_rayleigh
    real(wp),    intent(in) :: rayleigh_xsec(nnu)
    logical(wl), intent(in) :: do_clrsky, do_allsky
    real(wp), dimension(ncol, nlay), intent(in) :: clwp, ciwp, reliq, reice
    integer,     intent(in) :: liq_nsteps, ice_nsteps
    real(wp),    intent(in) :: liq_step_size, radliq_lwr, ice_step_size, radice_lwr
    real(wp),    intent(in) :: lut_liq(liq_nsteps, nnu, 3), lut_ice(ice_nsteps, nnu, 3)
    integer,     intent(in) :: nfrac
    real(wp),    intent(in) :: sfc_frac(nnu, max(nfrac,1))
    real(wp), dimension(ncol, max(nfrac,1)), intent(inout) :: sfc_dn_dir, sfc_dn_dif, sfc_up
    real(wp), dimension(ncol, nlay  ), intent(inout) :: tau, ssa, Rdif, Tdif, source_dn, source_up
    real(wp), dimension(ncol, nlay+1), intent(inout) :: albedo, src
    real(wp), dimension(ncol, nlay+1), intent(inout) :: flux_up_clrsky, flux_dn_clrsky, flux_up, flux_dn
    logical(wl), intent(in) :: do_dir                     ! direct beam flux, for the last sky solved
    real(wp), dimension(ncol, nlay+1), intent(inout) :: flux_dir
    ! -----------------
    integer  :: icol, ilay, ilev, inu, k, isky, top_level, sfc_level, lay_index, inc_index, trans_index, ifrac
    real(wp) :: tau_abs, t_r, tau_s, w0_s, g_s, inc_flux, flux_dn_dir, flux_dif, flux_upw, denom
    real(wp) :: Rdif_s, Tdif_s, source_dn_s, source_up_s
    real(wp) :: ltau, ltaussa, ltaussag, itau, itaussa, itaussag, tau2, ssa2, g2, taussa, f, wf, tau12, tauscat12

    if(top_at_1) then
      top_level = 1
      sfc_level = nlay+1
    else
      top_level = nlay+1
      sfc_level = 1
    end if

    !$acc        data copyin(nus, dnus, solar_source, fax_a, fax_b, fax_T0, fax_c, fax_sigma0, xsec_p) &
    !$acc             copyin(mtckd_cself, mtckd_cfrgn, mtckd_n, rayleigh_xsec, lut_liq, lut_ice, sfc_frac)
    !$omp target data map(to:nus, dnus, solar_source, fax_a, fax_b, fax_T0, fax_c, fax_sigma0, xsec_p) &
    !$omp             map(to:mtckd_cself, mtckd_cfrgn, mtckd_n, rayleigh_xsec, lut_liq, lut_ice, sfc_frac)
    !$acc                         parallel loop    gang vector default(present) &
    !$acc   private(tau_abs, t_r, tau_s, w0_s, g_s, inc_flux, flux_dn_dir, flux_dif, flux_upw, denom) &
    !$acc   private(Rdif_s, Tdif_s, source_dn_s, source_up_s) &
    !$acc   private(ltau, ltaussa, ltaussag, itau, itaussa, itaussag, tau2, ssa2, g2, taussa, f, wf, tau12, tauscat12) &
    !$acc   private(lay_index, inc_index, trans_index)
    !$omp target teams distribute parallel do simd &
    !$omp   private(tau_abs, t_r, tau_s, w0_s, g_s, inc_flux, flux_dn_dir, flux_dif, flux_upw, denom) &
    !$omp   private(Rdif_s, Tdif_s, source_dn_s, source_up_s) &
    !$omp   private(ltau, ltaussa, ltaussag, itau, itaussa, itaussag, tau2, ssa2, g2, taussa, f, wf, tau12, tauscat12) &
    !$omp   private(lay_index, inc_index, trans_index)
    do icol = 1, ncol
      do ilev = 1, nlay+1
        if(do_clrsky) flux_up_clrsky(icol,ilev) = 0._wp
        if(do_clrsky) flux_dn_clrsky(icol,ilev) = 0._wp
        if(do_allsky) flux_up       (icol,ilev) = 0._wp
        if(do_allsky) flux_dn       (icol,ilev) = 0._wp
        if(do_dir)    flux_dir      (icol,ilev) = 0._wp
      end do
      do ifrac = 1, nfrac
        sfc_dn_dir(icol,ifrac) = 0._wp
        sfc_dn_dif(icol,ifrac) = 0._wp
        sfc_up    (icol,ifrac) = 0._wp
      end do
      do inu = 1, nnu
        !
        ! Gas optical properties, as tau_absorption_from_fits and add_tau_rayleigh
        !
        do ilay = 1, nlay
          tau_abs = tau_absorption_from_fits(icol, ilay, inu, ncol, nlay, nnu, nus(inu), play, tlay, &
                      fax_ngas, fax_a, fax_b, fax_T0, fax_c, fax_sigma0, fax_rs, fax_x, fax_w,      &
                      xsec_ngas, xsec_p, xsec_w, mtckd_ngas, mtckd_cself, mtckd_cfrgn, mtckd_n,     &
                      mt_logr, mt_sw, mt_fw, inv_2kT)
          if(has_rayleigh) then
            t_r = dry_num(icol,ilay) * rayleigh_xsec(inu)
            tau(icol,ilay) = tau_abs + t_r
            ssa(icol,ilay) = t_r/(tau_abs + t_r)
          else
            tau(icol,ilay) = tau_abs
            ssa(icol,ilay) = 0._wp
          end if
        end do
        do isky = 1, 2
          if(isky == 1 .and. .not. do_clrsky) cycle
          if(isky == 2 .and. .not. do_allsky) cycle
          !
          ! Direct beam, layer reflectance and transmittance and sources, as sw_dif_and_source
          !
          inc_flux = solar_source(inu) * dnus(inu) * toa_factor(icol)
          flux_dn_dir = inc_flux * mu0(icol)
          if(isky == 1) flux_dn_clrsky(icol,top_level) = flux_dn_clrsky(icol,top_level) + flux_dn_dir
          if(isky == 2) flux_dn       (icol,top_level) = flux_dn       (icol,top_level) + flux_dn_dir
          if(do_dir .and. (isky == 2 .or. .not. do_allsky)) &
            flux_dir(icol,top_level) = flux_dir(icol,top_level) + flux_dn_dir
          do k = 1, nlay
            if(top_at_1) then
              lay_index   = k
              inc_index   = lay_index
              trans_index = lay_index+1
            else
              lay_index   = nlay-k+1
              inc_index   = lay_index+1
              trans_index = lay_index
            end if
            tau_s = tau(icol,lay_index)
            w0_s  = ssa(icol,lay_index)
            g_s   = 0._wp
            if(isky == 2) then
              !
              ! Cloud optical properties as in the cloud optics, delta_scale_2str_k,
              !   and increment_2stream_by_2stream of the gas (1) by the clouds (2)
              !
              call compute_all_from_table(liq_nsteps, nnu, inu, liq_step_size, radliq_lwr, lut_liq, &
                                          clwp(icol,lay_index), reliq(icol,lay_index), ltau, ltaussa, ltaussag)
              call compute_all_from_table(ice_nsteps, nnu, inu, ice_step_size, radice_lwr, lut_ice, &
                                          ciwp(icol,lay_index), reice(icol,lay_index), itau, itaussa, itaussag)
              tau2   = ltau    + itau
              taussa = ltaussa + itaussa
              g2     = (ltaussag + itaussag) / max(epsilon(tau2), taussa)
              ssa2   = taussa / max(epsilon(tau2), tau2)
              if(tau2 > eps) then
                f  = g2  * g2
                wf = ssa2 * f
                tau2 = (1._wp - wf) * tau2
                ssa2 = (ssa2 - wf) / (1.0_wp - wf)
                g2   = (g2   -  f) / (1.0_wp -  f)
              end if
              tau12     = tau_s + tau2
              tauscat12 = tau_s * w0_s + tau2 * ssa2
              if(tauscat12 > eps) then
                g_s   = (tau_s * w0_s * g_s + tau2 * ssa2 * g2) / tauscat12
                w0_s  = tauscat12 / tau12
                tau_s = tau12
              end if
            end if
            call sw_dif_and_source(mu0(icol), tau_s, w0_s, g_s, flux_dn_dir, Rdif_s, Tdif_s, source_dn_s, source_up_s)
            Rdif     (icol,lay_index) = Rdif_s
            Tdif     (icol,lay_index) = Tdif_s
            source_dn(icol,lay_index) = source_dn_s
            source_up(icol,lay_index) = source_up_s
            if(isky == 1) flux_dn_clrsky(icol,trans_index) = flux_dn_clrsky(icol,trans_index) + flux_dn_dir
            if(isky == 2) flux_dn       (icol,trans_index) = flux_dn       (icol,trans_index) + flux_dn_dir
            if(do_dir .and. (isky == 2 .or. .not. do_allsky)) &
              flux_dir(icol,trans_index) = flux_dir(icol,trans_index) + flux_dn_dir
          end do
          !
          ! Transport of diffuse radiation, as adding: albedo and source from the surface up, ...
          !
          albedo(icol,sfc_level) = sfc_alb_dif(inu,icol)
          src   (icol,sfc_level) = flux_dn_dir*sfc_alb_dir(inu,icol)
          do k = nlay, 1, -1
            if(top_at_1) then
              lay_index   = k
              inc_index   = lay_index
              trans_index = lay_index+1
            else
              lay_index   = nlay-k+1
              inc_index   = lay_index+1
              trans_index = lay_index
            end if
            denom = 1._wp/(1._wp - Rdif(icol,lay_index)*albedo(icol,trans_index))                 ! Eq 10
            albedo(icol,inc_index) = Rdif(icol,lay_index) + &
                  Tdif(icol,lay_index)*Tdif(icol,lay_index) * albedo(icol,trans_index) * denom    ! Equation 9
            src(icol,inc_index) = source_up(icol,lay_index) + &
                  Tdif(icol,lay_index) * denom *                                                 &
                  (src(icol,trans_index) + albedo(icol,trans_index)*source_dn(icol,lay_index))   ! Equation 11
          end do
          !
          ! ... then fluxes from the top down, with denom computed again instead of stored
          !
          flux_dif = 0._wp
          flux_upw = flux_dif * albedo(icol,top_level) + src(icol,top_level)                     ! Eq 12
          if(isky == 1) flux_up_clrsky(icol,top_level) = flux_up_clrsky(icol,top_level) + flux_upw
          if(isky == 2) flux_up       (icol,top_level) = flux_up       (icol,top_level) + flux_upw
          do k = 1, nlay
            if(top_at_1) then
              lay_index   = k
              inc_index   = lay_index
              trans_index = lay_index+1
            else
              lay_index   = nlay-k+1
              inc_index   = lay_index+1
              trans_index = lay_index
            end if
            denom = 1._wp/(1._wp - Rdif(icol,lay_index)*albedo(icol,trans_index))
            flux_dif = (Tdif(icol,lay_index)*flux_dif +                                          &  ! Equation 13
                        Rdif(icol,lay_index)*src(icol,trans_index) + source_dn(icol,lay_index)) * denom
            flux_upw = flux_dif * albedo(icol,trans_index) + src(icol,trans_index)               ! Equation 12
            if(isky == 1) then
              flux_up_clrsky(icol,trans_index) = flux_up_clrsky(icol,trans_index) + flux_upw
              flux_dn_clrsky(icol,trans_index) = flux_dn_clrsky(icol,trans_index) + flux_dif
            else
              flux_up       (icol,trans_index) = flux_up       (icol,trans_index) + flux_upw
              flux_dn       (icol,trans_index) = flux_dn       (icol,trans_index) + flux_dif
            end if
          end do
          !
          ! Surface fluxes by spectral range, for the last sky solved
          !
          if(isky == 2 .or. .not. do_allsky) then
            do ifrac = 1, nfrac
              sfc_dn_dir(icol,ifrac) = sfc_dn_dir(icol,ifrac) + sfc_frac(inu,ifrac) * flux_dn_dir
              sfc_dn_dif(icol,ifrac) = sfc_dn_dif(icol,ifrac) + sfc_frac(inu,ifrac) * flux_dif
              sfc_up    (icol,ifrac) = sfc_up    (icol,ifrac) + sfc_frac(inu,ifrac) * flux_upw
            end do
          end if
        end do
      end do
    end do
    !$acc end data
    !$omp end target data
  end subroutine sw_solver_2stream_ddq
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Absorption optical depth of one layer at one spectral point, the spectral loop of tau_absorption_from_fits
  !
  pure function tau_absorption_from_fits(icol, ilay, inu, ncol, nlay, nnu, nu, play, tlay,   &
                  fax_ngas, fax_a, fax_b, fax_T0, fax_c, fax_sigma0, fax_rs, fax_x, fax_w,   &
                  xsec_ngas, xsec_p, xsec_w, mtckd_ngas, mtckd_cself, mtckd_cfrgn, mtckd_n,  &
                  mt_logr, mt_sw, mt_fw, inv_2kT) result(tau)
    !$acc routine seq
    !$omp declare target
    integer,  intent(in) :: icol, ilay, inu, ncol, nlay, nnu
    real(wp), intent(in) :: nu
    real(wp), dimension(ncol, nlay), intent(in) :: play, tlay
    integer,  intent(in) :: fax_ngas, xsec_ngas, mtckd_ngas
    real(wp), dimension(0:2, fax_ngas, nnu),     intent(in) :: fax_a, fax_b, fax_c
    real(wp), dimension(fax_ngas),               intent(in) :: fax_T0
    real(wp), dimension(fax_ngas, nnu),          intent(in) :: fax_sigma0, fax_rs
    real(wp), dimension(ncol, nlay, fax_ngas),   intent(in) :: fax_x, fax_w
    real(wp), dimension(0:3, xsec_ngas, nnu),    intent(in) :: xsec_p
    real(wp), dimension(ncol, nlay, xsec_ngas),  intent(in) :: xsec_w
    real(wp), dimension(mtckd_ngas, nnu),        intent(in) :: mtckd_cself, mtckd_cfrgn, mtckd_n
    real(wp), dimension(ncol, nlay, mtckd_ngas), intent(in) :: mt_logr, mt_sw, mt_fw
    real(wp), dimension(ncol, nlay),             intent(in) :: inv_2kT
    real(wp) :: tau
    ! -----------------
    integer  :: igas
    real(wp) :: x, w, q, num, den, delta_T, c0, lin, rs, acc
    real(wp) :: cself, cfrgn, R ! MT_CKD

    acc = 0._wp
    !
    ! Functional approximation to cross-sections
    !
    do igas = 1, fax_ngas
      x       = fax_x (icol, ilay, igas)
      delta_T = tlay(icol, ilay) - fax_T0(igas)
      c0      = fax_c(0, igas, inu)
      lin     = fax_c(2, igas, inu)
      rs      = fax_rs(igas, inu)
      w       = rs * x + (1._wp - rs)
      ! c0/w + (1-c0)*w = q/w
      q       = c0 + (1._wp - c0) * w * w
      ! P_scale*T_scale = [lin*x*q + (1-lin)*w] * num / (q * den)
      num     = (lin * x * q + (1._wp - lin) * w) &
              * (fax_a(0, igas, inu) + delta_T * (fax_a(1, igas, inu) + fax_a(2, igas, inu) * delta_T))
      den     = q &
              * (fax_b(0, igas, inu) + delta_T * (fax_b(1, igas, inu) + fax_b(2, igas, inu) * delta_T))
      acc = acc &
        + fax_sigma0(igas, inu) * max(num / den, 0._wp) & ! cross-section [m**2/mol]
        * fax_w(icol, ilay, igas)
    end do
    !
    ! Cross-sections pressure and temperature dependence following doi:10.1029/2022MS003239
    !
    do igas = 1, xsec_ngas
      acc = acc &
        + (xsec_p(0, igas, inu)                                                         &
           + tlay(icol, ilay) * (xsec_p(1, igas, inu) + xsec_p(2, igas, inu) * tlay(icol, ilay)) &
           + xsec_p(3, igas, inu) * play(icol, ilay))                                   &
        * xsec_w(icol, ilay, igas)
    end do
    !
    ! MT_CKD continuum
    !
    do igas = 1, mtckd_ngas
      ! (T0/T)**(1+n) written as exp((1+n)*log(T0/T)) with the log hoisted
      cself = exp((1._wp + mtckd_n(igas, inu)) * mt_logr(icol, ilay, igas)) &
            * mt_sw(icol, ilay, igas) * mtckd_cself(igas, inu)
      cfrgn = mt_fw(icol, ilay, igas) * mtckd_cfrgn(igas, inu)
      R = nu * tanh(nu * inv_2kT(icol, ilay))
      acc = acc + R * (cself + cfrgn)
    end do
    tau = max(0._wp, acc)
  end function tau_absorption_from_fits
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Cloud optical depth, tau*ssa and tau*ssa*g at one spectral point, as compute_all_from_table
  !   in the cloud optics of ICON; 0 without condensate
  !
  pure subroutine compute_all_from_table(nsteps, nbnd, ibnd, step_size, offset, lut, lwp, re, &
                                         tau, taussa, taussag)
    !$acc routine seq
    !$omp declare target
    integer,  intent(in)  :: nsteps, nbnd, ibnd
    real(wp), intent(in)  :: step_size, offset
    real(wp), intent(in)  :: lut(nsteps, nbnd, 3)      ! tau_table, ssa_table, asy_table
    real(wp), intent(in)  :: lwp, re
    real(wp), intent(out) :: tau, taussa, taussag
    integer  :: index
    real(wp) :: fint, t, ts

    if(lwp > 0._wp) then
      index = min(floor((re - offset)/step_size)+1, nsteps-1)
      fint = (re - offset)/step_size - (index-1)
      t   = lwp * &
            (lut(index,  ibnd,1) + fint * (lut(index+1,ibnd,1) - lut(index,ibnd,1)))
      ts  = t   * &
            (lut(index,  ibnd,2) + fint * (lut(index+1,ibnd,2) - lut(index,ibnd,2)))
      taussag = ts * &
            (lut(index,  ibnd,3) + fint * (lut(index+1,ibnd,3) - lut(index,ibnd,3)))
      taussa  = ts
      tau     = t
    else
      tau     = 0._wp
      taussa  = 0._wp
      taussag = 0._wp
    end if
  end subroutine compute_all_from_table
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! As lw_source_noscat in mo_rte_solver_kernels
  !
  pure subroutine lw_source_noscat(top_at_1, lay_source, lev_source, levp1_source, tau, trans, &
                                   source_dn, source_up)
    !$acc routine seq
    !$omp declare target
    logical(wl), intent(in) :: top_at_1
    real(wp),    intent(in) :: lay_source,   & ! Planck source at layer center
                               lev_source,   & ! Planck source at levels (layer edges),
                               levp1_source, & ! Planck source at level +1 (layer edges),
                               tau,          & ! Optical path (tau/mu)
                               trans           ! Transmissivity (exp(-tau))
    real(wp), intent(inout):: source_dn, source_up
    real(wp), parameter  :: tau_thresh = sqrt(sqrt(epsilon(tau)))
    real(wp)             :: fact, source_inc, source_dec

    if(tau > tau_thresh) then
      fact = (1._wp - trans)/tau - trans
    else
      fact = tau * (0.5_wp + tau * (-1._wp/3._wp + tau * 1._wp/8._wp) )
    end if
    source_inc = (1._wp - trans) * levp1_source + 2._wp * fact * (lay_source - levp1_source)
    source_dec = (1._wp - trans) * lev_source   + 2._wp * fact * (lay_source - lev_source)
    if (top_at_1) then
      source_dn = source_inc
      source_up = source_dec
    else
      source_up = source_inc
      source_dn = source_dec
    end if
  end subroutine lw_source_noscat
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! One layer of sw_dif_and_source in mo_rte_solver_kernels. flux_dn_dir is the incident direct
  !   flux on entry and the transmitted one on return.
  !
  pure subroutine sw_dif_and_source(mu0_s, tau_s, w0_s, g_s, flux_dn_dir, Rdif, Tdif, source_dn, source_up)
    !$acc routine seq
    !$omp declare target
    real(wp), intent(in   ) :: mu0_s, tau_s, w0_s, g_s
    real(wp), intent(inout) :: flux_dn_dir
    real(wp), intent(  out) :: Rdif, Tdif, source_dn, source_up
    ! Variables used in Meador and Weaver
    real(wp) :: gamma1, gamma2, gamma3, gamma4, alpha1, alpha2
    ! Ancillary variables
    real(wp), parameter :: min_k = 1.e4_wp * epsilon(1._wp) ! Suggestion from Chiel van Heerwaarden
    real(wp) :: k, exp_minusktau, k_mu, k_gamma3, k_gamma4
    real(wp) :: RT_term, exp_minus2ktau
    real(wp) :: Rdir, Tdir, Tnoscat, inc_flux

    inc_flux = flux_dn_dir
    ! Zdunkowski Practical Improved Flux Method "PIFM"
    !  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    gamma1 = (8._wp - w0_s * (5._wp + 3._wp * g_s)) * .25_wp
    gamma2 =  3._wp *(w0_s * (1._wp -         g_s)) * .25_wp
    k = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), min_k))
    k_mu     = k * mu0_s
    exp_minusktau = exp(-tau_s*k)
    exp_minus2ktau = exp_minusktau * exp_minusktau
    RT_term = 1._wp / (k      * (1._wp + exp_minus2ktau)  + &
                       gamma1 * (1._wp - exp_minus2ktau) )
    ! Equation 25
    Rdif = RT_term * gamma2 * (1._wp - exp_minus2ktau)
    ! Equation 26
    Tdif = RT_term * 2._wp * k * exp_minusktau
    if(mu0_s > 0._wp) then
      ! Equation 14, multiplying top and bottom by exp(-k*tau)
      !   and rearranging to avoid div by 0.
      RT_term =  w0_s * RT_term/merge(1._wp - k_mu*k_mu, &
                                      epsilon(1._wp),    &
                                      abs(1._wp - k_mu*k_mu) >= epsilon(1._wp))
      gamma3 = (2._wp - 3._wp * mu0_s *         g_s ) * .25_wp
      gamma4 =  1._wp - gamma3
      alpha1 = gamma1 * gamma4 + gamma2 * gamma3           ! Eq. 16
      alpha2 = gamma1 * gamma3 + gamma2 * gamma4           ! Eq. 17
      ! Transmittance of direct, unscattered beam.
      k_gamma3 = k * gamma3
      k_gamma4 = k * gamma4
      Tnoscat = exp(-tau_s/mu0_s)
      Rdir = RT_term  *                                            &
          ((1._wp - k_mu) * (alpha2 + k_gamma3)                  - &
           (1._wp + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau - &
           2.0_wp * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau * Tnoscat)
      Tdir = -RT_term *                                                             &
            ((1._wp + k_mu) * (alpha1 + k_gamma4)                  * Tnoscat - &
             (1._wp - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * Tnoscat - &
             2.0_wp * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau)
      ! Final check that energy is not spuriously created
      Rdir    = max(0.0_wp, min(Rdir, (1.0_wp - Tnoscat       ) ))
      Tdir    = max(0.0_wp, min(Tdir, (1.0_wp - Tnoscat - Rdir) ))
      source_up   =    Rdir * inc_flux
      source_dn   =    Tdir * inc_flux
      flux_dn_dir = Tnoscat * inc_flux
    else
      source_up   = 0._wp
      source_dn   = 0._wp
      flux_dn_dir = 0._wp
    end if
  end subroutine sw_dif_and_source
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Planck function, as B_nu in mo_gas_optics_utils
  !
  pure function B_nu(T, nu)
    !$acc routine seq
    !$omp declare target
    real(wp), intent(in) :: T, nu
    real(wp)             :: B_nu
    B_nu = 100._wp*2._wp*planck_h*((nu*100._wp)**3)*(lightspeed**2) / &
         (exp((planck_h * lightspeed * nu * 100._wp) / (boltzmann_k * T)) - 1._wp)
  end function B_nu
end module mo_gas_optics_ddq_fused_kernels
