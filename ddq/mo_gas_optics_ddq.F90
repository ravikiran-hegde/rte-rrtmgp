module mo_gas_optics_ddq
! This code implements data-driven quadrature
!
! Copyright 2026-, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> ## Class implementing data-driven quadrature
!>
!> Implements a class for computing spectrally-resolved gas optical properties and source functions
!>   given atmopsheric physical properties (profiles of temperature, pressure, and gas concentrations)
!>   The class must be initialized with data (provided as a netCDF file) before being used.
!>
!> Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
!>   external stellar radiation (shortwave radiation in the Earth's atmosphere).
!>   The variant is chosen based on whether a solar_source is supplied during initialization.
!
! -------------------------------------------------------------------------------------------------module gas_optics_ddq
  use mo_rte_kind,           only: wp, wl
  use mo_rte_config,         only: check_extents, check_values
  use mo_rte_util_array,     only : zero_array
  use mo_rte_util_array_validation, &
                             only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props,      only: ty_optical_props
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_optics_constants,   only: boltzmann_k, lightspeed, planck_h
  use mo_gas_optics_util_string, only: lower_case, string_in_array, string_loc_in_array
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_arry, &
                                   ty_optical_props_1scl, &
                                   ty_optical_props_2str, &
                                   ty_optical_props_nstr
  use mo_gas_optics,         only: ty_gas_optics
  use mo_gas_optics_utils,   only: compute_Planck_source, get_layer_number, interp_tlev_from_tlay
  use mo_gas_optics_ddq_kernels, &
                             only: tau_absorption_from_fits, add_tau_rayleigh, &
                                   fax_norder, fax_nterms, xsec_nterms
  use mo_gas_optics_ddq_fused_kernels, &
                             only: get_layer_number_ddq => get_layer_number, tau_absorption_from_fits_invariants, &
                                   lw_solver_noscat_ddq, sw_solver_2stream_ddq
  implicit none
  private

  integer, parameter :: gas_name_len = 8
  ! -------------------------------------------------------------------------------------------------
  type, extends(ty_gas_optics), public :: ty_gas_optics_ddq
    private
    ! Spectral discretization
    ! -------------------------------------
    real(wp), allocatable :: nus(:), weights(:)  ! (nnu) Wavenumbers and weights for the DDQ scheme
    ! -------------------------------------
    ! Functional approximations to cross-sections (fax) for line absorption
    ! -------------------------------------
    ! Which gases?
    character(len=gas_name_len), &
             allocatable :: fax_species_names(:)
    !
    ! Temperature dependence of absorption coefficients
    !
    real(wp), allocatable :: fax_a(:,:,:), fax_b(:,:,:) ! (0:2, nspecies, nnu)
    real(wp), allocatable :: fax_T0(:) ! (nspecies)
    !
    ! Pressure dependence
    !
    real(wp), allocatable :: fax_c(:,:,:)    ! (0:2, nspecies, nnu)
    real(wp), allocatable :: fax_p0(:)       ! (nspecies)
    !
    ! Reference cross-section, self-broadening factor
    !
    real(wp), allocatable :: fax_sigma0(:,:) ! (     nspecies, nnu), reference absorption coefficient at p_0, T_0
    real(wp), allocatable :: fax_S(:)        ! (     nspecies), self-broadening coefficients
    real(wp), allocatable :: fax_vmr0(:)     ! (     nspecies), reference volume mixing ratios
    ! -------------------------------------
    ! cross-section fits (xsec)
    ! -------------------------------------
    character(len=gas_name_len), &
              allocatable :: xsec_species_names(:)
    real(wp), allocatable :: xsec_p(:,:,:) ! (0:3, nspecies, nnu), (const, T, T^2, p scaling of cross-section)
    ! -------------------------------------
    ! MT_CKD continuum absorption (mtckd)
    ! -------------------------------------
    character(len=gas_name_len), &
              allocatable :: mtckd_species_names(:)
    real(wp), allocatable :: mtckd_cself(:, :), mtckd_cfrgn(:, :), mtckd_n(:, :) ! self- and foreign continuua
    real(wp)              :: mtckd_T0, mtckd_p0
    ! -------------------------------------
    ! Rayleigh
    ! -------------------------------------
    real(wp), allocatable :: rayleigh_xsec(:) ! (nnu)
    ! -------------------------------------
    ! Solar source
    ! -------------------------------------
    real(wp), allocatable :: solar_source(:) ! (nnu)
    ! -------------------------------------
    ! Unique list of all gases known to the class
    character(len=gas_name_len), &
              allocatable :: species_names(:)
  contains
    procedure, public :: source_is_internal
    procedure, public :: source_is_external
    procedure, public :: is_loaded
    procedure, public :: get_press_min
    procedure, public :: get_press_max
    procedure, public :: get_temp_min
    procedure, public :: get_temp_max
    procedure, public :: get_ngas
    procedure, public :: get_gases

    procedure, public :: load
    procedure, public :: set_tsi
    procedure, public :: finalize

    procedure, public  :: gas_optics_int
    procedure, public  :: gas_optics_ext

    ! Gas optics and the RTE solvers in one kernel
    procedure, public  :: rte_lw => rte_lw_ddq
    procedure, public  :: rte_sw => rte_sw_ddq
  end type
  !
  ! The nu-independent quantities of tau_absorption_from_fits, for rte_lw_ddq and rte_sw_ddq
  !
  type :: ty_nu_independent
    integer :: fax_ngas = 0, xsec_ngas = 0, mtckd_ngas = 0
    integer,  dimension(:),     allocatable :: fax_num_index, xsec_num_index, mtckd_num_index
    real(wp), dimension(:,:),   allocatable :: dry_num, fax_rs, inv_2kT
    real(wp), dimension(:,:,:), allocatable :: fax_x, fax_w, xsec_w, mt_logr, mt_sw, mt_fw
    real(wp), dimension(:,:,:), allocatable :: no_lut            ! stands in for cloud tables without clouds
  end type ty_nu_independent
contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute gas optical depth given temperature, pressure, and composition
  !
  function gas_optics_ext(this,                         &
                          play, plev, tlay, gas_desc,   & ! mandatory inputs
                          optical_props, toa_src,       & ! mandatory outputs
                          col_dry) result(error_msg)
    class(ty_gas_optics_ddq), intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props !!
    real(wp), dimension(:,:), intent(  out) :: toa_src   !! Incoming solar irradiance(ncol,ngpt)
    character(len=128)                      :: error_msg !! Empty if successful
    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry !! Column dry amount (molecules/cm^2); dim(ncol,nlay)
    ! --------------
    integer :: ncol, nlay, nnu
    integer :: icol, inu
    real(wp), dimension(size(play,1), size(play,2), this%get_ngpt()) :: tau_rayleigh
    real(wp), dimension(size(play,1), size(play,2)) :: col_dry_used

    ! ----------------------------------------------------------
    error_msg = ""
    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = this%get_ngpt()

    call optical_props%set_top_at_1(play(1,1) < play(1, nlay))

    !$acc        data copyout( optical_props%tau) create(   col_dry_used)
    !$omp target data map(from:optical_props%tau) map(alloc:col_dry_used)
    ! Absoption optical depth
    error_msg =  compute_tau_absorption(this,   &
                    play, plev, tlay, gas_desc, &
                    optical_props%tau, col_dry, &
                    col_dry_out = col_dry_used)

    if (error_msg == "") then
      if (.not. allocated(this%rayleigh_xsec)) then
        select type(optical_props)
          type is (ty_optical_props_2str)
            call zero_array(ncol, nlay, nnu, optical_props%ssa)
          type is (ty_optical_props_nstr)
            call zero_array(ncol, nlay, nnu, optical_props%ssa)
        end select
      else
        select type(optical_props)
          type is (ty_optical_props_2str)
            call add_tau_rayleigh(ncol, nlay, nnu,  &
                                  col_dry_used,          &
                                  this%rayleigh_xsec,    &
                                  optical_props%tau, optical_props%ssa)
          type is (ty_optical_props_nstr)
            call add_tau_rayleigh(ncol, nlay, nnu,  &
                                  col_dry_used,          &
                                  this%rayleigh_xsec,    &
                                  optical_props%tau, optical_props%ssa)
        end select
      end if
    end if
    !$acc end data
    !$omp end target data
    if (error_msg /= "") return

    select type(optical_props)
      type is (ty_optical_props_2str)
        call zero_array(ncol, nlay, nnu, optical_props%g)
      type is (ty_optical_props_nstr)
        call zero_array(optical_props%get_nmom(), &
                      ncol, nlay, nnu, optical_props%p)
    end select

    ! ----------------------------------------------------------
    !
    ! External source function is constant
    !
    if(check_extents) then
      if(.not. extents_are(toa_src, ncol, nnu)) &
        error_msg = "gas_optics(): array toa_src has wrong size"
    end if
    if(error_msg  /= '') return

    !$acc enter data create(toa_src)
    !$omp target enter data map(alloc:toa_src)
    !$acc parallel loop collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do inu = 1,nnu
       do icol = 1,ncol
          toa_src(icol,inu) = this%solar_source(inu) * this%weights(inu)
       end do
    end do
    !$acc exit data copyout(toa_src)
    !$omp target exit data map(from:toa_src)

  end function gas_optics_ext
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute gas optical depth and Planck source functions,
  !  given temperature, pressure, and composition
  !
  function gas_optics_int(this,                             &
                          play, plev, tlay, tsfc, gas_desc, &
                          optical_props, sources,           &
                          col_dry, tlev) result(error_msg)
    class(ty_gas_optics_ddq), intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    real(wp), dimension(:),   intent(in   ) :: tsfc      !! surface skin temperatures [K]; (ncol)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props !! Optical properties
    class(ty_source_func_lw    ),  &
                              intent(inout) :: sources    !! Planck sources
    character(len=128)                      :: error_msg  !! Empty if successful
    real(wp), dimension(:,:), intent(in   ), &         !!
                          optional, target :: col_dry, &  !! Column dry amount (molecules/cm^2); dim(ncol,nlay)
                                                 tlev     !! level temperatures [K]l (ncol,nlay+1)
    ! --------------
    integer :: ncol, nlay, nnu
    !
    ! Source function needs temperature at interfaces/levels and at layer centers
    !   Allocate small local array for tlev unconditionally
    !
    real(wp), dimension(size(plev,1),size(plev,2)), target  :: tlev_arr
    real(wp), dimension(:,:),                       pointer :: tlev_wk

    ! ----------------------------------------------------------
    error_msg = ""
    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = this%get_ngpt()
    ! --------------
    call optical_props%set_top_at_1(play(1,1) < play(1, nlay))

    ! Absoption optical depth
    error_msg = compute_tau_absorption(this,   &
                    play, plev, tlay, gas_desc, &
                    optical_props%tau, col_dry)
    if(error_msg /= "") return

    select type(optical_props)
      type is (ty_optical_props_2str)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(ncol, nlay, nnu, optical_props%g)
      type is (ty_optical_props_nstr)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(optical_props%get_nmom(), &
                      ncol, nlay, nnu, optical_props%p)
    end select

    !
    ! Planck function sources, interpolating tlev from tlay if necessary
    !
    if (present(tlev)) then
      tlev_wk => tlev
    else
      tlev_arr = interp_tlev_from_tlay(ncol, nlay, tlay, play, plev)
      tlev_wk => tlev_arr
    end if
    call compute_Planck_source(ncol, nlay, nnu,              &
                               this%nus, this%weights, tlay, &
                               sources%lay_source)
    call compute_Planck_source(ncol, nlay+1, nnu,            &
                               this%nus, this%weights, tlev_wk, &
                               sources%lev_source)
    call compute_Planck_source(ncol,         nnu,            &
                               this%nus, this%weights, tsfc, &
                               sources%sfc_source)

    call zero_array(ncol, nnu, sources%sfc_source_Jac)

  end function gas_optics_int
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Absorption optical depth
  !
  function compute_tau_absorption(this,             &
                          play, plev, tlay, gas_desc, &
                          tau_abs, col_dry, col_dry_out) result(error_msg)
    class(ty_gas_optics_ddq),   intent(in ) :: this
    real(wp), dimension(:,:),   intent(in ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),         intent(in ) :: gas_desc  !! Gas volume mixing ratios
    real(wp), dimension(:,:,:), intent(out) :: tau_abs   !! Cabsorption optical depth; (ncol,nlay,nnu)
    real(wp), dimension(:,:),   intent(in   ), &
                          optional, target :: col_dry     !! Column dry amount (molecules/cm^2); dim(ncol,nlay)
    real(wp), dimension(:,:),   intent(out   ), &
                        optional :: col_dry_out     !! Column dry amount (molecules/cm^2); dim(ncol,nlay)

    character(len=128)                      :: error_msg
    ! -----------------------
    real(wp), dimension(size(play, 1),size(play, 2)), &
                                    target  :: dry_num_arr
    real(wp), dimension(:,:),       pointer :: dry_num
    real(wp), dimension(:,:,:), allocatable :: vmrs
    integer,  dimension(:),     allocatable :: temp
    integer,  dimension(size(this%fax_species_names))   ::   fax_num_index
    integer,  dimension(size(this%xsec_species_names))  ::  xsec_num_index
    integer,  dimension(size(this%mtckd_species_names)) :: mtckd_num_index
    character(len=gas_name_len), &
              dimension(:), allocatable     :: gases_to_use
    integer :: ncol, nlay, ngas, nnu
    integer :: icol, ilay, igas, idx_h2o
    ! -----------------------------------------------------------
    error_msg = ""
    ncol = size(play, 1)
    nlay = size(play, 2)
    nnu = this%get_ngpt()
    ! -----------------------------------------------------------
    ! Data validation
    !
    ! initialization
    if (.not. this%is_loaded()) then
      error_msg = 'Spectral configuration not loaded'
      return
    end if
    !
    ! Check input data sizes and values
    !
    if(check_extents) then
      if(.not. extents_are(play, ncol, nlay  )) &
        error_msg = "gas_optics(): array play has wrong size"
      if(.not. extents_are(tlay, ncol, nlay  )) &
        error_msg = "gas_optics(): array tlay has wrong size"
      if(.not. extents_are(plev, ncol, nlay+1)) &
        error_msg = "gas_optics(): array plev has wrong size"
      if(.not. extents_are(tau_abs, ncol, nlay, nnu)) &
        error_msg = "gas_optics(): optical depth have the wrong extents"
      if(present(col_dry)) then
        if(.not. extents_are(col_dry, ncol, nlay)) &
          error_msg = "gas_optics(): array col_dry has wrong size"
      end if
    end if

    if(error_msg == '') then
      if(check_values) then
        if(any_vals_outside(play, this%get_press_min(),this%get_press_max())) &
          error_msg = "gas_optics(): array play has values outside range"
        if(any_vals_less_than(plev, 0._wp)) &
          error_msg = "gas_optics(): array plev has values outside range"
        if(any_vals_outside(tlay, this%get_temp_min(),  this%get_temp_max())) &
          error_msg = "gas_optics(): array tlay has values outside range"
        if(present(col_dry)) then
          if(any_vals_less_than(col_dry, 0._wp)) &
            error_msg = "gas_optics(): array col_dry has values outside range"
        end if
      end if
    end if

    if(error_msg /= '') return
    ! -----------------------------------------------------------

    call match_gases(this, gas_desc, gases_to_use, fax_num_index, xsec_num_index, mtckd_num_index)

    ! vmr array, gas index 0 stands in for gases users haven't supplied (vmr = 0);
    ! column dimension first for unit-stride access in the kernel
    allocate(vmrs(ncol, nlay, 0:size(gases_to_use)))
    !$acc        data copyin(play, plev, tlay) create(   vmrs)
    !$omp target data map(to:play, plev, tlay) map(alloc:vmrs)
    call zero_array(ncol, nlay, vmrs(:,:,0))
    do igas = 1, size(gases_to_use)
      error_msg = gas_desc%get_vmr(gases_to_use(igas), vmrs(:,:,igas))
      if (error_msg /= "") exit
    end do

    if (error_msg == "") then
      if (present(col_dry)) then
        dry_num => col_dry
      else
        dry_num => dry_num_arr
        idx_h2o = max(string_loc_in_array("h2o", gases_to_use), 0)
        dry_num_arr = get_layer_number(ncol, nlay,       &
                                      vmrs(:,:,idx_h2o), &
                                      plev)
                                       ! dry air column amounts computation
      end if
      call tau_absorption_from_fits(ncol, nlay, this%get_ngpt(), size(gases_to_use), &
                    this%nus, &
                    play, tlay, dry_num, vmrs, &
                    size(fax_num_index), fax_num_index,   &
                    this%fax_a, this%fax_b, this%fax_T0, this%fax_c, this%fax_p0, this%fax_sigma0, this%fax_S, this%fax_vmr0, &
                    size(xsec_num_index), xsec_num_index,  &
                    this%xsec_p,                           &
                    size(mtckd_num_index), mtckd_num_index, &
                    this%mtckd_cself, this%mtckd_cfrgn, this%mtckd_n, this%mtckd_T0, this%mtckd_p0, &
                    tau_abs)

      if (present(col_dry_out)) then
        !$acc                         parallel loop    collapse(2) copyin(dry_num)  copyout(col_dry_out)
        !$omp target teams distribute parallel do simd collapse(2) map(to:dry_num) map(from:col_dry_out)
        do ilay = 1, nlay
          do icol = 1, ncol
            col_dry_out(icol, ilay) = dry_num(icol, ilay)
          end do
        end do
      end if
    end if
    !$acc end data
    !$omp end target data
  end function compute_tau_absorption
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Gases provided by the user that the scheme knows, and where each set of fits finds them (0 if missing)
  !
  subroutine match_gases(this, gas_desc, gases_to_use, fax_num_index, xsec_num_index, mtckd_num_index)
    class(ty_gas_optics_ddq),    intent(in ) :: this
    type(ty_gas_concs),          intent(in ) :: gas_desc
    character(len=gas_name_len), dimension(:), allocatable, &
                                 intent(out) :: gases_to_use
    integer, dimension(:),       intent(out) :: fax_num_index, xsec_num_index, mtckd_num_index
    ! -----------------------
    character(len=32),           dimension(:), allocatable :: temp_gas_names
    character(len=gas_name_len), dimension(:), allocatable :: provided_gases
    integer :: igas

    allocate(temp_gas_names(gas_desc%get_num_gases()), &
             provided_gases(gas_desc%get_num_gases()))
    ! temp_gas_names is of uniform string length
    temp_gas_names(:) = gas_desc%get_gas_names()
    do igas = 1, size(temp_gas_names)
      provided_gases(igas) = trim(temp_gas_names(igas))
    end do
    allocate(gases_to_use(count([(string_in_array(provided_gases(igas), this%species_names), &
                                  igas = 1, size(provided_gases))])))
    gases_to_use = pack(provided_gases, &
                        mask = [(string_in_array(provided_gases(igas), this%species_names), &
                                 igas = 1, size(provided_gases))])
    fax_num_index   = [(max(string_loc_in_array(this%fax_species_names(igas), gases_to_use), 0), &
                        igas = 1, size(this%fax_species_names))]
    xsec_num_index  = [(max(string_loc_in_array(this%xsec_species_names(igas), gases_to_use), 0), &
                        igas = 1, size(this%xsec_species_names))]
    mtckd_num_index = [(max(string_loc_in_array(this%mtckd_species_names(igas), gases_to_use), 0), &
                        igas = 1, size(this%mtckd_species_names))]
  end subroutine match_gases
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Gas optics and rte_lw in one kernel: broadband longwave fluxes without scattering and with one
  !   quadrature angle, as gas_optics followed by rte_lw, without fields of spectral dimension.
  !   With clouds (condensate paths, effective radii and the lookup tables of the cloud optics,
  !   (nsteps, nnu, 3) for extinction, single-scattering albedo and asymmetry), flux_up and flux_dn
  !   are all-sky fluxes and the optional *_clrsky arguments clear-sky ones.
  !
  function rte_lw_ddq(this, play, plev, tlay, tsfc, gas_desc, sfc_emis, flux_up, flux_dn, tlev, &
                      flux_up_clrsky, flux_dn_clrsky,                                           &
                      clwp, ciwp, reliq, reice, lut_liq, radliq_lwr, liq_step_size,             &
                      lut_ice, radice_lwr, ice_step_size) result(error_msg)
    class(ty_gas_optics_ddq), intent(in ) :: this
    real(wp), dimension(:,:), intent(in ) :: play, &     !! layer pressures [Pa, mb]; (ncol,nlay)
                                             plev, &     !! level pressures [Pa, mb]; (ncol,nlay+1)
                                             tlay        !! layer temperatures [K]; (ncol,nlay)
    real(wp), dimension(:),   intent(in ) :: tsfc        !! surface skin temperatures [K]; (ncol)
    type(ty_gas_concs),       intent(in ) :: gas_desc    !! Gas volume mixing ratios
    real(wp), dimension(:,:), intent(in ) :: sfc_emis    !! emissivity at surface []; (nnu, ncol)
    real(wp), dimension(:,:), target, intent(inout) :: flux_up, flux_dn  !! broadband fluxes [W/m2]; (ncol,nlay+1)
    real(wp), dimension(:,:), intent(in ) :: tlev        !! level temperatures [K]; (ncol,nlay+1)
    real(wp), dimension(:,:), target, optional, intent(inout) :: flux_up_clrsky, flux_dn_clrsky
    real(wp), dimension(:,:),   optional, intent(in ) :: clwp, ciwp, reliq, reice !! [g/m2], [microns]; (ncol,nlay)
    real(wp), dimension(:,:,:), optional, intent(in ) :: lut_liq, lut_ice
    real(wp),                   optional, intent(in ) :: radliq_lwr, liq_step_size, radice_lwr, ice_step_size
    character(len=128)                    :: error_msg
    ! -----------------------
    ! Secant and weight of the one-angle Gaussian quadrature in rte_lw
    real(wp), parameter :: lw_Ds = 1._wp/0.6096748751_wp, lw_wts = 1._wp
    type(ty_nu_independent) :: st
    real(wp), dimension(:,:,:), allocatable :: trans, source_up
    real(wp), dimension(:,:),   pointer     :: flux_up_clr, flux_dn_clr
    integer     :: ncol, nlay, nnu
    logical(wl) :: top_at_1, do_clrsky, do_allsky
    ! -----------------------
    ncol = size(play, 1)
    nlay = size(play, 2)
    nnu  = this%get_ngpt()
    error_msg = check_fused_inputs(this, ncol, nlay, plev, tlay, flux_up, flux_dn, &
                                   clwp, ciwp, reliq, reice, lut_liq, lut_ice,     &
                                   radliq_lwr, liq_step_size, radice_lwr, ice_step_size)
    if (check_extents .and. error_msg == '') then
      if(.not. extents_are(tlev, ncol, nlay+1)) error_msg = "rte_lw_ddq(): array tlev has wrong size"
      if(.not. extents_are(sfc_emis, nnu, ncol)) error_msg = "rte_lw_ddq(): array sfc_emis has wrong size"
    end if
    if (error_msg /= '') return

    do_allsky = present(clwp)
    do_clrsky = .not. do_allsky .or. present(flux_up_clrsky)
    if (do_allsky .and. do_clrsky) then
      flux_up_clr => flux_up_clrsky
      flux_dn_clr => flux_dn_clrsky
    else
      flux_up_clr => flux_up
      flux_dn_clr => flux_dn
    end if
    top_at_1 = play(1, 1) < play(1, nlay)
    allocate(trans(ncol, nlay, 2), source_up(ncol, nlay, 2))
    !$acc        data copyin(play, plev, tlay, tlev, tsfc, sfc_emis) copy(flux_up, flux_dn) create(trans, source_up)
    !$omp target data map(to:play, plev, tlay, tlev, tsfc, sfc_emis) map(tofrom:flux_up, flux_dn) map(alloc:trans, source_up)
    !$acc        data copy(flux_up_clr, flux_dn_clr) if(do_allsky .and. do_clrsky)
    !$omp target data map(tofrom:flux_up_clr, flux_dn_clr) if(do_allsky .and. do_clrsky)
    error_msg = tau_absorption_invariants(this, play, plev, tlay, gas_desc, st)
    if (error_msg == "") then
      if (do_allsky) then
        !$acc        data copyin(clwp, ciwp, reliq, reice)
        !$omp target data map(to:clwp, ciwp, reliq, reice)
        call lw_solver_noscat_ddq(ncol, nlay, nnu, logical(top_at_1, wl), lw_Ds, lw_wts, this%nus, this%weights, &
                                  play, tlay, tlev, tsfc, sfc_emis,                                         &
                                  st%fax_ngas, this%fax_a, this%fax_b, this%fax_T0, this%fax_c, this%fax_sigma0, &
                                  st%xsec_ngas, this%xsec_p, st%mtckd_ngas, this%mtckd_cself, this%mtckd_cfrgn, &
                                  this%mtckd_n, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w,                    &
                                  st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT,                                &
                                  do_clrsky, do_allsky, clwp, ciwp, reliq, reice,                            &
                                  size(lut_liq, 1), liq_step_size, radliq_lwr, lut_liq,                      &
                                  size(lut_ice, 1), ice_step_size, radice_lwr, lut_ice,                      &
                                  trans, source_up, flux_up_clr, flux_dn_clr, flux_up, flux_dn)
        !$acc end data
        !$omp end target data
      else
        ! Without clouds the cloud arguments are not referenced
        call lw_solver_noscat_ddq(ncol, nlay, nnu, logical(top_at_1, wl), lw_Ds, lw_wts, this%nus, this%weights, &
                                  play, tlay, tlev, tsfc, sfc_emis,                                         &
                                  st%fax_ngas, this%fax_a, this%fax_b, this%fax_T0, this%fax_c, this%fax_sigma0, &
                                  st%xsec_ngas, this%xsec_p, st%mtckd_ngas, this%mtckd_cself, this%mtckd_cfrgn, &
                                  this%mtckd_n, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w,                    &
                                  st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT,                                &
                                  do_clrsky, do_allsky, play, play, play, play,                              &
                                  2, 1._wp, 0._wp, st%no_lut, 2, 1._wp, 0._wp, st%no_lut,                    &
                                  trans, source_up, flux_up_clr, flux_dn_clr, flux_up, flux_dn)
      end if
    end if
    call finalize_invariants(st)
    !$acc end data
    !$omp end target data
    !$acc end data
    !$omp end target data
  end function rte_lw_ddq
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Gas optics and rte_sw in one kernel: broadband shortwave fluxes with two-stream scattering, as
  !   gas_optics followed by rte_sw; flux_dn is total. The solar source of the data is scaled to tsi
  !   if given, and per column by toa_factor (for example 0 at night). Clouds as for rte_lw_ddq.
  !   Optionally, surface fluxes weighted per spectral point by sfc_frac(nnu, nfrac), for the
  !   all-sky fluxes, or the clear-sky ones without clouds.
  !
  function rte_sw_ddq(this, play, plev, tlay, gas_desc, mu0, sfc_alb_dir, sfc_alb_dif, flux_up, flux_dn, &
                      flux_up_clrsky, flux_dn_clrsky, tsi, toa_factor,                                    &
                      clwp, ciwp, reliq, reice, lut_liq, radliq_lwr, liq_step_size,                       &
                      lut_ice, radice_lwr, ice_step_size,                                                 &
                      sfc_frac, sfc_dn_dir, sfc_dn_dif, sfc_up, flux_dn_dir) result(error_msg)
    class(ty_gas_optics_ddq), intent(in ) :: this
    real(wp), dimension(:,:), intent(in ) :: play, &     !! layer pressures [Pa, mb]; (ncol,nlay)
                                             plev, &     !! level pressures [Pa, mb]; (ncol,nlay+1)
                                             tlay        !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in ) :: gas_desc    !! Gas volume mixing ratios
    real(wp), dimension(:),   intent(in ) :: mu0         !! cosine of solar zenith angle, > 0; (ncol)
    real(wp), dimension(:,:), intent(in ) :: sfc_alb_dir, sfc_alb_dif !! surface albedo; (nnu,ncol)
    real(wp), dimension(:,:), target, intent(inout) :: flux_up, flux_dn  !! broadband fluxes [W/m2]; (ncol,nlay+1)
    real(wp), dimension(:,:), target, optional, intent(inout) :: flux_up_clrsky, flux_dn_clrsky
    real(wp),                   optional, intent(in ) :: tsi
    real(wp), dimension(:),     optional, intent(in ) :: toa_factor !! (ncol)
    real(wp), dimension(:,:),   optional, intent(in ) :: clwp, ciwp, reliq, reice
    real(wp), dimension(:,:,:), optional, intent(in ) :: lut_liq, lut_ice
    real(wp),                   optional, intent(in ) :: radliq_lwr, liq_step_size, radice_lwr, ice_step_size
    real(wp), dimension(:,:),   optional, intent(in ) :: sfc_frac   !! (nnu,nfrac)
    real(wp), dimension(:,:),   optional, intent(inout) :: sfc_dn_dir, sfc_dn_dif, sfc_up !! (ncol,nfrac)
    real(wp), dimension(:,:), target, optional, intent(inout) :: flux_dn_dir !! direct beam [W/m2]; (ncol,nlay+1)
    character(len=128)                    :: error_msg
    ! -----------------------
    type(ty_nu_independent) :: st
    real(wp), dimension(:,:), pointer     :: flux_dir
    real(wp), dimension(:,:), allocatable :: tau, ssa, Rdif, Tdif, source_dn, source_up, albedo, src
    real(wp), dimension(:,:), pointer     :: flux_up_clr, flux_dn_clr
    real(wp), dimension(size(this%nus))    :: solar_source
    real(wp), dimension(size(play, 1))     :: factor
    real(wp), dimension(size(this%nus), 1) :: no_frac
    real(wp), dimension(:), allocatable    :: rayleigh_xsec
    integer     :: ncol, nlay, nnu, nfrac, icol
    logical(wl) :: top_at_1, do_clrsky, do_allsky, has_rayleigh
    ! -----------------------
    ncol = size(play, 1)
    nlay = size(play, 2)
    nnu  = this%get_ngpt()
    error_msg = check_fused_inputs(this, ncol, nlay, plev, tlay, flux_up, flux_dn, &
                                   clwp, ciwp, reliq, reice, lut_liq, lut_ice,     &
                                   radliq_lwr, liq_step_size, radice_lwr, ice_step_size)
    if (error_msg == '' .and. .not. allocated(this%solar_source)) error_msg = "rte_sw_ddq(): no solar source"
    if (check_extents .and. error_msg == '') then
      if(.not. extents_are(mu0, ncol)) error_msg = "rte_sw_ddq(): array mu0 has wrong size"
      if(.not. extents_are(sfc_alb_dir, nnu, ncol) .or. .not. extents_are(sfc_alb_dif, nnu, ncol)) &
        error_msg = "rte_sw_ddq(): surface albedo arrays have wrong size"
    end if
    if (error_msg /= '') return

    do_allsky = present(clwp)
    do_clrsky = .not. do_allsky .or. present(flux_up_clrsky)
    if (do_allsky .and. do_clrsky) then
      flux_up_clr => flux_up_clrsky
      flux_dn_clr => flux_dn_clrsky
    else
      flux_up_clr => flux_up
      flux_dn_clr => flux_dn
    end if
    if (present(flux_dn_dir)) then
      flux_dir => flux_dn_dir
    else
      flux_dir => flux_up
    end if
    solar_source(:) = this%solar_source(:)
    if (present(tsi)) solar_source(:) = solar_source(:) * tsi / sum(this%solar_source(:) * this%weights(:))
    has_rayleigh = allocated(this%rayleigh_xsec)
    if (has_rayleigh) then
      rayleigh_xsec = this%rayleigh_xsec
    else
      allocate(rayleigh_xsec(nnu))
      rayleigh_xsec(:) = 0._wp
    end if
    nfrac = 0
    if (present(sfc_frac)) nfrac = size(sfc_frac, 2)
    no_frac(:,:) = 0._wp
    top_at_1 = play(1, 1) < play(1, nlay)
    allocate(tau(ncol, nlay), ssa(ncol, nlay), Rdif(ncol, nlay), Tdif(ncol, nlay), &
             source_dn(ncol, nlay), source_up(ncol, nlay), albedo(ncol, nlay+1), src(ncol, nlay+1))
    !$acc        data copyin(play, plev, tlay, mu0, sfc_alb_dir, sfc_alb_dif) copy(flux_up, flux_dn) &
    !$acc             create(factor, tau, ssa, Rdif, Tdif, source_dn, source_up, albedo, src)
    !$omp target data map(to:play, plev, tlay, mu0, sfc_alb_dir, sfc_alb_dif) map(tofrom:flux_up, flux_dn) &
    !$omp             map(alloc:factor, tau, ssa, Rdif, Tdif, source_dn, source_up, albedo, src)
    !$acc        data copy(flux_up_clr, flux_dn_clr) if(do_allsky .and. do_clrsky)
    !$omp target data map(tofrom:flux_up_clr, flux_dn_clr) if(do_allsky .and. do_clrsky)
    !$acc        data copy(flux_dir) if(present(flux_dn_dir))
    !$omp target data map(tofrom:flux_dir) if(present(flux_dn_dir))
    ! toa_factor may exist on the device only
    if (present(toa_factor)) then
      !$acc                         parallel loop    gang vector copyin(toa_factor) default(present)
      !$omp target teams distribute parallel do simd map(to:toa_factor)
      do icol = 1, ncol
        factor(icol) = toa_factor(icol)
      end do
    else
      !$acc                         parallel loop    gang vector default(present)
      !$omp target teams distribute parallel do simd
      do icol = 1, ncol
        factor(icol) = 1._wp
      end do
    end if
    error_msg = tau_absorption_invariants(this, play, plev, tlay, gas_desc, st)
    if (error_msg == "") then
      if (do_allsky .and. nfrac > 0) then
        !$acc        data copyin(clwp, ciwp, reliq, reice) copy(sfc_dn_dir, sfc_dn_dif, sfc_up)
        !$omp target data map(to:clwp, ciwp, reliq, reice) map(tofrom:sfc_dn_dir, sfc_dn_dif, sfc_up)
        call sw_solver(clwp, ciwp, reliq, reice, size(lut_liq, 1), liq_step_size, radliq_lwr, lut_liq, &
                       size(lut_ice, 1), ice_step_size, radice_lwr, lut_ice, nfrac, sfc_frac, sfc_dn_dir, sfc_dn_dif, sfc_up)
        !$acc end data
        !$omp end target data
      else if (do_allsky) then
        !$acc        data copyin(clwp, ciwp, reliq, reice)
        !$omp target data map(to:clwp, ciwp, reliq, reice)
        call sw_solver(clwp, ciwp, reliq, reice, size(lut_liq, 1), liq_step_size, radliq_lwr, lut_liq, &
                       size(lut_ice, 1), ice_step_size, radice_lwr, lut_ice, 0, no_frac, flux_up, flux_up, flux_up)
        !$acc end data
        !$omp end target data
      else if (nfrac > 0) then
        !$acc        data copy(sfc_dn_dir, sfc_dn_dif, sfc_up)
        !$omp target data map(tofrom:sfc_dn_dir, sfc_dn_dif, sfc_up)
        call sw_solver(play, play, play, play, 2, 1._wp, 0._wp, st%no_lut, 2, 1._wp, 0._wp, st%no_lut, &
                       nfrac, sfc_frac, sfc_dn_dir, sfc_dn_dif, sfc_up)
        !$acc end data
        !$omp end target data
      else
        call sw_solver(play, play, play, play, 2, 1._wp, 0._wp, st%no_lut, 2, 1._wp, 0._wp, st%no_lut, &
                       0, no_frac, flux_up, flux_up, flux_up)
      end if
    end if
    call finalize_invariants(st)
    !$acc end data
    !$omp end target data
    !$acc end data
    !$omp end target data
    !$acc end data
    !$omp end target data
  contains
    ! The kernel with the arguments that are the same in all cases above
    subroutine sw_solver(c_lwp, c_iwp, c_rel, c_rei, n_liq, step_liq, lwr_liq, t_liq, &
                         n_ice, step_ice, lwr_ice, t_ice, n_frac, frac, f_dir, f_dif, f_up)
      real(wp), dimension(:,:),   intent(in)    :: c_lwp, c_iwp, c_rel, c_rei, frac
      integer,                    intent(in)    :: n_liq, n_ice, n_frac
      real(wp),                   intent(in)    :: step_liq, lwr_liq, step_ice, lwr_ice
      real(wp), dimension(:,:,:), intent(in)    :: t_liq, t_ice
      real(wp), dimension(:,:),   intent(inout) :: f_dir, f_dif, f_up
      call sw_solver_2stream_ddq(ncol, nlay, nnu, logical(top_at_1, wl), this%nus, this%weights, solar_source, &
                                 mu0, factor, sfc_alb_dir, sfc_alb_dif, play, tlay, st%dry_num,               &
                                 st%fax_ngas, this%fax_a, this%fax_b, this%fax_T0, this%fax_c, this%fax_sigma0, &
                                 st%xsec_ngas, this%xsec_p, st%mtckd_ngas, this%mtckd_cself, this%mtckd_cfrgn, &
                                 this%mtckd_n, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w,                       &
                                 st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT, has_rayleigh, rayleigh_xsec,      &
                                 do_clrsky, do_allsky, c_lwp, c_iwp, c_rel, c_rei,                             &
                                 n_liq, step_liq, lwr_liq, t_liq, n_ice, step_ice, lwr_ice, t_ice,             &
                                 n_frac, frac, f_dir, f_dif, f_up,                                             &
                                 tau, ssa, Rdif, Tdif, source_dn, source_up, albedo, src,                      &
                                 flux_up_clr, flux_dn_clr, flux_up, flux_dn,                                   &
                                 logical(present(flux_dn_dir), wl), flux_dir)
    end subroutine sw_solver
  end function rte_sw_ddq
  !--------------------------------------------------------------------------------------------------------------------
  function check_fused_inputs(this, ncol, nlay, plev, tlay, flux_up, flux_dn, &
                              clwp, ciwp, reliq, reice, lut_liq, lut_ice,     &
                              radliq_lwr, liq_step_size, radice_lwr, ice_step_size) result(error_msg)
    class(ty_gas_optics_ddq),   intent(in) :: this
    integer,                    intent(in) :: ncol, nlay
    real(wp), dimension(:,:),   intent(in) :: plev, tlay, flux_up, flux_dn
    real(wp), dimension(:,:),   optional, intent(in) :: clwp, ciwp, reliq, reice
    real(wp), dimension(:,:,:), optional, intent(in) :: lut_liq, lut_ice
    real(wp),                   optional, intent(in) :: radliq_lwr, liq_step_size, radice_lwr, ice_step_size
    character(len=128)                     :: error_msg

    error_msg = ""
    if (.not. this%is_loaded()) then
      error_msg = 'Spectral configuration not loaded'
      return
    end if
    if (present(clwp) .neqv. (present(ciwp) .and. present(reliq) .and. present(reice) .and.          &
                              present(lut_liq) .and. present(radliq_lwr) .and. present(liq_step_size) .and. &
                              present(lut_ice) .and. present(radice_lwr) .and. present(ice_step_size))) then
      error_msg = "fused fluxes: provide all cloud arguments or none"
      return
    end if
    if (check_extents) then
      if(.not. extents_are(plev, ncol, nlay+1) .or. .not. extents_are(tlay, ncol, nlay)) &
        error_msg = "fused fluxes: atmospheric arrays have wrong size"
      if(.not. extents_are(flux_up, ncol, nlay+1) .or. .not. extents_are(flux_dn, ncol, nlay+1)) &
        error_msg = "fused fluxes: flux arrays have wrong size"
      if (present(clwp)) then
        if(.not. extents_are(clwp, ncol, nlay) .or. .not. extents_are(ciwp, ncol, nlay) .or. &
           .not. extents_are(reliq, ncol, nlay) .or. .not. extents_are(reice, ncol, nlay))   &
          error_msg = "fused fluxes: cloud arrays have wrong size"
        if(size(lut_liq, 2) /= this%get_ngpt() .or. size(lut_ice, 2) /= this%get_ngpt() .or. &
           size(lut_liq, 3) /= 3 .or. size(lut_ice, 3) /= 3) &
          error_msg = "fused fluxes: cloud tables have wrong size"
      end if
    end if
  end function check_fused_inputs
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! The nu-independent quantities of tau_absorption_from_fits on the device, one gas at a time;
  !   play, plev and tlay must be present on the device
  !
  function tau_absorption_invariants(this, play, plev, tlay, gas_desc, st) result(error_msg)
    class(ty_gas_optics_ddq), intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, plev, tlay
    type(ty_gas_concs),       intent(in   ) :: gas_desc
    type(ty_nu_independent),  intent(inout) :: st
    character(len=128)                      :: error_msg
    character(len=gas_name_len), dimension(:), allocatable :: gases_to_use
    real(wp), dimension(size(play, 1), size(play, 2)) :: vmr
    integer :: ncol, nlay, nnu, igas, idx_h2o

    error_msg = ""
    ncol = size(play, 1)
    nlay = size(play, 2)
    nnu  = this%get_ngpt()
    allocate(st%fax_num_index(size(this%fax_species_names)), st%xsec_num_index(size(this%xsec_species_names)), &
             st%mtckd_num_index(size(this%mtckd_species_names)))
    call match_gases(this, gas_desc, gases_to_use, st%fax_num_index, st%xsec_num_index, st%mtckd_num_index)
    st%fax_ngas   = size(st%fax_num_index)
    st%xsec_ngas  = size(st%xsec_num_index)
    st%mtckd_ngas = size(st%mtckd_num_index)
    idx_h2o = max(string_loc_in_array("h2o", gases_to_use), 0)
    allocate(st%dry_num(ncol, nlay), st%fax_x(ncol, nlay, st%fax_ngas), st%fax_w(ncol, nlay, st%fax_ngas),  &
             st%fax_rs(st%fax_ngas, nnu), st%xsec_w(ncol, nlay, st%xsec_ngas),                               &
             st%mt_logr(ncol, nlay, st%mtckd_ngas), st%mt_sw(ncol, nlay, st%mtckd_ngas),                     &
             st%mt_fw(ncol, nlay, st%mtckd_ngas), st%inv_2kT(ncol, nlay), st%no_lut(2, nnu, 3))
    st%no_lut(:,:,:) = 0._wp
    !$acc enter data copyin(st)
    !$omp target enter data map(to:st)
    !$acc enter data create(st%dry_num, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w) &
    !$acc            create(st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT) copyin(st%no_lut)
    !$omp target enter data map(alloc:st%dry_num, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w) &
    !$omp            map(alloc:st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT) map(to:st%no_lut)
    !$acc        data create(vmr)
    !$omp target data map(alloc:vmr)
    ! Dry air needs water vapor first; gases users haven't supplied (index 0) have vmr = 0
    call zero_array(ncol, nlay, vmr)
    if (idx_h2o > 0) error_msg = gas_desc%get_vmr(gases_to_use(idx_h2o), vmr)
    if (error_msg == "") call get_layer_number_ddq(ncol, nlay, vmr, plev, st%dry_num)
    if (error_msg == "") call zero_array(ncol, nlay, vmr)
    do igas = 0, size(gases_to_use)
      if (error_msg /= "") exit
      if (igas > 0) error_msg = gas_desc%get_vmr(gases_to_use(igas), vmr)
      if (error_msg == "") &
        call tau_absorption_from_fits_invariants(ncol, nlay, nnu, igas, play, tlay, st%dry_num, vmr, &
                  st%fax_ngas, st%fax_num_index, this%fax_c, this%fax_p0, this%fax_S, this%fax_vmr0,  &
                  st%xsec_ngas, st%xsec_num_index,                                                    &
                  st%mtckd_ngas, st%mtckd_num_index, this%mtckd_T0, this%mtckd_p0,                    &
                  st%fax_x, st%fax_w, st%fax_rs, st%xsec_w, st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT)
    end do
    !$acc end data
    !$omp end target data
  end function tau_absorption_invariants
  !--------------------------------------------------------------------------------------------------------------------
  subroutine finalize_invariants(st)
    type(ty_nu_independent), intent(inout) :: st
    if (.not. allocated(st%dry_num)) return
    !$acc exit data delete(st%dry_num, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w) &
    !$acc           delete(st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT, st%no_lut)
    !$omp target exit data map(release:st%dry_num, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w) &
    !$omp           map(release:st%mt_logr, st%mt_sw, st%mt_fw, st%inv_2kT, st%no_lut)
    !$acc exit data delete(st)
    !$omp target exit data map(release:st)
    deallocate(st%dry_num, st%fax_x, st%fax_w, st%fax_rs, st%xsec_w, st%mt_logr, st%mt_sw, st%mt_fw, &
               st%inv_2kT, st%no_lut)
  end subroutine finalize_invariants
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Initialiation
  !
  function load(this,                       &
                nus, weights, &
                fax_species_names, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, fax_vmr0, &
                xsec_species_names, xsec_p, &
                mtckd_species_names, mtckd_cself, mtckd_cfrgn, mtckd_n, mtckd_T0, mtckd_p0, &
                rayleigh_xsec, solar_source) &
                result(error_msg)
    class(ty_gas_optics_ddq), intent(inout) :: this
    ! -------------------------------------
    ! Spectral discretization
    real(wp), intent(in) :: nus(:), weights(:)  ! (nnu) Wavenumbers and weights for the DDQ scheme
    ! -------------------------------------
    ! Functional approximations to cross-sections (fax) for line absorption
    character(len=*), intent(in) :: fax_species_names(:)
    real(wp), intent(in) :: fax_a(:,:,:), fax_b(:,:,:) ! (0:2, nspecies, nnu)
    real(wp), intent(in) :: fax_T0(:)       ! (nspecies)
    real(wp), intent(in) :: fax_c(:,:,:)    ! (0:2, nspecies, nnu)
    real(wp), intent(in) :: fax_p0(:)       ! (nspecies)
    real(wp), intent(in) :: fax_sigma0(:,:) ! (     nspecies, nnu), reference absorption coefficient at p_0, T_0
    real(wp), intent(in) :: fax_S(:)        ! (     nspecies), self-broadening coefficients
    real(wp), intent(in) :: fax_vmr0(:)     ! (     nspecies), reference volume mixing ratios
    ! -------------------------------------
    ! cross-section fits (xsec)
    character(len=*), intent(in) :: xsec_species_names(:)
    real(wp),         intent(in) :: xsec_p(:,:,:) ! (0:3, nspecies, nnu), (const, T, T^2, p scaling of cross-section)
    ! -------------------------------------
    ! MT_CKD continuum absorption (mtckd)
    character(len=*), intent(in) :: mtckd_species_names(:)
    real(wp),         intent(in) :: mtckd_cself(:, :), mtckd_cfrgn(:, :), mtckd_n(:, :) ! self- and foreign continuua
    real(wp),         intent(in) :: mtckd_T0, mtckd_p0
    ! -------------------------------------
    ! Rayleigh cross-sections
    real(wp), optional, &
                      intent(in) :: rayleigh_xsec(:) ! (nnu)
    ! solar source
    real(wp), optional, &
                      intent(in) :: solar_source(:) ! (nnu)
    ! -------------------------------------
    character(len=gas_name_len), dimension(:), &
      allocatable :: all_names
    character(len = 128) :: error_msg

    ! -------------------------------------
    integer :: nnu, fax_nspecies, xsec_nspecies, mtckd_nspecies, i
    ! -------------------------------------
    !$acc enter data copyin(this)
    call this%finalize()
    error_msg = this%ty_optical_props%init(band_lims_wvn = &
                                             spread(nus, dim=1, ncopies=2))
    if(len_trim(error_msg) /= 0) return

    nnu = size(nus)
    fax_nspecies   = size(  fax_species_names)
    xsec_nspecies  = size( xsec_species_names)
    mtckd_nspecies = size(mtckd_species_names)

    !
    ! Check dimensions of input data
    !
    if (.not. extents_are(weights, nnu)) &
      error_msg = "Need the same number of weights, nus"
    if (.not. extents_are(fax_a, fax_norder+1, fax_nspecies, nnu) .or. &
        .not. extents_are(fax_b, fax_norder+1, fax_nspecies, nnu)) &
      error_msg = "Wrong dimensions for fax_a and/or fax_b"
    if (.not. extents_are(fax_c, fax_nterms+1, fax_nspecies, nnu)) &
      error_msg = "Wrong dimensions for fax_c"
    if (.not. extents_are(fax_sigma0,          fax_nspecies, nnu)) &
      error_msg = "Wrong dimensions for fax_sigma0"
    if (.not. extents_are(fax_T0, fax_nspecies) .or. &
        .not. extents_are(fax_p0, fax_nspecies) .or. &
        .not. extents_are(fax_vmr0, fax_nspecies) .or. &
        .not. extents_are(fax_S,  fax_nspecies)) &
      error_msg = "fax_T0, fax_p0, fax_vmr0, fax_S depend only on species"

    if(xsec_nspecies > 0) then
      if (.not. extents_are(xsec_p, xsec_nterms+1, xsec_nspecies, nnu)) &
        error_msg = "Wrong dimensions for xsec_p"
    end if

    if(mtckd_nspecies > 0) then
      if (.not. extents_are(mtckd_cself, mtckd_nspecies, nnu) .or.  &
          .not. extents_are(mtckd_cfrgn, mtckd_nspecies, nnu) .or. &
          .not. extents_are(mtckd_n,     mtckd_nspecies, nnu)) &
        error_msg = "Wrong dimensions for mtckd arrays"
    end if

    if(present(rayleigh_xsec)) then
      if (.not. extents_are(rayleigh_xsec, nnu)) &
        error_msg = "rayleigh_xsec should have same extents as nus"
    end if
    if(present(solar_source)) then
      if (.not. extents_are(solar_source, nnu)) &
        error_msg = "solar_source should have same extents as nus"
    end if
    if(error_msg /= "") return
    ! ----------------------------------------------------------
    allocate(this%nus(nnu), this%weights(nnu))
    this%nus = nus
    this%weights = weights

    allocate(this%fax_species_names(  fax_nspecies),      &
             this%fax_a(0:fax_norder, fax_nspecies, nnu), &
             this%fax_b(0:fax_norder, fax_nspecies, nnu), &
             this%fax_c(0:fax_nterms, fax_nspecies, nnu), &
             this%fax_sigma0(         fax_nspecies, nnu), &
             this%fax_S (             fax_nspecies),      &
             this%fax_vmr0(           fax_nspecies))
    this%fax_species_names = fax_species_names
    this%fax_a = fax_a
    this%fax_b = fax_b
    this%fax_c = fax_c
    this%fax_sigma0 = fax_sigma0
    this%fax_S  = fax_S
    this%fax_vmr0 = fax_vmr0
    this%fax_T0 = fax_T0
    this%fax_p0 = fax_p0

    allocate(this%xsec_species_names(   xsec_nspecies),    &
             this%xsec_p(0:xsec_nterms, xsec_nspecies, nnu))
    this%xsec_species_names = xsec_species_names
    this%xsec_p = xsec_p

    allocate(this%mtckd_species_names(mtckd_nspecies),    &
             this%mtckd_cself(mtckd_nspecies, nnu), &
             this%mtckd_cfrgn(mtckd_nspecies, nnu), &
             this%mtckd_n    (mtckd_nspecies, nnu))
    this%mtckd_species_names = mtckd_species_names
    this%mtckd_cself = mtckd_cself
    this%mtckd_cfrgn = mtckd_cfrgn
    this%mtckd_n     = mtckd_n
    this%mtckd_T0    = mtckd_T0
    this%mtckd_p0    = mtckd_p0

    !$acc        enter data copyin(this%nus, this%weights)                                   &
    !$acc                   copyin(this%fax_a, this%fax_b, this%fax_c, this%fax_sigma0)          &
    !$acc                   copyin(this%fax_S, this%fax_vmr0, this%fax_T0, this%fax_p0)          &
    !$acc                   copyin(this%xsec_p, this%mtckd_cself, this%mtckd_cfrgn, this%mtckd_n)
    !$omp target enter data map(to:this%nus, this%weights)                                       &
    !$omp                   map(to:this%fax_a, this%fax_b, this%fax_c, this%fax_sigma0)          &
    !$omp                   map(to:this%fax_S, this%fax_vmr0, this%fax_T0, this%fax_p0)          &
    !$omp                   map(to:this%xsec_p, this%mtckd_cself, this%mtckd_cfrgn, this%mtckd_n)

    if (present(rayleigh_xsec)) then
      allocate(this%rayleigh_xsec(nnu))
      this%rayleigh_xsec = rayleigh_xsec
      !$acc        enter data copyin(this%rayleigh_xsec)
      !$omp target enter data map(to:this%rayleigh_xsec)
    end if
    if (present(solar_source)) then
      allocate(this%solar_source(nnu))
      this%solar_source = solar_source
      !$acc        enter data copyin(this%solar_source)
      !$omp target enter data map(to:this%solar_source)
    end if

    ! Make a consolidated list of all gases
    allocate(all_names(fax_nspecies + xsec_nspecies + mtckd_nspecies))
    all_names(:fax_nspecies              ) = &
      fax_species_names(:)
    all_names(fax_nspecies+1:fax_nspecies+xsec_nspecies) = &
      xsec_species_names(:)
    all_names(fax_nspecies+xsec_nspecies+1:) = &
      mtckd_species_names(:)
    do i = 2, size(all_names)
      if (string_in_array(all_names(i), all_names(:i-1))) all_names(i) = ""
    end do

    allocate(this%species_names(count(all_names /= "")))
    this%species_names(:) = pack(all_names, mask = all_names /= "")
  end function load
  !------------------------------------------------------------------------------------------
  ! Total solar irradiance
  !------------------------------------------------------------------------------------------
  function set_tsi(this, tsi) result(error_msg)
    !
    !> Scale the solar source function without changing the spectral distribution
    !
    class(ty_gas_optics_ddq), intent(inout) :: this
    real(wp),                 intent(in   ) :: tsi !! user-specified total solar irradiance;
    character(len=128)                      :: error_msg !! Empty if successful

    real(wp) :: norm
    integer  :: igpt, length
    ! ----------------------------------------------------------
    error_msg = ""
    if(tsi < 0._wp) then
      error_msg = 'tsi out of range'
    else if(.not. allocated(this%solar_source)) then
      error_msg = 'no solar source to rescale with TSI'
    else
      !
      ! Scale the solar source function to the input tsi
      !
      norm = 0._wp
      length = size(this%solar_source)
      !$acc parallel loop gang vector reduction(+:norm)
      !$omp target teams distribute parallel do simd reduction(+:norm)
      do igpt = 1, length
         norm = norm + this%solar_source(igpt)
      end do

      norm = 1._wp/norm

      !$acc parallel loop gang vector
      !$omp target teams distribute parallel do simd
      do igpt = 1, length
         this%solar_source(igpt) = this%solar_source(igpt) * tsi * norm
      end do
    end if

  end function set_tsi
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Release device copies and deallocate memory
  !
  subroutine finalize(this)
    class(ty_gas_optics_ddq), intent(inout) :: this

    if(this%is_loaded()) then
      !$acc        exit data delete(this%nus, this%weights)                                   &
      !$acc                  delete(this%fax_a, this%fax_b, this%fax_c, this%fax_sigma0)      &
      !$acc                  delete(this%fax_S, this%fax_vmr0, this%fax_T0, this%fax_p0)      &
      !$acc                  delete(this%xsec_p, this%mtckd_cself, this%mtckd_cfrgn, this%mtckd_n)
      !$omp target exit data map(release:this%nus, this%weights)                              &
      !$omp                  map(release:this%fax_a, this%fax_b, this%fax_c, this%fax_sigma0) &
      !$omp                  map(release:this%fax_S, this%fax_vmr0, this%fax_T0, this%fax_p0) &
      !$omp                  map(release:this%xsec_p, this%mtckd_cself, this%mtckd_cfrgn, this%mtckd_n)
      deallocate(this%nus, this%weights)
      deallocate(this%fax_species_names, this%fax_a, this%fax_b, this%fax_c, this%fax_sigma0)
      deallocate(this%fax_S, this%fax_vmr0, this%fax_T0, this%fax_p0)
      deallocate(this%xsec_species_names, this%xsec_p)
      deallocate(this%mtckd_species_names, this%mtckd_cself, this%mtckd_cfrgn, this%mtckd_n)
      deallocate(this%species_names)

      if(allocated(this%rayleigh_xsec)) then
        !$acc        exit data delete(this%rayleigh_xsec)
        !$omp target exit data map(release:this%rayleigh_xsec)
        deallocate(this%rayleigh_xsec)
      end if
      if(allocated(this%solar_source)) then
        !$acc        exit data delete(this%solar_source)
        !$omp target exit data map(release:this%solar_source)
        deallocate(this%solar_source)
      end if
      !$acc exit data delete(this)
      call this%ty_optical_props%finalize()
    end if
  end subroutine finalize
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if class has been initialized
  !
  pure function is_loaded(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    logical                              :: is_loaded
    is_loaded = allocated(this%nus)
  end function is_loaded
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if initialized for internal sources/longwave, false otherwise
  !
  pure function source_is_internal(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    logical                          :: source_is_internal
    source_is_internal = .not. allocated(this%solar_source)
  end function source_is_internal
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if initialized for external sources/shortwave, false otherwise
  !
  pure function source_is_external(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    logical                          :: source_is_external
    source_is_external = allocated(this%solar_source)
  end function source_is_external

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> How man gases are known to the optics?
  !
  pure function get_ngas(this)
    class(ty_gas_optics_ddq), intent(in)  :: this
    integer :: get_ngas !! names of the gases for which absorption coefficents

    get_ngas = size(this%species_names)
  end function get_ngas
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> For which species are absorption coefficients available?
  !
  pure function get_gases(this)
    class(ty_gas_optics_ddq), intent(in)  :: this
    character(32), dimension(get_ngas(this)) :: get_gases !! names of the gases for which absorption coefficents

    get_gases = this%species_names
  end function get_gases
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the minimum pressure for which absorption coefficient fits are valid
  !
  pure function get_press_min(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_press_min !! minimum pressure for which fits are valid

    get_press_min = 0.01_wp
  end function get_press_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the maximum pressure for which absorption coefficient fits are valid
  !
  pure function get_press_max(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_press_max !! maximum pressure for which fits are valid

    get_press_max = 1100.e2_wp
  end function get_press_max

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the minimum temparature for which absorption coefficient fits are valid
  !
  pure function get_temp_min(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_temp_min !! minimum temperature for which fits are valid

    get_temp_min = 120._wp
  end function get_temp_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the maximum temparature for which absorption coefficient fits are valid
  !
  pure function get_temp_max(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_temp_max !! maximum temperature for which fits are valid

    get_temp_max = 360._wp
  end function get_temp_max
  !--------------------------------------------------------------------------------------------------------------------

end module mo_gas_optics_ddq
