!---------------------------------------------------------------------------------
! This code is to test the routines in the generic_TOPAZ module of
! MOM5       
! Author: Moon,B.-K. (moonbk@jbnu.ac.kr)
!         Jung H. C. (dhtyphoon@naver.com)
! History
!       - 2016-11-08: initial version
!       - 2018-03-08: modular version
!---------------------------------------------------------------------------------

      MODULE generic_TOPAZ_model

      USE ocean_types_mod,   only: ocean_types_init, ocean_grid_type,&
                                   ocean_thickness_type, ocean_domain_type,&
                                   ocean_velocity_type, ocean_time_type,&
                                   ocean_time_steps_type, ocean_options_type,&
                                   ocean_external_mode_type, ocean_prog_tracer_type,&
                                   ocean_diag_tracer_type

      use g_tracer_utils,        only: g_tracer_type, g_diag_type, g_tracer_init,&
                                       g_tracer_set_values, g_tracer_get_pointer,&
                                       g_tracer_get_values, g_tracer_set_common4column

      use time_manager_mod,      only: NOLEAP, JULIAN, set_calendar_type,&
                                       set_date, get_date, set_time, get_time,& 
                                       time_type

      USE coupler_types_mod,     only: coupler_2d_bc_type

      USE netcdf_var
      USE netcdf_io

      IMPLICIT NONE

      !Linked Lists of all prog and diag tracers in this module
      !Ensure these pointers are "save"d between the calls
      type(g_tracer_type), save, pointer         :: tracer_list => NULL()
      type(ocean_prog_tracer_type), dimension(:), pointer, save :: T_prog =>NULL()
      type(ocean_diag_tracer_type), dimension(:), pointer, save :: T_diag =>NULL()
      type(coupler_2d_bc_type),       target, save   :: IOB_struc
      type(time_type),private,        target, save   :: Time_in
      type(ocean_grid_type),          target, save   :: Grid
      type(ocean_thickness_type),     target, save   :: Thickness
      type(ocean_domain_type),        target, save   :: Domain
      type(ocean_time_type),private,  target, save   :: TP_Time
      type(ocean_time_steps_type),private,target, save   :: Time_steps
      type(ocean_options_type),       target, save   :: Ocean_options
      type(ocean_velocity_type),      target, save   :: Velocity
      type(ocean_external_mode_type), target, save   :: Ext_mode

      integer :: index_temp=-1
      integer :: index_salt=-1
      integer :: index_con_temp=-1
      integer :: index_frazil=-1
      integer :: index_irr=-1

      integer :: T_prog_num
      integer :: T_diag_num
      integer, allocatable, dimension(:) :: dim_nc_var_id
      integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, ni, nj, nk
      integer :: si, ei, sj, ej

      real, dimension(:,:,:,:), allocatable      :: diff_cbt
      real, dimension(:,:,:), allocatable        :: opacity
      real, dimension(:,:), allocatable          :: hblt_depth
      real, dimension(:,:), allocatable          :: sw_pen
      real, dimension(:,:), allocatable          :: ST,SS  !!sst, sss
      real, dimension(:,:,:,:), allocatable      :: tp_rho    !!density
      real, dimension(:,:,:), allocatable        :: rad    !!irradiance,radiance flux(w/m2)

      real, dimension(:,:,:,:), allocatable      :: dummy4d
      real, dimension(:,:,:), allocatable        :: dummy3d
      real, dimension(:,:), allocatable          :: dummy2d
      real, dimension(:,:,:), allocatable        :: field3d
      real, dimension(:,:,:), allocatable        :: field2d

      character(len=32)                          :: nc_varname
      real, dimension(:,:), allocatable          :: co2_sc_no !!CO2 Schmitd number
      real, dimension(:,:), allocatable          :: co2_alpha !!CO2 Solubility
      real, dimension(:,:), allocatable          :: co2_csurf !!CO2 sea-surf con
      real, dimension(:,:), allocatable          :: o2_sc_no  !!O2 Schmitd number
      real, dimension(:,:), allocatable          :: o2_alpha  !!O2 Solubility
      real, dimension(:,:), allocatable          :: o2_csurf  !!O2 sea-surf con

      ! setting the number of prognostic and diagnostic tracers  
      ! both determined by reading field_table 
      integer :: num_prog_tracers=-1 ! (e.g., temp, salt, age)
      integer :: num_diag_tracers=-1 ! (e.g., frazil, pH)

      real, dimension(:,:,:), allocatable        :: climf_fed
      real, dimension(:,:,:), allocatable        :: climf_lith
      real, dimension(:,:,:), allocatable        :: climf_no3_wet
      real, dimension(:,:,:), allocatable        :: climf_no3_dry
      real, dimension(:,:,:), allocatable        :: climf_nh4_wet
      real, dimension(:,:,:), allocatable        :: climf_nh4_dry

      real, dimension(:,:,:,:), allocatable      :: clim_chl

      real, public, parameter                    :: epsln = 1.0e-40

      real, dimension(2)                         :: hist_o2

      CONTAINS

!!==================================================================================      
!!==================================================================================
      SUBROUTINE generic_TOPAZ_init(nlev, start_time, dt)

      use fms_io_mod,            only: set_domain

      use ocean_barotropic_mod,  only: ocean_barotropic_init
      use ocean_grids_mod,       only: ocean_grids_init, set_ocean_grid_size
      use ocean_domains_mod,     only: ocean_domain_init, set_ocean_domain,& 
                                       get_local_indices
      use ocean_thickness_mod,   only: ocean_thickness_init
      use ocean_workspace_mod,   only: ocean_workspace_init
      use ocean_grids_mod,       only: set_ocean_hgrid_arrays, set_ocean_vgrid_arrays,& 
                                       init_grids_diag
      use ocean_topog_mod,       only: ocean_topog_init
      use ocean_util_mod,        only: ocean_util_init
      use ocean_obc_mod,         only: ocean_obc_init, ocean_obc_end
      use ocean_operators_mod,   only: ocean_operators_init
      use ocean_coriolis_mod,    only: ocean_coriolis_init
      use ocean_velocity_mod,    only: ocean_velocity_init

      use diag_manager_mod,      only: diag_manager_init, diag_manager_end

      use ocean_tracer_mod,      only: ocean_prog_tracer_init, ocean_diag_tracer_init
      use ocean_tracer_util_mod, only: ocean_tracer_util_init

      use ocean_generic_mod,     only: ocean_generic_init, ocean_generic_flux_init
      use generic_tracer,        only: generic_tracer_get_list, generic_tracer_get_diag_list

      use generic_TOPAZ,         only: user_allocate_arrays, user_deallocate_arrays


      IMPLICIT NONE

      character(len=*), intent(in) :: start_time
      integer,          intent(in) :: nlev
      real,             intent(in) :: dt
      character :: dummy
      integer :: mon, day
      integer :: layout(2)=(/1,1/)
      integer :: io_layout(2)=(/0,0/)
      logical :: debug = .false.
      logical :: have_obc              =.false.
      integer :: horz_grid=1
      logical :: use_blobs             =.false.
      logical :: use_velocity_override =.false.
      logical :: cmip_units            =.false.
      integer :: vert_coordinate
      integer :: vert_coordinate_class
      integer :: vert_coordinate_type
      integer :: n, ntau, i, j, k, t
      integer :: axes(3)
      type(time_type) :: init_time
      real    :: dtime_t
      integer, dimension(:,:), allocatable    :: grid_kmt, mask_coast
      real, dimension(:,:,:), allocatable  :: grid_tmask

      real, dimension(:,:,:), allocatable   :: get_3dvar
      real, dimension(:,:,:,:), allocatable :: get_4dvar

      vert_coordinate = 2; vert_coordinate_class = 1;
      vert_coordinate_type = 1

      call set_calendar_type(JULIAN)
      TP_Time%calendar  = JULIAN; TP_Time%taum1 = 3; TP_Time%tau = 3; TP_Time%taup1= 1
      TP_Time%init       = .true.
      TP_Time%Time_step = set_time(INT(dt));Time_steps%dtts=dt
      ntau=3

      read(start_time(6:7),*) mon
      read(start_time(9:10),*) day
      Time_in=set_date(2,mon,day,0)

      call ocean_grids_init(vert_coordinate, vert_coordinate_class, horz_grid, debug=debug)
      call set_ocean_grid_size(Grid, 'INPUT/grid_spec.nc') !!here Grid set
      call ocean_domain_init()
      call set_ocean_domain(Domain, Grid, layout=layout, io_layout=io_layout) !!here Domain set
      call set_domain(Domain%domain2d)
      call ocean_workspace_init(Domain, Grid)
      call set_ocean_hgrid_arrays(Domain, Grid)
      call ocean_topog_init(Domain, Grid, 'INPUT/grid_spec.nc', vert_coordinate_type)
      call ocean_obc_init(have_obc, Time_steps, Domain, Grid, Ocean_options, vert_coordinate_type, debug=debug)
      call set_ocean_vgrid_arrays(Domain, Grid, have_obc)
      call ocean_util_init(Domain, Grid)
      call diag_manager_init()
      call init_grids_diag(Grid, TP_Time)

      call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)

      call ocean_types_init()
      call ocean_tracer_util_init(Grid, Domain, use_blobs)
      call ocean_coriolis_init(Grid, Domain, TP_Time, Time_steps, Ocean_options, horz_grid, debug=debug)
      call ocean_velocity_init(Grid, Domain, TP_Time, Time_steps, Ocean_options, Velocity, &
                               horz_grid, have_obc, use_blobs, use_velocity_override, debug=debug)
      call ocean_barotropic_init(Grid, Domain, TP_Time, Time_steps, Ocean_options, Ext_mode, have_obc, &
                                 vert_coordinate, vert_coordinate_class, horz_grid, cmip_units, use_blobs,&
                                 use_velocity_override, debug=debug)
      call ocean_thickness_init(TP_Time, Time_steps, Domain, Grid, Ext_mode, Thickness,&
                                vert_coordinate, vert_coordinate_class, vert_coordinate_type, &
                                use_blobs, dtime_t, debug=debug)
      call ocean_operators_init(Grid, Domain, Thickness, horz_grid)

      isc=1;iec=1;Domain%isc=isc;Domain%iec=iec
      jsc=1;jec=1;Domain%jsc=jsc;Domain%jec=jec
      nk=nlev;Grid%nk=nk
      isd=isc-1;ied=iec+1;Domain%isd=isd;Domain%ied=ied
      jsd=jsc-1;jed=jec+1;Domain%jsd=jsd;Domain%jed=jed

      ! initialize prognostic tracers 
      T_prog => ocean_prog_tracer_init(Grid, Thickness, Ocean_options, Domain, TP_Time, Time_steps, &
                num_prog_tracers, vert_coordinate_type, have_obc, &
                cmip_units, use_blobs, debug=debug)
      ! initialize diagnostic tracers 
      T_diag => ocean_diag_tracer_init(TP_Time, Thickness, vert_coordinate_type, num_diag_tracers,&
                use_blobs)

      allocate(diff_cbt(isd:ied,jsd:jed,nk,2)); diff_cbt = 0.0
      allocate(opacity(isd:ied,jsd:jed,nk)); opacity = 0.0
      allocate(hblt_depth(isd:ied,jsd:jed)); hblt_depth = 0.0
      allocate(sw_pen(isd:ied,jsd:jed));     sw_pen = 0.0
      allocate(ST(isd:ied,jsd:jed)); ST = 0.0
      allocate(SS(isd:ied,jsd:jed)); SS = 0.0
      allocate(tp_rho(isd:ied,jsd:jed,nk,3)); tp_rho = 0.0
      allocate(rad(isd:ied,jsd:jed,nk)); rad = 0.0

      allocate(co2_sc_no(isd:ied,jsd:jed)); co2_sc_no=0.0
      allocate(co2_alpha(isd:ied,jsd:jed)); co2_alpha=0.0
      allocate(co2_csurf(isd:ied,jsd:jed)); co2_csurf=0.0
      allocate(o2_sc_no(isd:ied,jsd:jed)); o2_sc_no=0.0
      allocate(o2_alpha(isd:ied,jsd:jed)); o2_alpha=0.0
      allocate(o2_csurf(isd:ied,jsd:jed)); o2_csurf=0.0

      allocate(dummy4d(isd:ied,jsd:jed,nk,3)); dummy4d=0.0
      allocate(dummy3d(isd:ied,jsd:jed,nk)); dummy3d = 0.0
      allocate(dummy2d(isd:ied,jsd:jed)); dummy2d= 0.0
      allocate(field3d(isc:iec,jsc:jec,nk)); field3d = 0.0
      allocate(field2d(isc:iec,jsc:jec,1)); field2d=0.0

      allocate(get_3dvar(isc:iec,jsc:jec,1)); get_3dvar=0.0
      allocate(get_4dvar(isc:iec,jsc:jec,nk,1)); get_4dvar=0.0

      allocate(climf_fed(isc:iec,jsc:jec,12)); climf_fed=0.0
      allocate(climf_lith(isc:iec,jsc:jec,12)); climf_lith=0.0
      allocate(climf_no3_wet(isc:iec,jsc:jec,12)); climf_no3_wet=0.0
      allocate(climf_no3_dry(isc:iec,jsc:jec,12)); climf_no3_dry=0.0
      allocate(climf_nh4_wet(isc:iec,jsc:jec,12)); climf_nh4_wet=0.0
      allocate(climf_nh4_dry(isc:iec,jsc:jec,12)); climf_nh4_dry=0.0
      allocate(clim_chl(isc:iec,jsc:jec,nk,12));clim_chl=0.0
!!===================================INIT TOPAZ======================================================!!
      nc_fname=TRIM('./INPUT/'//'ocean_topaz_column.res.nc')
      call NF_OPEN

      T_prog_num=0
      do n= 1, num_prog_tracers

          if (T_prog(n)%name == 'temp') then
            index_temp = n
            cycle
          end if
          if (T_prog(n)%name == 'salt') then
            index_salt = n
            cycle
          end if

          get_4dvar=0.0 
          CALL GET_VAR4D(T_prog(n)%name, get_4dvar)

          if (allocated(T_prog(n)%field)) deallocate(T_prog(n)%field);allocate(T_prog(n)%field(isd:ied,jsd:jed,nk,3));T_prog(n)%field=0.0
          T_prog(n)%field(isc:iec,jsc:jec,:,TP_Time%taup1)=get_4dvar(isc:iec,jsc:jec,:,1)

          if (allocated(T_prog(n)%stf)) deallocate(T_prog(n)%stf);allocate(T_prog(n)%stf(isd:ied,jsd:jed)); T_prog(n)%stf=0.0
          if (allocated(T_prog(n)%th_tendency))deallocate(T_prog(n)%th_tendency);allocate(T_prog(n)%th_tendency(isd:ied,jsd:jed,nk));T_prog(n)%th_tendency=0.0

          T_prog_num=T_prog_num+1
      enddo

      T_diag_num=0
      do n= 1, num_diag_tracers
      
         if (T_diag(n)%name == 'con_temp') then
           index_con_temp = n
           cycle
         end if
         if (T_diag(n)%name == 'frazil') then
           index_frazil = n
           cycle
         end if
         if (T_diag(n)%name == 'irr' ) then
           index_irr = n
           cycle
         end if

         get_4dvar=0.0
         CALL GET_VAR4D(T_diag(n)%name, get_4dvar)

         deallocate(T_diag(n)%field);allocate(T_diag(n)%field(isd:ied,jsd:jed,nk));T_diag(n)%field=0.0
         T_diag(n)%field(isc:iec,jsc:jec,:)=get_4dvar(isc:iec,jsc:jec,:,1)

         T_diag_num=T_diag_num+1
      enddo

      CALL NF_CLOSE

      allocate(dim_nc_var_id(T_diag_num+T_prog_num+4)) !T_diag + T_prog + irrad + hblt_depth + diff_cbt + tp_rho

      if (allocated(Thickness%rho_dzt)) deallocate(Thickness%rho_dzt);allocate(Thickness%rho_dzt(isd:ied,jsd:jed,nk,3));Thickness%rho_dzt=0.0
      if (allocated(Thickness%rho_dztr)) deallocate(Thickness%rho_dztr);allocate(Thickness%rho_dztr(isd:ied,jsd:jed,nk));Thickness%rho_dztr=0.0
      if (allocated(Thickness%dzt)) deallocate(Thickness%dzt);allocate(Thickness%dzt(isd:ied,jsd:jed,nk)); Thickness%dzt=0.0
      if (allocated(Thickness%dzwt)) deallocate(Thickness%dzwt);allocate(Thickness%dzwt(isd:ied,jsd:jed,0:nk)); Thickness%dzwt=0.0
      if (allocated(Velocity%current_wave_stress))deallocate(Velocity%current_wave_stress);allocate(Velocity%current_wave_stress(isd:ied,jsd:jed));Thickness%dzwt=0.0

!      dummy2d=0.
!      dummy2d(isc:iec,jsc:jec)=Grid%dat(si:ei,sj:ej)
!      if (allocated(Grid%dat)) deallocate(Grid%dat);allocate(Grid%dat(isd:ied,jsd:jed)); Grid%dat=dummy2d
      allocate(mask_coast(isd:ied,jsd:jed)); mask_coast=0
      allocate(grid_kmt(isd:ied,jsd:jed));grid_kmt=Grid%nk
      allocate(grid_tmask(isd:ied,jsd:jed,nk));grid_tmask=1.

      call g_tracer_set_common4column(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,mask_coast,grid_kmt,init_time)
      call user_deallocate_arrays
      call user_allocate_arrays

      call generic_tracer_get_list(tracer_list)

      !!====read init. flux coef of co2, o2
      get_3dvar=0.0
      nc_fname=TRIM('./INPUT/flux_clim/'//&
                    'ocean_topaz_airsea_flux_column.res.nc')
      call NF_OPEN

      nc_varname='o2_flux_alpha_ocn'
      call GET_VAR3D(nc_varname,get_3dvar)
      o2_alpha(isc:iec,jsc:jec)=get_3dvar(isc:iec,jsc:jec,1)

      nc_varname='o2_flux_sc_no_ocn'
      call GET_VAR3D(nc_varname,get_3dvar)
      o2_sc_no(isc:iec,jsc:jec)=get_3dvar(isc:iec,jsc:jec,1)

      nc_varname='o2_flux_csurf_ocn'
      call GET_VAR3D(nc_varname,get_3dvar)
      o2_csurf(isc:iec,jsc:jec)=get_3dvar(isc:iec,jsc:jec,1)

      nc_varname='co2_flux_alpha_ocn'
      call GET_VAR3D(nc_varname,get_3dvar)
      co2_alpha(isc:iec,jsc:jec)=get_3dvar(isc:iec,jsc:jec,1)

      nc_varname='co2_flux_sc_no_ocn'
      call GET_VAR3D(nc_varname,get_3dvar)
      co2_sc_no(isc:iec,jsc:jec)=get_3dvar(isc:iec,jsc:jec,1)

      nc_varname='co2_flux_csurf_ocn'
      call GET_VAR3D(nc_varname,get_3dvar)
      co2_csurf(isc:iec,jsc:jec)=get_3dvar(isc:iec,jsc:jec,1)

      !fed(dissolved fe)
      nc_fname=TRIM('./INPUT/flux_clim/'//&
                    'Soluble_Fe_Flux_PI_column.nc')
      call NF_OPEN
      nc_varname='FLUX'
      call GET_VAR3D(nc_varname,climf_fed)
      call NF_CLOSE

      !lith
      nc_fname=TRIM('./INPUT/flux_clim/'//&
                    'Mineral_Fe_Flux_PI_column.nc')
      call NF_OPEN
      nc_varname='FLUX'
      call GET_VAR3D(nc_varname,climf_lith)
      call NF_CLOSE

      !no3, nh4, alk
      nc_fname=TRIM('./INPUT/flux_clim/'//&
                    'depflux_total.mean.1860_column.nc')
      call NF_OPEN
      nc_varname='NO3_WET_DEP'
      call GET_VAR3D(nc_varname,climf_no3_wet)
      nc_varname='NO3_DRY_DEP'
      call GET_VAR3D(nc_varname,climf_no3_dry)
      nc_varname='NH4_WET_DEP'
      call GET_VAR3D(nc_varname,climf_nh4_wet)
      nc_varname='NH4_DRY_DEP'
      call GET_VAR3D(nc_varname,climf_nh4_dry)    
      call NF_CLOSE

      !chlorophyll
      nc_fname=TRIM('./INPUT/ocean_topaz_chl_clim.nc')
      call NF_OPEN
      nc_varname='chl'
      call GET_VAR4D(nc_varname,clim_chl)
      call NF_CLOSE

      !open output file
      nc_fname=TRIM('./ocean_topaz_tracer.nc')
      call NF_CREATE

      nx=iec; ny=jec; nz=nk
      write(times,100) 'seconds', trim(start_time)
      100   format(A,' since ',A)
      call DEF_DIM
      call DEF_TIME
      call DEF_DEPTH

      dim_nc_var_id=0
      do n = 1, num_prog_tracers
        varname_in=T_prog(n)%name
        long_name=T_prog(n)%longname
        short_name=T_prog(n)%name
        units=T_prog(n)%units
        call DEF_VAR(dim_nc_var_id(n))
        call PUT_ATT(dim_nc_var_id(n))
      end do

      do n = 1, num_diag_tracers
        if ( n .ne. index_con_temp .and.&
             n .ne. index_frazil .and.&
             n .ne. index_irr ) then
           varname_in=T_diag(n)%name
           long_name=T_diag(n)%longname
           short_name=T_diag(n)%name
           units=T_diag(n)%units
           call DEF_VAR(dim_nc_var_id(n+T_prog_num))
           call PUT_ATT(dim_nc_var_id(n+T_prog_num))
        end if
      end do

      varname_in='PAR'
      long_name='PAR'
      short_name='par'
      units='W/m2'
      call DEF_VAR(dim_nc_var_id(1+T_prog_num+T_diag_num))
      call PUT_ATT(dim_nc_var_id(1+T_prog_num+T_diag_num))

      varname_in='hblt_depth'
      long_name='ocean mixed layer thickness dfined by mixing scheme'
      short_name='hblt_depth'
      units='m'
      call DEF_VAR(dim_nc_var_id(2+T_prog_num+T_diag_num))
      call PUT_ATT(dim_nc_var_id(2+T_prog_num+T_diag_num))

      varname_in='diff_cbt'
      long_name='vertical diffusivity for temperature'
      short_name='diff_cbt'
      units='m'
      call DEF_VAR(dim_nc_var_id(3+T_prog_num+T_diag_num))
      call PUT_ATT(dim_nc_var_id(3+T_prog_num+T_diag_num))

      varname_in='rho'
      long_name='density'
      short_name='rho'
      units='kg/m3'
      call DEF_VAR(dim_nc_var_id(4+T_prog_num+T_diag_num))
      call PUT_ATT(dim_nc_var_id(4+T_prog_num+T_diag_num))

      call END_DEF

      END SUBROUTINE generic_TOPAZ_init

!!===============================================================================================


!!===============================================================================================
      subroutine generic_TOPAZ_column_physics(itt, MaxN, nsave, z_r)

      use ocean_generic_mod,     only: ocean_generic_column_physics

      use generic_tracer, only: generic_tracer_get_list
      use generic_tracer, only: generic_tracer_get_diag_list

      use generic_TOPAZ,         only: generic_TOPAZ_update_from_coupler
      use generic_TOPAZ,         only: generic_TOPAZ_update_from_bottom
      use generic_TOPAZ,         only: generic_TOPAZ_set_boundary_values
      use generic_TOPAZ,         only: user_allocate_arrays, user_deallocate_arrays
      use generic_TOPAZ,         only: generic_TOPAZ_register_diag

      use time_manager_mod,      only: get_date
      use time_interp_mod

      use bio_var,               only: h, dzwt, t, s, rho, par, nuh, wind, bio_airp, sfl
      use kpp,                   only: zsbl

      implicit none

      integer, intent(in)             :: itt, MaxN, nsave
      real, dimension(:), intent(in)  :: z_r

      integer             :: i, j, k, n
      integer             :: in_itt
      real                :: w2,w1, secs
      integer             :: year1, year2, month1, month2

      real, dimension(:), allocatable          :: depth

      real, dimension(:,:), allocatable        :: flux_fed
      real, dimension(:,:), allocatable        :: flux_lith
      real, dimension(:,:), allocatable        :: flux_no3_wet, flux_no3_dry
      real, dimension(:,:), allocatable        :: flux_nh4_wet, flux_nh4_dry
      real, dimension(:,:), allocatable        :: flux_o2, flux_co2
      real, dimension(:,:), allocatable        :: cair, pvel, spres, wind10

      !--------------------------------------------------------------------------

      allocate(depth(Grid%nk));depth=0.0
!!=====================================calc flux=============================
      allocate(flux_fed(isd:ied,jsd:jed));flux_fed=0.0
      allocate(flux_lith(isd:ied,jsd:jed));flux_lith=0.0
      allocate(flux_no3_wet(isd:ied,jsd:jed));flux_no3_wet=0.0
      allocate(flux_no3_dry(isd:ied,jsd:jed));flux_no3_dry=0.0
      allocate(flux_nh4_wet(isd:ied,jsd:jed));flux_nh4_wet=0.0
      allocate(flux_nh4_dry(isd:ied,jsd:jed));flux_nh4_dry=0.0
      allocate(flux_o2(isd:ied,jsd:jed));flux_o2=0.0
      allocate(flux_co2(isd:ied,jsd:jed));flux_co2=0.0
      allocate(cair(isd:ied,jsd:jed));cair=0.0 ! (mol/mol) 
      allocate(pvel(isd:ied,jsd:jed));pvel=0.0 !piston velocity (m/s)
      allocate(spres(isd:ied,jsd:jed));spres=0.0 !surf. pressure (Pa)
      allocate(wind10(isd:ied,jsd:jed));wind10=0.0 !10m wind (m/s)


      !call generic_TOPAZ_update_from_coupler(tracer_list)

      !Time_in = Time_in + TP_Time%Time_step
      call time_interp(Time_in,w2,year1,year2,month1,month2)
      w1 = 1.0-w2
!      print*,'clim flux weight:: ', month2, w2,' / ',month1, w1
        
      flux_fed(isc:iec,jsc:jec)=climf_fed(isc:iec,jsc:jec,month1)*w1+climf_fed(isc:iec,jsc:jec,month2)*w2
      flux_lith(isc:iec,jsc:jec)=climf_lith(isc:iec,jsc:jec,month1)*w1+climf_lith(isc:iec,jsc:jec,month2)*w2
      flux_no3_wet(isc:iec,jsc:jec)=climf_no3_wet(isc:iec,jsc:jec,month1)*w1+climf_no3_wet(isc:iec,jsc:jec,month2)*w2
      flux_no3_dry(isc:iec,jsc:jec)=climf_no3_dry(isc:iec,jsc:jec,month1)*w1+climf_no3_dry(isc:iec,jsc:jec,month2)*w2
      flux_nh4_wet(isc:iec,jsc:jec)=climf_nh4_wet(isc:iec,jsc:jec,month1)*w1+climf_nh4_wet(isc:iec,jsc:jec,month2)*w2
      flux_nh4_dry(isc:iec,jsc:jec)=climf_nh4_dry(isc:iec,jsc:jec,month1)*w1+climf_nh4_dry(isc:iec,jsc:jec,month2)*w2

      call g_tracer_set_values(tracer_list,'lith','stf',flux_lith,isd,jsd)
      call g_tracer_set_values(tracer_list,'fed','stf',flux_fed,isd,jsd)
      call g_tracer_set_values(tracer_list,'nh4','stf',flux_nh4_wet+flux_nh4_dry,isd,jsd)
      call g_tracer_set_values(tracer_list,'no3','stf',flux_no3_wet+flux_no3_dry,isd,jsd)
      call g_tracer_set_values(tracer_list,'alk','stf',-flux_no3_wet-flux_no3_dry,isd,jsd) 

      if ( itt .eq. 1 ) then
       call g_tracer_set_values(tracer_list,'dic','alpha',co2_alpha,isd,jsd)   
       call g_tracer_set_values(tracer_list,'dic','csurf',co2_csurf,isd,jsd)
       call g_tracer_set_values(tracer_list,'dic','sc_no',co2_sc_no,isd,jsd)
       call g_tracer_set_values(tracer_list,'o2','alpha',o2_alpha,isd,jsd)
       call g_tracer_set_values(tracer_list,'o2','csurf',o2_csurf,isd,jsd)
       call g_tracer_set_values(tracer_list,'o2','sc_no',o2_sc_no,isd,jsd)
      end if

      call g_tracer_get_values(tracer_list,'dic','alpha',co2_alpha,isd,jsd)
      call g_tracer_get_values(tracer_list,'dic','csurf',co2_csurf,isd,jsd)
      call g_tracer_get_values(tracer_list,'dic','sc_no',co2_sc_no,isd,jsd)
      call g_tracer_get_values(tracer_list,'o2','alpha',o2_alpha,isd,jsd)
      call g_tracer_get_values(tracer_list,'o2','csurf',o2_csurf,isd,jsd)
      call g_tracer_get_values(tracer_list,'o2','sc_no',o2_sc_no,isd,jsd)


      do j = jsc, jec
      do i = isc, iec
        wind10(i,j)=wind    !gotm 10m wind(m/s)
        spres(i,j)=bio_airp !gotm air pressure(pa)

        pvel(i,j)=9.36e-007*wind10(i,j)**2
        !o2
        cair(i,j)=o2_alpha(i,j)*0.21*spres(i,j)*9.7561e-006
        flux_o2(i,j)=pvel(i,j)*sqrt(660/(o2_sc_no(i,j)+epsln))* &
                                 (o2_csurf(i,j)-cair(i,j))
      end do
      end do
      flux_o2=-flux_o2
      call g_tracer_set_values(tracer_list,'o2','stf',flux_o2,isd,jsd)
     
      do j = jsc, jec
      do i = isc, iec
        !co2
        cair(i,j)=co2_alpha(i,j)*286.0e-6*spres(i,j)*9.7561e-006
        flux_co2(i,j)=pvel(i,j)*sqrt(660/(co2_sc_no(i,j)+epsln))* &
                                  (co2_csurf(i,j)-cair(i,j))
      end do
      end do
      flux_co2=-flux_co2
      call g_tracer_set_values(tracer_list,'dic','stf',flux_co2,isd,jsd)

!!=====================================calc flux=============================

!!==================================read ocean physics :: gotm================
        do k = 1, Grid%nk 
        do j = jsc, jec
        do i = isc, isc
           Thickness%dzt(i,j,k)=h(Grid%nk+1-k)
           T_prog(1)%field(i,j,k,TP_Time%taup1)=t(Grid%nk+1-k)
           T_prog(2)%field(i,j,k,TP_Time%taup1)=s(Grid%nk+1-k)
           tp_rho(i,j,k,TP_Time%taup1)=rho(Grid%nk+1-k)
           Thickness%rho_dzt(i,j,k,TP_Time%taup1)=rho(Grid%nk+1-k)*h(Grid%nk+1-k)
           diff_cbt(i,j,k,1)=nuh(Grid%nk+1-k-1)
           rad(i,j,k)=par(Grid%nk+1-k)
           Thickness%dzwt(i,j,k)=dzwt(Grid%nk+1-k)           
           hblt_depth(i,j)=zsbl*(-1.0)
        enddo 
        enddo
           depth(k)=z_r(Grid%nk+2-k)
        enddo
!!==================================read ocean physics :: gotm================
!        print*,'Thickness%rho_dzt...',Thickness%rho_dzt(1,1,:,TP_Time%taup1)
!        print*,'Thickness%dzt...',Thickness%dzt(1,1,:)
!        print*,'Thickness%dzwt...',Thickness%dzwt(1,1,:)
!        print*,'T_prog(1)%field...',T_prog(1)%field(1,1,:,TP_Time%taup1)
!        print*,'T_prog(2)%field...',T_prog(2)%field(1,1,:,TP_Time%taup1)
!        print*,'diff_cbt...',diff_cbt(1,1,:,1)
!        print*,'hblt_depth...',hblt_depth(1,1)
!        print*,'Dens%rho...',tp_rho(1,1,:,TP_Time%taup1)
!        print*,'rad...',rad(1,1,:)

      call g_tracer_get_values(tracer_list,'dic','alpha',co2_alpha,isd,jsd)   
      call g_tracer_get_values(tracer_list,'dic','csurf',co2_csurf,isd,jsd)
      call g_tracer_get_values(tracer_list,'dic','sc_no',co2_sc_no,isd,jsd)
      call g_tracer_get_values(tracer_list,'o2','alpha',o2_alpha,isd,jsd)
      call g_tracer_get_values(tracer_list,'o2','csurf',o2_csurf,isd,jsd)
      call g_tracer_get_values(tracer_list,'o2','sc_no',o2_sc_no,isd,jsd)


!!=================================TOPAZ MAIN ROUTINE================================================
      call ocean_generic_column_physics(Thickness, hblt_depth, TP_Time, &
           Grid, Time_steps%dtts, isd, jsd, T_prog, T_diag, sw_pen, opacity, diff_cbt, Velocity, rad)

       
      if ( itt == 1 ) then
        hist_o2(1)=T_prog(13)%field(1,1,1,TP_Time%taup1)
        call generic_TOPAZ_set_boundary_values(tracer_list,T_prog(1)%field(isd:ied,jsd:jed,1,TP_Time%taup1),&
                                               T_prog(2)%field(isd:ied,jsd:jed,1,TP_Time%taup1),tp_rho,isd,jsd,TP_Time%taup1,hist_o2(1))
      
      else if ( itt == 2 ) then
        hist_o2(2)=T_prog(13)%field(1,1,1,TP_Time%taup1)
        call  generic_TOPAZ_set_boundary_values(tracer_list,T_prog(1)%field(isd:ied,jsd:jed,1,TP_Time%taup1),&
                                               T_prog(2)%field(isd:ied,jsd:jed,1,TP_Time%taup1),tp_rho,isd,jsd,TP_Time%taum1,hist_o2(1))
      else if ( mod(itt,2) == 1 ) then
        call  generic_TOPAZ_set_boundary_values(tracer_list,T_prog(1)%field(isd:ied,jsd:jed,1,TP_Time%taup1),&
                                               T_prog(2)%field(isd:ied,jsd:jed,1,TP_Time%taup1),tp_rho,isd,jsd,TP_Time%taum1,hist_o2(1))
        hist_o2(1)=T_prog(13)%field(1,1,1,TP_Time%taup1)
      else if ( mod(itt,2) == 0 ) then
        call generic_TOPAZ_set_boundary_values(tracer_list,T_prog(1)%field(isd:ied,jsd:jed,1,TP_Time%taup1),&
                                               T_prog(2)%field(isd:ied,jsd:jed,1,TP_Time%taup1),tp_rho,isd,jsd,TP_Time%taum1,hist_o2(2))
        hist_o2(2)=T_prog(13)%field(1,1,1,TP_Time%taup1)
      end if
!!=================================TOPAZ MAIN ROUTINE================================================
      !------------------------------------------
      ! how to get values from topaz
      !------------------------------------------
      !call g_tracer_get_values(tracer_list,'o2','stf', dummy3d ,isd, jsd)
      !call g_tracer_get_pointer(tracer_list,'o2','stf',stf)


      if ( mod(itt,nsave) .eq. 0 ) then
      !------------------------------------------
      ! write topaz output
      !------------------------------------------
      secs= itt * Time_steps%dtts
      in_itt=itt/nsave

      call PUT_TIME(in_itt, secs)
      call PUT_DEPTH(depth)

      do n = 1, num_prog_tracers
        if ( n .ne. index_temp .and.&
             n .ne. index_salt ) then
          call g_tracer_get_values(tracer_list,T_prog(n)%name,'field', dummy3d ,isc,jsc)
          field3d(:,:,:) = dummy3d(1:iec,1:jec,:)
          call PUT_VAR(dim_nc_var_id(n), field3d, nk, in_itt)
        else !( .eq. index_temp & index_salt)
          field3d(:,:,:) = T_prog(n)%field(1:iec,1:jec,:,TP_Time%taup1)
          call PUT_VAR(dim_nc_var_id(n),field3d, nk, in_itt)
        end if
      enddo

      do n = 1, num_diag_tracers
        if ( n .ne. index_con_temp .and.&
             n .ne. index_frazil .and.&
             n .ne. index_irr ) then
          call g_tracer_get_values(tracer_list,T_diag(n)%name,'field', dummy3d ,isc,jsc)
          field3d(:,:,:) = dummy3d(1:iec,1:jec,:)
          call PUT_VAR(dim_nc_var_id(n+T_prog_num), field3d, nk, in_itt)    
        endif
      enddo

      field3d(1,1,:)=rad(1,1,:)
      call PUT_VAR(dim_nc_var_id(1+T_prog_num+T_diag_num),field3d, nk, in_itt)
      field2d(1,1,1)=hblt_depth(1,1) 
      call PUT_VAR(dim_nc_var_id(2+T_prog_num+T_diag_num),field2d, nk, in_itt)
      field3d(1,1,:)=diff_cbt(1,1,:,1)
      call PUT_VAR(dim_nc_var_id(3+T_prog_num+T_diag_num),field3d, nk, in_itt)
      field3d(1,1,:)=tp_rho(1,1,:,TP_Time%taup1)
      call PUT_VAR(dim_nc_var_id(4+T_prog_num+T_diag_num),field3d, nk, in_itt)

      end if !if ( mod(itt,nsave) .eq. 0 ) then

      if ( itt == MaxN ) call NF_CLOSE

      end subroutine generic_TOPAZ_column_physics

!===================================================================================================
      SUBROUTINE topaz_optic(z_w)

      USE bio_var,   only : h, rad, bioshade_, par, rho, topaz_climchl
      USE time_interp_mod
      USE observations, only: A

      implicit none

      real, dimension(:), intent(in)  :: z_w
      real, allocatable, dimension(:) :: depth_w

      real, allocatable, dimension(:) :: sw_frac_zt, sw_frac_zw, opacity
      real, allocatable, dimension(:) :: red, blue
      real :: sw_irr, f_vis
      integer :: n, k, i
      real :: w2, w1
      integer :: year1, year2, month1, month2
      real, allocatable, dimension(:) :: chl, tmp_chl

      if (.not. allocated(depth_w)) allocate(depth_w(Grid%nk));depth_w=0.0
      if (.not. allocated(chl)) allocate(chl(Grid%nk));chl=0.0
      if (.not. allocated(tmp_chl)) allocate(tmp_chl(Grid%nk));chl=0.0
      if (.not. allocated(sw_frac_zt)) allocate(sw_frac_zt(Grid%nk));sw_frac_zt=0.0
      if (.not. allocated(sw_frac_zw)) allocate(sw_frac_zw(Grid%nk));sw_frac_zw=0.0
      if (.not. allocated(opacity)) allocate(opacity(Grid%nk));opacity=0.0
      if (.not. allocated(red)) allocate(red(Grid%nk));red=0.0
      if (.not. allocated(blue)) allocate(blue(Grid%nk));blue=0.0

      Time_in = Time_in + TP_Time%Time_step
      
      if ( topaz_climchl ) then
        call time_interp(Time_in,w2,year1,year2,month1,month2)
        w1 = 1.0 - w2
!        print*,'chl clim weight:: ',month2, w2, ' / ', month1, w1
        do k = Grid%nk,1,-1
          chl(k)=clim_chl(isc,jsc,k,month1)*w1+clim_chl(isc,jsc,k,month2)*w2
          chl(k)=chl(k)/1000.*rho(k)        !!ug/kg -> mg/m3
        end do
      else
        do n= 1, num_diag_tracers
          if ( T_diag(n)%name .eq. 'chl' ) then
             do k = Grid%nk,1,-1
               chl(k)=T_diag(n)%field(isc,jsc,Grid%nk+1-k)/1000.*rho(k) !!ug/kg -> mg/m3
             end do
          end if
        end do
      end if


      f_vis=1.0-A
      i=Grid%nk
      sw_irr=rad(i)*f_vis
      red(i)=0.5*exp(-h(i)*(0.225+(0.037*chl(i)**0.629)))
      blue(i)=0.5*exp(-h(i)*(0.0232+(0.074*chl(i)**0.674)))

      sw_frac_zw(i)=f_vis*(red(i)+blue(i))
      sw_frac_zt(i)=0.5*(f_vis+sw_frac_zw(i))
      opacity(i)=-log(sw_frac_zw(i)/f_vis)/h(i)

      par(i)=sw_irr*exp(-opacity(i)*h(i)*0.5)
      sw_irr=sw_irr*exp(-opacity(i)*h(i))

      if ( rad(Grid%nk) .ne. 0.0 ) then
        bioshade_(i)=sw_irr/rad(Grid%nk)
      end if

      do i = Grid%nk-1,1,-1
       if ( z_w(i+1) <= 200.0 ) then
        red(i)=red(i+1)*exp(-h(i)*(0.2250+(0.037*chl(i)**0.629)))
        blue(i)=blue(i+1)*exp(-h(i)*(0.0232+(0.074*chl(i)**0.674)))

        sw_frac_zw(i)=f_vis*(red(i)+blue(i))
        sw_frac_zt(i)=0.5*(sw_frac_zw(i+1)+sw_frac_zw(i))
        opacity(i)=-log(sw_frac_zw(i)/sw_frac_zw(i+1))/h(i)

        par(i)=sw_irr*exp(-opacity(i)*h(i)*0.5)
        sw_irr=sw_irr*exp(-opacity(i)*h(i))
       else
        par(i)=0.0
        sw_irr=0.0
       end if

       if ( rad(Grid%nk) .ne. 0.0 ) then
         bioshade_(i)=sw_irr/rad(Grid%nk)
       end if

      end do

      END SUBROUTINE topaz_optic
!===================================================================================================

      SUBROUTINE topaz_w_adv(w_adv_discr)

      USE bio_var,    only: h, w
      USE util,       only: flux

      IMPLICIT NONE

      REAL, DIMENSION(:,:), ALLOCATABLE  :: prog
      INTEGER                            :: i, k
      INTEGER                            :: w_adv_discr
      INTEGER, PARAMETER                 :: adv_mode_0=0

      ALLOCATE(prog(num_prog_tracers,Grid%nk))

      do i = 1, num_prog_tracers
      do k = 1, Grid%nk
        prog(i,k) = T_prog(i)%field(1,1,Grid%nk+1-k,TP_Time%taup1)          
      end do
      end do

      do i = 1, num_prog_tracers
        CALL adv_center(Grid%nk,Time_steps%dtts,h,h,w,&
                        flux,flux,0.0,0.0,w_adv_discr,adv_mode_0,prog(i,:))
      end do

      do i =  1, num_prog_tracers
      do k = 1, Grid%nk
        T_prog(i)%field(1,1,k,TP_Time%taup1)=prog(i,Grid%nk+1-k)
      end do
      end do

      END SUBROUTINE topaz_w_adv



      END MODULE generic_TOPAZ_model
