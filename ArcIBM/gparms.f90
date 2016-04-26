Module gparms

!=======================================================================
! Fiscm/BISCM Global Parameters
!
! Description
!   Defines global params (constants) for the fiscm code
!
! Comments:
!    Originally from the OSCAR shallow water equation solver
!
! !REVISION HISTORY:
!  Original author(s): G. Cowles
!  2011/7/28       Xueming Zhu
!  2014/09/22      Zhixuan Feng
!  2015/11/13      Zhixuan Feng: Add "food_source" to implement various strategies for copepod food sources
!=======================================================================

    Implicit None

!single precision kind
!integer, parameter :: sp = selected_real_kind(6 , 37)

!double precision kind
    integer, parameter  :: sp = selected_real_kind(15,307)
    integer, parameter  :: dp = selected_real_kind(12,300)
!double precision constants
    real(sp), parameter :: an8th = 0.125_sp
    real(sp), parameter :: a4th  = 0.250_sp
    real(sp), parameter :: a3rd  = 1.0_sp/3.0_sp
    real(sp), parameter :: ahalf = 0.500_sp
    real(sp), parameter :: zero  = 0.000_sp
    real(sp), parameter :: one   = 1.000_sp
    real(sp), parameter :: eps   = 1.0E-6

! define parameters for N-D arrays (eg, 2-D,3-D,4-D) zfeng 10/07/2014
    integer, parameter  :: DIM1 = 1
    integer, parameter  :: DIM2 = 2
    integer, parameter  :: DIM3 = 3
    integer, parameter  :: DIM4 = 4

!max arrays sizes
    integer, parameter  :: max_state_vars = 200

!----------------------------------------------------------------
!string length
!    fstr:    filename length
!    vstr:    variable name length
!    sstr:    short string length
!    tstr:    text string
!    cstr:    string variable standard length
!    mstr:    message length
!----------------------------------------------------------------
    integer, parameter  :: fstr = 120
    integer, parameter  :: tstr = 120
    integer, parameter  :: vstr = 30
    integer, parameter  :: sstr = 15
    integer, parameter  :: cstr = 30
    integer, parameter  :: mstr = 120
!----------------------------------------------------------------
!time
!     day_2_sec:    convert days to seconds
!     sec_2_day:    convert seconds to days
!----------------------------------------------------------------

    real(sp), parameter :: day_2_sec = 86400.0_sp
    real(sp), parameter :: sec_2_day = one/day_2_sec

!----------------------------------------------------------------
!statistical
!      rvar:  variance of a uniform random walk [-1,1]
!----------------------------------------------------------------
    real(sp), parameter :: rvar = a3rd

!----------------------------------------------------------------
!trigonometric
!      pi:    pi
!     d2r:    convert degrees to radians
!     r2d:    convert radians to degrees
!----------------------------------------------------------------

    real(sp), parameter :: pi  = 3.141592653589793238_sp
    real(sp), parameter :: d2r = pi/180.0_sp
    real(sp), parameter :: r2d = 180.0_sp/pi

!----------------------------------------------------------------
!oceanic parameters
!    gacc        :  gravitational acceleration   [ms^-2]
!    omega_earth :  earths rotation rate         [s^-1]
!----------------------------------------------------------------
    real(sp), parameter :: gacc  = 9.8016_sp !default
    real(sp), parameter :: omega_earth = 7.292e-5_sp
    real(sp), parameter :: earth =   6371.004E03_sp

!----------------------------------------------------------------
!large and small numbers
!    hugenum = largest float
!    tinynum = smallest float
!----------------------------------------------------------------
    real(sp), parameter :: hugenum = huge(1.0_sp)
    real(sp), parameter :: tinynum = tiny(1.0_sp)

!----------------------------------------------------------------
! control netcdf output of variables
!    NETCDF_YES:  yes, output
!    NETCDF_NO:   no, do not output
!----------------------------------------------------------------
    integer, parameter  :: NETCDF_YES      = 1
    integer, parameter  :: NETCDF_NO       = 0
    integer, parameter  :: NCDO_HEADER     = 0
    integer, parameter  :: NCDO_ADD_STATES = 1
    integer, parameter  :: NCDO_OUTPUT     = 2

!----------------------------------------------------------------
! status of individuals
!    DEAD: -3
!    SETTLED: -2
!    EXITED DOMAIN: -1
!    UNKNOWN: 0
!    ACTIVE: 1
!----------------------------------------------------------------
    integer, parameter  :: DEAD    = -3
    integer, parameter  :: SETTLED = -2
    integer, parameter  :: EXITED  = -1
    integer, parameter  :: UNKNOWN = 0
    integer, parameter  :: ACTIVE  = 1

!----------------------------------------------------------------
! type of number
!   integer   : 1
!   float     : 2
!   logic     : 3
!   character : 4
!----------------------------------------------------------------
    integer, parameter  :: int_type = 1
    integer, parameter  :: flt_type = 2
    integer, parameter  :: log_type = 3
    integer, parameter  :: cha_type = 4

!----------------------------------------------------------------
! type of variables in NETCDF file
!  with time variable      : 1
!  without time variable   : 0
!----------------------------------------------------------------
    integer, parameter  :: TIME_VARYING = 1
    integer, parameter  :: CONSTANT = 0

!----------------------------------------------------------------
! frame status
!     setup  : 0
!     status : 1
!----------------------------------------------------------------
    integer, parameter  :: FRAME_SETUP = 0
    integer, parameter  :: FRAME_STATS = 1

!----------------------------------------------------------------
! diffusion type
!----------------------------------------------------------------
    integer, parameter  :: HDIFF_NONE     = 0
    integer, parameter  :: HDIFF_CONSTANT = 1
    integer, parameter  :: HDIFF_VARIABLE = 2     !unfinished
    integer, parameter  :: VDIFF_NONE     = 0
    integer, parameter  :: VDIFF_VARIABLE = 1
    integer, parameter  :: VDIFF_SPLINED  = 2     !unfinished
    integer, parameter  :: VDIFF_BINNED   = 3     !unfinished

!----------------------------------------------------------------
!-model setup control parameters
!----------------------------------------------------------------
    integer,  parameter :: max_nf = 250
    real(sp), parameter :: tpi  =3.14159265_sp/180.0_sp*6371.0_sp*1000.0_sp
    integer, parameter  :: iunit = 33
    real(sp),allocatable:: zpini(:),zptini(:)
    real(sp)            :: fbeg,fend,gt
    real(sp)            :: beg_time,end_time
    integer             :: sim_direction
    character(len=fstr) :: runcontrol
    logical             :: fexist

!----------------------------------------------------------------
!mesh params
!----------------------------------------------------------------
    integer, parameter  :: level_based = 1
    integer, parameter  :: layer_based = 2

!----------------------------------------------------------------
!Runge-Kutta integration coefficients
!----------------------------------------------------------------
    integer, parameter  :: nstage   = 1
    real(sp),parameter  :: alpha(1) = 1.0_sp
    integer, parameter  :: mstage   = 4
    real(sp),parameter  :: A_RK(4)  = (/0.0_sp,0.5_sp,0.5_sp,1.0_sp/)
    real(sp),parameter  :: B_RK(4)  = (/1.0_sp/6.0_sp,1.0_sp/3.0_sp,  &
                                        1.0_sp/3.0_sp,1.0_sp/6.0_sp/)
    real(sp),parameter  :: C_RK(4)  = (/0.0_sp,0.5_sp,0.5_sp,1.0_sp/)

!----------------------------------------------------------------
! version & var character name
!----------------------------------------------------------------
    character(len=fstr) :: ArcIBM_VERSION= "ArcIBM1.0"
    character(len=11)   :: nele_char,node_char,                    &
                           siglay_char,siglev_char,s_char,         &
                           x_char,y_char,z_char,h_char,nv_char,    &
                           a1u_char,a2u_char,                      &
                           aw0_char,awx_char,awy_char,             &
                           zeta_char,u_char,v_char,omega_char,     &
                           kh_char,ua_char,va_char,                &
                           wu_char,wv_char

!----------------------------------------------------------------
! namelist NML_FISCM : model control parameters
!----------------------------------------------------------------
    real(sp)            :: beg_time_days,end_time_days
    real(sp)            :: deltaT
    integer             :: ireport,ngroups
    integer             :: mjd_offset
    integer             :: year_cycle,nfiles_in
    character(len=fstr) :: forcing_file(max_nf)
!--------------------------------------------------------------------
! Variables below are added for structured forcing  ! (zfeng 09/07/2014)
    integer             :: spherical   ! 0 - x y(m)  ;1 - lon lat(deg)
    logical             :: Structured  ! T(rue)-structured; F(alse)-unstructured
    integer             :: TimeStepping! 1 - 1st-order forward Euler scheme for time stepping
                                       ! 2 - 2nd-order Runge-Kutta scheme for time stepping
                                       ! 3 - 3rd-order Runge-Kutta scheme for time stepping
                                       ! 4 - 4th-order Runge-Kutta scheme for time stepping
!   integer              :: grid_type  ! 1 - Arakawa A grid
                                       ! 2 - Arakawa B grid
                                       ! 3 - Arakawa C grid
                                       ! 4 - Arakawa D grid
    integer             :: nx          ! nx - number of grid points in x-dir (lon)
    integer             :: ny          ! ny - number of grid points in y-dir (lat)
    integer             :: nz          ! nz - number of vertical layers
    integer             :: ndays       ! ndays - number of days of output data
    integer             :: year        ! year - year of netcdf forcing
!--------------------------------------------------------------------
    logical             :: north_pole  ! zhuxm  angle_pre   ! to adjust the physical field coordinate
    integer             :: sz_cor      ! 0 - input s ;1 - input z
    real(sp)            :: Dtop, Dbot  ! top/bottom mixing layer depth
    integer             :: fix_dep     ! 0 - unfixed ;1 - fix(dep)
    integer             :: dvm_bio     ! 0 - nodvm   ;1 - dvm(bio)
    integer             :: wind_type   ! 0 - nowind  ;1 - wind
    integer             :: bio_fd      ! 0 - no food impat ;1 - food-limiting
    integer             :: n_extfile   ! 0 - no extfile ;1 - with extfile
    character(len=fstr) :: extfile_name(max_nf)
    real(sp)            :: dvmh_up,dvmh_dn  ! up-from surface;dn-from bottom
    integer             :: bcond_type  ! 0 - absorbing; 1 - non-absorbing
!    integer             :: velocity_scheme    ! 0 - instantanous velocity; 1 - whole-depth-average; 2 - average of active layers, use ActLayInd 
!    integer             :: temperature_scheme ! 0 - instantanous temp;     1 - whole-depth-average; 2 - average of active layers    

    namelist /NML_FISCM/ &
      beg_time_days,   &
      end_time_days,   &
      mjd_offset,      &
      deltaT,          &
      ireport,         &
      ngroups,         &
      year_cycle,      &
      nfiles_in,       &
      forcing_file,    &
      spherical,       &
      Structured,      & ! zfeng add for structured model input
      TimeStepping,    &
      nx,              & ! zfeng add
      ny,              & ! zfeng add
      nz,              & ! zfeng add
      ndays,           & ! zfeng add
      year,            & ! zfeng add
      north_pole,      &
      sz_cor,          &
      Dtop,            &
      Dbot,            &
      fix_dep,         &
      dvm_bio,         &
      wind_type,       &
      bio_fd,          &
      n_extfile,       &
      extfile_name,    &
      dvmh_up,         &
      dvmh_dn,         &
      bcond_type!,      &
!      velocity_scheme, &
!      temperature_scheme 

!----------------------------------------------------------------
! namelist NML_STATEVAR
!----------------------------------------------------------------
    character(len=fstr) :: state_varname
    character(len=fstr) :: state_longname
    character(len=fstr) :: state_units
    integer             :: state_netcdf_out
    integer             :: state_vartype
    integer             :: state_initval_int
    real(sp)            :: state_initval_flt
    character(len=fstr) :: state_from_ext_var

    Namelist /NML_STATEVAR/      &
      state_varname,             &
      state_longname,            &
      state_units,               &
      state_netcdf_out,          &
      state_vartype,             &
      state_initval_int,         &
      state_initval_flt,         &
      state_from_ext_var

!----------------------------------------------------------------
! namelist NML_GROUP
!----------------------------------------------------------------
    character(len=fstr) :: group_name
    integer             :: Tnind,space_dim
    integer             :: hdiff_type,vdiff_type
    real(sp)            :: hdiff_const_val,vdiff_const_val
    integer             :: vdiff_substeps
    logical             :: biology  ! Activate biology 
    logical             :: bio_food ! Activate food-limitation
    integer             :: food_source 
    ! implement different food sources: 1. ambient phytoplankton; 2. ambient combined phytoplankton and microzooplankton
    !                                   3. water column maximum phytoplankton 4. water column maximum combined phytop and microzoopl   
    real(sp)            :: LowFood   ! threshold value for low food starvation (say, 0.001 mmol-N m-3)
    integer             :: ActLayInd ! ActLayInd is the index of deepest depth that copepod can seek food maxima (use 13 for 100m; 15 for 160 m)
    integer             :: intvl_bio
    real(sp)            :: start_out
    integer             :: intvl_out,nstate_ud
    character(len=fstr) :: ini_file

    Namelist /NML_GROUP/  &
      group_name,         &
      Tnind,              &
      space_dim,          &
      hdiff_type,         &
      hdiff_const_val,    &
      vdiff_type,         &
      vdiff_const_val,    &
      vdiff_substeps,     &
      biology,            &
      bio_food,           &
      food_source,        &
      LowFood,            &
      ActLayInd,          &
      intvl_bio,          &
      start_out,          &
      intvl_out,          &
      nstate_ud,          &
      ini_file

!----------------------------------------------------------------
! namelist NML_COPEPOD
!----------------------------------------------------------------
  real(sp)           :: BEL_ALPHA ,PC0
  integer            :: BEL_A(12)

  Namelist /NML_COPEPOD/  &
    BEL_ALPHA,    &
    BEL_A,        &
    PC0

End Module gparms

