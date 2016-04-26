Program ArcIBM  ! ArcIBM version-1

! Arctic copepod Individual-Based Model (ArcIBM)
! Modified from FISCM (FVCOM i-state Configuration Model) code
! Dr. Zhixuan Feng
! created 09/02/2014
! Modified 11/13/2015
! Email: zfeng@whoi.edu
!
!=======================================================================
! ArcIBM Main Program
!
! Description
! - Main routine for ArcIBM - Setup and Time Looper
!
! - An individual-based model for Arctic copepods driven by
!    structured BIOMAS output in spherical curvilinear grid.
! -  The model is run in "off line" mode and
!    has lagrangian tracking & generic stage-based biological modules
!
! - Original code of Ji version was written for modeling Arctic copepods
!   forced by unstructured FVCOM outcome (velocity & temperature fields).
!
! - For details of previous applications refer to:
! -  Ji, R., Ashjian, C.J., Campbell, R.G., et al, 2012. Life history
!    and biogeography of Calanus copepods in the Arctic Ocean: An
!    individual-based modeling study. Prog. Oceanogr. 96,
!    40-56. doi:10.1016/j.pocean.2011.10.001
!
!
! !REVISION HISTORY:
!  Original author(s): G. Cowles
!  2011/7/28    Xueming Zhu
!  2014/09/02   Zhixuan Feng
!               1.Rename fiscm to ArcIBM because BIOMAS forcing is used
!                 instead of FVCOM.
!               2.Add new 'mod_setup' module and move 'setup' subroutine
!                   to this module
!               3.Add new 'mod_structured_compute' module to be driven
!                 structured model forcing.
!               4.Change 'gparms' module namelist to include input
!                 parameters for structured forcing.
!               5.Write new "mod_biomas_netcdf" module to convert
!                 BIOMAS model files (mesh + data output etc.) into
!                 a single netcdf forcing file, similar to FVCOM forcing.
!                 Generation of netcdf file should take place before IBM runs.
!               6.Add new 'mod_unstructured_compute' module to be driven
!                 unstructured model forcing and move the FVCOM
!                 computational part from main program to this module.
!               7.Add new 'biomas_driver' module which includes all BIOMAS
!                 structured forcing related algorithms.  
!=======================================================================

  use gparms
  use mod_setup
  use mod_igroup
  use mod_driver
  use mod_bio
  use mod_output
  use mod_forcing
  use mod_unstructured_compute
  use mod_structured_compute
  use utilities

  implicit none

  integer             :: n  ! group number
  real(sp)            :: t = 0.0
  integer             :: nvars ,uniqvars

!  write(*,'("Kind for single precision of delmar is",I2)') KIND(0.0)
!  write(*,'("Kind for double precision of delmar is",I2)') KIND(0.0D0)
   call drawline("-")
   write(*,*) 'The ArcIBM version is '//ArcIBM_VERSION
  !===================================================================
  !initialize the model (read namelists, allocate data, set ic's
  !===================================================================

  !----------------------------------------------------------
  ! read primary simulation data and setup groups
  !----------------------------------------------------------
  fbeg = -1.
  fend = -1.
  call setup   ! in a separate module "mod_setup"

  !----------------------------------------------------------
  ! initialize biology
  !  - set spawning times
  !  - set initial locations
  !  - set initial conditions for other states (e.g. length)
  ! added by zhuxm
  !----------------------------------------------------------
  do n=1,ngroups
       call init_bio(igroups(n),igroups(n)%Tnind)
  end do


  !----------------------------------------------------------
  ! define variables for output
  !----------------------------------------------------------
  call cdf_out(ngroups,igroups,0,t,NCDO_ADD_STATES)

  ! Determine whether structured or unstructured model forcing.
  ! Compute individual-based model with either one of them
  !
  ! zfeng changed 09/25/2014
  if(Structured) then ! structured forcing (BIOMAS)
    call structured_compute
  elseif(.NOT. Structured) then ! unstructured forcing (FVCOM)
    call unstructured_compute
  else     ! error
    write(*,*) 'ERROR! Must setup with either structured or unstructured forcing'
    stop  ! end program
  end if ! structured

  !===================================================================
  ! <=  end main loop over time
  !===================================================================

  !---------------------------------------------------------
  ! cleanup data
  !---------------------------------------------------------
  call drawline("-") ; write(*,*)'Finalizing Simulation' ; call drawline("-")

End Program ArcIBM
