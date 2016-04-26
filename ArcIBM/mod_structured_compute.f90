 Module mod_structured_compute

!=======================================================================
! BISCM BIOMAS Individual-Based Model (alfa version, in progress)
!
! BIOMAS: Pan-Arctic Biology/Ice/Ocean Modeling & Assimilation System
!         developed by Dr. Jinlun Zhang group (UW-APL)
!
! Detailed descriptions of BIOMAS refer to:
! - Zhang et al, 2010. Modeling the impact of declining sea ice on the
!   Arctic marine planktonic ecosystem. Journal of Geophysical Research,
!   115(C10), C10015. doi:10.1029/2009JC005387
!
! Description
!    - Read and setup structured grid
!    - Advect
!    - Diffuse
!    - Locate Particles
!
! Comments:
!    - This module is for STRUCTURED model only!
!    - Present code is tailored for BIOMAS forcing, but can be recoded to
!      model behavior of individuals, driven by any other models, eg. ROMS.
!    - BIOMAS is on Arakawa B-grid. Refer to Parallel Ocean Program
!      (POP) Reference Manual Chapter 3 Spatial Discretization
!    - Advection currently uses 1st Order Euler Step, update to a stiffer
!      solver (for example 4-stage RK used by C. Chen)
!    - Vertical diffusion uses Vissers modified random walk
!    - Vertical diffusion improvement possibilities:
!         a.) splines (previously used by J. Pringle, R. Ji, and M. Huret)
!         b.) tensioned splines to avoid needing to smooth diffusivity
!             before splining
!         c.) binned random walk (Thygesen and Adlandsvik, MEPS 347,2007)
!
! !REVISION HISTORY:
!
! 2014/09/03       Zhixuan Feng  (zfeng@whoi.edu)
!
!=======================================================================

 contains

  subroutine structured_compute

    use gparms
    use mod_igroup
    use mod_driver
    use mod_bio
    use mod_output
    use mod_forcing
    use biomas_driver
    use utilities
    implicit none

  integer             :: i,n,its,nits
  real(sp)            :: t = 0.0
  real(sp)            :: tsmin,tsmax
  integer             :: nvars ,uniqvars
  character(len=fstr) :: needvars(max_state_vars)

  call drawline("-")
  write(*,*) 'Lagrangian Tracking and i-state Development Using BIOMAS Forcing!'
  !----------------------------------------------------------
  ! setup connection with BIOMAS forcing
  !   determine which forcing vars needed from states
  !   determine which forcing vars needed for adv/diffusion
  !   read in the mesh from netcdf forcing file [ocean model dependent]
  !   setup frames to store data read from netcdf
  !----------------------------------------------------------
  nvars = 0
  call get_ext_varnames(ngroups,igroups,nvars,needvars)
  call structured_model_init(ngroups,igroups,nvars,needvars) ! module in 'biomas_driver.f90'
  call get_unique_strings(nvars,needvars,uniqvars);nvars = uniqvars
  if(nvars > 0 .and. maxval(igroups%space_dim) > 1) call setup_forcing(nvars,needvars)

  !----------------------------------------------------------
  ! NEW subroutine 'initial_indices' in module 'biomas_driver'
  ! --'initial_indices' updates indices (i,j,k) based on particle
  !    locations and add the indices to 'igroups' variable.
  !----------------------------------------------------------
  call initial_indices(ngroups,igroups)

  if(.not. checkstatus(ngroups,igroups,t,its))  then
     write(*,*) 'All dead!'
     stop
  endif
  !----------------------------------------------------------
  ! ensure times (spawning, simulation, forcing)
  ! are compatible
  !----------------------------------------------------------
   tsmin=beg_time
   tsmax=end_time
  !  call get_spawn_range(ngroups,igroups,tsmin,tsmax)
  call check_times(beg_time,end_time,tsmin,tsmax,fbeg,fend)

  !----------------------------------------------------------
  ! => begin main loop over time
  !----------------------------------------------------------
  call drawline("-")
  write(*,*)'Beginning Simulation'
  call drawline("-")

  t = beg_time
  gt=t
  nits = (end_time-beg_time)/deltaT+1


!===================================================================
! => set for beginning main loop over time
!===================================================================
  if(maxval(igroups%space_dim) > 1) then
    call update_forcing(t,3)
    call activate(ngroups,igroups,t,sim_direction)
    call interp_structured_forcing(ngroups,igroups,3) ! interp forcing variables (temp) for biology calculation
    call cdf_out(ngroups,igroups,its,t,NCDO_OUTPUT)
    call do_bio(ngroups,igroups,t,its)
  endif
!===================================================================
! => begin main loop over time
!===================================================================

  do its=1,nits
      t = t + deltaT
      gt= t; 
      
    !---------------------------------------------------------
    ! update forcing to time t
    !---------------------------------------------------------
    if(maxval(igroups%space_dim) > 1) call update_forcing(t,4); 

    !---------------------------------------------------------
    ! check for activation through spawning
    !---------------------------------------------------------
    call activate(ngroups,igroups,t,sim_direction); 

    !---------------------------------------------------------
    ! advect and diffuse particles
    !---------------------------------------------------------
    call do_adv_diff(ngroups,igroups,deltaT,t);! write(*,*) 'Marker-9'

    !---------------------------------------------------------
    ! interpolate external vars at new particle positions
    !---------------------------------------------------------
    call interp_structured_forcing(ngroups,igroups,4); !4 is new current field

    !---------------------------------------------------------
    ! report progress to screen
    !---------------------------------------------------------
    if(.not. checkstatus(ngroups,igroups,t,its))exit
    !---------------------------------------------------------
    ! advance the biology in time
    !---------------------------------------------------------
    call do_bio(ngroups,igroups,t,its)
    !---------------------------------------------------------
    ! output the particle states to a NetCDF file
    !---------------------------------------------------------
    call cdf_out(ngroups,igroups,its,t,NCDO_OUTPUT)
    !---------------------------------------------------------
    ! update the time step
    !---------------------------------------------------------
    call exchange_forcing
  end do

  deallocate(igroups)

  end subroutine structured_compute

 End Module mod_structured_compute
