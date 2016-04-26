module mod_unstructured_compute

 !----------------------------------------------------------
 ! Note: zfeng created this new module on 09/25/2014
 !  and moved fiscm particle computation to this new module
 !----------------------------------------------------------

 contains

 subroutine unstructured_compute

  use gparms
  use mod_igroup
  use mod_driver
  use mod_bio
  use mod_output
  use mod_forcing
  use utilities

  implicit none

  integer             :: i,n,its,nits
  real(sp)            :: t = 0.0
  real(sp)            :: tsmin,tsmax
  integer             :: nvars ,uniqvars
  character(len=fstr) :: needvars(max_state_vars)

  !----------------------------------------------------------
  ! setup connection with forcing
  !   determine which forcing vars needed from states
  !   determine which forcing vars needed for adv/diffusion
  !   read in the grid (determine if it is )
  !   read in the mesh [ocean model dependent]
  !   setup frames to store data read from netcdf
  !----------------------------------------------------------
  nvars = 0
  call get_ext_varnames(ngroups,igroups,nvars,needvars)
  call ocean_model_init(ngroups,igroups,nvars,needvars)
  call get_unique_strings(nvars,needvars,uniqvars);nvars = uniqvars

  if(nvars > 0 .and. maxval(igroups%space_dim) > 1) call setup_forcing(nvars,needvars)

  !----------------------------------------------------------
  ! find elements containing particles
  !----------------------------------------------------------
  call update_element(ngroups,igroups)

  if(.not. checkstatus(ngroups,igroups,t,its))  stop 'all dead'

  !----------------------------------------------------------
  ! summarize all states to the screen
  !----------------------------------------------------------
  do n=1,ngroups
    call print_group_summary(igroups(n))
  end do

  !----------------------------------------------------------
  ! ensure times (spawning, simulation, forcing)
  ! are compatible
  !----------------------------------------------------------
   tsmin=beg_time
   tsmax=end_time
  !call get_spawn_range(ngroups,igroups,tsmin,tsmax)
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
    call sz_ini(ngroups,igroups)
    call activate(ngroups,igroups,t,sim_direction)
    call interp_forcing(ngroups,igroups,3)
    call cdf_out(ngroups,igroups,its,t,NCDO_OUTPUT)
    call do_bio(ngroups,igroups,t,its)
  endif
!===================================================================
! => begin main loop over time
!===================================================================

  do its=1,nits
      t = t + deltaT
      gt= t
    !---------------------------------------------------------
    ! update forcing to time t
    !---------------------------------------------------------
    if(maxval(igroups%space_dim) > 1) call update_forcing(t,4)

    !---------------------------------------------------------
    ! check for activation through spawning
    !---------------------------------------------------------
    call activate(ngroups,igroups,t,sim_direction)
    !---------------------------------------------------------

    ! advect and diffuse particles
    !---------------------------------------------------------
    call do_adv_diff(ngroups,igroups,deltaT,t)
    !---------------------------------------------------------
    ! interpolate external vars at new particle positions

    !---------------------------------------------------------
    call interp_forcing(ngroups,igroups,4) !4 is new current field
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

 end subroutine unstructured_compute

end module mod_unstructured_compute
