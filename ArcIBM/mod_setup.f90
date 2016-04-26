module mod_setup

!----------------------------------------------------
! move from mod_driver.f90 to this new module
! by zfeng 09/25/2014
!
! setup the model with either:
! 1. FVCOM unstructured forcing
! or 2. BIOMAS structured forcing
!----------------------------------------------------

contains

Subroutine setup
  use gparms
  use utilities
  use mod_output
  use mod_forcing

  implicit none
  integer             :: n,ios, missfile
  real(sp)            :: uno = 1.0
  real(sp)            :: t   = 0.0
  character(len=fstr) :: buffer
  logical             :: ffexist ! True if forcing_file exists

  !--------------------------------------------------------
  ! get primary control file name from the command line
  !--------------------------------------------------------
  n = iargc()
  if(n /= 1)then
    write(*,*)'incorrect usage: fiscm file.nml'
    write(*,*)'where file.nml is the main parameter file'
    stop 'stopping'
  endif
  call getarg(1,buffer)
  read(buffer,*) runcontrol

  !--------------------------------------------------------
  ! check for existence of primary control file
  !--------------------------------------------------------
  call drawline("-")
  write(*,*)'reading main parameters from: ',trim(runcontrol)
  inquire(file=trim(runcontrol),exist=fexist)
  if( .not. fexist)then
    write(*,*)'fatal error: namelist file: ',trim(runcontrol),' does not exist, stopping...'
    stop
  endif
  !--------------------------------------------------------
  ! open and read global namelist:  nml_fiscm
  !--------------------------------------------------------
  open(unit=iunit,file=trim(runcontrol),form='formatted')
  read(unit=iunit,nml=nml_fiscm,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read namelist fiscm from ',trim(runcontrol)
    stop
  endif
  !--------------------------------------------------------
  !sanity on number of groups
  !--------------------------------------------------------
  if(ngroups < 1) stop 'Fatal error: number of groups < 1'
  !--------------------------------------------------------
  !set begin/end time in seconds
  !--------------------------------------------------------
  beg_time = beg_time_days*day_2_sec
  end_time = end_time_days*day_2_sec
  !--------------------------------------------------------
  !set simulation direction based on order of begin/end times
  !--------------------------------------------------------
  if(end_time /= beg_time)then
    sim_direction = sign(uno, end_time-beg_time)
  else
    stop 'fatal error: begin and end time are identical'
  endif
  !--------------------------------------------------------
  !revised by Xinyou Lin to read extend file
  ! zhuxm changed extfile and extf_in to n_extfile
  ! extfile for information other than forcing
  !--------------------------------------------------------
  if(n_extfile /= 0) call open_ext_file(extfile_name,n_extfile)

  !--------------------------------------------------------
  ! Check for existence and open forcing file (if needed)
  ! --use keyword 'structured' to determine whether structured
  !   or unstructured forcing is used.
  !   structured = 0; unstructured = 1.
  ! zfeng changed on 09/26/2014
  !--------------------------------------------------------
  if( (.NOT. structured) .AND. nfiles_in > 0) then ! open unstructured forcing
    call open_unstructured_forcing(forcing_file,nfiles_in,fbeg,fend)
  elseif( structured .AND. nfiles_in > 0) then ! open structured forcing
  !--------------------------------------------------------
  ! Check existence of BIOMAS netcdf forcing
  ! If so, directly open netcdf forcing file.
  ! If not, generate netcdf file using 'biomas_nc' executable
  ! zfeng 10/10/2014
  !-------------------------------------------------------

    missfile = 0 ! flag for number of missing forcing files

    do n = 1, nfiles_in
       inquire(FILE=forcing_file(n),EXIST=ffexist)
       if (.NOT. ffexist ) THEN
          missfile = missfile + 1 ! count total number of missing files
          call drawline("-")
          write(*,*) 'Forcing file: '//trim(forcing_file(n))// ' is NOT in directory!'
       endif
    enddo ! n

    call drawline("-")
    write(*,*) 'Number of missing files: ', missfile

    if (missfile == 0) then !
       call drawline("-")
       write(*,*) 'All forcing files exist'
       call open_structured_forcing(forcing_file,nfiles_in,fbeg,fend)
    else
      call drawline("-")
      write(*,*) "Missing netCDF forcing files!"
      write(*,*) "Use 'biomas_nc' executable to generate netCDF files first!"
      write(*,*) "ArcIBM STOPS!"
      call drawline("=")
      stop

 ! open structured forcing files, returned values of fbeg & fend
      call open_structured_forcing(forcing_file,nfiles_in,fbeg,fend)
    endif ! missifile

  elseif(nfiles_in <= 0) then
    call drawline("-")
    write(*,*) 'NO forcing file for simulation!!!'
  else
    write(*,*) 'STOP: Forcing files must be either structured or unstructured!'
    stop
  endif ! nfiles_in


  !report
  write(*,*)'begin time       :',beg_time
  write(*,*)'end time         :',end_time
  if(sim_direction ==  1)then
     write(*,*)'direction        :    forward'
  elseif(sim_direction == -1)then
     write(*,*)'direction        :    backward'
  endif
  write(*,*)'time step(s)     :',deltaT
  write(*,*)'num groups       :',ngroups
  write(*,*)'number of report :',ireport
  call drawline("-")
  !---------------------------------------------------
  !read and allocate individual groups
  !---------------------------------------------------
  allocate(igroups(ngroups))

  do n=1,ngroups
    igroups(n) = group_(iunit,n,deltaT)

  end do
  !---------------------------------------------------
  !add state variables
  !---------------------------------------------------
  do n=1,ngroups
    call group_addstates(igroups(n))
  end do

  !----------------------------------------------------------
  ! summarize all states to the screen
  !----------------------------------------------------------
  do n=1,ngroups
    call print_group_summary(igroups(n))
  end do

  !open output files and assign netcdf ids to each group
  call cdf_out(ngroups,igroups,0,t,NCDO_HEADER)
  call drawline("-")
  write(*,*)'setup complete'
  call drawline("-")

End Subroutine setup

!---------------------------------------------------------------

end module mod_setup
