!=======================================================================
! Fiscm NetCDF Forcing Routines 
!
! Description
!   Setup and access NetCDF forcing 
!    - use frames (time slices) built from containers
!    - current frame interpolated from bracketing frames
!    
! Comments: 
!    - Modify to use D. Stuebe's NetCDF package once binned random-walk
!      and NURBS field reconstruction tested and working
!    - Requires FORTRAN90 NetCDF 3x libraries
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!  2011/7/28       Xueming Zhu
!  09/26/2014      zfeng add new subroutine 'open_structured_forcing'
!  11/05/2014      zfeng add new subroutine 'get_forcing_f3'
!=======================================================================
module mod_forcing
    use gparms
    use mod_igroup
    use utilities,only : ncdchk
    implicit none

    logical               :: debug = .false.
    integer, private      :: ffile_id,nDim,nVar,nAtt,uDid,fNum,uD_extant
    character(len=fstr), private :: ffile(max_nf)

    real(sp)              :: forcing_beg_time,forcing_end_time,forcing_deltaT
    integer               :: nframes
    real(sp), allocatable :: ftimes(:),tmptime(:)
    integer , allocatable :: fid_ft(:)
    integer , allocatable :: rec_ft(:)
    integer , allocatable :: fnr_ft(:)

    integer               :: now = 3 !frame 3 is current data
    integer               :: nex = 4 !frame 4 is next data

!----------------------------------------------------------------
! container type
!----------------------------------------------------------------
    type container
      integer :: dtype
      integer :: ttype
      integer :: ndims
      integer, pointer  :: dims(:)
      integer :: varid
      character(len=fstr) :: vname
      !can do in with single 1-d array if we could reshape allocatable
      real(sp), pointer :: f1(:)
      real(sp), pointer :: f2(:,:)
      real(sp), pointer :: f3(:,:,:)
      real(sp), pointer :: f4(:,:,:,:)
      integer , pointer :: i1(:)
      integer , pointer :: i2(:,:)
      integer , pointer :: i3(:,:,:)
      integer , pointer :: i4(:,:,:,:)
    end type container

!----------------------------------------------------------------
! data frame type
!----------------------------------------------------------------
    type dataframe
      real(sp) :: time
      real(sp) :: ncfid
      real(sp) :: iframe
      integer  :: nvars
      integer  :: id
      logical  :: initial_read
      type(container), pointer :: fdata(:)   ![nvars]
    end type dataframe

    type(dataframe), save       :: frame(4)

!===========
    real(sp), allocatable :: extftime(:)
    real(sp), allocatable :: ext_data1(:,:),ext_data2(:,:)
    integer , allocatable :: fid_exft(:),rec_exft(:),fnr_exft(:)
    integer               :: nframex,extf1,extf2
!===========

!overload ">" to be able to compare the order of time frames in time
    interface operator(>)
      module procedure frame_order
    end interface

    interface get_forcing
      module procedure get_forcing_f1
      module procedure get_forcing_f2
    end interface

contains


!========================================================================
! open unstructured forcing file
!   - make sure it exists
!   - read time (convert to sec if necessary.)
!   - zfeng changed subroutine name 'open_forcing_file' to 'open_unstructured_forcing'
!========================================================================
subroutine open_unstructured_forcing(ffiles_in,nfls,fbeg,fend)
  use utilities

  implicit none
  integer, intent(in)  :: nfls
  character(len=fstr)  :: ffiles_in(nfls),ffile_in
  integer              :: nfsfid(nfls),nfsnrec(nfls)
  real(sp),intent(out) :: fbeg,fend
  integer              :: fid,varid,i,ierr,n
  character(len=mstr)  :: msg
  character(len=fstr)  :: dname,tunits
!Added on Jun 11 ,2009
  integer              :: nb,ne,tmpii,y
  real(sp)             :: tmpxx

!==============
  nframes=0
  do n=1,nfls
     ffile_in=ffiles_in(n) !maybe it goes wrong
     msg = "error opening unstructured forcing file: "//trim(ffile_in)
     call ncdchk( nf90_open(trim(ffile_in),nf90_nowrite,fid),msg)

     nfsfid(n)= fid
     ffile(n) = ffile_in

     msg = "reading number of dimensions from: "//trim(ffile_in)
     call ncdchk(nf90_inquire(fid, nDim,nVar,nAtt,uDid,fNum) ,msg)

 ! Following lines are modified in case NetCDF time is not unlimited
    write(*,*) 'uDid=',uDid

    if(uDid == -1) then ! CANNOT find unlimited dimension (time)
      call ncdchk(nf90_inq_dimid(fid,'time',uDid))
      dname = 'time'
    end if

     uD_extant = 1
     if(uDid /= -1)then
       call ncdchk(nf90_inquire_dimension(fid, uDid, dname, uD_extant ))
     endif
     nframes = uD_extant+nframes
     nfsnrec(n) = uD_extant
  end do

!=====================================
!Revised on Jun 11, 2009
!add year_cycle  by zhuxm
!=====================================
  !allocate space to hold times and read in
  if(year_cycle==0)then
     allocate(ftimes(nframes)) ; ftimes = 0.0_sp
     allocate(fid_ft(nframes)) ; fid_ft = 0
     allocate(rec_ft(nframes)) ; rec_ft = 0
     allocate(fnr_ft(nframes)) ; fnr_ft = 0
  else
     allocate(ftimes(nframes*year_cycle+1)) ; ftimes = 0.0_sp   !two years
     allocate(fid_ft(nframes*year_cycle+1)) ; fid_ft = 0
     allocate(rec_ft(nframes*year_cycle+1)) ; rec_ft = 0
     allocate(fnr_ft(nframes*year_cycle+1)) ; fnr_ft = 0
  endif
!=====================================
  msg = "error reading time variable from netcdf file"
  uD_extant=0
  do n=1,nfls
     uD_extant=uD_extant+nfsnrec(n)
     ffile_id=nfsfid(n)
     call ncdchk(nf90_inq_varid(ffile_id,'time',varid),msg )

     allocate(tmptime(nfsnrec(n))) ; tmptime = 0.0
     call ncdchk(nf90_get_var(ffile_id, varid,tmptime ),msg)

     if(nf90_get_att(ffile_id, varid, 'units', tunits) == nf90_noerr)then
        if(index(tunits,'day') /= 0) tmptime = tmptime*day_2_sec
     endif
     tmptime=tmptime-mjd_offset*day_2_sec

     ftimes(uD_extant-nfsnrec(n)+1 : uD_extant) =tmptime(1:nfsnrec(n))
     fid_ft(uD_extant-nfsnrec(n)+1 : uD_extant) = ffile_id
     rec_ft(uD_extant-nfsnrec(n)+1 : uD_extant) = (/(i,i=1,nfsnrec(n))/)
     fnr_ft(uD_extant-nfsnrec(n)+1 : uD_extant) = n
     deallocate(tmptime)
  enddo

!Revised on Jun 11, 2009
!To make new time variable
  nb=0
  ne=0
  if(nframes > 1)then
     do i=1,nframes
        if(ftimes(i)>= 0.0_sp .and. nb == 0)nb=i
        if(ftimes(i)<= 365*24*3600.0_sp )ne=i
     end do
  endif
!check the recycle
  if(year_cycle/=0)then
     if(ftimes(ne)-ftimes(nb) == 365*24*3600.0_sp)ne=ne-1
  endif
!========
  n=0
  do i=nb,ne
     n=n+1
     tmpxx=ftimes(i)
     ftimes(n)= tmpxx

     tmpii=fid_ft(i)
     fid_ft(n)= tmpii

     tmpii=rec_ft(i)
     rec_ft(n)= tmpii

     tmpii=fnr_ft(i)
     fnr_ft(n)= tmpii
  enddo

  if(year_cycle/=0)then
    do y=2,year_cycle
      do i=(y-1)*n+1,y*n
        ftimes(i)= ftimes(i-n)+365*24*3600.0_sp
        rec_ft(i)= rec_ft(i-n)
        fid_ft(i)= fid_ft(i-n)
        fnr_ft(i)= fnr_ft(i-n)
      enddo
!make it closed
      ftimes(y*n+1)= ftimes(1)+y*365*24*3600.0_sp
      rec_ft(y*n+1)= rec_ft(1)
      fid_ft(y*n+1)= fid_ft(1)
      fnr_ft(y*n+1)= fnr_ft(1)
    enddo
  endif

  if(nframes > 1)then
    do i=2,nframes
      if(ftimes(i)-ftimes(i-1) <= 0.0)then
         write(*,*)'netcdf time is not monotonically increasing'
         write(*,*)i,ftimes(i),ftimes(i-1)
      endif
    end do
  endif

  !set begin/end/deltaT from forcing
  !make sure time is sequential
  forcing_beg_time = ftimes(1)
  if(year_cycle/=0)then
    forcing_end_time = ftimes(year_cycle*n+1)
    nframes=year_cycle*n+1
  else
    forcing_end_time = ftimes(nframes)
  endif

  !set fbeg/fend to return val
  fbeg = forcing_beg_time
  fend = forcing_end_time

  call drawline('-')
  write(*,*)'Opened up unstructured model forcing file: '
  do n=1,nfls
     write(*,*)trim(ffiles_in(n))
  enddo
  write(*,*)'number of frames in file: ',nfsnrec
  write(*,*)'forcing begins at:',gettime(int(fbeg))
  write(*,*)'forcing ends   at:',gettime(int(fend))
  call drawline('-')

end subroutine open_unstructured_forcing

!========================================================================
! Open structured forcing file
!  - zfeng added 09/26/2014
!  - zfeng changed 10/10/2014
!  - Most lines were copied from 'open_unstructured_forcing'
!========================================================================
subroutine open_structured_forcing(ffiles_in,nfls,fbeg,fend)
  use utilities

  implicit none
  integer, intent(in)  :: nfls
  character(len=fstr)  :: ffiles_in(nfls),ffile_in
  integer              :: nfsfid(nfls),nfsnrec(nfls)
  real(sp),intent(out) :: fbeg,fend
  integer              :: fid,varid,i,ierr,n
  character(len=mstr)  :: msg
  character(len=fstr)  :: dname,tunits
!Added on Jun 11 ,2009
  integer              :: nb,ne,tmpii,y
  real(sp)             :: tmpxx



!==============
  nframes=0
  do n=1,nfls
     ffile_in=ffiles_in(n) !maybe it goes wrong
     msg = "error opening structured forcing file: "//trim(ffile_in)
     call ncdchk( nf90_open(trim(ffile_in),nf90_nowrite,fid),msg)

     nfsfid(n)= fid
     ffile(n) = ffile_in

     msg = "reading number of dimensions from: "//trim(ffile_in)
     call ncdchk(nf90_inquire(fid, nDim,nVar,nAtt,uDid,fNum) ,msg)

     uD_extant = 1
     if(uDid /= 0)then
       call ncdchk(nf90_inquire_dimension(fid, uDid, dname, uD_extant ))
     endif
     nframes = uD_extant+nframes
     nfsnrec(n) = uD_extant
  end do

!=====================================
!Revised on Jun 11, 2009
!add year_cycle  by zhuxm
!=====================================
  !allocate space to hold times and read in
  if(year_cycle==0)then
     allocate(ftimes(nframes)) ; ftimes = 0.0_sp
     allocate(fid_ft(nframes)) ; fid_ft = 0
     allocate(rec_ft(nframes)) ; rec_ft = 0
     allocate(fnr_ft(nframes)) ; fnr_ft = 0
  else
     allocate(ftimes(nframes*year_cycle+1)) ; ftimes = 0.0_sp   !two years
     allocate(fid_ft(nframes*year_cycle+1)) ; fid_ft = 0
     allocate(rec_ft(nframes*year_cycle+1)) ; rec_ft = 0
     allocate(fnr_ft(nframes*year_cycle+1)) ; fnr_ft = 0
  endif
!=====================================
  msg = "error reading time variable from netcdf file"
  uD_extant=0
  do n=1,nfls
     uD_extant=uD_extant+nfsnrec(n)
     ffile_id=nfsfid(n)
     call ncdchk(nf90_inq_varid(ffile_id,'time',varid),msg )

     allocate(tmptime(nfsnrec(n))) ; tmptime = 0.0
     call ncdchk(nf90_get_var(ffile_id, varid,tmptime ),msg)

     if(nf90_get_att(ffile_id, varid, 'units', tunits) == nf90_noerr)then
        if(index(tunits,'day') /= 0) tmptime = tmptime*day_2_sec
     endif
     tmptime=tmptime-mjd_offset*day_2_sec

     ftimes(uD_extant-nfsnrec(n)+1 : uD_extant) =tmptime(1:nfsnrec(n))
     fid_ft(uD_extant-nfsnrec(n)+1 : uD_extant) = ffile_id
     rec_ft(uD_extant-nfsnrec(n)+1 : uD_extant) = (/(i,i=1,nfsnrec(n))/)
     fnr_ft(uD_extant-nfsnrec(n)+1 : uD_extant) = n
     deallocate(tmptime)
  enddo

!Revised on Jun 11, 2009
!To make new time variable
  nb=0
  ne=0
  if(nframes > 1)then
     do i=1,nframes
        if(ftimes(i)>= 0.0_sp .and. nb == 0)nb=i
        if(ftimes(i)<= 365*24*3600.0_sp )ne=i
     end do
  endif
!check the recycle
  if(year_cycle/=0)then
     if(ftimes(ne)-ftimes(nb) == 365*24*3600.0_sp)ne=ne-1
  endif
!========
  n=0
  do i=nb,ne
     n=n+1
     tmpxx=ftimes(i)
     ftimes(n)= tmpxx

     tmpii=fid_ft(i)
     fid_ft(n)= tmpii

     tmpii=rec_ft(i)
     rec_ft(n)= tmpii

     tmpii=fnr_ft(i)
     fnr_ft(n)= tmpii
  enddo

  if(year_cycle/=0)then
    do y=2,year_cycle
      do i=(y-1)*n+1,y*n
        ftimes(i)= ftimes(i-n)+365*24*3600.0_sp
        rec_ft(i)= rec_ft(i-n)
        fid_ft(i)= fid_ft(i-n)
        fnr_ft(i)= fnr_ft(i-n)
      enddo
!make it closed
      ftimes(y*n+1)= ftimes(1)+y*365*24*3600.0_sp
      rec_ft(y*n+1)= rec_ft(1)
      fid_ft(y*n+1)= fid_ft(1)
      fnr_ft(y*n+1)= fnr_ft(1)
    enddo
  endif

  if(nframes > 1)then
    do i=2,nframes
      if(ftimes(i)-ftimes(i-1) <= 0.0)then
         write(*,*)'netcdf time is not monotonically increasing'
         write(*,*)i,ftimes(i),ftimes(i-1)
      endif
    end do
  endif

  !set begin/end/deltaT from forcing
  !make sure time is sequential
  forcing_beg_time = ftimes(1)
  if(year_cycle/=0)then
    forcing_end_time = ftimes(year_cycle*n+1)
    nframes=year_cycle*n+1
  else
    forcing_end_time = ftimes(nframes)
  endif

  !set fbeg/fend to return val
  fbeg = forcing_beg_time
  fend = forcing_end_time

  call drawline('-')
  write(*,*)'Opened up structured model forcing file: '
  do n=1,nfls
     write(*,*)trim(ffiles_in(n))
  enddo
  write(*,*)'number of frames in file: ',nfsnrec
  write(*,*)'forcing begins at:',gettime(int(fbeg))
  write(*,*)'forcing ends   at:',gettime(int(fend))
  call drawline('-')

end subroutine open_structured_forcing


!========================================================================
! update the model forcing to time (t)
!   using linear interpolation
!========================================================================
subroutine update_forcing(t,iframe0)
  implicit none
  real(sp), intent(in) :: t
  integer :: i1,i2,i1f,i2f,iframe0
  logical inew

  call bracket(t,i1,i2) 
  i1f = 0 ; i2f = 0


  !determine if either i1 or i2 are already loaded 
  if(frame(1)%iframe == i1) i1f = 1  
  if(frame(2)%iframe == i1) i1f = 2  
  if(frame(1)%iframe == i2) i2f = 1  
  if(frame(2)%iframe == i2) i2f = 2  


  !load i1 if necessary 
  if(i1f == 0)then
    if(i2f == 1) then
      call read_frame(frame(2),i1)
      i1f = 2
    else
      call read_frame(frame(1),i1)
      i1f = 1
    endif
  endif


  !load i2 if necessary 
  if(i2f == 0)then
    if(i1f == 1) then
      call read_frame(frame(2),i2)
    else
      call read_frame(frame(1),i2)
    endif
  endif

  call interp_two_frames(frame(iframe0),t,frame(1),frame(2))

end subroutine update_forcing

!========================================================================
! setup frames (set variable names, allocate, etc)
!   we need 3 frames, one to store two time frames
!   from the netcdf files, and one to store data
!   interpolated to time (t)
!========================================================================
subroutine setup_forcing(nvars,varlist)
  implicit none
  integer, intent(in) :: nvars
  character(len=*)    :: varlist(nvars)
  integer             :: i

  do i=1,4
     frame(i) = frame_(nvars,varlist,i)
  enddo

  call frame_info(frame(1),FRAME_SETUP)

end subroutine setup_forcing 

!========================================================================
! dump frame info to screen 
!========================================================================
subroutine frame_info(frame,istatus)
  use utilities, only : drawline
  implicit none
  type(dataframe), intent(in) :: frame
  integer                     :: istatus
  !----------------------------
  integer             :: i,j,ndims,dims(4)
  real(sp)            :: fmax,fmin,fave 
  character(len=fstr) :: vname
  logical             :: tv
 
  if(frame%nvars < 1)return

  !------------------------------------------------------------
  ! dump frame vars and basic information (dimensions and type)
  !------------------------------------------------------------
  if(istatus == FRAME_SETUP)then
    call drawline("-")
    write(*,*)' forcing var   |   | dims |   dim1  dim2  dim3  dim4' 
    call drawline("-")
    do i=1,frame%nvars 
      tv = .false.
      if(frame%fdata(i)%ttype == TIME_VARYING) tv = .true.
      ndims = frame%fdata(i)%ndims
      vname = frame%fdata(i)%vname
      dims  = 0
      dims(1:ndims)  = frame%fdata(i)%dims(1:ndims)
      write(*,'(A15,A3,L1,A3,I4,A3,4I6)')vname,' | ',tv,' | ',ndims,' | ',(dims(j),j=1,ndims)
    end do
  !------------------------------------------------------------
  ! dump frame vars and stats on variables 
  !------------------------------------------------------------
  elseif(istatus == FRAME_STATS) then
    call drawline("-")
    write(*,*)' forcing var   |    min |    max |    ave'  
    call drawline("-")
    do i=1,frame%nvars
      fmax = calc_stat(frame,frame%fdata(i)%vname,'max')
      fmin = calc_stat(frame,frame%fdata(i)%vname,'min')
      fave = calc_stat(frame,frame%fdata(i)%vname,'mean')
      vname = frame%fdata(i)%vname
      write(*,'(A15,A3,F6.2,A3,F6.2,A3,F6.2)')vname,' | ',fmin,' | ',fmax,' | ',fave
    end do
  else
    write(*,*)'fatal error in frame_info'
    write(*,*)'istatus must be: ',FRAME_SETUP,' or ',FRAME_STATS
    stop
  endif
  call drawline("-")
  
end subroutine frame_info

!========================================================================
! frame type constructor 
!========================================================================
function frame_(nvars,varlist,id) result(frame)

  implicit none
  integer, intent(in) :: nvars
  character(len=fstr) :: varlist(nvars)
  integer, intent(in) :: id
  type(dataframe)     :: frame
  !--------------------------------------
  integer             :: dim1,dim2,i,i1,i2,i3,j,jj
  integer             :: fid,nAtts,vid,ndims,xtype,ttype,dtype
  character(len=mstr) :: msg,msg2
  character(len=fstr) :: vname,dname
  integer             :: dimids(NF90_MAX_VAR_DIMS)
  integer             :: dims(NF90_MAX_VAR_DIMS)
  integer             :: alldims(NF90_MAX_VAR_DIMS)
  integer             :: varids(nvars)

  !set initial values for time and frame_num
  frame%time   = -hugenum
  frame%ncfid = ffile_id 
  frame%iframe = -1
  frame%nvars  = nvars
  frame%id     = id
  ffile_id=fid_ft(1)
  
  !return if nvars = 0
  if(nvars < 1)return
  
  !make sure all the variables exist
  do i=1,nvars
    if(debug) write(*,*)'making sure: ',trim(varlist(i)),' exists'
    msg  = "error: var ["//trim(varlist(i))//"]"
    msg2 = " not found in "//trim(ffile(1))
    call ncdchk( nf90_inq_varid(ffile_id,varlist(i),varids(i)),msg,msg2 ) 
  end do

  !allocate primary dataspace
  allocate(frame%fdata(nvars)) 

  !for each variable, set name, get dims, allocate
  do i=1,nvars
    vid = varids(i) 
    call ncdchk(nf90_inquire_variable(ffile_id,vid,vname,xtype,ndims,dimids,nAtts))
    alldims = 0
    dims = 0
    jj = 0
    !loop over dimensions
    do j=1,ndims
      jj = jj + 1
      call ncdchk(nf90_inquire_dimension(ffile_id, dimids(j), dname, alldims(j)) )
      if(dimids(j) == uDid)then
         ndims = ndims - 1
         ttype = TIME_VARYING
      else
         dims(jj) = alldims(j)
         ttype = CONSTANT
      end if
    end do

    !set vartype and check
    if(xtype == NF90_FLOAT .or. xtype == NF90_DOUBLE)then
      dtype = flt_type
    elseif(xtype == NF90_INT)then
      dtype = int_type
    else
      write(*,*)'fatal error in frame_: type: ',xtype,' not defined'
      stop
    endif

    !set dimensions and check
    if(ndims > 4)then
      write(*,*)'fatal error in frame_: number of var dimensions exceeds 4 '
      stop
    endif
  
    !set container
    frame%fdata(i) = container_(dtype,ttype,vid,vname,ndims,dims)  

    !set initial read to false
    frame%initial_read = .false.
  end do
end function frame_

!========================================================================
! read data into a frame from file frame i 
!========================================================================
subroutine read_frame(frame,f)
 
  implicit none
  type(dataframe),intent(inout) :: frame
  integer,intent(in)            :: f
  integer                       :: i,varid,ndims,ttype,stat,dtype
  integer, allocatable, dimension(:) :: start,vsize

  !make sure frame is in bounds
  if(f > nframes .or. f < 1)then
    write(*,*)'fatal error in read_frame'
    write(*,*)'cannot read frame: ',f,' from file: ',trim(ffile(fnr_ft(f)))
    write(*,*)'max frames in netcdf file is: ',nframes
    stop 
  endif

  frame%iframe = f
  frame%time   = ftimes(f)
  ffile_id     = fid_ft(f)    

  if(debug)write(*,*)'============data frame: ',frame%id,' reading frame ',f

  !loop over vars, read data into frame
  do i=1,frame%nvars
    if(frame%initial_read.and. (frame%fdata(i)%ttype == CONSTANT))cycle
    if(debug)write(*,*)'reading data into: ',trim(frame%fdata(i)%vname)

    !set dimensions, stride, start
    varid = frame%fdata(i)%varid
    ndims = frame%fdata(i)%ndims
    ttype = frame%fdata(i)%ttype
    dtype = frame%fdata(i)%dtype   !!zhuxm added

    !static var, read in entirety 
    if(ttype == CONSTANT .and. dtype == flt_type)then
!!zhuxm changed from if to select case
      select case(ndims)
        case(1)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f1) )
        case(2)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f2) )
        case(3)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f3) )
        case(4)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f4) )
        case default
          write(*,*)'variable dimensions incorrect, ndims=',ndims
      end select
    elseif(ttype == CONSTANT .and. dtype == int_type)then
      select case(ndims)
        case(1)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i1) )
        case(2)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i2) )
        case(3)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i3) )
        case(4)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i4) )
        case default
          write(*,*)'variable dimensions incorrect, ndims=',ndims
      end select
    elseif(ttype == TIME_VARYING .and. dtype == flt_type) then !time-varying variable
      allocate(start(ndims+ttype),vsize(ndims+ttype))
      start(1:ndims) = 1
      vsize(1:ndims) = frame%fdata(i)%dims(1:ndims)
      start(ndims+1) = rec_ft(f)
      vsize(ndims+1) = 1
      select case(ndims)
        case(1)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f1,start,vsize) )
        case(2)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f2,start,vsize) )
        case(3)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f3,start,vsize) )
        case(4)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f4,start,vsize) )
        case default
          write(*,*)'variable dimensions incorrect, ndims=',ndims
      end select
      deallocate(start,vsize)
    elseif(ttype == TIME_VARYING .and. dtype == int_type)then
      allocate(start(ndims+ttype),vsize(ndims+ttype))
      start(1:ndims) = 1
      vsize(1:ndims) = frame%fdata(i)%dims(1:ndims)
      start(ndims+1) = rec_ft(f)
      vsize(ndims+1) = 1
      select case(ndims)
        case(1)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i1,start,vsize) )
        case(2)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i2,start,vsize) )
        case(3)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i3,start,vsize) )
        case(4)
          call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%i4,start,vsize) )
        case default
          write(*,*)'variable dimensions incorrect, ndims=',ndims
      end select
      deallocate(start,vsize)
    endif
  end do
  !set initial_frame to .true. will no longer read non-time-varying vars
  frame%initial_read = .true.

end subroutine read_frame
  
!========================================================================
! linearly interpolate between two data frames 
!========================================================================
subroutine interp_two_frames(frame,t,frame1,frame2)
  implicit none
  type(dataframe), intent(inout) :: frame
  real(sp), intent(in)           :: t
  type(dataframe), intent(in)    :: frame1,frame2

  !----------------------------------------
  real(sp) :: t1,t2,c1,c2
  integer  :: ndims,dtype,i
  integer  :: d1,d2,d3,d4



!! zhuxm add to check frame1 and frame2 are identical
  if(frame1%nvars /= frame2%nvars) then
     write(*,*)'Fatal error: frame1%nvars=',frame1%nvars
     write(*,*)'             frame2%nvars=',frame2%nvars
     write(*,*)'They should be identical.'
     stop
  endif
  !set interpolation coefficients
  t1 = frame1%time
  t2 = frame2%time
  if(t2 > t1)then
    c2 = (t - t1)/(t2-t1)
    c1 = 1.0_sp-c2
  else
    c1 = (t - t2)/(t1-t2)
    c2 = 1.0_sp-c1
  endif

  !$OMP PARALLEL PRIVATE(dtype, ndims)
  !$OMP DO
  !loop over vars, interpolate to intermediate frame 
  do i=1,frame1%nvars
    dtype = frame1%fdata(i)%dtype
    if(frame1%fdata(i)%dtype /= frame2%fdata(i)%dtype) then
      write(*,*)'Fatal error: frame1%fdata(i)%dtype=',frame1%fdata(i)%dtype
      write(*,*)'             frame2%fdata(i)%dtype=',frame1%fdata(i)%dtype
      write(*,*)'They should be identical.'
      stop
    endif

    ndims = frame1%fdata(i)%ndims
    if(frame1%fdata(i)%ndims /= frame2%fdata(i)%ndims) then
      write(*,*)'Fatal error: frame1%fdata(i)%ndims=',frame1%fdata(i)%ndims
      write(*,*)'             frame2%fdata(i)%ndims=',frame1%fdata(i)%ndims
      write(*,*)'They should be identical.'
      stop
    endif
    if(dtype == flt_type)then
      if(ndims == 1)     then
        frame%fdata(i)%f1 = c1*frame1%fdata(i)%f1 + c2*frame2%fdata(i)%f1
      elseif(ndims == 2) then
            frame%fdata(i)%f2 = c1*frame1%fdata(i)%f2 + c2*frame2%fdata(i)%f2
      elseif(ndims == 3) then
        frame%fdata(i)%f3 = c1*frame1%fdata(i)%f3 + c2*frame2%fdata(i)%f3
      elseif(ndims == 4) then
        frame%fdata(i)%f4 = c1*frame1%fdata(i)%f4 + c2*frame2%fdata(i)%f4
      endif

    elseif(dtype == int_type)then
      if(ndims == 1)     then
        frame%fdata(i)%i1 = c1*frame1%fdata(i)%i1 + c2*frame2%fdata(i)%i1
      elseif(ndims == 2) then
        frame%fdata(i)%i2 = c1*frame1%fdata(i)%i2 + c2*frame2%fdata(i)%i2
      elseif(ndims == 3) then
        frame%fdata(i)%i3 = c1*frame1%fdata(i)%i3 + c2*frame2%fdata(i)%i3
      elseif(ndims == 4) then
        frame%fdata(i)%i4 = c1*frame1%fdata(i)%i4 + c2*frame2%fdata(i)%i4
      endif
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine interp_two_frames

!========================================================================
! check frame order  
! returns true if frame2 is after frame1 
!========================================================================
function frame_order(frame1,frame2) result(order)
  implicit none
  type(dataframe), intent(in) :: frame1
  type(dataframe), intent(in) :: frame2
  logical :: order
  order = .false.
  if(frame2%time > frame1%time) order = .true. 
end function frame_order

!========================================================================
! container constructor
!   set name
!   set varid
!   set dimensions
!   allocate space
!========================================================================
function container_(dtype,ttype,vid,vname,ndims,dims) result(c)
  implicit none
  integer, intent(in) :: dtype 
  integer, intent(in) :: ttype 
  integer, intent(in) :: vid
  character(len=*)    :: vname
  integer, intent(in) :: ndims
  integer, intent(in) :: dims(ndims)
  type(container)     :: c 
  integer             :: adims(4)

  !check dimensions
  if(minval(dims) < 1)then 
    write(*,*)'fatal error in container_'
    write(*,*)'input dimensions provided < 0'
    write(*,*)dims
    stop
  endif

  !set and check the number of dimensions
  allocate(c%dims(ndims)) ; c%dims = dims 
  adims = 1
  adims(1:ndims) = dims

  c%dtype = dtype
  c%ttype = ttype
  c%vname = vname
  c%ndims = ndims
  c%varid = vid

  !allocate dataspace - float and initialize to zero
  if(dtype == flt_type)then
    if(ndims==1)then
       allocate(c%f1(adims(1)))  
       c%f1 = 0.0
    elseif(ndims==2)then
       allocate(c%f2(adims(1),adims(2)))
       c%f2 = 0.0
    elseif(ndims==3)then
       allocate(c%f3(adims(1),adims(2),adims(3)))
       c%f3 = 0.0
    elseif(ndims==4)then
       allocate(c%f4(adims(1),adims(2),adims(3),adims(4)))
       c%f4 = 0.0
    endif
  elseif(dtype == int_type)then
    if(ndims==1)then
       allocate(c%i1(adims(1)))  
       c%i1 = 0.0
    elseif(ndims==2)then
       allocate(c%i2(adims(1),adims(2)))
       c%i2 = 0.0
    elseif(ndims==3)then
       allocate(c%i3(adims(1),adims(2),adims(3)))
       c%i3 = 0.0
    elseif(ndims==4)then
       allocate(c%i4(adims(1),adims(2),adims(3),adims(4)))
       c%i4 = 0.0
    endif
  else
    write(*,*)'fatal error in container_'
    write(*,*)'data type: ',dtype,' not a valid type'
    stop
  endif
end function container_

!========================================================================
! calculate statistic /stat_type/ from variable /vname/ in frame /frame/ 
!  options are:  min,max,mean
!========================================================================
function calc_stat(frame,vname,stat_type) result(stat)
  implicit none
  type(dataframe), intent(in) :: frame
  character(len=*) :: vname
  character(len=*) :: stat_type 
  real(sp) :: stat
  integer  :: i,ndims,dimt,j
  if(index(vname,'_ext') /= 0)return 
  !determine the var
  i = valindx(frame,vname) 

  !get dimension
  ndims = frame%fdata(i)%ndims

  !-gwc need to expand to integer type
  if(stat_type == 'max')then
    if(ndims==1)stat = maxval(frame%fdata(i)%f1)
    if(ndims==2)stat = maxval(frame%fdata(i)%f2)
    if(ndims==3)stat = maxval(frame%fdata(i)%f3)
    if(ndims==4)stat = maxval(frame%fdata(i)%f4)
  elseif(stat_type == 'min')then
    if(ndims==1)stat = minval(frame%fdata(i)%f1)
    if(ndims==2)stat = minval(frame%fdata(i)%f2)
    if(ndims==3)stat = minval(frame%fdata(i)%f3)
    if(ndims==4)stat = minval(frame%fdata(i)%f4)
  elseif(stat_type == 'mean')then
    dimt = 1 
    do j=1,ndims
      dimt = dimt*frame%fdata(i)%dims(j)
    end do
    if(ndims==1)stat = sum(frame%fdata(i)%f1)/float(dimt)
    if(ndims==2)stat = sum(frame%fdata(i)%f2)/float(dimt)
    if(ndims==3)stat = sum(frame%fdata(i)%f3)/float(dimt)
    if(ndims==4)stat = sum(frame%fdata(i)%f4)/float(dimt)
  else
    write(*,*)'fatal error in calc_stat'
    write(*,*)'do not know how to compute statistic: ',trim(stat_type)
    stop
  endif

end function calc_stat

function valindx(frame,vname) result(indx)

  implicit none
  type(dataframe), intent(in) :: frame
  character(len=*)            :: vname
  integer                     :: indx,i

  if(index(vname,'_ext') /= 0)return
  indx = 0
  do i=1,frame%nvars
    if(frame%fdata(i)%vname == trim(vname)) indx = i
  end do

  if(indx == 0)then
    write(*,*)'error in valindx'
    write(*,*)'frame does not contain variable: ',trim(vname)
    stop
  endif   
end function valindx

!========================================================================
! interface data to outside world through pointer
!  - 1D, float
!========================================================================
subroutine get_forcing_f1(vname,p,iframe) 
  implicit none
  character(len=*) :: vname
  real(sp), pointer :: p(:)
  integer :: i,ndims,iframe
!!zhuxm  integer, allocatable :: dims(:)

  !get variable index
  i = valindx(frame(iframe),vname)

  !get dimensions
  ndims = frame(iframe)%fdata(i)%ndims
  if(ndims /= 1)then
    write(*,*)'error in get_forcing_f1'
    write(*,*)'number of dimensions of variable: ',trim(vname)
    write(*,*)'is: ',ndims,' cannot put this into a 1D array'
    stop
  endif

!!zhuxm  allocate(dims(ndims)) ; dims = 0
!!zhuxm  dims = frame(iframe)%fdata(i)%dims

  !set pointer
  p => frame(iframe)%fdata(i)%f1
end subroutine get_forcing_f1
!========================================================================
! interface data to outside world through pointer
!  - 2D, float
!========================================================================
subroutine get_forcing_f2(vname,p,iframe) 

  implicit none
  character(len=*)     :: vname
  real(sp), pointer    :: p(:,:)
  integer              :: i,ndims,iframe
!!zhuxm  integer, allocatable :: dims(:)

  !get variable index
  i = valindx(frame(iframe),vname)

  !get dimensions
  ndims = frame(iframe)%fdata(i)%ndims
  if(ndims /= 2)then
    write(*,*)'error in get_forcing_f2'
    write(*,*)'number of dimensions of variable: ',trim(vname)
    write(*,*)'is: ',ndims,' cannot put this into a 2D array'
    stop
  endif

!!zhuxm  allocate(dims(ndims)) ; dims = 0
!!zhuxm  dims = frame(iframe)%fdata(i)%dims

  !set pointer
  p => frame(iframe)%fdata(i)%f2
end subroutine get_forcing_f2

!========================================================================
! interface data to outside world through pointer
!  - 3D, float
! zfeng 11/05/2014
!========================================================================
subroutine get_forcing_f3(vname,p,iframe)

  implicit none
  character(len=*)     :: vname
  real(sp), pointer    :: p(:,:,:)
  integer              :: i,ndims,iframe
!!zhuxm  integer, allocatable :: dims(:)

  !get variable index
  i = valindx(frame(iframe),vname)

  !get dimensions
  ndims = frame(iframe)%fdata(i)%ndims
  if(ndims /= 3)then
    write(*,*)'error in get_forcing_f3'
    write(*,*)'number of dimensions of variable: ',trim(vname)
    write(*,*)'is: ',ndims,' cannot put this into a 3D array'
    stop
  endif

!!zhuxm  allocate(dims(ndims)) ; dims = 0
!!zhuxm  dims = frame(iframe)%fdata(i)%dims

  !set pointer
  p => frame(iframe)%fdata(i)%f3
end subroutine get_forcing_f3

!--------------------------------------------------------------------
function get_ncfid() result(fid)
  integer :: fid
  fid = ffile_id
end function get_ncfid

subroutine bracket(t,i1,i2) 
  use utilities

  real(sp), intent(in) :: t
  integer, intent(out) :: i1,i2
  !------------------------
  integer              :: i

  if(nframes <= 1)then
    write(*,*)'error in bracket: netcdf forcing file only has 1 frame'
    stop
  endif

  if(t < ftimes(1))then
    write(*,*)'error in bracket:'
    write(*,*)'model time is not in the range of the netcdf dataset'
    write(*,*)'model time       : ',gettime(int(t))
    write(*,*)'netcdf begin time: ',gettime(int(ftimes(1)))
    write(*,*)'netcdf end   time: ',gettime(int(ftimes(nframes)))
    stop
!!  added by zhuxm
  elseif(t > ftimes(nframes))then
    write(*,*)'error in bracket:'
    write(*,*)'model time is not in the range of the netcdf dataset'
    write(*,*)'model time       : ',gettime(int(t))
    write(*,*)'netcdf begin time: ',gettime(int(ftimes(1)))
    write(*,*)'netcdf end   time: ',gettime(int(ftimes(nframes)))
    stop
!! added by zhuxm
  endif

  do i=1,nframes-1
    if(ftimes(i) <= t .and. t <= ftimes(i+1)) then
      i1 = i
      i2 = i+1
      exit  ! zfeng add 10/30/2014
    endif
  end do
end subroutine bracket

subroutine extbracket(t,i1,i2)
  use utilities
  real(sp), intent(in) :: t
  integer, intent(out) :: i1
  integer, intent(out) :: i2
  !------------------------
  integer :: i

  if(nframex <= 1)  &
    stop 'error in extbracket: have no external forcing files have been read'

  if(t < extftime(1))then
    write(*,*)'error in extbracket:'
    write(*,*)'model time is not in the range of the netcdf dataset'
    write(*,*)'model time       : ',gettime(int(t))
    write(*,*)'netcdf begin time: ',gettime(int(extftime(1)))
    write(*,*)'netcdf end   time: ',gettime(int(extftime(nframex)))
    stop
  endif

  do i=1,nframex-1
    if(extftime(i) <= t .and. t <= extftime(i+1))then
      i1 = i
      i2 = i+1
    endif
  end do

end subroutine extbracket
  
subroutine exchange_forcing
  implicit none

  !----------------------------------------
  integer  :: ndims,dtype,i
  integer  :: d1,d2,d3,d4

  !$OMP PARALLEL DO PRIVATE(dtype, ndims)
  !loop over vars, exchage frame
  do i=1,frame(now)%nvars
    dtype = frame(now)%fdata(i)%dtype
    ndims = frame(now)%fdata(i)%ndims
    if(dtype == flt_type)then
!!zhuxm changed from if to select case
      select case(ndims)
      case(1)
        frame(now)%fdata(i)%f1 = frame(nex)%fdata(i)%f1
      case(2) 
        frame(now)%fdata(i)%f2 = frame(nex)%fdata(i)%f2
      case(3)
        frame(now)%fdata(i)%f3 = frame(nex)%fdata(i)%f3
      case(4)
        frame(now)%fdata(i)%f4 = frame(nex)%fdata(i)%f4
      case default
        write(*,*)'variable dimensions incorrect in exchange_forcing:'
        write(*,*)' ndims=',ndims
      end select
    else
      select case(ndims)
      case(1)
        frame(now)%fdata(i)%i1 = frame(nex)%fdata(i)%i1
      case(2)
        frame(now)%fdata(i)%i2 = frame(nex)%fdata(i)%i2
      case(3)
        frame(now)%fdata(i)%i3 = frame(nex)%fdata(i)%i3
      case(4)
        frame(now)%fdata(i)%i4 = frame(nex)%fdata(i)%i4
      case default
        write(*,*)'variable dimensions incorrect in exchange_forcing:'
        write(*,*)' ndims=',ndims
      end select
    endif
  enddo
  !$OMP END PARALLEL DO

end subroutine exchange_forcing


!========================================================================
! open the ext files
!   - make sure it exists
!   - read time (convert to sec if necessary.)
!========================================================================
subroutine open_ext_file(ffiles_in,nfls) 
  use utilities

  implicit none
  integer, intent(in) :: nfls
  character(len=fstr) :: ffiles_in(nfls),ffile_in
  integer             :: nfsfid(nfls),nfsnrec(nfls)
  
  integer             :: fid,varid,i,ierr,n
  character(len=mstr) :: msg
  character(len=fstr) :: dname,tunits

!Added on Jun 11 ,2009
  integer             :: nb,ne,tmpii,y
  real(sp)            :: tmpxx

!==============

  nframex=0 
  do n=1,nfls
     ffile_in=ffiles_in(n) !maybe it goes wrong
     msg = "error opening forcing file: "//trim(ffile_in)
     call ncdchk( nf90_open(trim(ffile_in),nf90_nowrite,fid),msg)

     nfsfid(n) = fid
     ffile(n)  = ffile_in

     msg = "reading number of dimensions from: "//trim(ffile_in)
     call ncdchk(nf90_inquire(fid, nDim,nVar,nAtt,uDid,fNum) ,msg)

     uD_extant = 1
     if(uDid /= 0)then
       call ncdchk(nf90_inquire_dimension(fid, uDid, dname, uD_extant ))
     endif
     nframex = uD_extant+nframex
     nfsnrec(n) = uD_extant
  end do

!=====================================
!Revised on Jun 11, 2009
!add year_cycle  by zhuxm
!=====================================
  !allocate space to hold times and read in 
  if(year_cycle==0)then   
     allocate(extftime(nframex))  ; extftime  = 0.0_sp
     allocate(fid_exft(nframex)) ; fid_exft = 0
     allocate(rec_exft(nframex)) ; rec_exft = 0
     allocate(fnr_exft(nframex)) ; fnr_exft = 0
  else
     allocate(extftime(nframex*year_cycle+1)) ; extftime = 0.0_sp   !two years
     allocate(fid_exft(nframex*year_cycle+1)) ; fid_exft = 0
     allocate(rec_exft(nframex*year_cycle+1)) ; rec_exft = 0
     allocate(fnr_exft(nframex*year_cycle+1)) ; fnr_exft = 0
  endif
!=====================================
  msg = "error reading time variable from netcdf file"
  uD_extant=0
  do n=1,nfls
     uD_extant=uD_extant+nfsnrec(n)
     ffile_id=nfsfid(n)
     call ncdchk( nf90_inq_varid(ffile_id,'time',varid),msg )

     allocate(tmptime(nfsnrec(n))) ; tmptime = 0.0
     call ncdchk( nf90_get_var(ffile_id, varid,tmptime ),msg)

     if(nf90_get_att(ffile_id, varid, 'units', tunits) == nf90_noerr)then
        if(index(tunits,'day') /= 0) tmptime = tmptime*day_2_sec
     endif

     extftime(uD_extant-nfsnrec(n)+1 : uD_extant) = tmptime(1:nfsnrec(n))
     fid_exft(uD_extant-nfsnrec(n)+1 : uD_extant) = ffile_id
     rec_exft(uD_extant-nfsnrec(n)+1 : uD_extant) = (/(i,i=1,nfsnrec(n))/)
     fnr_exft(uD_extant-nfsnrec(n)+1 : uD_extant) = n
     deallocate(tmptime)
  enddo
!Revised on Jun 11, 2009
!To make new time variable
  nb=0
  ne=0
  if(nframex > 1)then
    do i=1,nframex
      if(extftime(i)>= 0.0_sp .and. nb == 0) nb=i
      if(extftime(i)<= 365*24*3600.0_sp) ne=i
    end do
  endif

! check the recycle
  if(year_cycle/=0)then
    if(extftime(ne)-extftime(nb) == 365*24*3600.0_sp) ne=ne-1
  endif
!==================
  n=0
  do i=nb,ne
     n=n+1
     tmpxx=extftime(i)  !! selected the time > 0.0
     extftime(n)= tmpxx

     tmpii=fid_exft(i)
     fid_exft(n)= tmpii

     tmpii=rec_exft(i)
     rec_exft(n)= tmpii

     tmpii=fnr_exft(i)
     fnr_exft(n)= tmpii
  enddo

  if(year_cycle/=0)then
    do y=2,year_cycle
      do i=(y-1)*n+1,y*n   !! set time for next year
         extftime(i)=  extftime(i-n)+365*24.0*3600.0_sp
         rec_exft(i)=  rec_exft(i-n)
         fid_exft(i)=  fid_exft(i-n)
         fnr_exft(i)=  fnr_exft(i-n)
      enddo
!make it closed
      extftime(2*n+1)= extftime(1)+2*365*24.0*3600.0_sp
      rec_exft(2*n+1)= rec_exft(1)
      fid_exft(2*n+1)= fid_exft(1)
      fnr_exft(2*n+1)= fnr_exft(1)
    enddo
  endif

  if(nframex > 1)then
    do i=2,nframex
      if(ftimes(i)-ftimes(i-1) <= 0.0)then
         write(*,*)'netcdf time is not monotonically increasing'
         write(*,*)i,ftimes(i),ftimes(i-1)
      endif
    end do
  endif

  if(year_cycle/=0) nframex=year_cycle*n+1

end subroutine open_ext_file

!========================================================================
! get data by varname
!========================================================================
subroutine  get_data(vvname,p)

  character(len=*) :: vvname
  real(sp), pointer :: p(:)
  
  real(sp)  :: t1,t2,c1,c2,t
  integer, allocatable, dimension(:) :: start,vsize  
  integer :: fid,nAtts,vid,ndims,xtype,ttype,dtype
  integer :: dimids(NF90_MAX_VAR_DIMS)
  integer :: alldims(NF90_MAX_VAR_DIMS)
  integer :: dims(NF90_MAX_VAR_DIMS)
  integer :: j,jd
  character(len=10) :: dname

  t=gt
  jd=index(vvname,'_ext')-1

  call extbracket(t,extf1,extf2)
  !set interpolation coefficients
  t1 = extftime(extf1)
  t2 = extftime(extf2)
  if(t2 > t1)then
    c2 = (t - t1)/(t2-t1)
    c1 = 1.-c2
  else
    c1 = (t - t2)/(t1-t2)
    c2 = 1.-c1
  endif

  fid=fid_exft(extf1)
  call ncdchk( nf90_inq_varid(fid,vvname(1:jd),vid),'erro1','erro2' )
  call ncdchk( nf90_inquire_variable(fid,vid,vvname(1:jd),xtype,ndims,dimids,nAtts))
  do j=1,ndims
    call ncdchk(nf90_inquire_dimension(fid, dimids(j), dname, alldims(j)) )
  enddo
  allocate(start(ndims),vsize(ndims))
  if(ndims==2)then
     allocate(ext_data1(alldims(1),1))
     allocate(ext_data2(alldims(1),1))
  endif
  start(1:ndims-1) = 1
  vsize(1:ndims-1) = alldims(1:ndims-1)

  start(ndims) = rec_exft(extf1)
  vsize(ndims) = 1
  call ncdchk(  nf90_get_var(fid,  vid, ext_data1,start,vsize) )
  start(ndims) = rec_exft(extf2)
  call ncdchk(  nf90_get_var(fid,  vid, ext_data2,start,vsize) )

     P = c1*ext_data1(:,1)+c2*ext_data2(:,1)
!      if(ndims==3) call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f3,start,vsize) )
!      if(ndims==4) call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f4,start,vsize) )
      deallocate(start)
      deallocate(vsize)
      deallocate(ext_data1)
      deallocate(ext_data2)
end subroutine get_data

End Module mod_forcing
