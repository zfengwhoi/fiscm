Module mod_output

!=======================================================================
! Fiscm NetCDF Output Routines 
!
! Description:
!   Routines to output fiscm data to NetCDF 
!
! Comments:  Requires fortran90 NetCDF 3.x libraries
!    
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!  2011/7/28       Xueming Zhu
!=======================================================================

	use gparms
	use mod_igroup
	use utilities,only : ncdchk

	implicit none

	integer, private  :: time_did
	integer, private  :: nlag_did
	integer, private  :: time_vid
	integer, private  :: dynm2d(2)
	integer, private  :: dynm1d(1)
    character(len=mstr), private  :: msg,msg2

	public

	interface add_cdfglobs
	   module procedure add_cdfglobs_int
	   module procedure add_cdfglobs_flt
	   module procedure add_cdfglobs_str
	end interface

	contains

!-------------------------------------------
! output driver
!-------------------------------------------
subroutine cdf_out(ng,g,itnum,time,otype)
  integer, intent(in) :: ng
  type(igroup),intent(inout), dimension(ng) :: g
  integer , intent(in) :: itnum
  real(sp), intent(in) :: time 
  integer,  intent(in) :: OTYPE
  integer :: n
 
  select case(otype)
  !write the netcdf header file
  case(NCDO_HEADER)             
    do n=1,ng
      call write_header(g(n))    
    end do
  !add the user-defined state variables to the output header
  case(NCDO_ADD_STATES)
    do n=1,ng
      call add_states_to_nc(g(n)) !
    end do
  !output the state variables 
  case(NCDO_OUTPUT)
    do n=1,ng
      if(mod(itnum-1,g(n)%intvl_out)==0)call output_group(g(n),time)
    end do
  end select
end subroutine cdf_out

!-------------------------------------------
! write netcdf header
!-------------------------------------------
subroutine write_header(g)

  type(igroup), intent(inout) :: g
  character(len=fstr)  :: fname,fnametmp
  character(len=3)     :: num,nid
  integer              :: ierr,fid,nc
 
  g%frame_out = 1
  write(nid,'(I3.3)') g%id

  !construct the name (eventually use some sim string/directory loc)
  fnametmp = 'fiscm_group_'//nid//'.nc'
  nc=0
  fname=fnametmp
  do
    inquire(FILE=trim(fname),EXIST=fexist)
	if(fexist)then
	   nc=nc+1
	   write(num,'(I3.3)')nc
	   fname=trim(fnametmp)//'-'//num
	else
	   exit
	endif
  enddo

  g%fname_out = trim(fname)
  !open the file - write access
  call ncdchk( nf90_create(trim(fname),nf90_clobber,fid) )
  g%fid_out = fid

  !global attributes
  call ncdchk( nf90_put_att(fid,nf90_global,"title", "FVCOM PARTICLE TRACKING")) 
  call ncdchk( nf90_put_att(fid,nf90_global,"Conventions", "CF-1.0")) 
  call ncdchk( nf90_put_att(fid,nf90_global,"source"      ,"Particle FVCOM_3.0"))  !trim(FISCM_VERSION)) )
  call ncdchk( nf90_put_att(fid,nf90_global,"group"     ,g%id) )
  call ncdchk( nf90_put_att(fid,nf90_global,"DT_bio"    ,g%DT_bio) )

  !define the time dimension   
  call ncdchk(nf90_def_dim(fid,"time",nf90_unlimited, time_did) )
  dynm1d = (/time_did/)

  !time variable
  call ncdchk( nf90_def_var(fid,"time",nf90_double,dynm1d, time_vid) )
  call ncdchk( nf90_put_att(fid, time_vid,"long_name","time") )
  call ncdchk( nf90_put_att(fid, time_vid,"units","days since 0.0") )
  call ncdchk( nf90_put_att(fid, time_vid,"time_zone","none") )

  !lag dimension
  call ncdchk(nf90_def_dim(fid,"nlag",g%Tnind, nlag_did) )
  dynm2d = (/nlag_did,time_did/)

  !close the file
  call ncdchk( nf90_close(g%fid_out) )
end subroutine write_header

!-----------------------------------------------------
! add user defined states variables header to netcdf
!-----------------------------------------------------
subroutine add_states_header_to_nc(plist,fid,dims)
  type(pvar_list) :: plist
  integer, intent(in) :: fid
  integer, intent(in) :: dims(2)
  type(pvar_node), pointer :: plst,pnew
  type(pvar), pointer :: dum
  integer :: ierr
  
  plst => plist%begin
  pnew => plst%next

  if(.not.associated(plst%next))then
    write(*,*)'plist has no nodes'
    stop
  endif
  do 
    if(.not.associated(plst%next)) exit
      !dump the header if user selected to dump
    if(pnew%v%output == NETCDF_YES)then
       dum => pnew%v
       select case(dum%istype)
       case(flt_type)
          ierr = nf90_def_var(fid,dum%varname,nf90_double,dims,dum%nc_vid)
          if(ierr /= nf90_noerr) then
            write(*,*)'error defining netcdf header for variable: ',dum%varname
            write(*,*)trim(nf90_strerror(ierr)) ; stop 
          end if
       case(int_type)
          ierr = nf90_def_var(fid,dum%varname,nf90_int  ,dims,dum%nc_vid)
          if(ierr /= nf90_noerr) then
            write(*,*)'error defining netcdf header for variable: ',dum%varname
            write(*,*)trim(nf90_strerror(ierr)) ; stop 
          end if
       case default
          write(*,*)'not setup for netcdf dumping state variable of type:',dum%istype
          stop
       end select
       ierr = nf90_put_att(fid,dum%nc_vid,"long_name",dum%longname)
       if(ierr /= nf90_noerr) then
          write(*,*)'error defining netcdf header for variable: ',dum%varname
          write(*,*)trim(nf90_strerror(ierr)) ; stop 
       end if

       ierr = nf90_put_att(fid,dum%nc_vid,"units",dum%units)
       if(ierr /= nf90_noerr) then
          write(*,*)'error defining netcdf header for variable: ',dum%varname
          write(*,*)trim(nf90_strerror(ierr)) ; stop 
       end if
    endif
    plst => pnew
    pnew => pnew%next
  end do
end subroutine add_states_header_to_nc

!-------------------------------------------
! write netcdf data
!-------------------------------------------
subroutine write_cdf_data(plist,fid,frame)
  type(pvar_list) :: plist
  integer, intent(in) :: fid
  integer, intent(in) :: frame
  type(pvar_node), pointer :: plst,pnew
  type(pvar), pointer :: dum
  integer :: ierr
  integer :: dims(2)

  dims(1) = 1
  dims(2) = frame
  
  plst => plist%begin
  pnew => plst%next

  if(.not.associated(plst%next))then
    stop 'plist has no nodes'
  endif
  do 
    if(.not.associated(plst%next)) exit
    !dump the variable if user selected to dump
    if(pnew%v%output == NETCDF_YES)then
       dum => pnew%v
       select case(dum%istype)
       case(flt_type)
          ierr = nf90_put_var(fid, dum%nc_vid, dum%fvar,START=dims)
          if(ierr /= nf90_noerr) then
            write(*,*)'error dumping frame for variable: ',dum%varname
            write(*,*)trim(nf90_strerror(ierr))
			stop 
          end if
       case(int_type)
          ierr = nf90_put_var(fid, dum%nc_vid, dum%ivar,START=dims)
          if(ierr /= nf90_noerr) then
             write(*,*)'error dumping frame for variable: ',dum%varname
             write(*,*)trim(nf90_strerror(ierr))
			 stop 
          end if
       case default
          write(*,*)'not setup for netcdf dumping state variable of type:',dum%istype
          stop
        end select

    endif
    plst => pnew
    pnew => pnew%next
  end do
end subroutine write_cdf_data

!-----------------------------------------------
! add state variable definitions to netcdf file
!-----------------------------------------------
subroutine add_states_to_nc(g)

   implicit none
   type(igroup), intent(inout) :: g
   integer              :: ierr

   !open the file
   msg ='error opening'//trim(g%fname_out)
   call ncdchk(nf90_open(g%fname_out,nf90_write,g%fid_out),msg)
   call ncdchk(nf90_redef(g%fid_out) )

   !write the state variables slated for output
   call add_states_header_to_nc(g%state,g%fid_out,dynm2d)

   !close the file
   call ncdchk( nf90_close(g%fid_out) )
end subroutine add_states_to_nc

!-----------------------------------------------
! add additional parameters to netcdf global
! parameter list -- FLOATS
!-----------------------------------------------
subroutine add_cdfglobs_flt(g,desc,val)

    implicit none

	type(igroup), intent(inout) :: g
	character(len=*)     :: desc
	real(sp), intent(in) :: val
	integer              :: ierr

    !open the file
	ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
	if(ierr /= nf90_noerr)then
	 write(*,*)'error opening', trim(g%fname_out)
	  write(*,*)trim(nf90_strerror(ierr))
	endif
	call ncdchk( nf90_redef(g%fid_out) )

    !write variable to global
	call ncdchk( nf90_put_att(g%fid_out,nf90_global,trim(desc),val) )

    !close the file
	call ncdchk( nf90_close(g%fid_out) )
end subroutine add_cdfglobs_flt

!-----------------------------------------------
! add additional parameters to netcdf global
! parameter list -- INTEGERS
!-----------------------------------------------
subroutine add_cdfglobs_int(g,desc,val)
	type(igroup), intent(inout) :: g
	character(len=*)     :: desc
	integer , intent(in) :: val
	integer              :: ierr

    !open the file
	ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
	if(ierr /= nf90_noerr)then
	  write(*,*)'error opening', trim(g%fname_out)
	  write(*,*)trim(nf90_strerror(ierr))
	  stop
	endif
	call ncdchk( nf90_redef(g%fid_out) )

    !write variable to global
	call ncdchk( nf90_put_att(g%fid_out,nf90_global,trim(desc),val) )

    !close the file
	call ncdchk( nf90_close(g%fid_out) )
end subroutine add_cdfglobs_int

!-----------------------------------------------
! add additional parameters to netcdf global
! parameter list -- STRINGS
!-----------------------------------------------
subroutine add_cdfglobs_str(g,desc,val)
   type(igroup), intent(inout) :: g
   character(len=*)     :: desc,val
   integer              :: ierr

   !open the file
   ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
   if(ierr /= nf90_noerr)then
      write(*,*)'error opening', trim(g%fname_out)
      write(*,*)trim(nf90_strerror(ierr))
   endif
   call ncdchk(nf90_redef(g%fid_out))

   !write variable to global
   call ncdchk(nf90_put_att(g%fid_out,nf90_global,trim(desc),val))

   !close the file
   call ncdchk(nf90_close(g%fid_out))
end subroutine add_cdfglobs_str

!-----------------------------------------------
! output time dependent vars (states and time) 
! to netcdf files
!-----------------------------------------------
subroutine output_group(g,time)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: time
  integer              :: ierr
  integer              :: dims(1)

  !set frame
  dims(1) = g%frame_out

  !debug
  !write(*,*)'dumping: ',g%id,' file: ',trim(g%fname_out),' frame: ',g%frame_out

  !open the file
  ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
  if(ierr /= nf90_noerr)then
    write(*,*)'error opening', trim(g%fname_out)
    write(*,*)trim(nf90_strerror(ierr))
	stop
  endif

  !dump time
  call ncdchk( nf90_put_var(g%fid_out,time_vid,time/(3600*24),START=dims) )

  !dump state variable data
  call write_cdf_data(g%state,g%fid_out,g%frame_out)

  !increment frame counter
  g%frame_out = g%frame_out + 1

  !close up
  call ncdchk( nf90_close(g%fid_out) )
end subroutine output_group

end module mod_output
