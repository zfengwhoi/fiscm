Module mod_igroup

!=======================================================================
! Fiscm Igroup 
!
! Description
!   Defines a group type for Fiscm 
!     -initialize the group
!     -add state variables to the group
!     -summarize the group
!    
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!  2011/7/28       Xueming Zhu
!  2014/10/17     Zhixuan Feng
!  -- Modify function 'group_' to include states for structured
!     forcing calculation (e.g., particle indices i,j,k).
!  2015/03/10     Zhixuan Feng
!  -- Modify function 'group_' to include a state for subsurface 
!     phytoplankton maximum
!=======================================================================

	use gparms
	use mod_pvar

	Implicit None
	integer             :: nind
	character(len=fstr) :: fname_out
	real(sp)            :: group_deltat,tlast_out

	interface add_state
	  module procedure add_state_fvec
	  module procedure add_state_ivec
	end interface
  
	interface get_state
	  module procedure get_state_fvec
	  module procedure get_state_ivec
	end interface

!----------------------------------------------------------------
! Group type
!----------------------------------------------------------------
	type igroup 
	  integer             :: id                         !unique group id
	  integer             :: Tnind                      !total number of group individuals
	  integer             :: nind                       !active subset
	  integer             :: space_dim                  !spatial dimension of problem 0,2,3
	  character(len=fstr) :: name                       !group name (e.g. species)
	  integer             :: hdiff_type                 !horizontal diffusion type (HDIFF_NONE,HDIFF_CONST,HDIFF_VISSER)
	  real(sp)            :: hdiff_const_val            !constant value  
	  integer             :: vdiff_type
	  real(sp)            :: vdiff_const_val
	  integer             :: vdiff_substeps
	  integer             :: intvl_bio                  
	  real(sp)            :: DT_bio
	  logical             :: biology
          logical             :: bio_food
          integer             :: food_source
          integer             :: ActLayInd
	  integer             :: intvl_out
	  real(sp)            :: tlast_out
	  real(sp)            :: start_out
	  integer             :: frame_out
	  integer             :: fid_out
	  character(len=fstr) :: fname_out
	  character(len=fstr) :: statefile
	  type(pvar_list)     :: state
	  integer             :: nstate
	  integer             :: next
	  character(len=fstr) :: ext_var(max_state_vars,2)   !! 200
	  integer             :: nstate_ud
	end type igroup 

	type(igroup), allocatable :: igroups(:)

contains
!---------------------------------------------------------
!Create the initial group type (passive particle)
!---------------------------------------------------------
function group_(fid,id,deltaT) result(g)
  type(igroup) :: g
  integer, intent(in) :: fid
  integer, intent(in) :: id
  real(sp),intent(in) :: deltaT
  integer :: ios,i,ns
  real(sp) :: fval
  integer  :: ival

  !read the group namelist from unit fid
  read(unit=fid,nml=nml_group,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read group namelist from '//trim(runcontrol)
    stop
  endif
  !check problem size
  if(Tnind < 1)then
    write(*,*) 'error creating group: size must be greater than 0'
    stop 
  endif
  !check problem dimension
  if(space_dim < 0 .or. space_dim > 3)then
    write(*,*)'fatal error: space dim must be 0,2,3'
    stop
  endif

  !set time step
  g%id              = id
  g%Tnind           = Tnind
  g%nind            = Tnind
  g%space_dim       = space_dim
  g%name            = trim(group_name)
  g%hdiff_type      = hdiff_type
  g%hdiff_const_val = hdiff_const_val
  g%vdiff_type      = vdiff_type
  g%vdiff_const_val = vdiff_const_val
  g%vdiff_substeps  = vdiff_substeps
  g%intvl_bio       = intvl_bio
  g%DT_bio          = intvl_bio*deltaT
  g%biology         = biology
  g%bio_food        = bio_food
  g%food_source     = food_source
  g%ActLayInd       = ActLayInd
  g%intvl_out       = intvl_out
  g%tlast_out       = 0.0_sp
  g%start_out       = start_out
  g%frame_out       = 0
  g%fid_out         = 0
  g%fname_out       = " "
  g%statefile       = trim(ini_file)
  g%state           = pvar_list_() 
  g%nstate          = 0 
  g%next            = 0
  g%nstate_ud       = nstate_ud

 if(.not. structured) then ! for unstructured computation
  !(x,y,z) - particle position
  if(g%space_dim > 1)then
    fval=0.0_sp;call add_state(g,'x','x location','m',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'y','y location','m',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'h','bathymetry','m',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'pathlength','integrated trajectory','m',NETCDF_YES,fval)
    ival=0 ;    call add_state(g,'cell','cell containing particle','-',NETCDF_YES,ival)
    fval=0.0_sp;call add_state(g,'u','x-velocity','m/s',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'v','y-velocity','m/s',NETCDF_YES,fval)
  endif
  if(g%space_dim == 3)then
    fval=0.0_sp;call add_state(g,'w','w-velocity'  ,'m/s',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'z','z location'  ,'m'  ,NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'s','s coordinate','-'  ,NETCDF_YES,fval) 
  endif
 endif

 if(structured)  then ! for structured computation
  !(x,y,z) - particle position
  if(g%space_dim > 1)then
    fval=0.0_sp;call add_state(g,'x','x location','m or deg',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'y','y location','m or deg',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'h','bathymetry','m',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'pathlength','integrated trajectory','m',NETCDF_YES,fval)
    ival=0;     call add_state(g,'i','i-index of nearest S-Grid point','-',NETCDF_YES,ival)
    ival=0;     call add_state(g,'j','j-index of nearest S-Grid point','-',NETCDF_YES,ival)
    ival=0;     call add_state(g,'mask','mask of nearest S-Grid point','-',NETCDF_YES,ival)
    fval=0.0_sp;call add_state(g,'u','x-velocity','cm/s',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'v','y-velocity','cm/s',NETCDF_YES,fval)
  endif
  if(g%space_dim == 3)then
    ival=0;     call add_state(g,'k','k-index of level below particle','-',NETCDF_YES,ival)
    fval=0.0_sp;call add_state(g,'w','w-velocity'  ,'m/s',NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'z','z location'  ,'m'  ,NETCDF_YES,fval)
    fval=0.0_sp;call add_state(g,'s','s coordinate','-'  ,NETCDF_YES,fval)
  endif
 endif


  !(status) - particle status
  ival=0;call add_state(g,'status','particle status','-',NETCDF_YES,ival)

  !(status) - spawn time
  if(g%biology) then
     fval=0.0_sp
     call add_state(g,'tspawn','time of spawning','seconds',NETCDF_YES,fval)
  endif

! zfeng add 03/10/2015 Subsurface Phytoplankton Maximum
!  if(g%biology .AND. bio_food .AND. bio_SPM) then
!     fval=0.0_sp
!     call add_state(g,'SubPhytoMax','Subsurface Phytoplankton Max','mg chla/m3',NETCDF_YES,fval)
!  endif  

end function group_

!---------------------------------------------------------
!add user defined states to the group
!---------------------------------------------------------
subroutine group_addstates(g) 
  type(igroup)       :: g
  integer            :: ns,ios

  if(g%nstate_ud < 1)return

  rewind(iunit)
  do ns = 1,g%nstate_ud

    !read the group namelist from unit fid
    read(unit=iunit,nml=nml_statevar,iostat=ios)
	if((index(trim(state_varname),'chl')/=0  .or.   &
	    index(trim(state_varname),'CHL')/=0) .and. bio_fd ==0) cycle
    if(state_vartype == int_type)then
      call add_state(g,trim(state_varname),trim(state_longname),trim(state_units), &
                   state_netcdf_out,state_initval_int,trim(state_from_ext_var))
    else if(state_vartype == flt_type)then
      call add_state(g,trim(state_varname),trim(state_longname),trim(state_units), &
                   state_netcdf_out,state_initval_flt,trim(state_from_ext_var))
    else
      write(*,*)'fatal error: not setup for variable type ',state_vartype
	  stop ' in group_addstates'
    endif
  end do
end subroutine group_addstates

!---------------------------------------------------------
!add real states to the group
!---------------------------------------------------------
subroutine add_state_fvec(g,varname,longname,units,output,init_val,ext_name)
   type(igroup)                 :: g 
   character(len=*)             :: varname
   character(len=*)             :: longname
   character(len=*)             :: units
   integer                      :: output
   real(sp)                     :: init_val
   character(len=*), optional   :: ext_name
   type(pvar)                   :: new
   integer                      :: i

   !need to make sure we arent duplicating state names!
   allocate(new%fvar(g%Tnind)) 
   do i=1,g%Tnind
      new%fvar(i) = init_val
   end do
    
   !check output flag
   if(output < 0 .or. output > 1)then
     write(*,*)'error adding state variable: ',varname
     write(*,*)'output argument must be 0 (no output), or 1 (output)'
     stop
   endif
 
   new%idim     = g%Tnind
   new%istype   = flt_type
   new%varname  = trim(varname)
   new%longname = trim(longname)
   new%units    = trim(units)
   if(present(ext_name))then
     new%from_ext_var = trim(ext_name)
     if(ext_name /= 'NONE' .and. ext_name /= 'none')then
	   if(index(ext_name,'_ext') /= 0 .and. n_extfile==0)then
	      write(*,*)'You will try to read variable ', trim(ext_name)
		  write(*,*)'you have not set external files with this variable!'
		  stop 'error in add_state_fvec'
	   endif
       g%next = g%next + 1
       g%ext_var(g%next,1) = trim(varname)
       g%ext_var(g%next,2) = trim(ext_name)
     endif
   else
     new%from_ext_var = "NONE"
   endif

   new%output   = output

   call add_pvar_node(g%state,new)
   g%nstate = g%nstate + 1;
end subroutine add_state_fvec

!---------------------------------------------------------
!add integer states to the group
!---------------------------------------------------------
subroutine add_state_ivec(g,varname,longname,units,output,init_val,ext_name)
   type(igroup)                 :: g 
   character(len=*)             :: varname
   character(len=*)             :: longname
   character(len=*)             :: units
   integer                      :: output
   integer                      :: init_val
   character(len=*), optional   :: ext_name
   type(pvar)                   :: new
   integer                      :: i

   !need to make sure we arent duplicating state names!
 
   allocate(new%ivar(g%Tnind)) 
   
   do i=1,g%Tnind
     new%ivar(i) = init_val
   end do
  
   !check output flag
   if(output < 0 .or. output > 1)then
     write(*,*)'error adding state variable: ',varname
     write(*,*)'output argument must be 0 (no output), or 1 (output)'
     stop
   endif

   new%idim     = g%Tnind
   new%istype   = int_type
   new%varname  = trim(varname)
   new%longname = trim(longname)
   new%units    = trim(units)
   if(present(ext_name))then
     new%from_ext_var = trim(ext_name)
     if(ext_name /= 'NONE' .and. ext_name /= 'none')then
	   if(index(ext_name,'_ext') /= 0 .and. n_extfile==0)then
	      write(*,*)'You will try to read variable ', trim(ext_name)
		  write(*,*)'you have not set external files with this variable!'
		  stop 'error in add_state_ivec'
	   endif
       g%next = g%next + 1
       g%ext_var(g%next,1) = trim(varname)
       g%ext_var(g%next,2) = trim(ext_name)
       write(*,*)'ADDING EXT VAR: ',trim(ext_name),'  to : ',trim(varname)
     endif
   else
     new%from_ext_var = "NONE"
   endif

   new%output   = output

   call add_pvar_node(g%state,new)
   g%nstate = g%nstate + 1;

end subroutine add_state_ivec

!---------------------------------------------------------
!get real states from the group by name
!---------------------------------------------------------
subroutine get_state_fvec(varname,g,v)
   character(len=*)    :: varname
   type(igroup)        :: g
   real(sp), pointer   :: v(:)
   type(pvar), pointer :: p
   
   p => get_pvar(g%state,varname)
   v => p%fvar
   
end subroutine get_state_fvec

!---------------------------------------------------------
!get integer states from the group by name
!---------------------------------------------------------
subroutine get_state_ivec(varname,g,v)
   character(len=*)    :: varname
   type(igroup)        :: g
   integer, pointer    :: v(:)
   type(pvar), pointer :: p

   p => get_pvar(g%state,varname)
   v => p%ivar
end subroutine get_state_ivec

!----------------------------------------------------
! determine which variables, if any are needed for
! forcing
!----------------------------------------------------
subroutine get_ext_varnames(ng,g,lsize,varlist)
  integer, intent(in)      :: ng
  type(igroup), intent(in) :: g(ng)
  integer, intent(inout)   :: lsize
  character(len=*)         :: varlist(max_state_vars)   !! 200

  integer                  :: i,j 
  character(len=fstr)      :: var

  do i=1,ng
     do j=1,g(i)%next
        var = g(i)%ext_var(j,2)
        if(index(var,'_ext') /= 0 )cycle
        if(lsize > 0)then
          if(.not. inlist(var,varlist,lsize))then
             lsize = lsize + 1
             varlist(lsize) = var
          endif
        else
          lsize = lsize + 1
          varlist(lsize) = var
        endif
     end do
  end do

end subroutine get_ext_varnames

!---------------------------------------------------------
!check the variable whether is in variables list by name
!---------------------------------------------------------
function inlist(var,varlist,ndim) result(check) 
  character(len=*)    :: var
  integer, intent(in) :: ndim
  character(len=*)    :: varlist(ndim)
  logical :: check
  integer :: i
  check = .false.
  do i=1,ndim 
    if(trim(var) == trim(varlist(i))) check = .true.
  enddo
end function inlist

!---------------------------------------------------------
!print group information
!---------------------------------------------------------
subroutine print_group_summary(g) 
   use utilities
   type(igroup)        :: g

   call drawline("-")
   write(*,*)'Group Summary: ',g%id
   call drawline("-")
   write(*,'(A21,I10)')'group id number   :: ',g%id
   write(*,'(A21,A20)')'group name        :: ',g%name
   write(*,'(A21,I10)')'total individuals :: ',g%Tnind
   write(*,'(A21,I10)')'assigned slots    :: ',g%nind
   write(*,'(A21,I10)')'spatial dimension :: ',g%space_dim
  
   select case(g%hdiff_type)
     case(HDIFF_NONE) 
       write(*,'(A21,A20)')'hdiff type        :: ','HDIFF_NONE'
     case(HDIFF_CONSTANT) 
       write(*,'(A21,A20)')'hdiff type        :: ','HDIFF_CONSTANT'
       write(*,'(A21,F10.6)')'hdiff constant  :: ',g%hdiff_const_val
     case(HDIFF_VARIABLE)
       write(*,'(A21,A20)')'hdiff type        :: ','HDIFF_VARIABLE'
   end select
   select case(g%vdiff_type)
     case(VDIFF_NONE) 
       write(*,'(A21,A20)')'vdiff type        :: ','VDIFF_NONE'
     case(VDIFF_VARIABLE)   
       write(*,'(A21,A20)')'vdiff type        :: ','VDIFF_VARIABLE'
       write(*,'(A21,I10)')'vdiff substeps    :: ',vdiff_substeps
     case(VDIFF_SPLINED)
       write(*,'(A21,A20)')'vdiff type        :: ','VDIFF_SPLINED'
       write(*,'(A21,I10)')'vdiff substeps    :: ',vdiff_substeps
     case(VDIFF_BINNED)
       write(*,'(A21,A20)')'vdiff type        :: ','VDIFF_BINNED'
       write(*,'(A21,I10)')'vdiff substeps    :: ',vdiff_substeps
   end select
   if(g%biology) then
     write(*,'(A21,A10)')'biology           :: ','ACTIVE'
     write(*,'(A21,F10.2)')'bio time step(s)  :: ',g%DT_bio
   else
     write(*,'(A21,A10)')'biology           :: ','INACTIVE'
   endif
  
   if(g%bio_food) then
     write(*,'(A21,A10)')'food_dependency   :: ','ACTIVE'
   else
     write(*,'(A21,A10)')'food_dependency   :: ','INACTIVE'
   endif

   if(g%food_source == 1) then
     write(*,'(A21,A10)')'food_source       :: ','Phyt'
   elseif(g%food_source == 2) then
     write(*,'(A21,A10)')'food_source       :: ','Phyt+MZ'
   elseif(g%food_source == 3) then
     write(*,'(A21,A10)')'food_source       :: ','Max Pht'
     write(*,'(A21,A10)')'Active Layers     :: ','0-',g%ActLayInd
   elseif(g%food_source == 4) then
     write(*,'(A21,A10)')'food_source       :: ','Max Pht+MZ'
     write(*,'(A21,A10)')'Active Layers     :: ','0-',g%ActLayInd
   else
     write(*,'(A21,A10)')'food_source       :: ','UNKNOWN'      
   endif

   write(*,'(A21,I10)')  'output intval(its):: ',g%intvl_out
   write(*,'(A21,F10.2)')'output starts at  :: ',g%start_out
   write(*,'(A21,I10)')'NO. of state vars   :: ',g%nstate

   call print_state_vars(g%state)

end subroutine print_group_summary

!---------------------------------------------------------
! check status of particles, report, return false if
! no active particles left
!---------------------------------------------------------
function checkstatus(ng,g,time,its) result(validsim)
   use utilities
   implicit none
   logical                 :: validsim 
   integer, intent(in)     :: ng,its
   type(igroup),intent(in) :: g(ng)
   real(sp), intent(in)    :: time
   integer                 :: nACTIVE,nSETTLED,nDEAD
   integer                 :: nUNKNOWN,nEXITED,nTOTAL
   integer, pointer        :: istatus(:)
   integer                 :: n,NI
   integer, allocatable    :: mark(:)
   integer, save           :: dumphead = 1

   !assume we still have some active particles
   validsim = .true.

   !initialize counters
   nACTIVE  = 0 ;   nSETTLED = 0 ;   nDEAD    = 0
   nUNKNOWN = 0 ;   nEXITED  = 0 ;   nTOTAL   = 0

   do n=1,ng
     NI = g(n)%Nind
     call get_state('status',g(n),istatus)
     allocate(mark(NI)) ; mark = 0
     where(istatus(1:NI) == SETTLED)mark = 1
       nSETTLED = nSETTLED + sum(mark)
       mark = 0
     where(istatus(1:NI) == DEAD)mark = 1
       nDEAD = nDEAD + sum(mark)
       mark = 0
     where(istatus(1:NI) == ACTIVE)mark = 1
       nACTIVE = nACTIVE + sum(mark)
       mark = 0
     where(istatus(1:NI) == EXITED)mark = 1
       nEXITED = nEXITED + sum(mark)
       mark = 0
     where(istatus(1:NI) == UNKNOWN)mark = 1
       nUNKNOWN = nUNKNOWN + sum(mark)
       mark = 0
     deallocate(mark)
     nTOTAL = nTOTAL + NI 
   end do
   if(mod(dumphead-1,10)==0)then
     write(*,*)
     write(*,102)
     dumphead = 1
   endif
   if(its>0 .and. mod(its,ireport)==0)then
     write(*,101)its,gettime(int(time)),nTOTAL,nSETTLED,nDEAD,nACTIVE,nEXITED,nUNKNOWN
     dumphead = dumphead + 1
   endif
   if(nACTIVE + nUNKNOWN == 0)then
      validsim = .false.
      call drawline('-')
      write(*,*)'no active particles left in the simulation: shutting down early'
   endif
  
 101 format(I8,2X,A13,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8) 
 102 format("   IINT     simtime   ",6x,"TOTAL  ",1x,"SETTLED",1x,"  DEAD  ",1X," ACTIVE ",1X," EXITED ",1x," UNKNOWN") 
end function checkstatus

!---------------------------------------------------------
!get spawn limits
!---------------------------------------------------------
subroutine get_spawn_range(ng,g,tspawn_min,tspawn_max)
  integer, intent(in)      :: ng
  type(igroup), intent(in) :: g(ng)
  real(sp), intent(out) :: tspawn_min,tspawn_max
  integer :: np ,n
  real(sp), pointer :: tspawn(:) 

  !initialize spawning times to extreme numbers
  tspawn_min =  hugenum
  tspawn_max = -hugenum
  

  !search through groups for earliest and latest spawning times
  do n=1,ng
    np = g(n)%Nind
    call get_state('tspawn',g(n),tspawn)
    tspawn_min = min( minval(tspawn(1:np)) , tspawn_min)
    tspawn_max = max( maxval(tspawn(1:np)) , tspawn_max)
  end do
end subroutine get_spawn_range

End Module mod_igroup
