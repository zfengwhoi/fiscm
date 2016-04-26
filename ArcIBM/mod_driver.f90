Module mod_driver

!=======================================================================
! FISCM/BISCM Drivers
!
! Description
!  - Main driver for calling ocean-model-specific advect/diffuse/find routines 
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!  2011/7/28       Xueming Zhu
!  2014/10/29      Zhixuan Feng
!=======================================================================

	use gparms
	use mod_igroup
	use fvcom_driver
	use biomas_driver
	use mod_bio

	implicit none
	contains

!  s from z initialization
	subroutine sz_ini(ng,g)

	  integer     , intent(in)                   :: ng
	  type(igroup), intent(inout), dimension(ng) :: g
	  integer       :: n,np

	  do n=1,ng
	     np = g(n)%nind
	     if(g(n)%space_dim == 3) call sz_trans(np,g(n))
	  end do
	end subroutine sz_ini

!----------------------------------------------------
! advection-diffusion driver
!   call the appropriate advection and diffusion
!   routines, group-dependent
!
! zfeng 10/29/2014
! add 'structured' to provide options using structured forcing
!----------------------------------------------------
	subroutine do_adv_diff(ng,g,dT,time)

	  implicit none
	  integer,     intent(in)                  :: ng
	  type(igroup),intent(inout),dimension(ng) :: g
	  real(sp),intent(in)                      :: dT,time
	  integer                                  :: n,np
	  real(sp),pointer                         :: x(:),y(:)

	  do n=1,ng
        !set group dimension
	    np = g(n)%Tnind
	   select case(g(n)%space_dim)
	    case(0)
	      cycle
	   case(1)
	      write(*,*)'error in driver'
	      write(*,*)'1D simulation not setup'
	      stop
    !--------------------------------------------------------------------
	   case(2) ! 2D Problem
    !--------------------------------------------------------------------
       if(.NOT. structured) then
	      call advect2D( g(n) , dT,g(n)%nind)
	      if(g(n)%hdiff_type == HDIFF_CONSTANT)then
		     call rw_hdiff_constant( g(n) ,dT) 
	      elseif(g(n)%hdiff_type == HDIFF_VARIABLE)then
		     call rw_hdiff_variable( g(n),dT)
		  endif
	   elseif(structured) then
		     call advect2dRK4( g(n), dT, g(n)%nind, time)
		! write(*,*) 'To add avection2D!'

	   endif
    !--------------------------------------------------------------------
	    case(3) ! 3D Problem
    !--------------------------------------------------------------------
	    if(.NOT. structured) then

	      call advect3D( g(n) ,dT ,g(n)%nind,time)
	      if(g(n)%hdiff_type == HDIFF_CONSTANT)then
		     call rw_hdiff_constant( g(n) ,dT) 
	      elseif(g(n)%hdiff_type == HDIFF_VARIABLE)then
		     call rw_hdiff_variable( g(n),dT)
		  endif

	      if(g(n)%vdiff_type == VDIFF_VARIABLE)then
		     call rw_vdiff(g(n), dT, g(n)%vdiff_substeps)
	      elseif(g(n)%vdiff_type == VDIFF_SPLINED )then
		     call rw_vdiff_splined(g(n), dT, g(n)%vdiff_substeps)  
	      elseif(g(n)%vdiff_type == VDIFF_BINNED  )then
		     call rw_vdiff_binned(g(n), dT, g(n)%vdiff_substeps)
		  endif

		elseif(structured) then
		  SELECT CASE(TimeStepping)
		  CASE(1)
		       call advect3dEuler( g(n), dT, g(n)%nind,time)
		  CASE(2)
		       write(*,*) "Add Rk2 scheme later!"
		       stop
		  CASE(4)
		       call advect3dRK4( g(n) ,dT ,g(n)%nind,time)
		  CASE DEFAULT
		  write(*,*) "Error in 'mod_driver': Choose a right time stepping scheme!"
		  write(*,*) "1: Euler; 2: RK2; 4: RK4."
		  END SELECT

		else
		   call drawline("-")
		   write(*,*) 'Error in advect-diffuse particles!'
		   call drawline("-")
		endif

	    case default
	      write(*,*)'space_dim must be [0,2,3]'
	      stop
	    end select

	  end do

     if(.NOT. structured) then
      !update elements
	  call update_element(ng,g)
	 elseif(structured) then
	  !update indices
	  call update_indices(ng,g)
	 else
	  write(*,*) 'Error in update particles in advection-diffusion!'
	 endif


	end subroutine do_adv_diff

!---------------------------------------------------------
! driver to advance the biology in time 
!---------------------------------------------------------
   subroutine do_bio(ng,g,t,its)
	  integer     , intent(in)    :: ng,its
	  type(igroup), intent(inout), dimension(ng) :: g
	  real(sp), intent(in)        :: t
	  integer                     :: n

	  do n=1,ng
	    if(g(n)%biology .and. mod(its,g(n)%intvl_bio) == 0) call advance_bio(g(n),t)
	  end do
   end subroutine do_bio

!---------------------------------------------------------
! driver to interpolate forcing onto particle positions 
!---------------------------------------------------------
    subroutine interp_forcing(ng,g,iframe)

      implicit none
      integer     , intent(in)                   :: ng,iframe
      type(igroup), intent(inout), dimension(ng) :: g

      integer,  pointer   :: cell(:),istatus(:),i(:)
      real(sp), pointer   :: x(:),y(:),s(:),f(:)
      integer             :: n,vtype,v,itm
      character(len=fstr) :: evar,svar
      type(pvar),pointer  :: p

      do n=1,ng

        !skip 0-d groups
        if(g(n)%space_dim < 2)cycle

        !get the horizontal particle positions for the group
        call get_state('x',g(n),x)
        call get_state('y',g(n),y)
        call get_state('cell',g(n),cell)
        call get_state('status',g(n),istatus)

        !get the vertical position if sim is 3-d
        if(g(n)%space_dim == 3)  call get_state('s',g(n),s)

        !loop over external variables, interpolate onto points
        do v=1,g(n)%next
          !get varname and type
          svar  = g(n)%ext_var(v,1)
!          write(*,*) 'svar=',svar
          evar  = g(n)%ext_var(v,2)
!          write(*,*) 'evar=',evar
          !!zhuxm change 1 to pointer
          p=>get_pvar(g(n)%state,svar)
          vtype = p%istype  !! 1 !istype(g(n),svar) gwc debug
!          write(*,*) 'vtype=',vtype

          !get pointer to var
          if(vtype == flt_type)then
            call get_state(svar,g(n),f)
          elseif(vtype == int_type)then
            call get_state(svar,g(n),i)  !! zhuxm changed from evar to svar
          else
            write(*,*)'error in interp_forcing'
            write(*,*)'cant deal with state variable of type: ',vtype
            stop
          endif
          if(index(evar,'_ext') /= 0)then
            itm=10
          else
            itm=iframe
          endif
      !interpolate - 2D      interp in fvcom_drivers.f90
!      if(g(n)%space_dim == 2)then
          if(index(svar,'2DIM') /= 0 .or. index(svar,'2dim') /= 0 .or. g(n)%space_dim == 2)then
            if(vtype == flt_type)then
              call interp(g(n)%nind,x,y,cell,istatus,evar,f,itm)
            elseif(vtype == int_type)then
           !  call interp(g(n)%nind,x,y,cell,istatus,evar,i)
            endif
          else !3D
            if(vtype == flt_type)then
              call interp(g(n)%nind,x,y,s,cell,istatus,evar,f,itm) 
            elseif(vtype == int_type)then
           !  call interp(g(n)%nind,x,y,s,cell,istatus,evar,i)
            endif
          endif   
        end do !loop over external vars
      end do !group loop
      nullify(x,y,s,istatus,cell)
end subroutine interp_forcing

!----------------------------------------------------
! determine in which element each particle in the  
! group resides.
!----------------------------------------------------
subroutine update_element(ng,g)

  implicit none
  integer, intent(in) :: ng
  type(igroup), intent(inout) :: g(ng)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  integer , pointer :: cell(:)
  integer , pointer :: istatus(:)
  integer :: np,n

  do n=1,ng
    if(g(n)%space_dim < 2)cycle

    !set dimensinos
    np = g(n)%nind

    !set pointers to states 
    call get_state('x',g(n),x)
    call get_state('y',g(n),y)
    call get_state('cell',g(n),cell)
    call get_state('status',g(n),istatus)

    !update element containing cell 
    call find_element(np,x(1:np),y(1:np),cell(1:np),istatus(1:np))
  end do
  nullify(x,y,cell,istatus)
end subroutine update_element

End Module mod_driver
