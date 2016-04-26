Module fvcom_driver

!=======================================================================
! FVCOM Driver
!
! Description
!    - Read and setup mesh
!    - Advect
!    - Diffuse
!    - Locate Particles
!
! Comments:
!    - Advection currently uses 1st Order Euler Step, update to a stiffer
!      solver (for example 4-stage RK used by C. Chen)
!    - Vertical diffusion uses Vissers modified random walk
!    - Vertical diffusion improvement possibilities:
!         a.) splines (previously used by J. Pringle, R. Ji, and M. Huret)
!         b.) tensioned splines to avoid needing to smooth diffusivity
!             before splining
!         c.) binned random walk (Thygesen and Adlandsvik, MEPS 347,2007)
!    - In theory fiscm could be driven by an alternate ocean model if
!      these routines (initialize,advect,diffuse,find) are recoded
!
! !REVISION HISTORY:
!  Original author(s): G. Cowles
!  2011/7/28       Xueming Zhu   modified subroutines ADVECT2D,ADVECT3D
!                                interp_float2D,interp_float3D
!                                add random mixed layer as top and bottom
!                                boundary conditions, and add cubic spline
!                                interpolation for KH
!
!=======================================================================

     use gparms
     use mod_igroup
     use mod_forcing
     use utilities

    Implicit None

!Added by Xinyou Lin
     integer, allocatable :: nbe(:,:)            !!INDICES OF ELMNT NEIGHBORS
     integer, allocatable :: isonb(:)            !!NODE MARKER = 0,1,2   ^M
     integer, allocatable :: isbce(:)
     integer, allocatable :: nbvt(:,:)
!dimensions
     integer :: KB ,KBM1     !! zhuxm changed in order to keep consistent with FVCOM
     integer :: MGL,NGL
     integer :: Max_Elems
!mesh
     logical           :: mesh_setup = .false.
     real(sp), pointer :: xm(:),xm0(:)
     real(sp), pointer :: ym(:),ym0(:)
     real(sp), pointer :: xc(:),xc0(:)
     real(sp), pointer :: yc(:)   !!zhuxm ,yc0(:)
     real(sp), pointer :: hm(:)
     real(sp), pointer :: aw0(:,:)
     real(sp), pointer :: awx(:,:)
     real(sp), pointer :: awy(:,:)
     integer,  pointer :: NV(:,:)

     integer, pointer  :: ntve(:)
     integer, pointer  :: nbve(:,:)
     real(sp), pointer :: siglay(:,:)
     real(sp), pointer :: siglev(:,:)
     real(sp), pointer :: esiglay(:,:)
     real(sp), pointer :: esiglev(:,:)

!added by Xinyou
     real(sp), pointer :: a1u(:,:)
     real(sp), pointer :: a2u(:,:)
     real(sp), pointer :: chl(:,:)

     logical           :: grid_metrics

interface interp
    module procedure interp_float2D
    module procedure interp_float3D
end interface

interface interp_from_nodes
    module procedure interp_flt_from_nodes
end interface

contains

!----------------------------------------------------
! Read the mesh and interpolation coefficients
! Add mesh to output files for viz
!----------------------------------------------------
subroutine ocean_model_init(ng,g,lsize,varlist)
  use utilities, only : drawline,ncdchk
  integer, intent(in)      :: ng
  type(igroup), intent(in) :: g(ng)
  integer, intent(inout)   :: lsize
  character(len=*)         :: varlist(max_state_vars)
  !----------------------------------------------------
  character(len=mstr)      :: msg
  integer                  :: dimid,varid,fid
  character(len=fstr)      :: dname
  integer                  :: subset(3)
  integer                  :: i,n,ierr,ofid,k
  integer                  :: x_vid,y_vid,h_vid,nv_vid
  integer                  :: nele_did,node_did,three_did,zeta_vid
  integer                  :: aw0_vid,awx_vid,awy_vid,a1u_vid,a2u_vid
  integer                  :: four_did

  !return if group spatial dimensions are all 0-d or 1-d
  if(maxval(g%space_dim) < 2)return

  !get the forcing file netcdf id
  fid = get_ncfid()

  !add required time dependent variables to the list
  !--------------------------------------------------------
  ! open and read  time dependent variables namelist:  nml_ncvar
  ! commeant out by zhuxm, they are the same for each case
  !--------------------------------------------------------

!!zhuxm  open(unit=iunit,file=trim(runcontrol),form='formatted')
!!zhuxm  read(unit=iunit,nml=nml_ncvar,iostat=ios)
!!zhuxm  if(ios /= 0)then
!!zhuxm     write(*,*)'fvcom:fatal error: could not read fiscm namelist from',trim(runcontrol)
!!zhuxm     stop
!!zhuxm  endif
        nele_char ='nele'
        node_char ='node'
      siglay_char ='siglay'
      siglev_char ='siglev'
   if(spherical==1)then
           x_char ='lon'
           y_char ='lat'
   else
           x_char ='x'
           y_char ='y'
   endif
           z_char ='z'
           s_char ='s'
           h_char ='h'
          nv_char ='nv'
         a1u_char ='a1u'
         a2u_char ='a2u'
         aw0_char ='aw0'
         awx_char ='awx'
         awy_char ='awy'
        zeta_char ='zeta'
           u_char ='u'
           v_char ='v'
       omega_char ='ww'
          kh_char ='kh'
          ua_char ='ua'
          va_char ='va'
          wu_char ='uuwind'   ! 'uwind_speed'
          wv_char ='vvwind'   ! 'vwind_speed'

  do n=1,ng
      lsize = lsize + 1 ; varlist(lsize) = zeta_char
      lsize = lsize + 1 ; varlist(lsize) = h_char
    if (wind_type == 1)then
      lsize = lsize + 1 ; varlist(lsize) = wu_char
      lsize = lsize + 1 ; varlist(lsize) = wv_char
    endif
    if(g(n)%space_dim == 2)then
      lsize = lsize + 1 ; varlist(lsize) = ua_char
      lsize = lsize + 1 ; varlist(lsize) = va_char
    elseif(g(n)%space_dim ==3)then
      lsize = lsize + 1 ; varlist(lsize) = u_char
      lsize = lsize + 1 ; varlist(lsize) = v_char
      lsize = lsize + 1 ; varlist(lsize) = omega_char
      lsize = lsize + 1 ; varlist(lsize) = siglev_char   !!zhuxm
      !rji: if vdiff_type==0, then do not load kh from fvcom output
      if(g(n)%vdiff_type .gt. 0)then
        lsize = lsize + 1 ; varlist(lsize) = kh_char
      endif
    endif
  end do

  !determine number of elements
  msg = "dimension 'nele' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, nele_char, dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, NGL ))

  !determine number of nodes
  msg = "dimension 'node' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, node_char, dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, MGL ))

  !determine number of layers
  msg = "dimension 'siglay' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, siglay_char, dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, KBM1 ))
  KB = KBM1 + 1

  !allocate dataspace and read in mesh
  allocate(xm(MGL))
  msg = "error reading x coordinate"
  call ncdchk( nf90_inq_varid(fid,x_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, xm),msg)

  allocate(ym(MGL))
  msg = "error reading y coordinate"
  call ncdchk( nf90_inq_varid(fid,y_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, ym),msg)

  allocate(hm(MGL))
  msg = "error reading bathymetry h"
  call ncdchk( nf90_inq_varid(fid,h_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, hm),msg)

  allocate(NV(NGL,3))
  msg = "error reading nv"
  call ncdchk( nf90_inq_varid(fid,nv_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, NV),msg)

  allocate(siglay(MGL,KBM1))
  msg = "error reading siglay"
  call ncdchk( nf90_inq_varid(fid,siglay_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglay),msg)

  allocate(siglev(MGL,KB))
  msg = "error reading siglev"
  call ncdchk( nf90_inq_varid(fid,siglev_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglev),msg)

  allocate(a1u(NGL,4)) ; a1u   = zero
!!zhuxm  msg = "error reading a1u"
!!zhuxm  call ncdchk( nf90_inq_varid(fid,a1u_char,varid),msg )
!!zhuxm  call ncdchk(nf90_get_var(fid, varid, a1u),msg)

  allocate(a2u(NGL,4)) ; a2u   = zero
!!zhuxm  msg = "error reading a2u"
!!zhuxm  call ncdchk( nf90_inq_varid(fid,a2u_char,varid),msg )
!!zhuxm  call ncdchk(nf90_get_var(fid, varid, a2u),msg)

  allocate(aw0(NGL,3)) ; aw0 = a3rd
!!zhuxm  msg = "error reading aw0"
!!zhuxm  call ncdchk( nf90_inq_varid(fid,aw0_char,varid),msg )
!!zhuxm  call ncdchk(nf90_get_var(fid, varid, aw0),msg)

  allocate(awx(NGL,3)) ; awx = zero
!!zhuxm  msg = "error reading awx"
!!zhuxm  call ncdchk( nf90_inq_varid(fid,awx_char,varid),msg )
!!zhuxm  call ncdchk(nf90_get_var(fid, varid, awx),msg)

  allocate(awy(NGL,3)) ; awy = zero
!!zhuxm  msg = "error reading awy"
!!zhuxm  call ncdchk( nf90_inq_varid(fid,awy_char,varid),msg )
!!zhuxm  call ncdchk(nf90_get_var(fid, varid, awy),msg)

  allocate(xm0(MGL),ym0(MGL))
  allocate(xc (NGL),yc (NGL))

  !calculate cell center coordinates
  if(spherical == 1)then
     xm0 = xm
     ym0 = ym
!! transfer 0-180 --> 180-360    180-360 --> 0-180
     if(north_pole)then
       allocate(xc0(NGL)) !!zhuxm don't need it   ,yc0(NGL)
       xm=xm+180.0_sp  !! angle_pre    !! default 180
     endif
     where(xm>360.0_sp) xm=xm-360.0_sp   !!zhuxm  reverse west to east
  endif

  allocate(esiglay(NGL,KBM1),esiglev(NGL,KB))
  do i=1,NGL
    subset = NV(i,1:3)
    xc (i)=a3rd*(sum(xm(subset)))   !!zhuxm
    yc (i)=a3rd*(sum(ym(subset)))   !!zhuxm

    if(spherical == 1)then
      if(north_pole .and. (maxval(xm(subset))-minval(xm(subset)) >180.0_sp))then
!!zhuxm        xc0(i)=0.0_sp;  yc0(i)=0.0_sp
        if(xm(subset(1))>180.0_sp)then
           xc(i)=xc(i)+xm(subset(1))-360.0_sp
        else
           xc(i)=xc(i)+xm(subset(1))
        endif

        if(xm(subset(2))>180.0_sp)then
           xc(i)=xc(i)+xm(subset(2))-360.0_sp
        else
           xc(i)=xc(i)+xm(subset(2))
        endif

        if(xm(subset(3))>180.0_sp)then
           xc(i)=xc(i)+xm(subset(3))-360.0_sp
        else
           xc(i)=xc(i)+xm(subset(3))
        endif
        xc(i)  = a3rd*xc(i)
        if(xc(i)<0.0_sp) xc(i)=xc(i)+360.0_sp

        xc0(i) = xc(i)
        if(xc0(i) >= 0.0_sp .and. xc0(i) <=180.0_sp)then
           xc0(i) = xc0(i) + 180.0_sp
        elseif( xc0(i) > 180.0_sp .and. xc0(i) <=360.0_sp)  then
           xc0(i) = xc0(i) - 180.0_sp
        endif
!!zhuxm don't need it  yc0(i) = yc(i)
      else
        if(xc(i)<0.0_sp) xc(i)=xc(i)+360.0_sp
      endif
!!zhuxm    else
!!zhuxm      xc(i)  = a3rd*(sum(xm(subset)))
!!zhuxm      yc(i)  = a3rd*(sum(ym(subset)))
    endif

!!zhuxm     if(spherical == 1) then
!!zhuxm        xc0(i) = xc(i)
!!zhuxm        yc0(i) = yc(i)
!!zhuxm        if(xc0(i) >= 0.0_sp .and. xc0(i) <=180.0_sp)then
!!zhuxm           xc0(i) = xc0(i) + 180.0_sp
!!zhuxm        elseif( xc0(i) > 180.0_sp .and. xc0(i) <=360.0_sp)  then
!!zhuxm           xc0(i) = xc0(i) - 180.0_sp
!!zhuxm     endif

  !calculate cell-center siglay/siglev
    do k=1,KBM1
      esiglay(i,k) = a3rd*(sum(siglay(subset,k)))
      esiglev(i,k) = a3rd*(sum(siglev(subset,k)))
    end do
    esiglev(i,KB)  = a3rd*(sum(siglev(subset,KB)))
  end do

  ! calculate secondary connectivity (nbve/ntve)
  !determine nbve/ntve - secondary connectivity, used
  !for searching element containing point

  !----------------Node, Boundary Condition, and Control Volume-----------------------!
  allocate(NBE(0:NGL,3)) ;  NBE      = 0  !!INDICES OF ELEMENT NEIGHBORS
  allocate(NTVE(0:MGL))  ;  NTVE     = 0
  allocate(ISONB(0:MGL)) ;  ISONB    = 0  !!NODE MARKER = 0,1,2
  allocate(ISBCE(0:NGL)) ;  ISBCE    = 0
  !mark boundary elements
  grid_metrics = .false.
  if(grid_metrics == .false.)then
     CALL TRIANGLE_GRID_EDGE
     CALL SHAPE_COEF_GCN   !!zhuxm They can be commented out and read from input files
     grid_metrics = .true.
  endif

  !calculate node-based interpolation coefficients

  call drawline('-')
  write(*,*)'FVCOM mesh stats '
  call drawline('-')
  write(*,*)'Number of elements:: ',NGL
  write(*,*)'Number of nodes   :: ',MGL
  write(*,*)'Number of sig levs:: ',KB
  write(*,*)'xmin              :: ',minval(xm)
  write(*,*)'xmax              :: ',maxval(xm)
  write(*,*)'ymin              :: ',minval(ym)
  write(*,*)'ymax              :: ',maxval(ym)
  !flag that mesh is setup
  mesh_setup = .true.

!! zhuxm delete mesh.nc output, don't needed it.

end subroutine ocean_model_init

!----------------------------------------------------
! Random-Walk horizontal diffusion with constant
!    turbulent eddy diffusivity
!
! m. huret use a Gauss dist. , this probably wont
! converge to the correct continuous form of diffusion
! here we will use a uniform random walk
!----------------------------------------------------
subroutine rw_hdiff_constant(g, dT)

   use utilities, only : unitrand
   implicit none

   type(igroup), intent(inout) :: g
   real(sp), intent(in)        :: dT

   real(sp), pointer           :: x(:),y(:)
   integer,  pointer           :: istatus(:)
   integer                     :: i,np
   real(sp)                    :: tscale
   real(sp),allocatable        :: rany(:)

  !set pointers to x,y particle positions
   call get_state('x',g,x)
   call get_state('y',g,y)
   call get_state('status',g,istatus)

  !set dimensions for loops and time step
   np = g%nind

  !set diffusive time scale
   tscale = sqrt(2.0_sp*dT*g%hdiff_const_val)

!!zhuxm added call ran_from_range
  allocate(rany(np))
  !horizontal random walk
  call ran_from_range(-1.0_sp,1.0_sp,np,rany)
   do i=1,np
      if(istatus(i) < 1)cycle
!     x(i) = x(i) + normal()*tscale
!     y(i) = y(i) + normal()*tscale
      if(spherical == 0 )then
         x(i) = x(i) + rany(i)*tscale
         y(i) = y(i) + rany(i)*tscale
      elseif (spherical == 1)then
         x(i) = x(i)  + rany(i)*tscale/(tpi*cosd(y(i)) + eps)
         y(i) = y(i)  + rany(i)*tscale/tpi
      endif
  end do

  if(spherical == 1)then
     where( x < 0.0_SP)   x = x + 360.0_SP
     where( x > 360.0_SP) x = x - 360.0_SP
     where( y > 90.0_SP)  y = 180.0_SP  - y
     where( y < -90.0_SP) y = -180.0_SP - y
  endif
  !nullify pointers
  nullify(x,y,istatus)
end subroutine rw_hdiff_constant

!----------------------------------------------------
! Random-Walk horizontal diffusion with spatially
!   variable turbulent eddy diffusivity
!
! Use Visser's modified random walk to compute step
!----------------------------------------------------
subroutine rw_hdiff_variable(g, dT)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT

  stop 'error in rw_hdiff_variable, it is not yet setup'

end subroutine rw_hdiff_variable

!----------------------------------------------------
! Random-Walk vertical diffusion
!
!   - use eddy diffusivity from the model (kh)
!   - use Vissers modified random walk to compute jump
!----------------------------------------------------
subroutine rw_vdiff(g, dT, nstep)

  use utilities, only : unitrand
  implicit none

  type(igroup), intent(inout) :: g
  real(sp), intent(in)        :: dT
  integer,  intent(in)        :: nstep
  !----------------------------
  integer               :: n,p,np
  integer,  pointer     :: istatus(:),cell(:)
  real(sp), pointer     :: x(:),y(:),s(:),h(:),z(:)  !! zhuxm add z
  real(sp), allocatable :: kh(:),kh2(:),s1(:),s_shift(:)
  real(sp), allocatable :: zeta(:),depth(:),randy(:),rany(:)
  real(sp), allocatable :: dkh_dz(:),dkh_ds(:)
  real(sp), parameter   :: delta_s = 1.0E-6_sp
  real(sp)              :: deltaT,fac,dz,ztmp

  !set pointers to particle positions and status
  call get_state('status',g,istatus)
  call get_state('cell',g,cell)
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('h',g,h)
  call get_state('z',g,z)

  !set problem size and time step
  np = g%nind

  !allocate local data
  allocate(s_shift(np)); s_shift = zero
  allocate(zeta(np))   ; zeta    = zero
  allocate(depth(np))  ; depth   = zero
  allocate(kh(np))     ; kh      = zero
  allocate(kh2(np))    ; kh2     = zero
  allocate(s1(np))     ; s1      = zero
  allocate(dkh_ds(np)) ; dkh_ds  = zero
  allocate(dkh_dz(np)) ; dkh_dz  = zero
  allocate(randy(np))  ; randy   = zero
  allocate(rany(np))   ; rany    = zero

  deltaT = dT/float(nstep)
  !set constants
  fac = 2.0_sp/rvar*deltaT  ! rvar = 1.0/3.0

  call interp(np,x,y,cell,istatus,h_char,h,3,zeta_char,zeta)
  depth  = h+zeta

  ! ==> loop over substeps
  do n=1,nstep

    !--------------------------------------------------
    ! calculate d(kh)/d(s) - brute force
    !--------------------------------------------------
!!--------------------------------zhuxm 2011/7/21-----------------------
    !set derivative step (backward near free surface)
!    s1 = s+delta_s
!    call ran_from_range(0.0_sp,1.0_sp,np,randy)
!    where(s1 > 0.0_sp) s1 = -Dtop*randy/depth
    s1=s-s*0.01_sp
    where(s==0.0_sp)s1=-delta_s
!!--------------------------------zhuxm 2011/7/21-----------------------

    !evaluate kh at both locations
    call interp_kh(np,x,y,s,cell,istatus,kh,3,s1,kh2)

    !form the derivative d(kh)/d(s)
    dkh_ds = (kh2-kh)/(s1-s)
    dkh_dz = dkh_ds/depth

    !function evaluation at [z + 0.5*dkh/dz*deltat] - Visser
    s_shift = s + 0.5_sp*dkh_ds*deltaT/(depth**2)   !!zhuxm added 0.5*
!!--------------------------------zhuxm 2011/7/21-----------------------
    call ran_from_range(0.0_sp,1.0_sp,np,randy)
    where(s_shift > 0.0_sp) s_shift = -Dtop*randy/depth
    call interp_kh(np,x,y,s_shift,cell,istatus,kh,3)
!!--------------------------------zhuxm 2011/7/21-----------------------

    ! => main loop over particles
    call ran_from_range(-1.0_sp,1.0_sp,np,rany )
    call ran_from_range( 0.0_sp,1.0_sp,np,randy)
    do p=1,np
      if(istatus(p) < 1)cycle
      !update particle position using Visser modified random walk
      dz   = dkh_dz(p)*deltaT + rany(p)*sqrt(fac*kh(p)) !Visser-modified
      s(p) = s(p) + dz/depth(p)
!!--------------------------------zhuxm 2011/7/21-----------------------
      if(depth(p)<=Dtop+Dbot)then
         s(p)=-randy(p)
      else
         ztmp=s(p)*depth(p)
         if(ztmp>-Dtop)then
            s(p) = -Dtop*randy(p)/depth(p)
         elseif(ztmp<(Dbot-depth(p)))then
            s(p) = (Dbot*randy(p)-depth(p))/depth(p)
         endif
      endif
!!--------------------------------zhuxm 2011/7/21-----------------------
!      dz     = unitrand()*sqrt(2.*kh(p))                 !naive
      !set boundary conditions at free surface and bottom
!!zhuxm      s(p) = max(s(p),-(2.0_sp+s(p))) !reflect off bottom
!!      s(p) = min(s(p),0.0)         !don't pierce free surface
!!zhuxm      s(p) = min(s(p),-s(p))
    end do
    ! <= end particle loop
  end do
  ! <== end loop over substeps
  z = s*(h+zeta)  !! + zeta   !! zhuxm added
  !deallocate workspace and nullify pointers
  deallocate(s1,kh,kh2,dkh_ds,zeta,s_shift)
  nullify(x,y,s,h,istatus)
end subroutine rw_vdiff

!----------------------------------------------------
! Random-Walk vertical diffusion using splines
!   - spline vertical diffusivity to smooth profile
!   - adjust at ends
!   - use vissers modified random walk to compute jump
!
! Follow M. Huret's Implementation for the spline
!----------------------------------------------------
subroutine rw_vdiff_splined(g, dT, nstep)

  type(igroup), intent(inout) :: g
  real(sp), intent(in)        :: dT
  integer,  intent(in)        :: nstep

  real(sp)                    :: kh(KB)
  integer, pointer            :: istatus(:)
  integer                     :: n,p,np

  stop 'Error in rw_vdiff_slined, it is not yet setup'

  call get_state('status',g,istatus)
  np = g%nind

  ! => main loop over particles
  do p=1,np
    if(istatus(p) < 1)cycle

    !interpolate eddy diffusivity at model control points in column

    !spline eddy diffusivity - assume constant during step

    !loop over substeps
    do n=1,nstep

      !evaluate eddy diffusivity value and gradient at x,y,s

      !calculate random step using Visser formulation

      !use a random mixed layer, 1m from top and bottom

    end do

  end do
  ! <= end main particle loop
end subroutine rw_vdiff_splined

!----------------------------------------------------
! Random-Walk vertical diffusion using bins
!
!  - Thygesen and Adlandsvik, MEPS, v347, 2007
!----------------------------------------------------
subroutine rw_vdiff_binned(g, dT, nstep)
    type(igroup), intent(inout) :: g
    real(sp), intent(in)        :: dT
    integer,  intent(in)        :: nstep

    stop 'Error in rw_vdiff_binned, it is not yet setup!'

end subroutine rw_vdiff_binned
!----------------------------------------------------
! unfinished yet
!----------------------------------------------------

subroutine sz_trans(np,g)

  implicit none

  integer, intent(in) :: np
  type(igroup), intent(inout) :: g
  real(sp), pointer           :: x(:),y(:),z(:),s(:),h(:)
  integer , pointer           :: cell(:),istatus(:)
  real(sp), dimension(np)     :: zeta
  integer                     :: ip

  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('z',g,z)
  call get_state('s',g,s)
  call get_state('h',g,h)
  call get_state('cell',g,cell)
  call get_state('status',g,istatus)

  call interp(np,x,y,cell,istatus,h_char,h,3,zeta_char,zeta)

  if(sz_cor == 1)then
     z  = -z  !!zhuxm + zeta
     do ip=1,np
        if(h(ip)+zeta(ip) /=0.0)then
           s(ip)=z(ip)/(h(ip)+zeta(ip))
!!zhuxm           s(ip)=(z(ip)-zeta(ip))/(h(ip)+zeta(ip))
        else
           s(ip)=0.0
        endif
     enddo
  elseif(sz_cor == 0)then
     z = s*(h + zeta) !!zhuxm + zeta
  endif
  nullify(x,y,s,z,h,cell,istatus)
end subroutine sz_trans

!---------------------------------------------------
! 2-D Advection
!----------------------------------------------------
subroutine advect2D(g,deltaT,np)

  implicit none

  integer, intent(in)         :: np
  type(igroup), intent(inout) :: g
  real(sp), intent(in)        :: deltaT
  integer                     :: k,i,ns
  real(sp), pointer           :: x(:),y(:),h(:),u(:),v(:),pathlength(:)
  integer , pointer           :: cell0(:),istatus0(:)

  integer                     :: cell(np),istatus(np)
  real(sp), dimension(np)     :: u1,u2,v1,v2
  real(sp), dimension(np)     :: pdx,pdy,pdxt,pdyt
!  real(sp), dimension(np,1:mstage) :: chix,chiy
  real(sp), dimension(np,0:mstage) :: chix,chiy !rji,starting from 1 cause outofbound error


  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('h',g,h)
  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('pathlength',g,pathlength)
  call get_state('cell',g,cell0)
  call get_state('status',g,istatus0)

  !--Initialize Stage Functional Evaluations
  chix = 0.0_sp
  chiy = 0.0_sp

  istatus=0;
  where( istatus0 == 1) istatus =1
  cell=cell0
  !--Loop over RK Stages
  do ns=1,mstage

     !!Particle Position at Stage N (x,y)
     pdx(:) = x(:) + a_rk(ns)*deltaT*chix(:,ns-1)*FLOAT(istatus(:))
     pdy(:) = y(:) + a_rk(ns)*deltaT*chiy(:,ns-1)*FLOAT(istatus(:))

     if(bcond_type == 0)then
        call find_element(np,pdxt,pdyt,cell,istatus)
     elseif(bcond_type == 1)then
        call find_element_lazy_(np, pdxt, pdyt, cell, istatus)
     end if
     !!Calculate Velocity Field for Stage N Using C_RK Coefficients
     !interpolate velocity field to particle position
     call interp(np,pdx,pdy,cell,istatus,ua_char,u1,3,va_char,v1)  !!zhuxm
     call interp(np,pdx,pdy,cell,istatus,ua_char,u2,4,va_char,v2)  !!zhuxm
     chix(:,ns) = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2
     chiy(:,ns) = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2
  end do

  !--Sum Stage Contributions to get Updated Particle Positions-------------------!
  pdxt(:)  = x(:)
  pdyt(:)  = y(:)
  istatus=0;                         !! added by zhuxm
  where( istatus0 == 1) istatus =1
  cell=cell0

  do ns=1,mstage
     pdxt(:) = pdxt(:) + b_rk(ns)*deltaT*chix(:,ns)*FLOAT(istatus(:))
     pdyt(:) = pdyt(:) + b_rk(ns)*deltaT*chiy(:,ns)*FLOAT(istatus(:))
     if(bcond_type == 0)then
        call find_element(np,pdxt,pdyt,cell,istatus)
     elseif(bcond_type == 1)then
        call find_element_lazy_(np, pdxt, pdyt, cell, istatus)
     end if
  end do

  !--Update only particles that are still in the water
  if(bcond_type == 0)then ! Absorbing
     where(istatus0 == 1) istatus0 = istatus   !!zhuxm
     where(istatus /= 1) istatus = 0
     cell0=cell

     pathlength(:) = pathlength(:) + &
          sqrt((pdxt(:) - x(:))**2 + (pdyt(:) - y(:))**2) * FLOAT(istatus(:))
     x(:) = x(:)*FLOAT(1 - istatus(:)) + pdxt(:)*FLOAT(istatus(:))
     y(:) = y(:)*FLOAT(1 - istatus(:)) + pdyt(:)*FLOAT(istatus(:))
     u(:) = u(:)*FLOAT(1 - istatus(:)) + chix(:,mstage)*FLOAT(istatus(:))
     v(:) = v(:)*FLOAT(1 - istatus(:)) + chiy(:,mstage)*FLOAT(istatus(:))
  elseif(bcond_type == 1)then ! Non-absorbing (istatus0 never updated)
     do i = 1,np
        if(istatus(i) == 1)then ! Particle is still active
           cell0(i) = cell(i)
           pathlength(i) = pathlength(i) + &
                sqrt((pdxt(i) - x(i))**2 + (pdyt(i) - y(i))**2)
           x(i) = pdxt(i)
           y(i) = pdyt(i)
           u(i) = chix(i,mstage)
           v(i) = chiy(i,mstage)
        else ! Particle left the grid, don't update pathlength, x, y
           u(i) = 0
           v(i) = 0
        endif
     enddo
  endif

  !disassociate pointers
  nullify(x,y,h,u,v,cell0,istatus0)

end subroutine advect2D
!---------------------------------------------------
! 3-D Advection
!----------------------------------------------------
subroutine advect3D(g,deltaT,np,time)

  implicit none
  integer, intent(in)         :: np
  type(igroup), intent(inout) :: g
  real(sp), intent(in)        :: deltaT,time

  real(sp), pointer           :: x(:),y(:),z(:),s(:),h(:)
  real(sp), pointer           :: u(:),v(:),w(:),pathlength(:)
  integer , pointer           :: cell0(:),istatus0(:)
  integer                     :: k,i,ns,ni,istatus(np),p,cell(np)
  real(sp), dimension(np)     :: u0,u1,u2,v0,v1,v2,w0,w1,w2,wm
  real(sp), dimension(np)     :: zeta,zeta1,zeta2,pdx,pdy,pdz
  real(sp), dimension(np)     :: wu,wu1,wu2,wv,wv1,wv2
  real(sp), dimension(np)     :: pdxt,pdyt,pdzt

!  real(sp), dimension(np,1:mstage) :: chix,chiy,chiz
  real(sp), dimension(np,0:mstage) :: chix,chiy,chiz !rji outofbound error
  real(sp)                         :: diel,ztmp,randy(np),depth(np)

  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('z',g,z)
  call get_state('s',g,s)
  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('w',g,w)
  call get_state('h',g,h)
  call get_state('pathlength',g,pathlength)
  call get_state('cell',g,cell0)
  call get_state('status',g,istatus0)

  istatus=0;
  where( istatus0 == 1) istatus =1
  cell=cell0
  !disassociate pointers
  !--Initialize Stage Functional Evaluations
  zeta = 0.0_sp;  zeta1= 0.0_sp;  zeta2= 0.0_sp
  chix = 0.0_sp;  chiy = 0.0_sp;  chiz = 0.0_sp
  pdx  = x    ;   pdy  = y     ;  pdz  = s

  diel=time/3600.0-int(time/3600.0/24.0)*24.0

  !--Loop over RK Stages
  do ns=1,mstage
     !!Particle Position at Stage N (x,y,sigma)
     if(spherical == 0)then
         pdx(:) = x(:) + a_rk(ns)*deltaT*chix(:,ns-1)*FLOAT(istatus(:))
         pdy(:) = y(:) + a_rk(ns)*deltaT*chiy(:,ns-1)*FLOAT(istatus(:))
         pdz(:) = s(:) + a_rk(ns)*deltaT*chiz(:,ns-1)*FLOAT(istatus(:))
     elseif(spherical == 1)then
         pdx(:) = x(:) + a_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COSD(pdy(:)) + eps)*FLOAT(istatus(:))  !! zhuxm y(:)--> pdy(:)
         pdy(:) = y(:) + a_rk(ns)*deltaT*chiy(:,ns-1)/tpi*FLOAT(istatus(:))
         pdz(:) = s(:) + a_rk(ns)*deltaT*chiz(:,ns-1)*FLOAT(istatus(:))
         where( pdx < 0.0_SP)    pdx = pdx + 360.0_SP    !! zhuxm move them to here from outside
         where( pdx >= 360.0_SP) pdx = pdx - 360.0_SP    !! of the do loop
         where( pdy > 90.0_SP)   pdy = 180.0_SP - pdy    !! across pole
         where( pdy < -90.0_SP)  pdy = -180.0_SP - pdy
     endif
!!zhuxm changed the order for using them in top/bottom boundary---------
     if(bcond_type == 0)then
        call find_element(np,pdx,pdy,cell,istatus)     !! added by zhuxm
     elseif(bcond_type == 1)then
        call find_element_lazy_(np, pdx, pdy, cell, istatus)
     endif
     call interp(np,pdx,pdy,cell,istatus,h_char,h,3,zeta_char,zeta1)
!!--------------------------------zhuxm 2011/7/21-----------------------

     !!Adjust Sigma Position to Reflect Off Bottom (Mirroring)
!!zhuxm     pdz = max(pdz,-(2.0+pdz))

     !!Adjust Sigma Position to Remain Below Free Surface
!!zhuxm     pdz = min(pdz,-pdz)  !! zhuxm 0.0_sp)

!!--------------------------------zhuxm 2011/7/21-----------------------
      depth = h+zeta1
      call ran_from_range(0.0_sp,1.0_sp,np,randy)
      do p=1,np
         if(depth(p)<=Dtop+Dbot)then
            pdz(p)=-randy(p)
         else
            ztmp = pdz(p)*depth(p)
            if(ztmp>-Dtop)then
               pdz(p) = -Dtop*randy(p)/depth(p)
            elseif(ztmp<(Dbot-depth(p)))then
               pdz(p) = (Dbot*randy(p)-depth(p))/depth(p)
            endif
         endif
      enddo
!!--------------------------------zhuxm 2011/7/21-----------------------

     !!Calculate Velocity Field for Stage N Using C_RK Coefficients
     !interpolate velocity field to particle position


     call interp(np,pdx,pdy,pdz,cell,istatus,u_char,u1,3,v_char,v1,omega_char,w1)
     call interp(np,pdx,pdy,pdz,cell,istatus,u_char,u2,4,v_char,v2,omega_char,w2)
     call interp(np,pdx,pdy,cell,istatus,zeta_char,zeta2,4)

     u0 = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2
     v0 = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2
     w0 = (1.0_sp-c_rk(ns))*w1 + c_rk(ns)*w2
     zeta=(1.0_sp-c_rk(ns))*zeta1  + c_rk(ns)*zeta2
     if(wind_type == 1)then
        call interp(np,pdx,pdy,cell,istatus,wu_char,wu1,3,wv_char,wv1)
        call interp(np,pdx,pdy,cell,istatus,wu_char,wu2,4,wv_char,wv2)
        wu = (1.0_sp-c_rk(ns))*wu1 + c_rk(ns)*wu2
        wv = (1.0_sp-c_rk(ns))*wv1 + c_rk(ns)*wv2
        u0 = u0 + wu*0.02
        v0 = v0 + wv*0.02
     endif
!! what's the reason for setting w=0.0  ????????????????????????????????????
!Added by Xinyou Lin for DVM modelling:WP=WP+WM
     if(dvm_bio == 1.0)then
        w0=0.0                  !! zhuxm
        do ni = 1, np
           if(6.0 <= diel .and. diel <18.0)then
              pdzt(ni) = (1 + pdz(ni))*(h(ni)+zeta(ni))
!        w(ni) = w(ni) - 0.0075*tanh((pdzt(ni) -  2)*3.14159)
!         w(ni) = w(ni) - 0.0075*tanh((pdzt(ni) -  10))
              w0(ni) = w0(ni) - 0.0075*tanh((pdzt(ni) -  dvmh_dn)) !from bottom
           else
              pdzt(ni) = -pdz(ni)*(h(ni)+zeta(ni))
              w0(ni) = w0(ni) + 0.0075*tanh((pdzt(ni) - dvmh_up))   !from surface
           endif
        enddo
     endif
!///////////////////////

     chix(:,ns)  = u0(:)
     chiy(:,ns)  = v0(:)
     chiz(:,ns)  = w0(:)/(h(:)+zeta(:))  !delta_sigma/deltaT = ww/D
     !!Limit vertical motion in very shallow water
     where( (h + zeta) < eps)  chiz(:,ns) = 0.0_sp
  end do
  !--Sum Stage Contributions to get Updated Particle Positions-------------------!
  pdxt(:)  = x(:)
  pdyt(:)  = y(:)
  pdzt(:)  = s(:)
!!--------------------------------zhuxm 2011/7/21-----------------------
  istatus=0;
  where(istatus0 == 1) istatus =1
  cell=cell0
!!--------------------------------zhuxm 2011/7/21-----------------------
  do ns=1,mstage
     if(spherical == 0)then
        pdxt(:) = pdxt(:) + b_rk(ns)*deltaT*chix(:,ns)*FLOAT(istatus(:))
        pdyt(:) = pdyt(:) + b_rk(ns)*deltaT*chiy(:,ns)*FLOAT(istatus(:))
        pdzt(:) = pdzt(:) + b_rk(ns)*deltaT*chiz(:,ns)*FLOAT(istatus(:))
     elseif(spherical == 1)then
        pdxt(:) = pdxt(:) + b_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COSD(pdyt(:)) + eps)*FLOAT(istatus(:))  !!!! zhuxm y(:)--> pdyt(:)
        pdyt(:) = pdyt(:) + b_rk(ns)*deltaT*chiy(:,ns-1)/tpi*FLOAT(istatus(:))
        pdzt(:) = pdzt(:) + b_rk(ns)*deltaT*chiz(:,ns-1)*FLOAT(istatus(:))
        where( pdxt < 0.0_SP)   pdxt = pdxt + 360.0_SP   !! zhuxm moved them to here from the outside
        where( pdxt > 360.0_SP) pdxt = pdxt - 360.0_SP   !! of do loop
        where( pdyt > 90.0_SP)  pdyt =  180.0_SP - pdyt
        where( pdyt < -90.0_SP) pdyt = -180.0_SP - pdyt
     endif
!!--------------------------------zhuxm 2011/7/21-----------------------
     if(bcond_type == 0)then
        call find_element(np,pdxt,pdyt,cell,istatus)  ! Ben debugging 10/07/2014 change 'istatus0' to 'istatus'
     elseif(bcond_type == 1)then
        call find_element_lazy_(np,pdxt,pdyt,cell,istatus) ! Ben debugging 10/07/2014
     endif
     call interp(np,pdxt,pdyt,cell,istatus,h_char,h,3,zeta_char,zeta1)   !! zhuxm
     depth = h+zeta1

     call ran_from_range(0.0_sp,1.0_sp,np,randy)
!--Adjust Depth of Updated Particle Positions----------------------------------!
     do p=1,np
        if(depth(p)<=Dtop+Dbot)then
           pdzt(p)=-randy(p)
        else
           ztmp = pdzt(p)*depth(p)
           if(ztmp>-Dtop)then
              pdzt(p) = -Dtop*randy(p)/depth(p)
           elseif(ztmp<(Dbot-depth(p)))then
              pdzt(p) = (Dbot*randy(p)-depth(p))/depth(P)
           endif
        endif
!!zhuxm  s = max(s,-(2.0+s))                 !Reflect off Bottom
!!zhuxm  s = min(s,-s) !! zhuxm s=min(s,0.0) !Don t Pierce Free Surface
     enddo
!!--------------------------------zhuxm 2011/7/21-----------------------
  end do

  where(istatus0  == 1) istatus0 =istatus
  where(istatus /= 1) istatus=0
  cell0=cell
  if(bcond_type == 0)then ! Absorbing
     pathlength(:) = pathlength(:) + &
          sqrt((pdxt(:) - x(:))**2 + (pdyt(:) - y(:))**2) * FLOAT(istatus(:))
     x(:) = x(:)*FLOAT(1 - istatus(:)) + pdxt(:)*FLOAT(istatus(:))
     y(:) = y(:)*FLOAT(1 - istatus(:)) + pdyt(:)*FLOAT(istatus(:))
     s(:) = s(:)*FLOAT(1 - istatus(:)) + pdzt(:)*FLOAT(istatus(:))
     u(:) = u(:)*FLOAT(1 - istatus(:)) + chix(:,mstage)*FLOAT(istatus(:))
     v(:) = v(:)*FLOAT(1 - istatus(:)) + chiy(:,mstage)*FLOAT(istatus(:))
     w(:) = w(:)*FLOAT(1 - istatus(:)) + chiz(:,mstage)*FLOAT(istatus(:))
  elseif(bcond_type == 1)then ! Non-absorbing
     do i = 1,np
        if(istatus(i) == 1)then ! Only update particles on the grid
           pathlength(i) = pathlength(i) + &
                sqrt((pdxt(i) - x(i))**2 + (pdyt(i) - y(i))**2)
           x(i) = pdxt(i)
           y(i) = pdyt(i)
           s(i) = pdzt(i)
           u(i) = chix(i, mstage)
           v(i) = chiy(i, mstage)
           w(i) = chiz(i, mstage)
        endif
     enddo
  endif

  !--Evaluate Bathymetry and Free Surface Height at Updated Particle Position----!
!!zhuxm  call find_element(np,x,y,cell,istatus)
  call interp(np,x,y,cell,istatus0,h_char,h,4,zeta_char,zeta)            !!zhuxm

  !--Sigma adjustment if fixed depth tracking------------------------------------!
  if(fix_dep == 1)then
      s = (-zpini)/max((h+zeta),eps)  !  THIS IS REALLY PDZN = ((-LAG%ZPIN+EP) - EP)/(HP+EP)
                             !  WHERE ZPINI IS THE SPECIFIED FIXED DEPTH RELATIVE TO THE SS
      s = max(s,-1.0_SP)     ! Depth can change though if particle goes into shallower areas
  endif

  !--Calculate Particle Location in Cartesian Vertical Coordinate----------------!
  z = s*(h+zeta)  !! + zeta
  !disassociate pointers
  nullify(x,y,z,s,h,u,v,w,cell0,istatus0)
end subroutine advect3D
!----------------------------------------------------
! interpolation routine
!   interpolate scalar field [vname] to vector [v]
!   at particle position [x,y] in cell [cell]
!   if [istatus] < 0 (inactive) do not interpolate
!----------------------------------------------------
subroutine interp_float2D(np,x,y,cell,istatus,vn1,v1,iframe,vn2,v2,vn3,v3)

  implicit none
  integer, intent(in)    :: np,iframe
  real(sp),intent(in)    :: x(np),y(np)
  integer, intent(in)    :: cell(np),istatus(np)
  character(len=*)       :: vn1
  real(sp),intent(inout) :: v1(np)
  character(len=*),intent(inout),optional  :: vn2,vn3
  real(sp),intent(inout),optional          :: v2(np),v3(np)

  !-------------------------------
  real(sp), pointer      :: field1(:),field2(:),field3(:)  !pointer to FVCOM field
  integer                :: i,vts(3),icell,d1
  !added by xinyou
  integer                :: e1,e2,e3,n1,n2,n3
  real(sp)               :: xoc,yoc,dvdx,dvdy
  real(sp)               :: e0,ex,ey,tmpx
  real(sp)               :: xv(3),yv(3)
  integer                :: pts(3)

  !point to forcing variable
  if(iframe <5)then
     call get_forcing(trim(vn1),field1,iframe)
     if(present(v2)) call get_forcing(trim(vn2),field2,iframe)
     if(present(v3)) call get_forcing(trim(vn3),field3,iframe)
  else
     call get_data(trim(vn1),field1)
     if(present(v2)) call get_data(trim(vn2),field2)
     if(present(v3)) call get_data(trim(vn3),field3)
  endif
  !determine dimension of forcing variable
  d1 = size(field1,1)
  !gwc debug = 0th order
  if(d1 == NGL)then  !interpolate to element-based quantities
    do i=1,np
      icell = cell(i)
      if(istatus(i) < 1 .or. icell == 0)cycle
      e1  = nbe(icell,1)
      e2  = nbe(icell,2)
      e3  = nbe(icell,3)
!!zhuxm      if(spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)

!!zhuxm        pts(1:3) = NV(icell,1:3)
      xv = xm(NV(icell,1:3))           !!zhuxm
      if(spherical == 1 .and. north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then
!!zhuxm        yv = ym(pts)
!!zhuxm        tmpx = x(i)
!!zhuxm        xoc = tmpx - xc(icell)
!!zhuxm        yoc = y(i) - yc(icell)
          if(x(i) >= 0.0_sp .and. x(i) <180.0_sp)then
            tmpx = x(i) + 180.0_sp
          elseif(x(i) >= 180.0_sp .and. x(i)<=360.0_sp)  then
            tmpx = x(i) - 180.0_sp
          endif
          xoc = tmpx - xc0(icell)
!!zhuxm don't need it           yoc = y(i) - yc0(icell)
      endif
      dvdx = a1u(icell,1)*field1(icell)+a1u(icell,2)*field1(e1)     &
           + a1u(icell,3)*field1(e2)   +a1u(icell,4)*field1(e3)
      dvdy = a2u(icell,1)*field1(icell)+a2u(icell,2)*field1(e1)     &
           + a2u(icell,3)*field1(e2)   +a2u(icell,4)*field1(e3)
      v1(i) = field1(icell) + dvdx*xoc + dvdy*yoc

      if(present(v2))then
        dvdx = a1u(icell,1)*field2(icell)+a1u(icell,2)*field2(e1)     &
             + a1u(icell,3)*field2(e2)   +a1u(icell,4)*field2(e3)
        dvdy = a2u(icell,1)*field2(icell)+a2u(icell,2)*field2(e1)     &
             + a2u(icell,3)*field2(e2)   +a2u(icell,4)*field2(e3)
        v2(i) = field2(icell) + dvdx*xoc + dvdy*yoc
      endif

      if(present(v3))then
        dvdx = a1u(icell,1)*field3(icell)+a1u(icell,2)*field3(e1)     &
             + a1u(icell,3)*field3(e2)   +a1u(icell,4)*field3(e3)
        dvdy = a2u(icell,1)*field3(icell)+a2u(icell,2)*field3(e1)     &
             + a2u(icell,3)*field3(e2)   +a2u(icell,4)*field3(e3)
        v3(i) = field3(icell) + dvdx*xoc + dvdy*yoc
      endif

    end do
  elseif(d1 == MGL)then !vertex-based 2-D array
    do i=1,np
      icell = cell(i)
      if(istatus(i) < 1 .or. icell == 0)cycle
      n1  = NV(icell,1)
      n2  = NV(icell,2)
      n3  = NV(icell,3)
!!zhuxm      if(spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)
!!zhuxm        pts(1:3) = NV(icell,1:3)
      xv = xm(NV(icell,1:3)) !!zhuxm
      if(spherical == 1 .and. north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then
!!zhuxm        yv = ym(pts)
!!zhuxm        tmpx = x(i)
!!zhuxm        xoc = tmpx - xc(icell)
!!zhuxm        yoc = y(i) - yc(icell)
          if(x(i) >= 0.0_sp .and. x(i) <180.0_sp)then
            tmpx = x(i) + 180.0_sp
          elseif( x(i) >= 180.0_sp .and. x(i)<=360.0_sp)  then
            tmpx = x(i) - 180.0_sp
          endif
          xoc = tmpx - xc0(icell)
!!zhuxm          yoc = y(i) - yc0(icell)
      endif
      !----Linear Interpolation of Free Surface Height---------------------------------!
      e0 = aw0(icell,1)*field1(n1)+aw0(icell,2)*field1(n2)+aw0(icell,3)*field1(n3)
      ex = awx(icell,1)*field1(n1)+awx(icell,2)*field1(n2)+awx(icell,3)*field1(n3)
      ey = awy(icell,1)*field1(n1)+awy(icell,2)*field1(n2)+awy(icell,3)*field1(n3)
      v1(i) = e0 + ex*xoc + ey*yoc

      if(present(v2))then
        e0 = aw0(icell,1)*field2(n1)+aw0(icell,2)*field2(n2)+aw0(icell,3)*field2(n3)
        ex = awx(icell,1)*field2(n1)+awx(icell,2)*field2(n2)+awx(icell,3)*field2(n3)
        ey = awy(icell,1)*field2(n1)+awy(icell,2)*field2(n2)+awy(icell,3)*field2(n3)
        v2(i) = e0 + ex*xoc + ey*yoc
      endif

      if(present(v3))then
        e0 = aw0(icell,1)*field3(n1)+aw0(icell,2)*field3(n2)+aw0(icell,3)*field3(n3)
        ex = awx(icell,1)*field3(n1)+awx(icell,2)*field3(n2)+awx(icell,3)*field3(n3)
        ey = awy(icell,1)*field3(n1)+awy(icell,2)*field3(n2)+awy(icell,3)*field3(n3)
        v3(i) = e0 + ex*xoc + ey*yoc
      endif
    end do
  else
    write(*,*)'field has horizontal dimensions that is not nodes or elements'
    write(*,*)'do not know how to interpolate'
    stop
  endif
  nullify(field1)
  if(present(v2)) nullify(field2)
  if(present(v3)) nullify(field3)
end subroutine interp_float2D

!----------------------------------------------------
! interpolation routine
!   interp a 3D scalar to particle positions in 3-space
!----------------------------------------------------
subroutine interp_float3D(np,x,y,s,cell,istatus,vn1,v1,iframe,vn2,v2,vn3,v3)

  implicit none
  integer, intent(in)         :: np,iframe
  real(sp),intent(in)         :: x(np),y(np),s(np)
  integer, intent(in)         :: cell(np),istatus(np)
  character(len=*),intent(in) :: vn1
  real(sp),intent(inout)      :: v1(np)
  character(len=*),intent(inout),optional  :: vn2,vn3
  real(sp),intent(inout),optional          :: v2(np),v3(np)
  !-------------------------------

  real(sp), pointer      :: field1(:,:),field2(:,:),field3(:,:)  !pointer to FVCOM field
  integer                :: i,d1,d2,icell,k1,k2,ktype
  real(sp)               :: f1,f2,tmpx
  integer                :: vts(3)
  !added by xinyou
  integer                :: e1,e2,e3,k,n1,n2,n3
  real(sp)               :: xoc,yoc,dvdx,dvdy,ve01,ve02
  real(sp)               :: e0,ex,ey,ep01,ep02
  real(sp)               :: xv(3),yv(3)
  integer                :: pts(3)
  !point to forcing variable
  !get_forcing is in forcing.f90

  call get_forcing(trim(vn1),field1,iframe)
  if(present(v2)) call get_forcing(trim(vn2),field2,iframe)
  if(present(v3)) call get_forcing(trim(vn3),field3,iframe)

  !determine the dimensions
  d1 = size(field1,1)
  d2 = size(field1,2)

  !set vertical type (layers or levels)
  if(d2 == KB)then
    ktype = level_based
  else if(d2 == KBM1)then
    ktype = layer_based
  else
    write(*,*)'error in interp_float3D'
    write(*,*)'second dimension of field is not layers or levels'
    write(*,*)'KB==: ',d2
    stop
  endif

  ve01=0.0_sp
  ve02=0.0_sp

  !------------------------------------------------
  !interpolate to element-based quantities
  !------------------------------------------------
  if(d1 == NGL)then
     do i=1,np
      !set cell containing particle
      icell = cell(i)

      !loop if particle dead or out of domain
      if(istatus(i) < 1 .or. icell == 0)cycle
      if(s(i)<-1.0_sp)then      !! zhuxm
        write(*,*)'22222',i,'ssss=',s(i),np
        pause
      endif


      !determine the vertical layer brackets and interp coeffs
      call get_vert_interpcoefs(icell,s(i),k1,k2,f1,f2,ktype)

      !xinyou
      e1 = nbe(icell,1)
      e2 = nbe(icell,2)
      e3 = nbe(icell,3)
!!zhuxm      if(spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)

!!zhuxm        pts(1:3) = NV(icell,1:3)
      xv = xm(NV(icell,1:3)) !!zhuxm
      if(spherical == 1 .and. north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then
!!zhuxm         pts(1:3) = NV(icell,1:3)
!!zhuxm         xv = xm(pts)
!!zhuxm         yv = ym(pts)
!!zhuxm         tmpx = x(i)
!!zhuxm         xoc = tmpx - xc(icell)
!!zhuxm         yoc = y(i) - yc(icell)
         if(x(i) >= 0.0_sp .and. x(i) <180.0_sp)then
            tmpx = x(i) + 180.0_sp
         elseif( x(i) >= 180.0_sp .and. x(i)<=360.0_sp)  then
            tmpx = x(i) - 180.0_sp
         endif
         xoc = tmpx - xc0(icell)
!!zhuxm don't need it           yoc = y(i) - yc0(icell)
      endif
      k = k1
      dvdx = a1u(icell,1)*field1(icell,k)+a1u(icell,2)*field1(e1,k)     &
           + a1u(icell,3)*field1(e2,k)   +a1u(icell,4)*field1(e3,k)
      dvdy = a2u(icell,1)*field1(icell,k)+a2u(icell,2)*field1(e1,k)     &
           + a2u(icell,3)*field1(e2,k)   +a2u(icell,4)*field1(e3,k)
      ve01 = field1(icell,k) + dvdx*xoc + dvdy*yoc

      if(f2 /= zero)then
        k = k2
        dvdx = a1u(icell,1)*field1(icell,k)+a1u(icell,2)*field1(e1,k)     &
             + a1u(icell,3)*field1(e2,k)   +a1u(icell,4)*field1(e3,k)
        dvdy = a2u(icell,1)*field1(icell,k)+a2u(icell,2)*field1(e1,k)     &
             + a2u(icell,3)*field1(e2,k)   +a2u(icell,4)*field1(e3,k)
        ve02 = field1(icell,k) + dvdx*xoc + dvdy*yoc
      endif

      v1(i) = f1*ve01 + f2*ve02

      if(present(v2))then
         k = k1
         dvdx = a1u(icell,1)*field2(icell,k)+a1u(icell,2)*field2(e1,k)     &
              + a1u(icell,3)*field2(e2,k)   +a1u(icell,4)*field2(e3,k)
         dvdy = a2u(icell,1)*field2(icell,k)+a2u(icell,2)*field2(e1,k)     &
              + a2u(icell,3)*field2(e2,k)   +a2u(icell,4)*field2(e3,k)
         ve01 = field2(icell,k) + dvdx*xoc + dvdy*yoc

         if(f2 /= zero)then
           k = k2
           dvdx = a1u(icell,1)*field2(icell,k)+a1u(icell,2)*field2(e1,k)     &
                + a1u(icell,3)*field2(e2,k)   +a1u(icell,4)*field2(e3,k)
           dvdy = a2u(icell,1)*field2(icell,k)+a2u(icell,2)*field2(e1,k)     &
                + a2u(icell,3)*field2(e2,k)   +a2u(icell,4)*field2(e3,k)
           ve02 = field2(icell,k) + dvdx*xoc + dvdy*yoc
         endif
         v2(i) = f1*ve01 + f2*ve02
      endif

      if(present(v3))then
         k = k1
         dvdx = a1u(icell,1)*field3(icell,k)+a1u(icell,2)*field3(e1,k)     &
              + a1u(icell,3)*field3(e2,k)   +a1u(icell,4)*field3(e3,k)
         dvdy = a2u(icell,1)*field3(icell,k)+a2u(icell,2)*field3(e1,k)     &
              + a2u(icell,3)*field3(e2,k)   +a2u(icell,4)*field3(e3,k)
         ve01 = field3(icell,k) + dvdx*xoc + dvdy*yoc

         if(f2 /= zero)then
           k = k2
           dvdx = a1u(icell,1)*field3(icell,k)+a1u(icell,2)*field3(e1,k)     &
                + a1u(icell,3)*field3(e2,k)   +a1u(icell,4)*field3(e3,k)
           dvdy = a2u(icell,1)*field3(icell,k)+a2u(icell,2)*field3(e1,k)     &
                + a2u(icell,3)*field3(e2,k)   +a2u(icell,4)*field3(e3,k)
           ve02 = field3(icell,k) + dvdx*xoc + dvdy*yoc
         endif
         v3(i) = f1*ve01 + f2*ve02
      endif
    end do
  !------------------------------------------------
  !interpolate to node-based quantities
  !------------------------------------------------
  elseif(d1 == MGL)then
    do i=1,np
      !set cell containing particle
      icell = cell(i)

      !loop if particle dead or out of domain
      if(istatus(i) < 1 .or. icell == 0)cycle
      if(s(i)<-1.0_sp)then      !! zhuxm
        write(*,*)'22222',i,'ssss=',s(i),np
        pause
      endif

      !determine the vertical layer brackets and interp coeffs
      call get_vert_interpcoefs(icell,s(i),k1,k2,f1,f2,ktype)

      n1 = NV(icell,1)
      n2 = NV(icell,2)
      n3 = NV(icell,3)
!!zhuxm      if(spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)
!!zhuxm        pts(1:3) = NV(icell,1:3)
      xv = xm(NV(icell,1:3)) !!zhuxm
      if(spherical == 1 .and. north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then
!!zhuxm         pts(1:3) = NV(icell,1:3)
!!zhuxm         xv   = xm(pts)
!!zhuxm         yv   = ym(pts)
!!zhuxm         tmpx = x(i)
!!zhuxm         xoc  = tmpx - xc(icell)
!!zhuxm         yoc  = y(i) - yc(icell)
         if(x(i) >= 0.0_sp .and. x(i) <180.0_sp)then
            tmpx = x(i) + 180.0_sp
         elseif( x(i) >= 180.0_sp .and. x(i)<=360.0_sp)  then
            tmpx = x(i) - 180.0_sp
         endif
         xoc = tmpx - xc0(icell)
!!zhuxm            yoc = y(i) - yc0(icell)
      endif

      k = k1
      e0 = aw0(icell,1)*field1(n1,k)+aw0(icell,2)*field1(n2,k)+aw0(icell,3)*field1(n3,k)
      ex = awx(icell,1)*field1(n1,k)+awx(icell,2)*field1(n2,k)+awx(icell,3)*field1(n3,k)
      ey = awy(icell,1)*field1(n1,k)+awy(icell,2)*field1(n2,k)+awy(icell,3)*field1(n3,k)
      ep01 = e0 + ex*xoc + ey*yoc
      if(f2 /= zero)then
        k = k2
        e0 = aw0(icell,1)*field1(n1,k)+aw0(icell,2)*field1(n2,k)+aw0(icell,3)*field1(n3,k)
        ex = awx(icell,1)*field1(n1,k)+awx(icell,2)*field1(n2,k)+awx(icell,3)*field1(n3,k)
        ey = awy(icell,1)*field1(n1,k)+awy(icell,2)*field1(n2,k)+awy(icell,3)*field1(n3,k)
        ep02 = e0 + ex*xoc + ey*yoc
      endif
      v1(i) = f1*ep01 + f2*ep02

      if(present(v2))then
        k = k1
        e0 = aw0(icell,1)*field2(n1,k)+aw0(icell,2)*field2(n2,k)+aw0(icell,3)*field2(n3,k)
        ex = awx(icell,1)*field2(n1,k)+awx(icell,2)*field2(n2,k)+awx(icell,3)*field2(n3,k)
        ey = awy(icell,1)*field2(n1,k)+awy(icell,2)*field2(n2,k)+awy(icell,3)*field2(n3,k)
        ep01 = e0 + ex*xoc + ey*yoc
        if(f2 /= zero)then
          k = k2
          e0 = aw0(icell,1)*field2(n1,k)+aw0(icell,2)*field2(n2,k)+aw0(icell,3)*field2(n3,k)
          ex = awx(icell,1)*field2(n1,k)+awx(icell,2)*field2(n2,k)+awx(icell,3)*field2(n3,k)
          ey = awy(icell,1)*field2(n1,k)+awy(icell,2)*field2(n2,k)+awy(icell,3)*field2(n3,k)
          ep02 = e0 + ex*xoc + ey*yoc
        endif
        v2(i) = f1*ep01 + f2*ep02
      endif

      if(present(v3))then
        k = k1
        e0 = aw0(icell,1)*field3(n1,k)+aw0(icell,2)*field3(n2,k)+aw0(icell,3)*field3(n3,k)
        ex = awx(icell,1)*field3(n1,k)+awx(icell,2)*field3(n2,k)+awx(icell,3)*field3(n3,k)
        ey = awy(icell,1)*field3(n1,k)+awy(icell,2)*field3(n2,k)+awy(icell,3)*field3(n3,k)
        ep01 = e0 + ex*xoc + ey*yoc
        if(f2 /= zero)then
          k = k2
          e0 = aw0(icell,1)*field3(n1,k)+aw0(icell,2)*field3(n2,k)+aw0(icell,3)*field3(n3,k)
          ex = awx(icell,1)*field3(n1,k)+awx(icell,2)*field3(n2,k)+awx(icell,3)*field3(n3,k)
          ey = awy(icell,1)*field3(n1,k)+awy(icell,2)*field3(n2,k)+awy(icell,3)*field3(n3,k)
          ep02 = e0 + ex*xoc + ey*yoc
        endif
        v3(i) = f1*ep01 + f2*ep02
      endif
    end do
  else
    write(*,*)'error in interp_float3D'
    write(*,*)'field has horizontal dimensions that is not nodes or elements'
    write(*,*)'do not know how to interpolate'
    stop
  endif
  nullify(field1)
  if(present(v2)) nullify(field2)
  if(present(v3)) nullify(field3)
end subroutine interp_float3D

!----------------------------------------------------
! interpolation routine
!   modified by zhuxm 2011.7.25 from interp3D
!   interp 3D kh variable to particle positions in 3-space
!----------------------------------------------------
subroutine interp_kh(np,x,y,s,cell,istatus,v,iframe,s1,v1)

  implicit none
  integer, intent(in)             :: np,iframe
  real(sp),intent(in)             :: x(np),y(np),s(np)
  integer, intent(in)             :: cell(np),istatus(np)
  real(sp),intent(inout)          :: v(np)
  real(sp),intent(in   ),optional :: s1(sp)
  real(sp),intent(inout),optional :: v1(sp)
  !-------------------------------

  real(sp), pointer      :: kh(:,:),sl(:,:)  !pointer to FVCOM field
  integer                :: i,icell,k1,k2,ktype
  real(sp)               :: f1,f2,tmpx
  integer                :: vts(3)
  !added by xinyou
  integer                :: e1,e2,e3,k,n1,n2,n3
  real(sp)               :: xoc,yoc,dvdx,dvdy,ve01,ve02
  real(sp),allocatable   :: e0(:),ex(:),ey(:)
  real(sp),allocatable   :: khh(:),sig(:),ysp(:)
  real(sp)               :: xv(3),yv(3)
  integer                :: pts(3)
  !point to forcing variable
  !get_forcing is in forcing.f90

  call get_forcing(trim(kh_char)    ,kh,iframe)
  call get_forcing(trim(siglev_char),sl,iframe)

  !set vertical type (layers or levels)
  allocate(e0(KB),ex(KB),ey(KB))
  allocate(khh(KB),sig(KB),ysp(KB))
  !------------------------------------------------
  !interpolate to node-based quantities
  !------------------------------------------------
  do i=1,np
     !set cell containing particle
     icell = cell(i)
     !loop if particle dead or out of domain
     if(istatus(i) < 1 .or. icell == 0)cycle
     if(s(i)<-1.0_sp)then      !! zhuxm
        write(*,*)'22222',i,'ssss=',s(i),np
        pause
     endif
     n1  = NV(icell,1)
     n2  = NV(icell,2)
     n3  = NV(icell,3)
!!zhuxm     if(spherical == 0)then
     xoc = x(i) - xc(icell)
     yoc = y(i) - yc(icell)
!!zhuxm        pts(1:3) = NV(icell,1:3)
     xv = xm(NV(icell,1:3)) !!zhuxm
     if(spherical == 1 .and. north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then
!!zhuxm        pts(1:3) = NV(icell,1:3)
!!zhuxm        xv   = xm(pts)
!!zhuxm        yv   = ym(pts)
!!zhuxm        tmpx = x(i)
!!zhuxm        xoc  = tmpx - xc(icell)
!!zhuxm        yoc  = y(i) - yc(icell)
        if(x(i) >= 0.0_sp .and. x(i) <180.0_sp)then
           tmpx = x(i) + 180.0_sp
        elseif( x(i) >= 180.0_sp .and. x(i)<=360.0_sp)  then
           tmpx = x(i) - 180.0_sp
        endif
        xoc = tmpx - xc0(icell)
!!zhuxm        yoc = y(i) - yc0(icell)
     endif
     e0(:) = aw0(icell,1)*kh(n1,:)+aw0(icell,2)*kh(n2,:)+aw0(icell,3)*kh(n3,:)
     ex(:) = awx(icell,1)*kh(n1,:)+awx(icell,2)*kh(n2,:)+awx(icell,3)*kh(n3,:)
     ey(:) = awy(icell,1)*kh(n1,:)+awy(icell,2)*kh(n2,:)+awy(icell,3)*kh(n3,:)
     khh(:) = e0 + ex*xoc + ey*yoc

     e0(:) = aw0(icell,1)*sl(n1,:)+aw0(icell,2)*sl(n2,:)+aw0(icell,3)*sl(n3,:)
     ex(:) = awx(icell,1)*sl(n1,:)+awx(icell,2)*sl(n2,:)+awx(icell,3)*sl(n3,:)
     ey(:) = awy(icell,1)*sl(n1,:)+awy(icell,2)*sl(n2,:)+awy(icell,3)*sl(n3,:)
     sig(:) = e0 + ex*xoc + ey*yoc

     if(present(s1))then
        call splint(KB,sig,khh,s(i),v(i),s1(i),v1(i))
     else
        call splint(KB,sig,khh,s(i),v(i))
     endif
  end do
  nullify(kh,sl)
end subroutine interp_kh

!----------------------------------------------------
! element search routine
!----------------------------------------------------
subroutine find_element(np,x,y,cell,istatus,option)

  integer, intent(in)    :: np
  real(sp),intent(in)    :: x(np),y(np)
  integer, intent(inout) :: cell(np),istatus(np)
  character(len=*), optional :: option
  integer :: p

  do p=1,np
    if(istatus(p) <= 0)cycle
    !try a quick find
    cell(p) = find_element_lazy(x(p),y(p),cell(p))
    if(cell(p) /= 0)cycle
    !failed, try a robust find
    cell(p) = find_element_robust(x(p),y(p))
    if(cell(p) /= 0)cycle
    !failed, update status to lost
    istatus(p) = EXITED
  end do

end subroutine find_element

!----------------------------------------------------
! lazy element search routine
! Note that this is different from the function
! find_element_lazy that is used by find_elemnt. This
! subroutine searches all elements, whereas the
! the function searches only one.
!
! Only the lazy search is used by the non-stick boundary
! conditions
!----------------------------------------------------
subroutine find_element_lazy_(np,x,y,cell,istatus,option)

  integer, intent(in)    :: np
  real(sp),intent(in)    :: x(np),y(np)
  integer, intent(inout) :: cell(np),istatus(np)
  character(len=*), optional :: option
  integer :: p

  do p=1,np
    if(istatus(p) <= 0)cycle
    !try a quick find
    cell(p) = find_element_lazy(x(p),y(p),cell(p))
    if(cell(p) /= 0)cycle
    !failed, update status to lost
    istatus(p) = EXITED
  end do

end subroutine find_element_lazy_

!----------------------------------------------------
! find the element in which a point resides: lazy
!    - look in last known element
!    - look in last known elements neighbors
!    - should be reasonably robust
!----------------------------------------------------
function find_element_lazy(xp,yp,last) result(elem)
   implicit none
   real(sp), intent(in) :: xp
   real(sp), intent(in) :: yp
   integer,  intent(in) :: last
   !---------------------------
   integer  :: elem,j,k,test,iney,shi,pts(3)
   real(sp) :: xv(3),yv(3),xpn,ypn

   xpn=xp; ypn=yp; shi=1
   !make sure we have a mesh
   if(.not. mesh_setup)then
      write(*,*)'error in find_element_lazy:'
      write(*,*)'mesh is not setup yet'
      stop
   endif

   !initialize
   elem = 0 ; if(last == 0)return

   !load vertices
   pts(1:3) = NV(last,1:3)

   if(spherical == 1)then
      if(shi==1)then
         xv = xm(pts)
         yv = ym(pts)
      elseif(shi==0)then
         xv = xm0(pts)
         yv = ym0(pts)
      endif
      if(north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then      !!zhuxm
         xv = xm0(pts)
         if(xp >= 0.0_sp .and. xp <180.0_sp)then
            xpn = xp + 180.0_sp
         elseif( xp >= 180.0_sp .and. xp<=360.0_sp)  then
            xpn = xp - 180.0_sp
         endif
         shi=0
      endif
   elseif (spherical == 0)then
      xv = xm(pts)
      yv = ym(pts)
   endif

   if(isintriangle(xpn,ypn,xv,yv))then
      elem = last
   else
      if(grid_metrics)then
         outer: do j=1,3
                   test = NV(last,j)
                   do k=1,ntve(test)
                      iney = nbve(test,k)
                      pts(1:3) = NV(iney,1:3)
                      if(spherical == 1)then
                         if(shi==1)then
                            xv = xm(pts)
                            yv = ym(pts)
                         elseif(shi==0)then
                            xv = xm0(pts)
                            yv = ym0(pts)
                         endif
                         if(north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then   !! zhuxm
                            xv = xm0(pts)
                            if(xp >= 0.0_sp .and. xp <180.0_sp)then
                               xpn = xp + 180.0_sp
                            elseif( xp >= 180.0_sp .and. xp<=360.0_sp)then
                               xpn = xp - 180.0_sp
                            endif
                         endif
                      elseif(spherical == 0)then
                         xv = xm(pts)
                         yv = ym(pts)
                      endif
                      if(isintriangle(xpn,ypn,xv,yv))then
                         elem  = iney
                         exit outer
                      end if
                   end do
                end do outer
      end if !grid_metrics on
   endif
!=== check points around 0(365)degree
   return
end function find_element_lazy

!----------------------------------------------------
! find the element in which a point resides: robust
!    - search outward in increasing radius
!    - search only the closest max_check points
!----------------------------------------------------
function find_element_robust(xp,yp) result(elem)
   use utilities
   implicit none
   real(sp), intent(in) :: xp,yp
   integer              :: elem
   !---------------------------
   real(sp)             :: radlist(NGL),xv(3),yv(3)
   real(sp)             :: radlast,xpn,ypn
   integer              :: pts(3),min_loc,locij(1),cnt,shi
   integer, parameter   :: max_check = 15

   !make sure we have a mesh
   if(.not. mesh_setup)then
      stop 'error in find_element_robust mesh is not setup yet'
   endif
   !initialize
   xpn  = xp;  ypn = yp;  shi=1
   elem = 0 ;  cnt = 0
   if(spherical == 1)then
      radlist(1:NGL) = xc(1:NGL)-xp
      where(north_pole .and. radlist >  180.0_sp)  radlist=radlist -360.0_sp
      where(north_pole .and. radlist < -180.0_sp)  radlist=radlist +360.0_sp
      radlist(1:NGL) = sqrt((radlist(1:NGL)*cosd(yp))**2 + (yc(1:NGL)-yp)**2)
   else
      radlist(1:NGL) = sqrt((xc(1:NGL)-xp)**2 + (yc(1:NGL)-yp)**2)
   endif
   radlast = -1.0_sp
   l1: do while(cnt < max_check)
          cnt = cnt+1
          locij   = minloc(radlist,radlist>radlast)
          min_loc = locij(1)
          if(min_loc == 0) exit l1

          pts(1:3) = NV(min_loc,1:3)
          if(spherical == 1)then
             if(shi==1)then
                xv = xm(pts)
                yv = ym(pts)
             elseif(shi==0)then
                xv = xm0(pts)
                yv = ym0(pts)
             endif
             if(north_pole .and. maxval(xv)-minval(xv) >180.0_sp)then
                xv = xm0(pts)
                if(xp >= 0.0_sp .and. xp <180.0_sp)then
                   xpn = xp + 180.0_sp
                elseif( xp >= 180.0_sp .and. xp<=360.0_sp)  then
                   xpn = xp - 180.0_sp
                endif
                shi=0
             endif
          elseif (spherical == 0)then
             xv = xm(pts)
             yv = ym(pts)
          endif
          if(isintriangle(xpn,ypn,xv,yv))then
             elem = min_loc
             exit l1
          end if
          radlast = radlist(min_loc)
       end do l1
   return
end function find_element_robust

function interp_flt_from_nodes(cell,x,y,fin) result(fout)
   integer,  intent(in) :: cell
   real(sp), intent(in) :: x
   real(sp), intent(in) :: y
   real(sp), intent(in) :: fin(3)
   real(sp) :: fout
   real(sp) :: f0,fx,fy,deltaX,deltaY

   f0 = aw0(cell,1)*fin(1) + aw0(cell,2)*fin(2) + aw0(cell,3)*fin(3)
   fx = awx(cell,1)*fin(1) + awx(cell,2)*fin(2) + awx(cell,3)*fin(3)
   fy = awy(cell,1)*fin(1) + awy(cell,2)*fin(2) + awy(cell,3)*fin(3)
   deltaX = x-xc(cell)
   deltaY = y-yc(cell)
   fout = f0 + fx*deltaX + fy*deltaY

end function interp_flt_from_nodes

subroutine get_vert_interpcoefs(cell,sloc,k1,k2,f1,f2,ktype)

   integer,  intent(in)  :: ktype,cell
   real(sp), intent(in)  :: sloc
   integer,  intent(out) :: k1,k2
   real(sp), intent(out) :: f1,f2
   !---------------------------
   real(sp)              :: ds,my_sloc
   integer               :: i,k

   my_sloc = max(sloc,-one)
   my_sloc = min(sloc,-tinynum)
   !level-based data
   if(ktype == level_based)then
     do k=1,KB
       if(my_sloc >= esiglev(cell,k))exit
     end do
     if(k==1)then
       k1 = 1    ;    k2 = 1
       f1 = one  ;    f2 = zero
     elseif(k >= KB)then
       k1 = KB   ;    k2 = KB
       f1 = one  ;    f2 = zero
     else
       k1 = k-1  ;    k2 = k
       ds = esiglev(cell,k1)-esiglev(cell,k2)
       f1 = (my_sloc-esiglev(cell,k2))/ds
       f2 = one-f1
     endif
   !layer-based data
   elseif(ktype == layer_based)then
     do k=1,KBM1
       if(my_sloc >= esiglay(cell,k))exit
     end do
     if(k==1)then
       k1 = 1    ;    k2 = 1
       f1 = one  ;    f2 = zero
     elseif(k >= KBM1)then
       k1 = KBM1 ;    k2 = KBM1
       f1 = one  ;    f2 = zero
     else
       k1 = k-1  ;    k2 = k
       ds = esiglay(cell,k1)-esiglay(cell,k2)
       f1 = (my_sloc-esiglay(cell,k2))/ds
       f2 = one-f1
     endif
   !error
   else
     write(*,*)'error in get_vert_interpcoefs'
     write(*,*)'ktype must be :',level_based,' or ',layer_based
     write(*,*)'but is: ',ktype
     stop
   endif
 end subroutine get_vert_interpcoefs

 subroutine   TRIANGLE_GRID_EDGE
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TEMP,NB_TMP,CELLS,NBET
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: CELLCNT
  INTEGER                              :: I,J,II,JJ,NTMP,NCNT,NFLAG,JJB
  INTEGER                              :: N1,N2,N3,J1,J2,J3,MX_NBR_ELEM

  !----------------------------INITIALIZE----------------------------------------!
  ISBCE = 0
  ISONB = 0
  NBE   = 0
  !----DETERMINE NBE(i=1:n,j=1:3): INDEX OF 1 to 3 NEIGHBORING ELEMENTS----------!

  ALLOCATE(NBET(NGL,3)) ; NBET = 0
  ALLOCATE(CELLS(MGL,50)) ; CELLS = 0
  ALLOCATE(CELLCNT(MGL))  ; CELLCNT = 0
  DO I=1,NGL
     N1 = NV(I,1) ; CELLCNT(N1) = CELLCNT(N1)+1
     N2 = NV(I,2) ; CELLCNT(N2) = CELLCNT(N2)+1
     N3 = NV(I,3) ; CELLCNT(N3) = CELLCNT(N3)+1
     CELLS(NV(I,1),CELLCNT(N1)) = I
     CELLS(NV(I,2),CELLCNT(N2)) = I
     CELLS(NV(I,3),CELLCNT(N3)) = I
  END DO
  DO I=1,NGL
     N1 = NV(I,1)
     N2 = NV(I,2)
     N3 = NV(I,3)
     DO J1 = 1,CELLCNT(N1)
        DO J2 = 1,CELLCNT(N2)
           IF((CELLS(N1,J1) == CELLS(N2,J2)).AND. CELLS(N1,J1) /= I)NBE(I,3) = CELLS(N1,J1)
        END DO
     END DO
     DO J2 = 1,CELLCNT(N2)
        DO J3 = 1,CELLCNT(N3)
           IF((CELLS(N2,J2) == CELLS(N3,J3)).AND. CELLS(N2,J2) /= I)NBE(I,1) = CELLS(N2,J2)
        END DO
     END DO
     DO J1 = 1,CELLCNT(N1)
        DO J3 = 1,CELLCNT(N3)
           IF((CELLS(N1,J1) == CELLS(N3,J3)).AND. CELLS(N1,J1) /= I)NBE(I,2) = CELLS(N3,J3)
        END DO
     END DO
  END DO
 DEALLOCATE(CELLS,CELLCNT)
  !   IF(MSR)WRITE(IPT,*)  '!  NEIGHBOR FINDING      :    COMPLETE'
  !
  !--ENSURE ALL ELEMENTS HAVE AT LEAST ONE NEIGHBOR------------------------------!
  !
  NFLAG = 0
  DO I=1,NGL
     IF(SUM(NBE(I,1:3))==0)THEN
        NFLAG = 1
        WRITE(*,*)'ELEMENT ',I,' AT ',XC(I),YC(I),' HAS NO NEIGHBORS'
        STOP
     END IF
  END DO
  IF(NFLAG == 1) STOP
  !
  !----IF ELEMENT ON BOUNDARY SET ISBCE(I)=1 AND ISONB(J)=1 FOR BOUNDARY NODES J-!
  !
  DO I=1,NGL
     IF(MIN(NBE(I,1),NBE(I,2),NBE(I,3))==0)THEN    !!ELEMENT ON BOUNDARY
        ISBCE(I) = 1
        IF(NBE(I,1) == 0)THEN
           ISONB(NV(I,2)) = 1 ; ISONB(NV(I,3)) = 1
        END IF
        IF(NBE(I,2) ==0) THEN
           ISONB(NV(I,1)) = 1 ; ISONB(NV(I,3)) = 1
        END IF
        IF(NBE(I,3) ==0) THEN
           ISONB(NV(I,1)) = 1 ; ISONB(NV(I,2)) = 1
        END IF
     END IF
  END DO

  !==============================================================================|
  !             DEFINE NTVE, NBVE, NBVT                                          !
  !                                                                              !
  ! ntve(1:m):           total number of the surrounding triangles               !
  !                      connected to the given node                             !
  ! nbve(1:m, 1:ntve+1): the identification number of surrounding                !
  !                      triangles with a common node (counted clockwise)        !
  ! nbvt(1:m,ntve(1:m)): the idenfication number of a given node over            !
  !                      each individual surrounding triangle(counted            !
  !                      clockwise)                                              !
  !==============================================================================|

  !
  !----DETERMINE MAX NUMBER OF SURROUNDING ELEMENTS------------------------------!
  !
  MX_NBR_ELEM = 0
  DO I=1,MGL
     NCNT = 0
     DO J=1,NGL
        IF( FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0_SP) &
             NCNT = NCNT + 1
     END DO
     MX_NBR_ELEM = MAX(MX_NBR_ELEM,NCNT)
  END DO

  !
  !----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
  !
  ALLOCATE(NBVE(MGL,MX_NBR_ELEM+1))
  ALLOCATE(NBVT(MGL,MX_NBR_ELEM+1))
  !
  !--DETERMINE NUMBER OF SURROUNDING ELEMENTS FOR NODE I = NTVE(I)---------------!
  !--DETERMINE NBVE - INDICES OF NEIGHBORING ELEMENTS OF NODE I------------------!
  !--DETERMINE NBVT - INDEX (1,2, or 3) OF NODE I IN NEIGHBORING ELEMENT---------!
  !
  DO I=1,MGL
     NCNT=0
     DO J=1,NGL
        IF (FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0_SP)THEN
           NCNT = NCNT+1
           NBVE(I,NCNT)=J
           IF((NV(J,1)-I) == 0) NBVT(I,NCNT)=1
           IF((NV(J,2)-I) == 0) NBVT(I,NCNT)=2
           IF((NV(J,3)-I) == 0) NBVT(I,NCNT)=3
        END IF
     ENDDO
     NTVE(I)=NCNT
  ENDDO

  ALLOCATE(NB_TMP(MGL,MX_NBR_ELEM+1))
  DO I=1,MGL
     IF(ISONB(I) == 0) THEN
        NB_TMP(1,1)=NBVE(I,1)
        NB_TMP(1,2)=NBVT(I,1)
        DO J=2,NTVE(I)+1
           II=NB_TMP(J-1,1)
           JJ=NB_TMP(J-1,2)
           NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
           JJ=NB_TMP(J,1)
           IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
           IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
           IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
        ENDDO

        DO J=2,NTVE(I)+1
           NBVE(I,J)=NB_TMP(J,1)
        ENDDO

        DO J=2,NTVE(I)+1
           NBVT(I,J)=NB_TMP(J,2)
        ENDDO

        NTMP=NTVE(I)+1
        IF(NBVE(I,1) /= NBVE(I,NTMP)) THEN
           PRINT*, I,'NBVE(I) NOT CORRECT!!'
           STOP
        ENDIF
        IF(NBVT(I,1) /= NBVT(I,NTMP)) THEN
           PRINT*, I,'NBVT(I) NOT CORRECT!!'
           STOP
        END IF
     ELSE
        JJB=0

        DO J=1,NTVE(I)
           JJ=NBVT(I,J)
           IF(NBE(NBVE(I,J),JJ+2-INT((JJ+2)/4)*3) == 0) THEN
              JJB=JJB+1
              NB_TMP(JJB,1)=NBVE(I,J)
              NB_TMP(JJB,2)=NBVT(I,J)
           END IF
        ENDDO

        IF(JJB /= 1) THEN
           PRINT*, 'ERROR IN ISONB !,I,J', I,J
!           PAUSE
        END IF

        DO J=2,NTVE(I)
           II=NB_TMP(J-1,1)
           JJ=NB_TMP(J-1,2)
           NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
           JJ=NB_TMP(J,1)
           IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
           IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
           IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
        ENDDO

        DO J=1,NTVE(I)
           NBVE(I,J)=NB_TMP(J,1)
           NBVT(I,J)=NB_TMP(J,2)
        ENDDO
        NBVE(I,NTVE(I)+1)=0

     END IF
  END DO
  DEALLOCATE(NB_TMP)

  RETURN


 end subroutine TRIANGLE_GRID_EDGE

   SUBROUTINE SHAPE_COEF_GCN
   IMPLICIT NONE
   REAL(sp) X1,X2,X3,Y1,Y2,Y3,DELT,AI1,AI2,AI3,BI1,BI2,BI3,CI1,CI2,CI3
   REAL(sp) DELTX,DELTY,TEMP1,ANG1,ANG2,B1,B2,ANGLE,ART
   INTEGER  I,II,J,JJ,J1,J2
! if defined (SPHERICAL)
   REAL(sp) XXC,YYC,XXC1,YYC1,XXC2,YYC2,XXC3,YYC3,SIDE,&
            TY1,TY2,X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE1,SIDE2,SIDE3
   REAL(sp) XTMP1,XTMP2,XTMP3
   REAL(sp) XTMP11,XTMP21,XTMP31
!!!!!!!!!!!!

   DO I=1,NGL
     IF(ISBCE(I) == 0)THEN
       Y1 = YC(NBE(I,1))-YC(I)
       Y2 = YC(NBE(I,2))-YC(I)
       Y3 = YC(NBE(I,3))-YC(I)
      if(spherical == 1)then
       X1_DP = XC(I)
       Y1_DP = YC(I)
       X2_DP = XC(NBE(I,1))
       Y2_DP = YC(NBE(I,1))
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       X1=SIDE

       X2_DP=XC(NBE(I,2))
       Y2_DP=YC(NBE(I,2))
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       X2=SIDE

       X2_DP=XC(NBE(I,3))
       Y2_DP=YC(NBE(I,3))
       CALL ARCX(X1_DP,Y1_DP,X2_DP,Y2_DP,SIDE)
       X3=SIDE

       Y1=TPI*Y1
       Y2=TPI*Y2
       Y3=TPI*Y3
      else
       X1=XC(NBE(I,1))-XC(I)
       X2=XC(NBE(I,2))-XC(I)
       X3=XC(NBE(I,3))-XC(I)
      endif
       X1=X1/1000.0_SP
       X2=X2/1000.0_SP
       X3=X3/1000.0_SP
       Y1=Y1/1000.0_SP
       Y2=Y2/1000.0_SP
       Y3=Y3/1000.0_SP
       delt=(x1*y2-x2*y1)**2+(x1*y3-x3*y1)**2+(x2*y3-x3*y2)**2
       delt=delt*1000.

       a1u(i,1)=(y1+y2+y3)*(x1*y1+x2*y2+x3*y3)- &
                (x1+x2+x3)*(y1**2+y2**2+y3**2)
       a1u(i,1)=a1u(i,1)/delt
       a1u(i,2)=(y1**2+y2**2+y3**2)*x1-(x1*y1+x2*y2+x3*y3)*y1
       a1u(i,2)=a1u(i,2)/delt
       a1u(i,3)=(y1**2+y2**2+y3**2)*x2-(x1*y1+x2*y2+x3*y3)*y2
       a1u(i,3)=a1u(i,3)/delt
       a1u(i,4)=(y1**2+y2**2+y3**2)*x3-(x1*y1+x2*y2+x3*y3)*y3
       a1u(i,4)=a1u(i,4)/delt

       a2u(i,1)=(x1+x2+x3)*(x1*y1+x2*y2+x3*y3)- &
                (y1+y2+y3)*(x1**2+x2**2+x3**2)
       a2u(i,1)=a2u(i,1)/delt
       a2u(i,2)=(x1**2+x2**2+x3**2)*y1-(x1*y1+x2*y2+x3*y3)*x1
       a2u(i,2)=a2u(i,2)/delt
       a2u(i,3)=(x1**2+x2**2+x3**2)*y2-(x1*y1+x2*y2+x3*y3)*x2
       a2u(i,3)=a2u(i,3)/delt
       a2u(i,4)=(x1**2+x2**2+x3**2)*y3-(x1*y1+x2*y2+x3*y3)*x3
       a2u(i,4)=a2u(i,4)/delt
     end if
 if (spherical == 1)then
     X1=XM(NV(I,1))
     X2=XM(NV(I,2))
     X3=XM(NV(I,3))
     Y1=YM(NV(I,1))
     Y2=YM(NV(I,2))
     Y3=YM(NV(I,3))

!---------------------------------------
     CALL ARC(X2,Y2,X3,Y3,side1)
     CALL ARC(X3,Y3,X1,Y1,side2)
     CALL ARC(X1,Y1,X2,Y2,side3)
     CALL AREA(side1,side2,side3,art)
!---------------------------------------

     AI1=TPI*(Y2-Y3)
     AI2=TPI*(Y3-Y1)
     AI3=TPI*(Y1-Y2)
     CALL ARCX(x2,y2,x3,y3,side)
     BI1=side
     CALL ARCX(x3,y3,x1,y1,side)
     BI2=side
     CALL ARCX(x1,y1,x2,y2,side)
     BI3=side

     x2_dp = xc(i)
     y2_dp = yc(i)
     call ARCC(x1,y1,x2_dp,y2_dp,xxc1,yyc1)
     call ARCC(x2,y2,x2_dp,y2_dp,xxc2,yyc2)
     call ARCC(x3,y3,x2_dp,y2_dp,xxc3,yyc3)

     XTMP1  = X1*TPI-XC(I)*TPI
     XTMP2  = X2*TPI-XC(I)*TPI
     XTMP3  = X3*TPI-XC(I)*TPI
     XTMP11 = X1-XC(I)
     XTMP21 = X2-XC(I)
     XTMP31 = X3-XC(I)

     IF(XTMP11 >  180.0_SP)THEN
       XTMP1 = -360.0_SP*TPI+XTMP1
     ELSE IF(XTMP11 < -180.0_SP)THEN
       XTMP1 =  360.0_SP*TPI+XTMP1
     END IF
     IF(XTMP21 >  180.0_SP)THEN
       XTMP2 = -360.0_SP*TPI+XTMP2
     ELSE IF(XTMP21 < -180.0_SP)THEN
       XTMP2 =  360.0_SP*TPI+XTMP2
     END IF
     IF(XTMP31 >  180.0_SP)THEN
       XTMP3 = -360.0_SP*TPI+XTMP3
     ELSE IF(XTMP31 < -180.0_SP)THEN
       XTMP3 =  360.0_SP*TPI+XTMP3
     END IF

     CI1=XTMP2*TPI*(Y3-YC(I))*cos(d2r*YYC2)-&
         XTMP3*TPI*(Y2-YC(I))*cos(d2r*YYC3)

     CI2=XTMP3*TPI*(Y1-YC(I))*cos(d2r*YYC3)-&
         XTMP1*TPI*(Y3-YC(I))*cos(d2r*YYC1)

     CI3=XTMP1*TPI*(Y2-YC(I))*cos(d2r*YYC1)-&
         XTMP2*TPI*(Y1-YC(I))*cos(d2r*YYC2)

 else
     x1=xm(NV(i,1))-xc(i)
     x2=xm(NV(i,2))-xc(i)
     x3=xm(NV(i,3))-xc(i)
     y1=ym(NV(i,1))-yc(i)
     y2=ym(NV(i,2))-yc(i)
     y3=ym(NV(i,3))-yc(i)

     ai1=y2-y3
     ai2=y3-y1
     ai3=y1-y2
     bi1=x3-x2
     bi2=x1-x3
     bi3=x2-x1
     ci1=x2*y3-x3*y2
     ci2=x3*y1-x1*y3
     ci3=x1*y2-x2*y1
     ART = (XM(NV(I,2)) - XM(NV(I,1))) * (YM(NV(I,3)) - YM(NV(I,1))) - &
           (XM(NV(I,3)) - XM(NV(I,1))) * (YM(NV(I,2)) - YM(NV(I,1)))
     ART = ABS(.5_SP*ART)
 endif

     aw0(i,1)=-ci1/2./art
     aw0(i,2)=-ci2/2./art
     aw0(i,3)=-ci3/2./art
     awx(i,1)=-ai1/2./art
     awx(i,2)=-ai2/2./art
     awx(i,3)=-ai3/2./art
     awy(i,1)=-bi1/2./art
     awy(i,2)=-bi2/2./art
     awy(i,3)=-bi3/2./art
   end do
!
!--------boundary cells------------------------------------------------!
!
   do i=1,NGL
     if(isbce(i) > 1) then
       do j=1,4
         a1u(i,j)=0.0_SP
         a2u(i,j)=0.0_SP
       end do
     else if(isbce(i) == 1) then
       do j=1,3
         if(nbe(i,j) == 0) jj=j
       end do
       j1=jj+1-int((jj+1)/4)*3
       j2=jj+2-int((jj+2)/4)*3
       x1=xm(NV(i,j1))-xc(i)
       x2=xm(NV(i,j2))-xc(i)
       y1=ym(NV(i,j1))-yc(i)
       y2=ym(NV(i,j2))-yc(i)

    if (spherical == 1)then
       TY1=0.5*(YM(NV(I,J1))+YC(I))
       TY2=0.5*(YM(NV(I,J2))+YC(I))

       XTMP1  = xm(NV(i,j1))*TPI-xc(i)*TPI
       XTMP2  = xm(NV(i,j2))*TPI-xc(i)*TPI
       XTMP11 = xm(NV(i,j1))-xc(i)
       XTMP21 = xm(NV(i,j2))-xc(i)
       IF(XTMP11 >  180.0_SP)THEN
         XTMP1 = -360.0_SP*TPI+XTMP1
       ELSE IF(XTMP11 < -180.0_SP)THEN
         XTMP1 =  360.0_SP*TPI+XTMP1
       END IF
       IF(XTMP21 >  180.0_SP)THEN
         XTMP2 = -360.0_SP*TPI+XTMP2
       ELSE IF(XTMP21 < -180.0_SP)THEN
         XTMP2 =  360.0_SP*TPI+XTMP2
       END IF

       X1=XTMP1*cos(d2r*TY1)
       X2=XTMP2*cos(d2r*TY2)
       Y1=TPI*Y1
       Y2=TPI*Y2
    endif
       delt=x1*y2-x2*y1
       b1=(y2-y1)/delt
       b2=(x1-x2)/delt
       deltx=xm(NV(i,j1))-xm(NV(i,j2))
       delty=ym(NV(i,j1))-ym(NV(i,j2))

    if (spherical == 1)then
       x1_dp=XM(NV(I,J1))
       y1_dp=YM(NV(I,J1))
       x2_dp=XM(NV(I,J2))
       y2_dp=YM(NV(I,J2))
       call ARCX(x2_dp,y2_dp,x1_dp,y1_dp,side)

       DELTX=side
       DELTY=TPI*DELTY

    endif

       x1=xc(nbe(i,j1))-xc(i)
       x2=xc(nbe(i,j2))-xc(i)
       y1=yc(nbe(i,j1))-yc(i)
       y2=yc(nbe(i,j2))-yc(i)
    if (spherical == 1)then
       ty1=0.5*(yc(nbe(i,j1))+yc(i))
       ty2=0.5*(yc(nbe(i,j2))+yc(i))

       XTMP1  = xc(nbe(i,j1))*TPI-xc(i)*TPI
       XTMP2  = xc(nbe(i,j2))*TPI-xc(i)*TPI
       XTMP11 = xc(nbe(i,j1))-xc(i)
       XTMP21 = xc(nbe(i,j2))-xc(i)
       IF(XTMP11 >  180.0_SP)THEN
         XTMP1 = -360.0_SP*TPI+XTMP1
       ELSE IF(XTMP11 < -180.0_SP)THEN
         XTMP1 =  360.0_SP*TPI+XTMP1
       END IF
       IF(XTMP21 >  180.0_SP)THEN
         XTMP2 = -360.0_SP*TPI+XTMP2
       ELSE IF(XTMP21 < -180.0_SP)THEN
         XTMP2 =  360.0_SP*TPI+XTMP2
       END IF

       X1=XTMP1*COS(D2R*TY1)
       X2=XTMP2*COS(D2R*TY2)
       Y1=TPI*Y1
       Y2=TPI*Y2

    endif
       temp1=x1*y2-x2*y1

       if(abs(temp1).lt.1.e-6_SP)  then
         print*, 'shape_f of solid b. c. temp1=0'
         print*, 'i,jj,j1,j2,x1,x2,y1,y2'
         print*, i,jj,j1,j2,x1,x2,y1,y2
         print*, 'x1*y2==',x1*y2
         print*, 'x2*y1==',x2*y1
         stop
       end if

       a1u(i,1)=0.0_SP
       a1u(i,jj+1)=0.0_SP
       a1u(i,j1+1)=0.0_SP
       a1u(i,j2+1)=0.0_SP

       a2u(i,1)=0.0_SP
       a2u(i,jj+1)=0.0_SP
       a2u(i,j1+1)=0.0_SP
       a2u(i,j2+1)=0.0_SP
     end if
   end do

   return

 end subroutine SHAPE_COEF_GCN

   SUBROUTINE ARCX(XX1,YY1,XX2,YY2,ARCX1)

   IMPLICIT NONE
   INTEGER I,NX
   PARAMETER(NX=500)
   REAL(sp) :: XX1,YY1,XX2,YY2,ARCX1
   REAL(sp) :: X1,Y1,X2,Y2,TY
   REAL(sp) :: XTMP

   IF(XX1 == XX2)THEN
     ARCX1=0.
   ELSE
     X1=XX1*d2r
     Y1=YY1*d2r

     X2=XX2*d2r
     Y2=YY2*d2r

     XTMP  = X2-X1
     IF(XTMP >  PI)THEN
       XTMP = -2*PI+XTMP
     ELSE IF(XTMP < -PI)THEN
       XTMP =  2*PI+XTMP
     END IF

     TY=0.5*(Y2+Y1)
     ARCX1=EARTH*COS(TY)*XTMP
   END IF

   RETURN
   END SUBROUTINE ARCX

   SUBROUTINE ARCC(XX1,YY1,XX2,YY2,XXC,YYC)
   IMPLICIT NONE
   REAL(sp) :: XXC,YYC,XX1,YY1,XX2,YY2
   REAL(sp) :: X1,Y1,X2,Y2

   X1=XX1*D2R
   Y1=YY1*D2R

   X2=XX2*D2R
   Y2=YY2*D2R

   XXC=COS(Y1)*SIN(X1)+COS(Y2)*SIN(X2)
!   XXC=XXC/(COS(Y1)*COS(X1)+COS(Y2)*COS(X2))
!   XXC=ATAN(XXC)
   XXC=ATAN2(XXC,(COS(Y1)*COS(X1)+COS(Y2)*COS(X2)))
   XXC=XXC/D2R

!   IF(XXC .LT. 0.0) XXC=180.0+XXC
   IF(XXC < 0.0) XXC=360.0+XXC

   YYC=COS(Y1)*COS(Y1)+COS(Y2)*COS(Y2)+2.*COS(Y1)*COS(Y2)*COS(X1-X2)
!   YYC=SQRT(YYC)/(SIN(Y1)+SIN(Y2))
   YYC=ATAN2(SQRT(YYC),(SIN(Y1)+SIN(Y2)))
!   YYC=ATAN(YYC)
   YYC=90.-YYC/D2R

   RETURN
   END SUBROUTINE ARCC

   SUBROUTINE AREA(SIDE1,SIDE2,SIDE3,AREA1)
!--------------------------------------------------------------------
!      function:
!           calculate the area of a triangle on a spherical plane
!      input:
!           side1,side2 and side3: are 3 arc lenth for one triangle
!      output:
!           areal: is area of a triangle on a spherical plane
!--------------------------------------------------------------------
   IMPLICIT NONE
   REAL(sp) :: SIDE1,SIDE2,SIDE3,AREA1
   REAL(sp) :: PSUM,PM,QMJC

   SIDE1=SIDE1/EARTH
   SIDE2=SIDE2/EARTH
   SIDE3=SIDE3/EARTH
   IF(SIDE1 == 0. .OR. SIDE2 == 0. .OR. SIDE3 == 0.)THEN
     AREA1=0.
   ELSE
     PSUM=0.5*(SIDE1+SIDE2+SIDE3)
     PM=SIN(PSUM)*SIN(PSUM-SIDE1)*SIN(PSUM-SIDE2)*SIN(PSUM-SIDE3)
     PM=SQRT(PM)/(2.0*COS(SIDE1*0.5)*COS(SIDE2*0.5)*COS(SIDE3*0.5))
     QMJC = 2.0*ASIN(PM)

     AREA1=EARTH*EARTH*QMJC

   END IF

   RETURN
   END SUBROUTINE AREA
!===================================================================================|

   SUBROUTINE ARC(XX1,YY1,XX2,YY2,ARCL)
!----------------------------------------------------------------------------
!      function:
!           calculate the arc lenth for given two point on the spherical plane
!      input:
!           xx1,yy1,xx2,yy2 :are longitude and latitude of two points
!      output:
!           arcl :  arc lenth of two points in spherical plane
!-----------------------------------------------------------------------------

!  solve the arc length through the earth center
   IMPLICIT NONE
   REAL(sp) :: X1,Y1,X2,Y2,XA,YA,ZA,XB,YB,ZB,AB,AOB,ARCL
   REAL(sp) :: XX1,YY1,XX2,YY2

   X1=XX1*D2R
   Y1=YY1*D2R

   X2=XX2*D2R
   Y2=YY2*D2R

   XA=COS(Y1)*COS(X1)
   YA=COS(Y1)*SIN(X1)
   ZA=SIN(Y1)

   XB=COS(Y2)*COS(X2)
   YB=COS(Y2)*SIN(X2)
   ZB=SIN(Y2)

   AB=SQRT((XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2)
   AOB=(2.-AB*AB)/2.
   AOB=ACOS(AOB)
   ARCL=EARTH*AOB

   RETURN
   END SUBROUTINE ARC

End Module fvcom_driver
