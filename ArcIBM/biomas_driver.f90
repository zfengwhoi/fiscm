Module biomas_driver

!=======================================================================
! BIOMAS Driver
!
! Description
!    - Read and setup mesh
!    - Advect
!    - Diffuse (not implemented yet)
!    - Locate Particles
!
! Comments:
! - Advection can use either 1st Order Euler Step or 4-stage RK4 scheme. Recommend RK4.
! - Current code is only compatible with pan-Arctic BIOMAS grid and forcing.
!   a) Arakawa B-grid (horizontal): scalar & vector fields in separate grids.
!   b) Z-grid (vertical): ocean levels are based on actual depths. 
!                         w is at bottom of each level.
!   c) "North pole" is dispaced on the land of Alaska.
!   d) Periodic boundary in x-direction (longitudinal or zonal). 
!      Open boundary at ~39N and closed boudary at "north pole".      
!
! zfeng 10/14/2014 Started this new module to handle structured BIOMAS forcing.
!
! Subroutines included (texts to be added): 
!
!
!=======================================================================

     use gparms
     use mod_igroup
     use mod_forcing
     use utilities

     Implicit None

     !dimensions of BIOMAS mesh
     integer :: NSLON, NSLAT, NVLON, NVLAT, NLVLS
     logical              :: mesh_setup = .false.
     real(sp),allocatable :: slon(:,:),slat(:,:)
     real(sp),allocatable :: vlon(:,:),vlat(:,:)
     real(sp),allocatable :: HTN(:,:),HTE(:,:),HUS(:,:),HUW(:,:)
     real(sp),allocatable :: zw(:),zt(:)
     integer, allocatable :: mask(:,:)
     real(sp),allocatable :: bathy(:,:)
     
     ! save grid variables below for future calculations
     save :: slon,slat,vlon,vlat,HTN,HTE,HUS,HUW,zw,zt,mask,bathy
     logical       :: grid_metrics

     real(sp),parameter:: NorPole = 89.5
     ! when particle is very close to north pole, linear interplotion
     ! won't work, so use an average value of four neighboring grid points    

contains

!----------------------------------------------------
! Read the mesh, mask, ocean level thickness, bathymetry
! Add mesh to output files for viz
!----------------------------------------------------
subroutine structured_model_init(ng,g,lsize,varlist)
  use utilities, only      :  drawline,ncdchk
  integer, intent(in)      :: ng
  type(igroup), intent(in) :: g(ng)
  integer, intent(inout)   :: lsize
  character(len=*)         :: varlist(max_state_vars)
  !----------------------------------------------------
  character(len=mstr)      :: msg
  integer                  :: dimid,varid,fid
  character(len=fstr)      :: dname
  integer                  :: i,j,k,n,ierr,ofid
  integer                  :: slon_dimid,slat_dimid,vlon_dimid,vlat_dimid
  integer                  :: htn_dimid,hte_dimid,hus_dimid,huw_dimid
  integer                  :: mask_dimid,thick_dimid,zw_dimid,zt_dimid,bathy_dimid
  integer                  :: slon_varid,slat_varid,vlon_varid,vlat_varid
  integer                  :: htn_varid,hte_varid,hus_varid,huw_varid
  integer                  :: mask_varid,thick_varid,zw_varid,zt_varid,bathy_varid
  
  character(len=fstr) :: slon_char,slat_char,vlon_char,vlat_char
  character(len=fstr) :: htn_char,hte_char,hus_char,huw_char
  character(len=fstr) :: zw_char,zt_char,thick_char,mask_char,h_char
  character(len=fstr) :: u_char,v_char,ww_char,phyto_char,temp_char

  !return if group spatial dimensions are all 0-d or 1-d
  if(maxval(g%space_dim) < 2) return

  !get the forcing file netcdf id
  fid = get_ncfid()

        slon_char = 'slon'
        slat_char = 'slat'
        vlon_char = 'vlon'
        vlat_char = 'vlat'
        htn_char  = 'HTN'
        hte_char  = 'HTE'
        hus_char  = 'HUS'
        huw_char  = 'HUW'
!        kh_char   = 'kh'


        if(spherical==1) then
           x_char = 'slon'
           y_char = 'slat'
!        else
           x_char = 'x'
           y_char = 'y'
        endif

           zw_char = 'zw'
           zt_char = 'zt'
       thick_char  = 'thickness'
        mask_char  = 'mask'
           h_char  = 'h'
!        zeta_char ='zeta'
           u_char = 'u'
           v_char = 'v'   
          ww_char = 'ww'
          ua_char = 'u2d'
          va_char = 'v2d'
      ! phyto_char  = 'phyto'
      !  temp_char = 'temp' 
!          wu_char ='uuwind'   ! 'uwind_speed'
!          wv_char ='vvwind'   ! 'vwind_speed'

  do n=1,ng ! groups 
!    if (wind_type == 1)then
!      lsize = lsize + 1 ; varlist(lsize) = wu_char
!      lsize = lsize + 1 ; varlist(lsize) = wv_char
!    endif

    ! if biology is activated, temperature is needed
!    if(biology) then 
!      lsize = lsize + 1;  varlist(lsize) = temp_char
!    endif

    ! if food is activated, phytoplankton is needed
!    if(g(n)%bio_food .AND. g(n)%bio_SPM) then
!      lsize = lsize + 1;  varlist(lsize) = phyto_char
!    endif

    if(g(n)%space_dim == 2) then
      lsize = lsize + 1 ; varlist(lsize) = ua_char
      lsize = lsize + 1 ; varlist(lsize) = va_char
    elseif(g(n)%space_dim ==3) then
      lsize = lsize + 1 ; varlist(lsize) = u_char
      lsize = lsize + 1 ; varlist(lsize) = v_char
      lsize = lsize + 1 ; varlist(lsize) = ww_char
      
!      if(g(n)%vdiff_type .gt. 0) then
!        lsize = lsize + 1 ; varlist(lsize) = kh_char
!      endif
    endif
  end do  ! groups

  !determine number of longitude points of scalar grid 
  msg = "dimension 'slon' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, slon_char, dimid ), msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, NSLON ))

  !determine number of latitude points of scalar grid
  msg = "dimension 'slat' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, slat_char, dimid ), msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, NSLAT ))

  !determine number of longitude points of vector grid
  msg = "dimension 'vlon' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, vlon_char, dimid ), msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, NVLON))

  !determine number of latitude points of vector grid
  msg = "dimension 'vlat' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, vlat_char, dimid ), msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, NVLAT))
  
  !determine number of vertical layers
  msg = "dimension 'z' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, zw_char, dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, NLVLS ))

  !allocate dataspace and read in mesh
  allocate(slon(NSLON,NSLAT))
  msg = "error reading longitude coordinate of scalar grid"
  call ncdchk( nf90_inq_varid(fid,slon_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, slon),msg)

  allocate(slat(NSLON,NSLAT))
  msg = "error reading latitude coordinate of scalar grid"
  call ncdchk( nf90_inq_varid(fid,slat_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, slat),msg)

  allocate(vlon(NVLON,NVLAT))
  msg = "error reading longitude coordinate of vector grid"
  call ncdchk( nf90_inq_varid(fid,vlon_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, vlon),msg)

  allocate(vlat(NVLON,NVLAT))
  msg = "error reading latitude coordinate of vector grid"
  call ncdchk( nf90_inq_varid(fid,vlat_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, vlat),msg)

  allocate(HTN(NVLON,NVLAT))
  msg = "error reading length of north side of scalar grid"
  call ncdchk( nf90_inq_varid(fid,htn_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, HTN),msg)

  allocate(HTE(NVLON,NVLAT))
  msg = "error reading length of east side of scalar grid"
  call ncdchk( nf90_inq_varid(fid,hte_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, HTE),msg)

  allocate(HUS(NVLON,NVLAT))
  msg = "error reading length of north side of scalar grid"
  call ncdchk( nf90_inq_varid(fid,hus_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, HUS),msg)

  allocate(HUW(NVLON,NVLAT))
  msg = "error reading length of north side of scalar grid"
  call ncdchk( nf90_inq_varid(fid,huw_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, HUW),msg)

  allocate(zw(NLVLS))
  msg = "error reading ocean level depth"
  call ncdchk( nf90_inq_varid(fid,zw_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, zw),msg)

  allocate(zt(NLVLS))
  msg = "error reading ocean level center depth"
  call ncdchk( nf90_inq_varid(fid,zt_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, zt),msg)

  allocate(bathy(NSLON,NSLAT))
  msg = "error reading bathymetry h"
  call ncdchk( nf90_inq_varid(fid,h_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, bathy),msg)

  allocate(mask(NSLON,NSLAT))
  msg = "error reading BIOMAS grid mask"
  call ncdchk( nf90_inq_varid(fid,mask_char,varid),msg )
  call ncdchk( nf90_get_var(fid, varid, mask),msg)

  call drawline('-')
  write(*,*)'BIOMAS mesh stats '
  call drawline('-')
  write(*,*) 'Number of scalar grid longitude points: ',NSLON
  write(*,*) 'Number of scalar grid latitude points:  ',NSLAT
  write(*,*) 'Number of vector grid longitude points: ',NVLON
  write(*,*) 'Number of vector grid latitude points:  ',NVLAT
  write(*,*) 'Number of vertical layers:              ',NLVLS
  write(*,*) 'xmin:                                   ',minval(slon)
  write(*,*) 'xmax:                                   ',maxval(slon)
  write(*,*) 'ymin:                                   ',minval(slat)
  write(*,*) 'ymax:                                   ',maxval(slat)

  !flag that mesh is setup
  mesh_setup = .true.

end subroutine structured_model_init
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine initial_indices(ng,g)
! -- Find indices (i,j,k) based on initial particle
!   locations and then store indices in g(ng) to update states.
! -- This subroutine uses linear search algorithm to search for the
!   whole domain because initial indices are unknown.

  implicit none
  integer, intent(in) :: ng
  type(igroup), intent(inout) :: g(ng)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: z(:)
  integer , pointer :: i_index(:)
  integer , pointer :: j_index(:)
  integer , pointer :: k_index(:)
!  integer , pointer :: mask_value(:)
  integer , pointer :: istatus(:)
!  character(len=*), optional :: option
  integer :: np,n,p

  do n=1,ng  ! groups
    if(g(n)%space_dim < 2)cycle

    ! set dimensions
    np = g(n)%nind

  if (g(n)%space_dim == 2) then ! 2-D case
    ! set pointers to states
    call get_state('x',g(n),x)
    call get_state('y',g(n),y)
    call get_state('i',g(n),i_index)
    call get_state('j',g(n),j_index)
!    call get_state('mask',g(n),mask_value)
    call get_state('status',g(n),istatus)

    ! find grid indices (i_index,j_index) containing particles
    do p = 1, np
     if(istatus(p) <= 0) cycle
     call linear_search(x(p),y(p),i_index(p),j_index(p))
     ! when new indices could not be found, print error msg and particle #
     if (i_index(p)==-999 .OR. j_index(p)==-999) then
        write(*,*)  "Error in 'initial_indices': cannot find i_index or j_index!" 
        write(*,*)  "Particle Number is: ", p
     endif
    enddo 

    ! make a neighbor_search to make sure special cases are taken care of
    call update_ij(np,x(:),y(:),i_index(:),j_index(:),istatus(:))

  elseif (g(n)%space_dim == 3) then ! 3-D case
    ! set pointers to states
    call get_state('x',g(n),x)
    call get_state('y',g(n),y)
    call get_state('z',g(n),z)
    call get_state('i',g(n),i_index)
    call get_state('j',g(n),j_index)
!    call get_state('mask',g(n),mask_value)
    call get_state('k',g(n),k_index)
    call get_state('status',g(n),istatus)

    ! find grid indices (i_index,j_index) containing particles
    do p = 1, np
     if(istatus(p) <= 0) cycle
     call linear_search(x(p),y(p),i_index(p),j_index(p))
    ! when new indices could not be found
     if (i_index(p)==-999 .OR. j_index(p)==-999) then
        write(*,*) "Error in 'initial_indices': cannot find i_index or j_index!"        
        write(*,*) 'Particle Number is: ',p
     endif
    enddo

    ! make a neighbor_search to make sure special cases are taken care of
    call update_ij(np,x(:),y(:),i_index(:),j_index(:),istatus(:))

    ! find grid index (k) in the vertical
    call update_k(np,z(:),k_index(:),istatus(:))

  endif

  enddo ! group loop
  nullify(x,y,z,i_index,j_index,k_index,istatus)
end subroutine initial_indices
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine update_indices(ng,g)
! -- 'update_indices' updates indices (i,j,k) based on present particle
!   locations and then replace indices in g(ng) to update states.
! -- Use 'neighbor_search' algorithm, which only searchs for
!  neighboring 8 grid points.

  implicit none
  integer, intent(in) :: ng
  type(igroup), intent(inout) :: g(ng)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: z(:)
  integer , pointer :: i_index(:)
  integer , pointer :: j_index(:)
  integer , pointer :: k_index(:)
  integer , pointer :: istatus(:)
!  character(len=*), optional :: option
  integer :: np,n

  do n=1,ng
    if(g(n)%space_dim < 2) cycle

    ! set dimensinos
    np = g(n)%nind
    
  if (g(n)%space_dim == 2) then ! 2-D case
    ! set pointers to states
    call get_state('x',g(n),x)
    call get_state('y',g(n),y)
    call get_state('i',g(n),i_index)
    call get_state('j',g(n),j_index)
    call get_state('status',g(n),istatus)
    ! update grid indice (i,j) containing particles
    call update_ij(np,x(:),y(:),i_index(:),j_index(:),istatus(:))

  elseif (g(n)%space_dim == 3) then ! 3-D case
    ! set pointers to states
    call get_state('x',g(n),x)
    call get_state('y',g(n),y)
    call get_state('z',g(n),z)
    call get_state('i',g(n),i_index)
    call get_state('j',g(n),j_index)
    call get_state('k',g(n),k_index)
    call get_state('status',g(n),istatus)
    ! update grid indice (i,j) containing particles
    call update_ij(np,x(:),y(:),i_index(:),j_index(:),istatus(:))
    ! update grid index (k) in the vertical
    call update_k(np,z(:),k_index(:),istatus(:))
  endif

  enddo ! group loop
  nullify(x,y,z,i_index,j_index,k_index,istatus)
end subroutine update_indices

!--------------------------------------------------------------------
 subroutine update_ij(np,x,y,i_index,j_index,istatus)
! 'update_ij' update indices (i_index,j_index) based on particle 
! location (x,y) and 'neighbor_search' algorithm.
!
!
   use gparms
   implicit none

   integer, intent(in)    ::  np
   real(sp),intent(in)    ::  x(np),y(np)
   integer, intent(inout) ::  i_index(np),j_index(np),istatus(np)
   integer                ::  p,i,j

   do p = 1, np
     if(istatus(p) <= 0) cycle
     call neighbor_search(x(p),y(p),i_index(p),j_index(p))

     ! when new indices could not be found
     if (i_index(p)==-999 .OR. j_index(p)==-999) then
        write(*,*) "Error in 'update_ij': cannot update ij"
        write(*,*) 'Particle Number ',p
     endif
   enddo

 end subroutine update_ij
!--------------------------------------------------------------------

!--------------------------------------------------------------------
 subroutine update_k(np,z,k_index,istatus,option)
! 'update_k' update index k based on particle depth z
!-------------------------------------------------------------------
   use gparms
   implicit none

   integer, intent(in)    ::  np
   real(sp),intent(in)    ::  z(np)
   integer, intent(inout) ::  k_index(np),istatus(np)
   character(len=*), optional :: option
   integer                :: p,k

   do p = 1, np

     if(istatus(p) <= 0) cycle
     call find_k(z(p),k)
     k_index(p) = k

     ! write error message if k_nrst NOT found
   if(k == -999) then
     write(*,*) "Error in 'update_k': cannot update k_index "
     write(*,*) "Particle Number ", p
   endif

   enddo

 end subroutine update_k

!--------------------------------------------------------------------
 subroutine linear_search(lon,lat,i_nrst,j_nrst)
 !--Find nearest grid indices (i_nrst,k_nrst) on a structured mesh
 !  to a particle location (lon,lat) using linear search algorithm
 !--This algorithm searches the whole mesh domain and is time-
 !  consuming, so is only used in the beginning when the indices of
 !  particle is unknown.
 ! zfeng 10/27/2014
 !-----------------------------------------------------------------
   use gparms
   use utilities
   implicit none

   real(sp), intent(in)   :: lon,lat
   integer,  intent(out)  :: i_nrst,j_nrst
   integer  :: i,j
   real(sp) :: dist(nx,ny)


   ! make sure mesh is already set up
   if (.not. mesh_setup) then
     write(*,*) 'Error in linear_search:'
     write(*,*) 'BIOMAS mesh is not setup correctly'
     stop
   endif

   ! initialize
   i_nrst = -999
   j_nrst = -999

   ! Calculate distance between particle location and scalar grid points
   do j = 1, ny
     do i = 1, nx
        dist(i,j) = gcdist( lon, lat, slon(i,j),slat(i,j) )
     enddo
   enddo

   ! find (i_nrst,j_nrst) based on minimum distance
   call nearest_ij(dist,i_nrst,j_nrst)

   ! write error message of either i_nrst or j_nrst NOT found
   if(i_nrst == -999 .OR. j_nrst == -999) then
     write(*,*) "Error in 'linear_search': indices (i_nrst,j_nrst) for particle NOT found!"
   endif

 end subroutine linear_search

 !--------------------------------------------------------------------
 subroutine neighbor_search(lon,lat,i_nrst,j_nrst)
 !-- Find nearst grid indices (i_nrst,k_nrst) on a structured mesh
 !   to a point (lon,lat) using "neighbor search" algorithm
 !-- To save computational time, only search in 8 neighboring points
 !   because particle location change in one time step is subtle.
 !-- The minimal length of two neiboring grid points is 1780 m. Taking
 !   a extremly high velocity of 1 m/s and time step of 600 s,
 !   particle only travels 600 m, less than 1/3 of the minimal legth.
 !
 ! zfeng 10/28/2014
 !--------------------------------------------------------------------
   use gparms
   use utilities
   implicit none

   real(sp), intent(in)   :: lon,lat
   integer,  intent(inout)  :: i_nrst,j_nrst
   integer  :: i,j, i_old, j_old, i_new, j_new
   real(sp) :: min_dist
   real(sp),allocatable :: dist(:,:)
   real(sp), parameter :: large_val = 1E7

   ! make sure mesh is already set up
   if (.not. mesh_setup) then
     write(*,*) "Error in 'neighbor_search':"
     write(*,*) "BIOMAS mesh is not setup correctly"
     stop
   endif

   ! Calculate distance
   i_old = i_nrst
   j_old = j_nrst
   allocate(dist(i_old-1:i_old+1,j_old-1:j_old+1))

   if(j_old==1 .AND. i_old==1) then ! corner-1

      dist(i_old-1,j_old-1) = large_val
      dist(i_old  ,j_old-1) = large_val ! give a large dist value to outside point
      dist(i_old+1,j_old-1) = large_val
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(nx,j_old       ),slat(nx     ,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(i_old+1,j_old  ),slat(i_old+1,j_old  ))
      dist(i_old-1,j_old+1) = gcdist(lon,lat,slon(nx     ,j_old+1),slat(nx     ,j_old+1))
      dist(i_old  ,j_old+1) = gcdist(lon,lat,slon(i_old  ,j_old+1),slat(i_old  ,j_old+1))
      dist(i_old+1,j_old+1) = gcdist(lon,lat,slon(i_old+1,j_old+1),slat(i_old+1,j_old+1))

   elseif(j_old==1 .AND. i_old==nx) then ! corner-2

      dist(i_old-1,j_old-1) = large_val
      dist(i_old  ,j_old-1) = large_val
      dist(i_old+1,j_old-1) = large_val
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(i_old-1,j_old  ),slat(i_old-1,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(1      ,j_old  ),slat(1      ,j_old  ))
      dist(i_old-1,j_old+1) = gcdist(lon,lat,slon(i_old-1,j_old+1),slat(i_old-1,j_old+1))
      dist(i_old  ,j_old+1) = gcdist(lon,lat,slon(i_old  ,j_old+1),slat(i_old  ,j_old+1))
      dist(i_old+1,j_old+1) = gcdist(lon,lat,slon(1      ,j_old+1),slat(1      ,j_old+1))

   elseif(j_old==ny .AND. i_old==nx) then ! corner-3

      dist(i_old-1,j_old-1) = gcdist(lon,lat,slon(i_old-1,j_old-1),slat(i_old-1,j_old-1))
      dist(i_old  ,j_old-1) = gcdist(lon,lat,slon(i_old  ,j_old-1),slat(i_old  ,j_old-1))
      dist(i_old+1,j_old-1) = gcdist(lon,lat,slon(1      ,j_old-1),slat(1      ,j_old-1))
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(i_old-1,j_old  ),slat(i_old-1,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(1      ,j_old  ),slat(1      ,j_old  ))
      dist(i_old-1,j_old+1) = large_val
      dist(i_old  ,j_old+1) = large_val
      dist(i_old+1,j_old+1) = large_val

   elseif(j_old==ny .AND. i_old==1) then ! corner-4

      dist(i_old-1,j_old-1) = gcdist(lon,lat,slon(nx     ,j_old-1),slat(nx     ,j_old-1))
      dist(i_old  ,j_old-1) = gcdist(lon,lat,slon(i_old  ,j_old-1),slat(i_old  ,j_old-1))
      dist(i_old+1,j_old-1) = gcdist(lon,lat,slon(i_old+1,j_old-1),slat(i_old+1,j_old-1))
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(nx     ,j_old  ),slat(nx     ,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(i_old+1,j_old  ),slat(i_old+1,j_old  ))
      dist(i_old-1,j_old+1) = large_val
      dist(i_old  ,j_old+1) = large_val
      dist(i_old+1,j_old+1) = large_val

   elseif(j_old==1 .AND. i_old>1 .AND. i_old<nx) then

      dist(i_old-1,j_old-1) = large_val
      dist(i_old  ,j_old-1) = large_val ! give a large dist value to outside point
      dist(i_old+1,j_old-1) = large_val
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(i_old-1,j_old  ),slat(i_old-1,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(i_old+1,j_old  ),slat(i_old+1,j_old  ))
      dist(i_old-1,j_old+1) = gcdist(lon,lat,slon(i_old-1,j_old+1),slat(i_old-1,j_old+1))
      dist(i_old  ,j_old+1) = gcdist(lon,lat,slon(i_old  ,j_old+1),slat(i_old  ,j_old+1))
      dist(i_old+1,j_old+1) = gcdist(lon,lat,slon(i_old+1,j_old+1),slat(i_old+1,j_old+1))

   elseif(j_old==ny .AND. i_old>1 .AND. i_old<nx) then

      dist(i_old-1,j_old-1) = gcdist(lon,lat,slon(i_old-1,j_old-1),slat(i_old-1,j_old-1))
      dist(i_old  ,j_old-1) = gcdist(lon,lat,slon(i_old  ,j_old-1),slat(i_old  ,j_old-1))
      dist(i_old+1,j_old-1) = gcdist(lon,lat,slon(i_old+1,j_old-1),slat(i_old+1,j_old-1))
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(i_old-1,j_old  ),slat(i_old-1,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(i_old+1,j_old  ),slat(i_old+1,j_old  ))
      dist(i_old-1,j_old+1) = large_val
      dist(i_old  ,j_old+1) = large_val
      dist(i_old+1,j_old+1) = large_val

   elseif(j_old>1 .AND. j_old<ny .AND. i_old==1) then

      dist(i_old-1,j_old-1) = gcdist(lon,lat,slon(nx     ,j_old-1),slat(nx     ,j_old-1))
      dist(i_old  ,j_old-1) = gcdist(lon,lat,slon(i_old  ,j_old-1),slat(i_old  ,j_old-1))
      dist(i_old+1,j_old-1) = gcdist(lon,lat,slon(i_old+1,j_old-1),slat(i_old+1,j_old-1))
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(nx     ,j_old  ),slat(nx     ,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(i_old+1,j_old  ),slat(i_old+1,j_old  ))
      dist(i_old-1,j_old+1) = gcdist(lon,lat,slon(nx     ,j_old+1),slat(nx     ,j_old+1))
      dist(i_old  ,j_old+1) = gcdist(lon,lat,slon(i_old  ,j_old+1),slat(i_old  ,j_old+1))
      dist(i_old+1,j_old+1) = gcdist(lon,lat,slon(i_old+1,j_old+1),slat(i_old+1,j_old+1))

   elseif(j_old>1 .AND. j_old<ny .AND. i_old==nx) then

      dist(i_old-1,j_old-1) = gcdist(lon,lat,slon(i_old-1,j_old-1),slat(i_old-1,j_old-1))
      dist(i_old  ,j_old-1) = gcdist(lon,lat,slon(i_old  ,j_old-1),slat(i_old  ,j_old-1))
      dist(i_old+1,j_old-1) = gcdist(lon,lat,slon(1      ,j_old-1),slat(1      ,j_old-1))
      dist(i_old-1,j_old  ) = gcdist(lon,lat,slon(i_old-1,j_old  ),slat(i_old-1,j_old  ))
      dist(i_old  ,j_old  ) = gcdist(lon,lat,slon(i_old  ,j_old  ),slat(i_old  ,j_old  ))
      dist(i_old+1,j_old  ) = gcdist(lon,lat,slon(1      ,j_old  ),slat(1      ,j_old  ))
      dist(i_old-1,j_old+1) = gcdist(lon,lat,slon(i_old-1,j_old+1),slat(i_old-1,j_old+1))
      dist(i_old  ,j_old+1) = gcdist(lon,lat,slon(i_old  ,j_old+1),slat(i_old  ,j_old+1))
      dist(i_old+1,j_old+1) = gcdist(lon,lat,slon(1      ,j_old+1),slat(1      ,j_old+1))

   else

      do j = j_old-1, j_old+1
         do i = i_old-1, i_old+1
            dist(i,j) = gcdist( lon, lat, slon(i,j),slat(i,j) )
         enddo
      enddo

   endif

   ! find (i_new,j_new) based on minimum distance
   min_dist = minval(dist)

   i_new = -999
   j_new = -999

   do j = j_old-1, j_old+1
      do i = i_old-1, i_old+1
        if( abs(dist(i,j)-min_dist) < 0.1 ) then
           i_new = i
           j_new = j
           if(i_new==0)    i_new = nx
           if(i_new==nx+1) i_new = 1
           exit
        endif
      enddo
   enddo

!   write(*,*) 'i_new=',i_new
!   write(*,*) 'j_new=',j_new

  if(i_new == -999 .OR. j_new == -999) then
    write(*,*) "Error in 'neighor_search': cannot find (i_new,j_new)!"
  endif

  ! subroutine return new indices
  i_nrst = i_new; j_nrst = j_new

  deallocate(dist)

 end subroutine neighbor_search

!------------------------------------------------------------------------------
!--'find_k' finds index k_below  where a particle is located between levels
!   [k_below-1, k_below], using particle depth and ocean level center depth 'zt'
!-- All quantities (u,v,temp ...) except vertical velocity(w) are stored at
!   ocean level centers.
! zfeng 10/27/2014
!------------------------------------------------------------------------------

 subroutine find_k(pdep,k_below)

   use gparms
   implicit none

   real(sp), intent(in)   :: pdep
   integer,  intent(out)  :: k_below
   integer                :: k

   ! initialize
   k_below = -999

   if( pdep <= zt(1) ) then  ! particle is above top level center
       k_below = 1
   else
     do k = 2, nz

       if( pdep > zt(k-1) .AND. pdep <= zt(k) ) then
         k_below = k; exit
       else

         k_below = k
!         write(*,*) "Warning in 'find_k': particle depth is below ocean depth!"
!         write(*,*) "Use forcing data at the bottom level for this particle!"
       endif

     enddo
   endif

! write error message if k_below NOT found
   if(k_below == -999) then
     write(*,*) "Error in 'find_k': index 'k_below' for particle NOT found!"
   endif

 end subroutine find_k

!----------------------------------------------------------------------------
!--'find_kw' finds index k_below  where a particle is located between levels
!   [k_below-1, k_below], using particle depth and ocean level bottom depth 'zw'
!-- Vertical velocity(ww) are stored at ocean level bottoms.
! zfeng 12/03/2014
!---------------------------------------------------------------------------
 subroutine find_kw(pdep,k_below)

   use gparms
   implicit none

   real(sp), intent(in)   :: pdep
   integer,  intent(out)  :: k_below
   integer                :: k

   ! initialize
   k_below = -999

   if( pdep <= zw(1) ) then  ! particle is above bottom of top level
       k_below = 1
   else

     do k = 2, nz

       if( pdep > zw(k-1) .AND. pdep <= zw(k) ) then
         k_below = k; exit
       else
         k_below = k
!         write(*,*) "Warning in 'find_k': particle depth is below ocean depth!"
!         write(*,*) "Use forcing data at the bottom level for this particle!"
       endif

     enddo
   endif

! write error message if k_below NOT found
   if(k_below == -999) then
     write(*,*) "Error in 'find_k': index 'k_below' for particle NOT found!"
   endif

 end subroutine find_kw


!--------------------------------------------------------------------
! - Search for nearest vector grid point to the particle location,
!   based on given information of particle lon & lat, and nearest
!   scalar grid point indices. Note current subroutine is only applied
!   to Arakawa B-Grid and periodic boundary in x-direction.
! - zfeng 11/04/2014
!
!
!                 V4(i-1,j)__________V1(i,j)
!                 |         |    P   |
!                 |         |        |
!                 |_________S(i,j)___|
!                 |         |        |
!                 |         |        |
!                 V3(i-1,j-1)________V2(i.j-1)
!
!
! P is particle location: P(lon,lat), and S is the nearest scalar grid point;
! S is scalar grid point at cell center;
! V is vector grid point at cell corner (B-Grid convention).
! zfeng 11/03/2014
!
subroutine search4corners(lon,lat,i_scalar,j_scalar,i_vector,j_vector)

     implicit none
     real(sp),     intent(in)          :: lon, lat
     integer,      intent(in)          :: i_scalar,j_scalar
     integer,      intent(out)         :: i_vector,j_vector

     real(sp)  :: x(4),y(4),dis(4), min_dis
     integer   :: i, icorner, i_temp, j_temp

     i_temp = i_scalar
     j_temp = j_scalar

   if(j_scalar == 1) then ! scalar point at the outside latitude boundary
       j_temp = 2       ! move grid point to j = 2
   endif

   if(i_scalar == 1) then ! periodic boundary in x-dir (longitude)
       x(1) = vlon(i_temp,j_temp)
       y(1) = vlat(i_temp,j_temp)
       x(2) = vlon(i_temp,j_temp-1)
       y(2) = vlat(i_temp,j_temp-1)
       x(3) = vlon(nx,j_temp-1)
       y(3) = vlat(nx,j_temp-1)
       x(4) = vlon(nx,j_temp)
       y(4) = vlat(nx,j_temp)
   else
      x(1) = vlon(i_temp,j_temp)     ! right-upper corner
      y(1) = vlat(i_temp,j_temp)     !
      x(2) = vlon(i_temp,j_temp-1)   ! right-lower corner
      y(2) = vlat(i_temp,j_temp-1)   !
      x(3) = vlon(i_temp-1,j_temp-1) ! left-lower corner
      y(3) = vlat(i_temp-1,j_temp-1) !
      x(4) = vlon(i_temp-1,j_temp)   ! left-upper corner
      y(4) = vlat(i_temp-1,j_temp)   !
   endif

   ! calculate great circle distance from particle to four corners
   do i = 1, 4
      dis(i) = gcdist(lon,lat,x(i),y(i))
   enddo

   min_dis = minval(dis)

   ! find which corner point is nearest
   icorner = -999 ! initialize
   do i = 1, 4
      if ( abs(dis(i)-min_dis)<=1 ) icorner = i
   enddo

   ! find corner indices
   select case(icorner)
   case(1)
       i_vector = i_temp; j_vector = j_temp
   case(2)
       i_vector = i_temp; j_vector = j_temp-1
   case(3)
       if (i_temp == 1) i_temp = nx+1
       i_vector = i_temp-1; j_vector = j_temp-1
   case(4)
       if (i_temp == 1) i_temp = nx+1
       i_vector = i_temp-1; j_vector = j_temp
   case default
     write(*,*) 'i_vector=',i_vector
     write(*,*) 'j_vector=',j_vector
     write(*,*) "Error in 'search4corners': cannot find nearest vector grid point!"
   end select

end subroutine search4corners

!-----------------------------------------------------------------
! driver to interpolate structured forcing onto particle positions
! zfeng 10/29/2014
!-----------------------------------------------------------------
subroutine interp_structured_forcing(ng,g,iframe)

      implicit none
      integer     , intent(in)                   :: ng,iframe
      type(igroup), intent(inout), dimension(ng) :: g

      integer,  pointer   :: istatus(:),i(:),j(:),k(:)
      real(sp), pointer   :: x(:),y(:),z(:),f(:)
      integer             :: n, v
      character(len=fstr) :: evar, svar
!      type(pvar),pointer  :: p

 do n=1,ng

        !skip 0-d groups
        if(g(n)%space_dim < 2)cycle

        !get the horizontal particle positions for the group
        call get_state('x',g(n),x)
        call get_state('y',g(n),y)
        call get_state('i',g(n),i)
        call get_state('j',g(n),j)
        call get_state('status',g(n),istatus)

        !get the vertical position if sim is 3-d
        if(g(n)%space_dim == 3)  then
           call get_state('z',g(n),z)
           call get_state('k',g(n),k)
        endif

    do v = 1, g(n)%next

       ! g(n)%next = 2

       ! get varname and type
       svar = g(n)%ext_var(v,1)
       ! 'svar' is T & phytoplankton  
       
       evar = g(n)%ext_var(v,2)
       ! 'evar' is temp & phyto
    
       ! get state variable data
       call get_state(svar,g(n),f)

     if(g(n)%space_dim == 2) then
         if (evar/='phyto' .AND. evar/='zoopl1') then 
          call bilinear_interp(g(n)%nind,x,y,i,j,istatus,evar,iframe,0,f)
         endif
       ! zfeng 12/09/2015: food-seeking in 2D mode
       if(g(n)%bio_food .AND. g(n)%food_source == 3) then
          call get_state('SPM',g(n),f)    ! subsurface phytoplankton maxima
          call interp_SPM_2D(g(n)%nind,x,y,i,j,istatus,'phyto',iframe,0,f)
       elseif (g(n)%bio_food .AND. g(n)%food_source ==4 ) then
          call get_state('SFM',g(n),f)    ! subsurface food (P+M) maxima
          call interp_SPM_2D(g(n)%nind,x,y,i,j,istatus,'zoopl1',iframe,0,f)
       endif

     elseif(g(n)%space_dim == 3) then
          call trilinear_interp(g(n)%nind,x,y,z,i,j,k,istatus,evar,iframe,0,f)
     
       ! zfeng 11/13/2015: food seeking in 3D mode
       if(g(n)%bio_food .AND. g(n)%food_source == 3) then
          call get_state('SPM',g(n),f)    ! subsurface phytoplankton maxima
          call interp_SPM(g(n)%nind,x,y,z,i,j,k,istatus,'phyto',iframe,0,f)
       elseif(g(n)%bio_food .AND. g(n)%food_source ==4 ) then
          call get_state('SFM',g(n),f)    ! subsurface food (P+M) maxima
          call interp_SPM(g(n)%nind,x,y,z,i,j,k,istatus,'zoopl1',iframe,0,f)
       endif
    
     endif

    enddo

  enddo ! group loop
      
 nullify(x,y,z,i,j,k,istatus)
end subroutine interp_structured_forcing


!---------------------------------------------------------------
subroutine bilinear_interp(np,x,y,ip,jp,istatus,var_char,iframe,flag,interpdata)

  implicit none
  integer, intent(in)     :: np,iframe
  real(sp),intent(in)     :: x(np),y(np)
  integer, intent(in)     :: ip(np),jp(np)
  integer, intent(inout)  :: istatus(np)
  character(len=*)        :: var_char
  integer, intent(in)     :: flag ! Variable flag if vector or scalar quantity to be interpolated
                                     ! scalar: = 0; vector: = 1.
  real(sp),intent(inout)  :: interpdata(np)    ! data from bilinear interpolation

  !----------------------------
  integer                 :: i,j,olvl
  integer                 :: n, iv, jv
  real(sp), pointer       :: field2d(:,:)
  real(sp), allocatable   :: fcell(:,:), xcell(:,:),ycell(:,:)
  real(sp)                :: dx1,dx2,dpx1,dpx2,xtemp,ytemp,ftemp,ew_val1,ew_val2,ytemp1,ytemp2

!  write(*,*) 'NOT in use!'
  
  ! Obtain forcing data at time 'iframe' with a pointer 'field2d'
  if (trim(var_char)/='phyto' .AND. trim(var_char)/='zoopl1') call get_forcing_f2(trim(var_char),field2d,iframe)
 
  ! Initialize interpolation value
  interpdata = 0.0_sp

 PARTICLE: DO n = 1, np ! particle loop
    
  olvl = mask(ip(n),jp(n))

  ACTIVE: IF(istatus(n) < 1) then
      cycle
  else

  select case(flag)
  !==========
  case(0) ! interp scalar quantity at ocean level center (e.g. temperature)
  !==========

  ! Before interpolationm find neareast [iv,jv] of 4 surrounding vector grid points
  call search4corners(x(n),y(n),ip(n),jp(n),iv,jv)
  
  ! check if indices are out of bounds; if so, change particle status to -1 (exited)
  if(iv<1 .OR. iv>nx .OR. jv<1 .OR. jv>ny) then
    istatus(n) = -1
    write(*,*) "Indices are Out of Bounds in 'bilinear_interp'. Change status of Particle#",n,'to EXITED!'
    cycle
  endif
  
  ! Next, obtain grid locations and values of 4 points of the cell
  if(jv == ny) jv = ny-1 ! particle at inner closed boundary, move it to ny-1

  ! Allocate arrays
  allocate(fcell(iv:iv+1, jv:jv+1))  
  allocate(xcell(iv:iv+1, jv:jv+1))
  allocate(ycell(iv:iv+1, jv:jv+1))

 ! Deal with boudary and non-boundary cases
  IF (iv == nx ) then
   DO j = jv, jv+1
     xcell(iv,  j) = slon(nx,j)
     xcell(iv+1,j) = slon(1, j)
     ycell(iv,  j) = slat(nx,j)
     ycell(iv+1,j) = slat(1, j)
     fcell(iv,  j) = field2d(nx,j)
     fcell(iv+1,j) = field2d(1, j)  
    ENDDO
  ELSE  
    DO j = jv, jv+1
     DO i = iv, iv+1
         xcell(i,j) = slon(i,j)
         ycell(i,j) = slat(i,j)
         fcell(i,j) = field2d(i,j)
     ENDDO
    ENDDO
  ENDIF 

  ! Treatment of particles when too close to north pole
   if(abs(y(n))> NorPole) then
     interpdata(n) = sum(fcell)*0.25_sp
     deallocate(xcell,ycell,fcell)
     cycle
   endif


  ! Interpolation in the horizontal directions (x-dir & y-dir)

    ! calculate differences
    dx1 = xcell(iv+1,jv) - xcell(iv,jv) ! difference in two x-dir corners
    if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
    if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

     dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
     if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
     if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

     if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases, rotate grid

         xtemp            = xcell(iv  ,jv  )
         xcell(iv  ,jv  ) = xcell(iv  ,jv+1)
         xcell(iv  ,jv+1) = xcell(iv+1,jv+1)
         xcell(iv+1,jv+1) = xcell(iv+1,jv  )
         xcell(iv+1,jv  ) = xtemp

         ytemp            = ycell(iv  ,jv  )
         ycell(iv  ,jv  ) = ycell(iv  ,jv+1)
         ycell(iv  ,jv+1) = ycell(iv+1,jv+1)
         ycell(iv+1,jv+1) = ycell(iv+1,jv  )
         ycell(iv+1,jv  ) = ytemp

         ftemp            = fcell(iv  ,jv  )
         fcell(iv  ,jv  ) = fcell(iv  ,jv+1)
         fcell(iv  ,jv+1) = fcell(iv+1,jv+1)
         fcell(iv+1,jv+1) = fcell(iv+1,jv  )
         fcell(iv+1,jv  ) = ftemp

         ! re-calculate differences in the rotated grid
         dx1 = xcell(iv+1,jv) - xcell(iv,jv) ! difference in two x-dir corners
         if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
         if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

         dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
         if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
         if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

         if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! if still not working, using average of 4 
            interpdata(n)= sum(fcell)*0.25_sp
            deallocate(xcell,ycell,fcell)
            cycle
         endif

      endif     
   
     dpx1 = xcell(iv+1,jv) - x(n)  ! difference between particle x-location to a corner  
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv,jv)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-1 by interp in x-dir 
     ew_val1 = fcell(iv,jv)*dpx1/dx1+fcell(iv+1,jv)*dpx2/dx1

     ! calculate intermediate longitudes (ytemp1)
     ytemp1  = ycell(iv,jv)*dpx1/dx1+ycell(iv+1,jv)*dpx2/dx1

     dpx1 = xcell(iv+1,jv+1) - x(n)  ! difference between particle x-location to a corner        
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv,jv+1)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-2 by interp in x-dir 
     ew_val2 = fcell(iv,jv+1)*dpx1/dx2 + fcell(iv+1,jv+1)*dpx2/dx2

     ! calculate intermediate longitudes (ytemp2)
     ytemp2  = ycell(iv,jv+1)*dpx1/dx2 + ycell(iv+1,jv+1)*dpx2/dx2

     ! calculate final interpolation value
     if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= sum(fcell)*0.25_sp
        deallocate(xcell,ycell,fcell)
        cycle
     else
        ! interpolate in y-dir
        interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                      + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2

      ! constrain final interpolated value using max and min values of fcell
       interpdata(n) = max(interpdata(n),minval(fcell))
       interpdata(n) = min(interpdata(n),maxval(fcell))

        ! deallocate arrays
        deallocate(xcell,ycell,fcell)
     endif

  !===========
  case(1) ! interp vector quantity at ocean level center (u & v)
  !===========
  
  ! Assign indices
  iv = ip(n); jv = jp(n)
 
  ! check if indices are out of bounds; if so, change particle status to -1 (exited)
  if(iv<1 .OR. iv>nx .OR. jv<1 .OR. jv>ny) then
    istatus(n) = -1
    write(*,*) "Indices are Out of Bounds in 'bilinear_interp'. Change status of Particle#",n,'to EXITED!'
    cycle
  endif
  
  ! Next, obtain grid locations and values of 4 points of the cell

  ! boundary case in y-dir
  if(jv == 1) then 
     istatus(n) = -1 ! set partcile status as 'EXITED'
     cycle
  endif

  ! Allocate arrays
  allocate(fcell(iv-1:iv, jv-1:jv))  
  allocate(xcell(iv-1:iv, jv-1:jv))
  allocate(ycell(iv-1:iv, jv-1:jv))

 ! Deal with boudary and non-boundary cases
  IF (iv == 1 ) then
   DO j = jv-1, jv
     xcell(iv-1,j) = vlon(nx,j)
     xcell(iv,  j) = vlon(1, j)
     ycell(iv-1,j) = vlat(nx,j)
     ycell(iv,  j) = vlat(1, j)
     fcell(iv-1,j) = field2d(nx,j)
     fcell(iv,  j) = field2d(1, j)
    ENDDO
  ELSE
    DO j = jv-1, jv
     DO i = iv-1, iv
         xcell(i,j) = vlon(i,j)
         ycell(i,j) = vlat(i,j)
         fcell(i,j) = field2d(i,j)
     ENDDO
    ENDDO
  ENDIF  

  ! Treatment of particles when too close to north pole
   if(abs(y(n))> NorPole) then
     interpdata(n) = sum(fcell)*0.25_sp
     deallocate(xcell,ycell,fcell)
     cycle
   endif


  ! Interpolation in the horizontal directions (x-dir & y-dir)

    ! calculate differences
    dx1 = xcell(iv,jv-1) - xcell(iv-1,jv-1) ! difference in two x-dir corners
    if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
    if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

     dx2 = xcell(iv,jv) - xcell(iv-1,jv)
     if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
     if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

     if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases, rotate grid

         xtemp            = xcell(iv-1,jv-1)
         xcell(iv-1,jv-1) = xcell(iv-1,jv  )
         xcell(iv-1,jv  ) = xcell(iv  ,jv  )
         xcell(iv  ,jv  ) = xcell(iv  ,jv-1)
         xcell(iv  ,jv-1) = xtemp

         ytemp            = ycell(iv-1,jv-1)
         ycell(iv-1,jv-1) = ycell(iv-1,jv  )
         ycell(iv-1,jv  ) = ycell(iv  ,jv  )
         ycell(iv  ,jv  ) = ycell(iv  ,jv-1)
         ycell(iv  ,jv-1) = ytemp
         
         ftemp            = fcell(iv-1,jv-1)
         fcell(iv-1,jv-1) = fcell(iv-1,jv  )
         fcell(iv-1,jv  ) = fcell(iv  ,jv  )
         fcell(iv  ,jv  ) = fcell(iv  ,jv-1)
         fcell(iv  ,jv-1) = ftemp

         ! re-calculate differences in the rotated grid
         dx1 = xcell(iv,jv-1) - xcell(iv-1,jv-1) ! difference in two x-dir corners
         if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
         if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

         dx2 = xcell(iv,jv) - xcell(iv-1,jv)
         if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
         if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

         if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! if still not working, using average of 4 
            interpdata(n)= sum(fcell)*0.25_sp
            deallocate(xcell,ycell,fcell)
            cycle
         endif

      endif     
   
     dpx1 = xcell(iv,jv-1) - x(n)  ! difference between particle x-location to a corner  
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv-1,jv-1)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-1 by interp in x-dir 
     ew_val1 = fcell(iv-1,jv-1)*dpx1/dx1+fcell(iv,jv-1)*dpx2/dx1

     ! calculate intermediate longitudes (ytemp1)
     ytemp1  = ycell(iv-1,jv-1)*dpx1/dx1+ycell(iv,jv-1)*dpx2/dx1

     dpx1 = xcell(iv,jv) - x(n)  ! difference between particle x-location to a corner        
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv-1,jv)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-2 by interp in x-dir 
     ew_val2 = fcell(iv-1,jv)*dpx1/dx2 + fcell(iv,jv)*dpx2/dx2

     ! calculate intermediate longitudes (ytemp2)
     ytemp2  = ycell(iv-1,jv)*dpx1/dx2 + ycell(iv,jv)*dpx2/dx2

     ! calculate final interpolation value
     if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= sum(fcell)*0.25_sp
        deallocate(xcell,ycell,fcell)
        cycle
     else
        ! interpolate in y-dir
        interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                      + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2

      ! constrain final interpolated value using max and min values of fcell
       interpdata(n) = max(interpdata(n),minval(fcell))
       interpdata(n) = min(interpdata(n),maxval(fcell))

        ! deallocate arrays
        deallocate(xcell,ycell,fcell)
     endif

  !===========
  case default
  !===========
    write(*,*) "Error in 'bilinear_interp': A wrong flag value: Scalar:flag=0; vector:flag=1; ww:flag=2."
  END select

  ENDIF ACTIVE
 ENDDO PARTICLE ! end particle loop

end subroutine bilinear_interp

!---------------------------------------------------------------
subroutine trilinear_interp(np,x,y,z,ip,jp,kp,istatus,var_char,iframe,flag,interpdata)

  implicit none
  integer, intent(in)       :: np,iframe
  real(sp),intent(in)       :: x(np),y(np),z(np)    ! particle locations
  integer, intent(in)       :: ip(np),jp(np),kp(np) ! scalar indices of particle
  integer, intent(inout)    :: istatus(np)          ! particle status can change
  character(len=*)          :: var_char
  integer, intent(in)       :: flag ! Variable flag if vector or scalar quantity to be interpolated
                                   ! scalar: = 0; vector: = 1, ww = 2.
  real(sp),intent(inout)    :: interpdata(np)    ! data after trilinear interpolation

  !----------------------------
  integer                   :: i,j,k
  integer                   :: iv, jv, kv  ! vector indices of particles
  integer                   :: n,olvl
  real(sp), pointer         :: field3d(:,:,:)
  real(sp)      :: coef1, coef2  ! coeffients for linear interpolation in the vertical (z-dir)
 ! variables for the particle-containing cell
  real(sp), allocatable  :: cell(:,:,:), xcell(:,:),ycell(:,:),fcell(:,:)
  real(sp)      :: dx1, dx2, dpx1, dpx2, xtemp,ytemp,ftemp
  real(sp)      :: ew_val1, ew_val2  ! intermediate values after interpolation in East-West (x-dir)
  real(sp)      :: ytemp1, ytemp2    ! intermediate y-positions

  ! Obtain forcing data at time 'iframe' with a pointer 'field3d'
  call get_forcing_f3(trim(var_char),field3d,iframe)

  ! Initialize interpolation value
  interpdata = 0.0_sp

 PARTICLE: DO n = 1, np ! particle loop

  olvl = mask(ip(n),jp(n)) ! ocean levels of particle-containing cell

 ACTIVE: if(istatus(n) < 1 ) then ! only do interpolation for active partciles
       cycle
  else

  select case(flag)

 !=======
  case(0)  ! interp scalar quantity at ocean level center (temperature, diatoms, flagellates)
 !=======

   ! Before interp, find nearest [iv,jv] of 4 surounding vector grid points
    call search4corners(x(n),y(n),ip(n),jp(n),iv,jv)
    kv = kp(n)

   ! Check if indices are out of bounds; if so, change particle status to -1 (i.e.EXITED)
   if( iv<1 .OR. iv>nx .OR. jv<1 .OR. jv>ny .OR. kv<1 .OR. kv>nz) then
       istatus(n) = -1
       write(*,*) "Indices are Out of Bounds in 'Trilinear_interpolation'." 
       write(*,*) 'Change status of Particle # ',n, ' to EXITED!'
       cycle
   endif

   ! Next, obtain grid locations and values of 8 points of the cubic
    if( jv == ny )  jv = ny-1 ! particle at inner closed boundary, move it to ny-1 
    
   ! Allocate arrays
    allocate(  cell(iv:iv+1, jv:jv+1, kv-1:kv) ) ! 8-point cubic containing particle
    allocate( xcell(iv:iv+1, jv:jv+1) )
    allocate( ycell(iv:iv+1, jv:jv+1) )
    allocate( fcell(iv:iv+1, jv:jv+1) )

    ! Deal with both boundary and non-boundary cases
    IF( kv == 1 .AND. iv == nx ) THEN   ! boundary case-1
       
       DO j = jv, jv+1
          xcell(iv,  j) = slon(nx,j)
          xcell(iv+1,j) = slon(1, j)
          ycell(iv,  j) = slat(nx,j)
          ycell(iv+1,j) = slat(1, j)
       ENDDO

       DO k = kv-1, kv
          DO j = jv, jv+1
             cell(iv,  j,k) = field3d(nx,j,1)
             cell(iv+1,j,k) = field3d(1, j,1)
          ENDDO
       ENDDO

    ELSEIF( kv == 1 .AND. iv < nx ) THEN ! boundary case-2

       DO j = jv, jv+1
         DO i = iv, iv+1
            xcell(i,j) = slon(i,j)
            ycell(i,j) = slat(i,j)
         ENDDO
       ENDDO

       DO k = kv-1, kv
         DO j = jv, jv+1
           DO i = iv, iv+1
              cell(i,j,k) = field3d(i,j,1)
           ENDDO
         ENDDO
       ENDDO

    ELSEIF( kv > 1 .AND. iv == nx ) THEN ! boundary case-3

       DO j = jv, jv+1
          xcell(iv,  j) = slon(nx,j)
          xcell(iv+1,j) = slon(1, j)
          ycell(iv,  j) = slat(nx,j)
          ycell(iv+1,j) = slat(1, j)
       ENDDO

       DO k = kv-1, kv
         DO j = jv, jv+1
            cell(iv,  j,k) = field3d(nx,j,k)
            cell(iv+1,j,k) = field3d(1, j,k)
         ENDDO
       ENDDO

    ELSEIF( kv > 1 .AND. iv < nx ) THEN ! non-boundary normal case

       DO j = jv, jv+1
         DO i = iv, iv+1
            xcell(i,j) = slon(i,j)
            ycell(i,j) = slat(i,j)
         ENDDO 
       ENDDO

      DO k = kv-1, kv
        DO j = jv, jv+1
         DO i = iv, iv+1
            cell(i,j,k) = field3d(i,j,k)
         ENDDO
       ENDDO
     ENDDO

    ELSE
        write(*,*) "Error in 'trilinear_interp': iv or kv indices NOT correct!"    
        write(*,*) "Error Particle #", n
    ENDIF  

   ! Calculate interp coefficients in the vertical (z-dir)
    if( kv == 1 ) then ! particle is located above center of top layer
         coef1 = 0.5_sp
         coef2 = 0.5_sp
    else
         coef1 = ( zt(kv) - z(n) )  /( zt(kv) - zt(kv-1) )
         coef2 = ( z(n) - zt(kv-1) )/( zt(kv) - zt(kv-1) )
    endif 

   ! Interpolate in the vertical direction
    fcell(:,:) = coef1*cell(:,:,kv-1) + coef2*cell(:,:,kv)           

   ! Treatment of particles when they are close to north pole
   ! Interpolation value is an average of 4 neighboring points
   if(abs(y(n))> NorPole) then
     interpdata(n) = sum(fcell)*0.25_sp
     deallocate(cell,xcell,ycell,fcell)
     cycle
   endif 


  ! Interpolation in the horizontal directions (x-dir & y-dir)

    ! calculate differences
    dx1 = xcell(iv+1,jv) - xcell(iv,jv) ! difference in two x-dir corners
    if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
    if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp
      
     dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
     if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
     if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

     if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases, rotate grid
         
         xtemp            = xcell(iv  ,jv  )
         xcell(iv  ,jv  ) = xcell(iv  ,jv+1)
         xcell(iv  ,jv+1) = xcell(iv+1,jv+1)
         xcell(iv+1,jv+1) = xcell(iv+1,jv  )
         xcell(iv+1,jv  ) = xtemp
 
         ytemp            = ycell(iv  ,jv  )
         ycell(iv  ,jv  ) = ycell(iv  ,jv+1)
         ycell(iv  ,jv+1) = ycell(iv+1,jv+1)
         ycell(iv+1,jv+1) = ycell(iv+1,jv  )
         ycell(iv+1,jv  ) = ytemp

         ftemp            = fcell(iv  ,jv  )
         fcell(iv  ,jv  ) = fcell(iv  ,jv+1)
         fcell(iv  ,jv+1) = fcell(iv+1,jv+1)
         fcell(iv+1,jv+1) = fcell(iv+1,jv  )
         fcell(iv+1,jv  ) = ftemp
     
         ! re-calculate differences in the rotated grid
         dx1 = xcell(iv+1,jv) - xcell(iv,jv) ! difference in two x-dir corners
         if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
         if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

         dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
         if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
         if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

         if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! if still not working, using average of 4 
            interpdata(n)= sum(fcell)*0.25_sp
            deallocate(cell,xcell,ycell,fcell)
            cycle
         endif
    
     endif

     dpx1 = xcell(iv+1,jv) - x(n)  ! difference between particle x-location to a corner        
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv,jv) 
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-1 by interp in x-dir 
     ew_val1 = fcell(iv,jv)*dpx1/dx1+fcell(iv+1,jv)*dpx2/dx1

     ! calculate intermediate longitudes (ytemp1)
     ytemp1  = ycell(iv,jv)*dpx1/dx1+ycell(iv+1,jv)*dpx2/dx1

     dpx1 = xcell(iv+1,jv+1) - x(n)  ! difference between particle x-location to a corner        
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv,jv+1)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-2 by interp in x-dir 
     ew_val2 = fcell(iv,jv+1)*dpx1/dx2 + fcell(iv+1,jv+1)*dpx2/dx2

     ! calculate intermediate longitudes (ytemp2)
     ytemp2  = ycell(iv,jv+1)*dpx1/dx2 + ycell(iv+1,jv+1)*dpx2/dx2

     ! calculate final interpolation value
     if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= sum(fcell)*0.25_sp
        deallocate(cell,xcell,ycell,fcell)
        cycle
     else
        ! interpolate in y-dir
        interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                      + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2

      ! constrain final interpolated value using max and min values of fcell
       interpdata(n) = max(interpdata(n),minval(fcell)) 
       interpdata(n) = min(interpdata(n),maxval(fcell))

        ! deallocate arrays
        deallocate(cell,xcell,ycell,fcell)
     endif

!========
  case(1) ! interp vector quantity at ocean level center (u & v)
!========
   
   ! Assign indices
   iv = ip(n); jv = jp(n); kv = kp(n)

   ! Check if indices are out of bounds; if so, change particle status to -1 (i.e. exited)
   if( iv<1 .OR. iv>nx .OR. jv<1 .OR. jv>ny .OR. kv<1 .OR. kv>nz) then
       istatus(n) = -1
       write(*,*) "Indices are Out of Bounds in 'Trilinear_interpolation'."
       write(*,*) 'Change status of Particle # ',n, ' to EXITED!'
       cycle
   endif

   ! Next, obtain grid locations and values of 8 points of the cubic
    
    ! boundary case in y-dir  
    if( jv == 1 )  then ! particle at outside open boundary
      istatus(n) = -1   ! set particle status as 'exited'
    !  write(*,*) 'Particle # ',n,' exits computational domain!'
      cycle
    endif
    
    ! Allocate arrays
    allocate(  cell(iv-1:iv, jv-1:jv, kv-1:kv) ) ! 8-point cubic containing particle
    allocate( xcell(iv-1:iv, jv-1:jv) )
    allocate( ycell(iv-1:iv, jv-1:jv) )
    allocate( fcell(iv-1:iv, jv-1:jv) )

    ! Deal with both boundary and non-boundary cases: always periodic bounary in x-dir
    IF( kv == 1 .AND. iv == 1 ) THEN   ! boundary case-1
    
        DO j = jv-1, jv
              xcell(iv-1,j) = vlon(nx,j)
              xcell(iv,  j) = vlon(1, j)
              ycell(iv-1,j) = vlat(nx,j)
              ycell(iv,  j) = vlat(1, j)
        ENDDO
  
        DO k = kv-1, kv
           DO j = jv-1, jv
              cell(iv-1,j,k) = field3d(nx,j,1)
              cell(iv,  j,k) = field3d(1, j,1)
           ENDDO
        ENDDO

    ELSEIF( kv == 1 .AND. iv > 1) THEN ! boundary case-2

        DO j = jv-1, jv
           DO i = iv-1, iv
              xcell(i,j) = vlon(i,j)
              ycell(i,j) = vlat(i,j)
           ENDDO
        ENDDO

        DO k = kv-1, kv
           DO j = jv-1, jv
              DO i = iv-1, iv
                 cell(i,j,k) = field3d(i,j,1)
              ENDDO
           ENDDO
        ENDDO 

    ELSEIF( kv > 1 .AND. iv == 1 ) THEN ! boundary case-3

        DO j = jv-1, jv
           xcell(iv-1,j) = vlon(nx,j)
           ycell(iv,  j) = vlat(1, j)
        ENDDO
      
        DO k = kv-1, kv
           DO j = jv-1, jv
              cell(iv-1,j,k) = field3d(nx,j,k)
              cell(iv,  j,k) = field3d(1, j,k)
           ENDDO
        ENDDO

    ELSEIF( kv > 1 .AND. iv > 1) THEN ! non-boundary normal case

        DO j = jv-1, jv
           DO i = iv-1, iv
              xcell(i,j) = vlon(i,j)
              ycell(i,j) = vlat(i,j)
           ENDDO
        ENDDO
 
        DO k = kv-1, kv
           DO j = jv-1, jv
              DO i = iv-1, iv
                 cell(i,j,k) = field3d(i,j,k)
              ENDDO
           ENDDO
        ENDDO

    ELSE
        write(*,*) "Error in 'trilinear_interp': iv or kv indices NOT correct!"
        write(*,*) "Error Particle #", n
    ENDIF

   ! Calculate interp coefficients in the vertical (z-dir)
    if( kv == 1 ) then ! particle is located above center of top layer
         coef1 = 0.5_sp
         coef2 = 0.5_sp
    else
         coef1 = ( zt(kv) - z(n) )  /( zt(kv) - zt(kv-1) )
         coef2 = ( z(n) - zt(kv-1) )/( zt(kv) - zt(kv-1) )
    endif

   ! Interpolate in the vertical direction
    fcell(:,:) = coef1*cell(:,:,kv-1) + coef2*cell(:,:,kv)

   ! Treatment of particles when they are close to north pole
   ! Interpolation value is an average of 4 neighboring points
   if(abs(y(n))> NorPole) then
     interpdata(n) = sum(fcell)*0.25_sp
     deallocate(cell,xcell,ycell,fcell)
     cycle
   endif

   ! Interpolation in the horizontal directions (x-dir & y-dir)

    ! calculate differences
    dx1 = xcell(iv,jv-1) - xcell(iv-1,jv-1) ! difference in two x-dir corners
    if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
    if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

     dx2 = xcell(iv,jv) - xcell(iv-1,jv)
     if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
     if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

     if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases, rotate grid

         xtemp            = xcell(iv-1,jv-1)
         xcell(iv-1,jv-1) = xcell(iv-1,jv  )
         xcell(iv-1,jv  ) = xcell(iv  ,jv  )
         xcell(iv  ,jv  ) = xcell(iv  ,jv-1)
         xcell(iv  ,jv-1) = xtemp

         ytemp            = ycell(iv-1,jv-1)
         ycell(iv-1,jv-1) = ycell(iv-1,jv  )
         ycell(iv-1,jv  ) = ycell(iv  ,jv  )
         ycell(iv  ,jv  ) = ycell(iv  ,jv-1)
         ycell(iv  ,jv-1) = ytemp

         ftemp            = fcell(iv-1,jv-1)
         fcell(iv-1,jv-1) = fcell(iv-1,jv  )
         fcell(iv-1,jv  ) = fcell(iv  ,jv  )
         fcell(iv  ,jv  ) = fcell(iv  ,jv-1)
         fcell(iv  ,jv-1) = ftemp
     
         ! re-calculate differences in the rotated grid
         dx1 = xcell(iv,jv-1) - xcell(iv-1,jv-1) ! difference in two x-dir corners
         if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
         if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

         dx2 = xcell(iv,jv) - xcell(iv-1,jv)
         if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
         if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

         if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! if still not working, using average of 4
            interpdata(n)= sum(fcell)*0.25_sp
            deallocate(cell,xcell,ycell,fcell)
            cycle
         endif

     endif

     dpx1 = xcell(iv,jv-1) - x(n)  ! difference between particle x-location to a corner
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv-1,jv-1)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-1 by interp in x-dir
     ew_val1 = fcell(iv-1,jv-1)*dpx1/dx1+fcell(iv,jv-1)*dpx2/dx1

     ! calculate intermediate longitudes (ytemp1)
     ytemp1  = ycell(iv-1,jv-1)*dpx1/dx1+ycell(iv,jv-1)*dpx2/dx1

     dpx1 = xcell(iv,jv) - x(n)  ! difference between particle x-location to a corner
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv-1,jv)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-2 by interp in x-dir
     ew_val2 = fcell(iv-1,jv)*dpx1/dx2 + fcell(iv,jv)*dpx2/dx2

     ! calculate intermediate longitudes (ytemp2)
     ytemp2  = ycell(iv-1,jv)*dpx1/dx2 + ycell(iv,jv)*dpx2/dx2

     ! calculate final interpolation value
     if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= sum(fcell)*0.25_sp
        deallocate(cell,xcell,ycell,fcell)
        cycle
     else
       ! interpolate in y-dir
       interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                     + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2

      ! constrain final interpolated value using max and min values of fcell
       interpdata(n) = max(interpdata(n),minval(fcell))
       interpdata(n) = min(interpdata(n),maxval(fcell))

       ! deallocate arrays
       deallocate(cell,xcell,ycell,fcell)
     endif

 !========
  case(2)  ! interp scalar quantity at ocean level bottom (ww)
 !========

    ! Before interp, find nearest [iv,jv] of 4 surounding vector grid points
   call search4corners(x(n),y(n),ip(n),jp(n),iv,jv)

   ! Find kv using ocean level bottom depth 
   call find_kw(z(n),kv)
  
   ! Check if indices are out of bounds; if so, change particle status to -1 (i.e.EXITED)
   if( iv<1 .OR. iv>nx .OR. jv<1 .OR. jv>ny .OR. kv<1 .OR. kv>nz) then
       istatus(n) = -1
       write(*,*) "Indices are Out of Bounds in 'Trilinear_interpolation'."
       write(*,*) 'Change status of Particle # ',n, ' to EXITED!'
       cycle
   endif

   ! Next, obtain grid locations and values of 8 points of the cubic
    if( jv == ny )  jv = ny-1 ! particle at inner closed boundary, move it to ny-1 

   ! Allocate arrays
    allocate(  cell(iv:iv+1, jv:jv+1, kv-1:kv) ) ! 8-point cubic containing particle
    allocate( xcell(iv:iv+1, jv:jv+1) )
    allocate( ycell(iv:iv+1, jv:jv+1) )
    allocate( fcell(iv:iv+1, jv:jv+1) )
    
    ! Deal with both boundary and non-boundary cases
    IF( kv == 1 .AND. iv == nx ) THEN   ! boundary case-1

        DO j = jv, jv+1
           xcell(iv,  j) = slon(nx,j)
           xcell(iv+1,j) = slon(1, j)
           ycell(iv,  j) = slat(nx,j)
           ycell(iv+1,j) = slat(1, j)
        ENDDO

        DO k = kv-1, kv
           DO j = jv, jv+1
              cell(iv,  j,k) = field3d(iv,  j,1)
              cell(iv+1,j,k) = field3d(iv+1,j,1)
           ENDDO
        ENDDO

    ELSEIF( kv == 1 .AND. iv < nx ) THEN ! boundary case-2

       DO j = jv, jv+1
          DO i = iv, iv+1
             xcell(i,j) = slon(i,j)
             ycell(i,j) = slat(i,j)
          ENDDO
       ENDDO

       DO k = kv-1, kv
          DO j = jv, jv+1
             DO i = iv, iv+1
                cell(i,j,k) = field3d(i,j,1)
             ENDDO
          ENDDO
       ENDDO

    ELSEIF( kv > 1 .AND. iv == nx ) THEN ! boundary case-3

       DO j = jv, jv+1
          xcell(iv,  j) = slon(nx,j)
          xcell(iv+1,j) = slon(1, j)
          ycell(iv,  j) = slat(nx,j)
          ycell(iv+1,j) = slat(1, j)
       ENDDO

       DO k = kv-1, kv
          DO j = jv, jv+1
             cell(iv,  j,k) = field3d(nx,j,k)
             cell(iv+1,j,k) = field3d(1, j,k)
          ENDDO
       ENDDO

    ELSEIF( kv > 1 .AND. iv < nx) THEN ! non-boundary normal case

       DO j = jv, jv+1
          DO i = iv, iv+1
             xcell(i,j) = slon(i,j)
             ycell(i,j) = slat(i,j)
          ENDDO
       ENDDO

       DO k = kv-1, kv
          DO j = jv, jv+1
             DO i = iv, iv+1
                cell(i,j,k) = field3d(i,j,k)
             ENDDO
          ENDDO
       ENDDO

    ELSE
        write(*,*) "Error in 'trilinear_interp': iv or kv indices NOT correct!"
        write(*,*) "Error Particle #", n
    ENDIF

   ! Calculate interp coefficients in the vertical (z-dir)
    if( kv == 1 ) then ! particle is located above center of top layer
         coef1 = 0.5_sp
         coef2 = 0.5_sp
    else
         coef1 = ( zw(kv) - z(n) )  /( zw(kv) - zw(kv-1) )
         coef2 = ( z(n) - zw(kv-1) )/( zw(kv) - zw(kv-1) )
    endif

   ! Interpolate in the vertical direction
    fcell(:,:) = coef1*cell(:,:,kv-1) + coef2*cell(:,:,kv)

   ! Treatment of particles when they are close to north pole
   ! Interpolation value is an average of 4 neighboring points
   if(abs(y(n))> NorPole) then
     interpdata(n) = sum(fcell)*0.25_sp
     deallocate(cell,xcell,ycell,fcell)
     cycle
   endif

  ! Interpolation in the horizontal directions (x-dir & y-dir)
    ! calculate differences
    dx1 = xcell(iv+1,jv) - xcell(iv,jv) ! difference in two x-dir corners
    if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
    if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

     dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
     if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
     if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

     if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases, rotate grid
     
         xtemp            = xcell(iv  ,jv  )
         xcell(iv  ,jv  ) = xcell(iv  ,jv+1)
         xcell(iv  ,jv+1) = xcell(iv+1,jv+1)
         xcell(iv+1,jv+1) = xcell(iv+1,jv  )
         xcell(iv+1,jv  ) = xtemp

         ytemp            = ycell(iv  ,jv  )
         ycell(iv  ,jv  ) = ycell(iv  ,jv+1)
         ycell(iv  ,jv+1) = ycell(iv+1,jv+1)
         ycell(iv+1,jv+1) = ycell(iv+1,jv  )
         ycell(iv+1,jv  ) = ytemp

         ftemp            = fcell(iv  ,jv  )
         fcell(iv  ,jv  ) = fcell(iv  ,jv+1)
         fcell(iv  ,jv+1) = fcell(iv+1,jv+1)
         fcell(iv+1,jv+1) = fcell(iv+1,jv  )
         fcell(iv+1,jv  ) = ftemp

         ! re-calculate differences in the rotated grid
         dx1 = xcell(iv+1,jv) - xcell(iv,jv) ! difference in two x-dir corners
         if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
         if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

         dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
         if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
         if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

         if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! if still not working, using average of 4 grid
            interpdata(n)= sum(fcell)*0.25_sp
            deallocate(cell,xcell,ycell,fcell)
            cycle
         endif

     endif

     dpx1 = xcell(iv+1,jv) - x(n)  ! difference between particle x-location to a corner
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv,jv)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-1 by interp in x-dir
     ew_val1 = fcell(iv,jv)*dpx1/dx1+fcell(iv+1,jv)*dpx2/dx1

     ! calculate intermediate longitude (ytemp1)
     ytemp1  = ycell(iv,jv)*dpx1/dx1+ycell(iv+1,jv)*dpx2/dx1

     dpx1 = xcell(iv+1,jv+1) - x(n)  ! difference between particle x-location to a corner
     if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
     if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

     dpx2 = x(n)- xcell(iv,jv+1)
     if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
     if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

     ! calculate intermediate value-2 by interp in x-dir
     ew_val2 = fcell(iv,jv+1)*dpx1/dx2 + fcell(iv+1,jv+1)*dpx2/dx2

     ! calculate intermediate longitude (ytemp2)
     ytemp2  = ycell(iv,jv+1)*dpx1/dx2 + ycell(iv+1,jv+1)*dpx2/dx2

     if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= sum(fcell)*0.25_sp
        deallocate(cell,xcell,ycell,fcell)
        cycle
     else
       ! interpolate in y-dir
       interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                     + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2

      ! constrain final interpolated value using max and min values of fcell
       interpdata(n) = max(interpdata(n),minval(fcell))
       interpdata(n) = min(interpdata(n),maxval(fcell))

       ! deallocate arrays
       deallocate(cell,xcell,ycell,fcell)
    endif

 !===========
  case default
 !===========
    write(*,*) "Error in 'trilinear_interp': "
    write(*,*) "Given flag is wrong! Scalar: flag=0; vector: flag=1; ww: flag=2"
  End select

  endif ACTIVE

 enddo PARTICLE ! particle loop

 nullify(field3d)

end subroutine trilinear_interp


!----------------------------------------------------------------
! -- 'interp_SPM': obtain the subsurface phytoplankton or food maxima
!    for calculating food-dependent copepod development rate.
! -- When copepod migrate to layer of max food, also force them to corresponding depth
! -- zfeng 03/10/2015
! -- zfeng 11/25/2015
!---------------------------------------------------------------
subroutine interp_SPM(np,x,y,z,ip,jp,kp,istatus,var_char,iframe,flag,SPM)

 implicit none

 integer, intent(in)         :: np, iframe,flag
 real(sp),intent(in)         :: x(np),y(np)
 integer, intent(in)         :: ip(np),jp(np)
 integer, intent(inout)      :: kp(np) ! update index to the layer of max food
 real(sp),intent(inout)      :: z(np)  ! update depth to the layer of max food
 integer, intent(inout)      :: istatus(np)
 character(len=*)            :: var_char
 real(sp),intent(inout)      :: SPM(np)

 !--------------------
 integer                     :: i,j,k,kmax
 integer                     :: iv,jv,kv
 integer                     :: n, olvl
 real(sp), pointer           :: field3d(:,:,:)
 real(sp), allocatable       :: cell(:,:,:), xcell(:,:), ycell(:,:)
 real(sp)                    :: phytop(nz) ! vertical phytoplankton concentrations
 real(sp)                    :: phyto2(2,nz)
 real(sp)   :: dx1,dx2,dpx1,dpx2,ew_val1(nz),ew_val2(nz),ytemp1,ytemp2,mincell,maxcell 
! real(sp), parameter         :: LowFood = 1.0E-6_sp ! very low phyto threshold 


 ! Obtain forcing data at time 'iframe' with a pointer 'field3d'
 call get_forcing_f3(trim(var_char),field3d,iframe) 

 ! Initialize SPM (subsurface phytoplankton maximum)
 SPM = 0.0_sp

 PARTICLE: DO  n = 1, np  ! particle loop

   olvl = mask(ip(n),jp(n)) ! ocean levels of particle-containing cell

 ACTIVE: IF(istatus(n) < 1) THEN
         cycle      ! only do interpolation for active particles
   ELSE 
     
    ! before interp, find nearest [iv,jv] of 4 surrounding vector grid points
    call search4corners(x(n),y(n),ip(n),jp(n),iv,jv)
    kv = kp(n)

    ! check if indices are out of bounds; if so, change particle status to -1  
    if(iv<1 .OR. iv>nx .OR. jv<1 .OR. jv>ny .OR. kv<1 .OR. kv>nz) then
       istatus(n) = -1
       write(*,*) "Indices are Out of Bounds in 'interp_spm'!"
       write(*,*) "Change status of Particle #",n," to EXITED!"
       cycle
    endif 

    ! Next, obtain grid locations and values of 4*nz points 

    if(jv == ny) jv = ny-1 ! inner closed boundary 

    ! Allocate arrays
    allocate( cell(iv:iv+1, jv:jv+1, 1:nz)) 
    allocate(xcell(iv:iv+1, jv:jv+1))
    allocate(ycell(iv:iv+1, jv:jv+1))

    ! Deal with boundary and non-boundary cases
     IF (iv == nx)  then ! periodic boundary 
        
      DO k = 1, nz
        DO j = jv, jv+1
           cell(iv,  j, k) = field3d(nx, j, k) 
           cell(iv+1,j, k) = field3d(1,  j, k)
        ENDDO
      ENDDO 
  
      DO  j = jv, jv+1
          xcell(iv,  j) = slon(nx,j)
          xcell(iv+1,j) = slon(1, j)
          ycell(iv,  j) = slat(nx,j)
          ycell(iv+1,j) = slat(1, j)
      ENDDO

     ELSEIF(iv < nx) then ! normal case

       DO k = 1, nz
         DO j = jv, jv+1
           DO i = iv, iv+1
              cell(i,j,k) = field3d(i,j,k)
           ENDDO
         ENDDO
       ENDDO

       DO  j = jv, jv+1
           DO i = iv, iv+1
              xcell(i,j) = slon(i,j)
              ycell(i,j) = slat(i,j)
           ENDDO
       ENDDO

     ELSE 
         write(*,*) "Error in 'interp_spm': iv index is NOT correct!"
         write(*,*) "Error Particle #",n
     ENDIF

   ! calculate min & max values of the cell
     mincell = minval(cell)
     maxcell = maxval(cell)     

    ! Check if all values in cell are very small.
    ! If so, given zero
    IF(maxcell < LowFood) then
       SPM(n) = ZERO
       deallocate(cell,xcell,ycell)
       cycle
    else

   ! Special Case Treatment when particles are close to north pole.
   ! Interpolation value is an average of 4 neighboring points.
      if(abs(y(n))> NorPole) then
        phyto2 = sum(cell,1)
        phytop = sum(phyto2,1)*0.25_sp
        SPM(n) = maxval(phytop(1:ActLayInd)) ! find maximum food within active layer
        kmax   = maxloc(phytop(1:ActLayInd),DIM=1) ! find layer index of max phytop
        kp(n)  = kmax           ! update layer index
        z(n)   = zt(kmax)       ! update copepod depth to the depth of max phytop
        deallocate(cell,xcell,ycell)
        cycle
      endif

       ! Interpolate in horizontal directions, x-dir first
       dx1 = xcell(iv+1,jv) - xcell(iv,jv) 

       if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp  ! crossing 0 meridian
       if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

       dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
       if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
       if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

       if(abs(dx1)<0.0001_sp .OR. abs(dx2)<0.0001_sp) then
         phyto2 = sum(cell,1)
         phytop = sum(phyto2,1)*0.25_sp  ! average 4 horizontal grid points
         SPM(n) = maxval(phytop(1:ActLayInd))
         kmax   = maxloc(phytop(1:ActLayInd),DIM=1)
         kp(n)  = kmax
         z(n)   = zt(kmax)
         deallocate(cell,xcell,ycell)
         cycle
       endif

       dpx1 = xcell(iv+1,jv) - x(n)  ! difference between particle and a grid corner

       if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp  ! crossing 0 meridian
       if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

       dpx2 = x(n) - xcell(iv,jv)
       if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
       if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp       
      
       ! calculate intermediate value-1 by interp in x-dir
       do k = 1, nz
          ew_val1(k) = cell(iv,jv,k)*dpx1/dx1 + cell(iv+1,jv,k)*dpx2/dx1
       enddo

       ! calculate intermediate longitudes
       ytemp1 = ycell(iv,jv)*dpx1/dx1 + ycell(iv+1,jv)*dpx2/dx1

       !---------
       dpx1 = xcell(iv+1,jv+1) - x(n)  ! difference between particle and a grid corner

       if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp  ! crossing 0 meridian
       if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

       dpx2 = x(n) - xcell(iv,jv+1)
       if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
       if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

       ! calculate intermediate value-2 by interp in x-dir
       do k = 1, nz
          ew_val2(k) = cell(iv,jv+1,k)*dpx1/dx2 + cell(iv+1,jv+1,k)*dpx2/dx2
       enddo

       ! calculate intermediate longitudes
       ytemp2 = ycell(iv,jv+1)*dpx1/dx2 + ycell(iv+1,jv+1)*dpx2/dx2

       ! calculate final interpolated value
       if(abs(ytemp2-ytemp1)<0.0001_sp) then
         phyto2 = sum(cell,1)
         phytop = sum(phyto2,1)*0.25_sp  ! average 4 horizontal grid points
         SPM(n) = maxval(phytop(1:ActLayInd))
         kmax   = maxloc(phytop(1:ActLayInd),DIM=1)
         kp(n)  = kmax
         z(n)   = zt(kmax)
         deallocate(cell,xcell,ycell)
         cycle
       else 
         phytop = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2       
         
         SPM(n) = maxval(phytop(1:ActLayInd)) 
         kmax   = maxloc(phytop(1:ActLayInd),DIM=1)
         kp(n)  = kmax
         z(n)   = zt(kmax) 
         deallocate(cell,xcell,ycell)
       endif 
       
       ! constrain final interpolated value using max and min values of cell
       SPM(n) = max(SPM(n),mincell)
       SPM(n) = min(SPM(n),maxcell)             

    ENDIF

  ENDIF ACTIVE

 ENDDO PARTICLE       

end subroutine interp_spm


subroutine interp_SPM_2D(np,x,y,ip,jp,istatus,var_char,iframe,flag,SPM)
!----------------------------------------------------------------
! -- 'interp_SPM_2D': obtain the subsurface phytoplankton or food maxima
!    for calculating food-dependent copepod development rate.
! -- This is done in 2D mode, so depth information is not stored; otherwise, use 3D. 
! -- zfeng 12/09/2015
!---------------------------------------------------------------

 implicit none

 integer, intent(in)         :: np, iframe,flag
 real(sp),intent(in)         :: x(np),y(np)
 integer, intent(in)         :: ip(np),jp(np)
 integer, intent(inout)      :: istatus(np)
 character(len=*)            :: var_char
 real(sp),intent(inout)      :: SPM(np)

 !--------------------
 integer                     :: i,j,k
 integer                     :: iv,jv
 integer                     :: n, olvl
 real(sp), pointer           :: field3d(:,:,:)
 real(sp), allocatable       :: cell(:,:,:), xcell(:,:), ycell(:,:)
 real(sp)                    :: phytop(nz) ! vertical phytoplankton concentrations
 real(sp)                    :: phyto2(2,nz)
 real(sp)   :: dx1,dx2,dpx1,dpx2,ew_val1(nz),ew_val2(nz),ytemp1,ytemp2,mincell,maxcell 
! real(sp), parameter         :: LowFood = 1.0E-6_sp ! very low phyto threshold 


 ! Obtain forcing data at time 'iframe' with a pointer 'field3d'
 if (food_source == 3) then
  call get_forcing_f3('phyto',field3d,iframe) 
 elseif (food_source == 4) then
  call get_forcing_f3('zoopl1',field3d,iframe)
 else
  write(*,*) "Error in 'interp_SPM_2D: redefine 'food_source'!"
  stop 
 endif

 ! Initialize SPM (subsurface phytoplankton maximum)
 SPM = 0.0_sp

 PARTICLE: DO  n = 1, np  ! particle loop

   olvl = mask(ip(n),jp(n)) ! ocean levels of particle-containing cell

 ACTIVE: IF(istatus(n) < 1) THEN
         cycle      ! only do interpolation for active particles
   ELSE 
     
    ! before interp, find nearest [iv,jv] of 4 surrounding vector grid points
    call search4corners(x(n),y(n),ip(n),jp(n),iv,jv)

    ! check if indices are out of bounds; if so, change particle status to -1  
    if(iv<1 .OR. iv>nx .OR. jv<1 .OR. jv>ny ) then
       istatus(n) = -1
       write(*,*) "Indices are Out of Bounds in 'interp_spm_2D'!"
       write(*,*) "Change status of Particle #",n," to EXITED!"
       cycle
    endif 

    ! Next, obtain grid locations and values of 4*nz points 

    if(jv == ny) jv = ny-1 ! inner closed boundary 

    ! Allocate arrays
    allocate( cell(iv:iv+1, jv:jv+1, 1:nz)) 
    allocate(xcell(iv:iv+1, jv:jv+1))
    allocate(ycell(iv:iv+1, jv:jv+1))

    ! Deal with boundary and non-boundary cases
     IF (iv == nx)  then ! periodic boundary 
        
      DO k = 1, nz
        DO j = jv, jv+1
           cell(iv,  j, k) = field3d(nx, j, k) 
           cell(iv+1,j, k) = field3d(1,  j, k)
        ENDDO
      ENDDO 
  
      DO  j = jv, jv+1
          xcell(iv,  j) = slon(nx,j)
          xcell(iv+1,j) = slon(1, j)
          ycell(iv,  j) = slat(nx,j)
          ycell(iv+1,j) = slat(1, j)
      ENDDO

     ELSEIF(iv < nx) then ! normal case

       DO k = 1, nz
         DO j = jv, jv+1
           DO i = iv, iv+1
              cell(i,j,k) = field3d(i,j,k)
           ENDDO
         ENDDO
       ENDDO

       DO  j = jv, jv+1
           DO i = iv, iv+1
              xcell(i,j) = slon(i,j)
              ycell(i,j) = slat(i,j)
           ENDDO
       ENDDO

     ELSE 
         write(*,*) "Error in 'interp_spm_2D': iv index is NOT correct!"
         write(*,*) "Error Particle #",n
     ENDIF

   ! calculate min & max values of the cell
     mincell = minval(cell)
     maxcell = maxval(cell)     

    ! Check if all values in cell are very small.
    ! If so, given zero
    IF(maxcell < LowFood) then
       SPM(n) = ZERO
       deallocate(cell,xcell,ycell)
       cycle
    else

   ! Special Case Treatment when particles are close to north pole.
   ! Interpolation value is an average of 4 neighboring points.
      if(abs(y(n))> NorPole) then
        phyto2 = sum(cell,1)
        phytop = sum(phyto2,1)*0.25_sp
        SPM(n) = maxval(phytop(1:ActLayInd)) ! find maximum food within active layer
        deallocate(cell,xcell,ycell)
        cycle
      endif

       ! Interpolate in horizontal directions, x-dir first
       dx1 = xcell(iv+1,jv) - xcell(iv,jv) 

       if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp  ! crossing 0 meridian
       if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

       dx2 = xcell(iv+1,jv+1) - xcell(iv,jv+1)
       if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
       if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

       if(abs(dx1)<0.0001_sp .OR. abs(dx2)<0.0001_sp) then
         phyto2 = sum(cell,1)
         phytop = sum(phyto2,1)*0.25_sp  ! average 4 horizontal grid points
         SPM(n) = maxval(phytop(1:ActLayInd))
         deallocate(cell,xcell,ycell)
         cycle
       endif

       dpx1 = xcell(iv+1,jv) - x(n)  ! difference between particle and a grid corner

       if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp  ! crossing 0 meridian
       if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

       dpx2 = x(n) - xcell(iv,jv)
       if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
       if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp       
      
       ! calculate intermediate value-1 by interp in x-dir
       do k = 1, nz
          ew_val1(k) = cell(iv,jv,k)*dpx1/dx1 + cell(iv+1,jv,k)*dpx2/dx1
       enddo

       ! calculate intermediate longitudes
       ytemp1 = ycell(iv,jv)*dpx1/dx1 + ycell(iv+1,jv)*dpx2/dx1

       !---------
       dpx1 = xcell(iv+1,jv+1) - x(n)  ! difference between particle and a grid corner

       if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp  ! crossing 0 meridian
       if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

       dpx2 = x(n) - xcell(iv,jv+1)
       if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
       if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

       ! calculate intermediate value-2 by interp in x-dir
       do k = 1, nz
          ew_val2(k) = cell(iv,jv+1,k)*dpx1/dx2 + cell(iv+1,jv+1,k)*dpx2/dx2
       enddo

       ! calculate intermediate longitudes
       ytemp2 = ycell(iv,jv+1)*dpx1/dx2 + ycell(iv+1,jv+1)*dpx2/dx2

       ! calculate final interpolated value
       if(abs(ytemp2-ytemp1)<0.0001_sp) then
         phyto2 = sum(cell,1)
         phytop = sum(phyto2,1)*0.25_sp  ! average 4 horizontal grid points
         SPM(n) = maxval(phytop(1:ActLayInd))
         deallocate(cell,xcell,ycell)
         cycle
       else 
         phytop = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2       
         
         SPM(n) = maxval(phytop(1:ActLayInd)) 
         deallocate(cell,xcell,ycell)
       endif 
       
       ! constrain final interpolated value using max and min values of cell
       SPM(n) = max(SPM(n),mincell)
       SPM(n) = min(SPM(n),maxcell)             

    ENDIF

  ENDIF ACTIVE

 ENDDO PARTICLE       

end subroutine interp_spm_2D

!---------------------------------------------------
! 2-D Advection in structured grid forcing
! zfeng 12/05/2014
!----------------------------------------------------
subroutine advect2dRK4(g,deltaT,np,time)

  implicit none

  integer, intent(in)         :: np
  type(igroup), intent(inout) :: g
  real(sp), intent(in)        :: deltaT,time

  real(sp), pointer           :: x(:),y(:),h(:)
  real(sp), pointer           :: u(:),v(:),z(:),pathlength(:)
  integer , pointer           :: i0(:),j0(:),istatus0(:)
  integer,  dimension(np)     :: ip,jp,istatus
  integer                     :: i,ns,ni,n
  real(sp), dimension(np)     :: u0,u1,u2,v0,v1,v2
  real(sp), dimension(np)     :: zeta,zeta1,zeta2,pdx,pdy
  real(sp), dimension(np)     :: wu,wu1,wu2,wv,wv1,wv2
  real(sp), dimension(np)     :: pdxt,pdyt

  real(sp), dimension(np,0:mstage) :: chix,chiy !rji outofbound error

  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
!  call get_state('z',g,z)
  call get_state('i',g,i0)
  call get_state('j',g,j0)
  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('h',g,h)
  call get_state('pathlength',g,pathlength)
  call get_state('status',g,istatus0)

  istatus = 0   ! 
  where(istatus0 == 1) istatus = 1
  ip = i0; jp = j0

  ! Initialize stage functional evaluations
  pdx  = x;      pdy = y
  chix = 0.0_sp; chiy = 0.0_sp

 !--Loop over RK stages
 DO ns = 1, mstage
    !! particle position at stage ns (x,y,z)
    if(spherical==0) then
      pdx(:) = x(:) + a_rk(ns)*deltaT*chix(:,ns-1)
      pdy(:) = y(:) + a_rk(ns)*deltaT*chiy(:,ns-1)
    elseif(spherical==1) then
      pdx(:) = x(:) + a_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COSD(pdy(:))+eps)
      pdy(:) = y(:) + a_rk(ns)*deltaT*chiy(:,ns-1)/tpi
      where( pdy > 90.0_SP)     pdy = 180.0_SP - pdy
      where( pdy < -90.0_SP)    pdy = -180.0_SP - pdy
      where( pdx < 0.0_SP)      pdx = pdx + 360.0_SP
      where( pdx >= 360.0_SP)   pdx = pdx - 360.0_SP
    endif
  
  ! update indices i,j
  call update_ij(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:))

  ! calculate velocity for Stage-N using C_RK coefficients
  ! interpolate velocity to particle position
  ! Note that u & v are in cm/s, w is in m/s
  call bilinear_interp(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:),ua_char,3,1,u1)
  call bilinear_interp(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:),ua_char,4,1,u2)

  call bilinear_interp(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:),va_char,3,1,v1)
  call bilinear_interp(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:),va_char,4,1,v2)

  if(fix_dep == 0)  then  ! unfixed depth, interpolate vertical velocity 
  !   call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'ww',3,2,w1)
  !   call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'ww',4,2,w2)
  elseif(fix_dep == 1) then ! fix depth, impose zero vertical velocity
  !   w1 = 0.0_sp
  !   w2 = 0.0_sp
  else 
     write(*,*) "Error in 'advect2dRK4':fix_dep must be 0 (unfixed) or 1 (fix) !!"
     stop 
  endif

  u0 = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2
  v0 = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2

  chix(:,ns) = u0(:)*0.01_sp  ! cm/s to m/s
  chiy(:,ns) = v0(:)*0.01_sp

 ENDDO

 !--Sum stage contributions to get updated particle positions
 pdxt(:) = x(:)
 pdyt(:) = y(:)

 do ns=1,mstage
     if(spherical==0) then
      pdxt(:) = pdxt(:) + b_rk(ns)*deltaT*chix(:,ns-1)
      pdyt(:) = pdyt(:) + b_rk(ns)*deltaT*chiy(:,ns-1)
    elseif(spherical==1) then
      pdxt(:) = pdxt(:) + b_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COSD(pdyt(:))+eps)
      pdyt(:) = pdyt(:) + b_rk(ns)*deltaT*chiy(:,ns-1)/tpi
      where( pdyt > 90.0_SP)     pdyt = 180.0_SP - pdyt 
      where( pdyt < -90.0_SP)    pdyt = -180.0_SP - pdyt
      where( pdxt < 0.0_SP)      pdxt = pdxt + 360.0_SP    !
      where( pdxt >= 360.0_SP)   pdxt = pdxt - 360.0_SP    !
    endif
 enddo


 ! Calculate integrated trajectory and assign new locations

 if(bcond_type == 0) then ! absorbing boundary type     
     x(:) = x(:)*FLOAT(1 - istatus(:)) + pdxt(:)*FLOAT(istatus(:))
     y(:) = y(:)*FLOAT(1 - istatus(:)) + pdyt(:)*FLOAT(istatus(:))
     u(:) = u(:)*FLOAT(1 - istatus(:)) + chix(:,mstage)*FLOAT(istatus(:))
     v(:) = v(:)*FLOAT(1 - istatus(:)) + chiy(:,mstage)*FLOAT(istatus(:))

     if(spherical==0) then
        pathlength(:) = pathlength(:) + &
               sqrt((pdxt(:)-x(:))**2 + (pdyt(:)-y(:))**2) * FLOAT(istatus(:)) 
     elseif(spherical==1) then
        pathlength(:) = pathlength(:) + &
           tpi*sqrt( ((pdxt(:)-x(:))*COSD(pdyt(:)))**2 + (pdyt(:)-y(:))**2 ) *FLOAT(istatus(:))
     endif

 elseif(bcond_type == 1) then ! non-absorbing
     x(:) = pdxt(:)
     y(:) = pdyt(:)
     u(:) = chix(:,mstage)
     v(:) = chiy(:,mstage)

     if(spherical==0) then
        pathlength(:) = pathlength(:) + &
               sqrt((pdxt(:)-x(:))**2 + (pdyt(:)-y(:))**2)
     elseif(spherical==1) then
        pathlength(:) = pathlength(:) + &
           tpi*sqrt( ((pdxt(:)-x(:))*COSD(pdyt(:)))**2 + (pdyt(:)-y(:))**2 )
     endif

 endif 

 nullify(x,y,z,i0,j0,u,v,h,pathlength,istatus0)
 
end subroutine advect2dRK4


!----------------------------------------------------------
! 3-D Advection for structured grid forcing using RK4 scheme
! zfeng 11/15/2014
!----------------------------------------------------------
subroutine advect3dRK4(g,deltaT,np,time)

  implicit none
  integer, intent(in)         :: np
  type(igroup), intent(inout) :: g
  real(sp), intent(in)        :: deltaT,time

  real(sp), pointer           :: x(:),y(:),z(:),h(:)
  real(sp), pointer           :: u(:),v(:),w(:),pathlength(:)
  integer , pointer           :: i0(:),j0(:),k0(:),istatus0(:)
  integer,  dimension(np)     :: ip,jp,kp,istatus
  integer                     :: k,i,ns,ni,n
  real(sp), dimension(np)     :: u0,u1,u2,v0,v1,v2,w0,w1,w2
  real(sp), dimension(np)     :: zeta,zeta1,zeta2,pdx,pdy,pdz
  real(sp), dimension(np)     :: wu,wu1,wu2,wv,wv1,wv2
  real(sp), dimension(np)     :: pdxt,pdyt,pdzt

  real(sp), dimension(np,0:mstage) :: chix,chiy,chiz !rji outofbound error
  real(sp)                         :: diel,ztmp,randy(np),depth(np),max_pdepth

  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('z',g,z)
  call get_state('i',g,i0)
  call get_state('j',g,j0)
  call get_state('k',g,k0)
  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('w',g,w)
  call get_state('h',g,h)
  call get_state('pathlength',g,pathlength)
  call get_state('status',g,istatus0)

  istatus = 0   ! 
  where(istatus0 == 1) istatus = 1
  ip = i0; jp = j0; kp = k0

  ! Initialize stage functional evaluations
  pdx  = x;      pdy = y;       pdz = z
  chix = 0.0_sp; chiy = 0.0_sp; chiz = 0.0_sp

 !--Loop over RK stages
 DO ns = 1, mstage
    !! particle position at stage ns (x,y,z)
    if(spherical==0) then
      pdx(:) = x(:) + a_rk(ns)*deltaT*chix(:,ns-1)
      pdy(:) = y(:) + a_rk(ns)*deltaT*chiy(:,ns-1)
      pdz(:) = z(:) + a_rk(ns)*deltaT*chiz(:,ns-1)
    elseif(spherical==1) then
      pdx(:) = x(:) + a_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COSD(pdy(:))+eps)
      pdy(:) = y(:) + a_rk(ns)*deltaT*chiy(:,ns-1)/tpi
      pdz(:) = z(:) + a_rk(ns)*deltaT*chiz(:,ns-1)
      where( pdy > 90.0_SP)     pdy = 180.0_SP - pdy
      where( pdy < -90.0_SP)    pdy = -180.0_SP - pdy
      where( pdx < 0.0_SP)      pdx = pdx + 360.0_SP
      where( pdx >= 360.0_SP)   pdx = pdx - 360.0_SP
      where( pdz < 0.0_sp)      pdz = 0.0_sp
    endif

  ! Determine max depth a partcile can reach, or very bottom of the grid
  ! then constrain pdz
  DO n = 1, np
    max_pdepth = zw(mask(ip(n),jp(n)))
    if(pdz(n) > max_pdepth)  pdz(n) = max_pdepth
  ENDDO

  ! update indices i,j,k
  call update_ij(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:))
  call update_k(np,pdz(:),kp(:),istatus(:))

  ! calculate velocity for Stage-N using C_RK coefficients
  ! interpolate velocity to particle position
  ! Note that u & v are in cm/s, w is in m/s
  call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'u',3,1,u1)
  call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'u',4,1,u2)

  call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'v',3,1,v1)
  call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'v',4,1,v2)

  if(fix_dep == 0)  then  ! unfixed depth, interpolate vertical velocity
     call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'ww',3,2,w1)
     call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'ww',4,2,w2)
  elseif(fix_dep == 1) then ! fix depth, impose zero vertical velocity
     w1 = 0.0_sp
     w2 = 0.0_sp
  else 
     write(*,*) "Error in 'advect3dRK4':fix_dep must be 0 (unfixed) or 1 (fix) !!"
     stop 
  endif

  u0 = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2
  v0 = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2
  w0 = (1.0_sp-c_rk(ns))*w1 + c_rk(ns)*w2

  chix(:,ns) = u0(:)*0.01_sp  ! cm/s to m/s
  chiy(:,ns) = v0(:)*0.01_sp
  chiz(:,ns) = w0(:)

 ENDDO

 !--Sum stage contributions to get updated particle positions
 pdxt(:) = x(:)
 pdyt(:) = y(:)
 pdzt(:) = z(:)

 do ns=1,mstage
     if(spherical==0) then
      pdxt(:) = pdxt(:) + b_rk(ns)*deltaT*chix(:,ns-1)
      pdyt(:) = pdyt(:) + b_rk(ns)*deltaT*chiy(:,ns-1)
      pdzt(:) = pdzt(:) + b_rk(ns)*deltaT*chiz(:,ns-1)
    elseif(spherical==1) then
      pdxt(:) = pdxt(:) + b_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COSD(pdyt(:))+eps)
      pdyt(:) = pdyt(:) + b_rk(ns)*deltaT*chiy(:,ns-1)/tpi
      pdzt(:) = pdzt(:) + b_rk(ns)*deltaT*chiz(:,ns-1)
      where( pdyt > 90.0_SP)     pdyt = 180.0_SP - pdyt 
      where( pdyt < -90.0_SP)    pdyt = -180.0_SP - pdyt
      where( pdxt < 0.0_SP)      pdxt = pdxt + 360.0_SP    !
      where( pdxt >= 360.0_SP)   pdxt = pdxt - 360.0_SP    !
      where( pdzt < 0.0_sp)      pdzt = 0.0_sp
    endif
 enddo

  ! Determine max depth a partcile can reach, or very bottom of the grid
  ! then constrain pdzt
  DO n = 1, np
    if( mask(ip(n),jp(n)) == 0) then  ! land grid
        pdzt(n) = 0.0_sp
    else
        max_pdepth = zw(mask(ip(n),jp(n))) 
        if(pdzt(n) > max_pdepth)  pdzt(n) = max_pdepth

    endif
  ENDDO

 ! Calculate integrated trajectory and assign new locations

 if(bcond_type == 0) then ! absorbing boundary type     
     x(:) = x(:)*FLOAT(1 - istatus(:)) + pdxt(:)*FLOAT(istatus(:))
     y(:) = y(:)*FLOAT(1 - istatus(:)) + pdyt(:)*FLOAT(istatus(:))
     z(:) = z(:)*FLOAT(1 - istatus(:)) + pdzt(:)*FLOAT(istatus(:))
     u(:) = u(:)*FLOAT(1 - istatus(:)) + chix(:,mstage)*FLOAT(istatus(:))
     v(:) = v(:)*FLOAT(1 - istatus(:)) + chiy(:,mstage)*FLOAT(istatus(:))
     w(:) = w(:)*FLOAT(1 - istatus(:)) + chiz(:,mstage)*FLOAT(istatus(:))  

     if(spherical==0) then
        pathlength(:) = pathlength(:) + &
               sqrt((pdxt(:)-x(:))**2 + (pdyt(:)-y(:))**2) * FLOAT(istatus(:)) 
     elseif(spherical==1) then
        pathlength(:) = pathlength(:) + &
           tpi*sqrt( ((pdxt(:)-x(:))*COSD(pdyt(:)))**2 + (pdyt(:)-y(:))**2 ) *FLOAT(istatus(:))
     endif

 elseif(bcond_type == 1) then ! non-absorbing
     x(:) = pdxt(:)
     y(:) = pdyt(:)
     z(:) = pdzt(:)
     u(:) = chix(:,mstage)
     v(:) = chiy(:,mstage)
     w(:) = chiz(:,mstage)

     if(spherical==0) then
        pathlength(:) = pathlength(:) + &
               sqrt((pdxt(:)-x(:))**2 + (pdyt(:)-y(:))**2)
     elseif(spherical==1) then
        pathlength(:) = pathlength(:) + &
           tpi*sqrt( ((pdxt(:)-x(:))*COSD(pdyt(:)))**2 + (pdyt(:)-y(:))**2 )
     endif

 endif 

 nullify(x,y,z,i0,j0,k0,u,v,w,h,pathlength,istatus0)
end subroutine advect3dRK4

!----------------------------------------------------------------------
! 3-D Advection for structured grid forcing using forward Euler scheme
!----------------------------------------------------------------------
subroutine advect3DEuler(g,deltaT,np,time)

  implicit none
  integer, intent(in)         :: np
  type(igroup), intent(inout) :: g
  real(sp), intent(in)        :: deltaT,time

  real(sp), pointer           :: x(:),y(:),z(:),h(:)
  real(sp), pointer           :: u(:),v(:),w(:),pathlength(:)
  integer , pointer           :: i0(:),j0(:),k0(:),istatus0(:)
  integer,  dimension(np)     :: ip,jp,kp,istatus
  integer                     :: k,i,ns,ni,n
  real(sp), dimension(np)     :: u0,v0,w0
  real(sp), dimension(np)     :: pdx,pdy,pdz
  real(sp), dimension(np)     :: wu,wu1,wu2,wv,wv1,wv2

  real(sp), dimension(np)     :: chix,chiy,chiz
  real(sp)                    :: diel,ztmp,randy(np),depth(np),max_pdepth

  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('z',g,z)
  call get_state('i',g,i0)
  call get_state('j',g,j0)
  call get_state('k',g,k0)
  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('w',g,w)
  call get_state('h',g,h)
  call get_state('pathlength',g,pathlength)
  call get_state('status',g,istatus0)

  istatus = 0
  where(istatus0 == 1) istatus = 1
  ip = i0; jp = j0; kp = k0

  ! Initialize stage functional evaluations
  pdx  = x;      pdy = y;       pdz = z
  chix = 0.0_sp; chiy = 0.0_sp; chiz = 0.0_sp

  ! update indices i,j,k
!  call update_ij(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:))
!  call update_k(np,pdz(:),kp(:),istatus(:))

  ! calculate velocities
  ! interpolate velocities to particle position at previous time step (iframe=3)
  ! Note that u & v are in cm/s, w is in m/s
  call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'u',3,1,u0)
  call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'v',3,1,v0)

  if (fix_dep == 0) then
     call trilinear_interp(np,pdx(:),pdy(:),pdz(:),ip(:),jp(:),kp(:),istatus(:),'ww',3,2,w0)
  elseif (fix_dep == 1) then
     w0 = 0.0_sp
  else 
     write(*,*) "Error in 'advect3dEuler': fix_dep must be 0 (unfixed) or 1 (fix)!"
     stop 
  endif

  chix(:) = u0(:)*0.01_sp  ! cm/s to m/s
  chiy(:) = v0(:)*0.01_sp
  chiz(:) = w0(:)

 !--get updated particle positions using forward Euler scheme

 if(spherical==0) then
      pdx(:) = pdx(:) + deltaT*chix(:)
      pdy(:) = pdy(:) + deltaT*chiy(:)
      pdz(:) = pdz(:) + deltaT*chiz(:)
 elseif(spherical==1) then
      pdx(:) = pdx(:) + deltaT*chix(:)/(tpi*COSD(pdy(:))+eps)
      pdy(:) = pdy(:) + deltaT*chiy(:)/tpi
      pdz(:) = pdz(:) + deltaT*chiz(:)
      where( pdy > 90.0_SP)     pdy = 180.0_SP - pdy;  pdx = pdx + 180.0_sp   !
      where( pdy < -90.0_SP)    pdy = -180.0_SP - pdy; pdx = pdx + 180.0_sp
      where( pdx < 0.0_SP)      pdx = pdx + 360.0_SP    !
      where( pdx >= 360.0_SP)   pdx = pdx - 360.0_SP    !
      where( pdz < 0.0_sp)      pdz = 0.0_sp
  endif

  ! update indices i,j
  call update_ij(np,pdx(:),pdy(:),ip(:),jp(:),istatus(:))

  ! Determine max depth a partcile can reach, i.e. max depth of the grid; then constrain pdz
  ! If particle is located at a land grid, its depth is zero.
  DO n = 1, np
    if( mask(ip(n),jp(n)) == 0) then
        pdz(n) = 0.0_sp
    else
        max_pdepth = zw(mask(ip(n),jp(n)))
        if(pdz(n) >= max_pdepth)  pdz(n) = max_pdepth
    endif
  ENDDO
  
 ! update index k
 call update_k(np,pdz(:),kp(:),istatus(:))

  ! Calculate integrated trajectory and assign new locations
  if(bcond_type == 0) then ! absorbing boundary type
     x(:) = x(:)*FLOAT(1 - istatus(:)) + pdx(:)*FLOAT(istatus(:))
     y(:) = y(:)*FLOAT(1 - istatus(:)) + pdy(:)*FLOAT(istatus(:))
     z(:) = z(:)*FLOAT(1 - istatus(:)) + pdz(:)*FLOAT(istatus(:))
     u(:) = u(:)*FLOAT(1 - istatus(:)) + chix(:)*FLOAT(istatus(:))
     v(:) = v(:)*FLOAT(1 - istatus(:)) + chiy(:)*FLOAT(istatus(:))
     w(:) = w(:)*FLOAT(1 - istatus(:)) + chiz(:)*FLOAT(istatus(:))
     istatus0(:) = istatus(:) ! status may change here

     if(spherical==0) then
        pathlength(:) = pathlength(:) + &
               sqrt((pdx(:)-x(:))**2 + (pdy(:)-y(:))**2 + (pdz(:)-z(:))**2 ) * FLOAT(istatus(:))
     elseif(spherical==1) then
        pathlength(:) = pathlength(:) + &
           sqrt( tpi**2*((pdx(:)-x(:))*COSD(pdy(:)))**2 + tpi**2*(pdy(:)-y(:))**2 + (pdz(:)-z(:))**2 ) *FLOAT(istatus(:))
     endif

  elseif(bcond_type == 1) then ! non-absorbing boundary type
     x(:) = pdx(:)
     y(:) = pdy(:)
     z(:) = pdz(:)
     u(:) = chix(:)
     v(:) = chiy(:)
     w(:) = chiz(:)

     if(spherical==0) then
        pathlength(:) = pathlength(:) + &
               sqrt((pdx(:)-x(:))**2 + (pdy(:)-y(:))**2 + (pdz(:)-z(:))**2 )
     elseif(spherical==1) then
        pathlength(:) = pathlength(:) + &
           sqrt( tpi**2*((pdx(:)-x(:))*COSD(pdy(:)))**2 + tpi**2*(pdy(:)-y(:))**2 + (pdz(:)-z(:))**2 )
     endif

  endif

 nullify(x,y,z,i0,j0,k0,u,v,w,h,pathlength,istatus0)
end subroutine advect3DEuler

!==============================================================================
End Module biomas_driver
