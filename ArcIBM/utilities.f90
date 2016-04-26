Module utilities

!=======================================================================
! Fiscm Utilities 
!
! Description
!   Utilities for fiscm 
!    
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
! Fortran90 Random Number generators from the Web
!  2011/7/28       Xueming Zhu   modified function to subroutine to generate
!                                random array
!=======================================================================

  use gparms
  use netcdf

  implicit none

! default all is private.
  private
 
! public member functions:
  public gettime 
  public drawline 
  public isintriangle 
  public ncdchk 
  public get_unique_strings
  public ran_from_range 
  public unitrand 
  public normal 
  public splint
  public normal_circle
  public random_square
  public check_times
  public gcdist
  public nearest_ij

interface normal    
  module procedure normal 
  module procedure normal_zeromean 
end interface

contains

!------------------------------------------------------------------------------!
!   Return a Time String Days:Hours:Minutes:Seconds from Number of Seconds     !
!     input:                                                                   !
!       integer:  insecs [number of seconds]                                   !
!------------------------------------------------------------------------------!
function gettime(insecs) result(instring)

  implicit none
  character(len=13)   :: instring
  integer, intent(in) :: insecs 
  character(len=4)    :: s0
  character(len=2)    :: s1,s2,s3
  integer :: dtcp,htcp,mtcp,stcp
 
  dtcp = insecs/(3600*24)
  htcp = mod(insecs,(3600*24))/3600
  mtcp = mod(insecs,(3600))/60
  stcp = insecs - (dtcp*3600*24 + htcp*3600 + mtcp*60)
!! zxm  2010-10-26
  if(dtcp < 10000.)then
     write(s0,"(i4.4)")int(dtcp)
     write(s1,"(i2.2)")int(htcp)
     write(s2,"(i2.2)")int(mtcp)
     write(s3,"(i2.2)")int(stcp)
     instring = s0//":"//s1//":"//s2//":"//s3
  else
     instring = "> 10000 days"
  end if
!! zxm 2010-10-26
  return
  end function gettime

  !---------------------------------------------
  !draw a line 72 characters long repeating c 
  !optional: dump to unit iunit
  !---------------------------------------------
  subroutine drawline(c,iunit)
    character(len=1) :: c
    integer, intent(in), optional :: iunit
    character(len=72) :: line
    integer :: i
    line(1:1) = "!"
    do i=2,72
      line(i:i) = c
    end do
    if(present(iunit))then
      write(iunit,*)line
    else
      write(*,'(A72)')line
    endif
  end subroutine drawline
    
  !------------------------------------------------------------------------------
  !  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
  !  using algorithm used for scene rendering in computer graphics               |
  !  algorithm works well unless particle happens to lie in a line parallel      |
  !  to the edge of a triangle.                                                  |
  !  This can cause problems if you use a regular grid, say for idealized        |
  !  modelling and you happen to see particles right on edges or parallel to     |
  !  edges.                                                                      |
  !------------------------------------------------------------------------------
   logical function isintriangle(x0,y0,x,y)
   implicit none
   real(sp), intent(in) :: x0,y0
   real(sp), intent(in) :: x(3),y(3)
   !----------------------------------
   real(sp) :: f1,f2,f3
   !----------------------------------
   !revised by Xinyou Lin
   isintriangle = .true.

   f1 = (y0-y(1))*(x(2)-x(1)) - (x0-x(1))*(y(2)-y(1))
   f1 = f1*((y(3)-y(1))*(x(2)-x(1)) - (x(3)-x(1))*(y(2)-y(1)))  

   f2 = (y0-y(3))*(x(1)-x(3)) - (x0-x(3))*(y(1)-y(3))
   f2 = f2*((y(2)-y(3))*(x(1)-x(3)) - (x(2)-x(3))*(y(1)-y(3)))

   f3 = (y0-y(2))*(x(3)-x(2)) - (x0-x(2))*(y(3)-y(2))
   f3 =f3*((y(1)-y(2))*(x(3)-x(2)) - (x(1)-x(2))*(y(3)-y(2)))

   if(f1 <0.0_sp .or. f2 <0.0_sp .or.f3 <0.0_sp ) isintriangle = .false.
   return
   end function isintriangle



  !-----------------------------------------------
  ! runtime errors  - netcdf
  !-----------------------------------------------
subroutine ncdchk(status,message,message2)
  implicit none
  integer, intent ( in) :: status
  character(len=*), optional   :: message
  character(len=*), optional   :: message2

  if(status /= nf90_noerr) then
    if(present(message))print *, trim(message)
    if(present(message2))print *, trim(message2)
    print *, trim(nf90_strerror(status))
    stop
  end if
end subroutine ncdchk 

  !-----------------------------------------------
  ! search list of strings and return unique 
  !-----------------------------------------------
  subroutine get_unique_strings(ns,s,n_uniq)
    integer, intent(in)  :: ns
    character(len=*)     :: s(ns)
    integer, intent(out) :: n_uniq
    !------------------------------
    character(len=fstr)  :: stmp(ns)
    integer :: uniq(ns)
    integer :: i,j,ii

    if(ns < 1)then 
      write(*,*)'error in get_unique_strings:'
      stop 'number of strings to check == 0'
    endif

    uniq = 1 
    do i=1,ns    
      stmp(i) = s(i)
    end do
    
    do i=2,ns    
      do j=1,i-1
        if( s(j) == s(i)) uniq(i) = 0 
      end do
    end do

    !reinitialize s
    s(1:ns) = ""

    !count unique
    n_uniq = sum(uniq)

    !transfer uniq 
    ii = 0
    do i=1,ns
      if(uniq(i) == 1)then
        ii = ii + 1
        s(ii) = stmp(i)
      endif
    end do

  end subroutine get_unique_strings

  !-----------------------------------------------
  ! return a random number, precision "sp" using  
  ! fortran90 intrinsic "random_number"
  !-----------------------------------------------
  subroutine ran1(np,ran)  
    implicit none 
    integer,intent(in)   :: np
    real(sp),intent(out) :: ran(np)
    real(sp)             :: x
    integer              :: m
    integer,allocatable  :: k(:) 

    call random_seed   !! zhuxm
    call random_seed(size=m)
    allocate(k(m))
    call random_seed(get=k(1:m))
    call random_seed(put=k(1:m))
    call random_number(ran)
    deallocate(k)

  end subroutine ran1

  !-----------------------------------------------
  ! return a random number in specified interval 
  !-----------------------------------------------
  subroutine ran_from_range(fmin,fmax,np,rany)  
    implicit none 
    integer,intent(in)         :: np
    real(sp),intent(in)        :: fmin,fmax
    real(sp),intent(out)       :: rany(np)
    real(sp)                   :: ran(np)
    integer                    :: n

    call ran1(np,ran)
    do n=1,np
       rany(n)=(fmax - fmin) * ran(n) + fmin
    enddo
  end subroutine ran_from_range

  !-----------------------------------------------
  ! return unit random number (-1 or 1) 
  !-----------------------------------------------
  subroutine unitrand(np,unitr)  
    implicit none
    integer,intent(in)   :: np
    real(sp),intent(out) :: unitr(np)
    real(sp)             :: tmp,ran(np)
    integer              :: n
    call ran1(np,ran)
    do n=1,np
       tmp   = ran(n)
       unitr(n) = sign(1.0_sp,tmp-0.5_sp)
    enddo
  end subroutine unitrand 

  !-----------------------------------------------
  ! return random number from normal distribution 
  ! with mean ->  mean and standard dev -> sigma
  ! from?
  !-----------------------------------------------
  real(sp) function normal(mean,sigma) 
    implicit none 
    real(sp), intent(in) ::  mean
    real(sp), intent(in) ::  sigma 
    integer  :: flag 
    real(sp) :: fac,gsave,rsq,r1,r2,tmp,randy(1)
    save flag,gsave 
    data flag /0/ 

    
    if (flag == 0) then 
      rsq=2.0_sp 
      do while(rsq >= 1.0_sp.or.rsq == 0.0_sp) 
        call ran1(1,randy)
        r1=2.0_sp*randy(1)-1.0_sp
        call ran1(1,randy)
        r2=2.0_sp*randy(1)-1.0_sp 
        rsq=r1*r1+r2*r2 
      enddo 
      fac=sqrt(-2.0_sp*log(rsq)/rsq) 
      gsave=r1*fac 
      tmp=r2*fac 
      flag=1 
    else 
      tmp=gsave 
      flag=0 
    endif 
    
    normal=tmp*sigma+mean !tmp is distribution of standard normal
    return 
  end function normal


  !-----------------------------------------------
  ! return random number from normal distribution 
  ! with mean = 0.0 and dev = 1. 
  !-----------------------------------------------
  function normal_zeromean() result(normal)
    implicit none 
    integer  :: flag 
    real(sp) :: fac,gsave,rsq,r1,r2,randy(1)
    real(sp) :: normal
    save flag,gsave 
    data flag /0/ 
    if (flag == 0) then 
      rsq=2.0_sp 
      do while(rsq >= 1.0_sp.or.rsq == 0.0_sp)
        call ran1(1,randy)
        r1=2.0_sp*randy(1)-1.0_sp
        call ran1(1,randy)
        r2=2.0_sp*randy(1)-1.0_sp 
        rsq=r1*r1+r2*r2 
      enddo 
      fac=sqrt(-2.0_sp*log(rsq)/rsq) 
      gsave=r1*fac 
      normal=r2*fac 
      flag=1 
    else 
      normal=gsave 
      flag=0 
    endif 
    return 
  end function normal_zeromean

subroutine spline(n,x,y,yp1,ypn,ysp)

  !-----------------------------------------------
  ! fit a cubic splint (zero tension) to data 
  ! from numerical recipes
  ! in:  
  !   n:  dimension of data
  !   x:  independent variable 
  !   y:  dependent variable
  ! yp1:  boundary condition at i=1
  ! ypn:  boundary condition at i=n
  !  ys:  spline
  !-----------------------------------------------
    integer,  intent(in ) :: n
    real(sp), intent(in ) :: x(n),y(n)
    real(sp), intent(in ) :: yp1,ypn
    real(sp), intent(out) :: ysp(n)  
    !------------------------------
    integer               :: i,k
    real(sp)              :: p,qn,sig,un
	real(sp),allocatable  :: u(:)
  
    allocate(u(n))

	if (yp1 .gt. 0.99e30) then
        ysp(1)=0.
        u(1)=0.
    else
        ysp(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*ysp(i-1)+2.
       ysp(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    if (ypn .gt. 0.99e30) then
       qn=0.
       un=0.
    else
       qn=0.5
       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    ysp(n)=(un-qn*u(n-1))/(qn*ysp(n-1)+1.)
    do k=n-1,1,-1
       ysp(k)=ysp(k)*ysp(k+1)+u(k)
    end do

	deallocate(u)
    return
end subroutine spline
  
subroutine splint(nn,xx,yy,x1,y1,x2,y2)
  !-----------------------------------------------
  ! cubic spline interpolation
  !
  ! added by zhuxm 2011.7.25 from Numerical Recipes Fortran
  !
  ! in:  
  !   n:  dimension of data
  !  xa:  independent variable 
  !  ya:  dependent variable
  ! y2a:  the output from SPLINE subroutine
  !   x:  the given position to interpolate
  ! out:
  !   y:  the interpolated value returned
  !-----------------------------------------------
  implicit none

  integer, intent(in )  :: nn
  real(sp),intent(in )  :: xx(nn),yy(nn),x1
  real(sp),intent(out)  :: y1
  real(sp),intent(in ),optional  :: x2
  real(sp),intent(out),optional  :: y2

  integer               :: k,khi,klo,m,i,n,xn
  integer,parameter     :: nmax=100
  real(sp)              :: A,B,C,D,deltx,x,y
  real(sp),allocatable,dimension(:)  :: tmpx,tmpy,xa,ya,y2a

!!-------calculate nmax------------
  n=nn
  m=0

  do 
    if(n<nmax) then
	  n=n+n-1
	  m=m+1
	endif
	if(n>=nmax) exit
  enddo

  allocate(xa(n),ya(n),y2a(n),tmpx(n),tmpy(n))

  xa(1:nn)=xx(1:nn)
  ya(1:nn)=yy(1:nn)

!!-------make sure xa is increased from i=1,...,nn ------------
  if(xa(1)>xa(nn)) then
	 do k=1,nn
	    tmpx(k) =xa(nn-k+1)
		tmpy(k) =ya(nn-k+1)
	 enddo
	 xa(1:nn) =tmpx(1:nn)
	 ya(1:nn) =tmpy(1:nn)
  endif

  xa(1) =-1.0_sp
  xa(nn)= 0.0_sp

!!----------refine the grid-----------------------------------
  n=nn
  do k=1,m     
	 do i=1,n
	    tmpx(2*(i-1)+1)=xa(i)
		if(i/=n) tmpx(2*i)=(xa(i)+xa(i+1))/2.0_sp

	    tmpy(2*(i-1)+1)=ya(i)
		if(i/=n) tmpy(2*i)=(ya(i)+ya(i+1))/2.0_sp
	 enddo
	 n=2*n-1
	 xa(1:n)=tmpx(1:n)
	 ya(1:n)=tmpy(1:n)
  enddo
  call spline(n,xa,ya,0.0_sp,0.0_sp,y2a)

  xn=1
  if(present(x2)) xn=2

  do i=1,xn
     if(i==1) then
	    x=x1
	 elseif(present(x2) .and. i==2) then
	    x=x2
	 endif

!!----------find x location in array XA-----------------------
     klo=0  ;  khi=0
     if(x<xa(1))then
       write(*,*)'Error in SPLINT, x=', x,' < xa(1)=',xa(1)
	   stop
     elseif(x>xa(n))then
       write(*,*)'Error in SPLINT, x=', x,' > xa(n)=',xa(n)
  	   stop
     else
       do k=1,n-1
          if(xa(k)<=x .and. x<=xa(k+1))then
	         klo=k
	  	     khi=k+1
		     exit
	      endif
       enddo
     endif
!!----------Get the second derivatives of the interpolating function
!!----------at the tabulated points XA                -------------

print*,'khi klo',khi,klo

     deltx=xa(khi)-xa(klo)

     if(deltx == 0.0_sp) then
        write(*,*) k,x,xa
        write(*,*) deltx,xa(khi),khi,xa(klo),klo
        stop 'Bad xa input in splint, which must be increased!'
     endif
!!----------Get the cubic-spline interpolated value
     A=(xa(khi)-x)/deltx
     B=(x-xa(klo))/deltx
     C=(A**3-A)*(deltx**2)/6.0_sp
     D=(B**3-B)*(deltx**2)/6.0_sp

     y=A*ya(klo)+B*ya(khi)+C*y2a(klo)+D*y2a(khi)

     if(y<0.0_sp)then
!       write(*,*)y,' Error in SPLINT, kh < 0.0 happened!!'
	   y=1.0E-7_sp
     endif
     if(i==1) then
	    y1=y
	 elseif(present(x2) .and. i==2) then
	    y2=y
	 endif
  enddo
  deallocate(tmpx,tmpy)
  deallocate(xa,ya,y2a)

  return
end subroutine splint
     
  !-----------------------------------------------
  ! return x,y pseudo-normally distributed in 
  ! circle centered at (xc,yc), of radius r
  ! r will extend to two standard devs
  !-----------------------------------------------
  subroutine normal_circle(n,xc,yc,rad,x,y) 
    implicit none 
    integer,  intent(in)  :: n
    real(sp), intent(in)  :: xc,yc,rad 
    real(sp), intent(out) :: x(n),y(n)
    real(sp)              :: theta(n),rval
    integer               :: i 

    call ran_from_range(0.0_sp,2*pi,n,theta)
    do i=1,n
      rval  = normal(0.0_sp,rad/2.)
      x(i) = rval*cos(theta(i))
      y(i) = rval*sin(theta(i))
    end do
      
  end subroutine normal_circle 

  !-----------------------------------------------
  ! random distribution of [n] particles in a square 
  ! [xmin,xmax,ymin,ymax]
  !-----------------------------------------------
  subroutine random_square(n,xmin,xmax,ymin,ymax,x,y)
    implicit none
    integer,  intent(in)  :: n
    real(sp), intent(in)  :: xmin
    real(sp), intent(in)  :: xmax
    real(sp), intent(in)  :: ymin
    real(sp), intent(in)  :: ymax
    real(sp), intent(out) :: x(n)
    real(sp), intent(out) :: y(n)
    real(sp) :: theta,rval
    integer  :: i

    call ran_from_range(xmin,xmax,n,x)
    call ran_from_range(ymin,ymax,n,y)
  end subroutine random_square 

!------------------------------------------------------------------
! check times 
!------------------------------------------------------------------
 Subroutine check_times(mbeg,mend,smin,smax,fbeg,fend)

  implicit none
  real(sp), intent(in) :: mbeg,mend 
  real(sp), intent(in) :: smin,smax
  real(sp), intent(in) :: fbeg,fend
  !-----------------------------
  real(sp)             :: tmin,tmax

  call drawline('-')
  write(*,*)'checking coherency between: '
  write(*,*)'   spawning times '
  write(*,*)'   simulation begin/end time'
  write(*,*)'   model forcing begin/end time (if applicable)'
  call drawline('.')
  write(*,*)'time of earliest spawning: ',gettime(int(smin))
  write(*,*)'time of latest   spawning: ',gettime(int(smax))
  if(fbeg >= 0.0 .and. fend >= 0.0)then
    write(*,*)'netcdf forcing start time: ',gettime(int(fbeg))
    write(*,*)'netcdf forcing end time  : ',gettime(int(fend))
  endif
  write(*,*)'model simulation beg time: ',gettime(int(mbeg))
  write(*,*)'model simulation end time: ',gettime(int(mend))

  tmin = min(mend,mbeg)
  tmax = max(mend,mbeg)

  !make sure spawning range is contained inside the simulation range
  if(smin < tmin .or. smin > tmax .or. smax < tmin .or. smax > tmax)     &
     stop 'fatal error: spawning not contained within simulation bounds' 
  !make sure simulation range is contained within forcing data time range 
  !(if forcing active)
  if(fbeg >= 0.0 .and. fend >= 0.0)then
    if(tmin < fbeg .or. tmin > fend .or. tmax < fbeg .or. tmax > fend)   &
       stop 'fatal error: simulation not contained within forcing bounds' 
  endif
 End Subroutine check_times

!--------------------------------------------------------------------
! 'gcdist' function calculates great circle distance between
!   two points on the earth.
!   If lon/lat of two ponts are very close to each other, then
!   use flat-surface formula
! This function returns a real(sp) value in meters
! Input argument:
!  rlon1:  real longitude of point-1 in degrees
!  rlon1:  real latitude of point-1 in degrees
!  rlon2:  real longitude of point-2 in degrees
!  rlat2:  real latitude of point-2 in degrees
!  zfeng 10/20/2014
!--------------------------------------------------------------------
 real(sp) function gcdist(rlon1,rlat1,rlon2,rlat2)

 use gparms
 implicit none

 real(sp)  :: deltalon, deltalat, meanlat
 real(sp)  :: rlon1, rlat1, rlon2, rlat2
 real(sp)  :: clat1, clat2, slat1, slat2, cdlon, crd
 real(sp), parameter :: mindeg = 0.03_sp
 ! minimum degree difference considered as in a flat surface

 if (abs(rlat1-rlat2)<mindeg .AND. abs(rlon1-rlon2)<mindeg) then
    deltalon = (rlon1-rlon2)*d2r
    deltalat = (rlat1-rlat2)*d2r
    meanlat  = (rlat1+rlat2)*0.5_sp*d2r
    gcdist = earth*sqrt(deltalat**2+(deltalon*cos(meanlat))**2)
 else
    clat1 = cos(rlat1*d2r)
    slat1 = sin(rlat1*d2r)
    clat2 = cos(rlat2*d2r)
    slat2 = sin(rlat2*d2r)
    cdlon = cos((rlon1-rlon2)*d2r)
    crd   = slat1*slat2+clat1*clat2*cdlon
    gcdist = earth*acos(crd)

 endif

 end function gcdist

 !-------------------------------------------------------------------
 ! 'nearest_ij' finds indice (i_nrst,j_nrst) of minimum distance
 ! 'dist' is a 2-D array of size [nx,ny]
 ! zfeng 10/27/2014

 subroutine nearest_ij(dist,i_nrst,j_nrst)

   use gparms
   implicit none

   real(sp), intent(in)   :: dist(nx,ny)
   integer,  intent(out)  :: i_nrst,j_nrst
   integer                :: i,j
   real(sp)               :: min_dist
   real(sp), parameter    :: mdval = 0.1 ! 0.1 m

   ! initialize
   i_nrst = -999
   j_nrst = -999


   min_dist = minval(dist)

 !  write(*,*) 'min_dist ='
 !  write(*,*) min_dist

   do j = 1, ny
     do i = 1, nx
        if( abs(dist(i,j)-min_dist) < mdval ) then
           i_nrst = i
           j_nrst = j
           exit
        endif
     enddo
   enddo

! write error message of either i_nrst or j_nrst NOT found
   if(i_nrst == -999) then
     write(*,*) "Error: index 'i_nrst' for particle NOT found!"
   elseif (j_nrst == -999) then
     write(*,*) "Error: index 'j_nrst' for particle NOT found!"
   endif

 end subroutine nearest_ij


!-------------------------------------------------------------------

!====================================================================



end module utilities
