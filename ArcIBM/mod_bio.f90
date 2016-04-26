MODULE mod_bio

!====================================================================
!-----Calanus model, version 1.0------------------
!------Jess McNally, WHOI--------------------------
!  2011/7/28       Xueming Zhu  modified ini_bio
!  2014/12/15      Zhixuan Feng modified advance_bio to add constraints to 
!                               stage-based copepod development rates 
!====================================================================
     
	Use gparms
	Use mod_igroup

    IMPLICIT NONE 
  
! REAL BEL_BETA ! Beta for Belehradek equation, common for all calanus

!------STATE VARIABLES declared in copepod.nml-----

! ---STAGE 
! Stage number keeps track of stage, gender and maturity. 

! INTEGER iStage(N_IND)
! Keeps track of stage as an integer

! REAL pasd(N_IND) 
!'Proportional Accumulative Stage Days', to keep track of progress in development

! Stages = 1 Egg
!        = 2 N1
!        = 3 N2
!        = 4 N3
!        = 5 N4
!        = 6 N5
!        = 7 N6
!        = 8 C1
!        = 9 C2
!        =10 C3
!        =11 C4
!        =12 C5
!        =13 Male
!        =14 Immature Female
!        =15 Mature Female 


!parameters
  integer, parameter :: nstages  = 15
  real(sp)           :: BEL_BETA = -2.05 
!=====Passed from species specific vars module====
! For Belehradek equation Duration = a*(temperature+alpha)^beta
! will just declare them here now, for testing. These are values for C. fin
 
  contains

!----------------------------------------------------
!  initinalize status and biology spawning time
!----------------------------------------------------
  subroutine init_bio(g,Nind_start)
       use mod_output
       use utilities 

       implicit none
       type(igroup), intent(inout) :: g
       integer,intent(in)          :: Nind_start
       integer                     :: ios, i, np,npoint
       real(sp) , pointer          :: x(:),y(:),z(:),s(:),tspawn(:)
       integer  , pointer          :: istatus(:)
       real(sp)                    :: tmp
       
       !set current problem size
       if(Nind_start < 0)then
          write(*,*)'must specify more than 0 individuals in init_bio'
          stop
       elseif(Nind_start > g%Tnind)then
          write(*,*)'number of individuals in init_bio (',Nind_start,') exceeds problem size',g%Tnind
          stop
       endif
       g%Nind = Nind_start
       !set problem size
       np = g%Nind
       
       !read the namelist (if any) 
       if(g%biology)then
          rewind(iunit)
          read(unit=iunit,nml=nml_copepod,iostat=ios)
          write(*, nml=nml_copepod)
          close(iunit)
          !------------------------------------
          ! set the spawning time
          !------------------------------------
          call get_state('tspawn',g,tspawn)
          !        tspawn(1:np) = 0.0_sp*day_2_sec
          nullify(tspawn)
          allocate(zptini(np))
       endif
       !-----------------------------------
       ! 2D,3D -> initialize x,y 
       !-----------------------------------
       if(g%space_dim > 1)then
          call get_state('x',g,x)
          call get_state('y',g,y)
          call get_state('status',g,istatus)
          if(g%space_dim == 3)then
             call get_state('s',g,s)
             call get_state('z',g,z)
          endif
          !    call random_square(np,-2500._sp,2500._sp,-2500._sp,2500._sp,x,y)
          !    Revised by Xinyou Lin
          inquire(file=trim(g%statefile),exist=fexist)
          if(.not. fexist)then
             write(*,*)trim(g%statefile),' does not exist, stopping...'
             stop
          endif
          write(*,*)'Initializing x,y,z/s,status and spawn time,'
          write(*,*)' using file: ',trim(g%statefile)
          open(unit=iunit,file=trim(g%statefile),form='formatted')
          read(iunit,*)npoint
          if(npoint /= np)then
             write(*,*)'The number of points are wrong in',trim(g%statefile)
             stop
          endif
          
          if(g%space_dim == 2)then
             if(g%biology)then
                do i=1,np
                   read(iunit,*) npoint,x(i),y(i),zptini(i)
                enddo
                zptini = zptini*3600.0_sp  ! hours to seconds
             else
                do i=1,np
                   read(iunit,*) npoint,x(i),y(i)
                enddo
             endif
          elseif(g%space_dim == 3)then
             allocate(zpini(np))
             if(g%biology)then
                do i=1,np
                   read(iunit,*) npoint,x(i),y(i),zpini(i),zptini(i)
                enddo
                ! hours to seconds
                zptini=zptini*3600.0_sp
             else
                do i=1,np
                   read(iunit,*) npoint,x(i),y(i),zpini(i)
                enddo
             endif
             
             if(sz_cor == 0)then
                s(1:np) = zpini(1:np)
                z       = 0.0
             elseif(sz_cor == 1)then
                z(1:np) = zpini(1:np)
                s       = 0.0
             endif
             nullify(s,z)  
          endif
          
          !---- transform x from 0~360 to -180~180------- 
          if(spherical == 1) then
             if(north_pole) x = x + 180.0_sp  !! zhuxm angle_pre 
             where(x >= 360.0_sp) x=x-360.0_sp
          endif
          !---- transform x from 0~360 to -180~180------- 
!!!!!!!!!!cnt>0
          !       if(x(i) >= 0.0_sp .and. x(i) <=180.0_sp)then
          !          x(i) = x(i) + 180.0_sp
          !       elseif( x(i) > 180.0_sp .and. x(i) <=360.0_sp)  then
          !          x(i) = x(i) - 180.0_sp
          !       endif
!!!!!!!!!!!
          close(iunit)
          istatus=1;
          
          nullify(x,y,istatus)
       endif
       !-----------------------------------
       ! 3D -> initialize s-coordinate
       !-----------------------------------
       !add parameters to netcdf file header 
       call add_cdfglobs(g,"info","some kind of info")
     end subroutine init_bio

!------------------------------------------------------------------
! advance the biology (this routine is called from the main loop at
!                      a time interval of DT_bio )
!------------------------------------------------------------------
subroutine advance_bio(g,mtime)

	implicit none

	type(igroup), intent(inout) :: g
    real(sp),     intent(in)    :: mtime

    real(sp),     pointer :: PASD(:),T(:), Phyto(:)
    integer ,     pointer :: stage(:),istatus(:),diapause(:)
    integer               :: i,N
    real(sp)              :: dt, stg_duri, stg_durj, stg_durdiff
    real(sp)              :: develop_rate ! [stage/day] zfeng add 
    real(sp)              :: food_multiplier, food_saturation, food_starvation
    real(sp), parameter   :: max_rate = 1.0_sp ! define a maximum development rate to contrain calculated rate 
    ! This number is given by best-case estimation but will be adjusted to a literature value later  

    ! construct pointers to access and modify state variables for the group
    call get_state('PASD',g,PASD)
    call get_state('T',g,T)
    call get_state('stage',g,stage)
    call get_state('diapause',g,diapause) 
    call get_state('status',g,istatus)
!    if(bio_fd ==1)call get_state('CHL_2DIM',g,FOOD)
    
    ! Use 'bio_food' to activate food-dependency
    ! Use 'food_source' to allow select ;
    if(bio_food .AND. food_source == 1)   call get_state('phytoplankton',g,Phyto) ! zfeng 03/02/2015
    if(bio_food .AND. food_source == 2)   call get_state('food',g,Phyto)          ! Food (P+M) zfeng 11/13/2015 
    if(bio_food .AND. food_source == 3)   call get_state('SPM',g,Phyto)           ! Subsurface Phytoplankton Maxima
    if(bio_food .AND. food_source == 4)   call get_state('SFM',g,Phyto)           ! Subsurface food (P+M) Maxima    

  !set problem size
    N  = g%nind
    dt = g%DT_bio*sec_2_day !time step in days 

!------------MAIN LOOP-------------------
    DO i=1, N
!----check for diapause---needs to be filled in
       if(istatus(i) /= ACTIVE) cycle
!increment development
       SELECT CASE (stage(i))
!Development from egg to C5
       CASE(1:12)
!          stg_duri = BEL_A(stage(i))*((T(i)+BEL_ALPHA)**BEL_BETA) 
!          stg_durj = BEL_A(stage(i)+1)*((T(i)+BEL_ALPHA)**BEL_BETA)
!          stg_durdiff = stg_durj - stg_duri
           
           ! Temperature-dependent development days 
           stg_durdiff =  BEL_A(stage(i))*((T(i)+BEL_ALPHA)**BEL_BETA)

           ! Development rate (zfeng added constraint conditions to avoid extremes)
           if (T(i) <= -1.0_sp * BEL_ALPHA) then
               develop_rate = 0.0_sp  ! no growth due to very low temperature
           else 
               develop_rate = min(1.0_sp/stg_durdiff,max_rate) ! constrain development rate no larger than 1.0 stage/day
           endif
          
           !
!           if(bio_fd == 0) then  ! only temperature-dependent case
!              PASD(i) = PASD(i) + dt * develop_rate   
!           elseif(bio_fd == 1) then  ! temperature- & food-dependent case
!              PASD(i) = PASD(i) + dt*develop_rate*(1.0_sp - exp(-FOOD(i)/PC0))  ! old line        
            
           if(.NOT. bio_food) then ! only temperature-dependent, assuming unlimited food
              PASD(i) = PASD(i) + dt * develop_rate 
           elseif(bio_food) then   ! food-limitation factor (zfeng added conditions to avoid extremes)            
           ! Before Stage N3, food is not required for development
           ! So only applying food limitation after N3
             if(PASD(i)>=4.0_sp) then
                food_multiplier = 1.0_sp - exp(-1.0_sp * Phyto(i) / PC0)
                food_saturation = -1.0_sp * log(0.01_sp) * PC0  ! food saturated condition
               ! food_starvation = -1.0_sp * log(0.99_sp) * PC0  ! food starved condition
                food_starvation = LowFood  !       
                if( Phyto(i) >= food_saturation) then
                  develop_rate = develop_rate   ! food-independent
                elseif( Phyto(i) <= food_starvation) then
                  develop_rate = 0.0_sp         ! food is negligible, no development
                else
                  develop_rate = develop_rate * food_multiplier  ! food-dependent development   
                endif              
             endif
   
           ! Update stage
           PASD(i) = PASD(i) + dt * develop_rate 
       
           endif ! food-dependency

          
!If stage has surpassed 13, need to make half of those go to 14 instead (half male, half female)
!Simple solution is that 'odd individuals' stay male, and 'even' individuals become female
!Can maybe put a random generator in here later?
           if(stage(i)>13 .AND. MOD(i,2)==0) THEN
              stage(i) = stage(i)+1
           endif
!Males, eventually can put some function in to model their death rate.						
       CASE(13)
           CYCLE
!Females, immature and mature
       CASE(14:)
           CYCLE
       END SELECT
!truncate PASD to nearest whole number, to keep track of stage as int
       stage(i) = AINT(PASD(i))
    ENDDO
! nullify pointers 
    nullify(PASD,T,istatus,stage,diapause)
    if(bio_food) nullify(Phyto)
!   if(bio_fd ==1)nullify(FOOD)
End Subroutine advance_bio

!---------------------------------------------------------
! activate particles after spawning time reached 
!   for each particle:
!     if the status is still initial (status = UNKNOWN) and
!     we have exceeded the spawning time (or < for backwards)
!     set status to active 
!---------------------------------------------------------
subroutine activate(ng,g,t,direction)

  implicit none
  integer     , intent(in)    :: ng
  type(igroup), intent(inout), dimension(ng) :: g
  real(sp), intent(in) :: t
  integer,  intent(in) :: direction
  !---------------------------------
  real(sp), pointer    :: tspawn(:)
  integer , pointer    :: istatus(:)
  integer              :: n,p,np
  
  !forward model
  if(direction > 0)then
    do n=1,ng
	if(g(n)%biology) then
         call get_state('tspawn',g(n),tspawn)
         call get_state('status',g(n),istatus)
         np = g(n)%nind

          do p=1,np
 !          if(structured .AND. jp(p)==1)  istatus(p) = EXITED
           if(istatus(p)==UNKNOWN .and. t>=tspawn(p)) istatus(p) = ACTIVE
           if(t<=zptini(p)) istatus(p) = UNKNOWN
          end do
	   
        else
          call get_state('status',g(n),istatus)
          np = g(n)%nind
          do p=1,np
 !           if(structured .AND. jp(p)==1) istatus(p) = EXITED
            if(istatus(p)==UNKNOWN) istatus(p) = ACTIVE
          end do
	endif
    end do
  !reverse model
  else
    do n=1,ng
	if(g(n)%biology) then
          call get_state('tspawn',g(n),tspawn)
          call get_state('status',g(n),istatus)
          np = g(n)%nind
          do p=1,np
            if(istatus(p)==UNKNOWN .and. t>=tspawn(p)) istatus(p) = ACTIVE ! bug fix on t<=tspawn(p)
            if(t>=zptini(p)) istatus(p) = UNKNOWN
          end do
	else
          call get_state('status',g(n),istatus)
          np = g(n)%nind
          do p=1,np
             if(istatus(p)==UNKNOWN) istatus(p) = ACTIVE
          end do
	endif
    end do
  endif
end subroutine activate

End Module  mod_bio 
