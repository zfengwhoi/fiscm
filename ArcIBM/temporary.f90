module temporary
    implicit none

subroutine trilinear_interp(np,x,y,z,is_ind,js_ind,kz_ind,istatus,var_char,iframe,flag,interpdata)

  implicit none
  integer, intent(in)       :: np,iframe
  real(sp),intent(in)       :: x(np),y(np),z(np) ! particle locations
  integer, intent(in)       :: is_ind(np),js_ind(np),kz_ind(np),istatus(np) ! scalar indices
  character(len=*)          :: var_char
  integer, intent(in)       :: flag ! Variable flag if vector or scalar quantity to be interpolated
                                   ! scalar: = 0; vector: = 1, ww = 2.
  real(sp),intent(inout)    :: interpdata(np)    ! data after bilinear interpolation

  !----------------------------
  integer                   :: iv, jv  ! indices of containing grid cell
  integer                   :: i,j,n,olvl
  real(sp), pointer         :: field3d(:,:,:)
  real(sp)      :: field2d(nx,ny)
  real(sp)      :: coef1, coef2  ! coefs for linear interpolation in the vertical (z-dir)
 ! variables for the particle-containing cell
  real(sp)      :: xcell(2,2),ycell(2,2),fcell(2,2)
  real(sp)      :: dx1, dx2, dpx1,dpx2,dy
  real(sp)      :: ew_val1, ew_val2 ! intermediate field values after interpolation in East-West (x-dir)
  real(sp)      :: ytemp1,ytemp2    ! intermediate y-positions

  ! obtain forcing data with pointer 'field3d'
  call get_forcing_f3(trim(var_char),field3d,iframe)
!  write(*,*) size(field3d(1,:,1))
!  write(*,*) 'field3d(1,:,1) =', field3d(1,:,1)

  field2d = 0.0_sp
  interpdata = 0.0_sp

 PARTICLE: do n = 1, np ! particle loop

 !  ocean levels of the grid where particle is
  olvl = mask(is_ind(n),js_ind(n))
!  write(*,*) 'olvl=',olvl

!  write(*,*) 'istatus(',n,')=',istatus(n)
 ACTIVE: if(istatus(n) < 1 ) then ! only do interpolation for active partciles
       cycle
  else

  select case(flag)
  case(0)  ! interp scalar quantity at ocean level center (temperature, diatoms, flagellates)

      ! First, linear intepolation in the vertical (Z-layer)
      if(z(n) <= zt(1)) then ! particle is located above top layer center
         field2d = field3d(:,:,1) ! use data of top layer
      elseif (z(n) >= zt(olvl)) then ! particle is located below ocean bottom of containing grid
         field2d = field3d(:,:,olvl) ! use data of bottom layer
      else  ! conduct linear interpolation
         coef1 = (zt(kz_ind(n)) - z(n))/(zt(kz_ind(n)) - zt(kz_ind(n)-1))
         coef2 = (z(n) - zt(kz_ind(n)-1))/(zt(kz_ind(n)) - zt(kz_ind(n)-1))
         field2d = coef1 * field3d(:,:,kz_ind(n)-1) + coef2 * field3d(:,:,kz_ind(n))
      endif

      ! Second, bilinear interpolation in the horizontal
      ! Before interp, find nearest [iv,jv] of 4 surounding vector grid points
      call search4corners(x(n),y(n),is_ind(n),js_ind(n),iv,jv)

      if(jv==ny) jv = ny-1 ! particle at inner boundary, move it to ny-1

      ! obtain positions and data at four surrouding scalar points
      if(iv<nx) then ! normal case, no boundary points involved
       xcell(1,1) = slon(iv,jv);     xcell(1,2) = slon(iv,jv+1)
       xcell(2,1) = slon(iv+1,jv);   xcell(2,2) = slon(iv+1,jv+1)
       ycell(1,1) = slat(iv,jv);     ycell(1,2) = slat(iv,jv+1)
       ycell(2,1) = slat(iv+1,jv);   ycell(2,2) = slat(iv+1,jv+1)
       fcell(1,1) = field2d(iv,jv);  fcell(1,2) = field2d(iv,jv+1)
       fcell(2,1) = field2d(iv+1,jv);fcell(2,2) = field2d(iv+1,jv+1)
      elseif(iv==nx)  then! particle at lateral periodic boundary
       xcell(1,1) = slon(nx,jv);     xcell(1,2) = slon(nx,jv+1)
       xcell(2,1) = slon(1,jv);      xcell(2,2) = slon(1,jv+1)
       ycell(1,1) = slat(nx,jv);     ycell(1,2) = slat(nx,jv+1)
       ycell(2,1) = slat(1,jv);      ycell(2,2) = slat(1,jv+1)
       fcell(1,1) = field2d(nx,jv);  fcell(1,2) = field2d(nx,jv+1)
       fcell(2,1) = field2d(1,jv);   fcell(2,2) = field2d(1,jv+1)
      else
       write(*,*) "Error in 'trilinear_interp': Particle #", n
      endif

      ! Now, do interpolation

      ! calculate differences
      dx1 = xcell(2,1)-xcell(1,1) ! difference in two x-dir corners
      if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
      if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

      dx2 = xcell(2,2)-xcell(1,2)
      if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
      if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

      if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases, use average of 4 corners
         interpdata(n)= (fcell(1,1)+fcell(2,1)+fcell(1,2)+fcell(2,2))*0.25_sp
         exit
      endif

      dpx1 = xcell(2,1) - x(n)  ! difference between particle x-location to a corner
      if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
      if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

      dpx2 = x(n)- xcell(1,1)
      if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
      if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

      ! calculate intermediate value-1 by interp in x-dir
      ew_val1 = fcell(1,1)*dpx1/dx1+fcell(2,1)*dpx2/dx1

      ! calculate intermediate longitudes (ytemp1)
      ytemp1  = ycell(1,1)*dpx1/dx1+ycell(2,1)*dpx2/dx1

      dpx1 = xcell(2,2) - x(n)  ! difference between particle x-location to a corner
      if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
      if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

      dpx2 = x(n)- xcell(1,2)
      if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
      if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

      ! calculate intermediate value-2 by interp in x-dir
      ew_val2 = fcell(1,2)*dpx1/dx2 + fcell(2,2)*dpx2/dx2

      ! calculate intermediate longitudes (ytemp2)
      ytemp2  = ycell(1,2)*dpx1/dx2 + ycell(2,2)*dpx2/dx2

      if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= (fcell(1,1)+fcell(2,1)+fcell(1,2)+fcell(2,2))*0.25_sp

      else
      ! interpolate in y-dir
      interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2
      endif

  case(1) ! interp vector quantity at ocean level center (u & v)

      ! First, linear intepolation in the vertical (Z-layer)
      if(z(n) <= zt(1)) then ! particle is located above top layer center
         field2d = field3d(:,:,1) ! use data of top layer
      elseif (z(n) >= zt(olvl)) then ! particle is located below ocean bottom of containing grid
         field2d = field3d(:,:,olvl) ! use data of bottom layer
      else  ! conduct linear interpolation
         coef1 = (zt(kz_ind(n)) - z(n))/(zt(kz_ind(n)) - zt(kz_ind(n)-1))
         coef2 = (z(n) - zt(kz_ind(n)-1))/(zt(kz_ind(n)) - zt(kz_ind(n)-1))
         field2d = coef1 * field3d(:,:,kz_ind(n)-1) + coef2 * field3d(:,:,kz_ind(n))
      endif

      ! Second, bilinear interpolation in the horizontal
      iv = is_ind(n); jv = js_ind(n)

      if(jv==1) jv = 2 ! particle at inner boundary, move it to 2

      ! obtain positions and data at four surrouding vector points
      if(iv>1) then ! normal case, no boundary points involved
       xcell(1,1) = vlon(iv-1,jv-1); xcell(1,2) = vlon(iv-1,jv)
       xcell(2,1) = vlon(iv,jv-1);   xcell(2,2) = vlon(iv,jv)
       ycell(1,1) = vlat(iv-1,jv-1); ycell(1,2) = vlat(iv-1,jv)
       ycell(2,1) = vlat(iv,jv-1);   ycell(2,2) = vlat(iv,jv)
       fcell(1,1) = field2d(iv-1,jv-1);fcell(1,2) = field2d(iv-1,jv)
       fcell(2,1) = field2d(iv,jv-1);fcell(2,2) = field2d(iv,jv)
      elseif(iv==1) then ! particle at lateral periodic boundary
       xcell(1,1) = vlon(nx,jv-1);  xcell(1,2) = vlon(nx,jv)
       xcell(2,1) = vlon(1,jv-1);   xcell(2,2) = vlon(1,jv)
       ycell(1,1) = vlat(nx,jv-1);  ycell(1,2) = vlat(nx,jv)
       ycell(2,1) = vlat(1,jv-1);   ycell(2,2) = vlat(1,jv)
       fcell(1,1) = field2d(nx,jv-1);fcell(1,2) = field2d(nx,jv)
       fcell(2,1) = field2d(1,jv-1);fcell(2,2) = field2d(1,jv)
      else
       write(*,*) "Error in 'trilinear_interp': Particle #", n
      endif

      ! Now, do interpolation

      ! calculate differences
      dx1 = xcell(2,1)-xcell(1,1) ! difference in two x-dir corners
      if(dx1 < -180.0_sp)  dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
      if(dx1 >  180.0_sp)  dx1 = dx1 - 360.0_sp

      dx2 = xcell(2,2)-xcell(1,2)
      if(dx2 < -180.0_sp)  dx2 = dx2 + 360.0_sp
      if(dx2 >  180.0_sp)  dx2 = dx2 - 360.0_sp

      if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases
         interpdata(n)= (fcell(1,1)+fcell(2,1)+fcell(1,2)+fcell(2,2))*0.25_sp
         exit
      endif

      dpx1 = xcell(2,1) - x(n)  ! difference between particle x-location to a corner
      if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
      if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

      dpx2 = x(n)- xcell(1,1)
      if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
      if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

      ! calculate intermediate value-1 by interp in x-dir
      ew_val1 = fcell(1,1)*dpx1/dx1+fcell(2,1)*dpx2/dx1

      ! calculate intermediate longitudes (ytemp1)
      ytemp1  = ycell(1,1)*dpx1/dx1+ycell(2,1)*dpx2/dx1

      dpx1 = xcell(2,2) - x(n)  ! difference between particle x-location to a corner
      if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
      if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

      dpx2 = x(n)- xcell(1,2)
      if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
      if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

      ! calculate intermediate value-2 by interp in x-dir
      ew_val2 = fcell(1,2)*dpx1/dx2+fcell(2,2)*dpx2/dx2

      ! calculate intermediate longitudes (ytemp2)
      ytemp2  = ycell(1,2)*dpx1/dx2+ycell(2,2)*dpx2/dx2

      if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= (fcell(1,1)+fcell(2,1)+fcell(1,2)+fcell(2,2))*0.25_sp

      else
      ! interpolate in y-dir
      interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
                + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2
     endif

  case(2)  ! interp scalar quantity at ocean level bottom (ww)

      ! First, linear intepolation in the vertical (Z-layer)
      if(z(n) <= zw(1)) then ! particle is located above top layer center
         field2d = field3d(:,:,1) ! use data of top layer
      elseif (z(n) >= zw(olvl)) then ! particle is located below ocean bottom of containing grid
         field2d = field3d(:,:,olvl) ! use data of bottom layer
      else  ! conduct linear interpolation
         coef1 = (zw(kz_ind(n)) - z(n))/(zw(kz_ind(n)) - zw(kz_ind(n)-1))
         coef2 = (z(n) - zw(kz_ind(n)-1))/(zw(kz_ind(n)) - zw(kz_ind(n)-1))
         field2d = coef1 * field3d(:,:,kz_ind(n)-1) + coef2 * field3d(:,:,kz_ind(n))
      endif

      ! Second, bilinear interpolation in the horizontal
      ! Before interp, find nearest [iv,jv] of 4 surounding vector grid points
      call search4corners(x(n),y(n),is_ind(n),js_ind(n),iv,jv)

      if(jv==ny) jv = ny-1 ! particle at inner boundary, move it to ny-1

      ! obtain positions and data at four surrouding scalar points
      if(iv<nx) then ! normal case, no boundary points involved
       xcell(1,1) = slon(iv,jv);     xcell(1,2) = slon(iv,jv+1)
       xcell(2,1) = slon(iv+1,jv);   xcell(2,2) = slon(iv+1,jv+1)
       ycell(1,1) = slat(iv,jv);     ycell(1,2) = slat(iv,jv+1)
       ycell(2,1) = slat(iv+1,jv);   ycell(2,2) = slat(iv+1,jv+1)
       fcell(1,1) = field2d(iv,jv);  fcell(1,2) = field2d(iv,jv+1)
       fcell(2,1) = field2d(iv+1,jv);fcell(2,2) = field2d(iv+1,jv+1)
      elseif(iv==nx) then ! particle at lateral periodic boundary
       xcell(1,1) = slon(nx,jv);     xcell(1,2) = slon(nx,jv+1)
       xcell(2,1) = slon(1,jv);      xcell(2,2) = slon(1,jv+1)
       ycell(1,1) = slat(nx,jv);     ycell(1,2) = slat(nx,jv+1)
       ycell(2,1) = slat(1,jv);      ycell(2,2) = slat(1,jv+1)
       fcell(1,1) = field2d(nx,jv);  fcell(1,2) = field2d(nx,jv+1)
       fcell(2,1) = field2d(1,jv);   fcell(2,2) = field2d(1,jv+1)
      else
       write(*,*) "Error in 'trilinear_interp': Particle #", n
      endif

      ! Now, do interpolation

      ! calculate differences
      dx1 = xcell(2,1)-xcell(1,1) ! difference in two x-dir corners
      if(dx1 < -180.0_sp) dx1 = dx1 + 360.0_sp ! cell crossing 0 meridian
      if(dx1 >  180.0_sp) dx1 = dx1 - 360.0_sp

      dx2 = xcell(2,2)-xcell(1,2)
      if(dx2 < -180.0_sp) dx2 = dx2 + 360.0_sp
      if(dx2 >  180.0_sp) dx2 = dx2 - 360.0_sp

      if(abs(dx1)<0.0001_sp .OR. abs(dx2) <0.0001_sp) then ! special cases
         interpdata(n)= (fcell(1,1)+fcell(2,1)+fcell(1,2)+fcell(2,2))*0.25_sp
         exit
      endif

      dpx1 = xcell(2,1) - x(n)  ! difference between particle x-location to a corner
      if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
      if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

      dpx2 = x(n)- xcell(1,1)
      if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
      if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

      ! calculate intermediate value-1 by interp in x-dir
      ew_val1 = fcell(1,1)*dpx1/dx1+fcell(2,1)*dpx2/dx1

      ! calculate intermediate longitudes (ytemp1)
      ytemp1  = ycell(1,1)*dpx1/dx1+ycell(2,1)*dpx2/dx1

      dpx1 = xcell(2,2) - x(n)  ! difference between particle x-location to a corner
      if(dpx1 < -180.0_sp) dpx1 = dpx1 + 360.0_sp
      if(dpx1 >  180.0_sp) dpx1 = dpx1 - 360.0_sp

      dpx2 = x(n)- xcell(1,2)
      if(dpx2 < -180.0_sp) dpx2 = dpx2 + 360.0_sp
      if(dpx2 >  180.0_sp) dpx2 = dpx2 - 360.0_sp

      ! calculate intermediate value-2 by interp in x-dir
      ew_val2 = fcell(1,2)*dpx1/dx2+fcell(2,2)*dpx2/dx2

      ! calculate intermediate longitudes (ytemp2)
      ytemp2  = ycell(1,2)*dpx1/dx2+ycell(2,2)*dpx2/dx2

      if(abs(ytemp2-ytemp1)<0.0001_sp) then
        interpdata(n)= (fcell(1,1)+fcell(2,1)+fcell(1,2)+fcell(2,2))*0.25_sp
      else
      ! interpolate in y-dir
      interpdata(n) = (ytemp2-y(n))/(ytemp2-ytemp1)*ew_val1 &
           + (y(n)-ytemp1)/(ytemp2-ytemp1)*ew_val2
      endif

  case default
    write(*,*) "Error in 'trilinear_interp': "
    write(*,*) "Given flag is wrong! Scalar: flag=0; vector: flag=1; ww: flag=2"
  End select

  endif ACTIVE

 enddo PARTICLE ! particle loop

end subroutine trilinear_interp

end module temporary
