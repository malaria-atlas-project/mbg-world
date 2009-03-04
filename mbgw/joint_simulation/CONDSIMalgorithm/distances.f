! Copyright (c) Anand Patil, 2007

      SUBROUTINE euclidean(D,x,y,nx,ny,ndx,ndy,symm)

cf2py double precision dimension(nx,ny), intent(out)::D
cf2py double precision dimension(nx,ndx), intent(in)::x
cf2py double precision dimension(ny,ndy), intent(in)::y 
cf2py logical intent(in), optional:: symm=0
cf2py integer intent(hide), depend(x)::nx=shape(x,0)
cf2py integer intent(hide), depend(y)::ny=shape(y,0)
cf2py integer intent(hide), depend(x)::ndx=shape(x,1)
cf2py integer intent(hide), depend(y)::ndy=shape(y,1)

      DOUBLE PRECISION D(nx,ny), x(nx,ndx), y(ny,ndy)
      integer nx,ny,ndx,ndy,i,j,k
      LOGICAL symm
      DOUBLE PRECISION dist, dev
      
      do i=1,nx

        if(symm) then            
          do j=i+1,ny
            dist = 0.0D0
            do k=1,ndx
              dev=(x(i,k) - y(j,k))
              dist = dist + dev*dev
            end do
            D(i,j) = dsqrt(dist)
            D(j,i) = D(i,j)
          end do

        else
          do j=1,ny
            dist = 0.0D0
            do k=1,ndx
              dev=(x(i,k) - y(j,k))
              dist = dist + dev*dev
            end do
            D(i,j) = dsqrt(dist)
          end do    
        end if
          
      end do
      RETURN
      END



      SUBROUTINE geographic(D,x,y,nx,ny,symm)
! First coordinate is longitude, second is latitude.
! Assumes r=1.

cf2py double precision intent(out), dimension(nx,ny) :: D
cf2py double precision intent(in), dimension(nx,2) :: x
cf2py double precision intent(in), dimension(ny,2) :: y
cf2py logical intent(in), optional :: symm = 0
cf2py integer intent(hide), depend(x)::nx=shape(x,0)
cf2py integer intent(hide), depend(y)::ny=shape(y,0)

      DOUBLE PRECISION D(nx,ny), x(nx,2), y(ny,2)
      integer nx,ny,i,j,j_lo
      LOGICAL symm
      DOUBLE PRECISION clat1, clat2, dlat, dlon, a, sterm, cterm
      
      do i=1,nx
        clat1 = dcos(x(i,2))
        if(symm) then
            D(i,i)=0.0D0            
            j_lo = i+1
        else 
            j_lo = 1
        end if
        
        do j=j_lo,ny
            clat2 = dcos(y(j,2))
            dlat = (x(i,2)-y(j,2))*0.5D0
            dlon = (x(i,1)-y(j,1))*0.5D0
            a=dsin(dlat)**2 + clat1*clat2*dsin(dlon)**2
            sterm = dsqrt(a)
            cterm = dsqrt(1.0D0-a)
            D(i,j) = 2.0D0*DATAN2(sterm,cterm)    
            if(symm) then                  
                D(j,i) = D(i,j)
            end if
        end do          
      end do
      RETURN
      END

      FUNCTION ag(A,B)

! Finds angle of inclination of vector pointing from A to B on sphere
! First coordinate is longitude, second is latitude.
! Assumes r=1.
! From http://www.movable-type.co.uk/scripts/latlong.html (Thankyouthankyouthankyou)

cf2py intent(out) angle

      DOUBLE PRECISION A(2), B(2), sterm, cterm
      DOUBLE PRECISION dlon, lat2, lat1, ag
      PARAMETER (pi=3.141592653589793238462643d0)
      if (D.EQ.0.0D0) then
          ag=-999999D0
          RETURN
      end if
      
      dlon = B(1)-A(1)
      lat1 = A(2)
      lat2 = B(2)
      
      sterm = dsin(dlon)*dcos(lat2)
      cterm = dcos(lat1)*dsin(lat2)-dsin(lat1)*dcos(lat2)*dcos(dlon)
      ag = datan2(sterm, cterm)
      
      RETURN
      END

      SUBROUTINE aniso_geo_rad(D,x,y,nx,ny,inc,ecc,symm)
! First coordinate is longitude, second is latitude.
! Assumes r=1.

cf2py intent(out) D
cf2py logical intent(optional) :: symm = 0
cf2py intent(hide) nx
cf2py intent(hide) ny

      DOUBLE PRECISION D(nx,ny), x(nx,2), y(ny,2)
      integer nx,ny,i,j,j_lo
      LOGICAL symm
      DOUBLE PRECISION clat1, clat2, dlat, dlon, a, sterm, cterm
      DOUBLE PRECISION slat1, slat2, C(2), B(2), inc, ecc,q,dtheta

      do i=1,nx
        clat1 = dcos(x(i,2))
        slat1 = dsin(x(i,2))
        
        if(symm) then
            D(i,i)=0.0D0            
            j_lo = i+1
        else 
            j_lo = 1
        end if
        do j=j_lo,ny
            clat2 = dcos(y(j,2))
            slat2 = dsin(y(j,2))
            dlat = (x(i,2)-y(j,2))
            dlon = (x(i,1)-y(j,1))
            a=dsin(dlat*0.5D0)**2 + clat1*clat2*dsin(dlon*0.5D0)**2
            sterm = dsqrt(a)
            cterm = dsqrt(1.0D0-a)
            D(i,j) = 2.0D0*DATAN2(sterm,cterm)


            if (D(i,j).GT.0.0D0) then
                
                dtheta = theta-inc                
                dtheta = dcos(dtheta)                          
                dtheta=ecc*ecc*dtheta*dtheta
                D(i,j) = D(i,j) * dsqrt(1.0D0 - dtheta)


            end if
            
            if(symm) then                  
                D(j,i) = D(i,j)
            end if
        end do          
      end do
      RETURN
	END
