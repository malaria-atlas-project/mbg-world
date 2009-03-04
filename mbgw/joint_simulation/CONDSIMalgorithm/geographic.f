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