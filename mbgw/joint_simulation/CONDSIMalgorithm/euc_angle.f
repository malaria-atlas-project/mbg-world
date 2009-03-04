      SUBROUTINE euc_angle(theta,x,y,nx,ny,symm)
! First coordinate is longitude, second is latitude.
! Assumes r=1.


      DOUBLE PRECISION theta(nx,ny), x(nx,2), y(ny,2),dlat,dlon
      integer nx,ny,i,j,j_lo
      LOGICAL symm

      do i=1,nx

        if(symm) then
            theta(i,i)=0.0D0
            j_lo = i+1
        else
            j_lo = 1
        end if
        do j=j_lo,ny
            dlat = (x(i,2)-y(j,2))
            dlon = (x(i,1)-y(j,1))
            theta(i,j) = DATAN2(dlat,dlon)
        end do
      end do
      RETURN
	END
