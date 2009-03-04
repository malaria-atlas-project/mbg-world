      SUBROUTINE dist_to_aniso(out,D,theta,nx,ny,inc,ecc,symm)
! First coordinate is longitude, second is latitude.
! Assumes r=1.

      DOUBLE PRECISION D(nx,ny),theta(nx,ny),out(nx,ny)
      DOUBLE PRECISION ecc, inc, eccsq, dtheta
      integer nx,ny,i,j,j_lo
      LOGICAL symm

      eccsq = ecc*ecc

      do i=1,nx
        if(symm) then
            j_lo = i+1
        else
            j_lo = 1
        end if

        do j=j_lo,ny

            if (D(i,j).GT.0.0D0) then

                dtheta = theta(i,j)-inc
                dtheta = dcos(dtheta)
                dtheta= eccsq*dtheta*dtheta
                out(i,j) = D(i,j) * dsqrt(1.0D0 - dtheta)

            end if

            if(symm) then
                out(j,i) = D(i,j)
            end if
        end do
      end do
      RETURN
	END