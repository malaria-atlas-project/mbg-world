      subroutine logit(theta,n,ltheta)
c Maps (0,1) -> R.
cf2py intent(hide) n
cf2py intent(out) ltheta
      DOUBLE PRECISION theta(n), ltheta(n)
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)
      INTEGER n, i
      do i=1,n
          if (theta(i).LE.0) then
              ltheta(i) = -inf
          else if (theta(i).GE.1) then
              ltheta(i) = inf
          else
              ltheta(i) = dlog(theta(i) / (1.0D0-theta(i)))
          endif
      end do
      RETURN
      END
c


      subroutine invlogit(ltheta,n,theta) c Maps R -> (0,1).
cf2py intent(hide) n
cf2py intent(out) theta
      DOUBLE PRECISION theta(n), ltheta(n)
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)
      INTEGER n, i
      do i=1,n
          theta(i) = 1.0D0 / (1.0D0 + dexp(-ltheta(i)))
      end do
      RETURN
      END

c

      subroutine stukel_logit(theta,n,ltheta,a1,a2,na1,na2)
!
! Reference: Therese A. Stukel, 'Generalized Logistic Models', ! JASA vol 83 no 402, pp.426-431 (June 1988) !
cf2py intent(hide) n, na1, na2
cf2py intent(out) ltheta
cf2py intent(copy) theta
      DOUBLE PRECISION theta(n), ltheta(n)
      DOUBLE PRECISION a1(na1), a2(na2), a1t, a2t
      LOGICAL a1_isscalar, a2_isscalar
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)
      INTEGER n, i, na1, na2

      a1t = a1(1)
      a2t = a2(1)

      CALL logit(theta,n,ltheta)

      a1_isscalar = (na1.LT.n)
      a2_isscalar = (na2.LT.n)

      do i=1,n

          if (ltheta(i).GT.0.0D0) then
              if (.NOT.a1_isscalar) then
                  a1t = a1(i)
              end if
              if (a1t.GT.0.0D0) then
                  ltheta(i)=dlog(ltheta(i)*a1t+1.0D0)/a1t
              else if (a1t.LT.0.0D0) then
                  ltheta(i) = (1.0D0-dexp(-ltheta(i)*a1t))/a1t
              end if

          else if (ltheta(i).LT.0.0D0) then
              if (.NOT.a2_isscalar) then
                  a2t = a2(i)
              end if
              if (a2t.GT.0.0D0) then
                  ltheta(i)=-dlog(-ltheta(i)*a2t+1.0D0)/a2t
              else if (a2t.LT.0.0D0) then
                  ltheta(i)=-(1.0D0-dexp(ltheta(i)*a2t))/a2t
              end if

          else
              ltheta(i) = 0.0D0
          end if

      end do

      RETURN
      END

c

      subroutine stukel_invlogit(ltheta,n,theta,a1,a2,na1,na2)
!
! Reference: Therese A. Stukel, 'Generalized Logistic Models', ! JASA vol 83 no 402, pp.426-431 (June 1988) !
cf2py intent(hide) n, na1, na2
cf2py intent(out) theta
cf2py intent(copy) ltheta
      DOUBLE PRECISION theta(n), ltheta(n)
      DOUBLE PRECISION a1(na1), a2(na2), a1t, a2t
      LOGICAL a1_isscalar, a2_isscalar
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)
      INTEGER n, i, na1, na2

      a1t = a1(1)
      a2t = a2(1)

      a1_isscalar = (na1.LT.n)
      a2_isscalar = (na2.LT.n)

      do i=1,n
          if (ltheta(i).GT.0.0D0) then
              if (.NOT.a1_isscalar) then
                  a1t = a1(i)
              end if
              if (a1t.GT.0.0D0) then
                  ltheta(i) = (dexp(a1t*ltheta(i))-1.0D0)/a1t
              else if (a1t.LT.0.0D0) then
                  ltheta(i) = -dlog(1.0D0-a1t*ltheta(i))/a1t
              end if

          else if (ltheta(i).LT.0.0D0) then
              if (.NOT.a2_isscalar) then
                  a2t = a2(i)
              end if
              if (a2t.GT.0.0D0) then
                  ltheta(i) = -(dexp(-a2t*ltheta(i))-1.0D0)/a2t
              else if (a2t.LT.0.0D0) then
                  ltheta(i) = dlog(1.0D0+a2t*ltheta(i))/a2t
              end if

          else
              ltheta(i) = 0.5D0
          end if

      end do


      CALL invlogit(ltheta,n,theta)

      RETURN
      END


