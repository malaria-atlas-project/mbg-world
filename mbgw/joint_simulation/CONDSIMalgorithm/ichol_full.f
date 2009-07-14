!Argument definitions:

!c : Input covariance matrix, will not be destroyed
!n : size of c
!sig : Output matrix. Upper triangle will be overwritten with Cholesky
!factor. Initialize to zero.
!m : Output integer
!p : Output integer array of length n



      subroutine ichol_full(c,n,sig,m,p)
c
c Incomplete cholesky factorization
c Author: Anand Patil
c Date: May 6, 2007
c Port of mex function chol_incomplete.c by Matthias Seeger
c http://www.kyb.tuebingen.mpg.de/bs/people/seeger/
c
cf2py double precision dimension(n,n), intent(in)::c
cf2py double precision dimension(n,n), intent(out)::sig
cf2py integer dimension(n), intent(out)::p
cf2py double precision dimension(n), intent(hide)::rowvec
cf2py double precision dimension(n), intent(hide)::diag
cf2py integer intent(hide), depend(c):: n = shape(c,0)
cf2py integer intent(out)::m
cf2py double precision intent(in) :: reltol
cf2py threadsafe

      DOUBLE PRECISION c(n,n), sig(n,n), diag(n)
      DOUBLE PRECISION rowvec(n)
      integer p(n), n, m, i, j
      DOUBLE PRECISION maxdiag, tol, dtemp

      EXTERNAL DGEMV
      DOUBLE PRECISION ZERO, ONE, RELTOL, NEGONE
      PARAMETER (zero=0.0D0)
      PARAMETER (one=1.0D0)
      PARAMETER (negone = -1.0D0)
      PARAMETER (reltol = 0.000001

* DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
      EXTERNAL DSWAP
* DSWAP(N,DX,INCX,DY,INCY)

      EXTERNAL IDAMAX
* IDAMAX(N,DX,INCX)

!       Make diagonal and index vectors
      do i=1,n
        diag(i) = c(i,i)
        p(i)=i
      enddo

      maxdiag = diag(idamax(n,diag,1))

      tol = maxdiag * reltol
      m = n
!       Main loop
      do i=1,n

!         Find maximum remaining pivot
        l = idamax(n-i+1,diag(i),1)+i-1
        maxdiag = diag(l)


!         Early return if there are no big pivots left
        if (maxdiag .LE. tol) then
          do j=1,n
            p(j) = p(j)-1
          enddo
          m = i-1
          return
        endif

        if (i .NE. l) then
!         Swap p and diag's elements i and l

          itemp = p(i)
          p(i) = p(l)
          p(l) = itemp

          dtemp = diag(i)
          diag(i) = diag(l)
          diag(l) = dtemp

!         Swap the i and lth columns of sig
          CALL DSWAP(i,sig(1,i),1,sig(1,l),1)
        endif

!       Write diagonal element
        sig(i,i) = dsqrt(diag(i))

!       Assemble the row vector
        if (i.LT.n) then
            do j=i+1,n
              rowvec(j) = c(p(i),p(j))
            enddo
        endif

          if (i.GT.1) then

!               BLAS-less DGEMV might be useful if you ever do the
sparse version.
!               do j=i+1,n
!                 do k=1,i-1
!                   rowvec(j)=rowvec(j)-sig(k,j)*sig(k,i)
!                 enddo
!               enddo

!         Implement Cholesky algorithm.
            CALL DGEMV('T',i-1,n-i,negone,sig(1,i+1),
     1                  n,
     2                  sig(1,i),
     3                  1,
     4                  one,rowvec(i+1),1)
          endif

        if (i.LT.n) then
          do j=i+1,n
            sig(i,j) = rowvec(j) / sig(i,i)
            diag(j) = diag(j) - sig(i,j)*sig(i,j)
          enddo
        endif

      enddo

      do i=1,n
        p(i)=p(i)-1
      enddo

      return
      end
