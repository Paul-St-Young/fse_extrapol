      subroutine dposl(a,lda,n,b)
c     LINPACK Cholesky solve A*x = b using factors from dpoco/dpofa.
c     b is overwritten with the solution.
      integer lda,n
      double precision a(lda,*),b(*)
      double precision ddot,t
      integer k,kb

      do 10 k = 1, n
         t = ddot(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k) - t)/a(k,k)
   10 continue

      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call daxpy(k-1,t,a(1,k),1,b(1),1)
   20 continue
      return
      end
