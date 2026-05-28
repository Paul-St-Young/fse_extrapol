      subroutine basis(m,maxm,s)
      implicit none
      integer m,maxm
      integer ialpha,ien,i,j,ibeta,k,info
      integer fact(0:2*m+1)
      integer ipiv(m+1)
      double precision a(m+1,m+1),b(m+1,m+1)
      double precision x(m+1,m+1)
      integer delta(0:m,0:m)
      double precision s(0:maxm,0:2*maxm+1)
      integer nd

      nd=m+1

      fact(0)=1
      do i=1,2*m+1
         fact(i)=fact(i-1)*i
      enddo

      do ien=1,nd
         do ialpha=1,nd
            a(ien,ialpha)=dble(fact(ien+m))/dble(fact(ien+m-ialpha+1))
         enddo
      enddo

      do ialpha=1,nd
         do ibeta=1,nd
            if (ialpha.gt.ibeta) then
               b(ibeta,ialpha)=0.d0
            else
               b(ibeta,ialpha)=-1.d0/dble(fact(ibeta-ialpha))
            endif
         enddo
      enddo

c     Invert A in place via Gauss-Jordan elimination.
      call invmat(a,nd,info)
      if (info.ne.0) then
         write(6,*) 'basis: singular system, info=',info
         stop
      end if

c     X = B * A^{-1}
      do i=1,nd
         do j=1,nd
            x(i,j)=0.d0
            do k=1,nd
               x(i,j)=x(i,j)+b(i,k)*a(k,j)
            enddo
         enddo
      enddo

      do i=0,m
         do j=0,m
            delta(i,j)=0
         enddo
         delta(i,i)=1
      enddo

      do i=0,m
         do j=0,m
            s(i,j)=dble(delta(i,j))/dble(fact(i))
            s(i,j+m+1)=x(i+1,j+1)
         enddo
      enddo

      return
      end

      subroutine invmat(a,n,info)
c     Replace a(n,n) with its inverse via Gauss-Jordan with partial
c     pivoting.  Sized for tiny matrices (m+1 <= 7 in production).
      implicit none
      integer n,info,i,j,k,ip
      double precision a(n,n),tmp,pivot
      double precision aug(n,2*n)

      info=0
      do i=1,n
         do j=1,n
            aug(i,j)=a(i,j)
            aug(i,n+j)=0.d0
         enddo
         aug(i,n+i)=1.d0
      enddo

      do k=1,n
c        partial pivoting on column k
         ip=k
         pivot=abs(aug(k,k))
         do i=k+1,n
            if (abs(aug(i,k)).gt.pivot) then
               pivot=abs(aug(i,k))
               ip=i
            end if
         enddo
         if (pivot.eq.0.d0) then
            info=k
            return
         end if
         if (ip.ne.k) then
            do j=1,2*n
               tmp=aug(k,j)
               aug(k,j)=aug(ip,j)
               aug(ip,j)=tmp
            enddo
         end if
c        normalize pivot row
         tmp=aug(k,k)
         do j=1,2*n
            aug(k,j)=aug(k,j)/tmp
         enddo
c        eliminate other rows
         do i=1,n
            if (i.ne.k) then
               tmp=aug(i,k)
               if (tmp.ne.0.d0) then
                  do j=1,2*n
                     aug(i,j)=aug(i,j)-tmp*aug(k,j)
                  enddo
               end if
            end if
         enddo
      enddo

      do i=1,n
         do j=1,n
            a(i,j)=aug(i,n+j)
         enddo
      enddo
      return
      end
