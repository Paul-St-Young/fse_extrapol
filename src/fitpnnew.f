       subroutine fitpnnew(v,y,t,rk,wt,nk,nf,m,mp,rad,nknots,t0,vt0,t1,
     a  vt1,coul,vmad,vol,s,rknot,delta,rkbig,vkbig,vklrbig,nbig)
       implicit none
       integer m,nalphaknot,nknots,k,nk,i,alpha,ialpha,maxn,n,j,istart
       integer iskip,nalc,jrow,icol,nf,ld,info,ialc,mp
       double precision delta,vt0,vt1,vmad,vmad1,vmad2,chisq,vol,r,rad
       double precision s(0:mp,0:2*mp+1),ddplus(0:2*m+1)
     .                 ,ddminus(0:2*m+1)
       double precision c(0:nknots*(m+1)-1,0:nk)
       double precision v(0:nk),wt(0:nk),y(0:nk)
       double precision a(0:nknots*(m+1)-1,0:nknots*(m+1)-1)
       double precision w(0:nk),t(0:nknots*(m+1)-1),b(0:nknots*(m+1)-1)
       double precision bb(nknots*(m+1)),x(nknots*(m+1))
       double precision aa(nknots*(m+1),nknots*(m+1)),wknot(0:nknots)
       double precision rk(0:nk)
       logical coul,t0,t1
       double precision rknot(0:nknots),rcond,z(nknots*(m+1))
       integer nbig
       real*8 rkbig(nbig),vkbig(nbig),vklrbig(nbig)
       include 'geometry.cm'

c      3D only -- guaranteed by parameter (ndim=3) in geometry.cm.

       ld=nknots*(m+1)
       wknot(0)=0.d0
       do i=1,nknots
         wknot(i)=1.d0
       enddo

       delta=rad/dble(nknots)

       if(rk(nk)*delta.le.10.d0*2.d0*pi) then
            write(*,*)'ATTENTION, too much not points'
            write(*,*)'could lead to trouble'
       end if

c      Hermite-style polynomial basis coefficients on [0,1].
       call basis(m,mp,s)

       nalphaknot=nknots*(m+1)-1
       maxn=2*m+1

c      Fourier transform of all polynomials at every shell vector.
       do k=0,nk
         ialpha=-1
         do i=0,nknots-1
           r=i*delta
           call splint3D(maxn,r,delta,rk(k),vol,coul,ddplus,ddminus)
           do alpha=0,m
             ialpha=ialpha+1
             c(ialpha,k)=0.d0
             do n=0,maxn
               c(ialpha,k)=c(ialpha,k)+s(alpha,n)
     &         *(ddplus(n)+wknot(i)*ddminus(n)*(-1)**(n+alpha))
             enddo
             c(ialpha,k)=delta**alpha*c(ialpha,k)
           enddo
         enddo
       enddo

c      Build normal-equation matrices for the constrained least-squares
c      problem in (nf+1..nk) (high-k tail).
       do j=0,nalphaknot
         b(j)=0.d0
         do k=nf+1,nk
           b(j)=b(j)+c(j,k)*v(k)*wt(k)
         enddo
         do i=0,nalphaknot
           a(j,i)=0.d0
           do k=nf+1,nk
             a(j,i)=a(j,i)+c(j,k)*c(i,k)*wt(k)
           enddo
         enddo
       enddo

c      Apply constraints t(0)=vt0 (value) and/or t(1)=vt1 (derivative).
       if (t0.and..not.t1) then
         istart=1
         iskip=0
         nalc=nalphaknot
         do j=0,nalphaknot
           b(j)=b(j)-vt0*a(j,0)
         enddo
         t(0)=vt0
       elseif (t1.and..not.t0) then
         istart=0
         iskip=1
         nalc=nalphaknot
         do j=0,nalphaknot
           b(j)=b(j)-vt1*a(j,1)
         enddo
         t(1)=vt1
       elseif (t0.and.t1) then
         istart=2
         iskip=0
         nalc=nalphaknot-1
         do j=0,nalphaknot
           b(j)=b(j)-vt0*a(j,0)
           b(j)=b(j)-vt1*a(j,1)
         enddo
         t(0)=vt0
         t(1)=vt1
       else
         istart=0
         iskip=-1
         nalc=nalphaknot+1
       endif

c      Load reduced matrices.
       jrow=0
       do j=istart,nalphaknot
         if (j.ne.iskip) then
           jrow=jrow+1
           bb(jrow)=b(j)
           icol=0
           do i=istart,nalphaknot
             if (i.ne.iskip) then
               icol=icol+1
               aa(jrow,icol)=a(j,i)
             endif
           enddo
         endif
       enddo

c      Solve with LINPACK Cholesky.
       call dpoco(aa,ld,nalc,rcond,z,info)
       call dposl(aa,ld,nalc,bb)
       do jrow=1,nalc
         x(jrow)=bb(jrow)
       enddo

c      Scatter solution into t(0:nalphaknot).
       jrow=0
       do j=istart,nalphaknot
         if (j.ne.iskip) then
           jrow=jrow+1
           t(j)=x(jrow)
         endif
       enddo

c      Knot positions.
       rknot(0)=0.d0
       do i=1,nknots
         rknot(i)=rknot(i-1)+delta
       enddo

c      Long-range Fourier component on the shell grid.
       do k=0,nk
         w(k)=0.d0
         do i=0,nalphaknot
           w(k)=w(k)+t(i)*c(i,k)
         enddo
         y(k)=v(k)-w(k)
       enddo

c      Chi^2 diagnostic.
       chisq=0.d0
       do k=nf+1,nk
         chisq=chisq+y(k)**2*wt(k)
       enddo
       write(*,*) ' fitpn: Chi^2 =',chisq
       write(*,*)

c      Madelung-like constant.
       vmad1=t(1)
       vmad2=0.d0
       do k=0,nf
         vmad2=vmad2+wt(k)*y(k)
       enddo
       if (coul) then
         vmad=vmad1
       else
         vmad=vmad2+t(0)
       endif

c      Long-range continuation on the big-k grid for the diagnostic
c      integration in the calling routine.
       do k=1,nbig
         ialpha=-1
         do i=0,nknots-1
           r=i*delta
           call splint3D(maxn,r,delta,rkbig(k),vol,coul,ddplus,ddminus)
           do alpha=0,m
             ialpha=ialpha+1
             c(ialpha,1)=0.d0
             do n=0,maxn
               c(ialpha,1)=c(ialpha,1)+s(alpha,n)
     &         *(ddplus(n)+wknot(i)*ddminus(n)*(-1)**(n+alpha))
             enddo
             c(ialpha,1)=delta**alpha*c(ialpha,1)
           enddo
         enddo
         w(1)=0.d0
         do i=0,nalphaknot
           w(1)=w(1)+t(i)*c(i,1)
         enddo
         vklrbig(k)=(vkbig(k)-w(1))*vol
       enddo

       return
       end
