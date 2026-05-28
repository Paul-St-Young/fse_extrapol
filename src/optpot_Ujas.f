      subroutine optpot_Ujas(ifs0,nkact,cutr,iwrite)
      implicit none
      include 'geometry.cm'
      include 'system.cm'
      include 'kshell.cm'
      integer nkact
      real*8 wtk(0:Mnshlls),vk(0:Mnshlls),ulr(0:Mnshlls)
      integer k,iwrite,is
      real*8 cutr
      integer maxm,maxknots
      parameter(maxm=6,maxknots=80)
      integer nknots,m,maxn,nk1,ifs0
      logical t0,t1,coul
      real*8 sm(0:maxm,0:2*maxm+1)
     .,t(0:maxknots*(maxm+1)-1),rknot(0:maxknots)
     .,vt0,vt1,vmad,delta
      real*8 vlrinf,dv,vlr,x,y

      call fillk(rknorm,wtk,kmult,nshlls,mnshex,vol)
      nk1=nkact

      vk(0)=0.d0
      ulr(0)=0.d0
      do k=1,mnshex
        call ujasee(rknorm(k),vk(k),zeta,ifs0,vol,rho,hbs2m,e2)
        ulr(k)=vk(k)
      end do
      do k=1,iksbig
        call ujasee(rkbig(k),vkbig(k),zeta,ifs0,vol,rho,hbs2m,e2)
        vklrbig(k)=vkbig(k)
      end do

      m=2          ! cubic spline
      maxn=2*m+1
      nknots=16    ! number of knots

      t0=.false.
      vt0=0.d0
      t1=.true.    ! cusp constraint U'(0) = -e2/(2(d-1) hbs2m)
      vt1=-e2/((ndim-1.d0)*2.d0*hbs2m)
      coul=.false. ! spline is NOT divided by 1/r

      call fitpnnew(vk,ulr,t,rknorm,wtk,mnshex,nk1,m,maxm,cutr,nknots
     . ,t0,vt0,t1,vt1,coul,vmad,vol,sm,rknot,delta
     . ,rkbig,vkbig,vklrbig,iksbig)

      vlr=0.d0
      do is=1,nk1
        x=rknorm(is)
        call splint(rksk(0),sk(0),sk2(0),nsk+1,x,y)
        vlr=vlr+2.d0*(kmult(is)-kmult(is-1))
     .    *hbs2m*x**2*rho*(2.d0*vk(is)-ulr(is))*ulr(is)*vol*y
      end do
      vlrinf=0.d0
      do is=1,iksbig
        x=rkbig(is)
        call splint(rksk(0),sk(0),sk2(0),nsk+1,x,y)
        vlrinf=vlrinf+2.d0*x**(ndim+1)*hbs2m*rho
     .    *(2.d0*vkbig(is)*vol-vklrbig(is))*vklrbig(is)*y
      end do
      vlrinf=(ndim-1.d0)*pi/(2.d0*pi)**ndim*vlrinf*del
      dv=(vlrinf-vlr)
      write(6,*)'***************************************'
      write(6,*)'  Kinetic energy correction Ujas  '
      write(6,*)'     using optimized potential with ik=',iksbig
      if(electrongas) then
        write(6,*)' dT/N=',dv,' Ry'
      else
        write(6,*)' dT/N=',dv,' Ha'
      end if
      write(6,*)'***************************************'
      if(iwrite.gt.0) then
        do k=1,nk1
          write(2,*) rknorm(k),ulr(k),ulr(k)/vk(k)
        end do
      end if

      end
