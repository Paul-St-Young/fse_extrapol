      subroutine ewald(kc,alpha)
      implicit none
      real*8 kc,alpha,vlrn
      real*8 vewald_lr_k,x,y
      real*8 vlrinf,dv
      integer is
      include 'geometry.cm'
      include 'system.cm'
      include 'kshell.cm'

c     Discrete sum over k-shells inside (kmax, rksk(nsk)).
      vlrn=0.d0
      do is=1,Nshlls
        x=rknorm(is)
        if(x.ge.kmax) goto 88
        if(x.ge.rksk(nsk)) goto 88
        call splint(rksk(0),sk(0),sk2(0),nsk+1,x,y)
        vlrn=vlrn+e2*vewald_lr_k(x,kc,alpha)
     .         *y*(kmult(is)-kmult(is-1))/vol
      end do
88    continue

c     Continuous-k integral on the dense rkbig grid.
      vlrinf=0.d0
      do is=1,iksbig
        x=rkbig(is)
        call splint(rksk(0),sk(0),sk2(0),nsk+1,x,y)
        vlrinf=vlrinf+x**(ndim-1)*vewald_lr_k(x,kc,alpha)*y
      end do
      vlrinf=e2*(ndim-1.d0)*pi/(2.d0*pi)**ndim*vlrinf*del

      dv=(vlrinf-vlrn)
      write(6,*)'***************************************'
      write(6,*)'  Potential energy correction  '
      write(6,*)'     using Ewald potential with ik=',iksbig
      if(electrongas) then
        write(6,*)' dV/N=',dv,' Ry'
      else
        write(6,*)' dV/N=',dv,' Ha'
      end if
      write(6,*)'***************************************'
      end
