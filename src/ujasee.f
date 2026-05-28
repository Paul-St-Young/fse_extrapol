      subroutine ujasee(rk,uk,zeta,ifs0,vol,rho,hbs2m,e2)
      implicit none
      integer ifs0
      real*8 rk,uk,zeta,kf(2)
      real*8 s0,vk,arg,q,tn,s,vol,rho,hbs2m,e2,sk0
      include 'geometry.cm'

      kf(1)=2.d0*pi*( (1.d0+zeta)* (rho*ndim)
     .                 /(4.d0*pi*(ndim-1)))**(1.d0/ndim)
      kf(2)=2.d0*pi*( (1.d0-zeta) * (rho*ndim)
     .             /(4.d0*pi*(ndim-1)))**(1.d0/ndim)

      if(ifs0.eq.0) then
         s0=1.d0
      else
         s0=0.5d0*(1.d0+zeta)*sk0(rk,kf(1))
     .    +0.5d0*(1.d0-zeta)*sk0(rk,kf(2))
      end if
      vk=(2.d0*pi*(ndim-1)*e2)/(rk**(ndim-1))
      arg=2.d0*vk/(hbs2m*rk**2)*rho
      arg=1.d0+arg*s0**2
      q=sqrt(arg)
      tn=1.d0
      s=q*(1+q*tn)/(q+tn)
      s=(-1.d0+s)/s0
      uk=s/(2.d0*rho*vol)
      end
