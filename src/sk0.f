      real*8 function sk0(k,kf)
      implicit none
      include 'geometry.cm'
      real*8 k,kf,x
      sk0=0.d0
      if(kf.le.0.d0) return

      x=0.5d0*k/kf
      if(x.lt.1.d0) then
        if(ndim.eq.3) then
          sk0=0.5d0*x*(3.d0-x*x)
        else if(ndim.eq.2) then
          sk0=(2.d0/pi)*(asin(x)+x*sqrt(1.d0-x**2))
        end if
      else
        sk0=1.d0
      end if
      end
