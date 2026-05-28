      real*8 function vewald_lr_k(k,kcut,alpha)
c     Long-range part of the Ewald potential in k-space.
c     3D: gaussian decay  v_lr(k) = (2(d-1)pi/k^{d-1}) exp(-k^2/(4 alpha^2))
c     2D: erfc decay      v_lr(k) = (2(d-1)pi/k^{d-1}) erfc(k/(2 alpha))
      implicit none
      include 'geometry.cm'
      real*8 k,kcut,alpha,x
      vewald_lr_k=0.d0
      if(k.eq.0.d0) then
        if(ndim.eq.3) then
          vewald_lr_k=-pi/alpha**2
        elseif(ndim.eq.2) then
          vewald_lr_k=-sqrt(pi)/alpha
        end if
      else
        x=k/(2.d0*alpha)
        if(ndim.eq.3) then
          vewald_lr_k=2.d0*(ndim-1.d0)*pi/k**(ndim-1)*exp(-x**2)
        else if(ndim.eq.2) then
          vewald_lr_k=2.d0*(ndim-1.d0)*pi/k**(ndim-1)*erfc(x)
        end if
      end if
      end
