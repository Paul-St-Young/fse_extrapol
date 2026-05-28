      subroutine fillk(rknorm,wtk,kmult,nshlls,nshex,vol)
      implicit none
      include 'geometry.cm'
      integer nshlls,nshex
      integer kmult(0:nshlls)
      real*8 rknorm(0:nshex),wtk(0:nshex),vol
      integer k
      real*8 dk,cutk,con,vdown,vup

      wtk(0)=1.d0
      do k=1,nshlls
        wtk(k)=2.d0*(kmult(k)-kmult(k-1))
      enddo

      dk=rknorm(1)*0.5d0
      cutk=rknorm(nshlls)
      con=vol*2.d0*pi*(ndim-1)/(ndim*(2*pi)**ndim)
      vdown=1+2*kmult(nshlls)

      do k=nshlls+1,nshex
           vup=con*(cutk+dk*(k-nshlls))**ndim
           rknorm(k)=cutk+dk*(k-nshlls-0.5d0)
           wtk(k)=vup-vdown
           vdown=vup
      enddo

      return
      end
