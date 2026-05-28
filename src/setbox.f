      subroutine setbox(iwrite)
c     Compute vol, el2, elli, tpiell from ell and enumerate the
c     reciprocal-lattice shells with cut = 15 * tpiell(1).
      implicit none
      include 'geometry.cm'
      include 'system.cm'
      include 'kshell.cm'
      real*8 kcut
      integer la,k,iwrite

      vol=1.d0
      do la=1,ndim
        vol=vol*ell(la)
        el2(la)=ell(la)*0.5d0
        elli(la)=1.d0/ell(la)
        tpiell(la)=2.d0*pi/ell(la)
      end do
      kcut=15.d0*tpiell(1)
      call shells(ndim,tpiell,kcut,Nshlls,rkcomp,rknorm,kmult
     +   ,Nvects,Mnkv,Mnshlls,ndim,iwrite)
      if(iwrite.gt.0) then
        write(6,*)' k   kmult(k)   rknorm(k)'
        do k=0,Nshlls
          write (6,'(2i5,f12.4)') k,kmult(k),rknorm(k)
        enddo
      end if
      end
