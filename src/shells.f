      subroutine shells(ndim,a,cut,nshlls,rkcomp,rknorm,kmult
     +,nvects,mnkv,mnsh,mdim,iwrite)
      implicit none
      integer ndim,nshlls,mnsh,kmult(0:mnsh),nvects,mnkv,npts,l
     +  ,nkspan(3),i,icount(3),nzero,j,kj,jp,mdim,iwrite
      real*8 a(mdim),rkcomp(mdim,mnkv),rknorm(0:mnsh)
     +,x(3),cut,c2,rsq,rks,rnow

      kmult(0)=0
      rknorm(0)=0.d0
      c2=cut**2
      npts=1
      do l=1,ndim
        nkspan(l)=(0.00001+abs(cut/a(l)))
        icount(l)=-nkspan(l)
        npts=(2*nkspan(l)+1)*npts
      enddo

      nvects=0
      do 3 i=1,npts
        rsq=0.0d0
        nzero=0
        do 4 l=1,ndim
          if(nzero.eq.0) nzero=icount(l)
          x(l)=icount(l)*a(l)
          rsq=rsq+x(l)**2
          if(rsq.gt.c2) go to 30
4       continue
        if(nzero.le.0) go to 30
        if(nvects.gt.mnkv) then
          write (6,*)' shells: mnkv too small ',mnkv, nvects
          stop
        endif
        nvects=nvects+1
        do 6 j=1,nvects-1
          kj=j
          rks=0.d0
          do 66 l=1,ndim
66          rks=rks+rkcomp(l,j)**2
6         if(rks.ge.rsq*(1.d0+1.d-5)) go to 7
        kj=nvects
7       do 8 jp=nvects,kj+1,-1
          do 8 l=1,ndim
8           rkcomp(l,jp)=rkcomp(l,jp-1)
        do l=1,ndim
          rkcomp(l,kj)=x(l)
        enddo
30      do 31 l=1,ndim
          icount(l)=icount(l)+1
          if(icount(l).le.nkspan(l)) go to 3
31        icount(l)=-nkspan(l)
3     continue

      nshlls=0
      rnow=0.d0
      do i=1,nvects
        rsq=0.d0
        do l=1,ndim
          rsq=rsq+rkcomp(l,i)**2
        enddo
        if(rsq-rnow.gt.0.0001d0*rnow) nshlls=nshlls+1
        if(nshlls.gt.mnsh) then
          write (6,*)'shells: mnsh too small ',mnsh,nshlls
          stop
        endif
        rnow=rsq
        rknorm(nshlls)=sqrt(rnow)
        kmult(nshlls)=i
      enddo
      if(iwrite.gt.0) then
        open(92,file='kshells.dat')
        write(92,*) 'nshlls=',nshlls
        do i=0,nshlls
          write (92,'(2i5,f12.4)') i,kmult(i),rknorm(i)
        enddo
        close(92)
      end if

      return
      end
