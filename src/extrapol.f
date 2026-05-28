      program extrapol
      implicit none
      real*8 alpha,kc,rc,kc_conv,dim_cutoff
      real*8 dummy,yp1,ypn
      character*128 qid
      logical ife
      integer nkact,la,n,lpx
      include 'geometry.cm'
      include 'system.cm'
      include 'kshell.cm'

      electrongas=.false.       ! atomic units only in this build
      kmax=5.d0

      write(*,*) 'rc*kc'
      read(*,'(f8.6)') dim_cutoff
      write(*,*) 'nkact'
      read(*,'(i3)') nkact
      write(*,*) 'rs'
      read(*,'(f8.6)') rs
      write(*,*) 'nelec'
      read(*,'(i4)') Nparts
      e2=1.d0
      rho=ndim/(2.d0*(ndim-1)*pi*rs**ndim)
      write(*,*) 'ell'
      vol=1.d0
      do la=1,ndim
        read(*,'(3f18.8)') ell(la)
        vol=vol*ell(la)
      end do

737   format((a8,f10.4))
738   format((a8,i4))
      write(*,737) 'rc*kc = ', dim_cutoff
      write(*,738) 'nkact = ', nkact
      write(*,737) 'rs = ', rs
      write(*,738) 'nparts = ', Nparts
      write(*,'(a8,3f18.8)') 'ell = ', ell
      write(*,737) 'vol =', vol
      write(*,737) 'rho = ', Nparts/vol
      write(*,737) 'rho_rs = ', rho

      rc=max(ell(1),ell(2),ell(3))/2.d0
      kc_conv=dim_cutoff/rc
      write(*,737) 'rc = ', rc
      write(*,737) 'rec. kc = ', kc_conv
      if (abs(rho-Nparts/vol) > 0.0001d0) then
        write(*,*) 'check density and box'
        stop
      endif
      call setbox(0)
      kc=rknorm(nkact+1)
      kmax=kc
      write(*,738) 'nkact = ', nkact
      write(*,737) 'kc = ', kc
      write(*,737) 'kmax = ', kmax
      if (kc.lt.kc_conv) then
        write(*,*) 'increase nkact'
        stop
      endif
      iksbig=1000

      write(*,*) ' file qid or stop ?'
      read(*,'(a)') qid
      lpx=index(qid,' ')-1
      if(qid(1:lpx).eq.'stop') stop
      inquire(file=qid(1:lpx)//'.sk',exist=ife)
      if(.not.ife) then
        write(*,*) 'file sk does not exist'
        stop
      end if

      open(10,file=qid(1:lpx)//'.sk',form='formatted',status='unknown')
      nsk=0
      rksk(0)=0.d0
      sk(0)=0.d0
67    continue
      read(10,*,end=68) rksk(nsk+1),sk(nsk+1),dummy
      nsk=nsk+1
      goto 67
68    continue
      close(10)
      if (rksk(nsk).lt.kmax) then
        write(6,*) ' kmax  '
        stop
      end if

      yp1=0.d0
      ypn=.999d30
      write(*,*) ' nsk= ',nsk
      call spline(rksk(0),sk(0),nsk+1,yp1,ypn,sk2(0))
      write(*,*) ' kc kmax ',kc,kmax

      del=kmax/iksbig
      do n=1,iksbig
        rkbig(n)=(n-0.5d0)*del
      end do

      alpha=sqrt(kc/2.d0/rc)   ! Ewald cutoff parameter

      call ewald(kc,alpha)

      end
