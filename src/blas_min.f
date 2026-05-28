c     Minimal BLAS routines needed by the vendored LINPACK shims
c     dpoco / dposl.  Hand-rolled to avoid a system BLAS dependency.

      double precision function ddot(n,dx,incx,dy,incy)
      implicit none
      integer n,incx,incy,i,ix,iy
      double precision dx(*),dy(*),s
      ddot=0.d0
      if(n.le.0) return
      s=0.d0
      if(incx.eq.1 .and. incy.eq.1) then
        do i=1,n
          s=s+dx(i)*dy(i)
        end do
      else
        ix=1
        iy=1
        if(incx.lt.0) ix=(-n+1)*incx+1
        if(incy.lt.0) iy=(-n+1)*incy+1
        do i=1,n
          s=s+dx(ix)*dy(iy)
          ix=ix+incx
          iy=iy+incy
        end do
      end if
      ddot=s
      return
      end

      double precision function dasum(n,dx,incx)
      implicit none
      integer n,incx,i,ix
      double precision dx(*),s
      dasum=0.d0
      if(n.le.0) return
      s=0.d0
      if(incx.eq.1) then
        do i=1,n
          s=s+abs(dx(i))
        end do
      else
        ix=1
        if(incx.lt.0) ix=(-n+1)*incx+1
        do i=1,n
          s=s+abs(dx(ix))
          ix=ix+incx
        end do
      end if
      dasum=s
      return
      end

      subroutine dscal(n,da,dx,incx)
      implicit none
      integer n,incx,i,ix
      double precision da,dx(*)
      if(n.le.0) return
      if(incx.eq.1) then
        do i=1,n
          dx(i)=da*dx(i)
        end do
      else
        ix=1
        if(incx.lt.0) ix=(-n+1)*incx+1
        do i=1,n
          dx(ix)=da*dx(ix)
          ix=ix+incx
        end do
      end if
      return
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
      implicit none
      integer n,incx,incy,i,ix,iy
      double precision da,dx(*),dy(*)
      if(n.le.0) return
      if(da.eq.0.d0) return
      if(incx.eq.1 .and. incy.eq.1) then
        do i=1,n
          dy(i)=dy(i)+da*dx(i)
        end do
      else
        ix=1
        iy=1
        if(incx.lt.0) ix=(-n+1)*incx+1
        if(incy.lt.0) iy=(-n+1)*incy+1
        do i=1,n
          dy(iy)=dy(iy)+da*dx(ix)
          ix=ix+incx
          iy=iy+incy
        end do
      end if
      return
      end
