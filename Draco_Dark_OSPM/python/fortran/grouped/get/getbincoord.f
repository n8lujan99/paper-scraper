
      parameter(nmax=100)
      parameter(radtodeg=57.29578)
      real r(nmax),v(nmax),rmid(nmax),vmid(nmax),xid(nmax)
      real rmid2(nmax)
      character fname*20,cvel*30

      rmax=1.
      rmax=400.
      rmax=1.2
      slitsize=4.
      open(unit=1,file='bin_r.out',status='old')
      read(1,*)
      read(1,*)
      nr=0
      do i=1,nmax
         read(1,*,end=666) i1,x2,x3,x4
         nr=nr+1
         r(nr)=x4
         rmid(nr)=x3
         rmid2(nr)=0.3*rmid(nr)
         xid(nr)=float(nr)
         irmax=i1
      enddo
 666  continue
      close(1)
      open(unit=1,file='bin_v.out',status='old')
      read(1,*)
      read(1,*)
      nv=0
      do i=1,nmax
         read(1,*,end=667) i1,x2,x3,x4
         nv=nv+1
         v(nv)=x4
         vmid(nv)=x3
         ivmax=i1
      enddo
 667  continue
      close(1)

      xs=0.25
      xfac=16.
      xoff=8.
      xnp=8.
      nb=0
      ymin=0.05
      ymax=0.95
      do i=1,9
         if(i.eq.1) then
            nb=nb+1
            write(*,1001) nb,0.5,ymin,0.,0.
         else
            do j=1,nv
               nb=nb+1
               x=rmid(i)*cos(vmid(j)/radtodeg)
               y=rmid(i)*sin(vmid(j)/radtodeg)
               xp=(x/xs+xoff)/xfac
               yp=(y/xs)/xfac
c               ip=nint(x/xs)
c               jp=nint(y/xs)
c               call xlinint(x,nr,rmid2,xid,xid0)
c               ip=nint(xid0)
c               call xlinint(y,nr,rmid2,xid,xid0)
c               jp=nint(xid0)
               xp=ymin+(ymax-ymin)*(float(i-1)/xnp+0.5)
               yp=ymin+(ymax-ymin)*float(j-1)/float(nv-1)
               write(*,1001) nb,xp,yp,x,y
            enddo
            do j=nv,1,-1
               nb=nb+1
               x=-rmid(i)*cos(vmid(j)/radtodeg)
               y=rmid(i)*sin(vmid(j)/radtodeg)
               xp=(x/xs+xoff)/xfac
               yp=(y/xs)/xfac
c               call xlinint(-x,nr,rmid2,xid,xid0)
c               ip=-nint(xid0)
c               call xlinint(y,nr,rmid2,xid,xid0)
c               jp=nint(xid0)
               xp=0.5-float(i-1)/xnp
               xp=ymin+(ymax-ymin)*(0.5-float(i-1)/xnp)
               yp=float(j-1)/float(nv-1)
               yp=ymin+(ymax-ymin)*float(j-1)/float(nv-1)
               write(*,1001) nb,xp,yp,x,y
            enddo
         endif
      enddo

 1001 format(i4,4(1x,f8.3))
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end
