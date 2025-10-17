
      parameter(nmax=100000)
      real x(nmax),y(nmax),xa(nmax),ya(nmax)
      real xl(2),yl(2),yb(nmax),yc(nmax),yd(nmax)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xadd=2017.
      xmin=0.+xadd
      xmax=8.0+xadd
      ymin=0.
      ymax=490.
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Year',
     $     'IFU Fields Remaining (1000s)','')
c      call pglabel('Years since 1 Jan 2017',
c     $     'IFU Fields Remaining (1000s)','')

c      call pgsci(2)
      call pgsch(1.4)
      call pgptxt(7.5+xadd,400.,270.,0.,'Spring Complete')
      call pgptxt(6.8+xadd,400.,270.,0.,'Fall Complete')
      call pgsci(1)
      call pgsch(1.5)
      call pgsls(4)
      yl(1)=ymin
      yl(2)=ymax
      do i=1,7
         xl(1)=i+xadd
         xl(2)=i+xadd
         call pgline(2,xl,yl)
      enddo
      xl(1)=xmin
      xl(2)=xmax
      yl(1)=468.
      yl(2)=468.
      call pgsls(1)
      call pgline(2,xl,yl)

      open(unit=1,file='inyear',status='old')
      n=0
      do i=1,10000
         read(1,*,end=666) x1,x2
         n=n+1
         x(n)=x1+xadd
         y(n)=x2
      enddo
 666  continue
      close(1)

      call pgslw(6)
      call pgline(n,x,y)

      open(unit=1,file='outyear',status='old')
      n=0
      do i=1,10000
         read(1,*,end=667) i1,x2,x3
         n=n+1
         x(n)=x2+1.+xadd
         y(n)=x3-1.
      enddo
 667  continue
      close(1)

      call pgsci(2)
      call pgline(n,x,y)

      call pgend

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
