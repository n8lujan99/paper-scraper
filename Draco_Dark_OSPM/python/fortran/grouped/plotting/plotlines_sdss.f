
      parameter(nmax=10000)
      real x(nmax),y(nmax),yl(nmax),yh(nmax)

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=3500.
      xmax=5500.
      ymin=-0.4
      ymax=0.4
      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgenv(xmin,xmax,ymin,ymax,0,0)

      open(unit=1,file='out2',status='old')

      n=0
      do i=1,8000
         read(1,*,end=667) x1,x2,x3
         n=n+1
         x(n)=x1
         y(n)=x2-1.
         yl(n)=x2-x3-1.
         yh(n)=x2+x3-1.
      enddo
 667  continue
      close(1)

      call pgslw(5)
      call pgline(n,x,y)
      call pgsci(15)
      call pgslw(3)
      call pgline(n,x,yl)
      call pgline(n,x,yh)

      call pgend

      end
