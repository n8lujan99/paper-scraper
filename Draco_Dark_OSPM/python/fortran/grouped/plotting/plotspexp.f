
      parameter(nmax=10000)
      real x(nmax),y(nmax),y1(nmax),y2(nmax),y3(nmax)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=3500.
      xmax=5500.
      ymin=3000.
      ymax=32000.
      ymin=0.5
      ymax=1.5
      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgsch(1.5)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength','Counts','')
      call pgslw(4)
      call pgsch(1.)

      open(unit=1,file='j1',status='old')

      nl=0
      ic=0
      do il=1,10000
         read(1,*,end=667) x1,x2,x3,x4,x5
         n=n+1
         x(n)=x1
         y(n)=x2
         y1(n)=x3/x2
         y2(n)=x4/x2
         y3(n)=x5/x2
      enddo
 667  continue
      close(1)
c      call pgline(n,x,y)
      call pgsci(2)
      call pgline(n,x,y1)
      call pgsci(3)
      call pgline(n,x,y2)
      call pgsci(5)
      call pgline(n,x,y3)

      call pgend

      end
