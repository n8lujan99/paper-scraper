
      parameter(nmax=10000)
      real x(nmax),yl(nmax),yo(nmax)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=4.8
      xmax=7.0
      ymin=0.
      ymax=5.0
      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgsch(1.5)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('S/N','Detections per IFU','')
      call pgslw(5)

      open(unit=1,file='ndet.txt',status='old')

      n=0
      do il=1,10000
         read(1,*,end=666) x1,x2,x3
         n=n+1
         x(n)=x1
         yl(n)=x3
         yo(n)=x2
      enddo
 666  continue
      close(1)

      call pgsci(1)
      call pgline(n,x,yl)
      call pgptxt(6.7,0.4,0.,0.,'LAE')
      call pgsci(2)
      call pgline(n,x,yo)
      call pgptxt(6.7,1.5,0.,0.,'ALL')
      x(1)=xmin
      x(2)=xmax
      yl(1)=2.5
      yl(2)=2.5
      yo(1)=4.1
      yo(2)=4.1
      call pgsls(4)
      call pgsci(1)
      call pgline(2,x,yl)
      call pgsci(2)
      call pgline(2,x,yo)

      call pgend

      end
