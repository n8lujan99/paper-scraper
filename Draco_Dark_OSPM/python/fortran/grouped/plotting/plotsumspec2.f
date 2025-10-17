
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax)
      real yall(1000,nmax),ynl(nmax),ynu(nmax),xrat(100,1000)
      character file1*80,file2*80,c1*18

      yoff=0.3
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.7)
      call pgscf(2)
      call pgslw(2)

c      xcen=1050
      xcen=1215.66
      xrange=70.
      xrange=80.
      xmin=xcen-xrange
      xmax=xcen+xrange
c      xmin=775.
c      xmax=1500.
c      xmin=3500.
c      xmax=3900.
      ymin=-0.18
      ymax=0.6
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      x(1)=xmin
      x(2)=xmax
      y(1)=0.
      y(2)=0.
      call pgline(2,x,y)
      call pglabel("Wavelength","flux","")

      call pgslw(6)

      open(unit=1,file='listin',status='old')

      do il=1,10000
         read(1,*,end=666) file1,ic,iw
         open(unit=2,file=file1,status='old',err=668)
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3
            n=n+1
            x(n)=x1
c            y(n)=x2
            y(n)=x3+yoff
         enddo
 667     continue
         call pgsci(ic)
         call pgslw(iw)
         call pgline(n,x,y)
 668     close(2)
      enddo
 666  continue
      close(1)

      call pgend

      end
