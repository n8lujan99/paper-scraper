
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax)
      real yall(1000,nmax),ynl(nmax),ynu(nmax),xrat(100,1000)
      character file1*80,file2*80,c1*18

      call pgbegin(0,'?',1,2)
      call pgpap(0.,1.)
      call pgsch(1.7)
      call pgscf(2)
      call pgslw(2)

c      xmin=42.3
c      xmax=44.
c      ymin=log10(0.7)
c      ymax=log10(2.e3)
c      call pgenv(xmin,xmax,ymin,ymax,0,20)
c      call pglabel("log\D10\UL\DLy\ga","sn4.8/sn7.0","")
      xmin=9.
      xmax=31.
      ymin=-1.0
      ymax=3.0
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("f\Din","f\Din\U-f\Dout","")

      call pgslw(6)

      open(unit=1,file='list',status='old')

      icomp=9
      ic=0
      do il=1,10000
         read(1,*,end=666) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2
            n=n+1
            x(n)=x1
c            y(n)=log10(x2)
            y(n)=x2
         enddo
 667     continue
         close(2)
         call pgsci(ic)
         call pgline(n,x,y)
      enddo
 666  continue
      close(1)

      call pgslw(2)
      call pgsci(1)

      xmin=9.
      xmax=31.
      ymin=0.1
      ymax=0.3
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("f\Din","\gs/f\Dout","")

      call pgslw(6)

      open(unit=1,file='list',status='old')

      icomp=9
      ic=0
      do il=1,10000
         read(1,*,end=766) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,nmax
            read(2,*,end=767) x1,x2,x3
            n=n+1
            x(n)=x1
c            y(n)=log10(x2)
            y(n)=x3
         enddo
 767     continue
         close(2)
         call pgsci(ic)
         call pgline(n,x,y)
      enddo
 766  continue
      close(1)

      call pgend

      end
