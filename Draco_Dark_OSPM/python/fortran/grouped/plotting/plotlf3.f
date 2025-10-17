
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),x1a(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax),ycomp(nmax)
      real yall(1000,nmax),ynl(nmax),ynu(nmax),xrat(100,1000)
      real x1ab(nmax),yc(nmax)
      real xcompb(nmax),ycompb(nmax)
      real xratb(100,1000)
      character file1*80,file2*80,c1*18

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=42.3
      xmax=44.
      ymin=log10(60.)
      ymax=log10(1.5e6)

      call pgenv(xmin,xmax,ymin,ymax,0,20)
      call pglabel("log\D10\UL\DLy\ga","N","")

      open(unit=1,file='listin',status='old')

      call pgslw(5)
c      call pgsls(4)

      do il=1,10000
         read(1,*,end=766) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         nb=0
         do i=1,nmax
            read(2,*,end=767) x1,x2,x3,x4
            n=n+1
            x(n)=x1
            x2=max(0.01,x2)
            yc(n)=log10(x2)
            y(n)=log10(x4)
         enddo
 767     continue
         close(2)
         call pgsci(ic)
         call pgline(n,x,y)
      enddo
 766  continue
      rewind(1)
      call pgslw(15)
      call pgsci(1)
      call pgsls(1)
      call pgline(n,x,yc)

      call pgslw(5)
      call pgsls(1)
      do il=1,10000
         read(1,*,end=768) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         nb=0
         do i=1,nmax
            read(2,*,end=769) x1,x2,x3,x4
            n=n+1
            x(n)=x1
            rat=x2/x4
            x2=max(0.01,x2)
            y(n)=yc(n)-log10(rat)
         enddo
 769     continue
         close(2)
         call pgsci(ic)
c         call pgline(n,x,y)
      enddo
 768  continue
      close(1)

 1001 format(7(1x,f10.3))

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
