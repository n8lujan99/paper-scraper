
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax),xp(10),yp(10)
      character file1*80,file2*80,c1*18

      xp(1)=30.
      yp(1)=35.
      xp(2)=175.
      yp(2)=35.
      xp(3)=220.
      yp(3)=42.

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=42.3
      xmax=43.5
      ymin=1.
      ymax=9.
      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgsch(1.5)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('L','vol_5.5/6.5','')
      call pgslw(4)
      call pgsch(1.)

      open(unit=1,file='df1',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3
         n=n+1
         x(n)=x1
         y(n)=x3
      enddo
 666  continue
      close(1)
      call pgline(n,x,y)

      open(unit=1,file='df2',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3
         n=n+1
         x(n)=x1
         y(n)=x3
      enddo
 667  continue
      close(1)
      call pgsci(2)
      call pgline(n,x,y)

      call pgend

      end
