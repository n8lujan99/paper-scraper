
      parameter(nmax=700000)
      real x(nmax),y(nmax),xl(2),yl(2),xl2(2),yl2(2)
      character file1*40,c1*3

      open(unit=1,file='listin',status='old')
      
      xmin=-2.5
      xmax=2.5
      xl(1)=xmin
      xl(2)=xmax
      yl(1)=0.
      yl(2)=0.
      xl2(1)=0.
      xl2(2)=0.
      yl2(1)=xmin
      yl2(2)=xmax
      
      call pgbegin(0,'?',3,3)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      do i=1,1000
         read(1,*,end=667) file1
         open(unit=2,file=file1,status='old')
         c1=file1(1:3)
         n=0
         do j=1,nmax
            read(2,*,end=666) x1,x2
            n=n+1
            x(n)=x1
            y(n)=x2
         enddo
 666     continue
         close(2)
         call pgsci(1)
         call pgsch(1.5)
         call pgenv(xmin,xmax,xmin,xmax,0,0)
         call pgsch(1.0)
         call pgpt(n,x,y,17)
         call pgsci(2)
         call pgline(2,xl,yl)
         call pgline(2,xl2,yl2)
         call biwgt(x,n,xb,xs)
         call biwgt(y,n,yb,xs)
         call pgsch(5.0)
         call pgslw(4)
         call pgsci(4)
         call pgpt1(xb,yb,5)
         call pgslw(2)
         call pgsci(1)
         call pgsch(2.5)
         call pgmtxt('T',0.8,0.5,0.5,c1)
      enddo
 667  continue
      close(1)
      
      call pgend
      end
