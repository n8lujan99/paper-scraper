
      parameter(nmax=10000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(1000)
      character file1*80,file2*80,c1*3

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=0.
      xmax=113.
      ymin=0.96
      ymax=1.1

      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Fiber','noise ratio','')

      open(unit=1,file='list',status='old')

      na=0
      ic=0
      do ia=1,1000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         n=0
         nin=0
         do i=1,8000
            read(2,*,end=667) xd1,x1,x2
            if(x2.gt.0.5.and.x2.lt.2.0) then
               n=n+1
               x(n)=x1
               y(n)=x2
               if(x1.gt.15.and.x1.lt.90) then
                  nin=nin+1
                  xin(nin)=x2
               endif
            endif
         enddo
 667     continue
         close(2)
         call biwgt(xin,nin,xb,xs)
         ymax=0.
         do i=1,n
            y(i)=y(i)/xb
            if(x(i).gt.15.and.x(i).lt.90) ymax=max(ymax,y(i))
         enddo
         print *,file1(1:8),xb,xs,ymax
         if(ymax.lt.1.06) then
            ic=ic+1
            if(ic.eq.13) ic=1
            call pgsci(ic)
            call pgline(n,x,y)
            call pgslw(1)
         endif
      enddo
 666  continue
      close(1)

      call pgend

      end
