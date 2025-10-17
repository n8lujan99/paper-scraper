
      real x(1000),y(1000),yl(1000),yh(1000)
      character file1*40,cx*6,cy*6

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=1,file='pos.dat',status='old')
      do i=1,1000
         read(1,*,end=666) file1,xp,yp,xg,yg
         open(unit=2,file=file1,status='old')
         n=0
         do j=1,1000
            read(2,*,end=667) x1,x2,x3,x4
            n=n+1
            x(n)=x1
            y(n)=x2
            yl(n)=x3
            yh(n)=x4
         enddo
 667     continue
         close(2)

         xmin=xp-0.05
         xmax=xmin+0.1
         ymin=yp/1.5
         ymax=ymin+0.1
         call pgsvp(xmin,xmax,ymin,ymax)
         call pgswin(-1200.,1200.,0.,0.15)
         if(i.lt.5) then
            call pgsci(2)
         else
            call pgsci(1)
         endif
         call pgbox('bc',0.,0,'bc',0.,0)
         call pgline(n,x,y)
         call pgsls(4)
         call pgline(n,x,yl)
         call pgline(n,x,yh)
         call pgsls(1)
         write(cx,1001) xg
         write(cy,1001) yg
         call pgsch(0.5)
         call pgmtxt('T',-1.3,0.2,0.5,cx)
         call pgmtxt('T',-1.3,0.6,0.5,cy)
         call pgsch(1.5)
         
      enddo
 666  continue
      close(1)

      call pgend
 1001 format(f6.2)
      end
