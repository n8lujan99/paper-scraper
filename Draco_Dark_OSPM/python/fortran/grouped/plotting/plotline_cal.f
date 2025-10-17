
      parameter(nmax=10000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(100)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      real we(nmax),ee(nmax),y2(nmax),ys(nmax),xl(2),yl(2)
      character file1*80,file2*80,c1*3

      call pgbegin(0,'?',1,2)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=3500.
      xmax=5500.
      ymin=0.94
      ymax=1.06
      call pgsls(1)
      call pgslw(3)

      open(unit=1,file='extinction',status='old')
      ne=0
      do i=1,nmax
         read(1,*,end=888) x1,x2
         ne=ne+1
         we(ne)=x1
         ee(ne)=x2
      enddo
 888  continue
      close(1)

c      open(unit=2,file="s1b",status='old')
      open(unit=2,file="sdsscomp.dat",status='old')
      open(unit=11,file="out",status='unknown')
      call pgsci(1)
      call pgsch(2.0)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength (\(2078))','DEX/SDSS','')
      call pgsch(1.5)
      xl(1)=xmin
      xl(2)=xmax
      yl(1)=1.
      yl(2)=1.
      call pgslw(2)
      call pgsls(4)
      call pgline(2,xl,yl)
      call pgsls(1)
      n=0
      do i=1,8000
         read(2,*,end=667) x1,x2,x3
         n=n+1
         x(n)=x1
c         y(n)=1.+x2
         y(n)=x2
         ys(n)=x3
         call xlinint(x(n),ne,we,ee,ev)
         y2(n)=y(n)*ev
c         write(11,1101) x(n),y2(n)
      enddo
 667  continue
      close(2)
      
      do i=1,1036
         wave=3470.+float(i-1)*2.
         call xlinint(wave,n,x,y2,yv)
         write(11,1101) wave,yv
      enddo
      close(11)
 1101 format(f6.1,1x,f6.3)
      call pgslw(6)
      call pgline(n,x,y)
      call pgsci(2)
c      call pgline(n,x,y2)
      call pgslw(3)

      ymin=0.
      ymax=0.2
c      open(unit=2,file="r1",status='old')
      call pgsci(1)
      call pgsch(2.0)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength (\(2078))','Scatter','')
      call pgsch(1.5)
c      n=0
c      do i=1,8000
c         read(2,*,end=668) x1,x2
c         n=n+1
c         x(n)=x1
c         y(n)=x2
c      enddo
c 668  continue
c      close(2)
      call pgslw(6)
      call pgline(n,x,ys)

      open(unit=2,file="ball2",status='old')
      call pgsch(1.5)
      n=0
      do i=1,8000
         read(2,*,end=669) x1,x2
         n=n+1
         x(n)=x1
         y(n)=x2
      enddo
 669  continue
      close(2)
      call pgsci(2)
      call pgline(n,x,y)

      open(unit=2,file="repeat.dat",status='old')
      call pgsch(1.5)
      n=0
      do i=1,8000
         read(2,*,end=670) x1,x2,x3
         n=n+1
         x(n)=x1
         y(n)=x3
      enddo
 670  continue
      close(2)
      call pgsci(4)
      call pgline(n,x,y)

      n=2
      x(1)=xmin
      x(2)=xmax
      y(1)=0.134
      y(2)=y(1)
      call pgsci(4)
c      call pgline(2,x,y)

      close(1)

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
