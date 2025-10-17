
      parameter(nmax=500000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(100)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      real xl(2),yl(2)
      character file1*80,file2*80,c1*3

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xl(1)=50.
      xl(2)=190.
      yl(1)=4.
      yl(2)=0.005

      yl(1)=0.14
      yl(2)=0.01

      xmin=75.
      xmax=220.
      xmin=75.
      xmax=190.
      ymin=-8.
      ymax=35.
      ymax=90.

      xmin=80.
      xmax=180.
      xmin=60.
      xmax=230.
      ymin=0.
c      ymax=0.15
      ymax=0.013

      ymin=-7.
      ymax=29.
c      ymin=-200.
c     ymax=400.

      open(unit=1,file='limits.dat',status='old',err=555)
      read(1,*) xmin,xmax,ymin,ymax
 555  continue
      close(1)
      
      call pgsls(1)
      call pgslw(2)

      open(unit=1,file='linelist',status='old')

      nl=0
      ic=0
      do ia=1,700
         read(1,*,end=666) file1,ic
         itype=1
         if(ic.lt.0) then
            ic=-ic
            itype=4
         endif
         open(unit=2,file=file1,status='old')
         if(ia.eq.1) then
            call pgsci(1)
            call pgenv(xmin,xmax,ymin,ymax,0,0)
c            call pglabel('Mpc','s\U2\DX','')
c            call pglabel('Mpc','DD/RR, DD/DR','')
            call pglabel('Mpc','DD/DR-1','')
c            call pglabel('wavelength','S/N cut','')
c            call pglabel('Mpc','','')
            call pgsch(2.0)
            call pgsch(1.5)
         endif
         n=0
         do i=1,500000
            read(2,*,end=667) x1,x2
            n=n+1
            x(n)=x1
            y(n)=x2
            if(ia.ne.3) then
               call xlinint(x1,2,xl,yl,yc)
c               y(n)=y(n)*yc
c               y(n)=y(n)-yc
            endif
         enddo
 667     continue
         close(2)
         call pgslw(10)
c         call pgslw(3)
         call pgsci(ic)
         call pgsls(itype)
         call pgline(n,x,y)
         call pgsls(1)
         call pgslw(1)
      enddo
 666  continue
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
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
