
      parameter(nmax=12000)
      real x(nmax),y(nmax),w(nmax),ys(nmax,nmax)
      real yp(nmax),yin(nmax),yp2(nmax),y2(nmax),xnorm(nmax)
      real yloc(nmax),ylocs(nmax,nmax),xl(2),yl(2),xp(nmax)
      real ylocp(nmax),ylocin(nmax),yloc2(nmax)
      character file1*60

      fluxfac=2.
      cont=0.

      do i=1,1001
         w(i)=3500.+2.*float(i-1)
      enddo
      
      open(unit=1,file='list',status='old')

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)
      
      xmin=3500.
      xmax=5500.
      ymin=0.
      ymax=1.3

      n=0
      do i=1,nmax
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old',err=667)
         read(2,*)
         ns=0
         n=n+1
         do j=1,nmax
            read(2,*,end=667) x1,x2
            read(2,*,end=667) x1b,x2b
            ns=ns+1
            x(ns)=(x1+x1b)/2.
            y(ns)=(x2+x2b)/2.
c            y(ns)=(x2+x2b)/2.*x(ns)
         enddo
 667     continue
         close(2)
         do iw=1,1001
            call xlinint0(w(iw),ns,x,y,y0)
            ys(iw,n)=y0
         enddo
      enddo
 666  continue
      close(1)

      wmin=3500.
      wmax=4200.
      do i=1,n
         nin=0
         do iw=1,ns
            if(x(iw).gt.wmin.and.x(iw).lt.wmax.and.ys(iw,i).ne.0.) then
               nin=nin+1
               yin(nin)=ys(iw,i)
            endif
         enddo
         call biwgt(yin,nin,xb,xs)
         xnorm(i)=xb
      enddo

      open(unit=11,file='out',status='unknown')
      np=0
      do iw=1,ns
         nin=0
         do i=1,n
            if(ys(iw,i).ne.0) then
               nin=nin+1
               yin(nin)=ys(iw,i)/xnorm(i)
            endif
         enddo
         call biwgt(yin,nin,xb,xs)
         if(nin.gt.0) then
            np=np+1
            xp(np)=x(iw)
            yp(np)=xb
            write(11,*) x(iw),yp(np),xs,nin
         endif
      enddo
      close(11)

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("Wavelength","Relative Flux","")

      call pgsci(1)
      call pgline(np,xp,yp)

      end

      subroutine xlinint0(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.le.x(1)) yp=0.
      if(xp.ge.x(n)) yp=0.
      return
      end
