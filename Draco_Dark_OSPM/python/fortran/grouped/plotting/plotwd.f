
      parameter(nmax=12000)
      real x(nmax),y(nmax),w(nmax),ys(nmax,nmax)
      real yp(nmax),yin(nmax),yp2(nmax),y2(nmax),xnorm(nmax)
      real yloc(nmax),ylocs(nmax,nmax),xl(2),yl(2)
      real ylocp(nmax),ylocin(nmax),yloc2(nmax),wfix(nmax),ffix(nmax)
      character file1*60

      fluxfac=2.
      cont=0.
      
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

      open(unit=2,file='fix1.dat',status='old')
      nfix=0
      do i=1,nmax
         read(2,*,end=555) x1,x2
         nfix=nfix+1
         wfix(nfix)=x1
         ffix(nfix)=x2
      enddo
 555  continue
      close(2)

      n=0
      do i=1,nmax
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old',err=667)
         ns=0
         n=n+1
         do j=1,nmax
            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9
c            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
            call xlinint(x1,nfix,wfix,ffix,f0)
            ns=ns+1
            x(ns)=x1
            y(ns)=x2/fluxfac/x9
            ys(ns,i)=y(ns)
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

      wmin=3800.
      wmax=4200.
      do i=1,n
         nin=0
         do iw=1,ns
            if(x(iw).gt.wmin.and.x(iw).lt.wmax) then
               nin=nin+1
               yin(nin)=ys(iw,i)
            endif
         enddo
         call biwgt(yin,nin,xb,xs)
         xnorm(i)=xb
      enddo


      open(unit=11,file='out',status='unknown')
      do iw=1,ns
         do i=1,n
            yin(i)=ys(iw,i)/xnorm(i)
         enddo
         call biwgt(yin,n,xb,xs)
         call xlinint(x(iw),nfix,wfix,ffix,f0)
         yp(iw)=xb
         yp2(iw)=xb/f0
         write(11,*) x(iw),yp(iw),xs,n
      enddo
      close(11)

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("Wavelength","Relative Flux","")

      call pgsci(1)
      call pgline(ns,x,yp)
      call pgsci(2)
      call pgline(ns,x,yp2)

      open(unit=1,file='out_WD',status='old')
      ns=0
      do i=1,nmax
         read(1,*,end=755) x1,x2
         ns=ns+1
         x(ns)=x1
         yp(ns)=x2
      enddo
 755  continue
      close(1)
      call pgsci(4)
      call pgline(ns,x,yp)

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
