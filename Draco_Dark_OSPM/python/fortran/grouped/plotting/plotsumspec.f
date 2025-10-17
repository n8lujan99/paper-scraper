
      parameter(nmax=12000)
      real x(nmax),y(nmax),w(nmax),ys(nmax,nmax)
      real yp(nmax),yin(nmax),yp2(nmax),y2(nmax)
      real yloc(nmax),ylocs(nmax,nmax),xl(2),yl(2)
      real ylocp(nmax),ylocin(nmax),yloc2(nmax),wfix(nmax),ffix(nmax)
      character file1*60

      icont=0
      fluxfac=2.
      if(icont.eq.0) then
         wla=1215.67
         wc=100.
         wbin=1.
         wmin=wla-wc
         wmax=wmin+wc*2.
c         wmin=775.
c         wmax=1500.
         xntry=1.0
c         xntry=0.85
c         wla=4500.
c         wmin=3500.
c         wmax=5500.
c         wbin=2.
      else
         wla=4500.
         wc=1000.
         wbin=2.
         wmin=3500.
         wmax=5500.
         xntry=0.85
c         xntry=1.
      endif

      wmino=370.
      wmaxo=40000.
c      read *,wmino,wmaxo
      nw=nint((wmax-wmin)/wbin)+1
      do i=1,nw
         w(i)=wmin+float(i-1)*(wmax-wmin)/float(nw-1)
      enddo
      
      open(unit=1,file='list',status='old')

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)
      
      xmin=wmin
      xmax=wmax
      ymin=-0.1
      ymax=0.5
      ymin=-1
      ymax=2.5

      open(unit=2,file='sky_all.use',status='old')
      nfix=0
      do i=1,nmax
         read(2,*,end=555) x1,x2
         nfix=nfix+1
         wfix(nfix)=x1
         ffix(nfix)=x2
      enddo
 555  continue
      close(2)

      wsum=4.

      n=0
      do i=1,nmax
         read(1,*,end=666) file1,wave
c         wave=4500.
         open(unit=2,file=file1,status='old',err=777)
         ns=0
         sumo=0.
         sumr=0.
         n=n+1
         z=wave/wla-1.
         do j=1,nmax
c            read(2,*,end=667) x1,x2,x3
            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
c            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11
            w0=x1*wla/wave
c            if(w0.ge.wmin.and.w0.le.wmax) then
            if(w0.ge.wmin.and.w0.le.wmax.and.
     $           x1.gt.wmino.and.x1.lt.wmaxo) then
               if(x9.le.0.5) then
                  ns=1
                  goto 667
               endif
               call xlinint(x1,nfix,wfix,ffix,f0)
               ns=ns+1
               x(ns)=w0
               y(ns)=x12/fluxfac/x9
c               y(ns)=y(ns)/fluxfac/x9/f0
               yloc(ns)=x2/fluxfac/x9
               if(ns.gt.1) then
                  sumr=sumr+x(ns)-x(ns-1)
                  sumo=sumo+x1-x1old
               endif
               x1old=x1
c               print *,ns,x(ns),y(ns)
            endif
         enddo
 667     continue
         wsumr=wsum*wla/wave
         if(ns.gt.10) then
            sumr=sumr/float(ns-1)
            sumo=sumo/float(ns-1)
            frac=wave/wla
            frac=frac*wbin/sumr
            do j=1,nw
               call xlinint0(w(j),ns,x,y,y0)
               call xlinint0(w(j),ns,x,yloc,yloc0)
c               ys(j,n)=y0*frac
               ys(j,n)=y0
               y2(j)=ys(j,n)
               ylocs(j,n)=yloc0
               yloc2(j)=ylocs(j,n)
            enddo
c            call pgsci(1)
c            call pgenv(xmin,xmax,ymin,ymax,0,0)
c            call pgline(ns,x,y)
c            call pgsci(2)
c            call pgline(ns,x,yloc)
         else
            n=n-1
         endif
 777     continue
         close(2)
      enddo
 666  continue
      close(1)

      print *,n
      open(unit=11,file='out',status='unknown')
      do iw=1,nw
         nin=0
         do i=1,n
            if(ys(iw,i).ne.0.) then
               nin=nin+1
               yin(nin)=ys(iw,i)
               ylocin(nin)=ylocs(iw,i)
            endif
         enddo
         call biwgt(yin,nin,xb,xs)
         call biwgt(ylocin,nin,xlocb,xlocs)
         ntry=nint(xntry*float(nin))
         call biwgt(yin,ntry,xb,xs)
         call biwgt(ylocin,ntry,xlocb,xlocs)
         yp(iw)=xb
         ylocp(iw)=xlocb
         yp2(iw)=sum
         write(11,*) w(iw),yp(iw),ylocp(iw),xs,n
      enddo
      close(11)
c      close(12)

      call pgenv(xmin,xmax,ymin,ymax,0,0)
c      call pgenv(xmin,3900.,ymin,ymax,0,0)
      call pglabel("Wavelength","Flux","")
      xl(1)=xmin
      xl(2)=xmax
      yl(1)=0.
      yl(2)=0.
      call pgline(2,xl,yl)

      call pgslw(6)
      call pgsci(2)
      call pgline(nw,w,ylocp)
      call pgsci(1)
      call pgline(nw,w,yp)

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
