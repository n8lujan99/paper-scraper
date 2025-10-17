
      real xin(1000),yin(1000),yin2(1000)
      real wadc(5),adc(5),xl(2),yl(2)
      character c1*5,title*10
      parameter(pi=3.141593e0)

      imoff=0
      fwhm=1.34
      parang=292.
      az0=270.+parang
      az0=az0+90.
      az0=parang
      if(az0.gt.360.) az0=az0-360.
      if(az0.lt.0.) az0=az0+360.
      print *,az0
c      az0=318.5
c      az0=90.
      rsig0=fwhm/2.355
      dtr=180./pi
      rfib=0.75

      wadc(1)=3500.
      wadc(2)=4000.
      wadc(3)=4500.
      wadc(4)=5000.
      wadc(5)=5500.
      adc(1)=-0.71
      adc(2)=-0.34
      adc(3)=-0.085
      adc(4)=0.08
      adc(5)=0.20

      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.1)
      call pgslw(2)

      xmin=-4.
      xmax=4.

      xrs=0.
      xds=0.
      nstep=10
      xstep=2.*rfib/float(nstep-1)
      deltx=pi*rfib*rfib
      area=1.*deltx/(2.*rsig0*rsig0*pi)
c      area=area/0.90

      open(unit=2,file='coords.in',status='old')
      do iaoff=1,5,2
      ia=iaoff
      call pgsfs(1)
      call pgsci(1)
      call pgsch(1.1)
      call pgenv(xmin,xmax,xmin,xmax,0,0)
      call pgsch(1.3)
      if(ia.eq.1) title="3500"
      if(ia.eq.3) title="4500"
      if(ia.eq.5) title="5500"
      call pglabel("\gDRA","\gDDec","")
      call pgsch(1.6)
      xl(1)=xmin
      xl(2)=xmax
      yl(1)=0.
      yl(2)=0.
      call pgline(2,xl,yl)
      yl(1)=xmin
      yl(2)=xmax
      xl(1)=0.
      xl(2)=0.
      call pgline(2,xl,yl)
      call pgsfs(2)
      xaoff=adc(ia)*sin(az0/dtr)
      yaoff=adc(ia)*cos(az0/dtr)
      nin=0
      sumgw=0.
      do j=1,100
         read(2,*,end=667) x1,x2,i3
         call pgslw(4)
         call pgsci(i3)
         call pgcirc(x1+xaoff,x2+yaoff,rfib)
         gaus=0.
         xmoff=0.
         nsum=0
         xs=x1-rfib+xaoff
         ys=x2-rfib+yaoff
         do ix=1,nstep
            xp=xs+xstep*float(ix-1)
            do iy=1,nstep
               yp=ys+xstep*float(iy-1)
               dist0=sqrt((xp-x1-xaoff)**2+(yp-x2-yaoff)**2)
c               dist0=sqrt((xp-x1)**2+(yp-x2)**2)
               if(dist0.lt.rfib) then
                  dist=sqrt((xp-xrs)**2+(yp-xds)**2)
                  g=dist/rsig0
                  gaus=gaus+exp(-g*g/2.)*area
                  nsum=nsum+1
               endif
            enddo
         enddo
         gaus=gaus/float(nsum)
         xmoff=xmoff/float(nsum)
         if(imoff.eq.1) gaus=xmoff
         sumgw=sumgw+gaus
c         print *,gaus,xaoff,yaoff
         write(c1,1001) gaus
         call pgptxt(x1+xaoff,x2+yaoff,0.,0.5,c1)
         call pgsci(1)
         call pgslw(2)
         call pgptxt(-3.5,3.3,0.,0.,title)
      enddo
 667  continue
      rewind(2)
      print *,sumgw
      enddo
      close(2)

      call pgend
 1001 format(f5.2)

      end
