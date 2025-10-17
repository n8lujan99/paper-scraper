Marked
      parameter(nstot=10,nftot=10,nmax=100,ncmax=1000000)
      real fin(nstot,nftot,nmax),hin(nstot,nftot,nmax),sn(nstot)
      real hin1(nstot,nftot,nmax),hin2(nstot,nftot,nmax)
      real hin3(nstot,nftot,nmax),hin4(nstot,nftot,nmax)
      real xlhist(nmax),flhist(nmax),fout(nstot,nftot)
      real zl(10000),xl(10000),xfc(ncmax),xsc(ncmax),xzc(ncmax)
      real sncon(100),rcon(100),xfce(ncmax),x5c(ncmax),yin(10),xin(10)
      real xwc(200),xic(200),xhc(200),xcc(200,100),xccu(200)
      real we(nmax),xe(nmax),wave1(nmax),fcor1(nmax),ratflux(ncmax)
      real wave2(10000),fcor2(10000),xf1sa(100),xif1s(100),xf1s(ncmax)
      real*8 dx2,dflhist(nmax),df5hist(nmax),dfnc(nmax)
      integer*8 i6,i7,iddet(ncmax)
      character file1*40

      read *,sncut

c- set to 0 to apply correction, 1 for no correction
      inocor=0

c- get the flux to luminosity as a function of z
      open(unit=1,file='FtoLz2.dat',status='old')
      nz=0
      do i=1,10000
         read(1,*,end=668) x1,dx2
         nz=nz+1
         zl(nz)=x1
         xl(nz)=sngl(dx2/1.d40)
c         xl(nz)=(1.+x1)*sngl(dx2/1.d40)
      enddo
 668  continue
      close(1)

c- get the confirmed rate as a funtion of s/n
      open(unit=1,file='conrate.dat',status='old')
      read(1,*)
      ncon=0
      do i=1,10000
         read(1,*,end=670) x1,x2
         ncon=ncon+1
         sncon(ncon)=x1
c         x2=1.
         rcon(ncon)=x2
         if(inocor.eq.1) rcon(ncon)=1.
      enddo
 670  continue
      close(1)

c- get the correction to the 1-sigma noise (divide by this value)
      open(unit=1,file='fcor1sig_hdr3',status='old')
      nfcor1=0
      do i=1,10000
         read(1,*,end=671) x1,x2
         nfcor1=nfcor1+1
         wave1(nfcor1)=x1
         fcor1(nfcor1)=x2
c         fcor1(nfcor1)=x2/1.1
      enddo
 671  continue
      close(1)

c- get the correction to the original extinction
      open(unit=1,file='extinction',status='old')
      nfcor2=0
      do i=1,10000
         read(1,*,end=672) x1,x2
         nfcor2=nfcor2+1
         wave2(nfcor2)=x1
c         fcor2(nfcor2)=x2
         fcor2(nfcor2)=1.0
      enddo
 672  continue
      close(1)

c- get the completeness curves for the sn cut
      cc0=0.
      cc1=0.
      open(unit=1,file='sn_all.use',status='old')
      read(1,*) x1,xhc0
      do i=1,100
         read(1,*,end=555) x1,x2
         xic(i)=x1
         xcc(i,1)=x2
c         xcc(i,1)=x2/2.
      enddo
 555  continue
      close(1)

c- get the histograms of fout to fin, in both flux and s/n bins
      xf1sa(1)=1.4
      xf1sa(2)=1.8
      xf1sa(3)=3.0
      xf1sa(4)=10.0
      xf1sa(1)=1.4
      xf1sa(2)=2.1
      xf1sa(3)=5.0
      xf1sa(4)=10.0
      xif1s(1)=1.
      xif1s(2)=2.
      xif1s(3)=3.
      xif1s(4)=4.
c      fachist=1.7
      open(unit=2,file='hist_frac',status='old')
      read(2,*) fachist
      close(2)
      open(unit=2,file='histlist',status='old')
      nall=0
      do iall=1,100
         read(2,*,end=667) file1
         open(unit=1,file=file1,status='old')
         nall=nall+1
         nf=0
         do j=1,100
            read(1,*,end=666) i1,x2,x3,x4,i5,x6
            nf=nf+1
            x4o=x4
            sn(nall)=x6
            fout(nall,nf)=x4
            nh=i5
            do i=1,nh
               read(1,*) x1,x2,x3
c               read(1,*) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
               fin(nall,nf,i)=x1/fachist
c               fin(nall,nf,i)=x1
c               hin(nall,nf,i)=x3 
               hin1(nall,nf,i)=x3
               hin2(nall,nf,i)=x3
               hin3(nall,nf,i)=x3
               hin4(nall,nf,i)=x3
c               hin2(nall,nf,i)=x6
c               hin3(nall,nf,i)=x9
c               hin4(nall,nf,i)=x12
               if(i1.lt.5.and.nf.gt.1) then
                  hin1(nall,nf,i)=hin1(nall,nf-1,i)
               endif
            enddo
c- this is for no correction
            if(inocor.eq.1) then
               do i=1,nh
c                  hin(nall,nf,i)=0.
                  hin1(nall,nf,i)=0.
                  hin2(nall,nf,i)=0.
                  hin3(nall,nf,i)=0.
                  hin4(nall,nf,i)=0.
c                  if(fin(nall,nf,i).eq.x4) hin(nall,nf,i)=1.
                  if(fin(nall,nf,i).eq.x4) hin1(nall,nf,i)=1.
                  if(fin(nall,nf,i).eq.x4) hin2(nall,nf,i)=1.
                  if(fin(nall,nf,i).eq.x4) hin3(nall,nf,i)=1.
                  if(fin(nall,nf,i).eq.x4) hin4(nall,nf,i)=1.
               enddo
            endif
         enddo
 666     continue
         close(1)
      enddo
 667  continue
      close(2)

c- define the L bins
      xlmin=42.2
      xlmax=44.3
      nlbin=21
      do i=1,nlbin
         xlhist(i)=xlmin+(xlmax-xlmin)*float(i-1)/float(nlbin-1)
         dflhist(i)=0.d0
         df5hist(i)=0.d0
         dfnc(i)=0.
      enddo

c- read in the data
      xin(1)=4.8
      xin(2)=5.0
      xin(3)=5.5
      xin(4)=6.0
      xin(5)=6.5
      xin(6)=7.0
      fluxc=0.9
      open(unit=1,file='in',status='old')
      n=0
      ntot=0
      ncut=0
      fcut=15000.
      apcut=0.45
c      open(unit=12,file="outdf",status='unknown')
      do i=1,ncmax
         read(1,*,end=669) x1,x2,x3,x4,x5,i6,i7,x8,x9,x10,x11,x12,x13
         ntot=ntot+1
         call xlinint(x11,nfcor1,wave1,fcor1,fcor01)
         x8=1.
         f50=x8*fcor01*sncut
         x12orig=x12
         x12=min(x12,12.)
         rfit=sqrt(x12/2.3)*
     $        (1.64-(sncut-4.8)/(9.5-4.8))/(1.+0.07*(12.-x12))
         rfit=max(1.,rfit)
         f50=f50*rfit
         rat0=x1/x13
         if(x3.ge.sncut.and.x5.le.apcut) ncut=ncut+1
         if(x3.ge.sncut.and.f50.le.fcut.and.x5.gt.apcut
     $        .and.rat0.lt.2..and.x12.lt.12.) then
            n=n+1
            x5=max(x5,0.6)
c            xfc(n)=x13*x5
            xfc(n)=x13
            xfce(n)=max(x2*x5,x1*x5*0.20)
            xsc(n)=x3
            xzc(n)=x4
            x5c(n)=f50
c            ratflux(n)=x1/(x13*x5)
            ratflux(n)=1.
            iddet(n)=i7
            xf1s(n)=x8
         endif
      enddo
 669  continue
      close(1)
c      close(12)

c- add each individual point to the L bins
      open(unit=13,file='out0',status='unknown')
      do itot=1,n
         fo=xfc(itot)
         sno=xsc(itot)
         z=xzc(itot)
         ratf=ratflux(itot)
         foe=1.0*xfce(itot)

c         x5d=x5c(itot)*0.6
c         x5d=x5c(itot)*0.7
         x5d=x5c(itot)

         wave=1215.67*(1.+z)
         call gethin(xf1s(itot),xf1sa,xif1s,hin1,hin2,hin3,hin4,hin)
         call getcc(x5d,wave,xwc,xhc0,xic,xcc,xccu)
         call getfin(fo,foe,ratf,sno,z,nall,nf,nh,sn,fout,fin,hin,
     $        nlbin,xlhist,dflhist,df5hist,dfnc,nz,zl,xl,ncon,sncon,
     $        rcon,sum,inocor,xic,xccu,x5d,cc0,cc1)
         call xlinint(z,nz,zl,xl,xconv)
         xlum=40.+log10(fo*1.e-17*xconv*ratf)
         write(13,*) z,fo,xlum,sno,iddet(itot)
      enddo
      close(13)

c- skip this normalization
      dsum=0.
      do i=1,nlbin
         dsum=dsum+dflhist(i)
      enddo
      if(n.gt.0.and.ntot.gt.0) print *,n,ncut,ntot,
     $     float(ncut)/float(n),float(ncut)/float(ntot)
      fac=float(ntot)/float(n)
      fac=1.

      open(unit=11,file='out',status='unknown')
      ratmin=1e10
      do i=1,nlbin-1
         xbcen=(xlhist(i)+xlhist(i+1))/2.
         rat=sngl(dflhist(i)/df5hist(i))
         ratmin=min(ratmin,rat)
         xout1=fac*sngl(dflhist(i))
         xout2=fac*sngl(df5hist(i))
         if(xbcen.gt.43.3.and.rat.gt.ratmin) then
c            xout1=xout2*ratmin
         endif
         write(11,*) xbcen,xout1,xout2,sngl(dfnc(i))
      enddo
      close(11)

      end

      subroutine getfin(fo,foe,ratf,sno,z,nall,nf,nh,sn,fout,fin,
     $     hin,nlbin,xlhist,dflhist,df5hist,dfnc,nz,zl,xl,ncon,sncon,
     $     rcon,sum,inocor,xic,xccu,x5d,cc0,cc1)
      parameter(nstot=10,nftot=10,nmax=100)
      real fin(nstot,nftot,nmax),hin(nstot,nftot,nmax),sn(nstot)
      real xlhist(nmax),fout(nstot,nftot),sncon(ncon),rcon(ncon)
      real xi(nmax),xin(nmax),yin(nmax),zl(nz),xl(nz)
      real fuse(nmax),foina(3),xina(nmax),xic(200),xccu(200),xccor(200)
      real yin0(nmax),xinanc(nmax),yinnc(nmax)
      real*8 dflhist(nmax),df5hist(nmax),dfnc(nmax)

c      ncc=199
      ncc=100

      do i=1,nall
         xi(i)=float(i)
      enddo

c- find the nth element for the s/n histogram
      call xlinint(sno,nall,sn,xi,xi0)
      is=nint(xi0)
c- find the flux to L conversion
      call xlinint(z,nz,zl,xl,xconv)
c- find the confirmed rate
      call xlinint(sno,ncon,sncon,rcon,rcon0)
      do i=1,nf
         xi(i)=float(i)
         xin(i)=fout(is,i)
      enddo

c- explore at +-1 sigma of the input flux
      foina(1)=fo
      foina(2)=fo-foe
      foina(3)=fo+foe
c      isigmax=3
      isigmax=1
c      rcon0=rcon0/3.
      if(inocor.eq.1) then
         rcon0=1.
         isigmax=1
      endif
      do isig=1,isigmax

c- find the proper histogram
         foin=foina(isig)
         call xlinintf(foin,nf,xin,xi,i1,i2,frac)

c- if outside the bounds, scale from the last s/n histogram
         facf=1.0
c         if(is.eq.nall) then
         if(is.eq.nall.or.i1.eq.i2) then
            facf=fo/fout(is,i1)
            do i=1,nh
               fuse(i)=fin(is,i1,i)*facf
            enddo
         else
            do i=1,nh
               fuse(i)=fin(is,i1,i)*facf
            enddo
         endif

         do i=1,nh
c            call xlinint(fuse(i),ncc,xic,xccu,xccor(i))
            call xnear(fuse(i),ncc,xic,xccu,xccor(i))
         enddo

c- find the L for each of the histogram flux bins
         sum=0.
         do i=1,nh
            h1=hin(is,i1,i)
            h2=hin(is,i2,i)
            if(frac.eq.0.) then
               yin(i)=h1
            else
               yin(i)=h1+(h2-h1)*frac
            endif
c            xina(i)=40.+log10(fuse(i)*1.e-17*xconv)
            xina(i)=40.+log10(fuse(i)*1.e-17*xconv*ratf)

c- get the completeness array
            yin0(i)=yin(i)
            if(xccor(i).gt.cc1) then
               if(xccor(i).gt.0.) yin(i)=yin(i)/xccor(i)
            else
               yin(i)=yin(i)/cc1
            endif
            if(xccor(i).le.cc0) yin(i)=0.
            if(xccor(i).eq.0.) yin(i)=0.
c            if(fuse(i).lt.13.) yin(i)=0.
c            print *,i,fo,is,i1,fuse(i),yin(i),xccor(i),fout(is,i1)
c            if(fo/fuse(i).gt.2.) yin(i)=0.
c            if(xccor(i).lt.0.001.and.yin(i).gt.0.) 
c     $           print *,xccor(i),fuse(i),fo,yin(i),x5d
            sum=sum+yin(i)
         enddo

c         if(inocor.eq.1) then
            xinanc(1)=40.+log10(fo*1.e-17*xconv)
            yinnc(1)=1.
            do i=2,nh
               xinanc(i)=0.
               yinnc(i)=0.
            enddo
c         endif

c- add into the L binning
         do i=1,nh
            do j=1,nlbin-1
               if(xina(i).ge.xlhist(j).and.xina(i).lt.xlhist(j+1)) then
                  dflhist(j)=dflhist(j)+dble(yin(i)*rcon0)
                  df5hist(j)=df5hist(j)+dble(yin0(i))
               endif
               if(xinanc(i).ge.xlhist(j).and.
     $              xinanc(i).lt.xlhist(j+1)) then
                  dfnc(j)=dfnc(j)+dble(yinnc(i))
               endif
            enddo
         enddo
      enddo

      return
      end

      subroutine getcc(x5d,wave,xwc,xhc0,xic,xcc,xccu)
      parameter(nmax=200)
      real xwc(nmax),xic(nmax),xcc(200,100),xccu(nmax)
      real xin(nmax),yin(nmax)

      ncc=100
c- scale to the 50% value
      scale=x5d/xhc0
c      print *,scale,x5d
      do i=1,ncc
         xin(i)=scale*xic(i)
         yin(i)=xcc(i,1)
      enddo
      do i=1,ncc
c         call xlinint0(xic(i),ncc,xin,yin,y0)
         call xnear0(xic(i),ncc,xin,yin,y0)
         xccu(i)=y0
      enddo

      return
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
      if(xp.ge.x(n)) yp=y(n)
      return
      end

      subroutine xnear(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
c            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            yp=y(j)
            if(abs(xp-x(j+1)).le.abs(xp-x(j))) yp=y(j+1)
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end

      subroutine xlinint0(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            if(y(j).eq.0.) yp=0.
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end

      subroutine xnear0(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
c            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            yp=y(j)
            if(abs(xp-x(j+1)).le.abs(xp-x(j))) yp=y(j+1)
            if(y(j).eq.0.) yp=0.
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end

      subroutine xlinintf(xp,n,x,y,i1,i2,frac)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            i1=j
            i2=j+1
            frac=(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      frac=0.
      if(xp.lt.x(1)) then
         i1=1
         i2=1
      else
         i1=n
         i2=n
      endif
      return
      end

      subroutine gethin(xf1s,xf1sa,xif1s,hin1,hin2,hin3,hin4,hin)
      parameter(nstot=10,nftot=10,nmax=100)
      real hin(nstot,nftot,nmax)
      real hin1(nstot,nftot,nmax),hin2(nstot,nftot,nmax)
      real hin3(nstot,nftot,nmax),hin4(nstot,nftot,nmax)
      real xf1sa(100),xif1s(100)
      call xlinint(xf1s,4,xf1sa,xif1s,xiout)
      x1=0.
      x2=0.
      x3=0.
      x4=0.
      if(xiout.le.1) x1=1.
      if(xiout.ge.4) x4=1.
      if(xiout.gt.1.and.xiout.le.2) then
         x1=2.-xiout
         x2=1.-x1
      endif
      if(xiout.gt.2.and.xiout.le.3) then
         x2=3.-xiout
         x3=1.-x2
      endif
      if(xiout.gt.3.and.xiout.lt.4) then
         x3=4.-xiout
         x4=1.-x3
      endif

c      print *,xf1s,x1,x2,x3,x4
      do i=1,nstot
         do j=1,nftot
            do k=1,nmax
               hin(i,j,k)=x1*hin1(i,j,k)+x2*hin2(i,j,k)+
     $              x3*hin3(i,j,k)+x4*hin4(i,j,k)
            enddo
         enddo
      enddo

      return
      end
