Marked
      parameter(nstot=10,nftot=10,nmax=100,ncmax=1000000)
      real fin(nstot,nftot,nmax),hin(nstot,nftot,nmax),sn(nstot)
      real xlhist(nmax),flhist(nmax),fout(nstot,nftot),finm(nstot,nftot)
      real zl(10000),xl(10000),xfc(ncmax),xsc(ncmax),xzc(ncmax)
      real sncon(100),rcon(100),xfce(ncmax),x5c(ncmax),yin(10),xin(10)
      real xwc(200),xic(200),xhc(200),xcc(200,100),xccu(200)
      real*8 dx2,dflhist(nmax),df5hist(nmax),dfnc(nmax)
      integer*8 i6,i7
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
         rcon(ncon)=x2
         if(inocor.eq.1) rcon(ncon)=1.
      enddo
 670  continue
      close(1)

c- get the completeness curves for the sn cut
      open(unit=1,file='snlist',status='old')
      read(1,*)
      smin=1e10
      do i=1,1000
         read(1,*,end=777) file1,x1,x2,x3
         if(abs(x1-sncut).lt.smin) then
            smin=abs(x1-sncut)
            imin=i
         endif
      enddo
 777  continue
      rewind(1)
      read(1,*)
      do i=1,imin
         read(1,*) file1,x1,cc0,cc1,xf50
      enddo
      close(1)
      open(unit=1,file=file1,status='old')
      read(1,*) x1,xwc(1),xwc(2),xwc(3),xwc(4),xwc(5),xwc(6),xwc(7)
      read(1,*) x1,xhc(1),xhc(2),xhc(3),xhc(4),xhc(5),xhc(6),xhc(7)
      do i=1,199
         read(1,*,end=555) x1,x2,x3,x4,x5,x6,x7,x8
         xic(i)=x1
         xcc(i,1)=x2
         xcc(i,2)=x3
         xcc(i,3)=x4
         xcc(i,4)=x5
         xcc(i,5)=x6
         xcc(i,6)=x7
         xcc(i,7)=x8
      enddo
 555  continue
      close(1)

c- get the histograms of fout to fin, in both flux and s/n bins
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
            sn(nall)=x6
            fout(nall,nf)=x4
            finm(nall,nf)=x2
            nh=i5
            do i=1,nh
               read(1,*) x1,x2,x3
               fin(nall,nf,i)=x1
               hin(nall,nf,i)=x3
            enddo
c- this is for no correction
            if(inocor.eq.1) then
               do i=1,nh
                  hin(nall,nf,i)=0.
                  if(fin(nall,nf,i).eq.x4) hin(nall,nf,i)=1.
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
      xlmax=46.0
      nlbin=21
      do i=1,nlbin
         xlhist(i)=xlmin+(xlmax-xlmin)*float(i-1)/float(nlbin-1)
         dflhist(i)=0.d0
         df5hist(i)=0.d0
      enddo

c- read in the data (truncated it to 5 cols for quicker read)
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
      do i=1,ncmax
c         read(1,*,end=669) x1,x2,x3,x4,x5,i6,i7,x8,x9,x10,x11,x12,x13
         read(1,*,end=669) x1,x2,x3,x4,x5,i6,i7,x8
         ntot=ntot+1
         f50=x8*xf50
c         if(x3.ge.sncut.and.f50.gt.fcut) ncut=ncut+1
         if(x3.ge.sncut.and.x5.le.apcut) ncut=ncut+1
c         if(x3.ge.sncut.and.f50.gt.fcut) print *,i7,f50,x1,x1/f50,x8,x5
         if(x3.ge.sncut.and.f50.le.fcut.and.x5.gt.apcut) then
            n=n+1
c            x5=0.91
c            x5=1.
            x5=max(x5,0.6)
            xfc(n)=x1*x5
            xfce(n)=max(x2*x5,x1*x5*0.20)
            xsc(n)=x3
            xzc(n)=x4
c            yin(1)=x8
c            yin(2)=x9
c            yin(3)=x10
c            yin(4)=x11
c            yin(5)=x12
c            yin(6)=x13
c            call xlinint(sncut,6,xin,yin,y0)
c            x5c(n)=y0
            x5c(n)=f50
         endif
      enddo
 669  continue
      close(1)

c- add each individual point to the L bins
c      do itot=1,1
      do itot=1,n
c      do itot=110,110
         fo=xfc(itot)
         sno=xsc(itot)
         z=xzc(itot)
c         foe=1.3*xfce(itot)
         foe=1.0*xfce(itot)
         x5d=x5c(itot)
         wave=1215.67*(1.+z)
c         fo=100.
c         sno=47.0
c         z=2.5
c         print *,fo,foe,sno,x5d
         call getcc(x5d,wave,xwc,xhc,xic,xcc,xccu)
         call getfin(fo,foe,sno,z,nall,nf,nh,sn,fout,finm,fin,hin,
     $        nlbin,xlhist,dflhist,df5hist,dfnc,nz,zl,xl,ncon,sncon,
     $        rcon,sum,inocor,xic,xccu,x5d,cc0,cc1)
      enddo

c- skip this normalization
      dsum=0.
      do i=1,nlbin
         dsum=dsum+dflhist(i)
      enddo
      print *,n,ncut,ntot,float(ncut)/float(n),float(ncut)/float(ntot)
      fac=float(ntot)/float(n)
      fac=1.

      open(unit=11,file='out',status='unknown')
      do i=1,nlbin-1
         xbcen=(xlhist(i)+xlhist(i+1))/2.
         write(11,*) xbcen,fac*sngl(dflhist(i)),fac*sngl(df5hist(i)),
     $        sngl(dfnc(i))
      enddo
      close(11)

      end

      subroutine getfin(fo,foe,sno,z,nall,nf,nh,sn,fout,finm,fin,hin,
     $     nlbin,xlhist,dflhist,df5hist,dfnc,nz,zl,xl,ncon,sncon,rcon,
     $     sum,inocor,xic,xccu,x5d,cc0,cc1)
      parameter(nstot=10,nftot=10,nmax=100)
      real fin(nstot,nftot,nmax),hin(nstot,nftot,nmax),sn(nstot)
      real xlhist(nmax),fout(nstot,nftot),sncon(ncon),rcon(ncon)
      real xi(nmax),xin(nmax),yin(nmax),zl(nz),xl(nz),finm(nstot,nftot)
      real fuse(nmax),foina(3),xina(nmax),xic(200),xccu(200),xccor(200)
      real yin0(nmax),xinanc(nmax),yinnc(nmax)
      real*8 dflhist(nmax),df5hist(nmax),dfnc(nmax)

      ncc=199

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
c            facf=fo/finm(is,i1)
            facf=fo/fout(is,i1)
            do i=1,nh
               fuse(i)=fin(is,i1,i)*facf
            enddo
         else
            do i=1,nh
               fuse(i)=fin(is,i1,i)*facf
c               print *,fin(is,i1,i),fin(is+1,i1,i),xi0,is,i1
            enddo
         endif

         do i=1,nh
            call xlinint(fuse(i),ncc,xic,xccu,xccor(i))
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
            xina(i)=40.+log10(fuse(i)*1.e-17*xconv)

c- get the completeness array
            yin0(i)=yin(i)
            if(xccor(i).gt.cc1) then
               yin(i)=yin(i)/xccor(i)
            else
               yin(i)=yin(i)/cc1
            endif
            if(xccor(i).lt.cc0) yin(i)=0.
c            if(fuse(i).lt.13.) yin(i)=0.
c            print *,i,fo,is,i1,fuse(i),yin(i),xccor(i),fout(is,i1)
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

      subroutine getcc(x5d,wave,xwc,xhc,xic,xcc,xccu)
      parameter(nmax=200)
      real xwc(nmax),xic(nmax),xhc(nmax),xcc(200,100),xccu(nmax)
      real xin(nmax),yin(nmax)

      ncc=199
c- get closest wavelength array
      xmin=1e10
      do i=1,7
         if(abs(xwc(i)-wave).lt.xmin) then
            xmin=abs(xwc(i)-wave)
            iwave=i
         endif
      enddo
c- scale to the 50% value
      scale=x5d/xhc(iwave)
c      print *,scale,x5d
      do i=1,ncc
         xin(i)=scale*xic(i)
         yin(i)=xcc(i,iwave)
      enddo
      do i=1,ncc
         call xlinint(xic(i),ncc,xin,yin,y0)
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
