
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),xn2(nmax),yn2(nmax)
      real xb(nmax),yb(nmax),xs(nmax),ys(nmax),wc(nmax),rat(nmax)
      real wl(nmax),wh(nmax),xwe(nmax),xecor(nmax),rata(100,nmax)
      real xin(nmax),xwp(nmax),xpcor(nmax)
      character file1*80,file2*80,c1*15,nullstr*1
      logical simple,extend,anyf

      ipall=1

      ibin1=7
      ibin2=17
      ibin1=2
      ibin2=6

      nwb=11
      nw=nwb-1
      wmin=3850.
      wmax=5500.
      do i=1,nw
         wl(i)=wmin+(wmax-wmin)*float(i-1)/float(nw)
         wh(i)=wmin+(wmax-wmin)*float(i)/float(nw)
         wc(i)=(wl(i)+wh(i))/2.
      enddo

      open(unit=1,file="extcor.dat",status='old')
      ne=0
      do i=1,nmax
         read(1,*,end=671) x1,x2
         ne=ne+1
         xwe(ne)=x1
         xecor(ne)=x2
      enddo
 671  continue
      close(1)

      open(unit=1,file="postcor.dat",status='old')
      np=0
      do i=1,nmax
         read(1,*,end=672) x1,x2
         np=np+1
         xwp(np)=x1
         xpcor(np)=x2
      enddo
 672  continue
      close(1)

      if(ipall.eq.1) call pgbegin(0,'?',3,3)
c      if(ipall.eq.1) call pgbegin(0,'?',1,1)
      if(ipall.eq.0) call pgbegin(0,'?',1,1)
c      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      if(ipall.eq.0) then
         call pgsci(1)
         call pgenv(3800.,5500.,0.5,1.5)
      endif

      xmin=3500.
      xmax=5500.

      open(unit=1,file='listcomp',status='old')
      open(unit=11,file='out',status='unknown')

      nl=0
      ic=0
      do il=1,20000
         read(1,*,end=666) file1,file2
         open(unit=2,file=file1,status='old',err=888)
         goto 889
 888     continue
c         print *,"Does not exist: ",file1
         goto 866
 889     continue
         ymin=1e10
         ymax=-ymin
         n=0
         x10=0.
         sumd=0.
         nd=0
         do i=1,nmax
c            read(2,*,end=667) x1,x2
            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9
            if(x2.ne.0) then
               n=n+1
               x(n)=x1
               call xlinint(x(n),ne,xwe,xecor,y0)
               y(n)=x2/2./x9/y0
               call xlinint(x(n),np,xwp,xpcor,y0)
               y(n)=y(n)/y0
               sumd=sumd+y(n)
               nd=nd+1
            endif
         enddo
 667     continue
         close(2)
         if(nd.gt.0) sumd=sumd/float(nd)
         ic=2
         if(x9.lt.0.6) ic=3
         if(sumd.lt.1) goto 866

c         if(sumd.lt.50) goto 866

         call smbin(ibin1,n,x,y,nbb,xn,yn,ymin,ymax)

         im1=51
         ier=0
         iread=0
c         call ftgiou(im1,ier)
c         print *,"ier= ",ier
         call ftopen(im1,file2,iread,iblock,ier)
         if(ier.gt.0) goto 866
         call ftmahd(im1,2,ihd,ier)

         nsa=3841
         do i=1,nsa
            call ftgcve(im1,1,i,1,1,0.,flux,anyf,ier)
            call ftgcve(im1,2,i,1,1,0.,wave,anyf,ier)
            if(10**(wave).gt.5520.) goto 333
            ns=i
            xs(i)=10**(wave)
            ys(i)=flux
         enddo
 333     continue

         call smbin(ibin2,ns,xs,ys,nbb2,xn2,yn2,ymin,ymax)
         call ftclos(im1,ier)

         call getrat(nbb,xn,yn,nbb2,xn2,yn2,nw,wc,wl,wh,rat,voff)

         do i=1,nw
            xin(i)=rat(i)
         enddo
         call biwgt(xin,nw,xb0,xs0)
         do i=1,nw
c            rat(i)=rat(i)/xb0
         enddo
         ymin=1e10
         ymax=-ymin
         do i=1,nbb
c            yn(i)=yn(i)/xb0
c            call xlinint(xn(i),nw,wc,rat,y0)
c            yn(i)=yn(i)/y0
            ymin=min(ymin,yn(i))
            ymax=max(ymax,yn(i))
         enddo
         do i=1,nbb2
            ymin=min(ymin,yn2(i))
            ymax=max(ymax,yn2(i))
         enddo

         call pgsls(1)
         call pgslw(2)
         call pgsci(1)
         call pgsch(1.8)
         ybit=(ymax-ymin)/10.
         ymin=ymin-ybit
         ymax=ymax+ybit

         c1=file2(12:26)
         print *,il,nl,xb0,xs0,voff," ",c1
         write(11,*) il,nl,xb0,xs0,voff," ",c1
         xblo=0.77
         xbhi=1.3
c         xblo=0.5
c         xbhi=1.5
         if(xb0.gt.xblo.and.xb0.lt.xbhi.and.xs0.lt.0.35
     $        .and.ic.ne.3) then
            ic=1
            nl=nl+1
            do i=1,nw
               rata(i,nl)=rat(i)
            enddo
            if(ipall.eq.0) call pgline(nw,wc,rat)
         else
            if(ic.ne.3) ic=4
         endif
         if(ipall.eq.1) then
            call pgsch(1.2)
            call pgenv(xmin,xmax,ymin,ymax,0,0)
            call pgsci(ic)
            call pgline(nbb,xn,yn)
            call pgsch(1.8)
            call pgsci(1)
c            call pglabel('Wavelength (\(2078))',
c     $           '10\U-17\D ergs/cm\U2\D/s/\(2078)','')
            call pgmtxt('B',2.0,0.5,0.5,'Wavelength (\(2078))')
            call pgmtxt('L',1.5,0.5,0.5,
     $           '10\U-17\D ergs/cm\U2\D/s/\(2078)')
            call pgsch(1.9)
            call pgmtxt('T',-1.5,0.5,0.5,c1)
            call pgsch(1.6)
            call pgsci(2)
            call pgline(nbb2,xn2,yn2)
            call pgsci(1)
         endif

 866     continue
         close(2)
      enddo
 666  continue
      close(1)
      close(11)

      open(unit=11,file='out2',status='unknown')
      do i=1,nw
         do j=1,nl
            xin(j)=rata(i,j)
         enddo
         call biwgt(xin,nl,xbr,xsr)
         write(11,1101) wc(i),xbr,xsr,nl
         write(*,1101) wc(i),xbr,xsr,nl
      enddo
      close(11)

      call pgend

 1101 format(1x,f8.2,2(1x,f6.3),1x,i5)
      end
      
      subroutine getrat(nbb,xn,yn,nbb2,xn2,yn2,nw,wc,wl,wh,rat,voff)
      real xn(nbb),yn(nbb),xn2(nbb2),yn2(nbb2),wc(nw),rat(nw)
      real wl(nw),wh(nw),xin1(1000),xin2(1000),yin1(1000),yin2(1000)
      real rmsa(1000),ync(10000),yn2c(10000)

      call fitcont(nbb,xn,yn,ync)
      call fitcont(nbb2,xn2,yn2,yn2c)
      call getveloff(nbb,xn,ync,nbb2,xn2,yn2c,voff)

      do i=1,nw
         sum1=0.
         n1=0
         do j=1,nbb
            if(xn(j).ge.wl(i).and.xn(j).lt.wh(i)) then
               sum1=sum1+yn(j)
               n1=n1+1
               xin1(n1)=xn(j)
               yin1(n1)=yn(j)
            endif
         enddo
         sum1=sum1/float(n1)
         sum2=0.
         n2=0
         do j=1,nbb2
            if(xn2(j).ge.wl(i).and.xn2(j).lt.wh(i)) then
               sum2=sum2+yn2(j)
               n2=n2+1
               xin2(n2)=xn2(j)
               yin2(n2)=yn2(j)
            endif
         enddo
         sum2=sum2/float(n2)

         rat(i)=sum1/sum2
         
         ktry=100
         rstep=0.005
         rmin=rat(i)-rstep*float(ktry)/2.
         rmax=rat(i)+rstep*float(ktry)/2.
         rmsl=1.e10
         do k=1,ktry
            rtry=rmin+(rmax-rmin)*float(k-1)/float(ktry-1)
            do j=1,n1
               call xlinint(xin1(j),nbb2,xn2,yn2,y0)
               ynew=yin1(j)/rtry
               rmsa(j)=(y0-ynew)**2
            enddo
            call biwgt(rmsa,n1,xbrms,xsrms)
            if(xbrms.lt.rmsl) then
               rmsl=xbrms
               rbest=rtry
            endif
         enddo
c         print *,i,rat(i),rbest
         rat(i)=rbest
      enddo

      return
      end

      subroutine smbin(ibin,n,x,y,nbb,xn,yn,ymin,ymax)
      real x(n),y(n),xb(10000),yb(10000),xn(n),yn(n)

      ib1=(ibin-1)/2
      nbb=0
      do j=1,n,ibin
         nbb=nbb+1
         istart=max(0,j-ib1)
         iend=istart+ibin-1
         if(iend.gt.n) then
            iend=n
            istart=n-ibin+1
         endif
         sum=0.
         nb=0
         do is=istart,iend
            sum=sum+y(is)
            nb=nb+1
            yb(nb)=y(is)
            xb(nb)=x(is)
         enddo
         call biwgt(yb,nb,xbb,xsb)
c         yn(nbb)=xbb
         yn(nbb)=sum/float(nb)
         call biwgt(xb,nb,xbb,xsb)
         xn(nbb)=xbb
         ymin=min(ymin,yn(nbb))
         ymax=max(ymax,yn(nbb))
      enddo

      return
      end

      subroutine getveloff(n,x,y,n2,x2,y2,voff)
      real x(n),y(n),x2(n),y2(n)

      xmin=3900.
      xmax=5500.

      vmin=-50.
      vmax=150.
      nvel=100

      rmin=1e10
      do ivel=1,nvel
         v=vmin+float(ivel-1)*(vmax-vmin)/float(nvel-1)
         rms=0.
         do i=1,n
            if(x(i).gt.xmin.and.x(i).lt.xmax) then
               woff=v*x(i)/2.99e5
               call xlinint(x(i)+woff,n2,x2,y2,y0)
               rms=rms+(y0-y(i))**2
            endif
         enddo
         if(rms.lt.rmin) then
            rmin=rms
            voff=v
         endif
      enddo

      return
      end

      subroutine fitcont(n,x,y,yc)
      real x(n),y(n),yc(10000),yin(10000),yin2(10000)
      real xcont(1000),ycont(1000)

      xcut=0.2
      xmin=3900.
      xmax=5500.
      nc=8
      xs=(xmax-xmin)/float(nc)
      do ic=1,nc
         xlo=xmin+float(ic-1)*xs
         xhi=xmin+float(ic)*xs
         nin=0
         do i=1,n
            if(x(i).ge.xlo.and.x(i).lt.xhi) then
               nin=nin+1
               yin(nin)=y(i)
            endif
         enddo
         if(nin.eq.0) return
         call sort(nin,yin)
         ns=nint(xcut*float(nin))
         nin2=0
         do i=nin,nin-ns,-1
            nin2=nin2+1
            yin2(nin2)=yin(i)
         enddo
         call biwgt(yin2,nin2,xb1,xs1)
         xcont(ic)=(xhi+xlo)/2.
         ycont(ic)=xb1
      enddo

      do i=1,n
         call xlinint(x(i),nc,xcont,ycont,y0)
         yc(i)=y(i)/y0
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
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
