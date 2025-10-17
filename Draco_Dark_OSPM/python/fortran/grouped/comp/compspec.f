
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),ye(nmax),yec(nmax)
      real xc(nmax),yc(nmax),xin(nmax),yfita(100000,7),xerra(nmax)
      real xb(nmax),yb(nmax),xfit(10),yfit(10),xnorma(nmax)
      character a1*12,a2*12,an1*6,an2*6,file1*80,file2*80

      wmin=4000.
      wmax=5000.
      nf=18
      ibin=13
      ib1=(ibin-1)/2
      xib=float(ibin)

c      nfit=5
c      xfit(1)=3700.
c      xfit(2)=4100.
c      xfit(3)=4500.
c      xfit(4)=4900.
c      xfit(5)=5300.
      nfit=7
      xfit(1)=3600.
      xfit(2)=3900.
      xfit(3)=4200.
      xfit(4)=4500.
      xfit(5)=4800.
      xfit(6)=5100.
      xfit(7)=5400.

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=3500.
      xmax=5450.

      open(unit=1,file='compspec.in',status='old')

      ymtot=3e5
      ymmin=1e3
c      ymtot=3e4
c      ymmin=1e2

      ntot=0
      nl=0
      ic=0
      do il=1,10000
         read(1,*,end=666) i1,i2,x3,xf1,xf2,a1,a2
         if(xf1.gt.ymtot.or.xf2.gt.ymtot) goto 866
         if(xf1.lt.ymmin.or.xf2.lt.ymmin) goto 866
         call getaname(i1,an1,in1)
         call getaname(i2,an2,in2)

c- first file
         file1="spec/"//a1//"_"//an1(1:in1)//".spec"
         file2="spec/"//a1//"_"//an1(1:in1)//".list"
         open(unit=2,file=file1,status='old')
         ymin=1e10
         ymax=-ymin
         n=0
         nin=0
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9
            if(x2.gt.0) then
               n=n+1
               x(n)=x1
               if(x4.ne.0) then
                  xap=x6*x2/x4
               else
                  xap=0
               endif
               y(n)=x2/x9
               y(n)=x2
c               y(n)=xap
               ye(n)=x3/x9
               ymin=min(ymin,y(n))
               ymax=max(ymax,y(n))
c               if(x(n).gt.4000.and.x(n).lt.5000.and.x(n).ne.0) then
               if(x(n).gt.wmin.and.x(n).lt.wmax.and.x(n).ne.0) then
                  nin=nin+1
                  xin(nin)=x9
               endif                  
            endif
         enddo
 667     continue
         close(2)
         call biwgt(xin,nin,xb1,xs1)

c- second file
         file1="spec/"//a2//"_"//an2(1:in2)//".spec"
         file2="spec/"//a2//"_"//an2(1:in2)//".list"
         open(unit=2,file=file1,status='old')
         nc=0
         nin=0
         do i=1,nmax
            read(2,*,end=668) x1,x2,x3,x4,x5,x6,x7,x8,x9
            if(x2.gt.0) then
               nc=nc+1
               xc(nc)=x1
               if(x4.ne.0) then
                  xap=x6*x2/x4
               else
                  xap=0
               endif
               yc(nc)=x2/x9
               yc(nc)=x2
c               yc(nc)=xap
               yec(nc)=x3/x9
               ymin=min(ymin,yc(nc))
               ymax=max(ymax,yc(nc))
c               if(xc(nc).gt.4000.and.xc(nc).lt.5000.and.
c     $              xc(nc).ne.0) then
               if(xc(nc).gt.wmin.and.xc(nc).lt.wmax.and.
     $              xc(nc).ne.0) then
                  nin=nin+1
                  xin(nin)=x9
               endif                  
            endif
         enddo
 668     continue
         close(2)
         call biwgt(xin,nin,xb2,xs2)
         if(xb1.lt.0.6.or.xb2.lt.0.6) goto 866

         call pgenv(3500.,5500.,ymin,ymax,0,0)
         call pglabel("Wavelength","Relative Flux","")
         call pgsls(1)
         call pgslw(2)
         call pgsch(1.8)
         call pgsci(1)
         call pgline(n,x,y)
         call pgsci(2)
         call pgline(nc,xc,yc)
         call pgsci(1)
         call getoffset(n,x,y,ye,nc,xc,yc,yec,xfit,yfit,xnorm,xerr)
         do j=1,nfit
            xin(j)=yfit(j)
         enddo
         call biwgt(xin,5,xbout,xsout)
         if(xbout.lt.0.75.or.xbout.gt.1.4) goto 866
         if(xsout.gt.0.1) goto 866
c         print *,i1,i2,xf1,xf2,xbout,xsout
         ntot=ntot+1
c         print *,ntot,xbout,xsout
         xnorma(ntot)=xnorm
         xerra(ntot)=xerr
         do j=1,nfit
            yfita(ntot,j)=yfit(j)
         enddo
         goto 867
 866     continue
c         print *,i1,i2,xf1,xf2,xbout,xsout
 867     continue
      enddo

 666  continue
      close(1)

      call biwgt(xnorma,ntot,xbouta,xsouta)
      call biwgt(xerra,ntot,xbea,xsea)
      do j=1,nfit
         do i=1,ntot
            xin(i)=yfita(i,j)
         enddo
         call biwgt(xin,ntot,xbout,xsout)
         write(*,1001) xfit(j),xbout,xsout,xsouta,xbea,xsea,ntot
      enddo

      call pgend
 1001 format(f6.1,5(1x,f7.3),1x,i6)

      end

      subroutine getoffset(n,x,y,ye,nc,xc,yc,yec,xfit,yfit,xnorm,xerr)
      real x(n),y(n),xc(nc),yc(nc),ye(n),yec(n)
      real xfit(10),yfit(10),xin(1000),xin2(1000)

c      nfit=5
c      xhalf=200.
      nfit=7
      xhalf=100.

c- first get overall
      ws=3800.
      we=5200.
      nin=0
      sum1=0.
      den1=0.
      sum2=0.
      den2=0.
      do i1=1,n
         if(x(i1).gt.ws.and.x(i1).lt.we) then
            if(y(i1).gt.0..and.ye(i1).gt.0) then
               nin=nin+1
c               xin(nin)=y(i1)
c               xin2(nin)=ye(i1)/y(i1)
               xin(nin)=y(i1)
               xin2(nin)=ye(i1)/y(i1)
               sum1=sum1+y(i1)/ye(i1)/ye(i1)
               den1=den1+1./ye(i1)/ye(i1)
            endif
         endif
      enddo
      call biwgt(xin,nin,xb1,xs)
      call biwgt(xin2,nin,xb1e,xs1e)
      nin=0
      do i1=1,nc
         if(xc(i1).gt.ws.and.xc(i1).lt.we) then
            if(yc(i1).gt.0..and.yec(i1).gt.0) then
               nin=nin+1
               xin(nin)=yc(i1)
               xin2(nin)=yec(i1)/yc(i1)
               sum2=sum2+yc(i1)/yec(i1)/yec(i1)
               den2=den2+1./yec(i1)/yec(i1)
            endif
         endif
      enddo
      call biwgt(xin,nin,xb2,xs)
      call biwgt(xin2,nin,xb2e,xs2e)
      xnorm=xb1/xb2
      xnorm=1.
      xerr=xb1e/xb2e
      sum1=sum1/den1
      sum2=sum2/den2
      xerr=sum1/sum2

c- now versus wavelength
      do i=1,nfit
         ws=xfit(i)-xhalf
         we=xfit(i)+xhalf
         nin=0
         do i1=1,n
            if(x(i1).gt.ws.and.x(i1).lt.we) then
c               call xlinint(x(i1),nc,xc,yc,yv)
c               if(y(i1).gt.0..and.yv.gt.0.) then
               if(y(i1).gt.0.) then
                  nin=nin+1
c                  xin(nin)=y(i1)/y
                  xin(nin)=y(i1)
               endif
            endif
         enddo
         call biwgt(xin,nin,xb1,xs)
         nin=0
         do i1=1,nc
            if(xc(i1).gt.ws.and.xc(i1).lt.we) then
               if(yc(i1).gt.0.) then
                  nin=nin+1
                  xin(nin)=yc(i1)
               endif
            endif
         enddo
         call biwgt(xin,nin,xb2,xs)
c         yfit(i)=xb
         yfit(i)=xb1/xb2/xnorm
      enddo      

      return
      end

      subroutine getaname(i1,an1,in1)
      character an1*6
      if(i1.lt.10) then
         write(an1(1:1),1001) i1
         in1=1
      endif
      if(i1.ge.10.and.i1.lt.100) then
         write(an1(1:2),1002) i1
         in1=2
      endif
      if(i1.ge.100.and.i1.lt.1000) then
         write(an1(1:3),1003) i1
         in1=3
      endif
      if(i1.ge.1000.and.i1.lt.10000) then
         write(an1(1:4),1004) i1
         in1=4
      endif
      if(i1.ge.10000.and.i1.lt.100000) then
         write(an1(1:5),1005) i1
         in1=5
      endif
      if(i1.ge.100000.and.i1.lt.1000000) then
         write(an1(1:6),1006) i1
         in1=6
      endif
 1001 format(i1)
 1002 format(i2)
 1003 format(i3)
 1004 format(i4)
 1005 format(i5)
 1006 format(i6)
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
      if(xp.gt.x(n)) yp=y(n)
      return
      end
