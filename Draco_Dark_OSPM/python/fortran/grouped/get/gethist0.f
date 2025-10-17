Marked
      parameter(nmax=20000000)
      real xi(nmax),xo(nmax),xw(nmax),w(nmax),xiw(nmax),xh0(1000)
      real xin(nmax),xinw(nmax),xbl(50),xbh(50),xh1(1000),xh2(1000)
      real x1sig(nmax),xin2(nmax),xhnew(nmax)
      integer is(nmax)

      fcutin=2.5

      wlo=3500.
      whi=4000.
      wlo=3500.
      whi=5600.

      ncut=10
      ncut=0

      goto 667
      open(unit=1,file='weightn.use',status='old')
      nw=0
      do i=1,nmax
         read(1,*,end=667) x1,x2
         nw=nw+1
         xw(nw)=x1
         w(nw)=x2
c         w(nw)=1.
      enddo
 667  continue
      close(1)

c - assumes input, output
      open(unit=1,file='t1',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,i3,x4,sn,x6
         if(x4.gt.wlo.and.x4.lt.whi.and.x1.ge.fcutin) then
            rat=x1/x2
            if(rat.lt.10) then
               n=n+1
c               xi(n)=x1-5.
               xi(n)=x1
               xo(n)=x2
               is(n)=i3
c     call xlinint(x1,nw,xw,w,w0)
               w0=1.
               xiw(n)=w0
               x1sig(n)=x6
            endif
         endif
      enddo
 666  continue
      close(1)

c      xbl(1)=4.0
c      xbh(1)=6.0
c      xbl(2)=6.0
c      xbh(2)=10.0
c      xbl(3)=10.0
c      xbh(3)=14.0
c      xbl(4)=14.0
c      xbh(4)=20.0
c      xbl(5)=20.0
c      xbh(5)=30.0
c      xbl(6)=30.0
c      xbh(6)=60.0
c      nbin=6
      xbl(1)=3.8
      xbh(1)=6.0
      xbl(2)=4.1
      xbh(2)=6.3
      xbl(3)=4.6
      xbh(3)=6.9
      xbl(4)=5.3
      xbh(4)=7.7
      xbl(5)=7.0
      xbh(5)=9.0
      xbl(6)=12.1
      xbh(6)=13.9
      xbl(7)=23.5
      xbh(7)=26.5
      xbl(8)=43.0
      xbh(8)=47.0
      nbin=8
      open(unit=11,file='out',status='unknown')
      do ibin=1,nbin
         nin=0
         do j=1,n
            if(xo(j).ge.xbl(ibin).and.xo(j).lt.xbh(ibin)) then
               nin=nin+1
               xin(nin)=xi(j)
               xinw(nin)=xiw(j)
               xin2(nin)=x1sig(j)
            endif
         enddo
         xm=0.
         xs=0.
         if(nin.gt.ncut) then
            call getstat(nin,xin,xinw,xm,xs,nh,xh0,xh1,xh2,
     $           x16,x84)
            xavg=(xbl(ibin)+xbh(ibin))/2.
            call biwgt(xin2,nin,xb1sig,xs1sig)
            write(*,1001) nin,xm,xs,xavg,x16,x84,xb1sig
            write(11,*) nin,xm,xs,xavg,nh,sn
c- try shifting the first bin
c            bshift=1.0
c            if(ibin.eq.1) bshift=0.8
c            do j=1,nh
c               xhnew(j)=xh0(j)*bshift
c            enddo
            do j=1,nh
c               call xlinint(xh0(j),nh,xhnew,xh1,xhnew0)
               write(11,*) xh0(j),xh1(j),xh2(j)
c               write(11,*) xh0(j),xhnew0,xhnew0
            enddo
         else
            nh=80
            write(11,*) nin,xm,xs,xavg,nh,sn
            do j=1,nh
               write(11,*) 0.,0.,0.
            enddo
         endif
      enddo
      close(11)

 1001 format(i7,6(1x,f8.2))

      end

      subroutine getstat(n,x,xw,xm,xs,nh,xh0,xh1,xh2,x16,x84)
      real x(n),xw(n),xh0(1000),xh1(1000),xh2(1000),xin(1000)

      xhmin=1.
c      xhmax=80.
c      dh=1.
      xhmax=40.
      dh=0.5
      nh=nint((xhmax-xhmin+dh)/dh)

      sum1=0.
      sum2=0.
      do i=1,nh
         xh1(i)=0.
         xh2(i)=0.
         xcen=xhmin+(xhmax-xhmin)*float(i-1)/float(nh-1)
         xl=xcen-dh/2.
         xh=xl+dh
         xh0(i)=xcen
         do j=1,n
            if(x(j).ge.xl.and.x(j).lt.xh) then
               xh1(i)=xh1(i)+1.
               xh2(i)=xh2(i)+xw(j)
            endif
         enddo
         sum1=sum1+xh1(i)
         sum2=sum2+xh2(i)
      enddo
      sumt1=0.
      sumt2=0.
      xm16=1e10
      xm84=1e10
      do i=1,nh
         xh1(i)=xh1(i)/sum1
         xh2(i)=xh2(i)/sum2
         sumt1=sumt1+xh1(i)
         sumt2=sumt2+xh2(i)
         xin(i)=sumt2
      enddo
      call xlinint(0.16,nh,xin,xh0,x16)
      call xlinint(0.84,nh,xin,xh0,x84)

      sum1=0.
      sum2=0.
      do i=1,n
         sum1=sum1+xw(i)*x(i)
         sum2=sum2+xw(i)
      enddo
      xm0=sum1/sum2

      sum1=0.
      sum2=0.
      do i=1,n
         sum1=sum1+xw(i)*(x(i)-xm0)**2
         sum2=sum2+xw(i)
      enddo
      xs=sqrt(sum1/sum2)
      xm=xm0

      return

      xcut=3.*xs
      sum1=0.
      sum2=0.
c      print *,xs
      do i=1,n
         diff=abs(x(i)-xm0)
         if(diff.lt.xcut) then
            sum1=sum1+xw(i)*x(i)
            sum2=sum2+xw(i)
         endif
      enddo
      xm=sum1/sum2

      sum1=0.
      sum2=0.
      do i=1,n
         diff=abs(x(i)-xm)
         if(diff.lt.xcut) then
            sum1=sum1+xw(i)*(x(i)-xm)**2
            sum2=sum2+xw(i)
         endif
      enddo
      xs=sqrt(sum1/sum2)

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

