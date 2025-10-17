Complete
      parameter(nmax=10000)
      real xoff(nmax),yoff(nmax),w(nmax),wo(nmax)
      real xin1(nmax),xin2(nmax),xa(nmax),ya(nmax)
      real wadc(10),adc(10),xar(nmax),yar(nmax),adca(10)
      real xar2(nmax),yar2(nmax)
      integer id(nmax)
      character c1*12
      parameter(radtodeg=57.29578)

      nadc=5
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

      chicut=5.

      open(unit=1,file='out1',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7,i8
         if(x4.gt.chicut) then
            n=n+1
            xoff(n)=3600.*cos(x2/radtodeg)*(x1-x6)
            yoff(n)=3600.*(x2-x7)
            w(n)=x3
            id(n)=i8
         endif
      enddo
 666  continue
      close(1)

      wo(1)=3600.
      wo(2)=3800.
      wo(3)=4000.
      wo(4)=4200.
      wo(5)=4400.
      wo(6)=4600.
      wo(7)=4800.
      wo(8)=5000.
      wo(9)=5200.
      wo(10)=5400.

      avgx=0.
      avgy=0.
      do iw=1,10
         nin=0
         do i=1,n
            if(abs(wo(iw)-w(i)).lt.2.) then
               nin=nin+1
               xin1(nin)=xoff(i)
               xin2(nin)=yoff(i)
            endif
         enddo
         call biwgt(xin1,nin,xb1,xs1)
         call biwgt(xin2,nin,xb2,xs2)
         xa(iw)=xb1
         ya(iw)=xb2
         avgx=avgx+xb1
         avgy=avgy+xb2
      enddo
      avgx=avgx/float(10)
      avgy=avgy/float(10)

      sumoff=0.
      do iw=1,10
         call xlinint(wo(iw),nadc,wadc,adc,adc0)
         adca(iw)=adc0
         sumoff=sumoff+adc0
      enddo
      sumoff=-sumoff/10.

      do iw=1,10
         xa(iw)=xa(iw)-avgx
         ya(iw)=ya(iw)-avgy
      enddo

      xmin=1e10
      ymin=1e10
      do i=1,360
         ang=float(i)
         ca=cos(ang/radtodeg)
         sa=sin(ang/radtodeg)
         sumx=0.
         sumy=0.
         do iw=1,10
            xar(iw)=ca*xa(iw)-sa*ya(iw)-sumoff
            yar(iw)=sa*xa(iw)+ca*ya(iw)
            sumx=sumx+abs(adca(iw)-xar(iw))
            sumy=sumy+abs(yar(iw))
         enddo
         if(sumx.lt.xmin) then
            xmin=sumx
            ang0=ang
         endif
         if(sumy.lt.ymin) then
            ymin=sumy
            ang0y=ang
         endif
      enddo
      ang=ang0y
      ca=cos(ang/radtodeg)
      sa=sin(ang/radtodeg)

      x1=a*xa(1)-sa*ya(1)-sumoff
      if(x1.gt.0.and.ang.lt.180.) ang0y=ang+180.
      if(x1.gt.0.and.ang.gt.180.) ang0y=ang-180.

      ang=ang0

      ca=cos(ang/radtodeg)
      sa=sin(ang/radtodeg)
      sumy=0.
      do iw=1,10
         xar(iw)=ca*xa(iw)-sa*ya(iw)-sumoff
         yar(iw)=sa*xa(iw)+ca*ya(iw)
         sumy=sumy+abs(yar(iw))
      enddo
      ca=cos(ang0y/radtodeg)
      sa=sin(ang0y/radtodeg)
      sumy2=0.
      do iw=1,10
         xar2(iw)=ca*xa(iw)-sa*ya(iw)-sumoff
         yar2(iw)=sa*xa(iw)+ca*ya(iw)
         sumy2=sumy2+abs(yar2(iw))
      enddo

      open(unit=11,file='fitadc.out',status='unknown')
      write(11,1101) avgx,avgy,ang,sumoff,sumy/10.,ang0,ang0y
      do iw=1,10
         write(11,1102) wo(iw),xar(iw),yar(iw),xar2(iw),yar2(iw)
      enddo
      close(11)


      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      call pgenv(3400.,5600.,-0.9,0.5,0,0)
      call pglabel('Wavelength (\(2078))','ADC (arcsec)','')
      call pgslw(5)
      call pgline(5,wadc,adc)
      call pgsci(2)
      call pgline(10,wo,xar)
      open(unit=1,file="title",status='old',err=669)
      read(1,*) c1
      call pgsci(1)
      call pgslw(3)
      call pgmtxt('T',-1.4,0.5,0.5,c1)
 669  continue
      close(1)

      call pgend

 1101 format(7(1x,f7.3))
 1102 format(1x,f6.1,4(2x,f7.3))

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
