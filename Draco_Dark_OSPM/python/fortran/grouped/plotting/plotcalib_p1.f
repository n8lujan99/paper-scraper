
      parameter(nmax=10000)
      real xin(nmax),yin(nmax),yin2(nmax),xp(15),xpos(15)
      real xv1(nmax),xv2(nmax),yc(nmax),xl(2),yl(2)
      character file1*30,c3*30,cdate*6,af(50)*6
      character a3*6,a4*3,a5*3,a6*2,aold*3
      character adate(nmax)*6,aspec(nmax)*3,acon(nmax)*3,aamp(nmax)*2

      af(1)='201701'
      af(2)='201702'
      af(3)='201703'
      af(4)='201704'
      af(5)='201705'
      af(6)='201706'
      af(7)='201707'
      af(8)='201708'
      af(9)='201709'
      af(10)='201710'
      af(11)='201711'
      af(12)='201712'
      af(13)='201801'
      af(14)='201802'
      af(15)='201803'
      af(16)='201804'
      af(17)='201805'
      af(18)='201806'
      af(19)='201807'
      af(20)='201808'
      af(21)='201809'
      af(22)='201810'
      af(23)='201811'
      af(24)='201812'
      af(25)='201901'
      af(26)='201902'
      af(27)='201903'
      af(28)='201904'
      af(29)='201905'
      af(30)='201906'
      af(31)='201907'
      af(32)='201908'
      af(33)='201909'
      af(34)='201910'
      af(35)='201911'
      af(36)='201912'
      af(37)='202001'
      nmth=37

      np=15
      do i=1,np
         xpos(i)=-7.+float(i-1)
      enddo

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      call pgenv(-7.,7.,0.,0.145,0,0)
      call pglabel("Pixels","","")
      call pgmtxt('T',1.0,0.5,0.5,"IFUSLOT= 084")

      open(unit=1,file='084.dat',status='old')
      ntot=0
      ic=0
      do i=1,10000
c         read(1,*,end=777) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,
c     $        x11,x12,x13,x14,x15,a3,a4,a5,a6
         read(1,*,end=777) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,
     $        x11,x12,x13,x14,x15,a3
         xp(1)=x1
         xp(2)=x2
         xp(3)=x3
         xp(4)=x4
         xp(5)=x5
         xp(6)=x6
         xp(7)=x7
         xp(8)=x8
         xp(9)=x9
         xp(10)=x10
         xp(11)=x11
         xp(12)=x12
         xp(13)=x13
         xp(14)=x14
         xp(15)=x15
c         call getfwhm(15,xpos,xp,0.5,xfwhm,pmax,ymax,v1,v2)
         ic=ic+1
         if(ic.eq.14) ic=1
         if(i.lt.26) il=1
         if(i.ge.26) il=4
         call pgsci(ic)
         call pgsls(il)
         call pgline(15,xpos,xp)
      enddo
 777  continue

      call pgend

      end

      subroutine getfwhm(n,x,y,frac,fwhm,xmax2,ymax2,v1,v2)
      real x(n),y(n),y2(10000)

      data big /1.e20/

c      call spline(x,y,n,0.,0.,y2)

      ymax=-big
      ymax2=-big
      do i=1,n-1
         do ia=1,9
            xp=x(i)+float(ia-1)/9.*(x(i+1)-x(i))
c            call splint(x,y,y2,n,xp,yp)
            if(yp.gt.ymax2) then
               ymax2=yp
               xmax2=xp
            endif
         enddo
         if(y(i).gt.ymax) then
            ymax=y(i)
            imax=i
         endif
      enddo

      ymax2=ymax
      xmax2=x(imax)
      yhalf=ymax2*frac

      diff=big
      x1=x(1)
      do i=1,imax-1
         if(yhalf.ge.y(i).and.yhalf.lt.y(i+1)) then
            x1=x(i)+(yhalf-y(i))/(y(i+1)-y(i))*(x(i+1)-x(i))
         endif
      enddo

      diff=big
      x2=x(n)
      do i=imax,n-1
         if(yhalf.ge.y(i+1).and.yhalf.lt.y(i)) then
            x2=x(i+1)+(yhalf-y(i+1))/(y(i)-y(i+1))*(x(i)-x(i+1))
         endif
      enddo

      fwhm=x2-x1

c - get 1st and second moment:

      sum1=0.
      sum2=0.
      do i=1,n
         sum1=sum1+x(i)*y(i)
         sum2=sum2+y(i)
      enddo
      v1=sum1/sum2
      sum1=0.
      sum2=0.
      do i=1,n
         sum1=sum1+y(i)*(x(i)-v1)**2
         sum2=sum2+y(i)
      enddo
      v2=sqrt(sum1/sum2)

      return
      end
