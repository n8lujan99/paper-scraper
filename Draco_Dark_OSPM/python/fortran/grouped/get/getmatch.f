
      parameter(nmaxs=10000)
      real fitp(nmaxs*5,8),fito(nmaxs,14)
      integer imatchd(nmaxs)
      character a15*12,as(nmaxs)*12

      radiusd=5.
      waved=3.5

      open(unit=1,file='in',status='old')
      np=0
      do i=1,nmaxs
         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,
     $        x11,x12,x13,x14,a15
         np=np+1
         fitp(np,1)=x13
         fitp(np,2)=x14
         fitp(np,3)=x1
         fitp(np,4)=x9
         fito(np,1)=x1
         fito(np,2)=x2
         fito(np,3)=x3
         fito(np,4)=x4
         fito(np,5)=x5
         fito(np,6)=x6
         fito(np,7)=x7
         fito(np,8)=x8
         fito(np,9)=x9
         fito(np,10)=x10
         fito(np,11)=x11
         fito(np,12)=x12
         fito(np,13)=x13
         fito(np,14)=x14
         as(np)=a15
      enddo
 666  continue
      close(1)

      call getmatch(np,fitp,radiusd,waved,imatchd)

      sncut=3.5
      chicut0=3.5
      sigcutl=1.4
      sigcuth=15.0
      open(unit=11,file='out',status='unknown')
      do i=1,np
         sns=fito(i,9)
         chis=fito(i,11)
         chicut=chicut0
         if(sns.gt.8.) chicut=6.
         if(imatchd(i).eq.1.and.fito(i,9).gt.sncut.and.
     $        fito(i,11).lt.chicut.and.
     $        fito(i,5).gt.sigcutl.and.fito(i,5).lt.sigcuth) then
         write(11,1101) fito(i,1),fito(i,2),fito(i,3),fito(i,4),
     $        fito(i,5),fito(i,6),fito(i,7),fito(i,8),fito(i,9),
     $        fito(i,10),fito(i,11),fito(i,12),fito(i,13),fito(i,14),
     $        as(i)
         endif
      enddo

 1101 format(12(1x,f8.2)2(1x,f10.6),1x,a12)
c 1101 format(i4,1x,12(1x,f8.2)2(1x,f10.6),1x,a12)

      end

      subroutine getmatch(np,fitp,rad,wc,imatchd)
      parameter(nmaxs=10000)
      real fitp(nmaxs*5,8)
      real w(nmaxs),xf(nmaxs),sn(nmaxs)
      real*8 dra(nmaxs),ddec(nmaxs),drad
      integer iok(nmaxs),iok2(nmaxs)
      integer imatchd(np)
      parameter(radtodeg=57.29578)

      do i=1,np
         dra(i)=dble(fitp(i,1))
         ddec(i)=dble(fitp(i,2))
         w(i)=fitp(i,3)
         sn(i)=fitp(i,4)
         iok(i)=1
         iok2(i)=1
         imatchd(i)=0
      enddo

      cosd=cos(sngl(ddec(1))/radtodeg)

      do i=1,np-1
         if(iok(i).eq.1) then
            imax=i
            snmax=sn(i)
            do j=i+1,np
               drad=dble(cosd)*(dra(i)-dra(j))**2+(ddec(i)-ddec(j))**2
               drad=3600.d0*dsqrt(drad)
               radc=sngl(drad)
               if(radc.lt.rad.and.abs(w(i)-w(j)).lt.wc) then
                  iok(j)=0
                  iok2(j)=0
                  iok2(i)=0
                  if(sn(j).ge.snmax) then
                     imax=j
                     snmax=sn(j)
                     iok(i)=0
                  endif
               endif
            enddo
            iok2(imax)=1
         endif
      enddo

      do i=1,np
         if(iok2(i).eq.1) imatchd(i)=1
      enddo

      return
      end
