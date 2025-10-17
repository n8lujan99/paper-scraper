
      character a1*12
      parameter(radtodeg=57.29578)

      xd=1.7
      yd=xd
      open(unit=1,file='coords.in',status='old')
      do i=1,100000
         read(1,*,end=666) a1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,
     $        x13,x14,x15,x16,x17
         if(abs(x2).gt.xd) goto 888
         if(abs(x3).gt.yd) goto 888
         cosd=cos(x10/radtodeg)
         di=sqrt(x16*x16+x17*x17)
         rdiff=3600.*cosd*(x14-x9)
         ddiff=3600.*(x15-x10)
         ds=sqrt(rdiff*rdiff+ddiff*ddiff)
c         ati=radtodeg*atan2(x14/di,x15/di)
         ati=radtodeg*atan2(x17/di,x16/di)
         ats=radtodeg*atan2(rdiff/ds,ddiff/ds)
         if(ati.lt.0) ati=ati+360.
         if(ats.lt.0) ats=ats+360.
         ang=ats+180+ati
         if(ang.gt.360.) ang=ang-360.
         if(ang.lt.0.) ang=ang+360.

         ang0=270.+x12
         if(ang0.gt.360.) ang0=ang0-360.
         if(ang0.lt.0.) ang0=ang0+360.
         xdum=abs(ang0-x8)
         if(xdum.gt.80..and.xdum.lt.270.) then
            angf=x8+180.
            if(angf.gt.360.) angf=angf-360.
            if(angf.lt.0.) angf=angf+360.
         else
            angf=x8
         endif
c         ang0=ang
         print *,ati,ats,ang0,angf,x8
 888     continue
      enddo
 666  continue
      close(1)

      end
