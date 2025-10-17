Marked
      character a1*12,file1*80,a2*6

      rad1=600.

      read *,ra,dec,rad2,wave

      open(unit=1,file="/data/00115/gebhardt/detect/radec.dat",
     $     status='old')

      cosd=cos(dec/57.3)

      open(unit=11,file='out',status='unknown')

      do i=1,100000
         read(1,*,end=666) a1,x2,x3,x4
         rad=sqrt( (cosd*(ra-x2))**2 + (dec-x3)**2 )
         rad=rad*3600.
         if(rad.lt.rad1) then
            file1="/data/00115/gebhardt/detect/"//a1//
     $           "/dithall.use"
            open(unit=2,file=file1,status='old')
            do j=1,100000
               read(2,*,end=667) x1,x2,a2,x4,x5
               rad=sqrt( (cosd*(ra-x1))**2 + (dec-x2)**2 )
               rad=rad*3600.
               if(rad.lt.rad2) then
c                  print *,a1,rad,x4,x5
                  write(11,1101) "rf1",ra,dec,3.,wave,50,1,a1,
     $                 1.7,3,3.5,0.2,11,102
                  goto 667
               endif
            enddo
 667        continue
            close(2)
         endif
      enddo
 666  continue
 1101 format(a3,2(1x,f10.5),1x,f5.2,1x,f6.1,1x,i2,1x,i1,1x,a12,
     $     1x,f3.1,1x,i1,1x,f3.1,1x,f3.1,1x,i2,1x,i3)
      close(1)
      close(11)

      end
