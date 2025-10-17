
      character a4*12,file1*80,file2*80,a1*8,a2*8,a4c*7

c- 0 means it was used
c  1 means it was rejected
c  2 means it never made it into the sample

      open(unit=1,file='listcheck',status='old')

      radm=4.
      do i=1,100000
         read(1,*,end=666) x1,x2,i3,a4
         ibad=1
         ra=x1
         dec=x2
         cosd=cos(x2/57.3)
         file1='/data/00115/gebhardt/detect/'//a4//"/runstar"
         file2='/data/00115/gebhardt/detect/'//a4//"/res/in2"
         open(unit=2,file=file1,status='old')
         do j=1,1000
            read(2,*,end=667) a1,a2,x3,x4
            rad=3600.*sqrt( (cosd*(x3-ra))**2 + (x4-dec)**2)
            if(rad.lt.radm) then
               close(2)
               goto 668
            endif
         enddo
 667     continue
         close(2)
         ibad=2
         goto 669
 668     continue
         open(unit=3,file=file2,status='old')
         do j=1,1000
            read(3,*,end=669) a4c,i2,x3,x4,i5c
            if(a2(1:5).eq.a4c(3:7)) then
               if(i5c.eq.0) ibad=0
               goto 669
            endif
         enddo
 669     continue
         close(3)
         write(*,1101) i3,a4,ibad
      enddo
 666  continue
      close(1)
 1101 format(1x,i6,1x,a12,1x,i1)
      end
