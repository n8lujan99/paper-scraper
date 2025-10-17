Skipped
      real xa1(1000,10),xa2(1000,10),xa3(1000,10)
      real xa4(1000,10),xa5(1000,10),xa6(1000,10),xa7(1000,10)
      real xin2(1000),xin3(1000),xin4(1000),xin5(1000),xin6(1000)
      real xin7(1000)
      character cfile*50

      nl=9

      open(unit=1,file='list',status='old')
      na=0
      do i=1,10000
         read(1,*,end=666) cfile
         na=na+1
         open(unit=2,file=cfile,status='old')
         do j=1,nl
            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7
            xa1(na,j)=x1
            xa2(na,j)=x2
            xa3(na,j)=x3
            xa4(na,j)=x4
            xa5(na,j)=x5
            xa6(na,j)=x6
            xa7(na,j)=x7
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

      open(unit=11,file='out',status='unknown')
      open(unit=12,file='out2',status='unknown')
      xn=0.
      do j=1,nl
         do i=1,na
            xin2(i)=xa2(i,j)
            xin3(i)=xa3(i,j)
            xin4(i)=xa4(i,j)
            xin5(i)=xa5(i,j)
            xin6(i)=xa6(i,j)
            xin7(i)=xa7(i,j)
         enddo
         call biwgt(xin2,na,xb2,xs2)
         call biwgt(xin3,na,xb3,xs3)
         call biwgt(xin4,na,xb4,xs4)
         call biwgt(xin5,na,xb5,xs5)
         call biwgt(xin6,na,xb6,xs6)
         call biwgt(xin7,na,xb7,xs7)
         xc2=2.5*xs2
         xc3=2.5*xs3
         xc4=2.5*xs4
         xc5=2.5*xs5
         xc6=2.5*xs6
         xc7=2.5*xs7
         x2=0.
         x3=0.
         x4=0.
         x5=0.
         x6=0.
         x7=0.
         do i=1,na
            if(xa2(i,j).gt.xb2-xc2.and.xa2(i,j).lt.xb2+xc2) then
               x2=x2+xa2(i,j)
               x3=x3+xa3(i,j)
               x4=x4+xa4(i,j)
               x5=x5+xa5(i,j)
               x6=x6+xa6(i,j)
               x7=x7+xa7(i,j)
            endif
         enddo
         rat=0.
         if(x2.gt.0.) rat=x7/x2
         write(11,1001) xa1(1,j),nint(x2),nint(x3),nint(x4),nint(x5),
     $        nint(x6),nint(x7),rat
c         write(12,2001) xa1(1,j),x2/x2,x2/x3,x2/x4,x2/x5,x2/x6,x2/x7,rat
         if(j.gt.nl-4) xn=xn+x2
      enddo
      close(11)
      xn=xn/4.
      print *,xn

      open(unit=1,file='out',status='old')
      do j=1,nl
         read(1,*,end=668) x1,x2,x3,x4,x5,x6,x7,x8
c         xn=x2
         write(12,2001) x1,xn/x2,xn/x3,xn/x4,xn/x5,xn/x6,xn/x7,x8
      enddo
 668  continue
      close(1)
      close(12)

 1001 format(f5.2,6(2x,i8),2x,f9.7)
 2001 format(f5.2,6(2x,f8.2),2x,f9.7)
      end
