
      real xin(7),xine(7)
      character file1*20,file2*20,file3*10

      idum=-1
      file1='sn4.8.use'
      file2='sn4.8.euse'
c      file1='sn7.0.use'
c      file2='sn7.0.euse'

      open(unit=1,file=file1,status='old')
      open(unit=2,file=file2,status='old')

      file3='sn4.8_0000'
c      file3='sn7.0_0000'
      do iall=1,1000
         if(iall.lt.10) write(file3(10:10),3001) iall
         if(iall.ge.10.and.iall.lt.100) write(file3(9:10),3002) iall
         if(iall.ge.100.and.iall.lt.1000) write(file3(8:10),3003) iall
         if(iall.ge.1000.and.iall.lt.10000) write(file3(7:10),3004) iall
         open(unit=11,file=file3,status='unknown')

         read(1,1102) i0,i1,i2,i3,i4,i5,i6,i7
         read(1,1101) x0,(xin(i),i=1,7)
         read(2,*)
         read(2,*)
         write(11,1102) i0,i1,i2,i3,i4,i5,i6,i7
         write(11,1101) x0,(xin(i),i=1,7)
         ns=199
         do j=1,ns
            read(1,1101) x0,(xin(i),i=1,7)
            read(2,1101) x0,(xine(i),i=1,7)
            do i=1,7
               xin(i)=xin(i)+1.2*xine(i)*gasdev(idum)
c               xin(i)=xin(i)+xine(i)*(2.*ran2(idum)-1.)
               xin(i)=max(0.,xin(i))
               xin(i)=min(1.,xin(i))
            enddo
            write(11,1101) x0,(xin(i),i=1,7)
         enddo
         rewind(1)
         rewind(2)
         close(11)
      enddo
      close(1)
      close(2)

 1101 format(8(f10.6,1x))
 1102 format(8(i7,1x))
 3001 format(i1)
 3002 format(i2)
 3003 format(i3)
 3004 format(i4)

      end
