Skipped
      real xl(10000),xa3(10000),xa4(10000),xa(10000,10000)
      real xin(10000)
      character file1*40
      
      open(unit=1,file='listlf',status='old')
      n=0
      do i=1,10000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old',err=667)
         n=n+1
         n2=0
         do j=1,1000
            read(2,*,end=667) x1,x2,x3,x4
            n2=n2+1
            xl(n2)=x1
            xa(n2,n)=x2
            xa3(n2)=x3
            xa4(n2)=x4
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

      do j=1,n2
         do i=1,n
            xin(i)=xa(j,i)
         enddo
         call biwgt(xin,n,xb,xs)
         print *,xl(j),xb,xa3(j),xa4(j),xs
      enddo

      end
