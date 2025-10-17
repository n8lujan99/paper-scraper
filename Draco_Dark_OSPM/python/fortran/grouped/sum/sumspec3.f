
      real w(10000),xa(10000,100),xin(10000)
      character file1*40

      open(unit=1,file='slist',status='old')
      nall=0
      do iall=1,100
         read(1,*,end=666) file1
         nall=nall+1
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,10000
            read(2,*,end=667) x1,x2
            n=n+1
            w(n)=x1
            xa(n,nall)=x2
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

      open(unit=11,file='out',status='unknown')
      do i=1,n
         do j=1,nall
            xin(j)=xa(i,j)
         enddo
         call biwgt(xin,nall,xb,xs)
         write(11,*) w(i),xb,xs
      enddo
      close(11)

      end

