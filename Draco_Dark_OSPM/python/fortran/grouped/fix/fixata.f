
      parameter(nmax=10000)
      real w(nmax),xf(nmax),xe(nmax)
      integer na(nmax),ibad(nmax),ibad2(nmax)
      character file1*100

      xcutl=0.3
      xcuth=1.8

      read *,file1
      open(unit=1,file=file1,status='old')
      open(unit=11,file='out',status='unknown')

      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,i4
         n=n+1
         w(n)=x1
         xf(n)=x2
         xe(n)=x3
         na(n)=i4
         ibad(n)=0
         if(x2.lt.xcutl) ibad(n)=1
         if(x2.gt.xcuth) ibad(n)=1
         ibad2(n)=ibad(n)
      enddo
 666  continue
      close(1)

      do i=1,n
         if(ibad(i).eq.1) then
            ibad2(max(1,i-1))=1
            ibad2(i)=1
            ibad2(min(n,i+1))=1
         endif
      enddo

      do i=1,n
         if(ibad2(i).eq.0) write(11,1101) w(i),xf(i),xe(i),na(i)
      enddo
      close(11)

 1101 format(f7.1,2(1x,f8.4),1x,i5)

      end
