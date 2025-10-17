Marked
      real s1(1000),s2(1000),c1(1000),c2(1000)
      character file1*40,file2*40

      chicut=1.4

      read *,file1
      read *,file2
      open(unit=1,file=file1,status='old')
      open(unit=2,file=file2,status='old')

      n1=0
      do i=1,10000
         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7,x8
         n1=n1+1
         s1(n1)=x6
         c1(n1)=x8
      enddo
 666  continue
      close(1)

      n2=0
      do i=1,10000
         read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8
         n2=n2+1
         s2(n2)=x6
         c2(n2)=x8
      enddo
 667  continue
      close(1)

      sc1=3.
      sc2=4.
      nc12=0
      nc22=0
      do i=1,n2
         if(c2(i).lt.chicut) then
            if(s2(i).gt.sc1) nc12=nc12+1
            if(s2(i).gt.sc2) nc22=nc22+1
         endif
      enddo
      nc11=0
      nc21=0
      do i=1,n1
         if(c1(i).lt.chicut) then
            if(s1(i).gt.sc1) nc11=nc11+1
            if(s1(i).gt.sc2) nc21=nc21+1
         endif
      enddo

      x1=float(nc12)/float(n2)*float(n1)
      x2=float(nc22)/float(n2)*float(n1)

      x1n=float(nc11)
      x2n=float(nc21)

      igood=0
      if(x1n.gt.2.*x1.and.x1n.ge.2.) igood=1

      open(unit=11,file="out_stat",status='unknown')
      write(11,1001) file1(1:10),nc11,nc21,x1,x2,igood
      write(*,1001) file1(1:10),nc11,nc21,x1,x2,igood
      close(11)

 1001 format(1x,a10,2x,2(1x,i4),2(1x,f6.2),1x,i3)
      end



