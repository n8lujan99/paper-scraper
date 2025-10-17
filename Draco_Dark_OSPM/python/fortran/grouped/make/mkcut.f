
      real xin(14),xin2(4)
      character file1*120,afield*12,aspec*24,aexp*5

      wl=3510.
      wh=5490.
      snl=4.7
      snh=1e10
      chil=0.4
      chih=1.2
      sigl=1.8
      sigh=6.0
      conl=-2.5
      conh=7.0
      chifib=4.0

      open(unit=1,file='listin',status='old')
      open(unit=11,file='out',status='unknown')
      do i=1,1000000
         read(1,*,end=666) file1
         file1="output/"//file1
         open(unit=2,file=file1,status='old')
         do j=1,20000
            read(2,*,end=667,err=668) (xin(k),k=1,14),
     $           afield,(xin2(i1),i1=1,4),
     $           id,aspec,aexp,xifu,yifu,ipix,jpix,weight
            if(xin(9).gt.snl.and.xin(9).lt.snh) then
            if(xin(11).gt.chil.and.xin(11).lt.chih) then
            if(xin(1).gt.wl.and.xin(1).lt.wh) then
            if(xin(5).gt.sigl.and.xin(5).lt.sigh) then
            if(xin(7).gt.conl.and.xin(7).lt.conh) then
            if(xin2(4).lt.chifib) then
            write(11,1104) (xin(k),k=1,14),afield,(xin2(i1),i1=1,4),
     $           id,aspec,aexp,xifu,yifu,ipix,jpix,weight
            endif
            endif
            endif
            endif
            endif
            endif
         enddo
 668     continue
         print *,"something wrong ",j," ",file1(1:40)
 667     continue
         close(2)
c         print *,i,j,xin(1),xin(2),xin(3)

      enddo
 666  continue
      close(1)
      close(11)

 1104 format(12(1x,f8.2),2(1x,f10.6),1x,a12,3(1x,f6.2),1x,f5.2,1x,i5,1x,
     $     a24,1x,a5,2(1x,f6.2),1x,i4,1x,i4,1x,f5.3)

      end
