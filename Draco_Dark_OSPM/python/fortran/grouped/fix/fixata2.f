
      parameter (narrm=10000)
      real xd(narrm),wd(narrm),sd(narrm)
      real xd2(narrm),wd2(narrm),sd2(narrm)
      integer id(narrm),id2(narrm)
      character file1*40

      diff0=0.3
      xcutl=0.6
      read *,file1

      open(unit=2,file=file1,status='old')
      open(unit=11,file='out',status='unknown')
      nd=0
      do i=1,narrm
         read(2,*,end=667) x1,x2,x3,i4
         nd=nd+1
         wd(nd)=x1
         xd(nd)=x2
         sd(nd)=x3
         id(nd)=i4
         wd2(nd)=x1
         xd2(nd)=x2
         sd2(nd)=x3
         id2(nd)=i4
      enddo
 667  continue
      close(2)
      
      isize=5
      do i=10,1026
         il=i-isize
         ih=i+isize
         avg=0.
         navg=0
         if(xd(i).gt.xcutl) then
            do k=il,ih
               if(k.ne.i.and.xd(k).gt.0.) then
                  avg=avg+xd(k)
                  navg=navg+1
               endif
            enddo
            if(navg.gt.0) then
               avg=avg/float(navg)
               diff=avg-xd(i)
               if(abs(diff).gt.diff0) xd2(i)=0.
            endif
         else
            xd2(i-1)=0.
            xd2(i)=0.
            xd2(i+1)=0.
         endif
      enddo

      do i=1,nd
         xd2(i)=max(0.,xd2(i))
         xd2(i)=min(6.66,xd2(i))
         sd2(i)=max(0.,sd2(i))
         sd2(i)=min(6.66,sd2(i))
         write(11,1001) wd2(i),xd2(i),sd2(i),id2(i)
      enddo
      close(11)

 1001 format(1x,f6.1,1x,f6.3,1x,f8.5,1x,i5)
      end
