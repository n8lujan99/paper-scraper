
      parameter (narrm=10000)
      real xd(narrm),wd(narrm)
      character file1*40

      diff0=0.4

      open(unit=1,file='listin',status='old')
      do iall=1,100000
         read(1,*,end=777) file1

         open(unit=2,file=file1,status='old')
         nd=0
         do i=1,narrm
            read(2,*,end=667) x1,x2
            nd=nd+1
            wd(nd)=x1
            xd(nd)=x2
         enddo
 667     continue
         close(2)
         
         isize=5
         do i=10,1022
            il=i-isize
            ih=i+isize
            avg=0.
            navg=0
            if(xd(i).gt.0.2) then
               do k=il,ih
                  if(k.ne.i.and.xd(k).gt.0.) then
                     avg=avg+xd(k)
                     navg=navg+1
                  endif
               enddo
               if(navg.gt.0) then
                  avg=avg/float(navg)
                  diff=avg-xd(i)
                  if(abs(diff).gt.diff0) then
                     write(*,1001) file1(1:25),wd(i),xd(i),avg
                  endif
               else
c     print *,"not enough"
               endif
            endif
         enddo

      enddo
 777  continue
      close(1)

 1001 format(a25,1x,f6.1,1x,f5.2,1x,f5.2)
 706  continue
      end
