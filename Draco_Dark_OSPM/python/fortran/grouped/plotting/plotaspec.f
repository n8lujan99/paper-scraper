
      real xw(1036),xf(1036,100),xe(1036,100)
      real xp(1036)
      character file1*40
      
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=1,file='alist',status='old')
      do iall=1,1000
         read(1,*,end=666) file1,xw0
         xmin=xw0-20.
         xmax=xw0+20.
         open(unit=2,file=file1,status='old',err=777)
         read(2,*) ntf
         ymin=1e10
         ymax=-ymin
         do j=1,1036
            read(2,*) xw(j),(xf(j,i),xe(j,i),i=1,ntf)
            if(xw(j).gt.xmin.and.xw(j).lt.xmax) then
               do i=1,ntf
                  ymin=min(ymin,xf(j,i))
                  ymax=max(ymax,xf(j,i))
               enddo
            endif
         enddo
         close(2)

c         ymin=-0.2
c         ymax=0.5

         call pgsci(1)
         call pgslw(2)
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgslw(4)

         do i=1,ntf
            do j=1,1036
               xp(j)=xf(j,i)
            enddo
            call pgsci(i)
            call pgline(1036,xw,xp)
         enddo
 777     continue
         close(2)
      enddo
 666  continue
      close(1)

      call pgend
      end
