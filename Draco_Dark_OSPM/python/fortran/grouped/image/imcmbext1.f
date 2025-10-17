
      parameter (nmax=2000,narrm=2000)
      real xd1(nmax,nmax),xd2(nmax,nmax),xin(nmax*nmax)
      real xall(nmax,112)
      integer naxes(2)
      character file1*120,filea(nmax)*120,comm*120
      logical simple,extend,anyf

      open(unit=1,file="listin",status='old',err=706)

      n=0
      do i=1,2000
         read(1,*,end=666) file1
         n=n+1
         filea(n)=file1
         if(i.eq.1) print *,file1
      enddo
 666  continue

      nuse=nint(0.2*float(n))
      print *,nuse
      print *

      do i=1,n
         im1=50
         ier=0
         iread=0
         call ftopen(im1,filea(i),iread,iblock,ier)
         if(ier.ne.0) then
            write(*,*) 'Error opening image : ',filea(1)
            goto 706
         endif
         call ftmahd(im1,15,ihd,ier)
         call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=naxes(1)
         nrow=naxes(2)
         call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd1,anyf,ier)
         call ftmahd(im1,16,ihd,ier)
         call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd2,anyf,ier)

         call ftclos(im1,ier)
         if(naxis.eq.1) naxes(2)=1
         if(naxes(1).gt.narrm.or.naxes(2).gt.narrm) then
            write(*,"('Arrays too small - make narrm bigger')")
            goto 706
         endif

         do j=1,112
            nin=0
            do ia=300,800
               if(xd1(ia,j).ne.0) then
                  nin=nin+1
                  xin(nin)=xd1(ia,j)
               endif
            enddo
            call biwgt(xin,nin,xb1,xs1)
            nin=0
            do ia=300,800
               if(xd1(ia,j).ne.0) then
                  nin=nin+1
                  xin(nin)=xd2(ia,j)
               endif
            enddo
            call biwgt(xin,nin,xb2,xs2)
            if(xb1.gt.0) then
               rat=xs2/xb1
            else
               rat=-1.
            endif
            xall(i,j)=rat
         enddo
      enddo

      open(unit=11,file='out',status='unknown')
      do j=1,112
         nin=0
         do i=1,n
            if(xall(i,j).gt.0) then
               nin=nin+1
               xin(nin)=xall(i,j)
            endif
         enddo
         call biwgt(xin,nin,xb,xs)
         write(11,*) nin,j,xb,xs
      enddo
      close(11)

 706  continue
      end
