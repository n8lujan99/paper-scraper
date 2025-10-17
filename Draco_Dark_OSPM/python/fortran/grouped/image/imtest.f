
      parameter (narrm1=15000,narrm2=15000)
      real xd(narrm1,narrm2),xd2(narrm1,narrm2)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      print *,"Image"
      read *,file1
c      print *,"Which extension"
c      read *,iext
      iext=1

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         call ftclos(im1,ier)
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      ier=0
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)

      do j=1,nrow
         do i=1,ncol
            xd2(i,j)=xd(i,j)
         enddo
      enddo
      do j=2,nrow
         do i=516,517
c            ynew=0.33333*xd(i,j)+0.66666*xd(i,j-1)
            ynew=0.5*xd(i,j)+0.5*xd(i,j-1)
            xd2(i,j)=ynew
         enddo
      enddo

      do ja=6,nrow-6
         sum1=0.
         sum2=0.
         js=ja-5
         je=ja+5
         do j=js,je
            sum1=sum1+(xd2(516,j)+xd2(517,j))/2.
            sum2=sum2+(xd2(514,j)+xd2(515,j)+xd2(518,j)+xd2(519,j))/4.
         enddo
         print *,ja,sum1/sum2,sum1
      enddo

      ier=0
c      call ftgiou(51,ier)
      call ftinit(51,'imars.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
c      call ftcopy(im1,51,0,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftp2de(51,igc,narrm1,naxes(1),naxes(2),xd2,ier)
      call ftclos(51,ier)
      call ftclos(im1,ier)

 706  continue
      end
