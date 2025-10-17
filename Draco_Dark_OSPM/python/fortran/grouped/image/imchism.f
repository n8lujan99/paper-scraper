
      parameter (narrm1=2000,narrm2=2000)
      real xd(narrm1,narrm2),xin(narrm1*narrm2)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      file1='in.fits'
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
c      print *,iext,ihd,ier
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      print *,naxes(1),naxes(2)
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)

      js=5
      je=105
      is=50
      ie=990
      nin=0
      do j=js,je
         do i=is,ie
            if(xd(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xd(i,j)
            endif
         enddo
      enddo
      call biwgt(xin,nin,xb,xs)
      do j=1,nrow
         do i=1,ncol
            xd(i,j)=xd(i,j)/xb
         enddo
      enddo

      ier=0
      call ftinit(51,'imars.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftp2de(51,igc,narrm1,naxes(1),naxes(2),xd,ier)
      call ftclos(51,ier)
      call ftclos(im1,ier)

 706  continue
      end
