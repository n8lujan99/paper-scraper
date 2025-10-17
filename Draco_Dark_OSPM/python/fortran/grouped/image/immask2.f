
      parameter (narrm=5000)
      real xd(narrm,narrm)
      integer naxes(2)
      character file1*80
      logical simple,extend,anyf

      read *,file1,val
      iext=1

      rad=546.
      xcen=651.
      ycen=639.

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
c         goto 1
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxes(1).gt.narrm.or.naxes(2).gt.narrm) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm,ncol,nrow,xd,anyf,ier)

      valm=0.
      do i=1,ncol
         do j=1,nrow
            if(xd(i,j).eq.valm) xd(i,j)=val
            rad0=sqrt((float(i)-xcen)**2+(float(j)-ycen)**2)
            if(rad0.gt.rad) xd(i,j)=-1000
         enddo
      enddo

      call ftclos(im1,ier)
      ier=0
      call ftinit(50,'immask.fits',iblock,ier)
      call ftphps(50,-32,naxis,naxes,ier)
c      call ftcopy(im1,50,0,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftp2de(50,igc,narrm,naxes(1),naxes(2),xd,ier)
      call ftclos(50,ier)

 706  continue
      end
