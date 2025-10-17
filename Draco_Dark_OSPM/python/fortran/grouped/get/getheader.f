Marked
      parameter (narrm=2000)
      integer naxes(2)
      character file1*40,cspec*100
      logical simple,extend,anyf

      file1='in.fits'
      iext=1

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      call ftgkye(im1,'SPECID',specid,file1,ier)
      print *,specid

      call ftclos(im1,ier)

 706  continue
      end

