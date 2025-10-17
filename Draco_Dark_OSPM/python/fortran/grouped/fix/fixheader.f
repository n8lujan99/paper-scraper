
      parameter (narrm=2000)
      integer naxes(2)
      character file1*80,cspec*100,cname*8,c0*3
      logical simple,extend,anyf

      read *,file1
      read *,c0
      iext=1

      im1=51
      ier=0
      iread=1
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      call ftdkey(51,"IFUID",ier)
      ier=0
      call ftpkys(51,'IFUID',c0,"IFU head ID",ier)
      call ftclos(im1,ier)

 706  continue
      end

