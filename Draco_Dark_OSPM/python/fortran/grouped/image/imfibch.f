
      parameter (narrm1=15000,narrm2=15000)
      real xd(narrm1,narrm2)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      read *,file1
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
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      icol=480
      icol=1032
      x1=xd(icol,1)
      do i=2,nrow
         x2=xd(icol,i)
         diff=x2-x1
         if(diff.lt.7.) then
c            print *,i,x2-x1," ",file1(1:50)
         endif
         x1=x2
         if(i.eq.38) print *,i,x2,diff
      enddo

 706  continue
      end
