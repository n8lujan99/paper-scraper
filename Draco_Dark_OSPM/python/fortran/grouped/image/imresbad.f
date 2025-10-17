
      parameter (narrm=10000)
      real xd(narrm,narrm),xin(narrm*narrm)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

c      print *,"Image"
      read *,file1
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
      naxis=2
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      if(naxes(1).gt.narrm.or.naxes(2).gt.narrm) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      xlo=-0.2
      xhi=0.2

      open(unit=11,file='out',status='unknown')
      ymin=1e10
      ymax=-ymin
      nh=0
      do j=1,nrow
         do i=1,ncol
            if(xd(i,j).lt.xlo.or.xd(i,j).gt.xhi) write(11,1101) 
     $           file1(1:8),i,j,xd(i,j)
         enddo
      enddo
      close(11)

 706  continue
 1101 format(a8,2(1x,i5),1x,f7.2)
      end
