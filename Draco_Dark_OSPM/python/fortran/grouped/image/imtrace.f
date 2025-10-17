
      parameter (narrm1=15000,narrm2=15000)
      real xd(narrm1,narrm2),xin(narrm1)
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

      open(unit=11,file='out',status='unknown')
      xdiff0=1.5
      do j=1,111
         nin=0
         do i=1,ncol
            nin=nin+1
            xin(nin)=xd(i,j+1)-xd(i,j)
         enddo
         call biwgt(xin,nin,xb,xs)
         do i=1,ncol
            xd1=xd(i,j+1)-xd(i,j)
            xdiff=xd1-xb
            if(abs(xdiff).gt.xdiff0) write(11,1001) i,j,xdiff,file1
         enddo
      enddo
      close(11)

 706  continue
 1001 format(i4,1x,i4,1x,f7.2,2x,a16)
      end
