
      parameter (narrm1=1036,narrm2=35000)
      real xd(narrm1,narrm2),xin(narrm2*500)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      frac=0.93

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
c      print *,naxes(1),naxes(2)
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      w0=3470.
      dw=2.0

      w1=3500.
      w2=3860.

      i1=nint((w1-w0)/dw)+1
      i2=nint((w2-w0)/dw)+1
      nin=0
      do i=i1,i2
         do j=1,nrow
            if(xd(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xd(i,j)
            endif
         enddo
      enddo

      call biwgt(xin,nin,xb,xs)
      nin=nint(frac*float(nin))
      call biwgt(xin,nin,xb,xs)
      b1=xb
      
      w1=5090
      w2=5500.

      i1=nint((w1-w0)/dw)+1
      i2=nint((w2-w0)/dw)+1
      nin=0
      do i=i1,i2
         do j=1,nrow
            if(xd(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xd(i,j)
            endif
         enddo
      enddo

      call biwgt(xin,nin,xb,xs)
      nin=nint(frac*float(nin))
      call biwgt(xin,nin,xb,xs)
      b2=xb
      
      print *,b1,b2," ",file1(2:18)

 706  continue
      end
