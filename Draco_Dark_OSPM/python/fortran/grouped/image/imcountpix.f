
      parameter (narrm=10000)
      real xd(narrm,narrm),xin(narrm*narrm),xin2(narrm*narrm)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

c 1    call qc1('Image ','imstat.def',file1)
c      call qi1('Which extension ','imstat.def',iext)
c      call savdef
c      print *,"Image"
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

c      print *,"X range"
c      read *,ix1,ix2
c      print *,"Y range"
c      read *,iy1,iy2

      ix1=1
      ix2=ncol
      iy1=1
      iy2=nrow

      n=0
      nc=0
      do j=iy1,iy2
         do i=ix1,ix2
            if(xd(i,j).gt.0.and.xd(i,j).lt.0.7) then
            nc=nc+1
            endif
         enddo
      enddo

      open(unit=11,file='imstat.out',status='unknown')
      write(11,1101) file1,nc
      close(11)
 1101 format(a20,1x,i10)

 706  continue
      end
