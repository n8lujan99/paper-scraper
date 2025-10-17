
      parameter (narr1=1036,narr2=70000)
      real xd(narr1,narr2),xin(narr2)
      integer naxes(2)
      character file1*40
      logical simple,extend,anyf

      file1='in.fits'
c      iext=2
      iext=3

      im1=50
      ier=0
c      call ftgiou(im1,ier)
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
      if(naxes(1).gt.narr1.or.naxes(2).gt.narr2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narr1,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      iy1=1
      iy2=nrow

      open(unit=11,file='imstat.out',status='unknown')
      open(unit=12,file='out',status='unknown')
      do i=1,ncol
         wave=3470.+float(i-1)*2.
         n=0
         sum=0.
         do j=iy1,iy2
            if(xd(i,j).ne.0) then
               n=n+1
               xin(n)=xd(i,j)
               sum=sum+xin(n)
               if(i.eq.500) write(12,*) xd(i,j)
            endif
         enddo
         if(n.gt.10) then
            call sort(n,xin)
            nc=nint(float(n)*0.95)
            call biwgt(xin,nc,xb,xs)
            write(11,1101) wave,n,xb,xs,sum
         endif
      enddo
      close(11)
      close(12)

 1101 format(f6.1,1x,i10,3(1x,1pe13.5))

 706  continue
      end
