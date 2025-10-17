
      parameter (narrm=10000)
      real xd(narrm,narrm),xin(narrm*narrm)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      read *,file1
      read *,iext
      read *,i1,i2,j1,j2
      read *,i1b,i2b,j1b,j2b

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


      n=0
      sum=0.
      do j=j1,j2
         do i=i1,i2
            if(xd(i,j).ne.0) then
               n=n+1
               xin(n)=xd(i,j)
               sum=sum+xin(n)
            endif
         enddo
      enddo
      call biwgt(xin,n,xb,xs)

      n=0
      sum=0.
      do j=j1b,j2b
         do i=i1b,i2b
            if(xd(i,j).ne.0) then
               n=n+1
               xin(n)=xd(i,j)
               sum=sum+xin(n)
            endif
         enddo
      enddo
      call biwgt(xin,n,xbb,xsb)

      rat=0.
      if(xb.gt.0) rat=xbb/xb

      open(unit=11,file='imstat.out',status='unknown')
      write(11,1101) n,rat,xb,xbb
      close(11)
 1101 format(i10,3(1x,1pe13.5))

 706  continue
      end
