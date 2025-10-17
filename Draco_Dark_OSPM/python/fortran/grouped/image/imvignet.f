
      parameter (narrm1=1500,narrm2=1500)
      real xd(narrm1,narrm2),xin(narrm1*narrm2)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      is=300
      ie=700

c      print *,"Image"
      read *,file1

c      print *,"Which extension"
      read *,iext

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

      n=0
      do j=1,nrow
         do i=is,ie
            if(xd(i,j).ne.0) then
               n=n+1
               xin(n)=xd(i,j)
            endif
         enddo
      enddo
      call biwgt(xin,n,xball,xs)
      cut=-2.*xs
      nt=0
      do j=1,nrow
         n=0
         do i=is,ie
            if(xd(i,j).ne.0) then
               n=n+1
               xin(n)=xd(i,j)-xball
            endif
         enddo
         call biwgt(xin,n,xb,xs)
c         print *,j,xb,xs
         if(n.gt.10.and.xb.lt.cut) nt=nt+1
      enddo
      write(*,1001) nt,xball,cut,file1

 1001 format(i3,2(1x,f8.1),1x,a102)

 706  continue
      end
