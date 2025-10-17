
      parameter (narrm1=2000,narrm2=2000)
      real xd(narrm1,narrm2),xin(narrm1),xin2(narrm1),xi(narrm1)
      real xy1(narrm1*narrm1),xy2(narrm1*narrm1)
      integer naxes(2),ixall(narrm1*narrm1),iyall(narrm1*narrm1)
      integer ixall2(narrm1*narrm1),iyall2(narrm1*narrm1)
      integer itrap(narrm1),jtrap(narrm1)
      integer itrap2(narrm1),jtrap2(narrm1)
      character file1*180
      logical simple,extend,anyf

      tcut=10.
      tcut2=-8.

      file1='in.fits'
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
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)

      ntrap=10
      istart=10
      iend=ncol-istart
      ihalf=9

      jstart=1
      
      do j=jstart,nrow
         nin2=0
         do i=istart,iend
            ilo=i-ihalf
            ihi=i+ihalf
            nin=0
            do is=ilo,ihi
               if(is.ne.i.and.xd(is,j).ne.0) then
                  nin=nin+1
                  xin(nin)=xd(is,j)
               endif
            enddo
            call biwgt(xin,nin,xb,xs)
            nin2=nin2+1
            if(xs.le.0) then
               xin2(nin2)=0.
            else
               xin2(nin2)=(xd(i,j)-xb)/xs
            endif
            xi(nin2)=float(i)
         enddo
         call sort2(nin2,xin2,xi)
         do ib=1,nin2
            if(xin2(ib).gt.tcut) then
               nall=nall+1
               ixall(nall)=nint(xi(ib))
               iyall(nall)=j
               xy1(nall)=xin2(ib)
            endif
            if(xin2(ib).lt.tcut2) then
               nall2=nall2+1
               ixall2(nall2)=nint(xi(ib))
               iyall2(nall2)=j
               xy2(nall2)=xin2(ib)
            endif
         enddo
      enddo

      nt=0
      nt2=0
      do i=1,ncol
         sum=0
         do iall=1,nall
            if(i.eq.ixall(iall)) sum=sum+1
         enddo
         if(sum.gt.ntrap) then
            nt=nt+1
            itrap(nt)=i
         endif
         sum2=0
         do iall=1,nall2
            if(i.eq.ixall2(iall)) sum2=sum2+1
         enddo
         if(sum2.gt.ntrap) then
            nt2=nt2+1
            itrap2(nt2)=i
         endif
      enddo

      open(unit=11,file='trap.out',status='unknown')
      do i=1,nt
         jmin=nrow
         jmax=1
         do k=1,nall
            if(ixall(k).eq.itrap(i)) jmin=min(jmin,iyall(k))
            if(ixall(k).eq.itrap(i)) jmax=max(jmax,iyall(k))
         enddo
         jmin=max(1,jmin-1)
         jmax=min(nrow,jmax+1)
         jhalf=(jmax+jmin)/2
         write(11,*) itrap(i),jmin,jmax,xy1(jmin),1
      enddo
      do i=1,nt2
         jmin=nrow
         jmax=1
         do k=1,nall2
            if(ixall2(k).eq.itrap2(i)) jmin=min(jmin,iyall2(k))
            if(ixall2(k).eq.itrap2(i)) jmax=max(jmax,iyall2(k))
         enddo
         jmin=max(1,jmin-1)
         jmax=min(nrow,jmax+1)
         jhalf=(jmax+jmin)/2
         write(11,*) itrap2(i),jmin,jmax,xy2(jmin),2
      enddo
      close(11)

c- no need to write out the new fits file
      goto 865
      do it=1,nt
         do j=jtrap(it),nrow
            xd(itrap(it),j)=0.
         enddo
      enddo

      ier=0
      call ftinit(51,'imtrap.fits',iblock,ier)
      call ftphps(51,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftp2de(51,igc,narrm1,naxes(1),naxes(2),xd,ier)
      call ftclos(51,ier)
 865  continue
      call ftclos(im1,ier)

 706  continue
      end
