
      parameter (narrm1=1036,narrm2=100000)
      real xd(narrm1,narrm2),xde(narrm1,narrm2),xb2a(narrm2)
      real sp(narrm1),spe(narrm1),xba(narrm2),xin(narrm2)
      integer nha(narrm2)
      integer*8 i1,id(narrm2)
      character file1*180,cid*2

      cid="02"
      read *,cid

      file1='spf'//cid//'.fits'
      call getfits(file1,xd)
      file1='spe'//cid//'.fits'
      call getfits(file1,xde)
      file1='spt'//cid//'.txt'

      open(unit=1,file=file1,status='old')
      n=0
      do i=1,narrm2
         read(1,*,end=666) i1
         n=n+1
         id(n)=i1
      enddo
 666  continue
      close(1)

      open(unit=11,file='out',status='unknown')
      do j=1,n
         do i=1,1036
            sp(i)=xd(i,j)
            spe(i)=xde(i,j)
         enddo
         call checkspec(1036,sp,spe,id(j),xb,xb2,na,np,nn)
         write(11,*) id(j),xb,xb2,na,np,nn
      enddo
      close(11)

      end

      subroutine checkspec(n,sp,spe,id,xb,xb2,na,np,nn)
      real sp(n),spe(n),xin(n)
      real*8 dsp(3000)
      integer*8 id

      nin=0
      do i=1,n
         if(spe(i).ne.0) then
            nin=nin+1
            xin(nin)=sp(i)
         endif
      enddo
      call biwgt(xin,nin,xb,xs)

      nin=0
      do i=1,n
         if(spe(i).ne.0) then
            nin=nin+1
            xin(nin)=((sp(i)-xb)/spe(i))**2
         endif
      enddo
      call biwgt(xin,nin,xb2,xs2)

      nin=0
      do i=2,n-1
         nin=nin+1
         xin(nin)=(sp(i-1)+sp(i)+sp(i+1))/3.
      enddo
      np=0
      nn=0
      na=0
      do i=2,nin-1
         diff1=xin(i)-xin(i-1)
         diff2=xin(i+1)-xin(i)
         if(diff1.gt.0.and.diff2.gt.0) then
            np=np+1
         elseif(diff1.lt.0.and.diff2.lt.0) then
            nn=nn+1
         else
            na=na+1
         endif
      enddo

      ifft=0
      if(ifft.eq.1) then
      do i=1,1032
         dsp(i)=dble(sp(i))
         dsp(i+1032)=0.d0
      enddo
      call drealft(dsp,1032,1)
      nhigh=0
      xhigh=100.
      do i=1,1032
         k=2*i-1
         yp=sngl(
     $        dsqrt(dsp(k)*dsp(k)+dsp(k+1)*dsp(k+1)))
         if(yp.gt.xhigh) nhigh=nhigh+1
      enddo
      endif

      return
      end

      subroutine getfits(file1,xd)
      parameter (narrm1=1036,narrm2=100000)
      real xd(narrm1,narrm2)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      iext=1
      im1=51
      ier=0
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
      print *,naxes(1),naxes(2)
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

 706  continue
      return
      end
