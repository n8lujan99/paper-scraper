
      parameter (narrm1=2000,narrm2=100000)
      real xd(narrm1,narrm2),xd2(narrm1,narrm2),ra(narrm2),dec(narrm2)
      integer naxes(2)
      character file1*180,ttype(3)*10
c      character nullstr*1,name*8,cname(narrm2)*18
      character nullstr*1,name*8,cname(narrm2)*24
      logical simple,extend,anyf

      file1="in.fits"

      im1=0
      ier=0
      iread=0
      call ftgiou(im1,ier)
      call ftopen(im1,file1,iread,iblock,ier)

      call ftmahd(im1,1,ihd,ier)
      call ftgkye(im1,'FWHM',fwhm,file1,ier)

      call ftmahd(im1,2,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)

      call ftmahd(im1,3,ihd,ier)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd2,anyf,ier)

      call ftmahd(im1,4,ihd,ier)
      call ftgkns(im1,'TTYPE',1,3,ttype,nfound,ier)
      open(unit=11,file='out',status='unknown')
      do i=1,nrow
         call ftgcve(im1,1,i,1,1,0.,ra(i),anyf,ier)
         call ftgcve(im1,2,i,1,1,0.,dec(i),anyf,ier)
         call ftgcvs(im1,3,i,1,1,nullstr,cname(i),anyf,ier)
         write(11,*) i,ra(i),dec(i),cname(i)
      enddo
      close(11)

      call ftclos(im1,ier)

 706  continue
      end
