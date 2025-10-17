
      parameter (narrm1=14000,narrm2=14000,nmax=1000)
      real xd(narrm1,narrm2),xd2(narrm1,narrm2)
      real xin1(narrm1),xin2(narrm1)
      real r(nmax),f1(nmax),f1s(nmax),f2(nmax),f2s(nmax)
      integer naxes(2)
      character file1*120
      logical simple,extend,anyf

      fmax=15.2
      fac=7.

      fmax=15.2
      fac=8.5

      file1='in.fits'
      iext1=1

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext1,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         print *,naxes(1),naxes(2)
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=max(1,naxes(2))
      call ftg2de(im1,igc,-666.,narrm1,ncol,nrow,xd,anyf,ier)

      file1='image.fits'
      iext2=1

      im2=0
      ier=0
      call ftgiou(im2,ier)
      iread=0
      call ftopen(im2,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im2,iext2,ihd,ier)
      call ftghpr(im2,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im2,igc,-666.,narrm1,ncol,nrow,xd2,anyf,ier)
      call ftclos(im2,ier)

      xcen=19.
      ycen=19.
      dx=0.27

      nr=15
      rmin=0.
      rmax=5.
      rd=(rmax-rmin)/float(nr-1)

      do i=1,nr
         r1=rmin+(rmax-rmin)*float(i-1)/float(nr-1)
         r2=r1+rd
         r1=r1/dx
         r2=r2/dx
         n=0
         do ix=1,ncol
            do iy=1,nrow
               rad=sqrt((float(ix)-xcen)**2+(float(iy)-ycen)**2)
               if(rad.gt.r1.and.rad.le.r2) then
                  n=n+1
                  xin1(n)=(xd(ix,iy)-fmax)/fac
                  xin2(n)=xd2(ix,iy)
               endif
            enddo
         enddo
         call biwgt(xin1,n,xb1,xs1)
         call biwgt(xin2,n,xb2,xs2)
         r(i)=dx*(r1+r2)/2.
         f1(i)=xb1
         f1s(i)=xs1
         f2(i)=xb2
         f2s(i)=xs2
         print *,r(i),f1(i),f1s(i),f2(i),f2s(i)
      enddo

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

c      call pgenv(0.,6.,-0.003,0.032,0,0)
      call pgenv(0.,6.,-0.003,0.11,0,0)
      call pglabel("radius (arcsec)",
     $     "flux (10\U-17\D ergs/cm\U2\D/s/\(2078))","")

      call pgslw(8)
      call pgsci(2)
      call pgline(nr,r,f1)
      call pgsci(4)
      call pgline(nr,r,f2)

      call pgend
      
 706  continue
      end
