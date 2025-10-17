
      parameter(nmax=1036,nfmax=100000,nmaxs=35000)
c- nmax is spectral; nfmax is fiber spectra; nmaxs is positions
      real xd(nmax,nfmax),xd2(nmax,nfmax),rat(nfmax),dect(nfmax)
      real az(nfmax),wave(nmax),ram(nfmax),decm(nfmax)
      real xb1a(nfmax),xb2a(nfmax)
      integer iexpnum(nfmax),idca(nfmax)

c- read in the fits file
      call getspec(nf,xd,xd2,rat,dect,az,iexpnum)

c- define wavelength array
      nw=1036
      do i=1,nw
         wave(i)=3470.+2.*float(i-1)
      enddo

c- get the positions with continuum
      call getposc(nf,nw,wave,xd,xd2,rat,dect,az,nc,idca,xb1a,xb2a)
      print *,nc
      open(unit=11,file='out',status='unknown')
      do i=1,nc
         ip=idca(i)
         write(11,*) rat(ip),dect(ip),xb1a(ip),xb2a(ip),ip
      enddo
      close(11)

      end

      subroutine getposc(nf,nw,wave,xd,xd2,rat,dect,az,nc,idca,
     $     xb1a,xb2a)
      parameter(nmax=1036,nfmax=100000)
c- nmax is spectral; nfmax is fiber spectra; nmaxs is positions
      real xd(nmax,nfmax),xd2(nmax,nfmax),rat(nfmax),dect(nfmax)
      real az(nfmax),wave(nmax)
      real xin1(nfmax),xin2(nfmax),xin3(nfmax),xin4(nfmax)
      real xb1a(nfmax),xb2a(nfmax),xs1a(nfmax),xs2a(nfmax)
      integer idc(nfmax),idca(nfmax)

      xsigb=4.0
      xsigr=4.0
      w1b=3700.
      w2b=4100.
      w1r=5100.
      w2r=5500.
      
      do i=1,nf
         n1=0
         n2=0
         do j=1,nw
            if(wave(j).gt.w1b.and.wave(j).lt.w2b.and.xd(j,i).ne.0.) then
               n1=n1+1
               xin1(n1)=xd(j,i)
            endif
            if(wave(j).gt.w1r.and.wave(j).lt.w2r.and.xd(j,i).ne.0.) then
               n2=n2+1
               xin2(n2)=xd(j,i)
            endif
         enddo
         call biwgt(xin1,n1,xb1,xs1)
         call biwgt(xin2,n2,xb2,xs2)
         xb1a(i)=xb1
         xb2a(i)=xb2
         xs1a(i)=xs1
         xs2a(i)=xs2
         xin3(i)=xb1
         xin4(i)=xb2
      enddo

      call biwgt(xin3,nf,xb1,xs1)
      call biwgt(xin4,nf,xb2,xs2)
      xcut1=xb1+xsigb*xs1
      xcut2=xb2+xsigr*xs2

      do i=1,nf
         idc(i)=0
         if(xb1a(i).gt.xcut1.or.xb2a(i).gt.xcut2) idc(i)=1
      enddo
      ic=0
      do i=1,nf
         if(idc(i).eq.1) then
            ic=ic+1
            idca(ic)=i
         endif
      enddo
      nc=ic

      return
      end

      subroutine getspec(nrow,xd,xd2,rat,dect,az,iexpnum)
      parameter(nmax=1036,nfmax=100000)
      real xd(nmax,nfmax),xd2(nmax,nfmax),rat(nfmax),dect(nfmax)
      real az(nfmax)
      integer naxes(2),iexpnum(nfmax)
      character ttype(3)*10,nullstr*1,cname(nfmax)*24,file1*120
      logical simple,extend,anyf

      im1=51
      ier=0
      iread=0
      call ftopen(im1,"in.fits",iread,iblock,ier)
      call ftmahd(im1,1,ihd,ier)
      call ftgkye(im1,'FWHM',fwhm,file1,ier)
      if(ier.gt.0) then
         ier=0
         fwhm=1.6
      endif
      call ftgkye(im1,'AZ',az0,file1,ier)
      if(ier.gt.0) then
         ier=0
         az0=200
      endif
      call ftmahd(im1,2,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd,anyf,ier)
      call ftmahd(im1,3,ihd,ier)
      call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd2,anyf,ier)
      call ftmahd(im1,4,ihd,ier)
      call ftgkns(im1,'TTYPE',1,3,ttype,nfound,ier)
      do i=1,nrow
         call ftgcve(im1,1,i,1,1,0.,rat(i),anyf,ier)
         call ftgcve(im1,2,i,1,1,0.,dect(i),anyf,ier)
         call ftgcvs(im1,3,i,1,1,nullstr,cname(i),anyf,ier)
c         call ftgcvj(im1,4,i,1,1,0,iexpnum(i),anyf,ier)
         iexpnum(i)=1
      enddo
      call ftclos(im1,ier)

      return
      end
