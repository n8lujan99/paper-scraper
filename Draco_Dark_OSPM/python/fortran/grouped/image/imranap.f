
      parameter (narrm=22000)
      real xd(narrm,narrm)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      rad=1.5
      pix=0.15
      prad=rad/pix

c- ab=-2.5*log10(flux) + 30

      ig=1
      if(ig.eq.1) then
         file1="D3.g.fits"
         off=18.
      else
         file1="D3.u.fits"
         off=8
      endif
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

      iseed=-1

      open(unit=11,file='out',status='unknown')
      nt=0
      n25=0
      n26=0
      n27=0
      n28=0
      n29=0
      do iall=1,100000
         ix1=nint(ran2(iseed)*ncol)
         iy1=nint(ran2(iseed)*nrow)

         ix1=max(15,ix1)
         ix1=min(ncol-15,ix1)
         iy1=max(15,iy1)
         iy1=min(nrow-15,iy1)

      ilo=nint(float(ix1)-prad)-1
      iup=ilo+nint(2.*prad)+1
      jlo=nint(float(iy1)-prad)-1
      jup=jlo+nint(2.*prad)+1

      sum=0.
      n=0
      do i=ilo,iup
         do j=jlo,jup
            if(xd(i,j).eq.0) goto 888
            rad=sqrt((float(i-ix1))**2+(float(j-iy1))**2)
            if(rad.lt.prad) then
               n=n+1
               sum=sum+xd(i,j)
            endif
         enddo
      enddo
      sum=sum+off
      nt=nt+1
      ab=30.
      if(sum.gt.1) ab=-2.5*log10(sum)+30.
      write(11,*) ix1,iy1,sum,ab,n
      if(ab.lt.25) n25=n25+1
      if(ab.lt.26) n26=n26+1
      if(ab.lt.27) n27=n27+1
      if(ab.lt.28) n28=n28+1
      if(ab.lt.29) n29=n29+1
 888  continue

      enddo
      close(11)

      xnt=float(nt)
      print *,float(n25)/xnt,float(n26)/xnt,
     $     float(n27)/xnt,float(n28)/xnt,float(n29)/xnt

 706  continue
      end
