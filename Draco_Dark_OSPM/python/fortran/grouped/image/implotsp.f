
      parameter (narrm0=10000,narrm=3040)
      real xd(narrm,narrm0)
      real w(narrm0),x(narrm0),xb(narrm0),y1(narrm0)
      integer naxes(2)
      character file1*160,cout*10
      logical simple,extend,anyf

      nsm=11
      nsmh=5
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      cout="spXXXX.txt"

      file1='imars.fits'
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
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxes(1).gt.narrm.or.naxes(2).gt.narrm0) then
         write(*,"('Arrays too small - make narrm bigger')")
         print *,naxes
         goto 706
      endif
      ncol=naxes(1)
      nrow=max(1,naxes(2))
      call ftg2de(im1,igc,0.,narrm,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      iout=1000
      nplot=0
      ic=0
c      do i=1,nrow
      do i=500,nrow
         nb=0
         do j=1,ncol
            w(j)=3470.+float(j-1)*2.
            x(j)=xd(j,i)
            if(w(j).gt.3900..and.w(j).lt.5300.) then
               nb=nb+1
               xb(nb)=x(j)
            endif
         enddo
         call biwgt(xb,nb,xb1,xs1)
         if(xb1.gt.100..and.xb1.lt.1000.) then
            print *,i,xb1
            do j=1,ncol
               x(j)=x(j)/xb1
            enddo
            if(nplot.eq.0) then
               call pgsci(1)
               call pgslw(2)
               call pgenv(3500.,5500.,0.,2.0,0,0)
               call pglabel("Wavelength","Relative Flux","")
            endif

            do ism=1,ncol
               sum=0.
               jstart=max(ism-nsmh,1)
               jend=min(jstart+nsm,ncol)
               ns=0
               do j=jstart,jend
                  sum=sum+x(j)
                  ns=ns+1
               enddo
               y1(ism)=sum/float(ns)
            enddo

            iout=iout+1
            write(cout(3:6),1101) iout
            open(unit=11,file=cout,status='unknown')
            do j=1,ncol
               write(11,*) w(j),y1(j)
            enddo
            close(11)
            ic=ic+1
            if(ic.eq.7) ic=ic+1
            if(ic.eq.17) ic=1
            call pgsci(ic)
            call pgslw(5)
            call pgline(ncol,w,y1)
            nplot=nplot+1
            if(nplot.eq.16) then
c            if(nplot.eq.1) then
               nplot=0
            endif
         endif
      enddo

 706  continue
      call pgend
 1101 format(i4)
      end
