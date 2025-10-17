Complete
      parameter (narrm1=2000,narrm2=2000,nmax2=2000,nspec=1036)
      real xd(narrm1,narrm2),xe(narrm1,narrm2),ra(nmax2),dec(nmax2)
      real spec(nspec),spece(nspec),wave(1000),diff(1000),wspec(nspec)
      real xcont(1000)
      integer irowa(1344)
      character file1*12,file2*11

      read *,file1
      read *,file2
      call getfits(file1,file2,xd,xe,nrow)
      call getpos(file1,file2,ra,dec)

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'/xwin',3,3)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)


      do i=1,1036
         wspec(i)=3470.+2.*float(i-1)
      enddo

      ntot=0
      do i=1,1344
         call getfib(nrow,ra,dec,ra(i),dec(i),3.0,nf,irowa)
         call sumfib(nf,irowa,xd,xe,spec,spece)

         call fitspabs(nspec,spec,spece,nout,wave,diff,xcont)
         do j=1,nout
            if(diff(j).gt.1.4.and.xcont(j).lt.2.) then
               ntot=ntot+1
               if(ntot.eq.10) goto 667
               iwave=nint((wave(j)-3470.)/2.)-1
               print *,iwave,i,ra(i),dec(i),wave(j),diff(j),file2
               call plotspec(wave(j),nspec,wspec,spec)
            endif
         enddo
      enddo
 667  continue
      call pgend()

      end

      subroutine plotspec(wave0,n,wave,spec)
      real wave(n),spec(n)
      xmin=wave0-150.
      xmax=wave0+150.
      ymin=1e10
      ymax=-1.e10
      do i=1,n
         if(wave(i).gt.xmin.and.wave(i).lt.xmax) then
            ymin=min(ymin,spec(i))
            ymax=max(ymax,spec(i))
         endif
      enddo

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pgline(n,wave,spec)
      return
      end

      subroutine fitspabs(n,spec,spece,nout,wave,diff,xcont)
      real spec(n),spece(n),xin(1000)
      real wave(1000),diff(1000),xcont(1000)

      nw=40
      nwh=20
      ncen=5
      ncon=nw
      nstart=nw
      nend=n-nw
      nout=0
      do i=nstart,nend
         i1=max(1,i-ncon)
         i2=i1+nwh
         i3=i2+nw
         i4=min(n,i3+nwh)
         nc=0
         do j=i1,i2
            if(spec(j).ne.0) then
               nc=nc+1
               xin(nc)=spec(j)
            endif
         enddo
         do j=i3,i4
            if(spec(j).ne.0) then
               nc=nc+1
               xin(nc)=spec(j)
            endif
         enddo
         call biwgt(xin,nc,xbc,xsc)

         na=0
         i2a=i2+18
         i2b=i2a+ncen
         do j=i2,i3
            if(spec(j).ne.0) then
               na=na+1
               xin(na)=spec(j)
            endif
         enddo
c         do j=i2b,i3
c            if(spec(j).ne.0) then
c               na=na+1
c               xin(na)=spec(j)
c            endif
c         enddo
         call biwgt(xin,na,xba,xsa)

         nout=nout+1
         wave(nout)=3470.+2.*float(i-1)
         diff1=(xbc-xba)/xsa
         diff2=(xbc-xba)/xsc
         diff(nout)=min(diff1,diff2)
         xcont(nout)=xbc
      enddo

      return
      end

      subroutine sumfib(nf,irowa,xd,xe,spec,spece)
      parameter (narrm1=2000,narrm2=2000,nmax2=2000,nspec=1036)
      real xd(narrm1,narrm2),xe(narrm1,narrm2)
      real spec(nspec),spece(nspec)
      integer irowa(nf)

      do i=1,nspec
         nin=0
         sum1=0.
         sum2=0.
         do j=1,nf
            irow=irowa(j)
            xd0=xd(i,irow)
            xe0=xe(i,irow)
            if(xd0.ne.0.and.xe0.ne.0) then
               nin=nin+1
               sum1=sum1+xd0/xe0/xe0
               sum2=sum2+1./xe0/xe0
            endif
         enddo
         if(nin.gt.0) then
            spec(i)=sum1/sum2
            spece(i)=sqrt(1./sum2)
         else
            spec(i)=0.
            spece(i)=0.
         endif
         wave=3470.+2.*float(i-1)
      enddo      

      return
      end

      subroutine getfib(n,ra,dec,ra1,dec1,rad,nf,irowa)
      real ra(n),dec(n)
      integer irowa(1344)
      parameter(pi=3.141593e0,radtodeg=57.29578)

      cosd=cos(dec1/radtodeg)
      nf=0
      do i=1,n
         r=3600.*sqrt((cosd*(ra(i)-ra1))**2 + (dec(i)-dec1)**2)
         if(r.lt.rad) then
            nf=nf+1
            irowa(nf)=i
         endif
      enddo

      return
      end      

      subroutine getpos(file1,file2,ra,dec)
      parameter (nmax2=2000)
      real ra(nmax2),dec(nmax2)
      character file1*12,file2*11,filet*80
      character a10*5,a8*28,a3*6,a9*17,amp*2

      do i=1,1344
         ra(i)=0.
         dec(i)=0.
      enddo

      filet="/scratch/03946/hetdex/detect/dithall/"//file1//".dithall"
      open(unit=1,file=filet,status='old')
      do i=1,1000000
         read(1,*,end=666) x1,x2,a3,x4,x5,x6,x7,a8,a9,a10
         if(a8(7:17).eq.file2) then
            if(a10.eq.'exp01') iexp=1
            if(a10.eq.'exp02') iexp=2
            if(a10.eq.'exp03') iexp=3
            amp=a8(19:20)
            if(amp.eq.'LL') iamp=1
            if(amp.eq.'LU') iamp=2
            if(amp.eq.'RL') iamp=3
            if(amp.eq.'RU') iamp=4
            read(a8(22:24),*) ifib
            na=(iexp-1)*448+(iamp-1)*112+ifib
            ra(na)=x1
            dec(na)=x2
         endif
      enddo
 666  continue
      close(1)

      return
      end

      subroutine getfits(file1,file2,xd,xe,nrow)
      parameter (narrm1=2000,narrm2=2000)
      real xd(narrm1,narrm2),xe(narrm1,narrm2)
      integer naxes(2)
      character file1*12,file2*11,filet*80
      logical simple,extend,anyf

      filet="cal_out/"//file1//"_"//file2//"_cal.fits"

      im1=50
      iread=0
      iext=1
      call ftopen(im1,filet,iread,iblock,ier)
      if(ier.ne.0) then
         call ftclos(im1,ier)
         write(*,*) 'Error opening image : ',filet
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)

      iext=2
      call ftmahd(im1,iext,ihd,ier)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xe,anyf,ier)

      call ftclos(im1,ier)
 706  continue
      return
      end
