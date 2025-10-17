
      parameter (narrm1=1000,narrm2=440000)
      real xd(narrm1,narrm2),x(narrm1),y(narrm1)
      real xl(10000),yl(10000),ylf(10000)
      real*8 dx2
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      
      file1='b0_chunk0.fits'
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
c      print *,iext,ihd,ier
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

      open(unit=1,file="lw_in",status="old")
      read(1,*)
      nl=0
      do i=1,10000
         read(1,*,end=669) x1,dx2
         if(x1.lt.1211..or.x1.gt.1219.) then
            nl=nl+1
            xl(nl)=x1
            yl(nl)=sngl(dx2/2.5d40)
            yl(nl)=(yl(nl)+5.)/5.5
         endif
      enddo
 669  continue
      close(1)

      open(unit=1,file="wavelength_array.txt",status="old")
      read(1,*)
      n=0
      do i=1,1000
         read(1,*) x1
         n=n+1
         x(n)=x1
         ylf(n)=1.
         if(x1.gt.1201..and.x1.lt.1231.) ylf(n)=0.
      enddo
      close(1)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)
      do j=1,nrow
         ymax=0.
         ymin=1e10
         do i=1,ncol
            y(i)=exp(-xd(i,j))
            ymax=max(ymax,y(i))
            ymin=min(ymin,y(i))
         enddo
         print *,ymin,ymax
         call pgenv(1162.,1265.,0.,1.05,0,0)
         call pgline(n,x,y)
c         call pgsci(2)
c         call pgline(nl,xl,yl)
         call pgsci(2)
         call pgline(n,x,ylf)
         call pgsci(1)
      enddo
      call pgend

 706  continue
      end
