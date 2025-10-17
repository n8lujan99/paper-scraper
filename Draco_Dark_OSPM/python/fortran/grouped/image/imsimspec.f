
      parameter (narrm1=1000,narrm2=500000)
      real xd(narrm1,narrm2),x(narrm1),y(narrm1)
      real xl(10000),yl(10000),ylf(10000)
      real xall(1000,10000),xin(10000)
      real xall2(1000,10000)
      real*8 dx2
      integer naxes(2),jrms(1000)
      character file1*180,filea(10)*20
      logical simple,extend,anyf

c      filea(1)='b0_chunk0.fits'
c      filea(2)='b0_chunk1.fits'
c      filea(3)='b0_chunk2.fits'
c      filea(4)='b0_chunk3.fits'
c      filea(5)='b0_chunk4.fits'
c      filea(6)='b0_chunk5.fits'
c      filea(7)='b0_chunk6.fits'
c      filea(8)='b0_chunk7.fits'
c      filea(9)='b0_chunk7.fits'
      filea(1)='impact0_chunk0.fits'
      filea(2)='impact0_chunk1.fits'
      filea(3)='impact0_chunk2.fits'
      filea(4)='impact0_chunk3.fits'
      filea(5)='impact0_chunk4.fits'
      filea(6)='impact0_chunk5.fits'
      filea(7)='impact0_chunk6.fits'
      filea(8)='impact0_chunk7.fits'
      filea(9)='impact0_chunk7.fits'
      iext=1

      open(unit=11,file='out',status='unknown')

      ip1=0
      if(ip1.eq.1) then
         call pgbegin(0,'?',3,3)
      else
         call pgbegin(0,'?',1,1)
      endif
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      nall=0
      do iall=1,8
      file1=filea(iall)
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
c         if(x1.lt.1211..or.x1.gt.1219.) then
            nl=nl+1
            xl(nl)=x1
            yl(nl)=sngl(dx2/2.5d40)
c            yl(nl)=(yl(nl)+5.)/5.5
            yl(nl)=(yl(nl)+5.)/10.5
c         endif
      enddo
 669  continue
      close(1)

      open(unit=1,file="wavelength.txt",status="old")
      n=0
      do i=1,1000
         read(1,*) x1
         n=n+1
         x(n)=x1
         ylf(n)=0.6
         if(x1.gt.1201..and.x1.lt.1231.) ylf(n)=0.
c         if(x1.gt.1206..and.x1.lt.1226.) ylf(n)=0.
      enddo
      close(1)

      rms0=1.e10
      nrms=0
      do j=1,nrow
         rms=0.
         sum1=0.
         sum2=0.
         sum3=0.
         n1=0
         n2=0
         n3=0
         do i=1,ncol
            y(i)=exp(-xd(i,j))
            if(x(i).gt.1190..and.x(i).lt.1240.) rms=rms+(y(i)-ylf(i))**2
            if(x(i).gt.1190..and.x(i).lt.1200.) then
               n1=n1+1
               sum1=sum1+y(i)
            endif
            if(x(i).gt.1230..and.x(i).lt.1240.) then
               n1=n1+1
               sum1=sum1+y(i)
            endif
            if(x(i).gt.1200..and.x(i).lt.1210.) then
               n2=n2+1
               sum2=sum2+y(i)
            endif
            if(x(i).gt.1220..and.x(i).lt.1230.) then
               n2=n2+1
               sum2=sum2+y(i)
            endif
            if(x(i).gt.1160..and.x(i).lt.1180.) then
               n3=n3+1
               sum3=sum3+y(i)
            endif
            if(x(i).gt.1240..and.x(i).lt.1260.) then
               n3=n3+1
               sum3=sum3+y(i)
            endif
         enddo
         sum1=sum1/float(n1)
         sum2=sum2/float(n2)
         sumd=sum1-sum2
         sum3=sum3/float(n3)
         rms=sqrt(rms)
         if(rms.lt.rms0) then
            rms0=rms
            jmin=j
         endif
c         if(rms.lt.4.0) then
         if(rms.lt.4.0.and.sumd.gt.0.4.and.sum3.lt.0.6) then
c         if(rms.lt.5.0.and.sumd.gt.0.38.and.sum3.lt.0.6
c     $        .and.sum2.lt.0.1) then
            nrms=nrms+1
            jrms(nrms)=j
            print *,rms,sumd,sum2,sum3
         endif
         if(iplot.eq.1) then
            call pgenv(1162.,1265.,0.,1.05,0,0)
            call pgline(n,x,y)
            call pgsci(2)
            call pgline(n,x,ylf)
            call pgsci(1)
            iplot=0
         endif
      enddo

      do irms=1,nrms
         nall=nall+1
         write(11,*) filea(iall),jrms(irms)
         do i=1,ncol
            y(i)=exp(-xd(i,jrms(irms)))
            xall2(i,nall)=y(i)
         enddo
      enddo

      do i=1,ncol
         y(i)=exp(-xd(i,jmin))
         xall(i,iall)=y(i)
      enddo
      if(ip1.eq.1) then
      call pgenv(1162.,1265.,0.,1.05,0,0)
      call pgslw(3)
      call pgline(n,x,y)
      call pgsci(2)
      call pgslw(5)
      call pgline(n,x,ylf)
      call pgsci(1)
      call pgslw(2)
      endif
      print *,iall,nall,jmin,rms0

      enddo
      print *,nall
      close(11)

      open(unit=11,file='out2',status='unknown')
      do j=1,1000
         sum=0.
         do i=1,nall
            xin(i)=xall2(j,i)
            sum=sum+xin(i)
         enddo
         avg=sum/float(nall)
         call biwgt(xin,nall,xb,xs)
         y(j)=xb
c         y(j)=xin(420)
c     y(j)=avg
         write(11,*) x(j),y(j)
      enddo
      close(11)
      call pgenv(1162.,1265.,0.,1.05,0,0)
      call pgsci(2)
c     call pgline(n,x,ylf)
      call pgslw(5)
      call pgline(nl,xl,yl)
      call pgsci(1)
      call pgslw(3)
      call pgline(n,x,y)
      call pgend

 706  continue
      end
