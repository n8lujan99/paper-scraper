
      parameter (narrm1=1036,narrm2=10000)
      real xd(narrm1,narrm2),xde(narrm1,narrm2),xdf(narrm1,narrm2)
      real xa(narrm1,narrm2),xin(narrm2),xres(100,narrm1)
      real xaf(narrm1,narrm2),xresf(100,narrm1)
      real wave(narrm1),ywave(narrm1)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      frac=0.94

      nc=0
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)
      call pgenv(3500.,4000.,-0.1,0.1,0,0)
      call pgslw(1)

      open(unit=1,file="list",status='old')
      do iall=1,100
         read(1,*,end=669) file1
         nall=nall+1

      ier=0
      iread=0
      call ftopen(51,file1,iread,iblock,ier)
      call ftmahd(51,1,ihd,ier)
      call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(51,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
      call ftmahd(51,2,ihd,ier)
      call ftg2de(51,igc,0.,narrm1,ncol,nrow,xde,anyf,ier)
      call ftmahd(51,5,ihd,ier)
      call ftg2de(51,igc,0.,narrm1,ncol,nrow,xdf,anyf,ier)
      call ftclos(51,ier)

      w0=3470.
      dw=2.0
      w1=3500.
      w2=5500.
      i1=nint((w1-w0)/dw)+1
      i2=nint((w2-w0)/dw)+1

      n=i2-i1+1
      do j=1,nrow
         n2=0
         do i=i1,i2
            n2=n2+1
            xv=xd(i,j)
            xvf=xdf(i,j)
            if(xde(i,j).eq.0.or.xd(i,j).eq.0.) then
               xv=0.
               xvf=0.
            endif
            xa(n2,j)=xv
            xaf(n2,j)=xvf
         enddo
      enddo
      
      do i=1,n
         nin=0
         do j=1,nrow
            if(xa(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xa(i,j)
            endif
         enddo
         call biwgt(xin,nin,xb,xs)
         nin=nint(frac*float(nin))
         call biwgt(xin,nin,xb,xs)
         xres(nall,i)=xb
      enddo

      do i=1,n
         nin=0
         do j=1,nrow
            if(xaf(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xaf(i,j)
            endif
         enddo
         call biwgt(xin,nin,xb,xs)
         nin=nint(frac*float(nin))
         call biwgt(xin,nin,xb,xs)
         xresf(nall,i)=xb
         wave(i)=3500.+2.*float(i)
         ywave(i)=xb
      enddo
      nc=nc+1
      call pgsci(nc)
      if(nc.gt.14) nc=0
c      call pgline(n,wave,ywave)

      enddo
 669  continue
      close(1)
      call pgend

      open(unit=11,file='out',status='unknown')
      do j=1,n
         nin=0
         do i=1,nall
            if(xres(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xres(i,j)
            endif
         enddo
         call biwgt(xin,nin,xb1,xs)
         nin=0
         do i=1,nall
            if(xresf(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xresf(i,j)
            endif
         enddo
         call biwgt(xin,nin,xb2,xs)
         write(11,*) 3500.+2.*float(j-1),xb1,xb2,xs
      enddo
      close(11)

      end
