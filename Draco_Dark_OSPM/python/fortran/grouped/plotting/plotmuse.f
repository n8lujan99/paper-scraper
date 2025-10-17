
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax),yb2(nmax),yn2(nmax),y2(nmax)
      character file1*12,file2*80,c1*18,cname*25

      xlya=1215.666

      idumb=1 ! 0 is normal, 1 if per A
      iff=0 ! 0 for one line, 1 for both lines
      nf=18
      ibin=7
      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      call pgbegin(0,'?',3,3)
c      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xw=30.
      xmin=xlya-xw
      xmax=xlya+xw

      open(unit=1,file='list',status='old')

      nl=0
      ic=0
      do il=1,20000
         read(1,*,end=666) file1,wave
         file2="spec/"//file1//".fits"
         z=wave/xlya-1.
         call getfits(file2,n,x,y)
c         print *,n,z,file2
         ymin=1e10
         ymax=-ymin
         do j=1,n
            x(j)=x(j)/(1.+z)
            if(x(j).gt.xmin.and.x(j).lt.xmax) then
               ymin=min(ymin,y(j))
               ymax=max(ymax,y(j))
            endif
         enddo
         nbb=0
         do j=1,n,ibin
            nbb=nbb+1
            istart=max(0,j-ib1)
            iend=istart+ibin-1
            if(iend.gt.n) then
               iend=n
               istart=n-ibin+1
            endif
            sum=0.
            nb=0
            do is=istart,iend
               sum=sum+y(is)
               nb=nb+1
               yb(nb)=y(is)
               xb(nb)=x(is)
            enddo
            call biwgt(yb,nb,xbb,xsb)
            yn(nbb)=xbb
            call biwgt(xb,nb,xbb,xsb)
            xn(nbb)=xbb
c            ymin=min(ymin,yn(nbb))
c            ymax=max(ymax,yn(nbb))
         enddo
         icp=1
         call pgsls(1)
         call pgslw(2)
         call pgsci(1)
         call pgsch(1.8)
         ybit=(ymax-ymin)/10.
         ymin=ymin-ybit/2.
         ymax=ymax+ybit

         call pgsch(1.2)
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgline(nbb,xn,yn)
         call pgsch(1.8)
         if(idumb.eq.0) then
            call pglabel('Wavelength (\(2078))',
     $           '1e-17 ergs/cm\U2\D/s','')
         else
c            call pglabel('Wavelength (\(2078))',
c     $           '10\U-17\D ergs/cm\U2\D/s/\(2078)','')
            call pgmtxt('B',2.0,0.5,0.5,'Wavelength (\(2078))')
            call pgmtxt('L',1.5,0.5,0.5,
     $           '10\U-17\D ergs/cm\U2\D/s/\(2078)')
         endif

         call pgsch(1.5)
         call pgsci(1)
 866     continue
         close(2)
      enddo
 666  continue
      close(1)

      call pgend

      end

      subroutine getfits(file1,n,x,y)
      parameter(nmax=10000)
      real x(nmax),y(nmax)
      integer naxes(2)
      character file1*80
      logical simple,extend,anyf

      im1=0
      ier=0
      iext=2
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)

      n=3681
      do i=1,n
         call ftgcve(im1,1,i,1,1,0.,wave,anyf,ier)
         call ftgcve(im1,3,i,1,1,0.,flux,anyf,ier)
         x(i)=wave
         y(i)=flux
      enddo

 706  continue
      call ftclos(im1,ier)

      return
      end
