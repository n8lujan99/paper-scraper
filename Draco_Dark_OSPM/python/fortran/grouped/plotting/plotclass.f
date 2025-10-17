
      parameter(nmax=50000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax)
      character file1*80,file2*80,c1*18

      ibin=153
c      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=1140.
      xmax=1270.

      ymin=-0.1
      ymax=2.8



      open(unit=1,file='list',status='old')

      nl=0
      ic=0
      do il=1,10000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2
            if(x2.ge.-1e10.or.x2.le.1e10) then
               n=n+1
               x(n)=x1
               y(n)=x2
            else
               goto 866
            endif
c            if(x(n).eq.4740.) print *,x2,file1
         enddo
 667     continue
         close(2)
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
         enddo
      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgsch(1.2)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength (\(2078))','Normalized Flux','')
      call pgslw(5)
      call pgsch(1.)
         c1=file1(1:17)
         ic=ic+1
         icp=int(float(ic)/4.)+1
         if(icp.eq.14) ic=1
         if(ic.eq.14) ic=1
c         call pgsci(ic)
         call pgline(nbb,xn,yn)
         do k=1,nbb
            if(xn(k).gt.1210.and.xn(k).lt.1240..and.yn(k).gt.0.9) then
               print *,il
               goto 866
            endif
         enddo
 866     continue
         close(2)
      enddo
 666  continue
      close(1)

      call pgend

      end
