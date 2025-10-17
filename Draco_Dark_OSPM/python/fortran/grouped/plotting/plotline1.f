
      parameter(nmax=10000)
      real x(nmax),y1(nmax),yb(nmax),xb(nmax)
      real xsm(nmax),ysm(nmax)
      
      open(unit=1,file='in',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         x(n)=x1
         y1(n)=x2/1.e-17
c         y2(n)=x3
      enddo
 666  continue
      close(1)

      ibin=5
      ib1=(ibin-1)/2
      xib=float(ibin)
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
            sum=sum+y1(is)
            nb=nb+1
            yb(nb)=y1(is)
            xb(nb)=x(is)
         enddo
         call biwgt(yb,nb,xbb,xsb)
         ysm(nbb)=xbb
c         ysm(nbb)=sum/float(nb)
         call biwgt(xb,nb,xbb,xsb)
         xsm(nbb)=xbb
c         print *,xbb
      enddo

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

c      call pgenv(3750.,4150.,-2.0,3.0,0,0)
      call pgenv(3850.,4050.,-0.01,0.01,0,0)
c      call pgline(n,x,y1)
c      call pgsci(2)
      call pgslw(7)
      call pgline(nbb,xsm,ysm)

      open(unit=1,file='1.txt',status='old',err=667)
      n=0
      do i=1,nmax
         read(1,*,end=668) x1,x2
         n=n+1
         x(n)=x1
         y1(n)=x2/200.
      enddo
 668  continue
      call pgsci(2)
      call pgline(n,x,y1)
 667  continue
      close(1)

      call pgend

      end
      
