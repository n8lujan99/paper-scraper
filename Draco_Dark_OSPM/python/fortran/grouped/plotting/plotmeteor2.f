
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax),yb2(nmax),yn2(nmax),y2(nmax)
      character file1*80,file2*80,c1*18,cname*25

      idumb=1 ! 0 is normal, 1 if per A
      iff=1 ! 0 for one line, 1 for both lines
      nf=18
c      ibin=7
      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=3500.
      xmax=5500.

      call pgsci(1)
      call pgsch(1.8)
      ymin=-50
      ymax=8000.
      
      call pgsch(1.2)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pgsch(1.8)
      call pgmtxt('B',2.0,0.5,0.5,'Wavelength (\(2078))')
      call pgmtxt('L',1.5,0.5,0.5,
     $     '10\U-17\D ergs/cm\U2\D/s/\(2078)')

      open(unit=1,file='splist',status='old')

      nl=0
      ic=0
      yoff=0.
      do il=1,20000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         ymin=1e10
         ymax=-ymin
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
            if(x2.ne.0) then
               n=n+1
               x(n)=x1
               y(n)=x6/2.+yoff
            endif
         enddo
 667     continue
         close(2)
         yoff=yoff+370.

         ic=ic+1
         if(ic.eq.12) ic=1
         call pgsls(1)
         call pgslw(2)
         call pgsci(ic)
         call pgsch(1.8)
         call pgline(n,x,y)

 866     continue
         close(2)
      enddo
 666  continue
      close(1)

      call pgend

      end
