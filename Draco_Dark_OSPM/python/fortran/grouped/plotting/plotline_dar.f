
      parameter(nmax=10000)
      real x(nmax),y(nmax),y2(nmax)
      character file1*80,file2*80,c1*3

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=3500.
      xmax=5500.
      ymin=0.3
      ymax=1.0
      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("Wavelength","Aperture Correction","")

      open(unit=1,file='list',status='old')

      nl=0
      ic=0
      do il=1,1000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,1036
            read(2,*,end=667) x1,x2,x3
            n=n+1
            x(n)=x1
            y(n)=x2*0.98
            y2(n)=x3
         enddo
 667     continue
         close(2)
         call pgslw(4)
         ic=ic+1
         if(ic.eq.13) ic=1
         call pgsci(ic)
         call pgslw(3)
         call pgline(n,x,y)
         call pgsls(4)
         call pgslw(6)
         call pgline(n,x,y2)
         call pgsls(1)
      enddo
 666  continue
      close(1)

      call pgend

      end
