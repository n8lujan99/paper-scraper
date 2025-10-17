
      parameter(nmax=10000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(100)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      character file1*80,file2*80,c1*3

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=1.
      xmax=1032.
      xmin=3500.
      xmax=5500.
      xmin=0.
      xmax=6.0
      ymin=-1.5
      ymax=1.5
      ymin=0.
      ymax=0.25
      ymin=-40
      ymax=40.
      ymin=0.
      ymax=1.1
      call pgsls(1)
      call pgslw(2)

      open(unit=1,file='linelist',status='old')

      nl=0
      ic=0
c      do il=1,1000
c         do ia=1,12
c         do ia=1,70
         do ia=1,700
            read(1,*,end=666) file1,ic
            open(unit=2,file=file1,status='old')
            if(ia.eq.1) then
               c1=file1(19:21)
               call pgsci(1)
               call pgenv(xmin,xmax,ymin,ymax,0,0)
c               call pglabel('Wavelength','Flux (1e-17 cgs)','')
c               call pglabel('Wavelength','Normalized Flux','')
c               call pglabel('Wavelength','Offset from linear','')
c               call pglabel('Contrast','CDF','')
               call pglabel('Flux Noise (1e-17 ergs/cm\U2\D/s)',
     $              'Normalized Histogram','')
               call pgsch(2.0)
c               call pgmtxt('B',-1.4,0.5,0.5,c1)
               call pgsch(1.5)
            endif
            n=0
            xmaxh=0.
            do i=1,8000
               read(2,*,end=667) x1,x2
               n=n+1
               x(n)=x1
               y(n)=x2
               xmaxh=max(xmaxh,y(n))
            enddo
 667        continue
            close(2)
            do i=1,n
               y(i)=y(i)/xmaxh
            enddo
c            if(y(1).gt.10.) goto 555
c            if(y(1).lt.-22.) goto 555
c            if(y(416).gt.6.) goto 555
c            if(y(1020).lt.-3.) goto 555
c            do i=1,n
c               y(i)=y(i)-y(16)-20.
c            enddo
            call pgslw(8)
c            ic=ic+1
c            if(ic.eq.13) ic=1
            call pgsci(ic)
            call pgline(n,x,y)
            call pgslw(1)
 555        continue
         enddo
c      enddo
 666  continue
      close(1)

      call pgend

      end
