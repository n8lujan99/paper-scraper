
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax),xin(20),xine(20),yin(nmax,20)
      character file1*80,file2*80,c1*18

      wavec=50.
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

      xmin=3500.
      xmax=5500.

      open(unit=1,file='splista',status='old')

      nl=0
      ic=0
      do il=1,20000
         read(1,*,end=666) file1,wave
c         read(1,*,end=666) file1,xnorm
         xmin=wave-wavec
         xmax=wave+wavec
         open(unit=2,file=file1,status='old',err=888)
         goto 889
 888     continue
         print *,"Does not exist: ",file1
         goto 866
 889     continue
         read(2,*) ifib
         ymin=1e10
         ymax=-ymin
         n=0
         x10=0.
         do i=1,nmax
            read(2,*,end=667) x1,(xin(j),xine(j),j=1,ifib)
            if(x1.gt.xmin.and.x1.lt.xmax) then
               n=n+1
               x(n)=x1
               do j=1,ifib
                  yin(n,j)=xin(j)
                  ymin=min(ymin,xin(j))
                  ymax=max(ymax,xin(j))
               enddo
            endif
         enddo
 667     continue
         close(2)
         c1=file1(1:17)
         icp=1
         call pgsls(1)
         call pgslw(2)
         call pgsci(1)
         call pgsch(1.8)
         ybit=(ymax-ymin)/10.
         ymin=ymin-ybit
         ymax=ymax+ybit

         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgsch(1.8)
         call pglabel('Wavelength',
     $        '1e-17 ergs/cm\U2\D/s','')
c         if(il.eq.1) call pglabel('Wavelength',
c     $        'Counts','')
         ic=0
         do j=1,ifib
            do k=1,n
               y(k)=yin(k,j)
            enddo
            ic=ic+1
            call pgsci(ic)
            call pgline(n,x,y)
         enddo
         call pgsci(1)
 866     continue
         close(2)
      enddo
 666  continue
      close(1)

      call pgend

      end
