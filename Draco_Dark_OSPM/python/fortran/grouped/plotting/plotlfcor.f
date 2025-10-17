
      parameter(nmax=10000)
      real x(nmax),y(nmax,100),yc(nmax,100),xp(100),yp(100)
      real xin(nmax),xin2(nmax),yp2(100)
      character file1*40

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=42.3
      xmax=44.
      ymin=0.
      ymax=1.6

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("log\D10\UL\DLy\ga","N_corr","")

      open(unit=1,file='listin',status='old')

      nl=0
      do il=1,10000
         read(1,*,end=766) file1,xnorm
         open(unit=2,file=file1,status='old',err=333)
         goto 334
 333     continue
         print *,"Not here: ",file1
         close(2)
         goto 767
 334     continue
         nl=nl+1
         n=0
         do i=1,nmax
            read(2,*,end=767) x1,x2,x3,x4,x5
            n=n+1
            x(n)=x1
            y(nl,n)=x2
            yc(nl,n)=x5
         enddo
 767     continue
         close(2)
      enddo
 766  continue
      close(1)

      np=0
      do i=1,n
         sum1=0.
         sum2=0.
         do j=1,nl
            xin(j)=y(j,i)
            xin2(j)=yc(j,i)
            sum1=sum1+xin(j)
            sum2=sum2+xin2(j)
         enddo
         call biwgt(xin,nl,xb,xs)
         xn=xb
         call biwgt(xin2,nl,xb,xs)
         xc=xb
         if(xn.gt.0..and.xc.gt.0.) then
            cor=xc/xn
            np=np+1
            xp(np)=x(i)
            yp(np)=cor
            yp2(np)=sum2/sum1
         else
            cor=0.
         endif
         print *,x(i),cor,xn,xc,sum1,sum2
      enddo

      call pgsci(2)
      call pgslw(8)
      call pgline(np,xp,yp)
      call pgsci(4)
      call pgline(np,xp,yp2)
      call pgend

      end
