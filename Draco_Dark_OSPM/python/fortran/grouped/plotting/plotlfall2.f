
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),x1a(nmax),xin(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax),ycomp(nmax)
      real yall(9000,nmax),ynl(nmax),ynu(nmax)
      real x1ab(nmax),xlhps(10),ylhps(10),ylmin(nmax,10),ylmax(nmax,10)
      real xcompb(nmax),ycompb(nmax),ypl(nmax),xp(nmax),yp(nmax)
      real xslf(nmax),yslf(nmax),yc(nmax),xnum(nmax)
      real yhi(nmax),ylo(nmax),xfit(100),yfit(100),xnump(nmax)
      real xin1(nmax)
      real*8 pstar,xlstar,alpha,dx,dy
      character file1*80,file2*80,c1*18,filea(nmax)*40

      iplota=0 ! 1 for individual, 0 for summary
      iplotn=0 ! 1 for N, 0 for density
      ratcut=10.
      ratcut=1000.
      ratcut0=1e10
      del_log=0.105
c   (1397-568)/41253*(50/3600)^2*75 = 2.907e-4, using 75 since that is our max
c   assuming 10.9 Gpc^3 in 540 sq deg for 1.88<z<3.52      
      del_log=2.907e-4          ! this is for a single field
      del_log=del_log/4.2      ! empirical to match 
c      del_log=del_log/1.08      ! 8% on average is lost, might be more
c      del_log=del_log/log10(1./0.105)  ! since L bins are in 0.105 and we need per dex

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=42.3
c      xmin=42.5
      xmax=44.
      ymin=log10(60.)
      ymax=log10(1.5e6)
      if(iplotn.eq.0) then
         ymin=log10(3.e-6)
         ymax=log10(2.e-2)
      endif

c- or just read the table from Silverrush

      open(unit=2,file='silver.txt',status='old')
      nsil=0
      do i=1,nmax
         read(2,*,end=669) x1,x2,x3
         nsil=nsil+1
         xslf(nsil)=x1
         yslf(nsil)=(x2+x3)/2.
      enddo
 669  continue
      close(2)
      nslf=nsil

      call pgenv(xmin,xmax,ymin,ymax,0,20)
      if(iplotn.eq.1) then
         call pglabel("log\D10\UL\DLy\ga","N","")
      else
         call pglabel("log\D10\UL\DLy\ga","dN/dlogL","")
      endif

      call pgsci(1)

      call pgslw(3)
      call pgsls(1)

      open(unit=1,file='listin',status='old')

      sumn=0.
      sumt=0.
      ntall=0
      ntall1=0
      ntall2=0
      ic=2
      nl=0
      do il=1,10000
         read(1,*,end=766) file1
         open(unit=2,file=file1,status='old',err=333)
         goto 334
 333     continue
c         print *,"Not here: ",file1
         close(2)
         goto 499
 334     continue
         nl=nl+1
         filea(nl)=file1
         n=0
         nb=0
         nsum=0
         do i=1,nmax
            read(2,*,end=767) x1,x2
            n=n+1
            xp(n)=x1
            yp(n)=x2
         enddo
 767     continue
         close(2)
         np=n

         ic=ic+1
         if(ic.eq.15) ic=2
         call pgsci(ic)
         call pgsls(1)
         call pgslw(4)
         call pgline(np,xp,yp)
 499     continue
      enddo
 766  continue
      close(1)

      print *,nslf
      call pgsci(1)
      call pgslw(6)
      call pgline(nslf,xslf,yslf)

 1001 format(7(1x,f10.3))

      call pgend

      end

      subroutine fitnorm(np,xp,yp,xnump,nfit,xfit,yfit,xnorm0)
      real xp(np),yp(np),xnump(np),xfit(nfit),yfit(nfit)

      xmin=42.7
      xmax=43.1
c      xmin=43.1
c      xmax=43.2
      xmin=42.9
      xmax=43.3
      sum=0.
      sum2=0.
      nd=0
      do i=1,np
         if(xp(i).ge.xmin.and.xp(i).le.xmax) then
            call xlinint(xp(i),nfit,xfit,yfit,y0)
            diff=yp(i)-y0
c            sum=sum+diff
            sum=sum+xnump(i)*diff
            sum2=sum2+xnump(i)
            nd=nd+1
         endif
      enddo
c      if(nd.gt.0) sum=sum/float(nd)
      if(sum2.gt.0) sum=sum/sum2
      xnorm0=sum

      return
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end
