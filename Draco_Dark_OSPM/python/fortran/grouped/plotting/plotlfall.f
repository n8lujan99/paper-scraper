
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

      open(unit=1,file='all.txt',status='old')
      nfit=0
      do i=1,100
         read(1,*,end=667) x1,x2
         nfit=nfit+1
         xfit(nfit)=x1
         yfit(nfit)=x2
      enddo
 667  continue
      close(1)      

      do i=1,nmax
         do ic=1,10
            ylmin(i,ic)=1e10
            ylmax(i,ic)=-1e10
         enddo
      enddo

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

c- from Konno et al.
      pstar=3.9d-4
      xlstar=0.849d43
      alpha=-1.8d0
c- from silverrush, z=2.2, Umeda et al.
c      pstar=6.9d-4
c      xlstar=0.631d43
c      alpha=-1.53d0
c- from silverrush, z=3.3, Umeda et al.
c      pstar=7.4d-3
c      xlstar=0.195d43
c      alpha=-1.19d0  !typo in paper?

      vol=1.0                   !Gpc^3
c      vol=1.2 !Gpc^3
c      vol=2. !Gpc^3
      xoff=log10(1.e9*vol)
      nslf=100
      do i=1,nslf
         xslf(i)=xmin+(xmax-xmin)*float(i-1)/float(nslf-1)
         dx=dble(xslf(i))
         dx=10.d0**dx
         dy=pstar*((dx/xlstar)**(alpha+1)) * dexp(-(dx/xlstar))
     $        *log(10.)
         yslf(i)=sngl(dlog10(dy))
         if(iplotn.eq.1) yslf(i)=sngl(dlog10(dy))+xoff
      enddo

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
      open(unit=11,file='outnorm',status='unknown')
      open(unit=12,file='out1',status='unknown')

      sumn=0.
      sumt=0.
      ntall=0
      ntall1=0
      ntall2=0
      ic=2
      nl=0
      do il=1,10000
         read(1,*,end=766) file1,xnorm
c         read(1,*,end=766) file1,xnorm,icp
         open(unit=2,file=file1,status='old',err=333)
         goto 334
 333     continue
         print *,"Not here: ",file1
         close(2)
         goto 499
 334     continue
         nl=nl+1
         filea(nl)=file1
         n=0
         nb=0
         nsum=0
         do i=1,nmax
            read(2,*,end=767) x1,x2,x3,x4
c            read(2,*,end=767) x1,x2,x3,x4,x5,x6
c            read(2,*,end=767) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
c            read(2,*,end=767) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,
c     $           x13,x14
c            x2=x14
            n=n+1
            x(n)=x1
            x2o=x2
            x2=max(0.001,x2)
c            x2=max(0.01,x6)
            x1a(i)=x1
            y(n)=x2*xnorm
            y(n)=log10(y(n)/del_log)
            yc(n)=x4
            if(iplotn.eq.0) y(n)=y(n)-xoff
            yall(nl,i)=y(n)
c            if(yc(n).le.1) yall(nl,i)=0.
            xnum(n)=x4
            nsum=nsum+nint(x4)
c            if(i.eq.6) print *,nl,x(n),y(n)
         enddo
 767     continue
         close(2)

         do i=1,n
c            if(yc(i).gt.-1) then
            if(yc(i).ge.0.5) then
               ilo=i
               goto 444
            else
               yall(nl,i)=0.
            endif
         enddo
 444     continue
         do i=n,1,-1
c            if(yc(i).ge.1.) then
            if(yc(i).ge.0.5) then
c            if(yc(i).gt.-1) then
               ihi=i
               goto 445
            else
               yall(nl,i)=0.
            endif
         enddo
 445     continue
c         ihi=n
c         ilo=1
c         print *,ilo,ihi
         np=0
c         fac1=1.5
c         fac2=1.2
c         fac3=1.1
         fac1=1.
         fac2=1.
         fac3=1.
         do i=ilo,ihi
            np=np+1
            xp(np)=x(i)
            yp(np)=y(i)
            if(np.eq.1) yp(np)=yp(np)+log10(fac1)
            if(np.eq.2) yp(np)=yp(np)+log10(fac2)
            if(np.eq.3) yp(np)=yp(np)+log10(fac3)
            if(np.eq.1) yall(nl,i)=yall(nl,i)+log10(fac1)
            if(np.eq.2) yall(nl,i)=yall(nl,i)+log10(fac2)
            if(np.eq.3) yall(nl,i)=yall(nl,i)+log10(fac3)
            xnump(np)=xnum(i)
c            if(i.lt.ihi.and.y(i).le.y(i+1)) goto 499
c            if(xp(np).gt.43.8) print *,nl,xp(np),yp(np)

c            if(xp(np).gt.42.4.and.xp(np).lt.42.5) 
c     $           print *,nl,xp(np),yp(np),filea(nl)
            if(xp(np).gt.42.95.and.xp(np).lt.43.05) 
     $           write(12,*) nl,xp(np),yp(np),filea(nl),xnorm,nsum
c            if(xp(np).gt.43.7.and.xp(np).lt.43.8) 
c     $           print *,nl,xp(np),yp(np),yall(nl,i)
         enddo
         call fitnorm(np,xp,yp,xnump,nfit,xfit,yfit,xnorm0)
         write(11,*) xnorm0

         ic=ic+1
         if(ic.eq.15) ic=2
         call pgsci(ic)
c         call pgsci(icp)
         call pgsls(1)
         call pgslw(1)
c         call pgslw(5)
         if(iplota.eq.1) call pgline(np,xp,yp)
         call pgslw(3)
 499     continue
      enddo
 766  continue
      close(1)
      close(11)
      close(12)

      call pgsci(1)
      call pgslw(5)
      call pgline(nslf,xslf,yslf)

      open(unit=11,file='out2',status='unknown')
      np=0
      do i=1,n
         nin=0
         sum=0.
         do j=1,nl
c     if(i.eq.15) print *,j,10**yall(j,i)
            if(yall(j,i).ne.0.) then
c            if(10**yall(j,i).gt.3.e-7) then
               nin=nin+1
               xin(nin)=yall(j,i)
               xin1(nin)=10**yall(j,i)
               sum=sum+xin1(nin)
            endif
         enddo
         if(i.eq.1) nin1=nin
         rat=float(nin)/float(nin1)
         sum=sum/float(nin)
         call biwgt(xin,nin,xb0,xs0)
         call biwgt(xin1,nin,xb0,xsl)
         if(x(i).lt.43.5) then
            y(i)=log10(xb0*rat)
         else
            y(i)=log10(sum)
c            y(i)=log10(sum*rat)
c            if(rat.gt.0.) y(i)=log10(sum/rat)
         endif
         if(nin.lt.15) xs0=max(0.25,xs0)
c         if(nin.gt.20) then
         if(nin.gt.2) then
            np=np+1
            xp(np)=x(i)
            yp(np)=y(i)
            yhi(np)=y(i)+xs0
            ylo(np)=y(i)-xs0
         endif
c         print *,i,nin,x(i),y(i),xs0
         write(11,*) x(i),y(i),xs0,nin
      enddo
      close(11)
      if(iplota.eq.0) then
         call pgsci(2)
         call pgslw(8)
         call pgline(np,xp,yp)
         call pgsls(4)
         call pgslw(7)
         call pgline(np,xp,yhi)
         call pgline(np,xp,ylo)
      endif

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
