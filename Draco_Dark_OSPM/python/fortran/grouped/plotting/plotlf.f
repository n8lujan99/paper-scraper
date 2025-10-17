
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),x1a(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax),ycomp(nmax)
      real yall(1000,nmax),ynl(nmax),ynu(nmax),xrat(100,1000)
      real x1ab(nmax),xlhps(10),ylhps(10),ylmin(nmax,10),ylmax(nmax,10)
      real xcompb(nmax),ycompb(nmax),ypl(nmax),yph(nmax),yp(nmax)
      real xratb(100,1000),xslf(nmax),yslf(nmax)
      real yall0(100,10,nmax),yall1(100,10,nmax),yall2(100,10,nmax)
      real*8 pstar,xlstar,alpha,dx,dy
      character file1*80,file2*80,c1*18

      iplotn=0 ! 1 for N, 0 for density
      ratcut=10.
      ratcut=1000.
      ratcut0=1e10
      del_log=0.105
c      del_log=1e-5  ! this is for a single field
      del_log=2.81e-5  ! this is for a single field

      xlhps(1)=42.61
      xlhps(2)=42.77
      xlhps(3)=42.88
      xlhps(4)=43.00
      xlhps(5)=43.12
      xlhps(6)=43.25
      xlhps(7)=43.38
      xlhps(8)=43.50
      ylhps(1)=1.0e-3
      ylhps(2)=6.0e-4
      ylhps(3)=3.3e-4
      ylhps(4)=2.0e-4
      ylhps(5)=1.7e-4
      ylhps(6)=1.1e-4
      ylhps(7)=2.2e-5
      ylhps(8)=5.0e-5
      do i=1,8
         ylhps(i)=log10(ylhps(i))
      enddo
      do i=1,nmax
         do ic=1,10
            ylmin(i,ic)=1e10
            ylmax(i,ic)=-1e10
         enddo
      enddo

      call pgbegin(0,'?',1,1)
c      call pgbegin(0,'/null',1,1)
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
         ymin=log10(1.e-6)
         ymax=log10(1.e-2)
c         ymax=log10(1.)
      endif

c- from HPS
      pstar=2.2d-4
      xlstar=1.2d43
      alpha=-1.7d0
c- from Ouchi
      pstar=3.9d-4
      xlstar=0.849d43
      alpha=-1.8d0
      vol=1.2 !Gpc^3
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

      call pgenv(xmin,xmax,ymin,ymax,0,20)
      if(iplotn.eq.1) then
         call pglabel("log\D10\UL\DLy\ga","N","")
      else
         call pglabel("log\D10\UL\DLy\ga","dN/dlogL","")
      endif

c      call pgsci(1)
c      call pgline(nslf,xslf,yslf)
c      call pgpt(8,xlhps,ylhps,17)
      call pgsci(1)

      call pgslw(3)
      call pgsls(4)

      open(unit=1,file='listin',status='old')

      sumn=0.
      sumt=0.
      ntot=0
      do il=1,10000
         read(1,*,end=666) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         nb=0
         ntot=ntot+1
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3
            rat=x2/x3
            sumt=sumt+x3
            if(rat.lt.ratcut0) then
               n=n+1
               x(n)=x1
               x2=max(0.01,x2)
               if(il.eq.1) then
                  x1a(i)=x1
                  xcomp(i)=x2
                  ycomp(i)=x3
                  n1=i
                  y(n)=x2
               else
                  call xlinint(x1,n1,x1a,ycomp,y0)
                  y(n)=x2
               endif
               y(n)=log10(y(n)/del_log)
               if(iplotn.eq.0) y(n)=y(n)-xoff
               xrat(il,i)=xcomp(i)/x2
            else
               sumn=sumn+x3
            endif
         enddo
 667     continue
         close(2)
         call pgsci(ic)
         call pgline(n,x,y)
      enddo
 666  continue
      rewind(1)

      if(ntot.lt.20) call pgslw(15)
      if(ntot.ge.20) ratcut=ratcut0
      call pgsls(1)

      sumn=0.
      sumt=0.
      ntall=0
      ntall1=0
      ntall2=0
      do il=1,10000
         read(1,*,end=766) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         nb=0
         if(ic.eq.1) ntall1=ntall1+1
         if(ic.eq.2) ntall2=ntall2+1
         do i=1,nmax
            read(2,*,end=767) x1,x2,x3,x4
            rat=x2/x3
            sumt=sumt+x3
            if(rat.lt.ratcut) then
               n=n+1
               x(n)=x1
               x2=max(0.01,x2)
               if(il.eq.1) then
                  x1a(i)=x1
                  xcomp(i)=x2
                  ycomp(i)=x3
                  n1=i
                  y(n)=x2
               else
                  call xlinint(x1,n1,x1a,ycomp,y0)
                  y(n)=x2
               endif
               y(n)=log10(y(n)/del_log)
               if(iplotn.eq.0) y(n)=y(n)-xoff
               xrat(il,i)=xcomp(i)/x2
               if(x4.gt.0) then
                  pfrac=sqrt(x4)/x4
               else
                  pfrac=0.9
               endif
               ypl(n)=log10(10**y(n)-pfrac*10**y(n))
               yph(n)=log10(10**y(n)+pfrac*10**y(n))
               ylmin(n,ic)=min(ylmin(n,ic),ypl(n))
               ylmax(n,ic)=max(ylmax(n,ic),yph(n))
               if(ic.eq.1) then
                  yall0(n,ic,ntall1)=y(n)
                  yall1(n,ic,ntall1)=ypl(n)
                  yall2(n,ic,ntall1)=yph(n)
               elseif(ic.eq.2) then
                  yall0(n,ic,ntall2)=y(n)
                  yall1(n,ic,ntall2)=ypl(n)
                  yall2(n,ic,ntall2)=yph(n)
               endif
            else
               sumn=sumn+x3
            endif
         enddo
 767     continue
         close(2)
         call pgsci(ic)
         if(ntot.lt.20) call pgline(n,x,y)
c         call pgline(n,x,ypl)
c         call pgline(n,x,yph)
      enddo
 766  continue
      close(1)

      if(ntot.lt.20) goto 888

      call pgslw(10)
      do ic=1,2
         call pgsci(ic)
         do i=1,n
            y(i)=ylmin(i,ic)
         enddo
c         call pgline(n,x,y)
         do i=1,n
            y(i)=ylmax(i,ic)
         enddo
c         call pgline(n,x,y)
      enddo

      do ic=1,2
         call pgsci(ic)
         if(ic.eq.1) ntall=ntall1
         if(ic.eq.2) ntall=ntall2
         do i=1,n
            nt=0
            do j=1,ntall
               nt=nt+1
               yin(nt)=yall0(i,ic,j)
               nt=nt+1
               yin(nt)=yall1(i,ic,j)
               nt=nt+1
               yin(nt)=yall2(i,ic,j)
            enddo
            call biwgt(yin,nt,xb1,xs1)
c            yup=xb1+5.*xs1
c            nup=0
c            do j=1,nt
c               if(yin(j).gt.yup) then
c                  nup=nup+1
c                  print *,ic,j,yup,yin(j)
c               endif
c            enddo
c            print *,nup
            n1=5
            n2=nt-5+1
c            n1=nint(float(nt)*0.01)
c            n2=nint(float(nt)*0.99)
c            print *,nt,n1,n2
            yp(i)=xb1
            ypl(i)=yin(n1)
            yph(i)=yin(n2)
         enddo
c         call pgline(n,x,yp)
c         call pgsls(4)
         call pgline(n,x,ypl)
         call pgline(n,x,yph)
         call pgsls(1)
      enddo
 888  continue

      nt=il-1
      do i=1,n
c         write(*,1001) x(i),(xrat(j,i),j=1,nt)
      enddo

      call pgsci(1)
      call pgslw(2)
      call pgline(nslf,xslf,yslf)

 1001 format(7(1x,f10.3))

      call pgend

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
