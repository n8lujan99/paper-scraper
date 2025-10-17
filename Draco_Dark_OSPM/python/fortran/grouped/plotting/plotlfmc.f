
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),x1a(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax),ycomp(nmax)
      real yall(1000,nmax),ynl(nmax),ynu(nmax),xrat(100,1000)
      real x1ab(nmax),xlhps(10),ylhps(10),ylmin(nmax,10),ylmax(nmax,10)
      real xcompb(nmax),ycompb(nmax),ypl(nmax),yph(nmax)
      real xp(nmax),yp(nmax)
      real xratb(100,1000),xslf(nmax),yslf(nmax),yl(nmax),yh(nmax)
      real yall0(100,10,nmax),yall1(100,10,nmax),yall2(100,10,nmax)
      real*8 pstar,xlstar,alpha,dx,dy
      character file1*80,file2*80,c1*18

      iplotn=0 ! 1 for N, 0 for density
      ratcut=10.
      ratcut=1000.
      ratcut0=1e10
      del_log=0.105

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

c      call pgbegin(0,'?',2,2)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.3)
      call pgscf(2)
      call pgslw(2)

      xmin=42.3
c      xmin=42.5
      xmax=43.8
      ymin=log10(60.)
      ymax=log10(1.5e6)
      if(iplotn.eq.0) then
         ymin=log10(1.e-6)
         ymax=log10(0.01)
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
      vol=1.8 !Gpc^3
      vol=2. !Gpc^3
      vol=4.5 !Gpc^3
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

      open(unit=1,file='listmc',status='old')

      sumn=0.
      sumt=0.
      ntot=0
      do il=1,10000
         read(1,*,end=666) file1,ic,c1
         call pgsci(1)
         call pgslw(2)
         call pgenv(xmin,xmax,ymin,ymax,0,20)
         call pgsch(1.3)
         if(iplotn.eq.1) then
            call pglabel("log\D10\UL\DLy\ga","N","")
         else
            call pglabel("log\D10\UL\DLy\ga","dN/dlogL","")
         endif
         call pgsch(1.8)
         call pgmtxt('T',-2.5,0.9,1.,c1)
         call pgsch(1.2)
         call pgsci(1)
         call pgslw(3)
         call pgsls(4)

         open(unit=2,file=file1,status='old')
         n=0
         nb=0
         ntot=ntot+1
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3,x4,x5
            rat=x2/x3
            sumt=sumt+x3
            if(rat.lt.ratcut0) then
               n=n+1
               x(n)=x1
               x2=max(0.01,x2)
               x1a(i)=x1
               xcomp(i)=x2
               ycomp(i)=x3
               n1=i
               y(n)=x2
               yl(n)=max(10.,x2-x5)
               yh(n)=x2+x5
               y(n)=log10(y(n)/del_log)
               yl(n)=log10(yl(n)/del_log)
               yh(n)=log10(yh(n)/del_log)
               if(iplotn.eq.0) y(n)=y(n)-xoff
               if(iplotn.eq.0) yl(n)=yl(n)-xoff
               if(iplotn.eq.0) yh(n)=yh(n)-xoff
               xrat(il,i)=xcomp(i)/x2
            else
               sumn=sumn+x3
            endif
         enddo
 667     continue
         close(2)
         call pgsci(ic)
         call pgsls(1)
         call pgslw(8)
         call pgline(n,x,y)
         call pgslw(2)
         ip=0
         do ip0=1,n
            ip=ip+1
            xp(ip)=x(ip0)
            yp(ip)=yh(ip0)
         enddo
         do ip0=n,1,-1
            ip=ip+1
            xp(ip)=x(ip0)
            yp(ip)=yl(ip0)
         enddo
         ip=ip+1
         xp(ip)=x(1)
         yp(ip)=yh(1)
         call pgpoly(ip,xp,yp)
         call pgsci(1)
         call pgslw(2)
         call pgline(nslf,xslf,yslf)
         call pgbox('',0.,0,'bclst',0.,0)
      enddo
 666  continue

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
