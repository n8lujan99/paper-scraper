
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),x1a(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax),ycomp(nmax)
      real yall(1000,nmax),ynl(nmax),ynu(nmax),xrat(100,1000)
      real x1ab(nmax),xlhps(10),ylhps(10),ylmin(nmax,10),ylmax(nmax,10)
      real xcompb(nmax),ycompb(nmax),ypl(nmax),yph(nmax),yp(nmax)
      real xratb(100,1000),xslf(nmax),yslf(nmax),xin(nmax)
      real yall0(100,10,nmax),yall1(100,10,nmax),yall2(100,10,nmax)
      real*8 pstar,xlstar,alpha,dx,dy
      character file1*80,file2*80,c1*18

      iplotn=0 ! 1 for N, 0 for density
      ratcut=10.
      ratcut=1000.
      ratcut0=1e10
      del_log=0.105

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
      vol=1.8 !Gpc^3
      vol=2. !Gpc^3
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

      call pgsci(1)

      call pgslw(3)

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
         nin=0
         do i=1,n
            if(x(i).gt.42.7.and.x(i).lt.43.5) then
               nin=nin+1
               call xlinint(x(i),nslf,xslf,yslf,y0)
               xin(nin)=y(i)-y0
            endif
         enddo
         call biwgt(xin,nin,xb0,xs0)
         do i=1,n
            y(i)=y(i)-xb0
         enddo
         call pgsci(ic)
         call pgline(n,x,y)
      enddo
 666  continue
      close(1)

      call pgsci(2)
      call pgslw(5)
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
