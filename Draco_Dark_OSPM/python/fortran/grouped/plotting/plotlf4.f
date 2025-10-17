
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),x1a(nmax)
      real xb(nmax),yb(nmax),yin(nmax),xcomp(nmax),ycomp(nmax)
      real yall(1000,nmax),ynl(nmax),ynu(nmax),xrat(100,1000)
      real x1ab(nmax),yoff(nmax)
      real xcompb(nmax),ycompb(nmax)
      real xratb(100,1000),xslf(nmax),yslf(nmax)
      real*8 pstar,xlstar,alpha,dx,dy
      character file1*80,file2*80,c1*18

      ratcut=25.
      ratcut0=1e10

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=42.3
      xmax=44.
      ymin=log10(60.)
      ymax=log10(1.5e6)

      pstar=2.2d-4
      xlstar=1.2d43
      alpha=-1.7d0
      vol=0.9 !Gpc^3
      xoff=log10(1.e9*vol**3)
      nslf=100
      do i=1,nslf
         xslf(i)=xmin+(xmax-xmin)*float(i-1)/float(nslf-1)
         dx=dble(xslf(i))
         dx=10.d0**dx
         dy=pstar * ((dx/xlstar)**alpha) * dexp(-(dx/xlstar))
         yslf(i)=sngl(dlog10(dy))+xoff
      enddo

      call pgenv(xmin,xmax,ymin,ymax,0,20)
      call pglabel("log\D10\UL\DLy\ga","N","")

      call pgline(nslf,xslf,yslf)

      call pgslw(2)
      call pgsls(4)

      open(unit=1,file='list2',status='old')

      sumn=0.
      sumt=0.
      do il=1,10000
         read(1,*,end=666) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         nb=0
         do i=1,nmax
c            read(2,*,end=667) x1b,x2b,x1,x4,x5,x6,x7,x2
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
               y(n)=log10(y(n))
               xrat(il,i)=xcomp(i)/x2
            else
               sumn=sumn+x3
            endif
         enddo
 667     continue
         close(2)
         if(ic.gt.13) ic=ic-10
         call pgsci(ic)
         call getoff(n,x,y,nslf,xslf,yslf,yoff(il))
         call pgline(n,x,y)
      enddo
 666  continue
      rewind(1)

      call pgslw(8)
      call pgsls(1)

      sumn=0.
      sumt=0.
      do il=1,10000
         read(1,*,end=766) file1,ic
         open(unit=2,file=file1,status='old')
         n=0
         nb=0
         do i=1,nmax
c            read(2,*,end=667) x1b,x2b,x1,x4,x5,x6,x7,x2
            read(2,*,end=767) x1,x2,x3
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
               y(n)=log10(y(n))
               xrat(il,i)=xcomp(i)/x2
            else
               sumn=sumn+x3
            endif
            y(n)=y(n)+yoff(il)
         enddo
 767     continue
         close(2)
         if(ic.gt.13) ic=ic-10
         call pgsci(ic)
         call pgline(n,x,y)
      enddo
 766  continue
      close(1)

 1001 format(7(1x,f10.3))

      call pgend

      end

      subroutine getoff(n,x,y,n2,x2,y2,yoff)
      real x(n),y(n),x2(n2),y2(n2),xin(10000)

      xlo=42.8
      xhi=43.5

      nin=0
      do i=1,n
         if(x(i).ge.xlo.and.x(i).le.xhi) then
            call xlinint(x(i),n2,x2,y2,yo)
            nin=nin+1
            xin(nin)=yo-y(i)
         endif
      enddo
      call biwgt(xin,nin,xb,xs)
      yuse=20.
      yoff=xb
      yoff=log10(yuse)
      print *,10**xb,yuse
      do i=1,n
         y(i)=y(i)+yoff
      enddo

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
      if(xp.ge.x(n)) yp=y(n)
      return
      end
