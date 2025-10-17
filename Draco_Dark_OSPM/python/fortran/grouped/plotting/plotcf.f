
      parameter(nmax=10000)
      real x(nmax),y1(nmax),y2(nmax),xl(2),yl(2)
      real xcf(nmax),ycf(nmax),xcfs(nmax),ycfs(nmax),xps2(nmax)
      real x1ps(nmax),y1ps(nmax),xt(nmax),yt(nmax),xps(nmax)
      real xpall(nmax,10),ypall(nmax,10),xin(nmax),yin(nmax)
      integer nall(10)
      character file1*40

      ilog=0
      
c      xc1=0.04
c      xc2=0.065
      xc1=0.045
      xc2=0.07
      xc3=0.085
      xc4=0.14
c      if(ilog.eq.1) then
         xc1=log10(xc1)
         xc2=log10(xc2)
         xc3=log10(xc3)
         xc4=log10(xc4)
c      endif
      
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.2)
      call pgscf(2)
      call pgslw(2)

      open(unit=2,file='list',status='old')

      ic=1
      nf=0
      do iall=1,10
         read(2,*,end=667) file1
         nf=nf+1

c      open(unit=1,file='in',status='old')
      open(unit=1,file=file1,status='old')
      n=0
      n2=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         x1=log10(x1)
         x2=log10(x2)
         xcf(n)=x1
         ycf(n)=x2
         x(n)=xcf(n)
         if(x1.lt.xc1) then
            n2=n2+1
            xcfs(n2)=x1
            ycfs(n2)=x2
         endif
         if(x1.gt.xc2.and.x1.lt.xc3) then
            n2=n2+1
            xcfs(n2)=x1
            ycfs(n2)=x2
         endif
         if(x1.gt.xc4) then
            n2=n2+1
            xcfs(n2)=x1
            ycfs(n2)=x2
         endif
      enddo
 666  continue
      close(1)

      xmin=0.015
      xmax=0.9
      if(ilog.eq.0) xmax=0.3
      if(ilog.eq.1) then
         xmin=log10(xmin)
         xmax=log10(xmax)
      endif
c      ymin=log10(2.2e3)
c      ymax=log10(8.e4)
      ymin=log10(7.2e2)
      ymax=log10(8.e5)

      ns=n
      call smooth(n2,xcfs,ycfs,ns,x,y1,y2)

      do i=1,ns
         y2(i)=10**ycf(i)/10**y1(i)
         if(ilog.eq.0) then
            xps(i)=10**xcf(i)
            xps2(i)=10**x(i)
         else
            xps(i)=xcf(i)
            xps2(i)=x(i)
         endif
      enddo

      call pgvport(0.15,0.85,0.5,0.85)
      call pgwindow(xmin,xmax,ymin,ymax)
      if(ilog.eq.0) then
         call pgbox('bcst',0.,0,'bcnstl',0.,0)
      else
         call pgbox('bcstl',0.,0,'bcnstl',0.,0)
      endif
      call pgsch(1.2)
      call pgmtxt('L',2.1,0.5,0.5,"P(k)")

      ic=ic+1
      call pgsci(ic)
      call pgslw(4)
      call pgline(n,xps,ycf)
      call pgsci(1)
c      call pgline(ns,x,y1)
      call pgslw(2)

      ymin=0.95
      ymax=1.18
      call pgvport(0.15,0.85,0.15,0.5)
      call pgwindow(xmin,xmax,ymin,ymax)
      if(ilog.eq.0) then
         call pgbox('bcnst',0.,0,'bcnst',0.,0)
      else
         call pgbox('bcnstl',0.,0,'bcnst',0.,0)
      endif
      call pgmtxt('L',2.1,0.5,0.5,"P(k)/P\Dsmooth")
      call pgmtxt('B',2.2,0.5,0.5,"Wavenumber")
      call pgsci(ic)
      call pgslw(2)
      call pgline(ns,xps2,y2)
      if(iall.eq.1) then
         ns1=ns
         do i=1,ns1
            x1ps(i)=xps2(i)
            y1ps(i)=y2(i)
         enddo
      endif
      call pgsci(1)

      nall(nf)=ns
      do i=1,ns
         xpall(i,nf)=xps2(i)
         ypall(i,nf)=y2(i)
      enddo

      enddo
 667  continue
      close(2)

      nt=500
      do i=1,nt
         xt(i)=xmin+(xmax-xmin)*float(i-1)/float(nt-1)
c         call xlinint(xt(i),ns1,x1ps,y1ps,y1v)
c     call xlinint(xt(i),ns,xps2,y2,y2v)
         sum=0.
         do j=1,nf
            do k=1,nall(j)
               xin(k)=xpall(k,j)
               yin(k)=ypall(k,j)
            enddo
            call xlinint(xt(i),nall(j),xin,yin,y2v)
            sum=sum+y2v
         enddo
c         yt(i)=(y1v+y2v)/2.
         yt(i)=sum/float(nf)
      enddo
      call pgslw(7)
      call pgsci(1)
      call pgline(nt,xt,yt)
      
      call pgend

      end

      subroutine smooth(n,x,y,n2,x2,y2,y3)
      parameter(nmax=20000,mm=2,nwk=nmax+6*(nmax*mm+1),mm2=mm*2)
      real x(n),y(n),x2(n2),y2(n2),y3(n)
      real*8 dx(nmax),dy(nmax),wx(nmax),cf(nmax),wk(nwk),val,splder
      real*8 q(mm2),wy(nmax)

      if(n.gt.nmax) print *,'make nmax bigger in smooth'

      read *,val
c      val=0.d0
      md=3
      if(val.eq.0.) md=2
      m=2

      do i=1,n
         dx(i)=dble(x(i))
         dy(i)=dble(y(i))
         wx(i)=1.d0
         wy(i)=1.d0
      enddo

      call gcvspl(dx,dy,nmax,wx,wy,m,n,1,md,val,cf,nmax,wk,ier)
c      call gcvspl(dx,dy,nmax,wx,1.d0,m,n,1,md,val,cf,nmax,wk,ier)
      if(ier.ne.0) print *,'ier= ',ier

      do i=1,n2
         in=i
         y2(i)=sngl(splder(0,m,n,dble(x2(i)),dx,cf,in,q))
      enddo

      do i=1,n
         in=i
         y3(i)=sngl(splder(0,m,n,dble(x(i)),dx,cf,in,q))
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
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
