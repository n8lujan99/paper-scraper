
      parameter(nmax=20000)
      real x(nmax),y(nmax),xs(nmax),ys(nmax),y3(nmax)
      real xsm(nmax),ysm(nmax),ynew(nmax)

c      print *,"xmin,xmax,ymin,ymax: "
c      read *,xmin,xmax,ymin,ymax

      open(unit=1,file='tmp2',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         x(n)=x1
         y(n)=x2
      enddo
 666  continue
      close(1)

      open(unit=1,file='smline.out',status='old')
      ns=0
      do i=1,nmax
         read(1,*,end=667) x1,x2
         ns=ns+1
         xs(ns)=x1
         ys(ns)=x2
      enddo
 667  continue
      close(1)

      open(unit=11,file='out',status='unknown')
      do i=1,n
         call xlinint(x(i),ns,xs,ys,y0)
         write(11,*) x(i),y(i),y0
      enddo
      close(11)
      
      end

      subroutine smooth(n,x,y,n2,x2,y2,y3)
      parameter(nmax=20000,mm=2,nwk=nmax+6*(nmax*mm+1),mm2=mm*2)
      real x(n),y(n),x2(n2),y2(n2),y3(n)
      real*8 dx(nmax),dy(nmax),wx(nmax),cf(nmax),wk(nwk),val,splder
      real*8 q(mm2),wy(nmax)

      if(n.gt.nmax) print *,'make nmax bigger in smooth'

      call qd1('Enter smoothing val ','smflat.def',val)
      call savdef
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
