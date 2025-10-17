
      parameter(nmax=1000)
      real xf(nmax),yf(nmax)
      parameter(cee=2.99e+5,pi=3.141593,ee=2.71828)

      vred=16700.
      
      open(unit=1,file='fixspec.dat',status='old')
      nf=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         nf=nf+1
         xf(nf)=x1
         yf(nf)=x2
      enddo
 666  continue
      close(1)

      open(unit=1,file='in',status='old')
      open(unit=11,file='out',status='unknown')
      do i=1,10000
         read(1,*,end=667) x1,x2
         x1in=x1/(1.+vred/cee)
         call xlinint(x1in,nf,xf,yf,y0)
         write(11,*) x1,x2*y0
      enddo
 667  continue
      close(1)
      close(11)

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
