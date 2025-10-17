
      real x(1000),xn(1000),y(1000)
      
      read *,f50
      read *,f0
      open(unit=1,file='cc_norm.txt',status='old')
      open(unit=11,file='outfix',status='unknown')
      write(11,*) 0.5,f50

      n=0
      do i=1,100
         read(1,*,end=666) x1,x2,x3
         n=n+1
         x(n)=x1
         xn(n)=x(n)*f50/f0
         y(n)=x2
      enddo
 666  continue
      close(1)

      do i=1,100
         x0=1.+0.5*float(i-1)
         call xlinint(x0,100,xn,y,y0)
         write(11,*) x0,y0,0.01
      enddo

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

      
