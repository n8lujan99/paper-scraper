
      real x(1000),y(1000)

      read *,snc
      open(unit=1,file='vol.dat',status='old')
      read(1,*)
      n=0
      do i=1,1000
         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7
         if(snc.eq.4.8) xuse=x2
         if(snc.eq.5.0) xuse=x3
         if(snc.eq.5.5) xuse=x4
         if(snc.eq.6.0) xuse=x5
         if(snc.eq.6.5) xuse=x6
         if(snc.eq.7.0) xuse=x7
         n=n+1
         x(n)=x1
         y(n)=xuse
      enddo
 666  continue
      close(1)
      
      open(unit=1,file='in_vol',status='old')
      open(unit=11,file='out',status='unknown')
      do i=1,1000
         read(1,*,end=667) x1,x2,x3,x4
         call xlinint(x1,n,x,y,y0)
         write(11,*) x1,x2*y0,x3,x4
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


