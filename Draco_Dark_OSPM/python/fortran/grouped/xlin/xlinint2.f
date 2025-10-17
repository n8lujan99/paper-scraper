      subroutine xlinint2(xp,n,x,y,yp,jin,jout)
      real x(n),y(n)
      do j=jin,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            jout=j
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      jout=1
      return
      end
