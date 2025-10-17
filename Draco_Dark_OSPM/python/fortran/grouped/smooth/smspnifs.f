
      parameter(nmax=100000000,nmax2=100000)
      real x(nmax),y(nmax),xs(nmax2),ys(nmax2),y3(nmax2)
      real yn(nmax2),yb(nmax2),xn(nmax2),xb(nmax2)

      ibin=5

      open(unit=1,file='in',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         x(n)=x1
         y(n)=x2
      enddo
 666  continue
      close(1)

      open(unit=11,file='out',status='unknown')

      ib1=(ibin-1)/2
      xib=float(ibin)
      nbb=0
      do j=1,n,ibin
         nbb=nbb+1
         istart=max(1,j-ib1)
         iend=istart+ibin-1
         if(iend.gt.n) then
            iend=n
            istart=n-ibin+1
         endif
         sum=0.
         nb=0
         do is=istart,iend
            sum=sum+y(is)
            nb=nb+1
            yb(nb)=y(is)
            xb(nb)=x(is)
         enddo
         sum=sum/float(nb)
         call biwgt(yb,nb,xbb,xsb)
         yn(nbb)=xbb
         call biwgt(xb,nb,xbb,xsb)
         xn(nbb)=xbb
      enddo

      do i=1,n
         call xlinint(x(i),nbb,xn,yn,y0)
         if(y(i).gt.0.) then
c            rat=y(i)/y0
            rat=(y(i)-y0)/y(i)
         else
            rat=0.
         endif
         write(11,*) x(i),y(i),rat,y0
      enddo
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

