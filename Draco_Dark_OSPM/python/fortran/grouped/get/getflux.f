
      real xslf(10000),yslf(10000),x(10000),y(10000)
      real*8 pstar,xlstar,alpha,dx,dy
      
c      call pgbegin(0,'?',1,1)
c      call pgpap(0.,1.)
c      call pgsch(1.5)
c      call pgscf(2)
c      call pgslw(2)

      xmin=41.7
      xmin=41.5
      xmax=44.
      ymin=log10(60.)
      ymax=log10(2.e6)

c- from Ouchi
      pstar=3.9d-4
      xlstar=0.849d43
      alpha=-1.8d0
      vol=1.e9
      nslf=1000
      do i=1,nslf
         xslf(i)=xmin+(xmax-xmin)*float(i-1)/float(nslf-1)
         dx=dble(xslf(i))
         dx=10.d0**dx  
         dy=pstar*((dx/xlstar)**(alpha+1)) * dexp(-(dx/xlstar))
     $        *log(10.)
         yslf(i)=sngl(dlog10(dy))+log10(vol)
         x(i)=10**(xslf(i)-43.)
c         y(i)=10**yslf(i)
         y(i)=1./x(i)*10**yslf(i)
c         print *,x(i)*19.1,y(i)
      enddo

      xmin=x(1)
      xmax=x(nslf)
      ymin=0.
      ymax=y(1)
      iseed=-1
      do i=1,20000000
         xran=xmin+(xmax-xmin)*ran2(iseed)
         yran=ymin+(ymax-ymin)*ran2(iseed)
         call xlinint(xran,nslf,x,y,ynew)
         xran=xran*19.1
         if(yran.le.ynew) print *,xran,yran
      enddo
c 1e43/5.226e58/1e-17 = 19.1
c      call pgenv(xmin,xmax,ymin,ymax,0,20)
c      call pglabel("log\D10\UL\DLy\ga","dN/dlogL","")
c      call pgline(nslf,xslf,yslf)

c      call pgend

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
      
