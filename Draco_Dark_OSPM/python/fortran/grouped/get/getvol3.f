
      parameter(nmax=10000)
      real f(nmax),c(nmax),z(nmax),fc(nmax)
      real zl(nmax),zfc(nmax),xlc(nmax),xcc(nmax)
      real*8 dx2

      flim=0.018
      flim=0.05
      
      open(unit=1,file='sn4.8.ccor',status='old')
      read(1,*)
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         f(n)=x1
         c(n)=x2
      enddo
 666  continue
      close(1)

      call xlinint(flim,n,c,f,f02)

      open(unit=1,file='FtoLz2.dat',status='old')
      n2=0
      do i=1,nmax
         read(1,*,end=667) x1,dx2
         n2=n2+1
         z(n2)=x1
         fc(n2)=sngl(dx2/1.d42)*1.e-17
      enddo
 667  continue
      close(1)

      zmin=1.88
      zmax=3.52
      nz=100
      do j=1,nz
         zl(j)=zmin+(zmax-zmin)*float(j-1)/float(nz-1)
         call xlinint(zl(j),n2,z,fc,fc0)
         zfc(j)=fc0
      enddo

      call xlinint(1.9,n2,z,fc,fc0)
      call xlinint(3.5,n2,z,fc,fc1)
      fc0=fc0*f02
      fc1=fc1*f02
c      print *,f02
c      print *,fc0,fc1,fc1/fc0

      xmin=0.
      xmax=1.5
      nl=200
      n3=0
      do i=1,nl
         xl=xmin+(xmax-xmin)*float(i-1)/float(nl-1)
         x=10**xl
         nt=0
         do j=1,nz
            ftry=x/zfc(j)
            if(ftry.gt.f02) nt=nt+1
         enddo
         if(nt.gt.0) then
            n3=n3+1
            xlc(n3)=xl+42.
            xcc(n3)=float(nz)/float(nt)
         endif
      enddo

      open(unit=1,file='out',status='old')
      open(unit=11,file='out2',status='unknown')
      do i=1,nmax
         read(1,*,end=668) x1,x2,x3,x4
         call xlinint(x1,n3,xlc,xcc,xc0)
         xc0=min(5.,xc0)
         write(11,*) x1,x2*xc0,x2,x3,x4,xc0
      enddo
 668  continue
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

      
