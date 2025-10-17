Skipped
      parameter(nmax=10000)
      real zl(nmax),xl(nmax),xlum(nmax),xlum0(nmax)
      real*8 dx2
      integer n1(nmax),n2(nmax),n3(nmax),n4(nmax),n5(nmax),n6(nmax)
      character cfile*50

      xly=1215.666

      nlum=9
      xlum(1)=42.3
      xlum(2)=42.4
      xlum(3)=42.6
      xlum(4)=42.8
      xlum(5)=43.0
      xlum(6)=43.2
      xlum(7)=43.4
      xlum(8)=43.6
      xlum(9)=43.8
      do i=1,nlum
         n1(i)=0
         n2(i)=0
         n3(i)=0
         n4(i)=0
         n5(i)=0
         n6(i)=0
         xlum0(i)=xlum(i)
         xlum(i)=10**(xlum(i)-40.)
      enddo

c- get the flux to luminosity as a function of z
      open(unit=1,file='FtoLz2.dat',status='old')
      nz=0
      do i=1,10000
         read(1,*,end=668) x1,dx2
         nz=nz+1
         zl(nz)=x1
         xl(nz)=sngl(dx2/1.d40)
      enddo
 668  continue
      close(1)

      read *,cfile
      open(unit=1,file=cfile,status='old')

      fac=0.8
      fac=0.65
      fac=0.45
      do i=1,10000000
         read(1,*,end=666) x1,x2,x3,x4,x5,x6
c         x4p8=4.8*x4/x5*fac
c         x7p0=7.0*x4/x5*fac
         x4p8=4.8*x4*fac
         x5p0=5.0*x4*fac
         x5p5=5.5*x4*fac
         x6p0=6.0*x4*fac
         x6p5=6.5*x4*fac
         x7p0=7.0*x4*fac
         if(x3.gt.3510..and.x3.lt.5490.) then
            if(x4.gt.0.and.x5.gt.0.5) then
               z=x3/xly-1.
               call xlinint(z,nz,zl,xl,x0)
               do il=1,nlum
                  x=xlum(il)/x0/1.e-17
                  if(x4p8.lt.x) n1(il)=n1(il)+1
                  if(x5p0.lt.x) n2(il)=n2(il)+1
                  if(x5p5.lt.x) n3(il)=n3(il)+1
                  if(x6p0.lt.x) n4(il)=n4(il)+1
                  if(x6p5.lt.x) n5(il)=n5(il)+1
                  if(x7p0.lt.x) n6(il)=n6(il)+1
               enddo
            endif
         endif
      enddo
 666  continue
      close(1)

      open(unit=11,file='out',status='unknown')
      do i=1,nlum
         rat=0.
         if(n1(i).gt.0.) rat=float(n2(i))/float(n1(i))
         write(11,1001) xlum0(i),n1(i),n2(i),n3(i),n4(i),n5(i),n6(i)
      enddo
      close(11)

 1001 format(f5.2,6(2x,i8))
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

