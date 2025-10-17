
      parameter(nmax=10000)
      real w(nmax),t1(nmax),t2(nmax),t3(nmax),t4(nmax)
      real t5(nmax),we(nmax),fe(nmax)

      open(unit=1,file='in',status='old')
      n=0
      do i=1,1000
c         read(1,*,end=666) x1,x2,x3,x4,x5,x6
         read(1,*,end=666) x1,x2,x3,x4,x5
         n=n+1
         w(n)=x1
         t1(n)=x2
         t2(n)=x3
         t3(n)=x4
         t4(n)=x5
c         t5(n)=x6
      enddo
 666  continue
      close(1)

c      open(unit=1,file='/scratch/00115/gebhardt/detect/extcor.dat',
      open(unit=1,file='../../extcor.dat',
     $     status='old')
      ne=0
      do i=1,1000
         read(1,*,end=667) x1,x2
         ne=ne+1
         we(ne)=x1
         fe(ne)=x2
      enddo
 667  continue
      close(1)

      open(unit=11,file='out',status='unknown')
      do i=1,n
         call xlinint(w(i),ne,we,fe,fcor)
c         write(11,1101) w(i),t1(i)*fcor,t2(i),t3(i),t4(i),t5(i)
         write(11,1101) w(i),t1(i)*fcor,t2(i)*fcor,t3(i)*fcor,t4(i)*fcor
      enddo
      close(11)

 1101 format(1x,f6.1,4(2x,f5.3))
c 1101 format(1x,f6.1,5(2x,f5.3))
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
