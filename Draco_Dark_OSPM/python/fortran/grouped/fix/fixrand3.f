
      parameter(nmax=16000000)
      real wave(nmax),rat(nmax)
      character a5*3,a6*8,a7*3,a9*13,amp*7
      
      open(unit=1,file='ratrand',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         wave(n)=x1
         rat(n)=x2
      enddo
 666  continue
      close(1)

      iseed=-1

c      open(unit=1,file='simm_all',status='old')
      open(unit=1,file='simm0_all',status='old')
      open(unit=11,file='sim_out',status='unknown')
      nsim=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3,x4,a5,a6,a7,x8,a9
         amp=a9(7:13)
         if(amp.eq."309_046".or.amp.eq."509_094".or.amp.eq."305_034")
     $        then
            ran=ran2(iseed)
            if(ran.lt.0.5) goto 668
         endif
         call xlinint(x3,n,wave,rat,r0)
         if(r0.gt.0.9) then
            write(11,1101) x1,x2,x3,x4,a5,a6,a7,x8,a9
         else
            ran=ran2(iseed)
            if(ran.lt.r0) write(11,1101) x1,x2,x3,x4,a5,a6,a7,x8,a9
         endif
 668     continue
      enddo
 667  continue
      close(1)
      close(11)

 1101 format(2(f9.5,1x),f8.3,1x,f6.2,1x,a3,1x,a8,1x,a3,1x,f8.3,1x,a13)
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
