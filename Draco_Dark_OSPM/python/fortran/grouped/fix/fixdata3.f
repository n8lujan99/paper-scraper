
      parameter(nmax=16000000)
      real wave(nmax),rat(nmax),wc(nmax),rc(nmax)
      character a5*11,a9*10
      
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

      open(unit=1,file='speccor.dat',status='old')
      ncor=0
      do i=1,nmax
         read(1,*,end=566) x1,x2
         ncor=ncor+1
         wc(ncor)=x1
         rc(ncor)=x2
      enddo
 566  continue
      close(1)

      iseed=-1

      open(unit=1,file='data.in',status='old')
      open(unit=11,file='data1_out',status='unknown')
      open(unit=12,file='data2_out',status='unknown')
      nsim=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3,x4,a5,x6,x7,x8,a9
         call xlinint(x3,n,wave,rat,r0)
         call xlinint(x3,ncor,wc,rc,r0b)
         sn=x4*r0b
         if(r0.lt.1.1) then
            write(11,1101) x1,x2,x3,x4,a5,x6,x7,x8,a9
         else
            ran=ran2(iseed)
            if(ran.lt.1./r0) write(11,1101) x1,x2,x3,x4,a5,x6,x7,x8,a9
         endif
         if(sn.ge.4.8) write(12,1101) x1,x2,x3,sn,a5,x6,x7,x8,a9
      enddo
 667  continue
      close(1)
      close(11)
      close(12)

 1101 format(2(f9.5,1x),f8.3,1x,f6.2,1x,a11,1x,f6.2,1x,f7.3
     $     ,1x,f6.3,1x,a10)
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
