
      parameter(nmax=16000000)
      real r(nmax),d(nmax),wave(nmax)
      
      open(unit=1,file='data.use',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3
         n=n+1
         r(n)=x1
         d(n)=x2
         wave(n)=x3
      enddo
 666  continue
      close(1)
      print *,n

      iseed=-1

      nsim=3000000
      open(unit=11,file='sim_out',status='unknown')
      do i=1,nsim
         ip=nint(ran2(iseed)*n)
         ip=max(ip,1)
         ip=min(ip,n)
         iw=nint(ran2(iseed)*n)
         iw=max(iw,1)
         iw=min(iw,n)
         write(11,*) r(ip),d(ip),wave(iw)
      enddo
      close(11)

      end
