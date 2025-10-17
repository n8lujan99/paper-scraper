
      parameter(nmax=16000000)
      real wave(nmax),sn(nmax)
      
      open(unit=1,file='data.use',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4
         n=n+1
         wave(n)=x3
         sn(n)=x4
      enddo
 666  continue
      close(1)
      print *,n

      iseed=-1

      open(unit=1,file='simm_all',status='old')
      open(unit=11,file='sim_out',status='unknown')
      nsim=0
      do i=1,nmax
         read(1,*,end=667) x1,x2
         iw=nint(ran2(iseed)*n)
         iw=max(iw,1)
         iw=min(iw,n)
         waveo=wave(iw)
         write(11,*) x1,x2,waveo,sn(iw)
      enddo
 667  continue
      close(1)
      close(11)

      end
