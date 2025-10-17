
      parameter(nmax=16000000)
      real wave(nmax),sn(nmax),ws(nmax),ss(nmax)
      
c      open(unit=1,file='data.use',status='old')
c      open(unit=1,file='data.bad',status='old')
      open(unit=1,file='skyn.dat',status='old')
      n=0
      do i=1,nmax
c         read(1,*,end=666) x1,x2,x3,x4
         read(1,*,end=666) x1,x2
         n=n+1
c         wave(n)=x3
c         sn(n)=x4
         ws(n)=x1
         ss(n)=x2/200.
      enddo
 666  continue
      close(1)
      print *,n

      iseed=-1

      open(unit=1,file='ifupos.all',status='old')
      open(unit=11,file='ifu_out',status='unknown')
      nsim=0
      do i=1,nmax
         read(1,*,end=667) x1,x2
c         iw=nint(ran2(iseed)*n)
c         iw=max(iw,1)
c         iw=min(iw,n)
c         waveo=wave(iw)
c     write(11,*) x1,x2,waveo,sn(iw)
         do j=1,100
            iw=nint(ran2(iseed)*n)
            iw=max(iw,1)
            iw=min(iw,n)
            waveo=ws(iw)
            rs=ran2(iseed)
c            if(rs.lt.ss(iw)) then
            if(rs.ge.ss(iw)) then
               write(11,*) x1,x2,waveo,5.0
               goto 888
            endif
         enddo
 888     continue
c         waveo=3500.+ran2(iseed)*2000.
c         write(11,*) x1,x2,waveo,5.0
      enddo
 667  continue
      close(1)
      close(11)

      end
