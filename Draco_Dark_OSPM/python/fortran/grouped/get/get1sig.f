Completed
      parameter(nmax=1036)
      real w(nmax),s(nmax),se(nmax)

      read *,wave
      f1sig=0.

      open(unit=1,file='spec.out',status='old')
      n=0
      do i=1,1036
         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7,x8,x9
         n=n+1
         w(n)=x1
         s(n)=x2
         if(x9.le.0) then
            se(n)=0.
         else
            se(n)=x3/x9
         endif
      enddo
 666  continue
      close(1)
      if(n.lt.10) goto 888

      nflimh=3
      nw=n
      diff=1e10
      do iw=1,nw
         dn=abs(w(iw)-wave)
         if(dn.lt.diff) then
            diff=dn
            iw0=iw
         endif
      enddo
      xnoise2=0.
      ilo=max(1,iw0-nflimh)
      ihi=min(nw,iw0+nflimh)
      do i=ilo,ihi
         xnoise2=xnoise2+se(i)**2
      enddo
      xnoise2=sqrt(xnoise2)
      f1sig=xnoise2
 888  continue
      open(unit=11,file='get1sig.out',status='unknown')
      write(11,*) f1sig
      close(11)

      end
