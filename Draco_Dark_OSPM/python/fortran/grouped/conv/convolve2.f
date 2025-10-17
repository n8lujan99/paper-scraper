
      parameter (nmax=50000,pi=3.141592e0)
      real r(nmax),s(nmax),yp(nmax)
      integer isum(nmax)
      character file1*80

      file1="j1"
      rsig=0.2

      h4=0.

      open(unit=1,file=file1,status='old')
      n=0
      sum=0.
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         r(n)=x1
         s(n)=x2
         sum=sum+s(n)
      enddo
 666  continue
      close(1)
      call sort2(n,r,s)

      open(unit=1,file='convolve.out',status='unknown')
      sum2=0.
      do i=1,n
         sig=r(i)*rsig
         den2=2.*sig*sig
         xp=r(i)
         yp(i)=0.
         isum(i)=0
         rlo=max(0.,xp-5.*sig)
         rup=xp+5.*sig
         rstep=(rup-rlo)/1000.
         do rnew=rlo,rup,rstep
            w=(xp-rnew)/sig
            gaus=exp(-w*w/2.)/sqrt(den2*pi)*rstep
            call xlinint(rnew,n,r,s,sg)
            yp(i)=yp(i)+sg*gaus*(1.+h4*fh4(w))
            isum(i)=isum(i)+1
         enddo
         sum2=sum2+yp(i)
      enddo
      do i=1,n
         write(1,*) r(i),yp(i)/sum2*sum
      enddo

      end

      function fh4(x)
      fh4=1./sqrt(24.)*(4.*x*x*x*x-12.*x*x+3.)
      return
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end
