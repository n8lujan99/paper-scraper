
      parameter(nmax=4000000)
      real r(nmax),d(nmax),b(1000),xi(1000),bs(1000)
      real rs(nmax),ds(nmax)
      real*8 corr(1000)
      integer*8 na,ncorr(1000)
      character cout*3,fout*7
      parameter(pi=3.141593e0,radtodeg=57.29578)

      read *,i1,i2,cout
      fout="xxx.out"
      write(fout(1:3),1001) cout
 1001 format(a3)

      nbin=51
      b1=0.
      b2=3.
      do i=1,nbin
         b(i)=b1+(b2-b1)*float(i-1)/float(nbin-1)
         bs(i)=b(i)*b(i)
         xi(i)=float(i)
         ncorr(i)=0
         corr(i)=0.d0
      enddo

      open(unit=1,file='c1a',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3
         n=n+1
         r(n)=x1
         d(n)=x2
      enddo
 666  continue
      close(1)
      print *,n

      open(unit=1,file='sim3',status='old')

      nsim=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3
         nsim=nsim+1
         rs(nsim)=x1
         ds(nsim)=x2
      enddo
 667  continue
      close(1)
      print *,nsim

c     do i=1,n
      i2=min(i2,nsim-1)
      ic=0
      do i=i1,i2
         ic=ic+1
         if(ic.eq.1000.or.ic.eq.10000.or.ic.eq.100000) then
            do ib=1,nbin
c               print *,i,ncorr(ib),b(ib)
               print *,i,sngl(corr(ib)),b(ib)
            enddo
         endif

         r1=rs(i)
         d1=ds(i)
         jstart=1
         do j=jstart,n
            cdec=cos((d1+d(j))/2./radtodeg)
            rcheck=(r(j)-r1)*cdec
            r2=rcheck*rcheck
            diff=r2+(d1-d(j))**2
            if(diff.le.bs(nbin)) then
               call xlinint(diff,nbin,bs,xi,x0)
               ibin=nint(x0)
c               ncorr(ibin)=ncorr(ibin)+1
               corr(ibin)=corr(ibin)+1.d0
            endif
         enddo
      enddo

      open(unit=11,file=fout,status='unknown')
      do i=1,nbin
c         write(11,*) b(i),ncorr(i)
         write(11,*) b(i),sngl(corr(i))
      enddo
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
      
