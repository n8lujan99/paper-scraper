
      parameter(nmax=1000)
      real r1(nmax),s1(nmax),r2(nmax),s2(nmax)
      real r(nmax),s(nmax),rp(10000),sp(10000)
      parameter(pi=3.141593e0)

      rfwhm=0.14
c      rfwhm=0.22
      rsig=rfwhm/2.35

      open(unit=1,file='phot.dat',status='old')
      n1=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,x5,x6
         if(x1.lt.150.) then
            n1=n1+1
            r1(n1)=log10(x1)
            s1(n1)=22.5-2.5*log10(x3+x6)
            r(n1)=x1
            s(n1)=x3+x6
         endif
      enddo
 666  continue
      close(1)

      open(unit=1,file='N4826Ks.pr',status='old')
      do i=1,6
         read(1,*)
      enddo
      n2=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3,x4,x5
         n2=n2+1
         r2(n2)=log10(x2)
         s2(n2)=x4
      enddo
 667  continue
      close(1)

c- convolve
      xs=r(1)
      xe=r(n1)
      n=1000
      xdiff=4.*rfwhm
      nr=1000
      do i=1,n
         rp(i)=xs+(xe-xs)/float(n-1)*float(i-1)
         xmin=rp(i)-xdiff
         xmax=rp(i)+xdiff
         gaus=0.
         sum=0.
         do j=1,nr
            xt=xmin+(xmax-xmin)/float(nr-1)*float(j-1)
            call xlinint(abs(xt),n1,r,s,st)
            if(xt.ge.0) then
               w=(xt-rp(i))/rsig
c               g0=exp(-w*w/2.)
               g0=exp(-w*w/2.)*xt*xt
               gaus=gaus+st*g0
               sum=sum+g0
            endif
         enddo
         sp(i)=gaus/sum
      enddo
      rat=1.0
      open(unit=11,file='out',status='unknown')
      do i=1,n
         rout=rp(i)
         call xlinint(rout,n1,r,s,st)
         rp(i)=log10(rp(i))
         sp(i)=22.5-2.5*log10(sp(i)*rat)
         st=22.5-2.5*log10(st)
         write(11,1101) rout,sp(i),st
c         print *,i,rp(i),sp(i)
      enddo         
      close(11)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

c      xmin=log10(0.05)
c      xmax=log10(30.)
      xmin=log10(0.25)
      xmax=log10(100.)
      ymin=17.
      ymax=9.5
      ymax=10.0

      call pgenv(xmin,xmax,ymin,ymax,0,10)
      call pglabel("R(arcsec)","SB_K","")

      call pgline(n2,r2,s2)
      call pgsci(2)
      call pgslw(6)
      call pgline(n1,r1,s1)
c      call pgsci(4)
c      call pgline(n,rp,sp)

      call pgend
      
 1101 format(f7.3,2(1x,f7.3))
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
