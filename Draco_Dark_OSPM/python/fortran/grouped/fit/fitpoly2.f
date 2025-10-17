
      parameter(nmax=1000)
      real x(nmax),y(nmax),xp(nmax),yp(nmax),c(5)
      real c1a(nmax),c2a(nmax),c3a(nmax),c4a(nmax)
      real c0a(nmax),rbc(6),ynew(nmax),xin(nmax),wg(nmax),tg(nmax)
      real xin2(nmax),yin2(nmax),xa(nmax),ya(nmax)
      character title*12

      open(unit=1,file='greg.tp',status='old')
      ng=0
      do i=1,nmax
         read(1,*,end=888) x1,x2
         ng=ng+1
         wg(ng)=x1
         tg(ng)=x2/1.515
      enddo
 888  continue
      close(1)
         
c      ifit0=1 ! 5th order poly
      ifit0=2 ! Robin's fit

c - Robin's coefficients
      rbc(1)=-0.01972573
      rbc(2)=+0.43904144
      rbc(3)=-3.86654830
      rbc(4)=+16.7814045
      rbc(5)=-35.7178268
      rbc(6)=+29.7167950

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=3500.
      xmax=5500.
      ymin=0.95
      ymax=1.05
      call pgsls(1)
      call pgslw(1)

      call pgsci(1)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength','kg/gz','')
c      call pglabel('Wavelength','data/fit',title)

      np=100
      do i=1,np
         xp(i)=xmin+(xmax-xmin)/float(np-1)*float(i-1)
         xrbc=xp(i)/1000.
         fit=0.
         do ifit=1,6
            fit=(fit*xrbc)+rbc(ifit)
         enddo
         ynew(i)=fit
         call xlinint(xp(i),ng,wg,tg,y0)
         ynew(i)=fit/y0
      enddo

      call pgline(np,xp,ynew)
c      call pgsci(2)
c      call pgline(ng,wg,tg)

      call pgend
 1101 format(1x,f6.1,5(1x,f7.4))

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
      if(xp.ge.x(n)) yp=y(n)
      return
      end
