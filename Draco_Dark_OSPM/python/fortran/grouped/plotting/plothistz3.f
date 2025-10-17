
      parameter(nmax=18000000)
      real x(nmax),wave(nmax),xp(1000),yp(1000),xn(nmax),sx(10)
      real wz(10),cz(10),yp2(nmax),xs(nmax),ys(nmax),xslo(10),xshi(10)
      real xwlo(10),xwhi(10),xsky(10000),ysky(10000),yskyp(10000)
      real xdata(nmax),ydata(nmax),xpo(1000),ypo(1000),xin(1000)
      real xcor(2000),ycor(2000),ypc(nmax),sna(100000)
      character file1*80

      xmin=3500.
      xmax=3800.
      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      open(unit=11,file='hout',status='unknown')

      open(unit=1,file='skyn.dat',status='old')
      nsky=0
      do i=1,nmax
         read(1,*,end=556) x1,x2
         nsky=nsky+1
         xsky(nsky)=x1
         ysky(nsky)=x2
      enddo
 556  continue
      close(1)

      open(unit=1,file='speccor.in',status='old')
      ncor=0
      do i=1,nmax
         read(1,*,end=456) x1,x2
         ncor=ncor+1
         xcor(ncor)=x1
c         ycor(ncor)=x2
         ycor(ncor)=1.
      enddo
 456  continue
      close(1)

      open(unit=1,file='data.in',status='old')
c      open(unit=1,file='data_out',status='old')
      ndata=0
      do i=1,nmax
         read(1,*,end=557) x1,x2,x3,x4
         call xlinint(x3,ncor,xcor,ycor,y0)
         sn=x4*y0
         if(sn.ge.4.8) then
            ndata=ndata+1
            xdata(ndata)=x3
            ydata(ndata)=x4
         endif
      enddo
 557  continue
      close(1)
      print *,ndata

c      open(unit=1,file='simm_all',status='old')
      open(unit=1,file='sim_out',status='old')
      ns=0
      do i=1,nmax
         read(1,*,end=555) x1,x2,x3,x4
         ns=ns+1
         xs(ns)=x3
         ys(ns)=x4
      enddo
 555  continue
      close(1)
      print *,ns

      xslo(1)=4.8
      xshi(1)=5.2
      xslo(2)=6.7
      xshi(2)=7.5

      xwlo(1)=3500.
      xwhi(1)=3800.
      xwlo(2)=3800.
      xwhi(2)=4100.
      xwlo(3)=4100.
      xwhi(3)=4400.
      xwlo(4)=4400.
      xwhi(4)=4700.
      xwlo(5)=4700.
      xwhi(5)=5000.
      xwlo(6)=5000.
      xwhi(6)=5300.
      xwlo(7)=5300.
      xwhi(7)=5500.

      do iwave=1,7
         xmin=xwlo(iwave)
         xmax=xwhi(iwave)

      call pgsci(1)
      call pgslw(2)
c      call pgenv(xmin,xmax,0.3,1.25,0,0)
      call pgenv(xmin,xmax,0.3,2.5,0,0)
      call pglabel('Wavelength','Normalized Counts','')
      call pgslw(6)

      do iall=1,1
         call pgsci(iall)
      
      n=0
      do i=1,ns
         if(ys(i).gt.xslo(iall).and.ys(i).lt.xshi(iall)) then
            n=n+1
            wave(n)=xs(i)
         endif
      enddo
 666  continue
      close(1)
      
      nbin=150
      xbin=(xmax-xmin)/float(nbin-1)
      ymax=0.
      sum=0.
      do i=1,nbin
         xlo=xmin+xbin*float(i-1)
         xhi=xlo+xbin
         nb=0
         do j=1,n
            if(wave(j).ge.xlo.and.wave(j).lt.xhi) nb=nb+1
         enddo
         xp(i)=(xhi+xlo)/2.
         yp(i)=float(nb)
         sum=sum+yp(i)
         ymax=max(ymax,yp(i))
      enddo
      sum=ymax
      sumt=0.
      do i=1,nbin
         yp(i)=yp(i)/sum
         sumt=sumt+yp(i)
         xpo(i)=xp(i)
         ypo(i)=yp(i)
      enddo
      call pgline(nbin,xp,yp)

      nbin=150
      xbin=(xmax-xmin)/float(nbin-1)
      ymax=0.
      sum=0.
      do i=1,nbin
         xlo=xmin+xbin*float(i-1)
         xhi=xlo+xbin
         nb=0
         do j=1,ndata
            if(xdata(j).ge.xlo.and.xdata(j).lt.xhi) then
               nb=nb+1
            endif
         enddo
         xp(i)=(xhi+xlo)/2.
         yp(i)=float(nb)
         ypc(i)=yp(i)
         sum=sum+yp(i)
         ymax=max(ymax,yp(i))
      enddo
      sum=ymax
      sumd=sum
      nb=0
      do i=1,nbin
         yp(i)=yp(i)/sum
c         write(11,*) xpo(i),xp(i),ypo(i),yp(i)
         if(yp(i).gt.0) then
            nb=nb+1
            xin(nb)=ypo(i)/yp(i)
         endif
      enddo
      call biwgt(xin,nb,xb0,xs0)
      do i=1,nbin
         yp(i)=yp(i)*xb0
         rat=0.
         if(yp(i).gt.0.) rat=ypo(i)/yp(i)
         ncut=nint((1.-rat)*ypc(i))
c         if(rat.lt.1.) print *,xpo(i),ypc(i),ncut,rat

         xlo=xmin+xbin*float(i-1)
         xhi=xlo+xbin
         nsn=0
         do j=1,ndata
            if(xdata(j).ge.xlo.and.xdata(j).lt.xhi) then
               nsn=nsn+1
               sna(nsn)=ydata(j)
            endif
         enddo
         if(nsn.gt.1) call sort(nsn,sna)
         if(rat.ge.0.95) ncut=1
         sncut=sna(ncut)
         if(sncut.lt.4.91) sncut=4.8
         if(ypo(i).lt.0.1) sncut=4.8
         write(11,*) xpo(i),ypo(i),yp(i),sncut,ncut
      enddo
      call pgsci(2)
      call pgline(nbin,xp,yp)

      smax=0.
      do isky=1,nsky
         if(xsky(isky).gt.xmin.and.xsky(isky).lt.xmax) then
            smax=max(smax,ysky(isky))
         endif
      enddo
      do isky=1,nsky
         yskyp(isky)=ysky(isky)/smax
      enddo
      call pgsci(4)
c      call pgline(nsky,xsky,yskyp)

      enddo
      enddo
      close(11)

      call pgend

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
