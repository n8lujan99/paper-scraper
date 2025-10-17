
      parameter(nmax=18000000)
      real x(nmax),wave(nmax),xp(1000),yp(1000),xn(nmax),sx(10)
      real wz(10),cz(10),yp2(nmax),xs(1000),ys(1000)
      character file1*80

      open(unit=1,file='smline.out',status='old')
      ns=0
      do i=1,nmax
         read(1,*,end=555) x1,x2
         ns=ns+1
         xs(ns)=x1
         ys(ns)=x2
      enddo
 555  continue
      close(1)

      file1="z1"
      open(unit=1,file=file1,status='old')

      nz=5
      wz(1)=3550.
      wz(2)=3800.
      wz(3)=4100.
      wz(4)=4500.
      wz(5)=4800.
      cz(1)=0.87
      cz(2)=0.87
      cz(3)=0.92
      cz(4)=0.98
      cz(5)=0.98

      n=0
      xmin=3520.
      xmax=5400.
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,x5,x6
c         if(x3.gt.3.0) x2=x2*0.8
c     if(x1.lt.4000.) x2=x2*0.95
c         call xlinint(x1,nz,wz,cz,cz0)
c         x2=x2*cz0
         x4=x4*cz0
         if(x2.gt.4.8.and.x6.gt.40.) then
            n=n+1
            wave(n)=x1
         endif
      enddo
 666  continue
      close(1)
      print *,n

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      call pgenv(xmin,xmax,0.,1.05,0,0)
      call pglabel('Wavelength','Normalized Counts','')
      call pgslw(6)
      
      nbin=50
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
      enddo
      call pgline(nbin,xp,yp)

      ic=1
      open(unit=1,file='zlist',status='old')
      do iall=1,1000
         read(1,*,end=667) file1,ymax2
         open(unit=2,file=file1,status='old')
         n=0
         do j=1,nmax
            read(2,*,end=668) x1,x2,x3
            if(x3.lt.80.0) then
               n=n+1
               wave(n)=x1
            endif
         enddo
 668     continue
         close(2)
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
            yp2(i)=float(nb)
            sum=sum+yp2(i)
            ymax=max(ymax,yp2(i))
         enddo
         sum=ymax
         sumt2=0.
         do i=1,nbin
            yp2(i)=yp2(i)/sum
            sumt2=sumt2+yp2(i)
         enddo
         print *,sumt,sumt2,sumt2/sumt
         ic=ic+1
         if(ic.eq.3) ic=ic+1
         call pgsci(ic)
c         call pgline(nbin,xp,yp2)
      enddo
 667  continue
      close(1)

      do i=1,nbin
         call xlinint(xp(i),ns,xs,ys,ys0)
         yp2(i)=yp2(i)*ys0
c         print *,i,xp(i),yp(i)/yp2(i)
      enddo
      call pgsci(2)
      call pgline(nbin,xp,yp2)

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
