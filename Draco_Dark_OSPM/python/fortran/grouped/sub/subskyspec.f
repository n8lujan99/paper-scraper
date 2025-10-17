
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),ragn(nmax),agnv(nmax)
      real xb(nmax),yb(nmax),wsky(nmax),snorm(nmax)
      character file1*80,file2*80,c1*18

      nf=18
      ibin=9
c      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      open(unit=1,file='skyres_norm',status='old')
      nsky=0
      do i=1,10000
         read(1,*,end=555) x1,x2
         nsky=nsky+1
         wsky(nsky)=x1
         snorm(nsky)=x2
      enddo
 555  continue
      close(1)

      open(unit=1,file='agn.dat',status='old')
      nagn=0
      do i=1,1000
         read(1,*,end=556) x1,x2
         nagn=nagn+1
         ragn(nagn)=x1
         agnv(nagn)=x2
      enddo
 556  continue
      close(1)

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=3810.
      xmax=3870.
      xmin=3500.
      xmax=5500.
      xmin=22000.
      xmax=24300.
      
c      xmin=4630.
c      xmax=4730.
c      xmin=6500.
c      xmax=6700.

c      xmin=8400.
c      xmax=8900.

      open(unit=1,file='splist',status='old')
      open(unit=11,file='out',status='unknown')

      wave0=22000.
      agn0=1150.
      wave1=24000.
      agn1=1550.
      slope=(agn1-agn0)/(wave1-wave0)

      nl=0
      ic=0
      do il=1,10000
c         read(1,*,end=666) file1
         read(1,*,end=666) file1,xnorm
         open(unit=2,file=file1,status='old')
         ymin=1e10
         ymax=-ymin
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3,x4,x5
c            read(2,*,end=667) x1,x2,x3,x4,x5,x6
c            x2=x6
c            if(x2.ne.0) then
               n=n+1
               x(n)=x1
               y(n)=x2
               radagn=x5
               ymin=min(ymin,y(n))
               ymax=max(ymax,y(n))
c            endif
         enddo
 667     continue
         close(2)
         call xlinint(radagn,nagn,ragn,agnv,a0)
         xnorm=a0
c         print *,a0,xnorm
         nbb=0
         do j=1,n,ibin
            nbb=nbb+1
            istart=max(0,j-ib1)
            iend=istart+ibin-1
            if(iend.gt.n) then
               iend=n
               istart=n-ibin+1
            endif
            sum=0.
            nb=0
            do is=istart,iend
               sum=sum+y(is)
               nb=nb+1
               yb(nb)=y(is)
               xb(nb)=x(is)
            enddo
            call biwgt(yb,nb,xbb,xsb)
            yn(nbb)=xbb
            call biwgt(xb,nb,xbb,xsb)
            xn(nbb)=xbb
         enddo
         c1=file1(1:17)
         icp=1
         call pgsls(1)
         call pgslw(2)
         call pgsci(1)
         call pgsch(1.8)
         call pgsch(1.2)
         ybit=(ymax-ymin)/10.
         ymin=ymin-ybit
         ymax=ymax+ybit
         ymin=0.65
         ymax=1.4
         if(il.eq.1) call pgenv(xmin,xmax,ymin,ymax,0,0)
c         call pgenv(xmin,xmax,ymin,ymax,0,0)
         ic=ic+1
         if(ic.eq.13) ic=1
         call pgsci(ic)
         call pgline(n,x,y)
c         call pgsci(ic+1)
c         call pgline(nbb,xn,yn)
         call pgsch(1.8)
         do k=1,n
            call xlinint(x(k),nbb,xn,yn,y0)
            call xlinint(x(k),nsky,wsky,snorm,s0)
            soff=y0*s0
            rat=y(k)-soff
            rat1=rat
            agn=agn0+(x(k)-wave0)*slope
            rat=rat-agn*xnorm
            write(11,*) x(k),rat,rat1
         enddo
c         if(il.eq.1) call pglabel('Wavelength',
c     $        '1e-17 ergs/cm\U2\D/s/AA','')
c         if(il.eq.1) call pglabel('Wavelength',
c     $        'RMS','')
c         if(il.eq.1) call pglabel('Wavelength',
c     $        'Counts','')
c         call pgmtxt('T',0.9,0.5,0.5,file1(1:nf))
         call pgsch(1.5)
         call pgsci(1)
 866     continue
         close(2)
      enddo
 666  continue
      close(1)

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
