
      parameter(nmax=10000)
      real x(nmax),y(nmax)
      real xa(nmax),ya(nmax)
      character file1*80

      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=log10(0.015)
      xmax=log10(120.)
      ymin=log10(0.9)
      ymax=log10(2.e6)
      call pgsls(1)
      call pgsch(1.5)
      call pgenv(xmin,xmax,ymin,ymax,0,30)
      call pglabel("Radius (pc)","density n_stars/pc\U3","")
      call pgsch(1.5)
      call pgslw(6)

      open(unit=1,file='dlist',status='old')

      nl=0
      ic=1
      do il=1,11
         read(1,*,end=666) file1,ilog
         nl=nl+1
         open(unit=2,file=file1,status='old')
         read(2,*)
         n=0
         na=0
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3
            if(i.eq.1) print *,il,x1,x2,x3
            if(il.eq.9) goto 866
            if(il.eq.10) goto 866
c            if(il.eq.11) goto 866
            if(il.eq.7) then
               x1=x1+log10(0.06)
               x2=x2+log10(0.06)
               x3=x3+log10(0.06)
            endif
            if(il.eq.3) then
               x2=x2-log10(10.)
               x3=x3-log10(10.)
            endif
c            x2=x3
            if(x2.gt.-666) then
               n=n+1
               if(ilog.eq.1) then
                  x(n)=x1
                  y(n)=x2
               elseif(ilog.eq.2) then
                  x(n)=x1+log10(0.04)
                  y(n)=log10(x2*1.e6)
               elseif(ilog.eq.0) then
                  x(n)=log10(x1)
                  y(n)=log10(x2)
               endif
            endif
         enddo
 667     continue
         close(2)
         ic=ic+1
         call pgsci(ic)
         call pgline(n,x,y)
 866     continue
         close(2)
      enddo
 666  continue
      close(1)

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
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end
