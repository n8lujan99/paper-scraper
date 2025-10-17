
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax),xbadl(nmax),xbadu(nmax)
      real ym(nmax),ysub(nmax),yrat(nmax)
      character file1*80,file2*80,c1*18

      ibin=5
      ibin=27
      ib1=(ibin-1)/2
      xib=float(ibin)

      nbad=0
      open(unit=1,file='baduse.dat',status='old',err=555)
      do i=1,nmax
         read(1,*,end=555) x1,x2
         nbad=nbad+1
         xbadl(nbad)=x1
         xbadu(nbad)=x2
      enddo
 555  continue
      close(1)
      nbad=0

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,2)
      call pgpap(0.,1.)
c      call pgpap(0.,0.5)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=21000.
      xmax=24100.
      
      open(unit=1,file='in',status='old')
      open(unit=11,file='out',status='unknown')

      ymin=1e10
      ymax=-ymin
      n=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3,x4,x5
         n=n+1
         x(n)=x1
         y(n)=x2
         ym(n)=x3
         ysub(n)=y(n)-ym(n)
         if(ym(n).gt.0.) then
            yrat(n)=y(n)/ym(n)
         else
            yrat(n)=1.
         endif
      enddo
 667  continue
      close(1)

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
            sum=sum+yrat(is)
            nb=nb+1
            yb(nb)=yrat(is)
            xb(nb)=x(is)
         enddo
         call biwgt(yb,nb,xbb,xsb)
         yn(nbb)=xbb
         call biwgt(xb,nb,xbb,xsb)
         xn(nbb)=xbb
      enddo

      ymin=-10.
      ymax=10.
      call pgsch(1.3)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength (\(2078))','Subtracted Flux','')
      call pgline(n,x,ysub)
      ymin=0.7
      ymax=1.3
      call pgsch(1.3)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Wavelength (\(2078))','Divided Flux','')
      call pgline(n,x,yrat)
      call pgsci(2)
      call pgline(nbb,xn,yn)
      call pgsci(1)

      do i=1,n
         if(yrat(i).gt.1.17.or.yrat(i).lt.0.83) then
            sub=ysub(i)
            rat=1.
         else
            call xlinint(x(i),nbb,xn,yn,y0)
            sub=0.
            yrat(i)=yrat(i)/y0
            rat=yrat(i)
         endif
         write(11,*) x(i),sub,rat,yrat(i)
      enddo
      call pgsci(3)
      call pgslw(1)
      call pgline(n,x,yrat)
      call pgsci(1)

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
