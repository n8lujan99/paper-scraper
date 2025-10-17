Skipped
      parameter(nmax=100000)
      real x(nmax),xb(nmax),yb(nmax),g(nmax)
      real xbf(nmax),ybf(nmax),xbl(nmax),xbh(nmax),xbfl(nmax),xbfh(nmax)
      character title*10

      sig=1.08
      chicut=1.3

      title='          '
      open(unit=1,file='title',status='old',err=555)
      read(1,*) title
 555  close(1)

      open(unit=1,file='in',status='old')

      xmax=0.
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         if(x2.lt.chicut) then
            n=n+1
            x(n)=x1
            xmax=max(xmax,x1)
         endif
      enddo
 666  continue
      close(1)

      xmax=xmax*1.1
c      xmin=-7.1
c      xmax=7.1
      xmin=-xmax
      nb=50

      ymax=-1e10
      ymin=0.
      sumt=0
      do i=1,nb-1
         xl=xmin+(xmax-xmin)*float(i-1)/float(nb-1)
         xh=xmin+(xmax-xmin)*float(i)/float(nb-1)
         xb(i)=(xl+xh)/2.
         xbl(i)=xl
         xbh(i)=xh
         sum=0.
         do j=1,n
            if(x(j).ge.xl.and.x(j).lt.xh) sum=sum+1
         enddo
         yb(i)=sum
         ymax=max(ymax,yb(i))
         sumt=sumt+yb(i)
      enddo

      con=0.
      amp=sumt*(xmax-xmin)/float(nb-1)
      nf=0
      sumc=0.
      do i=1,nb-1
         w=(xb(i)-0.)/sig
         gaus=exp(-w*w/2.)/sqrt(2.*3.1415*sig*sig)
         g(i)=con+amp*gaus
         if(xb(i).le.0) then
            nf=nf+1
            xbf(nf)=-xb(i)
            xbfl(nf)=-xbl(i)
            xbfh(nf)=-xbh(i)
            ybf(nf)=yb(i)
            sumc=sumc+ybf(nf)
c            print *,xbf(nf),sumc/sumt,sumt
         endif
      enddo
      ymax=ymax*1.05

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(3)

      igaus=0
      if(igaus.eq.1) then
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgslw(5)
         call pgline(nb-1,xb,yb)
         call pgsci(2)
         call pgline(nb-1,xb,g)
      else
         call pgenv(0.,xmax,ymin,ymax,0,0)
         call pglabel('S/N','N',title)
         call pgslw(5)
c         call pgline(nb-1,xb,yb)
         call kghist(nb-1,xbl,xbh,yb)
         call pgsci(2)
c         call pgline(nf,xbf,ybf)
         call kghist(nf,xbfl,xbfh,ybf)
      endif
         
      call pgend

      end

      subroutine kghist(n,xl,xh,y)
      real x(n),xl(n),xh(n),y(n)
      real xll(2),yll(2)

      do i=1,n
         y0=0.
         if(i.gt.1) y0=y(i-1)
         xll(1)=xl(i)
         xll(2)=xl(i)
         yll(1)=y0
         yll(2)=y(i)
         if(y(i).ne.0.) call pgline(2,xll,yll)
         xll(1)=xl(i)
         xll(2)=xh(i)
         yll(1)=y(i)
         yll(2)=y(i)
         if(y(i).ne.0.) call pgline(2,xll,yll)
         y0=0.
         if(i.lt.n) y0=y(i+1)
         xll(1)=xh(i)
         xll(2)=xh(i)
         yll(1)=y(i)
         yll(2)=y0
         if(y(i).ne.0.) call pgline(2,xll,yll)
      enddo
         
      return
      end
