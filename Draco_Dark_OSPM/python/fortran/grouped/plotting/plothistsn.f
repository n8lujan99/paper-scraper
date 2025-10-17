
      parameter(nmax=18000000)
      real x(nmax),sn(nmax),xp(1000),yp(1000),xn(nmax),sx(10)
      real snl(10),snh(10),yp1(1000),xsn(10),csn(10),ssn(10)
      character file1*80

      nsn=5
      xsn(1)=4.8
      csn(1)=0.7
      xsn(2)=5.0
      csn(2)=0.8
      xsn(3)=5.5
      csn(3)=0.9
      xsn(4)=6.0
      csn(4)=0.98
      xsn(5)=6.5
      csn(5)=1.0

      nss=8
      sx(1)=4.8
      sx(2)=5.0
      sx(3)=5.4
      sx(4)=5.8
      sx(5)=6.7
      sx(6)=7.0
      sx(7)=7.5
      sx(8)=8.5
      ssn(1)=1.0
      ssn(2)=0.96
      ssn(3)=0.88
      ssn(4)=0.87
      ssn(5)=0.86
      ssn(6)=0.70
      ssn(7)=0.65
      ssn(8)=0.60

      file1="in_hdr5_j1"
      open(unit=1,file=file1,status='old')

      n=0
      xmin=4.8
      xmax=7.5
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3
         n=n+1
         sn(n)=x3
      enddo
 666  continue
      close(1)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      call pgenv(xmin,xmax,0.,1.05,0,0)
      call pglabel('S/N','Normalized Counts','')
      call pgslw(5)
      
      nbin=30
      xbin=(xmax-xmin)/float(nbin-1)
      ymax=0.
      sum=0.
      do i=1,nbin
         xlo=xmin+xbin*float(i-1)
         xhi=xlo+xbin
         nb=0
         do j=1,n
            if(sn(j).ge.xlo.and.sn(j).lt.xhi) nb=nb+1
         enddo
         xp(i)=(xhi+xlo)/2.
         yp(i)=float(nb)
         call xlinint(xp(i),nsn,xsn,csn,cfrac)
         yp(i)=yp(i)*cfrac
         sum=sum+yp(i)
         ymax=max(ymax,yp(i))
      enddo
      sum=sum/5.
      sum=ymax
      sumt=0.
      do i=1,nbin
         yp(i)=yp(i)/sum
         sumt=sumt+yp(i)
      enddo
      call pgline(nbin,xp,yp)

      ic=1
      open(unit=1,file='hlist',status='old')
      do iall=1,1000
         read(1,*,end=667) file1,ymax2
         open(unit=2,file=file1,status='old')
         n=0
         do j=1,nmax
            read(2,*,end=668) x1
            n=n+1
            cfrac=1.0
            if(iall.eq.2) call xlinint(x1,nss,sx,ssn,cfrac)
            sn(n)=x1*cfrac
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
               if(sn(j).ge.xlo.and.sn(j).lt.xhi) nb=nb+1
            enddo
            xp(i)=(xhi+xlo)/2.
            yp(i)=float(nb)
            sum=sum+yp(i)
            ymax=max(ymax,yp(i))
         enddo
         sum=sum/5.
         sum=ymax2
         sumt2=0.
         do i=1,nbin
            yp(i)=yp(i)/sum
            sumt2=sumt2+yp(i)
         enddo
         print *,sumt,sumt2,sumt2/sumt
         ic=ic+1
         if(ic.eq.3) ic=ic+1
         call pgsci(ic)
         call pgline(nbin,xp,yp)
      enddo
 667  continue
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
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
