
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),xn2(nmax),yn2(nmax)
      real xb(nmax),yb(nmax),xs(nmax),ys(nmax)
      real wl(nmax),wh(nmax),xwe(nmax),xecor(nmax),rata(100,nmax)
      real xin(nmax),xwp(nmax),xpcor(nmax),rg(nmax)
      integer igg(nmax),ig(nmax)
      character file1*10,file2*15,c1*15,nullstr*1,f1*80,a1*10
      character cdg(nmax)*15
      logical simple,extend,anyf

      wc=50.
      ibin1=1
      ibin2=1

      open(unit=12,file='ofile',status='unknown')

      open(unit=1,file="zdex_desi.txt",status='old')
      ng=0
      do i=1,nmax
         read(1,*,end=555) a1,x2,x3,x4,i5,i6,x7
         ng=ng+1
         cdg(ng)=a1
         igg(ng)=1
         if(i6.gt.0) igg(ng)=0
c         if(x7.gt.0.4) igg(ng)=0
         ig(ng)=i6
         rg(ng)=x7
      enddo
 555  continue
      close(1)

      call pgbegin(0,'?',4,4)
c      call pgask(.false.)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      open(unit=1,file='listin',status='old')

      nl=0
      ic=0
      do il=1,20000
         read(1,*,end=666) file1,file2,wave
         f1="dexspec/"//file1//".spec"
         open(unit=2,file=f1,status='old',err=888)
         goto 889
 888     continue
         print *,"Does not exist: ",f1
         goto 866
 889     continue
         igood=1
         do j=1,ng
            if(file1.eq.cdg(j)) then
               igood=igg(j)
               rgw=rg(j)
            endif
         enddo
         
         ymin=1e10
         ymax=-ymin
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3
            if(x2.ne.0) then
               n=n+1
               x(n)=x1
               y(n)=x2/2.
            endif
         enddo
 667     continue
         close(2)

         f1="desispec/desi_"//file2//".spec"
         open(unit=2,file=f1,status='old',err=788)
         goto 789
 788     continue
         print *,"Does not exist: ",f1
         goto 866
 789     continue
         ns=0
         do i=1,nmax
c            read(2,*,end=567) x1,x2,x3
            read(2,*,end=567) x1,x2
            if(x2.ne.0) then
               ns=ns+1
               xs(ns)=x1
               ys(ns)=x2
            endif
         enddo
 567     continue
         close(2)
         if(ns.eq.0) then
            print *,"ns=0 ",f1
            goto 866
         endif

         call smbin(ibin1,n,x,y,nbb,xn,yn,ymin,ymax)

         call smbin(ibin2,ns,xs,ys,nbb2,xn2,yn2,ymin,ymax)

         call pgsls(1)
         call pgslw(2)
         call pgsci(1)
         call pgsch(1.8)

         xmin=wave-wc
         xmax=wave+wc

         ymin=1e10
         ymax=-ymin
         do i=1,n
            if(x(i).gt.xmin.and.x(i).lt.xmax) then
               ymin=min(ymin,y(i))
               ymax=max(ymax,y(i))
            endif
         enddo
         do i=1,ns
            if(xs(i).gt.xmin.and.xs(i).lt.xmax) then
               ymin=min(ymin,ys(i))
               ymax=max(ymax,ys(i))
            endif
         enddo
         call pgsci(1)
         if(igood.eq.0) call pgsci(4)
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgsci(1)
c         call pgline(nbb,xn,yn)
         call pgline(n,x,y)
         call pgsch(1.6)
         call pgsci(1)
         call pglabel('Wavelength \(2078)',
     $        '1e-17 ergs/cm\U2\D/s/\(2078)','')

         call pgsch(1.9)
         call pgsci(1)
         call pgmtxt('T',1.5,0.1,0.5,file1)
         call pgsci(2)
         call pgmtxt('T',1.5,0.75,0.5,file2)
         call pgsci(2)
c         call pgline(nbb2,xn2,yn2)
         call pgline(ns,xs,ys)

c         read *,iout
         iout=-666
         write(12,*) iout," ",file1," ",file2

         call pgsch(1.6)
         call pgsci(1)

 866     continue
         close(2)
      enddo
 666  continue
      close(1)
      close(12)

      call pgend
 1105 format(f5.3)

      end
      
      subroutine getrat(nbb,xn,yn,nbb2,xn2,yn2,nw,wc,wl,wh,rat)
      real xn(nbb),yn(nbb),xn2(nbb2),yn2(nbb2),wc(nw),rat(nw)
      real wl(nw),wh(nw),xin1(1000),xin2(1000),yin1(1000),yin2(1000)
      real rmsa(1000)

      do i=1,nw
         sum1=0.
         n1=0
         do j=1,nbb
            if(xn(j).ge.wl(i).and.xn(j).lt.wh(i)) then
               sum1=sum1+yn(j)
               n1=n1+1
               xin1(n1)=xn(j)
               yin1(n1)=yn(j)
            endif
         enddo
         sum1=sum1/float(n1)
         sum2=0.
         n2=0
         do j=1,nbb2
            if(xn2(j).ge.wl(i).and.xn2(j).lt.wh(i)) then
               sum2=sum2+yn2(j)
               n2=n2+1
               xin2(n2)=xn2(j)
               yin2(n2)=yn2(j)
            endif
         enddo
         sum2=sum2/float(n2)

         rat(i)=sum1/sum2
         
         ktry=100
         rstep=0.005
         rmin=rat(i)-rstep*float(ktry)/2.
         rmax=rat(i)+rstep*float(ktry)/2.
         rmsl=1.e10
         do k=1,ktry
            rtry=rmin+(rmax-rmin)*float(k-1)/float(ktry-1)
            do j=1,n1
               call xlinint(xin1(j),nbb2,xn2,yn2,y0)
               ynew=yin1(j)/rtry
               rmsa(j)=(y0-ynew)**2
            enddo
            call biwgt(rmsa,n1,xbrms,xsrms)
            if(xbrms.lt.rmsl) then
               rmsl=xbrms
               rbest=rtry
            endif
         enddo
c         print *,i,rat(i),rbest
         rat(i)=rbest
      enddo

      return
      end

      subroutine smbin(ibin,n,x,y,nbb,xn,yn,ymin,ymax)
      real x(n),y(n),xb(10000),yb(10000),xn(n),yn(n)

      ib1=(ibin-1)/2
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
c         yn(nbb)=xbb
         yn(nbb)=sum/float(nb)
         call biwgt(xb,nb,xbb,xsb)
         xn(nbb)=xbb
         ymin=min(ymin,yn(nbb))
         ymax=max(ymax,yn(nbb))
      enddo

      return
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
