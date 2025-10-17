
      parameter(nmax=300000)
      real xb(nmax),xg(nmax),sn(nmax),y(nmax),sn2(nmax)
      real xwl(10),xwh(10)

      snmax=200.
      wavemax=8000.

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      xmin=4.7
      xmax=6.7
      ymin=0.
      ymax=0.45

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("S/N","Cumulative Unconfirmed Rate","")
      y(1)=0.1
      y(2)=0.1
      sn(1)=xmin
      sn(2)=xmax
      call pgline(2,sn,y)
      call pgsci(4)
      call pgptxt(6.1,0.27,0.,0.,'350-390')
      call pgsci(3)
      call pgptxt(6.1,0.25,0.,0.,'390-450')
      call pgsci(2)
      call pgptxt(6.1,0.23,0.,0.,'450-550')

      nw=3
      xwl(1)=3500.
      xwh(1)=3900.
      xwl(2)=3900.
      xwh(2)=4500.
      xwl(3)=4500.
      xwh(3)=5500.
      do iw=1,3
         wmin=xwl(iw)
         wmax=xwh(iw)

      open(unit=1,file='b2',status='old')
      nb=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3 
         if(x1.lt.snmax.and.x3.ge.wmin.and.x3.le.wmax) then
            nb=nb+1
            xb(nb)=x1
         endif
      enddo
 666  continue
      close(1)
      open(unit=1,file='g2',status='old')
      ng=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3
         if(x1.lt.snmax.and.x3.ge.wmin.and.x3.le.wmax) then
            ng=ng+1
            xg(ng)=x1
         endif
      enddo
 667  continue
      close(1)

      n=15
      do i=1,n
         sn(i)=xmin+(xmax-xmin)*float(i-1)/float(n-1)
         n1=0
         do j=1,nb
            if(xb(j).gt.sn(i)) n1=n1+1
         enddo
         n2=0
         do j=1,ng
            if(xg(j).gt.sn(i)) n2=n2+1
         enddo
         y(i)=float(n1)/float(n2+n1)
c         print *,sn(i),y(i)
      enddo

      if(iw.eq.1) call pgsci(4)
      if(iw.eq.2) call pgsci(3)
      if(iw.eq.3) call pgsci(2)
      call pgslw(5)
      call pgline(n,sn,y)
      call pgslw(2)

      nn=5
      sn2(1)=4.8
      sn2(2)=5.0
      sn2(3)=5.2
      sn2(4)=5.4
      sn2(5)=10.0
      do i=2,nn
         snlo=sn2(i-1)
         snhi=sn2(i)
         sn(i-1)=(snlo+snhi)/2.
         n1=0
         do j=1,nb
            if(xb(j).ge.snlo.and.xb(j).lt.snhi) n1=n1+1
         enddo
         n2=0
         do j=1,ng
            if(xg(j).ge.snlo.and.xg(j).lt.snhi) n2=n2+1
         enddo
         if(float(n2+n1).gt.0) then
            y(i-1)=float(n1)/float(n2+n1)
         else
            y(i-1)=0.
         endif
         print *,iw,sn(i-1),y(i-1),n1,n2+n1
      enddo
      if(iw.eq.1) call pgsci(4)
      if(iw.eq.2) call pgsci(3)
      if(iw.eq.3) call pgsci(2)
      call pgslw(5)
      call pgsls(4)
      call pgline(nn-2,sn,y)
      call pgslw(2)
      call pgsls(1)

      enddo

      call pgend
      end
