
      parameter(nmax=300000)
      real xb(nmax),xg(nmax),sn(nmax),y(nmax),sn2(nmax)

      snmax=200.
      wavemax=8000.

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=1,file='b2',status='old')
      nb=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3
         if(x1.lt.snmax.and.x3.lt.wavemax) then
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
         if(x1.lt.snmax.and.x3.lt.wavemax) then
            ng=ng+1
            xg(ng)=x1
         endif
      enddo
 667  continue
      close(1)

      xmin=4.7
      xmax=6.7
      ymin=0.
      ymax=0.20

      n=20
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
         print *,sn(i),y(i)
      enddo

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("S/N","Unconfirmed Rate","")
      call pgslw(5)
      call pgline(n,sn,y)
      call pgslw(2)
      y(1)=0.1
      y(2)=0.1
      sn(1)=xmin
      sn(2)=xmax
      call pgsci(2)
      call pgline(2,sn,y)

      nn=9
      sn2(1)=4.7
      sn2(2)=4.8
      sn2(3)=4.9
      sn2(4)=5.0
      sn2(5)=5.1
      sn2(6)=5.2
      sn2(7)=5.4
      sn2(8)=7.0
      sn2(9)=10.0
c      sn2(10)=10.0
      do i=2,nn
         snlo=sn2(i-1)
         snhi=sn2(i)
         sn(i)=(snlo+snhi)/2.
         n1=0
         do j=1,nb
            if(xb(j).ge.snlo.and.xb(j).lt.snhi) n1=n1+1
         enddo
         n2=0
         do j=1,ng
            if(xg(j).ge.snlo.and.xg(j).lt.snhi) n2=n2+1
         enddo
         if(float(n2+n1).gt.0) then
            y(i)=float(n1)/float(n2+n1)
         else
            y(i)=0.
         endif
         print *,sn(i),y(i),n1,n2+n1
      enddo

      call pgend
      end
