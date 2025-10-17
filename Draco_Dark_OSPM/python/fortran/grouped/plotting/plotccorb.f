
      parameter(nmax=10000,nmax2=200)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),yin(nmax)
      real xb(nmax),yb(nmax),xw(nmax2),xf(nmax2),ysa(nmax2,nmax2)
      real xs(nmax2),ys(nmax2),ya(nmax2,nmax2),xin(nmax2),xcl(7)
      real yse(nmax2),ysae(nmax2,nmax2)
      character file1*80,file2*80,c1*18,a1*4,title(10)*15

      title(1)="3550-3820"
      title(2)="3820-4093"
      title(3)="4093-4364"
      title(4)="4364-4635"
      title(5)="4635-4907"
      title(6)="4907-5178"
      title(7)="5178-5450"
      ns=199
      xmin=1.0
      xmax=100.
      do i=1,ns
         xs(i)=xmin+float(i-1)*(xmax-xmin)/float(ns-1)
      enddo

      open(unit=1,file="bsn.in",status='old')
      read(1,*) xn1,xn2,xn3,xn4,xn5,xn6,xn7
      read(1,*) xnb1,xnb2,xnb3,xnb4,xnb5,xnb6,xnb7
      read(1,*) xcl(1),xcl(2),xcl(3),xcl(4),xcl(5),xcl(6),xcl(7)
c      read(1,*) ymin1,yset1,yset2
      close(1)
      ymin1=0.00001
      yset1=0.006
      yset2=0.0007

c      xn1=15.63
c      xn2=9.11
c      xn3=7.23
c      xn4=6.84
c      xn5=6.16
c      xn6=6.01
c      xn7=6.49

c      ibin=13
      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      call pgbegin(0,'?',3,2)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=0.0
      xmax=40.0
c      xmin=1.
c      xmax=9.0

      ymin=0.0
      ymax=1.05
c      ymax=0.02

      open(unit=1,file='flist',status='old')

c      do iw=1,7
      do iw=1,3
      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      call pgsch(1.9)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pgsch(2.1)
      call pglabel('Flux (1e-17 ergs/cm\U2\D/s)','Completeness','')
      call pgsch(2.4)
      call pgmtxt('T',1.0,0.5,0.5,title(iw))
      call pgmtxt('T',1.0,0.83,0.5,"\(2078")
      call pgsch(2.1)
      call pgslw(3)
      call pgsch(1.)

      ic=1
      ntot=0
      do il=1,10000
         read(1,*,end=666) file1
         ntot=ntot+1
         open(unit=2,file=file1,status='old')
         if(ntot.eq.1) then
            read(2,*) a1,i1,i2,i3
         else
            read(2,*)
         endif
         read(2,*) x0,x1,x2,x3
         xw(1)=x1/xn1
         xw(2)=x2/xn2
         xw(3)=x3/xn3
         xw(4)=x4/xn4
         xw(5)=x5/xn5
         xw(6)=x6/xn6
         xw(7)=x7/xn7
         scale=xw(iw)
c         scale=1.
         n=0
         do i=1,nmax
            read(2,*,end=667) x0,x1,x2,x3
            n=n+1
            x(n)=x0/scale
            xf(1)=x1/xnb1
            xf(2)=x2/xnb2
            xf(3)=x3/xnb3
            xf(4)=x4/xnb4
            xf(5)=x5/xnb5
            xf(6)=x6/xnb6
            xf(7)=x7/xnb7
            y(n)=xf(iw)
         enddo
 667     continue
         close(2)
         ic=ic+1
         if(ic.eq.14) ic=2
         call pgsci(ic)
         call pgline(n,x,y)
         do is=1,ns
            call xlinint(xs(is),n,x,y,yp)
            if(x(1).gt.xs(is)) yp=-666
            ya(is,ntot)=yp
         enddo
 866     continue
         close(2)
      enddo
 666  continue

      do is=1,ns
         nin=0
         do i=1,ntot
            if(ya(is,i).ge.0) then
               nin=nin+1
               xin(nin)=ya(is,i)
            endif
         enddo
         call biwgt(xin,nin,xbout,xsout)
         ys(is)=xbout
         yse(is)=xsout
c         yse(is)=(xin(nin)-xin(1))/2.
         ysa(iw,is)=min(xbout,1.0)
         ysae(iw,is)=yse(is)
      enddo

      do icheck=3,1,-1
         if(ys(icheck).lt.0) then
            ys(icheck)=ys(icheck+1)
            yse(icheck)=yse(icheck+1)
            ysa(iw,icheck)=ys(icheck)
            ysae(iw,icheck)=yse(icheck)
         endif
      enddo

c      goto 900
      xs0=10.
      do is=ns-1,1,-1
         if(xs(is).le.xs0) then
            if(ys(is).ge.ys(is+1)) ys(is)=ys(is+1)*0.999
            ysa(iw,is)=ys(is)
         endif
      enddo
      call xlinint(yset1,ns,ys,xs,xs0)
      do is=1,ns
         if(xs(is).le.xs0) then
            yval=(xs(is)/xs0)**4
c            yval=(xs(is)/xs0)**6
            yval=yval*yset1+ymin1
            ys(is)=yval
            ysa(iw,is)=ys(is)
         endif
      enddo
      call xlinint(yset2,ns,ys,xs,xs0)
      do is=1,ns
         if(xs(is).le.xs0) then
            yval=(xs(is)/xs0)**2
            yval=yval*yset2+ymin1
            ys(is)=yval
            ysa(iw,is)=ys(is)
         endif
      enddo
 900  continue

      do is=1,ns
c         ys(is)=max(ys(is),xcl(iw))
         if(ys(is).le.xcl(iw)) ys(is)=0.
         ys(is)=min(ys(is),1.0)
         ysa(iw,is)=ys(is)
      enddo

      call pgslw(9)
      call pgsci(1)
c      call pgline(ns,xs,ys)
      call pgslw(3)

      rewind(1)
      enddo

      close(1)

c      xn1=xn1*xnb1
c      xn2=xn2*xnb2
c      xn3=xn3*xnb3
c      xn4=xn4*xnb4
c      xn5=xn5*xnb5
c      xn6=xn6*xnb6
c      xn7=xn7*xnb7

      do iw=1,3
         do i=1,ns
            yin(i)=ysa(iw,i)
         enddo
         call xlinint(0.5,ns,yin,xs,f50)
         xin(iw)=f50
      enddo
c      write(*,1006) 0.5,(xin(i),i=1,7)
 1006 format(1x,f6.4,7(2x,f7.4))

      i0=0
      open(unit=11,file='out2',status='unknown')
      open(unit=12,file='out3',status='unknown')
      write(11,1102) i0,i1,i2,i3
      write(11,1101) 0.5,(xin(i),i=1,3)
      write(12,1102) i0,i1,i2,i3
      write(12,1101) 0.5,(xin(i),i=1,3)
c      write(11,1101) 0.5,xn1,xn2,xn3,xn4,xn5,xn6,xn7
      do i=1,ns
         write(11,1101) xs(i),(ysa(iw,i),iw=1,3)
         write(12,1101) xs(i),(ysae(iw,i),iw=1,3)
      enddo
      close(11)
      close(12)

      call pgend

 1101 format(8(f10.6,1x))
 1102 format(8(i7,1x))

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
      if(xp.ge.x(n)) yp=y(n)
      return
      end
