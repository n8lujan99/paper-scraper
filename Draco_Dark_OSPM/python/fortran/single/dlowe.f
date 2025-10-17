
c -- makes a lowess fit to data
      
      parameter(nmax=100000,nsimt=1000)

      real x(nmax),y(nmax),ys(nmax),rw(nmax),res(nmax),yinp(nmax)
      real yhigh(nmax),ylow(nmax),ya(nmax,nsimt),ymed(nmax)
      real yp(nmax)
      real wksp(nmax),ysold(nmax),yeold(nmax),yin(nsimt)
      integer iwksp(nmax),ip(nmax)
      character filein*30,title*40

      data big,vdlow /1.e10,0.1/

 10   call qc1('Velocity data ','dlowe.def',filein)
      open(unit=1,file=filein,status='old',err=10)

      call qi1('Number of simulations ','dlowe.def',nsim)

      call qr2('Input the min and max radius ','dlowe.def',rlow,rhigh)
      call qr1('Input mean velocity ','dlowe.def',vel)
      call qr1('Input F ','dlowe.def',f)
      call qc1('Plot Title ','dlowe.def',title)
      call savdef

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)

 1069 format(a)
c 1070 format(q,a:)

      xmin=big
      xmax=-big
      ymin=big
      ymax=-big

      n=0
      do i=1,nmax
c         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7,i8
         read(1,*,end=666) x1,x2,x3,x4
         rad=sqrt(x1*x1+x2*x2)
         if(rad.ge.rlow.and.rad.le.rhigh) then
            n=n+1
            x(n)=rad
c            x(n)=log10(rad)
            yp(n)=(x3-vel)*(x3-vel)
            yeold(n)=x4*x4
            ip(n)=17
c            if(i8.eq.577.or.i8.eq.251) ip(n)=24
            xmin=min(xmin,rad)
            xmax=max(xmax,rad)
         endif
      enddo
 666  continue
      close(1)

      print *,'Number of stars = ',n

      call sort3(n,x,yp,yeold,wksp,iwksp)

      do i=1,n
         yinp(i)=yp(i)-yeold(i)
      enddo

      call lowess2(x,yinp,yeold,n,f,1,0.,ys,rw,res)

      do i=1,n
         ysold(i)=sqrt(max(ys(i),vdlow))
         yp(i)=sqrt(yp(i))
         ymin=min(ymin,yp(i))
         ymax=max(ymax,yp(i))
      enddo

      idum=-1

      do isimt=1,nsim

         do i=1,n
            y(i)=(ysold(i)*gasdev(idum)+
     $           sqrt(yeold(i))*gasdev(idum))**2-
     $           yeold(i)
         enddo

         call lowess2(x,y,yeold,n,f,1,0.,ys,rw,res)
         
         do i=1,n
            ya(i,isimt)=sqrt(max(ys(i),vdlow))
         enddo

         write(*,'(a8,i4,1x,f7.1,a1,$)')
     $        'Did sim ',isimt,ya(n/2,isimt),char(13)
         call flush(6)

      enddo
      
      do i=1,n
         do j=1,nsim
            yin(j)=ya(i,j)
         enddo
         call biwgt(yin,nsim,xb,xs)
         bias=xb-ysold(i)
c         i5=max(1,nint(float(nsim)*.05))
c         i95=min(nsim,nint(float(nsim)*.95))
         i5=max(1,nint(float(nsim)*.16))
         i95=min(nsim,nint(float(nsim)*.84))
         ylow(i)=yin(i5)-bias
         yhigh(i)=yin(i95)-bias
         ymed(i)=xb-bias
c         print *,i,ylow(i),xb,yhigh(i)
      enddo

      xbit=(xmax-xmin)/10.
c      xmin=xmin-xbit
c      xmin1=xmin-xbit/100.
      xmin1=xmin
c      xmax=xmax+xbit/2.
      ybit=(ymax-ymin)/20.
c      ymin=ymin-ybit
c      ymax=ymax+ybit

c      xmax=log10(xmax)
      if(rlow.eq.0.) then
c         xmin=log10(xmin1)-.05
      else
c         xmin=log10(rlow)-.02
      endif
      print *,xmin,xmax

      call pgscf(2)
      call pgslw(1)
      call pgsch(1.25)
      call pgvport(.15,.85,.15,.85)
c      call pgwindow(xmin,xmax,-1.,ymax,0,0)
      call pgwindow(xmin,xmax,ymin,ymax,0,0)
c      call pgbox('bclnst',0.,0,'bcnst',0.,0)
      call pgbox('bcnst',0.,0,'bcnst',0.,0)
      call pglabel('R (arcsec)','',title)
c      call pglabel('','|\\gDV| (km s\\U-1\\D)','')
c      call pgsch(.5)
      call pgsch(1.25)
      do i=1,n
         call pgpoint(1,x(i),yp(i),ip(i))
      enddo
      call pgsch(1.25)
      call pgline(n,x,ymed)
      call pgsls(2)
      call pgline(n,x,ylow)
      call pgline(n,x,yhigh)
      call pgsls(4)

      open(unit=2,file='dlowe.out',status='unknown')
      do i=1,n
         write(2,*) 10**x(i),ymed(i),ylow(i),yhigh(i),yp(i),ysold(i)
      enddo
      close(2)

      call pgend

      end
