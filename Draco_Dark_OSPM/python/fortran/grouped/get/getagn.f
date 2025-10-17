
      parameter(nmax=10000)
      real x(nmax),y(nmax),xin(nmax),wg(nmax),sg(nmax),sgp(nmax)
      real yn0(nmax),xin2(nmax),xin3(nmax),xin4(nmax)
      real agna(nmax),ap(nmax),stot(nmax),rmsa(nmax)
      real xbadl(nmax),xbadh(nmax)
      real xn(nmax),yn(nmax),xb(nmax),yb(nmax)
      character file1*80,file2*80,c1*18

      wmin=22600.
      wmax=22800.
      w1=22500.
      w2=24000.

      idmax=3
c      idmax=0

      nbad=1
      xbadl(1)=1.
      xbadh(1)=1.

      open(unit=2,file='in',status='old')
      read(2,*) x1,x2,x3,i4,rad
      close(2)

      wave0=22000.
      agn0=1250.
      wave1=24000.
      agn1=1600.
      agn1=1750.
c      if(rad.gt.15.) agn1=1470.
      if(rad.gt.15.) agn1=1500.
      slope=(agn1-agn0)/(wave1-wave0)

      open(unit=1,file='good.spec',status='old')
      ng=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         ng=ng+1
         wg(ng)=x1
         sg(ng)=x2
         agna(ng)=agn0+(wg(ng)-wave0)*slope
         agna(ng)=agna(ng)/agn0
      enddo
 666  continue
      close(1)

      open(unit=2,file='in',status='old')
      open(unit=11,file='out',status='unknown')

      xmin=22000.
      xmax=24200.
      ymin=1.e10
      ymax=0.
      n=0
      do i=1,nmax
         read(2,*,end=667) x1,x2,x3
         n=n+1
         x(n)=x1
         y(n)=x2
         yn0(n)=x3
         if(x1.gt.xmin.and.x1.lt.xmax) then
            ymin=min(ymin,y(n))
            ymax=max(ymax,y(n))
         endif
      enddo
 667  continue
      close(2)

      ism=0
      if(ism.eq.1) then
         ibin=5
         ib1=(ibin-1)/2
         xib=float(ibin)
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
         ymin=1e10
         ymax=0.
         do i=1,n
            call xlinint(x(i),nbb,xn,yn,y0)
            y(i)=y0
            if(x(i).gt.xmin.and.x(i).lt.xmax) then
               ymin=min(ymin,y(i))
               ymax=max(ymax,y(i))
            endif
         enddo
      endif

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pgline(n,x,y)
c      call pgline(nbb,xb,yb)    

      idone=0
      ns=60
      smin=30.
      smax=80.
      if(rad.lt.5.) smax=160.
      na=50
      amin=0.
      amax=1500.
 888  continue
      diff=1.e30
      do ia=1,na
         atry=amin+float(ia-1)*(amax-amin)/float(na-1)
      do i=1,ns
         stry=smin+float(i-1)*(smax-smin)/float(ns-1)
         rms=0.
         nrms=0
         do j=1,n
            if(x(j).gt.xmin.and.x(j).lt.xmax) then
               call xlinint(x(j),ng,wg,sg,s0)
               call xlinint(x(j),ng,wg,agna,a0)
               s0=s0*stry
               a0=a0*atry
               rms=rms+(s0+a0-y(j))**2
               nrms=nrms+1
               rmsa(nrms)=(s0+a0-y(j))**2
            endif
         enddo
c         rms=sqrt(rms)
         call biwgt(rmsa,nrms,xb0,xs0)
         rms=sqrt(xb0)
         if(rms.lt.diff) then
            diff=rms
            sbest=stry
            abest=atry
         endif
      enddo
      enddo
      do i=1,ng
         sgp(i)=sg(i)*sbest
         ap(i)=agna(i)*abest
         stot(i)=sgp(i)+ap(i)
      enddo
      print *,idone,sbest,abest,diff
      if(idone.lt.idmax) then
         idone=idone+1
         ns=60
         smin=sbest*0.8
         smax=sbest*1.2
         na=50
         amin=abest*0.8
         amax=abest*1.2
         if(abest.lt.20.) then
            na=20
            amin=0.
            amax=20.
            if(abest.lt.9.) amax=10.
         endif
         goto 888
      endif

      call pgsci(4)
      call pgline(ng,wg,sgp)
      call pgsci(3)
      call pgline(ng,wg,ap)
      call pgsci(2)
      call pgline(ng,wg,stot)

      open(unit=11,file='out2',status='unknown')
      write(11,*) sbest,abest
      close(11)

      open(unit=11,file='out3',status='unknown')
      do i=1,ng
         call xlinint(wg(i),n,x,y,y0)
         write(11,*) wg(i),y0,stot(i),sgp(i),ap(i)
      enddo
      close(11)

      open(unit=11,file='out4',status='unknown')
      do i=1,n
         call xlinint(x(i),ng,wg,ap,a0)
         y(i)=y(i)-a0
      enddo
      nb=0
      do i=1,n
         if(x(i).gt.wmin.and.x(i).lt.wmax) then
            nb=nb+1
            xin(nb)=y(i)
         endif
      enddo
      call biwgt(xin,nb,xbnorm,xsnorm)
      do i=1,n
         write(11,*) x(i),y(i)/xbnorm
      enddo
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
