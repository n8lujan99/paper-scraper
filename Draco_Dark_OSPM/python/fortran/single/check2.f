
      parameter(nmax=1000000)
      real w(nmax),f(nmax),fe(nmax),sig(nmax),sn(nmax)
      real ch(nmax),wi(nmax),fi(nmax),wa(100),da(100),ca(100)
      real xin(nmax),xin2(nmax),xd50(100)
      real xwa(100,100),xca(100,100),xda(100,100)

c      sn0lo=4.8
      read *,sn0lo,sn0up
c      sn0up=1e10
c      sn0=6.0
c      sn0=5.0
      xlw0=1.8
      chi0=1.2
      xlw0=1.7
      chi0=2.2

      open(unit=1,file='a2',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7,x8
         igood=0
         if(x5.ge.sn0lo.and.x5.lt.sn0up.and.x4.ge.xlw0.and.
     $        x6.le.chi0) igood=1
         if(x1.eq.0) igood=1
         if(igood.eq.1) then
            n=n+1
            w(n)=x1
            f(n)=x2
            fe(n)=x3
            sig(n)=x4
            sn(n)=x5
            ch(n)=x6
            wi(n)=x7
            fi(n)=x8
         endif
      enddo
 666  continue
      close(1)

      ng=0
      nb=0
      dlo=25.
      dhi=35.
      wlo=5000.
      whi=5400.
      do i=1,n
         if(fi(i).gt.dlo.and.fi(i).le.dhi.and.
     $      wi(i).gt.wlo.and.wi(i).le.whi) then
            if(w(i).gt.0) then
               ng=ng+1
            else
               nb=nb+1
            endif
         endif
      enddo
      xc0=float(ng)/float(ng+nb)
      xmult=0.97/xc0
c      xmult=1.03
c      xmult=1.
      print *,xc0,xmult,ng,ng+nb

      nd=40
      dmin=2.
      dmax=40.
c      nd=30
c      dmin=4.
c      dmax=25.

      nw=8
      wmin=3550.
      wmax=5450.

      do k=1,nw-1
         wlo=wmin+float(k-1)*(wmax-wmin)/float(nw-1)
         whi=wmin+float(k)*(wmax-wmin)/float(nw-1)
         wv=(wlo+whi)/2.
         na=0
         do j=1,nd-1
            dlo=dmin+float(j-1)*(dmax-dmin)/float(nd-1)
            dhi=dmin+float(j)*(dmax-dmin)/float(nd-1)
            dv=(dlo+dhi)/2.
            ng=0
            nb=0
            do i=1,n
               if(fi(i).gt.dlo.and.fi(i).le.dhi.and.
     $            wi(i).gt.wlo.and.wi(i).le.whi) then
                  if(w(i).gt.0) then
                     ng=ng+1
                  else
                     nb=nb+1
                  endif
               endif
            enddo
            if(ng+nb.gt.0) then
               xc=xmult*float(ng)/float(ng+nb)
            else
               xc=0.
            endif
            iwv=nint(wv)
c            write(*,1001) iwv,xc,dv,ng,ng+nb
c            print *,wv,dv,xc
            xwa(k,j)=wv
            xca(k,j)=xc
            xda(k,j)=dv
            na=na+1
            wa(na)=wv
            da(na)=dv
            ca(na)=xc
         enddo
         do i=1,10
            xv=0.1+0.1*float(i-1)
            if(i.eq.10) xv=0.95
            call xlinint(xv,na,ca,da,dv)
            iwv=nint(wv)
c            write(*,1001) iwv,xv,dv
         enddo
         call xlinint(0.5,na,ca,da,dv)
         xd50(k)=dv
      enddo

      open(unit=11,file='out',status='unknown')
      write(11,1003) " flux ",(nint(xwa(k,1)),k=1,nw-1)
      write(11,1004) " 0.5  ",(xd50(k),k=1,nw-1)
      do j=1,nd-1
         write(11,1002) xda(1,j),(xca(k,j),k=1,nw-1)
      enddo
      close(11)

      s1=4.75
      do is=1,5-1
      s2=s1+0.5

      ns=0
      do i=1,n
         if(w(i).gt.4000..and.w(i).lt.5000.) then
c            if(sn(i).gt.4.8.and.sn(i).le.5.5.and.
            if(sn(i).gt.s1.and.sn(i).le.s2.and.
     $         f(i).gt.9..and.f(i).le.11.) then
               ns=ns+1
               xin(ns)=fi(i)
               xin2(ns)=f(i)
            endif
         endif
      enddo
      call biwgt(xin,ns,xbout,xsout)
      call biwgt(xin2,ns,xbout2,xsout2)
c      print *,(s2+s1)/2.,ns,xbout2,xbout,xsout
      s1=s2
      enddo

      print *
      print *,"S/N  N  Avg_err  RMS   avg_F   Ratio_err "
      filo=10.
      do ifa=1,4
      fihi=filo+10.
      s1=4.75
      do is=1,8-1
      s2=s1+1.0
      ns=0
      do i=1,n
         if(w(i).gt.4000..and.w(i).lt.5000.) then
            if(sn(i).gt.s1.and.sn(i).le.s2.and.
     $         f(i).gt.filo.and.f(i).le.fihi) then
               ns=ns+1
               xin(ns)=fe(i)
               xin2(ns)=f(i)
            endif
         endif
      enddo
      call biwgt(xin,ns,xbout,xsout)
      call biwgt(xin2,ns,xbout2,xsout2)
      print *,(s2+s1)/2.,ns,xbout,xsout2,xbout2,xsout2/xbout
      s1=s2
      enddo
      filo=fihi
      enddo

 1001 format(i4,2x,f5.3,1x,f6.2,2(1x,i5))
 1002 format(f6.2,7(1x,f7.4))
 1003 format(a6,7(3x,i4))
 1004 format(a6,7(1x,f7.4))
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
