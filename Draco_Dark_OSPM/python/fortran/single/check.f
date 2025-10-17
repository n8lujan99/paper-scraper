
      parameter(nmax=1000000)
      real w(nmax),f(nmax),fe(nmax),sig(nmax),sn(nmax)
      real ch(nmax),wi(nmax),fi(nmax),wa(1000),da(1000),ca(1000)
      real xin(nmax),xin2(nmax),xd50(1000)
      real xwa(1000,1000),xca(1000,1000),xda(1000,1000)
      character c1*15

      read *,sn0lo
      sn0up=1e10

      fwhm=0.
      open(unit=1,file='fwhm.use',status='old',err=444)
      read (1,*) c1,fwhm
 444  continue
      close(1)

      sn0=0.90
      if(fwhm.gt.0) then
c         sncorr=sn0+(1.-sn0)*(fwhm-1.2)/(2.5-1.2)
         sncorr=sn0+(1.-sn0)*(fwhm-1.2)/(2.5-1.2)
      else
         sncorr=1.
      endif
      sncorr=min(1.,sncorr)
      print *,fwhm,sncorr

      xlw0=1.8
      chi0=1.2
      xlw0=1.7
      chi0=2.2

      open(unit=1,file='a2',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666,err=755) x1,x2,x3,x4,x5,x6,x7,x8
         goto 756
 755     print *,"bad line at: ",i
         goto 757
 756     continue
         igood=0
         x5=x5*sncorr
c         if(x5.ge.sn0lo.and.x5.lt.sn0up.and.x4.ge.xlw0.and.
c     $        x6.le.chi0) igood=1
         if(x5.ge.sn0lo.and.x5.lt.sn0up.and.x4.ge.xlw0.and.
     $        x6.le.chi0) then
            igood=1
         else
            x1=0.
            igood=1
         endif
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
 757     continue
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
      xmult=1.
      print *,ng,nb

      nd=101
      dmin=0.5
      dmax=100.5

c      nw=8
      nw=4
      wmin=3550.
      wmax=5450.

      nall=0
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
            nall=nall+ng+nb
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
c         do i=1,10
c            xv=0.1+0.1*float(i-1)
c            if(i.eq.10) xv=0.95
c            call xlinint(xv,na,ca,da,dv)
c            iwv=nint(wv)
c            write(*,1001) iwv,xv,dv
c         enddo
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

c      print *,n,nall,float(nall)/float(n)
      goto 555

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

c      print *
c      print *,"S/N  N  Avg_err  RMS   avg_F   Ratio_err "
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
c      print *,(s2+s1)/2.,ns,xbout,xsout2,xbout2,xsout2/xbout
      s1=s2
      enddo
      filo=fihi
      enddo

 555  continue

 1001 format(i4,2x,f5.3,1x,f6.2,2(1x,i5))
 1002 format(f6.2,7(1x,f7.4))
 1003 format(a6,7(3x,i4))
 1004 format(a6,7(1x,f7.3))
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      if(xp.lt.x(1)) then
         yp=y(1)
         return
      endif
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.ge.x(n)) yp=y(n)
      return
      end
