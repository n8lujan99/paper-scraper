
      parameter(nmax=20000,nlmax=1000,nca=17,nm2=100000)
      real wave(nmax),flux(nmax),x(nmax),y(nmax),alpha(nca,nca)
      real wl(nmax),a(nca),wl0(nlmax),yin(nmax),covar(nca,nca)
      real fluxe(nmax),ye(nmax),wsky(nmax),fsky(nmax),fskye(nmax)
      real ydata(nmax),ysig(nmax),xd(1036,nm2),ra(nm2),dec(nm2)
      integer naxes(2)
      integer*8 idd(nm2)
      character file1*80,title*60,c1(nm2)*3,nullstr*1
      parameter(pi=3.141592e0,cee=2.99792458e5)
      logical simple,extend,anyf
      data big/1.e20/
      common/aval/ np,iforcew

      wave0=-5007.
c      wavec=50.
c      wavec=100.
      wavec=40.
      wavec=80.
      signsl=2.5
      fluxerr=0.01
      pixsize=2.

      nadd=5
      open(unit=1,file='listin',status='old')
      read(1,*) file1
      close(1)

c      file1="ifuspec_2020.fits"
      iext=5
      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         call ftclos(im1,ier)
         write(*,*) 'Error opening image : ',file1
      endif
      call ftmahd(im1,iext,ihd,ier)
      ier=0
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxis.eq.1) naxes(2)=1
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1036,ncol,nrow,xd,anyf,ier)

      call ftmahd(im1,2,ihd,ier)
      do i=1,nrow
         call ftgcvk(im1,1,i,1,1,0.,idd(i),anyf,ier)
         call ftgcvs(im1,2,i,1,1,nullstr,c1(i),anyf,ier)
         call ftgcve(im1,3,i,1,1,0.,ra(i),anyf,ier)
         call ftgcve(im1,4,i,1,1,0.,dec(i),anyf,ier)
      enddo

      call ftclos(im1,ier)

      iforcew=0
      if(wave0.lt.0) then
         wave0=-wave0
         iforcew=1
      endif
      nl=1
      wl0(1)=wave0
      wl(1)=wave0

      open(unit=11,file='lines.out',status='unknown')
      ido=0
      n=1036
      ntry=36826
      do iall=ntry,ntry
c      do iall=1,nrow
         do i=1,1036
            wave(i)=3470.+float(i-1)*2.
            fluxe(i)=fluxerr
            flux(i)=-666
            if(xd(i,iall).gt.-1e10.and.xd(i,iall).lt.1e10) 
     $           flux(i)=xd(i,iall)
            if(flux(i).eq.0) flux(i)=-666
         enddo

         do il=1,nl
            wlo=wl(il)-wavec
            wup=wl(il)+wavec
            nt=0
            ymin=big
            ymax=-big
            do i=1,n
               if(wave(i).gt.wlo.and.wave(i).lt.wup.and.
     $              nint(flux(i)).ne.-666) then
                  nt=nt+1
                  x(nt)=wave(i)
                  y(nt)=flux(i)
                  yin(nt)=y(nt)
                  ydata(nt)=y(nt)
                  ye(nt)=fluxe(i)
                  ymin=min(ymin,y(nt))
                  ymax=max(ymax,y(nt))
               endif
            enddo
            call biwgt(yin,nt,xb,xs)
            call biwgt(yin,12,xb,xs)

            if(ido.eq.0) then
c               call pgbegin(0,'?',3,3)
               call pgbegin(0,'?',1,1)
               call pgpap(0.,1.)
               call pgscf(2)
               call pgsch(1.4)
               call pgslw(2)
               ido=1
            endif
            call pgenv(wlo,wup,ymin,ymax*1.2,0,0)
            call pgpt(nt,x,y,17)
            do i=1,nt
               call pgslw(1)
               call pgerry(1,x(i),y(i)+ye(i),y(i)-ye(i),1.)
               call pgslw(2)
            enddo
            call pgsch(1.4)
            call pglabel('Wavelength','','')
            call pgsch(1.4)
            amp=(ymax-xb)*3.5
            a(1)=signsl
            a(2)=0.0
            a(3)=0.
            a(4)=xb
            a(5)=0.
            a(6)=wl(il)
            a(7)=amp
            a(8)=a(1)
            a(9)=0.
            na=8
            np=1

            call fitherms(nt,x,y,ye,a,na,covar,alpha,nca)

            h3=a(2)
            h4=a(3)
            con=a(4)
            xoff=a(5)
            rms=0.
            chi=0.
            do ia=1,nt
               yfit=con
               do i=1,np
                  sigg=a(i+nadd+np+np)
                  amp=a(i+nadd+np)
                  vel=a(i+nadd)+xoff
                  w=(x(ia)-vel)/sigg
                  gaus=exp(-w*w/2.)/sqrt(2.*pi*sigg**2)
                  yfit=yfit+amp*gaus*(1.+h3*fh3(w)+h4*fh4(w))
               enddo
               rms=rms+(y(ia)-yfit)**2
               chi=chi+((y(ia)-yfit)/ye(ia))**2
               y(ia)=yfit
               ysig(ia)=yfit-con
            enddo
            rms=sqrt(rms/float(nt))
            chi=chi/float(nt)
            call pgsci(2)
            call pgslw(6)
            call pgline(nt,x,y)
            call pgslw(2)
            call pgsch(1.1)
            call pgsci(1)
            call pgsls(1)
            wfit=a(6)+a(5)
            wfit0=wfit
            nfit=0
            do ii=1,nt
               if(x(ii).gt.(wfit0-6).and.x(ii).lt.(wfit0+6)) nfit=nfit+1
            enddo
         
            znew=(a(6)+a(5))/wl0(il)-1
            zerr=(a(6)+a(5)+sqrt(covar(5,5)))/wl0(il)-1
            zerr=zerr-znew
            iampch=0
            if(covar(7,7).gt.0.and.covar(7,7).lt.1e10) iampch=1
            if(iampch.eq.0) goto 766
            ampe=sqrt(covar(7,7))
            sigg=sqrt(sigg*sigg)
            xnp=4.*sigg
            xnp=xnp/pixsize
c         if(abs(wfit-wl0(il)).gt.20.or.amp.le.0.) goto 766

            if((wfit-x(1)).lt.6) goto 766
            if((x(nt)-wfit).lt.6) goto 766

c - get noise from the rms
            xnoise=rms*sqrt(xnp)
            ston=0.95*amp/xnoise/pixsize

c - get noise from the errors
            w1=wfit-xnp
            w2=wfit+xnp
            nerr=0
            xnoise2=0.
            xmaxn=0.
            signal=0
            do i=1,nt
               if(x(i).ge.w1.and.x(i).le.w2) then
                  nerr=nerr+1
                  xnoise2=xnoise2+ye(i)*ye(i)
                  signal=signal+ysig(i)
                  xmaxn=max(xmaxn,ye(i))
               endif
            enddo
            if(nerr.gt.0) then
               xnoise2=sqrt(xnoise2)
c            xnoise2=xmaxn*sqrt(float(nerr))
            else
               xnoise2=0.
            endif
            ston1=ston
            ston2=0.95*amp/xnoise2/pixsize
c            print *,'Noise_res, Noise_rms, sigma = ',xnoise2,xnoise,sigg
c            print *,'RMS, dAMP, S/N(res,rms) = ',rms,ampe/pixsize,
c     $           ston2,ston1,iall
c         ston=ston2
c            print *,"xnp, signal, noise, sn ",
c     $           nerr,xnp,signal,xnoise2,signal/xnoise2

            write(11,1101) wl0(il),wfit,amp/pixsize,
     $           a(8),sqrt(covar(8,8)),ston1,con,chi,ston2,iall,
     $           idd(iall),c1(iall),ra(iall),dec(iall)
            goto 7661
 766        continue
            write(11,1101) 0.,0.,0.,0.,0.,0.,0.,0.,0.,iall,0,"0",0.,0.
 7661       continue
         enddo
         goto 977
 966     continue
         write(11,1101) 0.,0.,0.,0.,0.,0.,0.,0.,0.,iall,0,"0",0.,0.
 977     continue

      enddo
      close(11)
 1101 format(f8.2,1x,f8.2,1x,f10.4,3(1x,f7.2),2(1x,f9.2),1x,f7.2,1x,i7,
     $     1x,i11,1x,a3,2(1x,f8.4))
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
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end

      subroutine fitherms(n,x,y,sig,a,na,covar,alpha,nca)

      parameter(ncaf=17)

      real x(n),y(n),a(na),covar(nca,nca)
      real alpha(nca,nca),sig(n)
      integer ia(ncaf)
      common/aval/ np,iforcew

      data tol,itermax/1.e-4,1000/

      nadd=5
      do i=1,na
         ia(i)=1
      enddo
      do i=nadd+1,np+nadd
         ia(i)=0
      enddo
      ia(1)=0
      ia(2)=0
      ia(3)=0
      ia(4)=1
      ia(5)=1
c     this is the sigma: 0 fix, 1 fit
      ia(8)=1
      if(iforcew.eq.1) ia(5)=0
      if(iforcew.eq.1) ia(8)=0

      alamda=-1
      alamo=1.e10
      cold=1.e10
      do iter=1,itermax
         call mrqminb(x,y,sig,n,a,ia,na,covar,alpha,nca,
     $        chisq,alamda)
         chirel=abs(cold-chisq)/chisq
         cold=chisq
         if(alamda.lt.alamo.and.chirel.lt.tol) goto 666
         if(alamda.gt.1.e9) goto 666
         alamo=alamda
      enddo
c      print *,'Hit max iteration'
 666  continue

      call mrqminb(x,y,sig,n,a,ia,na,covar,alpha,nca,
     $     chisq,0.)

      return
      end

      subroutine funcs(x,a,yfit,dyda,na)
      real a(na),dyda(na)
      parameter(pi=3.141593e0)
      common/aval/ np,iforcew

      h3=a(2)
      h4=a(3)
      con=a(4)
      xoff=a(5)
      nadd=5
      yfit=con
      do i=1,na
         dyda(i)=0.
      enddo
      dyda(4)=1.
      do i=1,np
         sig=a(i+nadd+np+np)
         amp=a(i+nadd+np)
         vel=a(i+nadd)+xoff
         w=(x-vel)/sig
         gaus=exp(-w*w/2.)/sqrt(2.*pi*sig*sig)
         yadd=amp*gaus*(1.+h3*fh3(w)+h4*fh4(w))
         yfit=yfit+yadd
         dyda(1)=dyda(1)+
     $        (-yadd)/sig+yadd*w*(x-vel)/sig/sig+amp*gaus*(
     $        -h3*dfh3(w)*(x-vel)/sig/sig-h4*dfh4(w)*(x-vel)/sig/sig)
         dyda(2)=dyda(2)+amp*gaus*fh3(w)
         dyda(3)=dyda(3)+amp*gaus*fh4(w)
         dyda(5)=dyda(5)+
     $        yadd*w/sig+amp*gaus*(-h3*dfh3(w)/sig-h4*dfh4(w)/sig)
         dyda(i+nadd+np)=yadd/amp
         dyda(i+nadd)=1.
         dyda(i+nadd+np+np)=
     $        (-yadd)/sig+yadd*w*(x-vel)/sig/sig+amp*gaus*(
     $        -h3*dfh3(w)*(x-vel)/sig/sig-h4*dfh4(w)*(x-vel)/sig/sig)
c         print *,i,sig,dyda(i+nadd+np+np),gaus,w,x,vel
c         read *
      enddo

      return
      end

      function fh3(x)
      fh3=1./sqrt(6.)*(2.*sqrt(2.)*x*x*x-3.*sqrt(2.)*x)
      return
      end
      function fh4(x)
      fh4=1./sqrt(24.)*(4.*x*x*x*x-12.*x*x+3.)
      return
      end
      function dfh3(x)
      dfh3=1./sqrt(6.)*(6.*sqrt(2.)*x*x-3.*sqrt(2.))
      return
      end
      function dfh4(x)
      dfh4=1./sqrt(24.)*(16.*x*x*x-24.*x)
      return
      end

      SUBROUTINE mrqminb(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *alamda)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=100)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,m,mfit
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
        call mrqcofb(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      j=0
      do 14 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).ne.0) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussjb(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcofb(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        j=0
        do 17 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END
      SUBROUTINE mrqcofb(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq)
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      PARAMETER (MMAX=100)
      INTEGER mfit,i,j,k,l,m
      REAL dy,sig2i,wt,ymod,dyda(MMAX)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      do 16 i=1,ndata
        call funcs(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END
      SUBROUTINE gaussjb(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=100)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
              
c                pause 'singular matrix in gaussj'
c                print *,'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
c        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
c        if (a(icol,icol).eq.0.) print *,'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL covar(npc,npc)
      INTEGER i,j,k
      REAL swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
