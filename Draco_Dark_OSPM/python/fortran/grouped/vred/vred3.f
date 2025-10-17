
      parameter(nmax=30000)
c      parameter(nmax=312*1032*112)
      parameter(nfmax=312)
      real arr1(1032,112),arr2(1032,112),xftf1(1032,112)
      real spec(1032*112*nfmax),wave(1032*112*nfmax),specloc(1032,112)
      real wmaster(nmax),xmaster(nmax),waves(1036),xata(1036,312)
      real speca(1032,nfmax*112),wavea(1032,nfmax*112),chia2(nfmax)
      real xftf(1032,nfmax*112),specs(1036,nfmax*112),sback(1032,112)
      real waveata(1036),xata1(1036),specsub(1032,112),wavesub(1032,112)
      real xin(1032*112),yin(1032*112),xp(10),sin(1032),win(1032)
      real xin2(1032*112),yin2(1032*112),yin3(10000),xres(1032,112)
      real xmw(10),wt1(1),xwh(10),xwbin(10),yplot(1036),chia(nfmax)
      integer naxes(2),ifibc(112),ipeaks(10),iwbin(10)
      character camp*2,cspecid*3,cifu*3,cifupos*3
      character camp2*2,cspecid2*3,cifu2*3,cifupos2*3
      character file1*120,cds9(312)*14,cmth*6

      iplot=0

      if(iplot.eq.1) then
         call pgbegin(0,'/xwin',2,2)
         call pgask(.false.)
         call pgpap(0.,1.)
         call pgscf(2)
         call pgsch(1.5)
         call pgslw(2)
      endif

c- set iuse to fit calibration from frame (0), or read in calibration (1)
      open(unit=1,file='vred2.in',status='old')
      read(1,*) imth,iuse
      close(1)
      write(cmth,2001) imth
 2001 format(i6)

      nwt=25
c      nwt=17
      wtmin=-1.0
      wtmax=1.0
      wt0=0
      whalf=8.
      nwt2=11
      wtmin2=-0.4/1000.
      wtmax2=0.4/1000.

c- read in all fits files
      ierc=0
      open(unit=1,file='list',status='old')
      nt=0
      nta=0
      wallmin=1e10
      wallmax=-wallmin
      do iall=1,nfmax
         read(1,*,end=666) file1
         call getfits(file1,8,arr2,cspecid2,camp2,cifu2,cifupos2,chi0,
     $        fchi2,ierc)
         call getfits(file1,11,arr1,cspecid,camp,cifu,cifupos,xdum1,
     $        xdum2,ierc)
         call getfits(file1,12,arr2,cspecid2,camp2,cifu2,cifupos2,xdum1,
     $        xdum2,ierc)
         if(cspecid(2:3).eq."  ") cspecid="00"//cspecid(1:1)
         if(cspecid(3:3).eq." ") cspecid="0"//cspecid(1:2)
         if(cifupos(2:3).eq."  ") cifupos="00"//cifupos(1:1)
         if(cifupos(3:3).eq." ") cifupos="0"//cifupos(1:2)

c- read in amp-to-amp from twilights
         call readcal2(xftf1,cifupos,camp,cspecid,
     $        nata,waveata,xata1,iata,cmth,xres,iuse)
         nta=nta+1
         cds9(nta)=cspecid//"_"//cifupos//"_"//cifu//"_"//camp
         chia(nta)=chi0
         chia2(nta)=fchi2
         do i=1,1032
            jin=1
            do j=1,112
               ntt=(nta-1)*112+j
               ftf0=xftf1(i,j)
c               res0=xres(i,j)
               res0=0.
               wavea(i,ntt)=arr2(i,j)
               if(iata.eq.1) arr1(i,j)=0.
               call xlinint2(wavea(i,ntt),nata,waveata,xata1,ata0,
     $              jin,jout)
               jin=max(1,jout-15)
               ftf1=ftf0
               ata1=ata0
               if(ftf0.eq.0) ftf1=1.
               if(ata0.eq.0) ata1=1.
               speca(i,ntt)=arr1(i,j)*(1.-res0)/ftf1/ata1
               if(speca(i,ntt).gt.-1e10.and.speca(i,ntt).lt.1e10) then
               else
                  speca(i,ntt)=0.
               endif
               if(ftf0.lt.0.3.or.ata0.lt.0.6) speca(i,ntt)=0.
               if(arr1(i,j).eq.0) speca(i,ntt)=0.
               if(arr1(i,j).ne.0) then
                  nt=nt+1
                  wave(nt)=wavea(i,ntt)
                  spec(nt)=speca(i,ntt)
                  wallmin=min(wallmin,wave(nt))
                  wallmax=max(wallmax,wave(nt))
               endif
            enddo
         enddo
      enddo
 666  continue
      close(1)
      print *,"Fibers, Amplifiers: ",ntt,nta

c- sort and compress with a brute force method (sort is too slow)
      wbin=0.1
      wallmin=3470.-wbin/2.
      wallmax=5530.+wbin/2.
      nwbin=nint((wallmax-wallmin)/wbin)-1
      nmaster=0
      do iwb=1,nwbin
         w1=wallmin+float(iwb-1)*wbin
         w2=w1+wbin
         nin=0
         do jwb=1,nt
            if(wave(jwb).gt.w1.and.wave(jwb).le.w2) then
               nin=nin+1
               yin3(nin)=spec(jwb)
            endif
         enddo
         call biwgt(yin3,nin,xb,xs)
         if(nin.eq.0) xb=0.
         nmaster=nmaster+1
         wmaster(nmaster)=(w2+w1)/2.
         xmaster(nmaster)=xb
      enddo
         
c- make master sky
c      call getmaster(nt,wave,spec,nmaster,wmaster,xmaster,1000,iuse,0)

c- find peaks for centering
      call getpeaks(nmaster,wmaster,xmaster,xp,ipeaks,iuse,0)
      call sort(4,xp)
      npeak=4
      xp(1)=4359.
      xp(2)=5084.
      xp(3)=5199.
      xp(4)=5461.

c- make new ftf and a2a
      call getftf(nmaster,wmaster,xmaster,nta,wavea,speca,xftf,
     $     waves,xata,iuse,0)
      if(iuse.eq.0) goto 555

c- iterate

c- subtract master sky
      open(unit=12,file='amp.out',status='unknown')
      write(12,"(a72)") "Spc_slt_iid_am Factor N_c Avg"//
     $     " Scale W0 W1 Nlo Avg_orig chi Frac_c2 Frac0"
      open(unit=13,file='ds9reg.dat',status='unknown')
      do iall=1,nta
         chi0=chia(iall)
         fchi2=chia2(iall)
c- first get combined spectra for an amplifier
         nin=0
         do j=1,112
            ntt=(iall-1)*112+j
            do i=1,1032
               if(speca(i,ntt).ne.0.) then
                  nin=nin+1
                  xin(nin)=wavea(i,ntt)
                  yin(nin)=speca(i,ntt)
               endif
            enddo
         enddo
         call getmaster(nin,xin,yin,nin2,xin2,yin2,21,iuse,0)

c- find the continuum sources
         nin3=0
         do j=1,112
            ntt=(iall-1)*112+j
            nin=0
            jin=1
            do i=1,1032
               if(speca(i,ntt).ne.0.) then
                  nin=nin+1
                  nin3=nin3+1
                  call xlinint2(wavea(i,ntt),nin2,xin2,yin2,wv,jin,jout)
                  jin=jout
                  xin(nin3)=speca(i,ntt)-wv
                  yin(nin)=speca(i,ntt)-wv
               endif
            enddo
            call biwgt(yin,nin,xb,xs)
            yin3(j)=xb
         enddo
         call biwgt(xin,nin3,xb,xs)
         xs=xs/sqrt(900.)
c         xcut=1.5*xs
         xcut=7.0*xs
         ncut=0
         do j=1,112
            ntt=(iall-1)*112+j
            ifibc(j)=0
            if(yin3(j).gt.xcut) ifibc(j)=1
            if(ifibc(j).eq.1) ncut=ncut+1
c            print *,ntt,ifibc(j),yin3(j),xs,xcut
         enddo

         nin=0
         do j=1,112
            if(ifibc(j).eq.0) then
               ntt=(iall-1)*112+j
               do i=1,1032
                  if(speca(i,ntt).ne.0.) then
                     nin=nin+1
                     xin(nin)=wavea(i,ntt)
                     yin(nin)=speca(i,ntt)
                  endif
               enddo
            endif
         enddo
         call getmaster(nin,xin,yin,nin2,xin2,yin2,21,iuse,0)
         wtmina=0.
         wtminb=0.

c         jin=1
c         nin2t=0
c         do i=1,nin2
c            if(xin2(i).gt.3550..and.xin2(i).lt.3650.) then
c              call xlinint2(xin2(i),nmaster,wmaster,xmaster,xv,jin,jout)
c              nin2t=nin2t+1
c              jin=jout
c              xin(nin2t)=yin2(i)/xv
c           endif
c         enddo
c         call biwgt(xin,nin2t,xbb,xsb)
c         cifupos=cds9(iall)(5:7)
c         camp=cds9(iall)(13:14)

         jin=1
         do i=1,nin2
            call xlinint2(xin2(i),nmaster,wmaster,xmaster,xv,jin,jout)
            jin=jout
            xin(i)=yin2(i)/xv
         enddo
         call biwgt(xin,nin2,xb,xs)
c         print *,xbb/xb,xbb,xb,nin2t,nin2,cifupos,camp
         if(ncut.ge.40) xb=1.
         if(xb.le.0.1) xb=1.
         if(nin2.eq.0) goto 888

c         goto 887
         wt0a=wt0
         ntry=0
 890     continue
         if(ntry.eq.8) goto 887
         xmin=1e10
         do iw2=1,nwt2
            wtry2=wtmin2+float(iw2-1)*(wtmax2-wtmin2)/float(nwt2-1)
            do iw=1,nwt
               wtry0=wt0a+wtmin+float(iw-1)*(wtmax-wtmin)/float(nwt-1)
               jin=1
               sum=0.
               sum2=0.
               do i=1,nin2
                  wtry=xin2(i)+wtry0
                  wtry=wtry+wtry2*(4500.-xin2(i))
                  stry=yin2(i)
                  call xlinint2(wtry,nmaster,wmaster,
     $                 xmaster,xv,jin,jout)
                  jin=jout
                  sum=sum+(stry-xv*xb)**2
                  if(wtry.gt.5455.and.wtry.lt.5472) 
     $                 sum2=sum2+(stry-xv*xb)**2
               enddo
               xbw=sum+5.*sum2
               if(xbw.lt.xmin) then
                  xmin=xbw
                  wtmina=wtry0
                  wtminb=wtry2
               endif
            enddo
         enddo
c         print *,cds9(iall)(1:3)," ",ntry,wtmina,wtminb*1000.
         if(wtmina.ge.wtmax+wt0a) then
            wt0a=wtmina
            ntry=ntry+1
            goto 890
         endif
         if(wtmina.le.wtmin+wt0a) then
            wt0a=wtmina
            ntry=ntry+1
            goto 890
         endif
 887     continue

         goto 889
         do i=1,nin2
            xin(i)=float(i)
         enddo
         nwbin=7
         xwbin(1)=3500
         xwbin(2)=3990.
         xwbin(3)=4500.
         xwbin(4)=4855.
         xwbin(5)=4865.
         xwbin(6)=5459.
         xwbin(7)=5470.
         do iwb=1,nwbin
            call xlinint2(xwbin(iwb),nin2,xin2,xin,xini,1,jout)
            iwbin(iwb)=xini
         enddo
         do iwb=1,nwbin-1
            xmw(iwb)=1e10
            ilo=iwbin(iwb)
            ihi=iwbin(iwb+1)
            xwh(iwb)=(xin2(ilo)+xin2(ihi))/2.
         enddo

         do iw=1,nwt
            wtry0=wtmin+float(iw-1)*(wtmax-wtmin)/float(nwt-1)
            do iwb=1,nwbin-1
               ilo=iwbin(iwb)
               ihi=iwbin(iwb+1)
               sum=0.
               jin=1
               do i=ilo,ihi
                  wtry=xin2(i)+wtry0
                  stry=yin2(i)
                  call xlinint2(wtry,nmaster,wmaster,
     $                 xmaster,xv,jin,jout)
                  jin=jout
                  sum=sum+(stry-xv*xb)**2
               enddo
               xbw=sum
               if(xbw.lt.xmw(iwb)) then
                  xmw(iwb)=xbw
                  wt1(iwb)=wtry0
               endif
            enddo
         enddo
c         do iwb=1,nwbin-1
c            print *,cds9(iall)(1:3)," ",iwb,xwh(iwb),wt1(iwb)
c         enddo
 889     continue
            
         goto 888
c- find the wavelength offset for the 4 lines
         xmin=1e10
         do iw2=1,nwt2
            wtry2=wtmin2+float(iw2-1)*(wtmax2-wtmin2)/float(nwt2-1)
         do iw=1,nwt
            wtry0=wtmin+float(iw-1)*(wtmax-wtmin)/float(nwt-1)
            nin=0
            do j=1,112
               ntt=(iall-1)*112+j
               jin=1
               if(ifibc(j).eq.0) then
                  do i=1,1032
                     wtry=wavea(i,ntt)+wtry0
                     wtry=wtry+wtry2*(516.-float(i))
                     stry=speca(i,ntt)
                     if(stry.ne.0.) then
                        do ip=1,npeak
                           if(wtry.gt.(xp(ip)-whalf).and.
     $                          wtry.lt.(xp(ip)+whalf)) then
                              call xlinint2(wtry,nmaster,wmaster,
     $                             xmaster,xv,jin,jout)
                              jin=jout
                              nin=nin+1
                              xin(nin)=(stry-xv*xb)**2
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
            call biwgt(xin,nin,xbw,xsw)
            if(xbw.lt.xmin) then
               xmin=xbw
               wtmina=wtry0
               wtminb=wtry2
            endif
         enddo
         enddo
 888     continue     

c- now subtract
         do j=1,112
            ntt=(iall-1)*112+j
            jin=1
            do i=1,1032
c               if(i.gt.200) then
               if(i.gt.0) then
                  xbuse=xb
               else
                  xbuse=xb*(1.+float(200-i)/199.*0.02)
               endif
c               call xlinint2b(wavea(i,ntt),nwbin,xwh,wt1,wtoff,1,jout)
c               wtry=wavea(i,ntt)+wtoff
               wtry=wavea(i,ntt)+wtmina
               wtry=wtry+wtminb*(4500.-wavea(i,ntt))
               call xlinint2(wtry,nmaster,wmaster,xmaster,
     $              xv,jin,jout)
               wavesub(i,j)=wtry
               if(speca(i,ntt).eq.0.) then
                  specsub(i,j)=0.
                  specloc(i,j)=0.
               else
                  specsub(i,j)=speca(i,ntt)-xv*xbuse
                  specloc(i,j)=speca(i,ntt)
               endif
               jin=jout
            enddo
         enddo
         nin=0
         nzero=0
         nztot=0
         do i=300,800
            do j=10,100
               nztot=nztot+1
c               if(ifibc(j).eq.0.and.specsub(i,j).ne.0.) then
               if(specsub(i,j).ne.0.) then
                  nin=nin+1
                  xin(nin)=specsub(i,j)
                  xin2(nin)=specloc(i,j)
               endif
               if(specsub(i,j).eq.0) nzero=nzero+1
            enddo
         enddo
         call biwgt(xin,nin,xbr,xsr)
         call biwgt(xin2,nin,xbr2,xsr2)

c         call getback(specsub,ifibc,sback,nlow)
         cifupos=cds9(iall)(5:7)
         camp=cds9(iall)(13:14)
         call getback2(specsub,ifibc,camp,cifupos,sback,nlow)
         do j=1,112
            ntt=(iall-1)*112+j
            do i=1,1032
               win(i)=wavesub(i,j)
               if(ncut.lt.40.and.specsub(i,j).ne.0) then
                  sin(i)=specsub(i,j)-sback(i,j)
               else
                  sin(i)=specsub(i,j)
               endif
c               sin(i)=specsub(i,j)
c               sin(i)=sback(i,j)
            enddo
            jin=1
            do i=1,1036
               call xlinint2(waves(i),1032,win,sin,xv,jin,jout)
               specs(i,ntt)=xv
               jin=jout
            enddo
         enddo
         if(iplot.eq.1) then
            do i=1,1036
               nin2=0
               do j=1,112
                  ntt=(iall-1)*112+j
                  if(specs(i,ntt).ne.0) then
                     nin2=nin2+1
                     xin2(nin2)=specs(i,ntt)
                  endif
               enddo
               call biwgt(xin2,nin2,xb0,xs0)
               yplot(i)=xb0
            enddo
            call pgenv(3500.,4500.,-5.,10.,0,0)
            call pgline(1036,waves,yplot)
         endif

c- write out the amplifier info
         xbr=max(-999.,xbr)
         xbr=min(999.,xbr)
         xbr2=max(-9999.,xbr2)
         xbr2=min(9999.,xbr2)
         xsr=min(999.,xsr)
         xb=min(9.99,xb)
         fzero=float(nzero)/float(nztot)
         chi0=min(666.,chi0)
         write(12,1201) cds9(iall),xb,ncut,xbr,xsr,
     $        wtmina,wtminb*1000.,nlow,xbr2,chi0,fchi2,fzero
         ix=10
         iy=iall*112-56
         ibad=0
         if(xb.lt.0.7) ibad=1
         if(xb.gt.1.1) ibad=1
         if(xsr.gt.20.) ibad=1
         if(xsr.lt.1.) ibad=1
         if(xsr.gt.16.and.abs(wtmina).gt.0.9) ibad=1
         if(ibad.eq.0) write(13,1301) ix,iy,cds9(iall)
         if(ibad.eq.1) write(13,1302) ix,iy,cds9(iall)
      enddo
      close(12)
      close(13)

 555  continue
c- write out each a2a from twilight
      if(iuse.eq.0) then
         open(unit=1,file='list',status='old')
         do iall=1,nta
            read(1,*,end=666) file1
            write(file1(26:32),1002) "ata.dat"
            open(unit=2,file=file1(1:32),status='unknown')
            do i=1,1036
               write(2,*) waves(i),xata(i,iall)
            enddo
            close(2)
         enddo
         close(1)
      endif

c- correct the subtracted image to the average offset
      ifix=1
      if(iuse.eq.1.and.ifix.eq.1) then
         open(unit=1,file=
c     $  'sky_res.use'
     $  '/scratch/03261/polonius/science_reductions/alldet/sky_res.use'
     $    ,status='old')
         ns=0
         do i=1,1036
            read(1,*,end=444) x1,x2
            ns=ns+1
            waveata(ns)=x1
            yplot(ns)=x2            
         enddo
 444     continue
         close(1)
         do i=1,1036
            call xlinint(waves(i),nmaster,wmaster,xmaster,sky0)
            call xlinint(waves(i),ns,waveata,yplot,res)
            sfix=res*sky0
            do j=1,nta*112
               specs(i,j)=specs(i,j)-sfix
            enddo
         enddo
      endif

c- write out the sky used
      if(iuse.eq.1) then
         open(unit=12,file='sky.out',status='unknown')
         do i=1,nmaster
            write(12,1202) wmaster(i),xmaster(i)
         enddo
         close(12)
      endif

 1202 format(1x,f7.2,2x,f11.3)

c- write out the sky-subtracted and rectified fits file
      if(iuse.eq.1) then
         naxis=2
         naxes(1)=1036
         naxes(2)=nta*112
         ier=0
         iblock=1
         igc=0
         call ftinit(51,'out.fits',iblock,ier)
         call ftphps(51,-32,naxis,naxes,ier)
         call ftp2de(51,igc,1036,naxes(1),naxes(2),specs,ier)
         call ftclos(51,ier)
      endif

 1001 format(a3)
 1002 format(a7)
 1201 format(a14,1x,f6.3,1x,i3,4(1x,f7.2),1x,i4,1x,f8.2,1x,f8.2,
     $     1x,f7.4,1x,f7.4)
 1301 format("# text("i2,", ",i5,") textangle=90 text={"a14"}")
 1302 format("# text("i2,", ",i5,
     $     ") textangle=90 text={"a14"} color=red")

      end

      subroutine getback(specsub,ifibc,sback,nlow)
      real specsub(1032,112),sback(1032,112),xin(1032*1000)
      real xin1(1032*10),xin2(1032*10),xin3(1032*10)
      real xla(3),xba(3)
      integer ifibc(112),ila(4)

      ifibh=3
      ispech=20
      ilo=300
      ihi=800

      ila(1)=100
      ila(2)=400
      ila(3)=700
      ila(4)=1000
      xla(1)=float(ila(1)+ila(2))/2.
      xla(2)=float(ila(2)+ila(3))/2.
      xla(3)=float(ila(3)+ila(4))/2.

c - get a smoothed background: iback=0 is by fiber, 1 is very local                                                                                                     
      iback=0
      if(iback.eq.0) then
         do j=1,112
            jlo=max(1,j-ifibh)
            jhi=min(112,j+ifibh)
            nin=0
            nin1=0
            nin2=0
            nin3=0
            do jt=jlo,jhi
               if(jt.ne.j) then
                  if(ifibc(jt).eq.0) then
c                     do i=ilo,ihi
c                        nin=nin+1
c                        xin(nin)=specsub(i,jt)
c                     enddo
                     do i=ila(1),ila(2)
                        nin1=nin1+1
                        xin1(nin1)=specsub(i,jt)
                     enddo
                     do i=ila(2),ila(3)
                        nin2=nin2+1
                        xin2(nin2)=specsub(i,jt)
                     enddo
                     do i=ila(3),ila(4)
                        nin3=nin3+1
                        xin3(nin3)=specsub(i,jt)
                     enddo
                  endif
               endif
            enddo
            if(nin1.gt.0.and.nin2.gt.0.and.nin3.gt.0) then
               call biwgt(xin1,nin1,xba(1),xs1)
               call biwgt(xin2,nin2,xba(2),xs2)
               call biwgt(xin3,nin3,xba(3),xs3)
               jin=1
c               print *,j,xba(1),xba(2),xba(3)
               do i=1,1032
                  call xlinint2b(float(i),3,xla,xba,xb0,jin,jout)
                  jin=jout
                  sback(i,j)=xb0
               enddo
            else
c            if(nin.ge.1) then
c               call biwgt(xin,nin,xb,xs)
c               do i=1,1032
c                  sback(i,j)=xb
c               enddo
               do i=1,1032
                  sback(i,j)=0.
               enddo
            endif
         enddo
      else
         do j=1,112
            jlo=max(1,j-ifibh)
            jhi=min(112,j+ifibh)
            do i=1,1032
               ilo2=max(1,i-ispech)
               ihi2=min(1032,i+ispech)
               nin=0
               do jt=jlo,jhi
                  if(ifibc(jt).eq.0) then
                     do it=ilo2,ihi2
                        nin=nin+1
                        xin(nin)=specsub(it,jt)
                     enddo
                  endif
               enddo
               if(nin.ge.1) then
                  call biwgt(xin,nin,xb,xs)
                  sback(i,j)=xb
               else
                  sback(i,j)=0.
               endif
            enddo
         enddo
      endif

c- get the number of fibers below the average value
      nt=0
      n=0
      do j=1,112
         do i=ilo,ihi
            if(specsub(i,j).ne.0) then
               n=n+1
               xin(n)=specsub(i,j)
            endif
         enddo
      enddo
      call biwgt(xin,n,xball,xs)
      cut=-2.*xs
      nlow=0
      do j=1,112
         n=0
         do i=ilo,ihi
            if(specsub(i,j).ne.0) then
               n=n+1
               xin(n)=specsub(i,j)-xball
            endif
         enddo
         call biwgt(xin,n,xb,xs)
         if(n.gt.10.and.xb.lt.cut) nlow=nlow+1
      enddo
c      print *,xball,xs,nlow

      return
      end

      subroutine getback2(xspecs,ifibc,amp,cifupos,xspecn,nlow)
      parameter(narrm=1032)
      real xspecs(narrm,112),xin(1032*1000)
      real xin1(1032*10),xin2(1032*10),xin3(1032*10),xin5i(112)
      real xspecn(narrm,112),xla(3),xba(3),xin4(1032),xin5(112)
      real xicol(10),xcon(10),xslope(10),xfix(1032,112),xntry(10)
      integer ifibc(112)
      integer iloa(10),ihia(10)
      character amp*2,cifupos*3

      nia=5
      iloa(1)=10
      ihia(1)=180
      iloa(2)=180
      ihia(2)=415
      iloa(3)=415
      ihia(3)=615
      iloa(4)=615
      ihia(4)=815
      iloa(5)=815
      ihia(5)=1015

c      xntry(1)=0.75
c      xntry(2)=0.69
c      xntry(3)=0.68
c      xntry(4)=0.666
c      xntry(5)=0.666

c      xntry(1)=0.42
      xntry(1)=0.50
      xntry(2)=0.45
      xntry(3)=0.45
      xntry(4)=0.45
      xntry(5)=0.48
      do ni=1,nia
         ilo=iloa(ni)
         ihi=ihia(ni)
         xicol(ni)=float(ilo+ihi)/2.
      enddo

      if(amp.eq.'LL'.or.amp.eq.'RU') then
         i1a=20
         i2a=112
         i1b=1
         i2b=20
         icuta=1
         icutb=15
      endif
      if(amp.eq.'LU'.or.amp.eq.'RL') then
         i1a=1
         i2a=93
         i1b=93
         i2b=112
         icuta=98
         icutb=112
      endif

      ifixs=0
      do ni=1,nia
         ilo=iloa(ni)
         ihi=ihia(ni)
         nin4=0
         nin5=0
         do j=1,112
            if(ifibc(j).eq.0) then
               nin=0
               do i=ilo,ihi
                  if(xspecs(i,j).ne.0.) then
                     nin=nin+1
                     xin(nin)=xspecs(i,j)
                  endif
               enddo
               if(nin.gt.0) then
                  call biwgt(xin,nin,xb,xs)
                  if(j.ge.i1a.and.j.le.i2a) then
                     nin4=nin4+1
                     xin4(nin4)=xb
                  endif
                  if(j.ge.i1b.and.j.le.i2b) then
                     nin5=nin5+1
                     xin5(nin5)=xb
                     xin5i(nin5)=j
                  endif
               else
                  xb=0.
               endif
            endif
         enddo
         call biwgt(xin4,nin4,xb4,xs4)
         ntry=nint(float(nin4)*xntry(ni))
         ntry=max(ntry,25)
         ntry=min(ntry,112)
         call biwgt(xin4,ntry,xb4,xs4)
         call gety0(nin5,xin5i,xin5,xb4,i1b,i2b,y0)
         if(nin4.eq.0) xb4=0.
         if(nin4.eq.0) y0=0.
         xcon(ni)=xb4
         xslope(ni)=y0
         if(xs4.gt.10..and.nin4.gt.0.) then
            ntry=min(nin4,30)
            call biwgt(xin4,ntry,xb4,xs4)
            xcon(ni)=xb4
            xslope(ni)=xb4
            ifixs=1
         endif
c         print *,xicol(ni),3470+2.*(xicol(ni)-1.),xcon(ni),xslope(ni)
c         print *,cifupos,amp," ",ni,ntry,xcon(ni),xslope(ni)
      enddo

      if(ifixs.eq.1) then
         do i=1,nia
            xslope(i)=xcon(i)
         enddo
      endif

      if(i1b.gt.10) then
         do i=1,1032
            call xlinint(float(i),nia,xicol,xcon,con)
            call xlinint(float(i),nia,xicol,xslope,y0)
            do j=i1a,i2a
               xspecn(i,j)=con
            enddo
            do j=i1b,i2b
               val=con+(y0-con)/float(i2b-i1b)*(float(j)-float(i1b))
               xspecn(i,j)=val
            enddo
         enddo
      else
         do i=1,1032
            call xlinint(float(i),nia,xicol,xcon,con)
            call xlinint(float(i),nia,xicol,xslope,y0)
            do j=i1a,i2a
               xspecn(i,j)=con
            enddo
            do j=i1b,i2b
              val=con+(y0-con)/float(i2b-i1b)*(float(20-i1b+1)-float(j))
               xspecn(i,j)=val
            enddo
         enddo
      endif

c- get the number of fibers below the average value
      ilo=300
      ihi=800
      n=0
      do j=1,112
         do i=ilo,ihi
            if(xspecs(i,j).ne.0) then
               n=n+1
               xin(n)=xspecs(i,j)
            endif
         enddo
      enddo
      call biwgt(xin,n,xball,xs)
      cut=-2.*xs
      nlow=0
      do j=1,112
         n=0
         do i=ilo,ihi
            if(xspecs(i,j).ne.0) then
               n=n+1
               xin(n)=xspecs(i,j)-xball
            endif
         enddo
         call biwgt(xin,n,xb,xs)
         if(n.gt.10.and.xb.lt.cut) nlow=nlow+1
      enddo

      return
      end

      subroutine gety0(n,x,y,xb4,i1b,i2b,y0)
      real x(n),y(n),xin(112)

      if(n.eq.0) then
         y0=0.
         return
      endif

      ny=20
      y0min=-10.
      y0max=10.
      ymin=1.e10
      imax=10
      do iy0=1,ny
         y0=y0min+(y0max-y0min)*float(iy0-1)/float(ny-1)
         nin=0
         do i=1,n
            nin=nin+1
            if(i1b.gt.10.) then
               val=xb4+(y0-xb4)/float(i2b-i1b)*(x(i)-float(i1b))
               xin(nin)=abs(val-y(i))
            else
               val=xb4+(y0-xb4)/float(i2b-i1b)*(float(i2b-i1b+1)-x(i))
               xin(nin)=abs(val-y(i))
            endif
         enddo
         call biwgt(xin,nin,xb,xs)
         if(xb.lt.ymin) then
            ymin=xb
            y0m=y0
         endif
      enddo
      y0=y0m

      return
      end


      subroutine getftf(nmaster,wmaster,xmaster,nta,wavea,speca,
     $     xftf,waves,xata,iuse,iplot)
      parameter(nmax=30000)
c      parameter(nmax=312*1032*112)
      parameter(nfmax=312)
      real wmaster(nmax),xmaster(nmax)
      real speca(1032,nfmax*112),wavea(1032,nfmax*112)
      real xftf(1032,nfmax*112),xin(nmax),yin(nmax)
      real yin2(nmax),yina(1032),ftfa(1036,112)
      real waves(1036),wave1(1032),spec1(1032),xata(1036,312)
      integer iftf(112)

      do i=1,1036
         waves(i)=3470.+float(i-1)*2.
      enddo

      if(iuse.eq.1) return

      do iall=1,nta
         if(iplot.eq.1) call pgenv(3470.,5540.,0.7,1.4,0,0)
         ic=0
         do j=1,112
            jt=(iall-1)*112+j
            iftf(j)=0
            nin=0
            do i=1,1032
               wave1(i)=wavea(i,jt)
               spec1(i)=speca(i,jt)
            enddo

c- first find which fibers to use
            jin=1
            jin2=1
            do i=1,1036
               call xlinint2b(waves(i),nmaster,wmaster,xmaster,xm1,
     $              jin,jout)
               jin=jout
               call xlinint2b(waves(i),1032,wave1,spec1,xs1,jin2,jout2)
               jin2=jout2
               ftfa(i,j)=xs1/xm1
               if(xs1.ne.0) then
                  nin=nin+1
                  xin(nin)=waves(i)
                  yin(nin)=xs1/xm1
               endif
            enddo
            do i=1,nin
               yin2(i)=yin(i)
            enddo
            call biwgt(yin2,nin,xb,xs)
c            print *,j,xb,xs
            if(xb.gt.0.65.and.xb.lt.1.3) then
               iftf(j)=1
               ic=ic+1
               if(ic.eq.14) ic=0
            else
               print *,j,xb,xs
            endif
         enddo

c- then compile to measure a2a
         do i=1,1036
            nin=0
            do j=1,112
               if(iftf(j).eq.1) then
                  nin=nin+1
                  yin(nin)=ftfa(i,j)
               endif
            enddo
            call biwgt(yin,nin,xb,xs)
            yin2(i)=xb
         enddo

         if(iplot.eq.1) then
            call pgslw(5)
            call pgsci(1)
            call pgline(1036,waves,yin2)
         endif

         do i=1,1036
            xata(i,iall)=yin2(i)
         enddo
      enddo

      return
      end

      subroutine getmaster(n,wave,spec,nmaster,wmaster,xmaster,ism,
     $     iuse,iplot)
      parameter(nmax=30000)
c      parameter(nmax=312*1032*112)
      real wave(n),spec(n)
      real wmaster(nmax),xmaster(nmax)

      call sort2(n,wave,spec)
      call smooth3(n,wave,spec,nmaster,wmaster,xmaster,ism)
c      print *,n,nmaster,float(nmaster)/1032.

      xmin=3500.
      xmax=5500.
      xmin=4035.
      xmax=4060.

      if(iplot.eq.1) then
         ymin=1e10
         ymax=-ymin
         do i=1,nmaster
            if(wmaster(i).gt.xmin.and.wmaster(i).lt.xmax) then
               ymin=min(ymin,xmaster(i))
               ymax=max(ymax,xmaster(i))
            endif
         enddo
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgline(nmaster,wmaster,xmaster)
         call pgsci(1)
      endif

      return
      end

      subroutine getfits(file1,iext,arr,cspecid,camp,cifu,cifupos,chi0,
     $     fchi2,ierc)
      real arr(1032,112),xin(112*1032)
      integer naxes(2)
      character file1*120,cname*15,cpos*1,chalf*1
      character camp*2,cspecid*3,cifu*3,cifupos*3
      logical simple,extend,anyf

      im1=0
      ier=0
c      iext=1
c      call ftgiou(im1,ier)
      im1=50
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         if(ierc.eq.0) write(*,*) 'Error opening image : ',file1
         ierc=1
         goto 706
      endif
      ier2=0
      call ftgkys(im1,"SPECID",cspecid,cname,ier2)
      call ftgkys(im1,"CCDPOS",cpos,cname,ier2)
      call ftgkys(im1,"CCDHALF",chalf,cname,ier2)
c      call ftgkys(im1,"AMP",camp,cname,ier2)
      cifu="   "
      call ftgkys(im1,"IFUID",cifu,cname,ier2)
      if(cifu(3:3).eq." ") then
         cifu(2:3)=cifu(1:2)
         cifu(1:1)="0"
      endif
      if(cifu(3:3).eq." ") then
         cifu(2:3)=cifu(1:2)
         cifu(1:1)="0"
      endif
      call ftgkys(im1,"IFUSLOT",cifupos,cname,ier2)
      camp=cpos//chalf

      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1032,ncol,nrow,arr,anyf,ier)

      call ftclos(im1,ier)
      if(iext.eq.8) then
         nin=0
         nc2=0
         do i=100,1000
            do j=10,100
               if(arr(i,j).gt.0) then
                  nin=nin+1
                  xin(nin)=arr(i,j)
                  if(arr(i,j).gt.2.) nc2=nc2+1
               endif
            enddo
         enddo
         call biwgt(xin,nin,xb,xs)
         chi0=666.
         if(nin.gt.0.) chi0=xb
         fchi2=1.
         if(nin.gt.0.) fchi2=float(nc2)/float(nin)
      endif
 706  continue
c      if(ier.ne.0) print *,"No file for: ",file1
      return
      end

      subroutine smooth3(nin,xin,yin,nin2,xin2,yin2,ns)
      parameter(nmax=30000)
c      parameter(nmax=312*1032*112)
      real xin(nin),yin(nin),xin2(nmax),yin2(nmax)
      real yin3(nmax)

      nin2=0
      do i=1,nin,ns
         n=0
         sum1=0.
         jmin=i
         jmax=min(nin,i+ns)
         do j=jmin,jmax
            n=n+1
            sum1=sum1+xin(j)
            yin3(n)=yin(j)
         enddo
         call biwgt(yin3,n,xb,xs)
         ntry=n
         call biwgt(yin3,ntry,xb,xs)
         sum1=sum1/float(n)
         nin2=nin2+1
         xin2(nin2)=sum1
         yin2(nin2)=xb
      enddo

      return
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end

      subroutine xlinint2(xp,n,x,y,yp,jin,jout)
      real x(n),y(n)
c      do j=1,n-1
      do j=jin,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            jout=j
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=0.
      if(xp.gt.x(n)) yp=0.
      jout=1
      return
      end

      subroutine xlinint2b(xp,n,x,y,yp,jin,jout)
      real x(n),y(n)
      do j=jin,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            jout=j
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      jout=1
      return
      end

      subroutine readcal2(xftf,cifupos,camp,cspecid,nata,waveata,xata1,
     $     iata,cmth,xres,iuse)

      real xftf(1032,112),waveata(1036),xata1(1036),xres(1032,112)
      integer naxes(2)
      character file1*100,file2*100
      character cifupos*3,camp*2,cmth*6,cspecid*3
      logical simple,extend,anyf

c- get fiber-to-fiber
c      file1="/work/03946/hetdex/maverick/virus_config/lib_calib/"//cmth
c     $     //"/i"//cifupos//"a"//camp//"cmbf.fits"
      file1="/scratch/projects/hetdex/lib_calib/"//cmth
     $     //"/i"//cifupos//"a"//camp//"cmbf.fits"

      im1=0
      ier=0
      iext=1
c      call ftgiou(im1,ier)
      im1=50
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 707
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1032,ncol,nrow,xftf,anyf,ier)
      call ftclos(im1,ier)

 707  continue
      if(ier.ne.0) print *,"No f2f for: ",file1

c- get amp-to-amp; if creating a2a (iuse=0) then return
      iata=0
      if(iuse.eq.0) return

c      file1="/work/03946/hetdex/maverick/virus_config/lib_calib/"//cmth
c     $     //"/i"//cifupos//"a"//camp//"ata.dat"
      file1="/scratch/projects/hetdex/lib_calib/"//cmth
     $     //"/i"//cifupos//"a"//camp//"ata.dat"

      open(unit=2,file=file1,status='old',err=708)
      aold=1.
      nata=0
      do i=1,1036
         read(2,*,end=711) x1,x2
         nata=nata+1
         waveata(nata)=x1
         if(nata.gt.10) then
            adiff=aold-x2
            if(adiff.gt.0.15) x2=aold
         endif
         if(x2.lt.0.7) x2=aold
         xata1(nata)=x2
         aold=xata1(nata)
      enddo
 711  continue
      close(2)
      goto 709
 708  continue
      print *,"No a2a for: ",file1
      iata=1
      do i=1,1036
         waveata(i)=3470.+float(i-1)*2.
         xata1(i)=1.
      enddo
 709  continue

c- get residual correction
c- skip for now since it is done in vred!!!
      goto 710

      file1="/scratch/projects/hetdex/lib_calib/reschi/"
     $     //"res"//cspecid//camp//".fits"

      im1=0
      ier=0
      iext=1
c      call ftgiou(im1,ier)
      im1=50
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 710
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1032,ncol,nrow,xres,anyf,ier)
      call ftclos(im1,ier)

 710  continue

      return
      end

      subroutine getpeaks(nm,wm,xm,xp,ipeaks,iuse,iplot)
      parameter(nmax=30000)
c      parameter(nmax=312*1032*112)
      real wm(nm),xm(nm),xp(10)
      real xin(nmax),yin(nmax),yin2(nmax)
      integer ipeaks(10)

      if(iuse.eq.0) return

c- find peaks in master                                                                                                                                                     
      npeak=4
      call smooth3(nm,wm,xm,nin,xin,yin,401)
      ymin=1e10
      ymax=-ymin
      do i=1,nm
         call xlinint(wm(i),nin,xin,yin,yv)
         yin2(i)=xm(i)-yv
         ymin=min(ymin,yin2(i))
         ymax=max(ymax,yin2(i))
      enddo

      if(iplot.eq.1) then
         call pgenv(wm(1),wm(nm),ymin,ymax,0,0)
c         call pgenv(wm(1),wm(nm),50.,300.,0,0)
         call pgline(nm,wm,yin2)
c         call pgline(nm,wm,xm)
c         call pgsci(2)
c         call pgline(nin,xin,yin)
c         call pgsci(1)
      endif

      do ip=1,npeak
         ymax=-1e10
         do i=1,nm
            if(ip.gt.1) then
               do ipc=1,ip-1
                  if(abs(xp(ipc)-wm(i)).lt.15) goto 666
               enddo
            endif
            if(yin2(i).gt.ymax) then
               ymax=yin2(i)
               xmax=wm(i)
               imax=i
            endif
 666        continue
         enddo
         xp(ip)=xmax
         ipeaks(ip)=imax
      enddo
      
      return
      end

