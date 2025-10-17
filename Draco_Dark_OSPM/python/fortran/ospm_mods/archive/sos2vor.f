c ********************************************************************
c * SOS2NUCLEI: transforms the s.o.s.-data into file, used by varea  *
c *             also prepares dL-dE plane                            *
c ********************************************************************
      program sos2nuclei
      INCLUDE 'libdefs.h'
      parameter(nmax=10000,nom=Norbit,ninter=40)
      parameter(nsos_max=55,idraw=1,nquadlength=10000,nsos_lim=1800)
      real xmin,xmax,ymin,ymax,fac1
      real x(nmax),y(nmax),t(nmax),xt(nmax),yt(nmax)
      real rs(nom,Nsos),vs(nom,Nsos),ts(nom,Nsos)
      real rperiseq,raposeq
      real xenv(nmax),yenv(nmax)
      real order(nmax),ts_seq(nmax),vs_seq(nmax),rs_seq(nmax)
      real dvs,drs,dvsmin,drsmin
      real rsmin,rsmax,vsmin,vsmax,dist,mindist,dist_lim
      real xshift,yshift,fracshift,minshift,maxshift,tooclose
      integer Ebin(nom),Lzbin(nom),nsosorb(nom),orb(nmax)
      integer Eseq(nom),Lzseq(nom),nseq,norbinseq(nom),seqoforb(nom)
      character cinc*40,lf*4
      integer nc,iml,false(nmax),idum,ntot,nsosmax

      integer iorb,iseq,itake,isosneu,j,nadd,ie0,iLz0,i
      integer iLz,iE,idata,imaxtemp,imaxtemplast,iaktorig
      integer iadd,ienv,i1,i2,ic,isos,imax,iakt,ix
      real xmaxorig,ymaxorig,xtemp,ymaxtemp,ysmall,xdiff,ydiff
      real xm,xn,xakt,xjunk

      real delta_nucx,delta_nucy

      idum = -1

c --- parameters
c     idraw = 1  : choose nsos_max points randomly
c     idraw != 1 : choose points by geometric condition
c     nsos_lim   : maximum allowed number of SOS-imprints per sequence (memory!)
      idum=-1
      tooclose=1.e-8         ! nuclei shouldn't be closer than tooclose
                             ! to avoid problems in varea

      ! for delta_I3
      dist_lim = 0.1          ! for dist < dist_lim, nuclei are mirrored
      fracshift = 0.001        ! fractional shift for edges of envelope
c      minshift = 0.001        ! lower limit for absolute shift of edges
c      maxshift = 0.01         ! upper limits for absolute shift of edges
      minshift = 0.0000001
      maxshift=100000.0
c --- end

c --- get library
      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call dataread()
      call galaxyread()
c      call haloread()
      call spacel()
      cden = dL(1,1)*ratML(1,1) / vol2d(1)
      core = 4.*pi*cden*rmin**3/(3.+cslope)
      call density()
      call tablesread()
      if(iquad.ne.0) call forcefix()

      fac1 = sqrt(GG*totlight*arcsec/distance/angrad)

      open(unit=68,file='norbit',status='old')
      read(68,*) Norb
      close(68)
c --- end


c      call pgbegin(0,'/cps',3,3)
c      call pgsch(1.2)
c      call pgscf(2)


c --- get orbital sos-information and sequences
      call readsos(rs,vs,ts,nsosorb,Ebin,Lzbin,angrad,fac1)
      call sosseq(Ebin,Lzbin,Eseq,Lzseq,nseq,Norb)
c --- end


c ===============================
c =                             =
c = create files for Voronoi    =
c =                             =
c ===============================

      open(unit=7,file='varea.script',status='unknown')

      ! preliminaries
      do iseq=1,nseq
         norbinseq(iseq)=0
         do iorb=1,Norb
            if(Ebin(iorb).eq.Eseq(iseq).and.Lzbin(iorb).eq.Lzseq(iseq)) 
     &           then
               norbinseq(iseq)=norbinseq(iseq)+1
               seqoforb(iorb)=iseq
            end if
         end do
      end do


c --- get nsos_max points out of the orbital SOS ...
      if(idraw.eq.0) then
         ! ... homogeneously in radius
         do iorb=1,Norb
            vsmin=1.e20
            vsmax=-1.e20
            do isos=1,nsosorb(iorb)
               rs_seq(isos)=10**rs(iorb,isos)
               vs_seq(isos)=vs(iorb,isos)
               vsmin=min(vsmin,vs_seq(isos))
               vsmax=max(vsmax,vs_seq(isos))
               ts_seq(isos)=ts(iorb,isos)
               order(isos)=float(isos)
            end do
            do isos=nsosorb(iorb)+1,nmax
               rs_seq(isos)=1.e12
               order(isos)=-99.
            end do
            call sort2(nmax,rs_seq,order)
            rsmin=rs_seq(1)
            rsmax=rs_seq(nsosorb(iorb))
            drsmin=(rsmax-rsmin)/float(Nsos)
            dvsmin=vsmax/float(Nsos)
            if(nsosorb(iorb).lt.60.or.norbinseq(seqoforb(iorb)).lt.4) 
     &           then
               itake=1
            else
               itake=0
            end if
            rs(iorb,1)=log10(rs_seq(1))
            vs(iorb,1)=vs_seq(int(order(1)))
            ts(iorb,1)=ts_seq(int(order(1)))
            isosneu=1
            do isos=2,nsosorb(iorb)
               dvs=abs(vs_seq(int(order(isos)))-vs(iorb,isosneu))
               drs=abs(rs_seq(isos)-10**rs(iorb,isosneu))
               if(itake.eq.1) then
                  isosneu=isosneu+1
                  rs(iorb,isosneu)=log10(rs_seq(isos))
                  vs(iorb,isosneu)=vs_seq(int(order(isos)))
                  ts(iorb,isosneu)=ts_seq(int(order(isos)))
               else
                  if(dvs.gt.dvsmin.or.drs.gt.drsmin) then
                     isosneu=isosneu+1
                     rs(iorb,isosneu)=log10(rs_seq(isos))
                     vs(iorb,isosneu)=vs_seq(int(order(isos)))
                     ts(iorb,isosneu)=ts_seq(int(order(isos)))
                  end if
               end if
            end do
            nsosorb(iorb)=isosneu
         end do
      else
         ! ... randomly choosen
         do iorb=1,Norb
            vsmin=1.e20
            vsmax=-1.e20
            if(nsosorb(iorb).lt.nsos_max.or.
     &           norbinseq(seqoforb(iorb)).lt.4) then
                 ! nothing to do: take all points of the SOS
            else
               if(nsos_max*norbinseq(seqoforb(iorb)).gt.nsos_lim) then
                  nsosmax=nint(float(nsos_lim)/
     &                 float(norbinseq(seqoforb(iorb))))
               else
                  nsosmax=nsos_max
               end if
               do isos=1,nsosorb(iorb)
                  rs_seq(isos)=rs(iorb,isos)
                  vs_seq(isos)=vs(iorb,isos)
                  ts_seq(isos)=ts(iorb,isos)
                  order(isos)=float(isos)
               end do
               do j=1,nsosmax
 9                continue
                  itake=nint(ran1(idum)*float(nsosorb(iorb)-j+1))
                  if(itake.le.0.or.itake.ge.nsosorb(iorb)-j+1) goto 9
                  rs(iorb,j)=rs_seq(itake)
                  vs(iorb,j)=vs_seq(itake)
                  ts(iorb,j)=ts_seq(itake)
                  if(itake.ne.nsosorb(iorb)-j+1) then
                     rs_seq(itake)=rs_seq(nsosorb(iorb)-j+1)
                     vs_seq(itake)=vs_seq(nsosorb(iorb)-j+1)
                     ts_seq(itake)=ts_seq(nsosorb(iorb)-j+1)
                  end if
               end do
               nsosorb(iorb)=nsosmax
            end if
         end do
      end if
c --- end 

      open(unit=29,file='sequences.out',status='unknown')

      do iseq=1,nseq
         ! combine all orbital SOS-imprints for the actual sequence
         ntot=0
         nadd=0
         xmax=-1.e20
         xmin=1.e20
         ymax=-1.e20
         ymin=1.e20
         ie0=Eseq(iseq)
         iLz0=Lzseq(iseq)
         do iorb=1,Norb
            if(Ebin(iorb).eq.ie0.and.Lzbin(iorb).eq.iLz0) then
               do isos=1,nsosorb(iorb)
                  ntot=ntot+1
                  x(ntot)=10**rs(iorb,isos)
                  y(ntot)=vs(iorb,isos)
                  xmax=max(x(ntot),xmax)
                  ymax=max(y(ntot),ymax)
                  xmin=min(x(ntot),xmin)
                  ymin=min(y(ntot),ymin)
                  t(ntot)=ts(iorb,isos)
                  orb(ntot)=iorb
               end do
               iLz = Lzbin(iorb)
               iE = Ebin(iorb)
            end if
         end do

         write(29,*) iseq,iLz,iE

         do i=1,ntot
            x(i)=x(i)/xmax
            y(i)=y(i)/ymax
         end do

         xmaxorig=xmax
         ymaxorig=ymax

         xmax=xmax/xmax
         xmin=xmin/xmax

         idata=ntot
         ! create an envelope for the sequence
         do i=1,idata
            if(y(i).eq.1.0) imax=i
         end do

         nadd = nadd + 1
         yshift = max(minshift,fracshift*y(imax))
         yshift = min(yshift,maxshift)
         xt(nadd)=x(imax)
         yt(nadd)=y(imax)+yshift

         ! from ymax to greater radii
         imaxtemp=imax
         imaxtemplast=-1
         do j=1,idata
            xtemp=x(imaxtemp)
            ymaxtemp=-1.e20
            do i=1,idata
               if(x(i).gt.xtemp) then
                  if(y(i).gt.ymaxtemp) then
                     imaxtemp=i
                     ymaxtemp=y(i)
                  end if
               end if
            end do
            if(imaxtemp.ne.imaxtemplast) then
               nadd = nadd + 1
               xshift=max(minshift,fracshift*x(imaxtemp))
               yshift=max(minshift,fracshift*y(imaxtemp))
               xshift=min(xshift,maxshift)
               yshift=min(yshift,maxshift)
               xt(nadd)=x(imaxtemp)+xshift
               yt(nadd)=y(imaxtemp)+yshift
            end if
            imaxtemplast=imaxtemp
            if(xtemp.eq.xmax) goto 909
         end do
         
 909     continue

         ! from ymax to smaller radii
         imaxtemp=imax
         imaxtemplast=-1
         do j=1,idata
            xtemp=x(imaxtemp)
            ymaxtemp=-1.e20
            do i=1,idata
               if(x(i).lt.xtemp) then
                  if(y(i).gt.ymaxtemp) then
                     imaxtemp=i
                     ymaxtemp=y(i)
                  end if
               end if
            end do
            if(imaxtemp.ne.imaxtemplast) then
               nadd = nadd + 1
               xt(nadd)=x(imaxtemp)*(1.-fracshift)
               yt(nadd)=y(imaxtemp)*(1.+fracshift)
            end if
            imaxtemplast=imaxtemp
            if(xtemp.eq.xmin) goto 910
         end do
         
 910     continue

         ! interpolate between the edges
         do i=nadd+1,nmax
            xt(i)=1.e12
            yt(i)=1.e12
         end do
         call sort2(nmax,xt,yt)
         iakt=nadd

         if(yt(1).ne.0.0) then
            nadd = nadd + 1
            xt(nadd) = xt(1)*0.9999
            yt(nadd) = 0.0
         end if
         if(yt(iakt).ne.0.0) then
            nadd = nadd + 1
            xt(nadd) = xt(iakt)*1.0001
            yt(nadd) = 0.0
         end if
         if(nadd.ne.iakt) then
            do i=nadd+1,nmax
               xt(i)=1.e12
               yt(i)=1.e12
            end do
            call sort2(nmax,xt,yt)
         end if

         xmin=xt(1)
         xmax=xt(nadd)

         ysmall=min((xt(nadd)-xt(1))/float(ninter),1.e-2)
         ysmall=max(ysmall,1.e-3)
         iaktorig=nadd

 1001    continue

         iakt=nadd
         
         do i=2,iakt
            ydiff=abs(yt(i)-yt(i-1))
            xdiff=abs(xt(i)-xt(i-1))
            ydiff=max(xdiff,ydiff)
            if(ydiff.gt.ysmall) then
               iadd=int(ydiff/ysmall)
               do j=1,iadd
                  xm=(yt(i)-yt(i-1))/(xt(i)-xt(i-1))
                  xn=yt(i-1)
                  xakt = xt(i-1)+(xt(i)-xt(i-1))/float(iadd+1)*float(j)
                  nadd = nadd +1 
                  xt(nadd) = xakt
                  yt(nadd) = xm * (xakt-xt(i-1)) + xn
               end do
            end if
         end do

         ienv = 0
         ! mirror points about vr=0
         do i=1,idata
            if(y(i).lt.dist_lim) then
               ienv = ienv + 1
               xenv(ienv)=1.001*x(i) ! ... to avoid problems with varea
               yenv(ienv)=-y(i)
            end if
         end do

         do i=1,nadd
            ienv=ienv+1
            xenv(ienv)=xt(i)
            yenv(ienv)=yt(i)
         end do

         if(nadd.gt.750.and.ysmall.lt.1000.0) then
            ysmall = ysmall*1.1
            nadd = iaktorig
            goto 1001
         end if

         do i=1,ienv
            ntot=ntot+1
            x(ntot)=xenv(i)
            y(ntot)=yenv(i)
c            x(ntot)=xenv(i)+fracshift*xenv(i)*gasdev(idum)
c            y(ntot)=yenv(i)+fracshift*yenv(i)*gasdev(idum)
            t(ntot)=0.
            orb(ntot)=0 ! orbit-number=0 identifies the envelope-points
         end do

         write(*,'("iseq=",i3,3x,"iE=",i3,3x,"iLz=",i3,
     &        3x,"sos=",i4,2x,"bound=",i4,2x,"tot=",i4,a1,$)') 
     &        iseq,ie0,iLz0,idata,ntot-idata,ntot,char(13)


         xmax=-1.e20
         xmin=1.e20
         ymax=-1.e20
         ymin=1.e20
         do i=1,ntot
            if(x(i).gt.xmax) xmax=x(i)
            if(y(i).gt.ymax) ymax=y(i)
            if(x(i).lt.xmin) xmin=x(i)
            if(y(i).lt.ymin) ymin=y(i)
         end do
         !plot SOS and envelope
         xmin=xmin-(xmax-xmin)/10.
         ymin=ymin-(ymax-ymin)/10.
         xmax=xmax+(xmax-xmin)/10.
         ymax=ymax+(ymax-ymin)/10.
c         call pgenv(xmin,xmax,ymin,ymax,0,0)
c         call pgpoint(idata,x,y,17)
c         call pgsci(2)
c         call pgslw(1)
c         call pgpoint(ienv,xenv,yenv,17)
c         call pgsci(1)
c
c         iml=nint(abs(iseq*1.))
c         call pgnumb(iml,0,0,cinc,nc)
c         call puttext(0.95,0.95,1.,'iseq = '//cinc(1:nc),
c     &        xmin,xmax,ymin,ymax)
c
c         call pglabel('x/xmax','y/ymax','') 
         
         if(iseq.lt.10) write(lf,'("000",i1)') iseq
         if(iseq.ge.10.and.iseq.lt.100) write(lf,'("00",i2)') iseq
         if(iseq.ge.100.and.iseq.lt.1000) write(lf,'("0",i3)') iseq
         if(iseq.ge.1000) write(lf,'(i4)') iseq
         
         ! final normalization
         xmax=-1.e20
         xmin=1.e20
         ymax=-1.e20
         ymin=1.e20
         do i=1,ntot
            x(i)=x(i)*xmaxorig
            y(i)=y(i)*ymaxorig
            if(x(i).gt.xmax) xmax=x(i)
            if(y(i).gt.ymax) ymax=y(i)
            if(x(i).lt.xmin) xmin=x(i)
            if(y(i).lt.ymin) ymin=y(i)
         end do

         xmax=xmax*1.25
         ymax=ymax*(1.2+ylim)

         xmax=xmax/float(nquadlength)
         ymax=ymax/float(nquadlength)

         do i=1,ntot
            x(i)=x(i)/xmax
            y(i)=(y(i)+1.1*abs(ymin))/ymax
         end do
         do i=1,ntot
            y(i)=y(i)+0.025
         end do
         xmaxorig=xmax
         ymaxorig=ymax

         delta_nucx = float(nquadlength)
         delta_nucy = float(nquadlength)
c         delta_nucx = (xmax/xmax - xmin/xmax)
c         delta_nucy = ((ymax + 1.1*abs(ymin))/ymax - 
c     &        (ymin+1.1*abs(ymin))/ymax)

         ! check if points fall beyond the edges 
         do i1=1,ntot
            false(i1)=0
            if(x(i1).ge.float(nquadlength)) false(i1)=1
            if(y(i1).ge.float(nquadlength)) false(i1)=1
            if(x(i1).le.0.0) false(i1)=1
            if(y(i1).le.0.0) false(i1)=1
         end do

         ! check if coordinates are real
         do i1=1,ntot
            ifx=1
            ify=1
            if(x(i1).gt.-1.e20.and.x(i1).lt.1.e20) ifx=0
            if(y(i1).gt.-1.e20.and.y(i1).lt.1.e20) ify=0
            ift=ifx+ify
            if(ift.ne.0) false(i1)=1 
         end do

         ! check if two points are too close to each other
         do i1=1,ntot
            do i2=i1+1,ntot
c               if(abs(x(i1)-x(i2))/x(i2).lt.tooclose) then
c                  if(abs(y(i1)-y(i2))/y(i2).lt.tooclose) then 
c                     false(i1)=1
c                  end if
c               end if
               if(abs(x(i1)-x(i2))/delta_nucx.lt.tooclose) then
                  if(abs(x(i1)-x(i2)).eq.0.0) then
                     false(i1)=1
                  else
                     if(abs(y(i1)-y(i2))/delta_nucy.lt.tooclose) then 
                        false(i1)=1
                     end if
                  end if
               end if
            end do
         end do
         
         open(unit=25,file=lf//'.nuclei',status='unknown')
         open(unit=27,file=lf//'.vor.cells',status='unknown')
         open(unit=26,file=lf//'.vor.sos',status='unknown')
         open(unit=28,file=lf//'.vor.plot',status='unknown')
         ic=0
         do i=1,ntot
            xjunk=0.30
            if(y(i).eq.0.0) y(i)=1.e-5
            if(false(i).ne.1) then
               ic=ic+1
               write(25,*) x(i),y(i),xjunk,ic-1
               write(26,'(2(f10.4,1x),(f5.2,1x),(i5,1x),(f12.6,1x),
     &              2(f13.6,1x))') x(i),y(i),xjunk,orb(i),t(i),
     &              xmaxorig,ymaxorig
               write(27,'("1 1 1 -1.00000")')
               write(28,*) x(i),y(i),orb(i),ymin
            end if
         end do
         close(26)
         close(25)
         close(27)
         close(28)
         write(7,'("./varea ",a," ",i6," 0")') lf,nquadlength

 9001    continue

      end do
      close(7)
      close(29)

      open(unit=8,file='nseq.out',status='unknown')
      write(8,*) nseq
      close(8)
c --- end 
         
c      call pgend

      print*

      end





      subroutine readsos(rs,vs,ts,nsosorb,E_bin,Lz_bin,rmax,fac1)
      INCLUDE 'libdefs.h'
      parameter(nom=Norbit)
      real rs(nom,Nsos),vs(nom,Nsos),ts(nom,Nsos)
      real rmax,d_Lz,d_E,rs_in,rvs,ts_in,mlnorm
      integer E_bin(nom),Lz_bin(nom),nsosorb(nom),iE,iLz
      integer iorblast,nseq,ntot,iorb

      open(unit=21,file="gden.norm",status='old')
      read(21,*) mlnorm
      close(21)

      open(unit=43,file='sos1.out',status='old')
      open(unit=44,file='sos2.out',status='old')
      ntot=0
      nseq=0
      iorblast=1

 700  continue
      read(44,*,end=777) iorb,iE,iLz,d_Lz,d_E
      read(43,*) rs_in,vs_in,ts_in
      if(iorb.ne.iorblast) then
         ntot=0
      end if
      ntot=ntot+1
      rs(iorb,ntot)=log10(rs_in*rmax)
      vs(iorb,ntot)=vs_in*fac1*mlnorm
      ts(iorb,ntot)=ts_in
      nsosorb(iorb)=ntot
      E_bin(iorb)=iE
      Lz_bin(iorb)=iLz
      iorblast=iorb
      goto 700

 777  continue
      close(43)
      close(44)
      return
      end



      subroutine sosseq(E_bin,Lz_bin,E_seq,Lz_seq,nseq,Norb)
      INCLUDE 'libdefs.h'
      parameter(nom=Norbit)
      integer E_bin(nom),Lz_bin(nom)
      integer E_seq(nom),Lz_seq(nom),nseq
      integer iLzlast,iElast,iLz,iE,ntot,iorb
      real d_Lz,d_E,xLz,xE,x3,x4
      open(unit=44,file='sos2.out',status='old')
      read(44,*) iorb,iE,iLz,d_Lz,d_E
      iLzlast=iLz
      iElast=iE
      close(44)
      nseq=0
      open(unit=40,file='integrals.out',status='old')
      ntot=0
      do i=1,Norb
         if(E_bin(i).eq.iElast.and.Lz_bin(i).eq.iLzlast) then 
            read(40,*) xLz,xE,x3,x4
         else
            nseq=nseq+1
            E_seq(nseq)=iElast
            Lz_seq(nseq)=iLzlast
            read(40,*) xLz,xE,x3,x4
            iLzlast=Lz_bin(i)
            iElast=E_bin(i)
         end if
      end do
      nseq=nseq+1
      E_seq(nseq)=iElast
      Lz_seq(nseq)=iLzlast
      close(40)
      return
      end


      subroutine plotELz()
      real E(50000),Lz(50000)
      integer n_int
      real x1,x2,x3,x4
      n_int=0
      open(unit=11,file='integrals.out',status='old')

 900  continue
      read(11,*,end=991) x1,x2,x3,x4
      E(2*i-1)=log10(-x2)
      E(2*i)=log10(-x2)
      Lz(2*i-1)=-x1
      Lz(2*i)=x1
      n_int=n_int+2
      goto 900

 991  continue
      close(11)
c      call gminmax(n_int,E,xmin,xmax)
c      call gminmax(n_int,Lz,ymin,ymax)
c      call pgenv(xmin,xmax,ymin,ymax,0,10)
c      call pgpoint(n_int,E,Lz,17)
c      call pglabel('log(-E)','L\\Dz\\U','')
      return
      end
