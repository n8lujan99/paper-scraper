
      INCLUDE 'moddefs.h'
      real y1(Nvel,Nvelbm),sumadph(Nvel),sigm(Nvelbm),sigd(Nvelbm)
      real y1p(Nvel),sumadp(Nvel),sumadpl(Nvel),rbinl(Nvelbm)
      real sigdh(Nvelbm),sigdl(Nvelbm),xl(2),yl(2),chip(Nvelbm)
      real rani(Nrani),vr(Nrani),vt(Nrani),vp(Nrani),vrot(Nrani)
      real veldp(Nvelbm),y12(Nvel,Nvelbm),y1p2(Nvel),velnsc(Nvelbm)
      real af(5),h3m(Nvelbm),h4m(Nvelbm),h3d(Nvelbm),h4d(Nvelbm)
      real rlmaj(Nvelbm),rlmin(Nvelbm),h3mmin(Nvelbm),h3mmaj(Nvelbm)
      real h3dmin(Nvelbm),h3dmaj(Nvelbm),h4mmin(Nvelbm),h4mmaj(Nvelbm)
      real h4dmin(Nvelbm),h4dmaj(Nvelbm),velmp(Nvelbm),signsc(Nvelbm)
      real rnsc(Nvelbm),vrt(Nrani),tani(Ntani)
      character filen(Nveld)*40,galname*40,cinc*40,cadd*2
      character clabel(20)*3

      ifitherm=1

      clabel(1)='mj'
      clabel(2)='a2'
      clabel(3)='a3'
      clabel(4)='a4'
      clabel(5)='a5'
      clabel(6)='a6'
      clabel(7)='a7'
      clabel(8)='a8'
      clabel(9)='a9'
      clabel(10)='a10'

      open(unit=11,file='ml.out',status='old')
      read(11,*) xmu
      read(11,*) alp,ent,chi
      read(11,*) x1,x2,ent,chi
      close(11)
      open(unit=11,file='ml2.out',status='unknown')
      write(11,*) xmu,chi,x1,x2
      close(11)
      xmu=1./sqrt(xmu)
      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call totmlread()
      call galreadm(galname)
      call phaseread()
      call velbin(filen)
      call smcoarse()
      call vlibread()
      call vdataread(filen)
      call weightread()
      call vaniread()
      call velmtod()

C---- fac1 is the factor needed to convert velocities to km/s
      fac1 = sqrt(GG*totlight*arcsec/distance/angrad)
   
      call pgbegin(0,'?',2,2)
      call pgsch(1.2)
      call pgscf(2)
      print *

c - internal moments

c - first get the radii and angles in the bin centers

      open(unit=11,file='bin_r.out',status='old')
      read(11,*)
      read(11,*)
      do i=1,Nrani
         read(11,*) i1,x2,x3
         rani(i)=log10(x3)
      enddo
      close(11)
      open(unit=11,file='bin_v.out',status='old')
      read(11,*)
      read(11,*)
      do i=1,Ntani
         read(11,*) i1,x2,x3
         tani(i)=x3
      enddo
      close(11)

      open(unit=11,file='intmom.out',status='unknown')
      write(11,*) '    R     Theta      V_r      V_theta    ',
     $     'V_rV_t      V_phi     v_phi     beta'

      ymin=big
      ymax=-big
      ymax2=-big
      do iv=1,Ntani
      do ir=1,Nrani
         sum1=0.
         sum2=0.
         sum3=0.
         sum4=0.
         sum5=0.
         sum6=0.
         do iorb=1,Norbit
            sum1=sum1+w(iorb)*vr2a(ir,iv,iorb)
            sum2=sum2+w(iorb)*vt2a(ir,iv,iorb)
            sum3=sum3+w(iorb)*vp1a(ir,iv,iorb)
            sum4=sum4+w(iorb)*vp2a(ir,iv,iorb)
            sum5=sum5+w(iorb)*v0a(ir,iv,iorb)
            sum6=sum6+w(iorb)*vrta(ir,iv,iorb)
         enddo
         vrp=sqrt(sum1/sum5)*fac1/xmu
         vtp=sqrt(sum2/sum5)*fac1/xmu
         vrtp=sum6/sum5*fac1/xmu*fac1/xmu
         vp1=sqrt(sum4/sum5-(sum3/sum5)**2)*fac1/xmu
         vp2=sqrt(sum4/sum5)*fac1/xmu
         vpp=vp1
         vrotp=abs(sum3/sum5*fac1)/xmu
         vtang=sqrt((vtp**2+vp2**2)/2.)
         beta=1.-vtang**2/vrp**2
         if(iv.eq.1) then
            vr(ir)=vrp
            vt(ir)=vtp
            vrt(ir)=vrtp
            vp(ir)=vpp
            vrot(ir)=vrotp
            write(*,1102) ir,10**rani(ir),vr(ir),vt(ir),
     $           vrt(ir),vp(ir),vrot(ir)
            ymax=max(ymax,vr(ir),vt(ir),vp(ir))
            ymax2=max(ymax2,vrot(ir))
         endif
         vrtp=0.
         write(11,1101) 10**rani(ir),tani(iv),vrp,vtp,vrtp,vpp,vrotp,
     $        beta
      enddo
      enddo
      close(11)
 1101 format(2(1x,f8.3),6(1x,f10.4))
 1102 format(1(1x,i3),1x,f7.3,5(1x,f11.3))

      call pgpage

      do ivel=1,Nvel
         velm(ivel)=velm(ivel)/xmu
      enddo

      call d3coarse()
      open(unit=21,file='moddat2.out',status='unknown')
      do i=1,Nbin
         write(21,*) i,d3c(i)
      enddo
      close(21)

      open(unit=21,file='moddat.out',status='unknown')
      do ivbin=1,Nvelb
         do ivel=1,Nvel
            y1(ivel,ivbin)=0.
            y12(ivel,ivbin)=0.
            do iorb=1,Norbit
               y1(ivel,ivbin)=y1(ivel,ivbin)+
     $              w(iorb)*v1lib(ivel,ivbin,iorb)
               y12(ivel,ivbin)=y12(ivel,ivbin)+
     $              w(iorb)*v1lib2(ivel,ivbin,iorb)
               write(21,*) ivbin,ivel,iorb,sumad(ivel,ivbin),
     $              sadfer(ivel,ivbin),w(iorb),v1lib(ivel,ivbin,iorb)
            enddo
         enddo
      enddo
      close(21)

      ymaxpd=0.
      yminpd=1.e10
      nmj=0
      do ivbin=1,Nvelb
         do ivel=1,Nvel
            y1p(ivel)=y1(ivel,ivbin)
            y1p2(ivel)=y12(ivel,ivbin)
            sumadp(ivel)=sumad(ivel,ivbin)
            sumadph(ivel)=sumadp(ivel)+sadfer(ivel,ivbin)
         enddo
         call getfwhm(Nvel,velm,y1p,0.5,fwm,xd1,xd2)
         sigm(ivbin)=fwm/2.35
         if(ifitherm.eq.1) then
            af(1)=xd2*sqrt(2.*pi)*fwm/2.35
            af(2)=xd1
            af(3)=fwm/2.35
            af(4)=0.
            af(5)=0.
            call fitherm(Nvel,velm,y1p,af,5)
            velmp(ivbin)=af(2)
            sigm(ivbin)=af(3)
            sigtot=sqrt(sigm(ivbin)**2+velmp(ivbin)**2)
            h3m(ivbin)=af(4)
            h4m(ivbin)=af(5)
         endif
         call getfwhm(Nvel,velm,y1p2,0.5,fwm,xd1,xd2)
         velnsc(ivbin)=xd1
         if(ifitherm.eq.1) then
            af(1)=xd2*sqrt(2.*pi)*fwm/2.35
            af(2)=xd1
            af(3)=fwm/2.35
            af(4)=0.
            af(5)=0.
            call fitherm(Nvel,velm,y1p2,af,5)
            if(iang(ivbin).eq.1.or.rbin(ivbin).eq.0) then
               nmj=nmj+1
               velnsc(nmj)=af(2)
               signsc(nmj)=af(3)
               if(rbin(ivbin).eq.0) then
                  rnsc(nmj)=rani(1)
               else
                  rnsc(nmj)=log10(rbin(ivbin)*angrad)
               endif
c               print *,10**rnsc(nmj),af(2),af(3)
            endif               
         endif
         call getfwhm(Nvel,velm,sumadp,0.5,fwm,xd1,xd2)
         sigd(ivbin)=fwm/2.35
         veldp(ivbin)=xd1
         if(ifitherm.eq.1) then
            af(1)=xd2*sqrt(2.*pi)*fwm/2.35
            af(2)=xd1
            af(3)=fwm/2.35
            af(4)=0.
            af(5)=0.
            call fitherm(Nvel,velm,sumadp,af,5)
            veldp(ivbin)=af(2)
            sigd(ivbin)=af(3)
            sigdtot=sqrt(sigd(ivbin)**2+veldp(ivbin)**2)
            h3d(ivbin)=af(4)
            h4d(ivbin)=af(5)
         endif
         call getfwhm(Nvel,velm,sumadph,0.5,fwm,xd1,xd2)
         diff=fwm/2.35-sigd(ivbin)
         if(ivbin.gt.1) diff=diff*10.
         sigdh(ivbin)=sigd(ivbin)+diff
         sigdl(ivbin)=sigd(ivbin)-diff
         ymaxpd=max(ymaxpd,sigm(ivbin),sigd(ivbin),velmp(ivbin),
     $        veldp(ivbin))
         yminpd=min(yminpd,veldp(ivbin))
         if(rbin(ivbin).eq.0) then
            rbinl(ivbin)=rani(1)
         else
            rbinl(ivbin)=log10(rbin(ivbin)*angrad)
         endif
         if(rbinl(ivbin).lt.rani(1)) rbinl(ivbin)=rani(1)
         write(*,2001) rbin(ivbin)*angrad,iang(ivbin),sigtot,sigdtot
      enddo
 2001 format(1x,f6.3,1x,i2,1x,f6.2,1x,f6.2)

c - first plot the projected dispersions

      rbit=(rani(Nrani)-rani(1))/20.
      yminpd=min(0.,yminpd)
      call pgenv(rani(1)-rbit,rani(Nrani),yminpd,1.05*ymaxpd,0,10)
      call pgsch(1.9)
      do i=1,Nvelb
         call pgsci(iang(i))
         if(isee(i).eq.1) call pgsci(8)
         call pgpt1(rbinl(i),sigm(i),17)
         call pgpt1(rbinl(i),sigd(i),23)
         call pgsci(12)
         if(isee(i).eq.1) call pgsci(8)
         call pgpt1(rbinl(i),velmp(i),13)
         call pgpt1(rbinl(i),veldp(i),7)
      enddo
      call pgsci(1)
      call pgsch(1.2)
      call pglabel('R','\\gs\\Dp\\U (km/s)','')

c - plot the 3d dispersions

      call pgenv(rani(1)-rbit,rani(Nrani),0.,1.05*ymax,0,10)

c - labels start
      call pgsch(2.5)
      call pgmtxt('T',15.,0.3,0.,galname)
      call pgsch(2.0)
      iml=nint(100./xmu**2)
      call pgnumb(iml,-2,0,cinc,nc)
      call pgmtxt('T',16.,0.3,0.,'M/L = '//cinc(1:nc))
      inclin=nint(xinclin*180./pi*10.)
      call pgnumb(inclin,-1,1,cinc,nc)
      call pgmtxt('T',14.,0.3,0.,'inclination = '//cinc(1:nc)//
     $     '\\Uo\\D')
      if(hole.gt.0.) then
         hole=hole/xmu**2
         ip=nint(log10(hole))
         ihole=nint(hole*10./10.**ip)
      else
         ip=0
         ihole=0
      endif
      call pgnumb(ihole,-1+ip,0,cinc,nc)
      call pgmtxt('T',12.,0.3,0.,'BH mass = '//cinc(1:nc))
      ip=nint(log10(abs(ent)))
      ient=nint(10.*ent/10.**ip)
      call pgnumb(ient,-1+ip,0,cinc,nc)
      call pgmtxt('T',8.,0.3,0.,'Entropy = '//cinc(1:nc))
      if(alp.gt.0) then
         ip=nint(log10(chi))
         ichi=nint(chi*1000./10.**ip)
         call pgnumb(ichi,-3+ip,0,cinc,nc)
         call pgmtxt('T',10.,0.3,0.,'\\Gx\\U2\\D = '//cinc(1:nc))
         ip=nint(log10(alp))
         ialp=nint(alp*100./10.**ip)
         call pgnumb(ialp,-2+ip,0,cinc,nc)
         call pgmtxt('T',6.,0.3,0.,'\\ga = '//cinc(1:nc))
      endif
c - end

      call pgsch(1.2)
      call pgline(Nrani,rani,vr)
      call pgpoint(Nrani,rani,vr,17)
      call pgsls(2)
      call pgline(Nrani,rani,vt)
      call pgsls(4)
      call pgline(Nrani,rani,vp)
      call pgsls(1)
      call pglabel('r','\\gs (km/s)','')

      call pgsch(1.5)
      xl(1)=rani(1)+0.8*(rani(Nrani)-rani(1))
      xl(2)=xl(1)+0.08*(rani(Nrani)-rani(1))
      yl(1)=ymax*0.98
      yl(2)=yl(1)
      xpc=xl(1)-0.5*(xl(2)-xl(1))
      ypc=yl(1)*0.99
      call pgptxt(xpc,ypc,0.,0.,'r ')
      call pgline(2,xl,yl)
      call pgsls(1)
      yl(1)=ymax*0.915
      yl(2)=yl(1)
      ypc=yl(1)*0.99
      call pgptxt(xpc,ypc,0.,0.,'\\(2185) ')
      call pgsls(2)
      call pgline(2,xl,yl)
      call pgsls(1)
      yl(1)=ymax*0.85
      yl(2)=yl(1)
      ypc=yl(1)*0.99
      call pgptxt(xpc,ypc,0.,0.,'\\(2186) ')
      call pgsls(4)
      call pgline(2,xl,yl)
      call pgsls(1)

      call pgsch(1.2)
      call pgenv(rani(1)-rbit,rani(Nrani),0.,1.05*ymax2,0,10)
      call pgline(Nrani,rani,vrot)
      call pglabel('r','V\\D\\gf\\U (km/s)','')

c -- plot the individual losvd's

      open(unit=12,file='losvd.out',status='unknown')
      vbit=(velm(Nvel)-velm(1))/15.
      schi2=0.
      chimin=big
      chimax=-big
      do ivbin=1,Nvelb
         ymax=0.
         ymin=big
         schi=0.
         do ivel=1,Nvel
            y1p(ivel)=y1(ivel,ivbin)
            y1p2(ivel)=y12(ivel,ivbin)
            sumadp(ivel)=sumad(ivel,ivbin)
c            sumadpl(ivel)=max(0.,sumadp(ivel)-sadfer(ivel,ivbin))
            sumadpl(ivel)=sumadp(ivel)-sadfer(ivel,ivbin)
            sumadph(ivel)=sumadp(ivel)+sadfer(ivel,ivbin)
            schi=schi+((sumadp(ivel)-y1p(ivel))/sadfer(ivel,ivbin))**2
            ymax=max(ymax,sumadph(ivel),y1p(ivel))
            ymin=min(ymin,sumadpl(ivel),y1p(ivel))
         enddo

         schi2=schi2+schi
         prob=gammq(float(Nvel-1)/2.,schi/2.)

         chip(ivbin)=schi
         chimin=min(chimin,chip(ivbin))
         chimax=max(chimax,chip(ivbin))

c         call pgenv(velm(1)-vbit,velm(Nvel)+vbit,0.,1.05*ymax,0,0)
         ymin=min(0.,ymin)
         call pgenv(velm(1)-vbit,velm(Nvel)+vbit,ymin,1.05*ymax,0,0)
         call pglabel('V (km/s)','f(V)','')
         call pgsls(1)
         call pgsch(1.9)
         call pgsci(2)
         call pgpoint(Nvel,velm,y1p,17)
         call pgsci(1)
         if(rbin(ivbin)*angrad.lt.2..and.isee(ivbin).ne.1) then
c            call pgpoint(Nvel,velm,y1p2,17)
c            call pgline(Nvel,velm,y1p2)
         endif
         call pgpoint(Nvel,velm,sumadp,23)
         call pgerry(Nvel,velm,sumadph,sumadpl,1.5)
         call pgsch(1.5)
         do ip=1,Nvel
            write(12,2002) ivbin,ip,velm(ip),sumadp(ip),sumadph(ip),
     $           sumadpl(ip),y1p(ip)
         enddo
         if(iminor(ivbin).eq.1) then
            cadd='mn'
         else
            cadd=clabel(iang(ivbin))
         endif
         if(rbin(ivbin).eq.0) cadd=''
         irb=nint(rbin(ivbin)*angrad*100.)
         call pgnumb(irb,-2,0,cinc,nc)
         call pgmtxt('T',-1.5,0.05,0.,'R\\D'//cadd//' \\U = '
     $        //cinc(1:nc)//'\\(2728)')
         if(prob.le.0.2) then
            if(prob.le..005) then
               call pgmtxt('T',-2.7,0.05,0.,'P < 0.01')
            else
               irb=nint(prob*100.)
               call pgnumb(irb,-2,0,cinc,nc)
               call pgmtxt('T',-2.7,0.05,0.,'P = '//cinc(1:nc))
            endif
         endif
         call pgsch(1.2)
      enddo
      close(12)
 2002 format(1x,i3,1x,i3,1x,f7.1,4(1x,f10.8))

      call pgenv(rani(1)-rbit,rani(Nrani),0.,chimax+chimax/10.,0,10)
      do i=1,Nvelb
         call pgpoint(1,rbinl(i),chip(i),17)
      enddo
      call pglabel('R','\\gx\\U2\\D','')

      if(ifitherm.eq.1) then
         n1=0
         n2=0
         do i=1,Nvelb
            if(iminor(i).eq.1.and.rbin(i).ne.0) then
               n1=n1+1
               h3mmin(n1)=h3m(i)
               h3dmin(n1)=h3d(i)
               h4mmin(n1)=h4m(i)
               h4dmin(n1)=h4d(i)
               rlmin(n1)=rbinl(i)
            elseif(iang(i).eq.1.or.rbin(i).eq.0) then
               n2=n2+1
               h3mmaj(n2)=h3m(i)
               h3dmaj(n2)=h3d(i)
               h4mmaj(n2)=h4m(i)
               h4dmaj(n2)=h4d(i)
               rlmaj(n2)=rbinl(i)
            endif
         enddo
         call gminmax(Nvelb,h3m,h3d,h3min,h3max)
         call pgenv(rani(1)-rbit,rani(Nrani),h3min,h3max,0,10)
         call kgline(n1,rlmin,h3mmin,4)
         call kgline(n2,rlmaj,h3mmaj,1)
         call kgline(n1,rlmin,h3dmin,3)
         call kgline(n2,rlmaj,h3dmaj,2)
         call pglabel('R','H3','')

         call gminmax(Nvelb,h4m,h4d,h4min,h4max)
         call pgenv(rani(1)-rbit,rani(Nrani),h4min,h4max,0,10)
         call kgline(n1,rlmin,h4mmin,4)
         call kgline(n2,rlmaj,h4mmaj,1)
         call kgline(n1,rlmin,h4dmin,3)
         call kgline(n2,rlmaj,h4dmaj,2)
         call pglabel('R','H4','')
      endif

      open(unit=1,file='gherm.out',status='unknown')
      do i=1,Nvelb
         write(1,*) iang(i),rbin(i)*angrad,velmp(i),sigm(i),
     $        h3m(i),h4m(i)
      enddo
      close(1)

      prob=gammq(float(Nvel*Nvelb-1)/2.,schi2/2.)
      print *,'Total Chi Prob = ',prob

      call pgend

      end

      subroutine gminmax(n,x1,x2,xmin,xmax)
      real x1(n),x2(n)
      data big/1.e30/
      xmin=big
      xmax=-big
      do i=1,n
         xmin=min(xmin,x1(i),x2(i))
         xmax=max(xmax,x1(i),x2(i))
      enddo
      xbit=(xmax-xmin)/15.
      xmin=xmin-xbit
      xmax=xmax+xbit
      return
      end

      subroutine kgline(n,x,y,il)
      real x(n),y(n)
      call pgsls(il)
      call pgline(n,x,y)
      call pgsls(1)
      return
      end
