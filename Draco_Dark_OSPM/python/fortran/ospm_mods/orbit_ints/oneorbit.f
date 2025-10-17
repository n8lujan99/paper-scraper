C=============================================================================D
C     program ONEORBIT runs a single orbit.  It is a condensed version of
C       librarian() and orbit().
C
C     This program is not used by library nor by model.
C
C=============================================================================D
      PROGRAM oneorbit
      INCLUDE 'libdefs.h'
      DIMENSION pos0(Nvar),pos(Nvar),dposdt(Nvar)
      PARAMETER (epsilon=0.025,nrmax=10000)
      CHARACTER dummy
      real D3T(Nrdat,Nvdat),xlibT(Nrdat,Nvdat),cont(10),cont2(10)
      real v0liba(Nrani,Ntani),vr2liba(Nrani,Ntani),vt2liba(Nrani,Ntani)
      real vp1liba(Nrani,Ntani),vp2liba(Nrani,Ntani)
      real v1libT(Nvel,Nrdat,Nvdat)
      real dEm(Nrani,Ntani),tm(Nrani,Ntani),Estart
      real rm(Nrani,Ntani),dr(Nrani,Ntani),dv(Nrani,Ntani)
      real xalt,yalt,xneu,yneu,rakt,vlt(Nvel)
      real rsosmean,rsossig,rsosdistmean,rsosdistsig
      
      integer idot(Nrani,Ntani),Nsos2,Npdat2,ivisit(Nrani,Ntani)
      integer ivalt,iralt,ipt,irneu,ivneu,Nskip2,iskip2

      character string*4
      external trans,ekin,derivs

      data big/1.e20/
      iseed = -1
      iskip=0
      iskip2=0

      sum1=0.
      sum2=0.
      nsum=0

      iplot=0
      if(iplot.eq.1) then
         call pgbegin(0,'/cps',2,2)
         call pgscf(2)
         call pgslw(2)
         call pgsch(1.1)
      endif

      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call dataread()
      call galaxyread()
      call spacel()
      cden = dL(1,1)*ratML(1,1) / vol2d(1)
      core = 4.*pi*cden*rmin**3/(3.+cslope)
c      print *,rmin,cden,core
      call density()
      call haloread()
      fac1 = sqrt(GG*totlight*arcsec/distance/angrad)

c      write (6,100)
c100   format ('Do you have a tables.out file you want to use?  y/n  ',$)
c      read (5,*)dummy
      dummy='y'
      if (dummy.eq.'y'.or.dummy.eq.'yes'.or.dummy.eq.'Y'.or.
     &   dummy.eq.'YES') then
        call tablesread()
      else
        call tables()
        call tableswrite()
      endif

      if(iquad.ne.0) call forcefix()
      call forcewrite()


c --- test segment
      open (unit=82,file='potential.file',status='unknown')
      do irr=1,1000
         r=10**(log10(perimin)+(log10(apomax)-log10(perimin))/999.*
     &           float(irr-1))
         do ivv=1,50
            v=float(ivv-1)/999.
            call potential(r,v,VV0)
            write (82,*) r*angrad/206.265*100.,v,VV0*fac1*fac1
         end do
      enddo
      close(82)
c --- end


c8     write (6,*)''

9     format ('      Nskip = ',$)
      write (6,9)
      read (5,*)Nskip
      Nskip2 = 4*Nskip

c 676  write(*,"('# crossings, Npdat, epsilon fac : '$)")
c      read(*,*) Nsos2,Npdat2,efac
      efac=1.
 676  write(*,"('# crossings, Npdat : '$)")
      read(*,*) Nsos2,Npdat2
c      if(Nsos2.gt.Nsos) then
c         write(*,*) 'Too many crossings'
c         goto 676
c     endif

      write (6,*)''
10    format (' launch iLz = ',$)
c      write (6,10)
c11    write(*,"('Input iLz, iE, and iv : '$)")
c      read *,iLz,iE,iv
c 11   write(*,"('Input iLz, iE, and v : '$)")
c      read *,iLz,iE,xv
      istart=0
 11   continue
      if(istart.eq.0) then
          write(*,"('Input rperi [arcsec], rapo [arcsec], and v : '$)")
         read *,xperi,xapo,xv
      else
c         write(*,"('Input rperi [arcsec],',
c     &        ' rapo [arcsec], and rstart [arcsec] : '$)")
c         read *,xperi,xapo,xv
         write(*,"('Input rperi [sec],',
     &        ' rapo [sec], rstart [sec],',
     &        ' vr [km/s] :'$)")
         read *,xperi,xapo,xr,xvr
         r = xr/angrad
         vr = xvr/fac1
         xv=1.0
      end if
      if(xv.lt.0.0) then
         istart=1
         goto 11
      end if

      rperi=xperi/angrad
      rapo=xapo/angrad
      if(istart.eq.0) then
         v = xv
         iv = IVneeV(v)
      else
         iv=99
      end if
      iE=IRneeR(rapo)
      iLz=IRneeR(rperi)

      print*,'launching iE iLz iv:',iE,iLz,iv

      call potential(rperi,0.,xphiperi)
      call potential(rapo,0.,xphiapo)
      if(xphiapo.le.xphiperi) goto 11
      xLz=sqrt(2.*(xphiapo-xphiperi)/(rapo**2-rperi**2))*rapo*rperi
      Etot=xphiperi+0.5*xLz*xLz/rperi/rperi

      
      open(unit=15,file='bin_r.out',status='old')
      read(15,*)
      read(15,*)
      do ira=1,Nrani
         read(15,*) i1,x1,x2,x3,x4
         do ita=1,Ntani
            dr(ira,ita)=(x3-x1)/angrad
         end do
      end do
      close(15)

      open(unit=15,file='bin_v.out',status='old')
      read(15,*)
      read(15,*)
      do ita=1,Ntani
         read(15,*) i1,x1,x2,x3,x4
         do ira=1,Nrani
            dv(ira,ita)=2.*pi*(x3-x1)/360.*RneeIRC(ira)
         end do
      end do
      close(15)

      open(unit=15,file='drdv.out',status='unknown')
      do ira=1,Nrani
         do ita=1,Ntani
            write(15,*) ira,ita,dr(ira,ita)/dv(ira,ita)
         end do
      end do
      close(15)


      Ecirc = energy(rperi,0.,0.,0.)
      
      dE = Etot-Ecirc
      iorblo=iorb

      iv0 = iv
      if(istart.eq.0) then
         do ir=Nrdat+1,2,-1
            rlo=RneeIR(ir-1)
            rup=RneeIR(ir)
            if (Ekin(rup,v).lt.0.0.and.Ekin(rlo,v).gt.0.0) goto 1888
         enddo

         print *,'Try Again'
c         goto 11
         goto 666

 1888    r = rtbis(ekin,rlo,rup,(rup-rlo)*1.e-6,v)
         
         vmin=1.e10
         vmax=-1.e10
         ivminp=100
         ivmaxp=-100

         th  = asin(v)
         vr  = 0.
         vth = 0.
      else
         th  = 0.
         if(2.*Ekin(r,0.).lt.vr**2) then
            print*,'Try Again!'
c            goto 11
            goto 666
         end if
         vth = sqrt(2.*Ekin(r,0.)-vr*vr)
      end if

c      open(unit=21,file='l2.out',status='unknown')
c      vphi = xLz/r/cos(th)
c      do i=1,nrmax
c         rakt=rperi+(rapo-rperi)/float(nrmax-1)*float(i-1)
c         if(Ekin(rakt,0.).ge.0.0) then
c            vth = sqrt(2.*Ekin(rakt,0.))
c            xL2 = rakt**2*(vth**2+vphi**2)
c            write(21,*) rperi*angrad,rapo*angrad,rakt*angrad,
c     &           (rapo-rakt)*angrad,xL2
c         end if
c      end do
c      close(21)
c
c      open(unit=21,file='l2sample.out',status='unknown')
c      vphi = xLz/r/cos(th)
c      do i=1,Nvdat,ivrat
c         v = VneeIV(i)
c         do ir=Nrdat+1,2,-1
c            rlo=RneeIR(ir-1)
c            rup=RneeIR(ir)
c            if (Ekin(rup,v).lt.0.0.and.Ekin(rlo,v).gt.0.0) goto 1919
c         enddo
c         goto 1920
c 1919    rakt = rtbis(ekin,rlo,rup,(rup-rlo)*1.e-6,v)
c         xL2 = rakt**2*vphi**2
c         write(21,*) i,xL2
c 1920    continue
c      end do
c      xL2alt=xL2
c      do i=1,Nrdat,irrat
c         rakt=RneeIR(i)
c         if(rakt.ge.rperi.and.rakt.le.rapo) then
c            if(Ekin(rakt,0.).ge.0.0) then
c               vth = sqrt(2.*Ekin(rakt,0.))
c               xL2 = rakt**2*(vth**2+vphi**2)
c               if(xL2.gt.xL2alt) then
c                  write(21,*) i,xL2
c               end if
c               xL2alt=xL2
c            end if
c         end if
c      end do
c      close(21)

 2113 continue

      pos0(1) = r
      pos0(2) = th
      pos0(3) = vr
      pos0(4) = vth

      vphi = xLz/r/cos(th)
      xL2 = r**2*(vth**2+vphi**2)
      
      iralt = -1
      ivalt = -1

      Estart = energy(r,th,vr,vth)
      Etot = Estart

      xelast = Etot
      xl2last = xL2

C---- initialization

      do iv=1,Nvdat
      do ir=1,Nrdat
        xlibT(ir,iv)=0.
        D3T(ir,iv)=0.
        do ivel=1,Nvel
           v1libT(ivel,ir,iv)=0.
        enddo
      enddo
      enddo

      do ira=1,Nrani
         do ita=1,Ntani
            v0liba(ira,ita)=0.
            vr2liba(ira,ita)=0.
            vt2liba(ira,ita)=0.
            vp1liba(ira,ita)=0.
            vp2liba(ira,ita)=0.
            dEm(ira,ita)=0.
            rm(ira,ita)=0.
            ivisit(ira,ita)=0
            tm(ira,ita)=0.
            idot(ira,ita)=0
         enddo
      enddo

      t=0.d0
      isos=0
      iesc=0
      thold = pos0(2)
      r_a = -1.e12
      r_p =  1.e12

      open(unit=50,file='orbit.out',status='unknown')
      open(unit=51,file='sos.out',status='unknown')
      open(unit=52,file='info.out',status='unknown')
      open(unit=53,file='mer.out',status='unknown')
      open(unit=54,file='denergy.out',status='unknown')
      open(unit=55,file='derivs.out',status='unknown')
      open(unit=56,file='dpos.out',status='unknown')
      open(unit=57,file='proorbit.out',status='unknown')
      open(unit=58,file='losvdorbit.out',status='unknown')
961   format (5(e12.6,1x))
962   format (4(e12.6,2x))
964   format (f13.7,2x,f13.7)

      write (52,*)' '
c      write (52,*)'      Nskip =',Nskip
c      write (52,*)' '
      write (52,*)' launch iLz =',iLz
      write (52,*)'  launch iE =',iE
      write (52,*)'  launch iv =',iv0
      write (52,*)' '

C---- TIME LOOP:
999     iskip=iskip+1
        iskip2=iskip2+1
        call derivs(pos0,dposdt)
        call step(pos0,xLz,dt,Nvar,dposdt,epsilon*efac)

 1212   format (f7.3,1x,f7.4,2x,4(f12.6,1x))
        write(55,1212) log10(pos(1)*angrad),pos(2),
     &       dposdt(1),dposdt(2),dposdt(3),dposdt(4)

        r = pos0(1)
        r_a = max(r_a,r)
        r_p = min(r_p,r)
        th = pos0(2)
        v = sin(th)
        cv = cos(th)
        vr = pos0(3)
        vth = pos0(4)
        vphi = xLz/r/cv

        if(IVneeV(v).lt.1) then
           ivneu = (-IVneeV(v)-1)/ivrat+1
        else
           ivneu = (IVneeV(v)-1)/ivrat+1
        end if
        irneu = IRCneeR(r)
        xneu=r*cos(th)
        yneu=r*sin(th)
c        if(iralt.ne.irneu.or.ivalt.ne.ivneu) then
           ivisit(irneu,ivneu)=ivisit(irneu,ivneu)+1
c        end if
        rm(irneu,ivneu)=rm(irneu,ivneu)+
     &       sqrt((xneu-xalt)**2+(yneu-yalt)**2)
        iralt=irneu
        ivalt=ivneu
        xalt = r*cos(th)
        yalt = r*sin(th)

        t = t + dt
        xde = energy(pos0(1),pos0(2),pos0(3),pos0(4))/Etot-1.
        xde2 = (energy(pos0(1),pos0(2),pos0(3),pos0(4))-xelast)/Etot
        xdl = r**2*(vth**2+vphi**2)/xL2-1.
        xdl2 = (r**2*(vth**2+vphi**2)-xl2last)/xL2

        xelast=energy(pos0(1),pos0(2),pos0(3),pos0(4))
        xl2last=r**2*(vth**2+vphi**2)

c        write(54,*) t,r*angrad,v,xde,xdl

        if(iv0.gt.0) then
           if (thold.le.0.0.and.th.ge.0.0) then
c           if(thold*th.le.0.0) then
              write (51,962)r,vr*fac1,vth*fac1,vphi*fac1
              isos=isos+1
              if(isos.le.Nsos2) then
                 rsos(isos) = r
                 Vrsos(isos) = abs(vr)
              endif
              write(*,'("isos = ",i4,a,$)') isos,char(13)
           endif
        else
           if(nint(th/2/pi).gt.nint(thold/2/pi)) then
              write (51,962)r,vr*fac1,vth*fac1,vphi*fac1
              isos=isos+1
              if(isos.le.Nsos2) then
                 rsos(isos) = r
                 Vrsos(isos) = abs(vr)
              endif
           endif
        end if
        ir3d=IRneeR(r)
        iv3d=IVneeV(v)
        if (ir3d.ge.1.and.ir3d.le.irmax)
     &    D3T(ir3d,abs(iv3d)) = D3T(ir3d,abs(iv3d)) + dt

C get the anisotropy stuff

c        ira=nint(r*float(Nrani-1)+1.)
        ira=IRCneeR(r)
        if(ira.le.Nrani.and.ira.ge.0) then
c           ita=nint(abs(2.*th/pi)*float(Ntani-1)+1.)
c           ita=IVCneeV(v)
           if(iv3d.ge.1) then
              ita=(iv3d-1)/ivrat+1
           else
              ita=(-iv3d-1)/ivrat+1
           end if

           if(ita.le.Ntani.and.ita.ge.1) then
              idot(ira,ita) = idot(ira,ita)+1
              v0liba(ira,ita)=v0liba(ira,ita)+dt
              vr2liba(ira,ita)=vr2liba(ira,ita)+dt*vr*vr
              vt2liba(ira,ita)=vt2liba(ira,ita)+dt*vth*vth
              vp1liba(ira,ita)=vp1liba(ira,ita)+dt*vphi
              vp2liba(ira,ita)=vp2liba(ira,ita)+dt*vphi*vphi
              dEm(ira,ita)=dEm(ira,ita)+
     &             (energy(r,th,vr,vth)-Estart)/Estart*dt
              tm(ira,ita)=tm(ira,ita)+dt
           endif
        endif

 2354   format (f11.5,1x,f7.3,1x,f7.4,1x,4(e9.3,1x),1x,f7.5)

        if (iskip.eq.Nskip) then

           write(54,2354) t,
     &          r*angrad,v,xde,xdl,xde2,xdl2,dt
           write (50,961)t,r,th*180./pi,xelast,xl2last
           write (53,962)r*cos(th),r*sin(th)
           iskip=0
           ex=energy(r,th,vr,vth)
           nsum=nsum+1
           sum1=sum1+ex*ex
           sum2=sum2+ex
        endif

C Un-comment if not edge-on:
C        tmp1=(vr*cv-vth*v)*sini
C        tmp2=vphi*sini
C        tmp3=(vr*v+vth*cv)*cosi
C        zt=r*v*sini
C and then comment these:
        tmp1=vr*cv-vth*v
        tmp2=vphi
        z=r*v
C to here.
        yt=r*cv

C un-comment if not edge-on:
C        do ireflect=1,2
C        tmp3=-tmp3

        do ip = 1,Npdat2

          phi = ran1(iseed)*pi
          sp = sin(phi)
          cp = cos(phi)

C un-comment if not edge-on:
C          z = -yt*cp*cosi + zt
          y = yt*sp

          if (iskip2.eq.Nskip2) then
             write (57,962) z,y
             iskip2=0
          endif

          rsky = sqrt(y*y+z*z)
          vsky = z/rsky
          irsky = IRneeR(rsky)
          ivsky = IVneeV(abs(vsky))

          if (irsky.le.irmax.and.irsky.gt.0) then
             xlibT(irsky,ivsky) = xlibT(irsky,ivsky) + dt

C pick one:  use first if edge-on, second if not
            tVc = tmp1*cp-tmp2*sp
C            tVc = tmp1*cp-tmp2*sp+tmp3
            ivel=nint(tVc*vmult+vadd)

            if(ivel.ge.1.and.ivel.le.Nvel) then
               v1libT(ivel,irsky,ivsky)=v1libT(ivel,irsky,ivsky)+dt
            endif

            vmax=max(vmax,tVc*fac1)
            vmin=min(vmin,tVc*fac1)
            ivminp=min(ivminp,ivel)
            ivmaxp=max(ivmaxp,ivel)
          endif
        enddo
c        enddo

        thold = pos0(2)
        call rk4(Nvar,pos0,dposdt,dt,pos,derivs)

        write(56,1212) log10(pos(1)*angrad),pos(2),
     &       pos(1)-pos0(1),pos(2)-pos0(2),
     &       pos(3)-pos0(3),pos(4)-pos0(4)

        do ii=1,Nvar
          pos0(ii)=pos(ii)
        enddo

        if(r.gt.10) goto 333

c      if (isos.lt.Nsos2.or.vr.gt.0.) goto 999
        if(isos.lt.Nsos2) goto 999
C---- END TIME LOOP

c      print *
c      print *,'ivmin, ivmax ',ivminp,ivmaxp
c      print *,'vmin vmax ',vmin,vmax

C---- scale maps

      tt = t*float(Npdat2)
c      tt = t*2.*float(Npdat2)

      xmax=-big
      xmax2=-big

      do iv=1,ivmax
      do ir=1,irmax
        xlibT(ir,iv) = xlibT(ir,iv) / tt
        D3T(ir,iv) = D3T(ir,iv) / t
        do ivel=1,Nvel
           v1libT(ivel,ir,iv)=v1libT(ivel,ir,iv)/tt
           vlt(ivel) = v1libT(ivel,ir,iv)
        enddo
        do ivel=1,Nvel
           v1libT(ivel,ir,iv)=(v1libT(ivel,ir,iv)+vlt(Nvel-ivel+1))/2.
           write(58,*) ir,iv,ivel,v1libT(ivel,ir,iv)
        end do
        xmax=max(D3T(ir,iv),xmax)
        xmax2=max(xlibT(ir,iv),xmax2)
      enddo
      enddo
      xmin=log(xmax/1000.)
      xmax=log(xmax)
      xmin2=log(xmax2/1000.)
      xmax2=log(xmax2)

      do ira=1,Nrani
         do ita=1,Ntani
            if(v0liba(ira,ita).ne.0.0) then
               v0liba(ira,ita)=v0liba(ira,ita)/t
               vr2liba(ira,ita)=vr2liba(ira,ita)/t
               vt2liba(ira,ita)=vt2liba(ira,ita)/t
               vp1liba(ira,ita)=vp1liba(ira,ita)/t
               vp2liba(ira,ita)=vp2liba(ira,ita)/t
               if(tm(ira,ita).gt.0.0) then
                  dEm(ira,ita)=dEm(ira,ita)/tm(ira,ita)
               end if
            end if
         end do
      end do

c      print*
c      print*
c      print*
c      print*,'internal moments: sr,st,sp'
c      print*

      open(unit=14,file='oneorbit.moments.out',status='unknown')
c      write(*,'(" ir  iv   r[arcsec]   sr    st    sp    #iter",
c     &     "     v0liba")')
      write(14,'(" ir  iv   r[arcsec]   sr    st    sp    #iter",
     &     "     v0liba")')
      do ira=1,Nrani
         do ita=1,Ntani
         if(v0liba(ira,ita).ne.0.0) then
            xsr=sqrt(vr2liba(ira,ita)/v0liba(ira,ita))*fac1
            xst=sqrt(vt2liba(ira,ita)/v0liba(ira,ita))*fac1
            xsp=sqrt(vp2liba(ira,ita)/v0liba(ira,ita))*fac1
         else
            xsr=0.0
            xst=0.0
            xsp=0.0
         end if
         if(idot(ira,ita).ne.0) then
            write(14,'(2(i3,1x),3x,f6.2,3x,3(f6.1,1x),2x,i5,3x,f8.5,
     &           2x,f8.5)') 
     &           ira,ita,
     &           RneeIRC(ira)*angrad,xsr,xst,xsp,idot(ira,ita),
     &           v0liba(ira,ita),vr2liba(ira,ita)
c            write(*,'(2(i3,1x),3x,f6.2,3x,3(f6.1,1x),2x,i5,3x,f8.5,
c     &           2x,f8.5)') 
c     &           ira,ita,
c     &           RneeIRC(ira)*angrad,xsr,xst,xsp,idot(ira,ita),
c     &           v0liba(ira,ita),vr2liba(ira,ita)
         end if
      end do
      end do
      close(14)


c      print*
c      print*
c      print*
c      print*,'energy:'
c      print*

      open(unit=14,file='oneorbit.energy.out',status='unknown')
c      write(*,'(" ir  iv   r[arcsec]     dE")')
      write(14,'(" ir  iv   r[arcsec]     dE")')
      do ira=1,Nrani
         do ita=1,Ntani
            if(dEm(ira,ita).ne.0.0) then
               xsr=sqrt(vr2liba(ira,ita)/v0liba(ira,ita))*fac1
               write(14,'(2(i3,1x),3x,f6.2,3x,(f11.8))') ira,ita,
     &              RneeIRC(ira)*angrad,dEm(ira,ita)
c               write(*,'(2(i3,1x),3x,f6.2,3x,(f11.8))') 
c     &              ira,ita,RneeIRC(ira)*angrad,dEm(ira,ita)
            end if
         end do
      end do
      close(14)


c      print*
c      print*
c      print*
c      print*,'bins:'
c      print*

      open(unit=14,file='oneorbit.bins.out',status='unknown')
c      write(*,'(" ir  iv   r[arcsec]    rm/dr  rm/dv",    
c     $     "  dr/dv   #visit  #iter/#visit")')
      write(14,'(" ir  iv   r[arcsec]    rm/dr  rm/dv",    
     $     "  dr/dv   #visit  #iter/#visit")')
      do ira=1,Nrani
         do ita=1,Ntani
            if(ivisit(ira,ita).ne.0) then
               xrm = rm(ira,ita)/float(ivisit(ira,ita))
               ipt = nint(float(idot(ira,ita))/float(ivisit(ira,ita)))
               write(14,'(2(i3,1x),3x,f6.2,5x,2(f6.3,1x),1x,
     &              f5.2,2x,i5,6x,i5)')
     &              ira,ita,
     &              RneeIRC(ira)*angrad,xrm/dr(ira,ita),xrm/dv(ira,ita),
     &              dr(ira,ita)/dv(ira,ita),ivisit(ira,ita),ipt
c               write(*,'(2(i3,1x),3x,f6.2,5x,2(f6.3,1x),1x,
c     &              f5.2,2x,i5,6x,i5)')
c     &              ira,ita,
c     &              RneeIRC(ira)*angrad,xrm/dr(ira,ita),xrm/dv(ira,ita),
c     &              dr(ira,ita)/dv(ira,ita),ivisit(ira,ita),ipt
            end if
         end do
      end do
      close(14)

      ncont=10

      if(icheck.ne.1) then
         do i=1,ncont
            j=ncont-i+1
            cont(j)=exp(xmin+(xmax-xmin)*(i-1)/float(ncont-1))
            cont2(j)=exp(xmin2+(xmax2-xmin2)*(i-1)/float(ncont-1))
c            print *,Nsos2,j,cont(j),cont2(j)
            icheck=1
         enddo
      else
         do i=1,ncont
            x1=exp(xmin+(xmax-xmin)*(i-1)/float(ncont-1))
            x2=exp(xmin2+(xmax2-xmin2)*(i-1)/float(ncont-1))
c            print *,Nsos2,x1,x2
         enddo
      endif

      if(iplot.eq.1) then
      call pgnumb(Nsos2,0,0,string,nc)

      call pgenv(0.,r_a,0.,r_a,1,0)
      call pglabel('r','z','')
      call pgsch(1.5)
      call pgmtext('T',-1.5,.85,.5,string)
      call pgsch(1.1)
      call pgconx(D3T,Nrdat,Nvdat,1,Nrdat,1,Nvdat,cont,10,trans)

      call pgenv(0.,r_a,0.,r_a,1,0)
      call pglabel('Y','Z','')
      call pgsch(1.5)
      call pgmtext('T',-1.5,.85,.5,string)
      call pgsch(1.1)
      call pgconx(xlibT,Nrdat,Nvdat,1,Nrdat,1,Nvdat,cont2,10,trans)
      endif

      tt = tt/2.

      Efinal = energy(pos0(1),pos0(2),pos0(3),pos0(4))
      error = Efinal/Etot - 1.

      call sort2(Nsos2,rsos,Vrsos)
      call sosmom(Nsos2,xm,rsos,Vrsos,isos)
      write (52,*)' kurtosis of ',isos,'-point orbital s.o.s. = ',xm
      write (52,*)' '
      write (52,*)' peri=',r_p,'    apo=',r_a
      write (52,*)' '
      write (52,*)' dE/E = ',error
      write (52,*)' '
      rsosmean = 0.
      rsossig = 0.
      do i=1,Nsos2
         rsos(i)=log10(rsos(i))
         rsosmean = rsosmean + rsos(i)
      end do
      rsosmean = rsosmean/float(Nsos2)
      rsosdistmean = 0.
      do i=1,Nsos2
         rsosdistmean = rsosdistmean + 
     &        sqrt((rsos(i)-rsosmean)**2+vrsos(i)**2)
      end do
      rsosdistmean = rsosdistmean/float(Nsos2)

 333  continue
      close(50)
      close(51)
      close(52)
      close(53)
      close(54)
      close(55)
      close(56)
      close(57)
      close(58)

c      write (6,*)' '
c      write (6,*)' kurtosis of ',isos,'-point orbital s.o.s. = ',xm
c      write (6,*)' '
c      write (6,*)' peri=',r_p,'    apo=',r_a
c      write (6,*)' '
c      write (6,*)' dE/E = ',error
c      write (6,*)' '

      eavg=sum2/float(nsum)
      esig=(sum1+float(nsum)*eavg*eavg-2.*eavg*sum2)/float(nsum)
      esig=sqrt(abs(esig))
c      print*,nsum,sum2,eavg,sum1
      write(6,*) ''
      print *,'Average E, Sig, Sig/E = ',eavg,esig,abs(esig/eavg)

      write (6,*)''

      if(iplot.eq.1) call pgend

c      print*
c      print*
c      print*
c      print*,'bin statistics'
c      print*

      open(unit=27,file='oneorbit.test',status='unknown')
      do irc=1,Nrani
         do ivc=1,Ntani
            if(idot(irc,ivc).ne.0) then
               write(27,*) irc,ivc,RneeIRC(irc)*angrad,idot(irc,ivc)
c               print*,irc,ivc,RneeIRC(irc)*angrad,idot(irc,ivc)
            end if
         end do
      end do
      close(27)
 666  continue

      END

      subroutine trans(visble,x,y,z)
      integer visble
      r=RneeIR(nint(x))
      sv=VneeIV(nint(y))
      xworld=r*sqrt(1.-sv*sv)
      yworld=r*sv
      if(visble.eq.1) then
c         call pgdraw(xworld,yworld)
      else
c         call pgmove(xworld,yworld)
      endif
c      print *,x,y,xworld,yworld,z,visble
      end
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
