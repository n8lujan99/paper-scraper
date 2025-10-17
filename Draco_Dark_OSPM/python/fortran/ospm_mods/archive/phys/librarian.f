C=============================================================================D
C     subroutine LIBRARIAN samples the phase space available to a star in
C       the potential.  It picks launch coordinates and feeds them to
C       subroutine ORBIT
C
C     USED BY LIBRARY
C
C=============================================================================D
      SUBROUTINE librarian()
      INCLUDE 'libdefs.h'
      PARAMETER (epsilon=0.05)
      real xlibT(Nrdat,Nvdat),D3T(Nrdat,Nvdat),v1libT(Nvel,Nrdat,Nvdat)
      real xlib(Nrlib,Nvlib),v1lib(Nvel,Nrlib,Nvlib),pos0(Nvar)
      real v0liba(Nrani,Ntani),vr2liba(Nrani,Ntani),vt2liba(Nrani,Ntani)
      real vp1liba(Nrani,Ntani),vp2liba(Nrani,Ntani),sumv2(Nrlib)
      real vrtliba(Nrani,Ntani)
      real vrpliba(Nrani,Ntani),vrliba(Nrani,Ntani),vtliba(Nrani,Ntani)
      real vrtpliba(Nrani,Ntani),vtpliba(Nrani,Ntani)
      real slush(Nrdat*Nvdat),slush2(Nrdat*Nvdat*Nvel)
      real slush3(Nrdat*Nvdat),D3(Nrlib,Nvlib)
      real Eo,deltaE_lim,Etot2,J_2

      real rstart,vrstart,thstart,vthstart,ttot
      real rsmax,rsmin,Esos,vrsmax
      integer is,nlaunch(norbsmax),iuptemp,ilotemp,iLzarr(1000)
      real rs_min(norbsmax),rs_max(norbsmax)
      real vs_min(norbsmax),vs_max(norbsmax)
      real phi_peri,phi_apo

      integer norbs,n_vs,n_rs,i_rs,n_rs0
      real vt_up,r_min,v_max
      real v_th,tooclose,v_r,vlim,vstep
      real d_rsmax,d_rsmin,vstep_min
      real r_prev,th_prev,vr_prev,vth_prev

      real v_max_last,vth_last,v_min,v_losvd,r_lo,r_up
      real vlo,vfac,r_frac,r_frac_lastorb,rshell,velfac
      integer irun,nsec,isec
      real xseclo,xsecup,rup,rlo,r2

      real radsos(Nsos,norbsmax),velsos(Nsos,norbsmax),rd,vd
      integer nsosorb(norbsmax),i_thorb

      real deltaE_lim_max,deltaE_last

      real rdist,vdist,vda,vabs

      real rorb,Estart
      real themax,zmax
      real facnorb

      external ekin

      facnorb=1.0
      n_rs0 = Nl
 3333 continue

      idum=-1
      velfac=sqrt(GG*totlight*arcsec/distance/angrad)

      tooclose=0.01
      efac1=efac0

      ! energy conservation
      efac=efac1

      v_losvd = 1./vmult
      vstep_min = v_losvd*abs(vfrac)

      ! launch check
      rdist=0.05
      vdist=0.05
      vabs = vstep_min

799   format (e12.6,' ',e12.6)

 804  format ('iorb=',i6,' iLz=',i4,' iE=',i4,' ir=',i4,' iv=',i4,$)
 805  format (' rp,ra=',f9.4,',',f9.4,$)
 806  format (' dE= ',f7.5,$)

800   format ('iorb=',i6,' iLz=',i4,' iE=',i4,' ir=',i4,' iv=',i4,$)
801   format (' deltaE= ',f7.5,$)
802   format (' rp,ra=',f9.6,',',f9.6,$)
861   format(1x,f4.2,$)
862   format(3x,f4.2)

 950  format(2(f10.6,1x),5(f11.6,1x))

901   format (' ESCAPE',$)

      open(unit=70,file='xlib.out',status='unknown')
      open(unit=71,file='v1lib.out',status='unknown')
      open(unit=72,file='vani.out',status='unknown')
      open(unit=74,file='sumv2.out',status='unknown')
      open(unit=75,file='d3.out',status='unknown')
      open(unit=76,file='v1lib2.out',status='unknown')
      open(unit=86,file='highervani.out',status='unknown')
      open(unit=43,file='sos1.out',status='unknown')
      open(unit=44,file='sos2.out',status='unknown')
      open(unit=45,file='periapo.out',status='unknown')
      open(unit=46,file='rorb.out',status='unknown')
      open(unit=47,file='i3orb.out',status='unknown')
      open(unit=91,file='nstep.out',status='unknown')

      open(unit=61,file='orbinit.out',status='unknown')

      inc=1
      if(xinclin.gt.1.48) inc=0

C---- blank out phase-space weight array
      do iorb=1,Norbit
        wphase(iorb)=0.
      enddo
 
      iseed=-1
      iorb=1

C---- fake the first dLz,dE
      rperi=RneeIR(1)
      rapo=RneeIR(3)
      call potential(rperi,0.,phi_peri)
      call potential(rapo,0.,phi_apo)
      xLz1=sqrt(2.*(phi_apo-phi_peri)/(rapo**2-rperi**2))
     &     *rapo*rperi
      Etot1=phi_peri+0.5*xLz*xLz/rperi/rperi
      rperi=RneeIR(1+iLzstep)
      rapo=RneeIR(1+iLzstep+1)
      call potential(rperi,0.,phi_peri)
      call potential(rapo,0.,phi_apo)
      xLz2=sqrt(2.*(phi_apo-phi_peri)/(rapo**2-rperi**2))
     &     *rapo*rperi
      Etot2=phi_peri+0.5*xLz*xLz/rperi/rperi
      dLz = xLz2 - xLz1
      xLz0 = xLz1 - dLz
      dE = Etot2 - Etot1
      Eo = Etot1 - dE
c --- end


      if(dLz.lt.1.e15.or.dLz.gt.-1.e15) then
         ! ok
      else
         dLz = 0.001
      end if


C---- pick out launch z-component of orbital angular momentum
C---- angular momentum binning determined by the coarse bins

      iLzup = Nrelz
      iLzlo = 1 + ideltaLz

      open(unit=94,file='launchpos.orbits',status='unknown')
      open(unit=95,file='orbit.initials',status='unknown')

      nLarr=0
      do iLz=iLzlo,10,1
         nLarr=nLarr+1
         iLzarr(nLarr)=iLz
      enddo
      do iLz=12,iLzup,iLzstep
         nLarr=nLarr+1
         iLzarr(nLarr)=iLz
      enddo
c      do iLz=iLzlo,iLzup,iLzstep
      do iLarr=1,nLarr
         iLz=iLzarr(iLarr)

         rperi = RneeIRperi(iLz)

C------ pick E ranging from the circular orbit to the most elongated orbit
         iElo = iLz + ideltaE
         iEup = Nrelz

C------ energy binning determined by the coarse bins
        do iE=iElo,iEup,iEstep

           vrsmax = -1.e20

           if(iE.eq.iLz) then
              rapo = (RneeIRapo(IRneeRapo(rperi))+
     &             RneeIRapo(IRneeRapo(rperi)+1))/2.
           else
              rapo = RneeIRapo(iE)
           end if

           if(rapo.lt.rmin) goto 912 !skip sequence: outside spatial grid ...
           if(rperi.gt.1.) goto 912 !skip sequence: outside spatial grid ...

           call potential(rperi,0.,phi_peri)
           call potential(rapo,0.,phi_apo)

           if(phi_peri.ge.phi_apo) then
              call potential(rperi,0.,phi_apo)
              call potential(rapo,0.,phi_peri)
           end if
           
           deltaE_lim = deltaE_max

           xLz=sqrt(2.*abs(phi_apo-phi_peri)/
     &          (rapo**2-rperi**2))*rapo*rperi
           Etot=phi_peri+0.5*xLz*xLz/rperi/rperi

           dE = Etot-Eo
           iorblo=iorb

           norbs = 0

           r_min = rperi
           r_lo = rperi
           r_up = rapo

           i_thorb = 0

           if(iE.eq.iLz) then
              n_rs = 3
           else
              n_rs = n_rs0
           end if

           r_frac = .0

c --------------------
c --- NEW SEQUENCE ---
c --------------------

           v_max = -1.e10
           do i=1,1000
              r=log10(rperi)+(log10(rapo)-log10(rperi))/
     &             float(1000-1)*float(i-1)
              r=10**r
              if(Ekin(r,0.).ge.0.0) then
                 v_max = max(v_max,sqrt(2.*Ekin(r,0.)))
              end if
           end do

           v_max_last = v_max
           irun = 1
           i_rs = 0 

           if(rperi.gt.1.0) goto 1611

c ************************
c * drop orbits from ZVC *
c ************************

           iangbin = -10
           nsec = 100

 1399      continue
           iangbin = iangbin + ivstep
	   iangbin = max(iangbin,1)
           if(iangbin.gt.Nvdat) goto 1400
           v = VneeIV(iangbin)
           xseclo=log10(relzmin/2.)
           xsecup=log10(2.5*relzmax)
           do isec=nsec,2,-1
              rup=10**(xseclo+(xsecup-xseclo)/
     &             float(nsec-1)*float(isec-1))
              rlo=10**(xseclo+(xsecup-xseclo)/
     &             float(nsec-1)*float(isec-2))
              if (Ekin(rup,v).lt.0.0.and.Ekin(rlo,v).gt.0.0) then 
                 goto 1888
              end if
           enddo
           if(iangbin.lt.Nvdat) goto 1399
           goto 1400
                
 1888      r = rtbis(ekin,rlo,rup,(rup-rlo)*1.e-6,v)
           ir = -1

c           if(Ekin(r,th).le.0.0) goto 1399

           Estart = energy(r,th,0.,0.)
           if(abs(Estart/Etot-1.).gt.deltaE_lim) then
              goto 1399
           end if

           vr = 0.
           vth = 0.
           th = asin(v)

           write(61,'(3(i3,1x),2(f11.6,1x)," 0.1 4",f11.6)') 
     &          iLz,iE,ir,r,vr,vth

           goto 1700
           
 1400      continue

c          **************************
c          * sample orbits from SOS *
c          **************************
           i_rs = i_rs + 1

           if(irun.eq.1) then
              r=log10(r_lo)+(log10(r_up)-log10(r_lo))/
     &             float(n_rs-1)*float(i_rs-1)
              r=10**r
              
              if(r.eq.rperi) then
                 r2=log10(r_lo)+(log10(r_up)-log10(r_lo))/
     &             float(n_rs-1)*float(i_rs+1-1)
                 r2=10**r2
                 r=(r2+r)/2.
              end if
              
              if(Ekin(r,0.).ge.0.0) then
                 vt_up = sqrt(2.*Ekin(r,0.))
              else
c                 write(*,'("irs=",i3,": no orbits")') i_rs
                 if(i_rs.lt.n_rs) then
                    goto 1400
                 else
                    goto 1610
                 end if
              end if
              
              inew_r = 1
           end if

 8001      continue

           ir = i_rs

           if(irun.eq.1) then
              if(inew_r.eq.1) then
	         vstep = min(vstep_min,v_max_last*abs(vfrac))*facnorb
                 vr = 0.999*sqrt(2.*Ekin(r,0.))
                 vr = 0.001
                 if(norbs.gt.1) then
                    do i_orb=1,norbs
                       do is=1,nsosorb(i_orb)
                          rd = abs(radsos(is,i_orb)-r)/r
                          if(rd.lt.0.01) then
                             vr=max(vr,velsos(is,i_orb),0.)
                          end if
                       end do
                    end do
                 end if
                 if(vr.eq.0.001) vr = 0.999*sqrt(2.*Ekin(r,0.))
                 inew_r = 0
              else
                 vr = vr_last - vstep
                 if(vr.lt.0.0.and.iangbin.gt.Nvdat) then
                    if(abs(vr_last)/vstep.lt.0.2.and.
     &                   i_thorb.eq.0) then
                       vr=0.001*sqrt(2.*Ekin(r,0.))
                       i_thorb = 1
                    else
                       i_thorb=0
                       goto 1610
                    end if
                 end if
              end if              
           else
              r = rshell
              vth = sqrt(max(2.*Ekin(rshell,0.),0.))
              vr = sqrt(max(2.*Ekin(r,0.)-vth*vth,0.))
           endif
           
           vth = sqrt(max(2.*Ekin(r,0.)-vr*vr,0.))
           th=0.
           vth_last = vth
           vr_last = vr


           write(61,'(3(i3,1x),3(f11.6,1x)," 1",f11.6)') iLz,iE,ir,r,
     &          vr,vstep,vth

c          ********************************
c          * check if already orbit there *
c          ********************************
           if(irun.eq.1) then
              inew_orb = 1
              if(norbs.gt.0) then
                 do i_orb=1,norbs
                    do is=1,nsosorb(i_orb)
                       rd = abs(radsos(is,i_orb)-r)/r
                       if(vr.gt.0.0) then
                          vd = abs(velsos(is,i_orb)-vr)/vr
                          vda = abs(velsos(is,i_orb)-vr)
                       else
                          vd = 0.0
                       end if
                       if(rd.lt.rdist.and.vd.lt.vdist) inew_orb = 0
                       if(rd.lt.rdist.and.vda.lt.vabs) inew_orb = 0
                    end do
                 end do
              end if

              if(inew_orb.eq.1) then
                 goto 1700
              else
                 goto 1609
              end if
           end if

 1700      continue

           rstart = r
           vrstart = vr
           thstart = th
           vthstart = vth

c -----------------
c --- NEW ORBIT ---
c -----------------

           pos0(1) = r
           pos0(2) = th
           pos0(3) = vr
           pos0(4) = vth

           if(r.eq.r_prev.and.th.eq.th_prev.and.vr.eq.vr_prev.and.
     &          vth.eq.vth_prev) then
              goto 1611
           else
              r_prev = r
              th_prev = th
              vr_prev = vr
              vth_prev = vth
           end if

           efac = efac1

           Estart = energy(pos0(1),pos0(2),pos0(3),pos0(4))

 866       continue
           if(inc.eq.0) then
              call orbit1(pos0,epsilon*efac,iesc,isos,r_a,r_p,xlibT,
     &             v1libT,D3T,v0liba,vr2liba,vt2liba,vrtliba,
     $             vp1liba,vp2liba,vrliba,vtliba,vrtpliba,vrpliba,
     $             vtpliba,sumv2,Esos,rorb,themax,zmax,ttot)
           else
              call orbit2(pos0,epsilon*efac,iesc,isos,r_a,r_p,xlibT,
     &             v1libT,D3T,v0liba,vr2liba,vt2liba,vrtliba,
     $             vp1liba,vp2liba,vrliba,vtliba,vrtpliba,vrpliba,
     $             vtpliba,sumv2,Esos,rorb,themax,zmax,ttot)
           endif

           deltaE = Esos/Estart -1.

           write(61,'(3(i3,1x),3(f11.6,1x)," 2",f11.6)') iLz,iE,ir,
     &	        rstart,vrstart,vstep,vth
            
c -- check energy conservation: if > 5% then re-run with new time step,
c    but only do it for the most radial orbits

           if(abs(deltaE).gt.deltaE_lim.and.efac.gt.efac_min) then
              if(efac.eq.efac1) then
                 deltaE_last=deltaE
                 efac=efac/2.
                 write(6,'("Trying new efac : ",f7.4,a1,$)') efac,
     &                   char(13)
                 pos0(1) = r
                 pos0(2) = th
                 pos0(3) = vr
                 pos0(4) = vth
                 goto 866
              else
                 ! check if smaller stepsize improves enough
                 if(abs((deltaE_last-deltaE)/deltaE).gt.0.2) then
                    deltaE_last=deltaE
                    efac=efac/2.
                    write(6,'("Trying new efac : ",f7.4,a1,$)') efac,
     &                   char(13)
                    pos0(1) = r
                    pos0(2) = th
                    pos0(3) = vr
                    pos0(4) = vth
                    goto 866
                 else
                    ! skip this orbit
                    write(*,'("bad energy conservation",2(f6.2,2x),
     &                   a1,$)') deltaE,deltaE_last,char(13)
                    goto 1609
                 end if
              end if
           endif
           if(abs(deltaE).ge.0.0.and.abs(deltaE).le.deltaE_lim) then
              !OK
           else
              ! skip this orbit
              write(*,'("bad energy conservation",
     &                   a1,$)') char(13)
              goto 1609
           end if

           if(iesc.eq.1) then
              write(*,'("orbit escaped",
     &                   a1,$)') char(13)
              goto 1609
           end if
           if(iesc.eq.-1) then
              write(*,'("too many steps",
     &                   a1,$)') char(13)
              write(91,*) iLz,iE,iv,r,th,vr,vth
              goto 1609
           end if

c          ***********************************************
c          * check if orbit is too close to other orbits *
c          *********************************************** 
           rsmin=1.e20
           rsmax=-1.e20
           vrsmax=-1.e20
           do is=1,isos
              rsmin=min(rsmin,rsos(is))
              rsmax=max(rsmax,rsos(is))
              vrsmax=max(vrsmax,Vrsos(is))
           end do 
           if(norbs.ne.0) then
              do i_orb=1,norbs
                 d_rsmin=abs(rsmin-rs_min(i_orb))/rsmin*100.0
                 d_rsmax=abs(rsmax-rs_max(i_orb))/rsmax*100.0
                 if(d_rsmin.lt.tooclose.and.d_rsmax.lt.tooclose) then
                    goto 1609
                 end if
              end do
           end if
           
           norbs = norbs + 1
           if(norbs.ge.norbsmax) then
              ! stop sequence
              norbs = norbs - 1
              goto 1611
           end if
           iv = norbs

	   write(61,'(3(i3,1x),3(f11.6,1x)," 3",f11.6)') iLz,iE,ir,
     &	        rstart,vrstart,vstep,vth

           write(46,*) iorb,rorb
           write(47,*) iorb,zmax,themax
           write(94,*) iLz,iE,log10(rstart*angrad*distance/arcsec/1.e3),
     &          vrstart*velfac
           write(95,950) rstart,vrstart,thstart,vthstart,xLz,Etot,ttot
           
           write (80,800)iorb,iLz,iE,ir,iv
           cv = cos(th)
           vphi = xLz/r/cv
           J_2=r**2*(vth**2+vphi**2)*angrad*angrad*velfac**2
           Etot2=Etot*velfac**2
              
           write(79,*) xLz,Etot,Etot2,J_2
           write (50,800)iorb,iLz,iE,ir,iv
           write ( 6,804)iorb,iLz,iE,ir,iv

c          *********************************
c          * store orbital SOS information *
c          *********************************

           v_max_last = -1.e10
           do is=1,isos
              rs = rsos(is)
              v_r = Vrsos(is)
              v_max_last = max(v_max_last,v_r)
              radsos(is,norbs) = rs
              velsos(is,norbs) = v_r
           end do
           nsosorb(norbs) = isos

           if(dE.lt.1.e10.and.dE.gt.-1.e10) then
              ! dE is a number
           else
              dE=1.0
           end if

           do is=1,isos
              write(43,*) rsos(is),Vrsos(is),tsos(is)
              write(44,*) iorb,iE,iLz,dLz,dE
           end do
           
           write (80,801)deltaE
           write (50,801)deltaE
           write ( 6,806)deltaE
           write (80,802)r_p,r_a
           write (50,802)r_p,r_a
           write ( 6,805)r_p*angrad,r_a*angrad
           call area(iorb,isos,xcurv)
           call sosmom(Nsos,xm,rsos,Vrsos,isos)
           
           if (iesc.eq.1) then
              write (80,901)
              write (50,901)
              write ( 6,901)
           endif
           write (80,*)''
           write (50,*)''
           write ( 6,*)''

C---------- compress xlibT() by binning 4x4 and assign to library at iorb
           do ivS=1,Nvlib
              do irS=1,Nrlib
                 xlib(irS,ivS)=0.
                 D3(irS,ivS)=0.
                 do ivel=1,Nvel
                    v1lib(ivel,irS,ivS)=0.
                 enddo
                 do k=1,ivrat
                    do j=1,irrat
                       xlib(irS,ivS)=xlib(irS,ivS)+
     $                      xlibT((irS-1)*irrat+j,(ivS-1)*ivrat+k)
                       D3(irS,ivS)=D3(irS,ivS)+
     $                      D3T((irS-1)*irrat+j,(ivS-1)*ivrat+k)
                       do ivel=1,Nvel
                          v1lib(ivel,irS,ivS)=v1lib(ivel,irS,ivS)+
     $                         v1libT(ivel,(irS-1)*irrat+j,
     $                         (ivS-1)*ivrat+k)
                       enddo
                    enddo
                 enddo
              enddo
           enddo

C-write out the libraries
 
           it=0
           it2=0
           do irt=1,Nrlib
              do ivt=1,Nvlib
                 it=it+1
                 slush(it)=xlib(irt,ivt)
                 slush3(it)=D3(irt,ivt)
                 do ivel=1,Nvel
                    it2=it2+1
                    slush2(it2)=v1lib(ivel,irt,ivt)
                 enddo
              enddo
           enddo
           
           write(70,*) (slush(i),i=1,it)
           write(71,*) (slush2(i),i=1,it2)
           write(72,*) v0liba,vr2liba,vt2liba,vrtliba,vp1liba,vp2liba
           write(74,*) sumv2
           write(75,*) (slush3(i),i=1,it)
           
           write(86,*) vrliba,vtliba,vrtpliba,vrpliba,vtpliba
           
           sumt=0.
           do ivel=1,Nvel
              sum=0.
              do irt=1,Nrlib
                 do ivt=1,Nvlib
                    sum=sum+v1lib(ivel,irt,ivt)
                 enddo
              enddo
              sumt=sumt+sum
              write(76,861) sum
           enddo
           write(76,862) sumt
           
           iorb=iorb+1
           if(iorb.gt.Norbit) then
              iorbup=iorb-1
              call phasvol(iorblo,iorbup,dLz,dE)
              write (50,*)' '
              write ( 6,*)' '

              if(facnorb.lt.50.0) then
                 close(80)
                 close(70)
                 close(71)
                 close(72)
                 close(74)
                 close(75)
                 close(86)
                 close(43)
                 close(44)
                 close(45)
                 close(94)
                 close(95)
                 close(46)
                 close(47)
                 close(61)
                 open (unit=79,file='integrals.out',status='unknown')
                 open (unit=80,file='librarian.out',status='unknown')
                 facnorb = facnorb * 1.5
                 n_rs0 = max(n_rs0-10,20)
                 goto 3333      ! try again
              else
                 goto 2001
              end if
           endif
           
           write(45,*) iLz,iE,rperi,rapo,rsmin,rsmax
           
           rs_min(norbs)=rsmin
           rs_max(norbs)=rsmax
           vs_max(norbs)=vrsmax
           nlaunch(norbs)=iv
           
           do is=1,isos
              if(rsmin.eq.rsos(is)) then
                 vs_min(norbs)=Vrsos(is)/vrsmax
              end if
           end do

c          ***************************************
c          * check if we can finish the sequence *
c          ***************************************
           ilomax=-1
           do i_orb=1,norbs
              ilotemp=IRneeR(rs_min(i_orb))
              if(ilotemp.ge.ilomax) then
                 ilomax=ilotemp
                 iupmax=IRneeR(rs_max(i_orb))
                 ishell=i_orb
              end if
           end do
           r_frac = rs_min(ishell)/rs_max(ishell)
           if(r_frac.gt.rsfrac) then
              if(abs(vs_max(ishell)/v_max).lt.0.1.or.
     &             abs(vs_max(ishell)*velfac).lt.10.0) then
                 ! stop sequence
                 goto 1611
              else
c                 if(iangbin.ge.Nvdat) then
c                    irun=3
c                 end if
              end if
           end if
           if(irun.eq.2) then
              ! stop sequence
              goto 1611
           end if
           if(irun.eq.3) then
              ! do only the orbit with vth = max
              goto 8989
           end if
           
 1609      continue

           if(iangbin.le.Nvdat) goto 1399
           goto 8001 
           
 1610      continue
              
           ! end of i_rs loop
           if(iangbin.le.Nvdat) goto 1399
           if(i_rs.lt.n_rs) goto 1400

 8989      continue

c           goto 1611
           if(iE.eq.iLz) goto 1611

c          **************************************
c          * check if shell orbits are included *
c          **************************************
           ilomax=-1
           do i_orb=1,norbs
              ilotemp=IRneeR(rs_min(i_orb))
              if(ilotemp.ge.ilomax) then
                 ilomax=ilotemp
                 iupmax=IRneeR(rs_max(i_orb))
                 ishell=i_orb
              end if
           end do
           r_frac = rs_min(ishell)/rs_max(ishell)
           if(r_frac.lt.rsfrac.and.irun.ne.2) then
              ! add shell orbits
              rshell = (rs_min(ishell)+rs_max(ishell))/2.
              if(rshell/rs_min(ishell).gt.1.2) then
                 irun = 3 
                 rshell = (rs_min(ishell)+rshell)/2.
              else
                 irun = 2
              end if
              i_rs = 0
              goto 1400
           end if

 1611  continue

c ************************
c * sequence is finished *
c ************************

         iorbup=iorb-1
         call phasvol(iorblo,iorbup,dLz,dE)
         write (50,*)' '
         write ( 6,*)' '

         open (unit=73,file='phase.out',status='unknown')
         write (73,*)(1./wphase(i),i=1,iorb-1)
         close(73)

 912     continue

         Eo = Etot
C ------ end of energy loop

 2000 enddo

      xLz0 = xLz

C---- end of angular momentum loop

      enddo

 2001 Norbtot = iorb-1

C---- convert phase space volumes to densities
      do iorb=1,Norbtot
        wphase(iorb) = 1./wphase(iorb)
      enddo

      close(80)
      close(70)
      close(71)
      close(72)
      close(74)
      close(75)
      close(86)
      close(43)
      close(44)
      close(45)
      close(94)
      close(46)
      close(47)
      close(91)

      close(61)

      RETURN
      END



