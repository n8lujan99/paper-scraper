C=============================================================================D
C     subroutine ORBIT integrates an orbit given a set of launch parameters
C       from subroutine LIBRARY and bins the positions over Nsos periods
C       ORBIT1 is for edge-on systems (theta > 85)
C
C     USED BY LIBRARY
C
C=============================================================================D
      SUBROUTINE orbit1(pos0,epsilon,iesc,isos,r_a,r_p,xlibT,
     $     v1libT,D3T,v0liba,vr2liba,vt2liba,vrtliba,vp1liba,
     $     vp2liba,vrliba,vtliba,vrtpliba,vrpliba,vtpliba,sumv2,Esos,
     $     rorb,themax,zmax,t)
      INCLUDE 'libdefs.h'
      real pos0(Nvar),pos(Nvar),dposdt(Nvar)
      real xlibT(Nrdat,Nvdat),v1libT(Nvel,Nrdat,Nvdat),D3T(Nrdat,Nvdat)
      real v0liba(Nrani,Ntani),vr2liba(Nrani,Ntani),vt2liba(Nrani,Ntani)
      real vp1liba(Nrani,Ntani),vp2liba(Nrani,Ntani),sumv2(Nrlib)
      real vrtliba(Nrani,Ntani)
      real vrpliba(Nrani,Ntani),vrliba(Nrani,Ntani),vtliba(Nrani,Ntani)
      real vrtpliba(Nrani,Ntani),vtpliba(Nrani,Ntani)

      real v0_last(Nrani,Ntani),maxdiff,rorb
      real zmax,themax

      real Esos,Eakt,ts
      external derivs

      real xtemp(Nsos),ytemp(Nsos)
      integer Nsosnew,istart

C---- initialization

      Nsosnew = Nsosmin

      istart=1
      ts=0.
      Esos=0.
      rorb=0.

      themax=0.
      zmax=0.
      
      fac1 = sqrt(GG*totlight*arcsec/distance/angrad)
      do i=1,Nrlib
         sumv2(i)=0.
      enddo
      nstep=0
      t=0.
      isos=0
      iesc=0
      thold = pos0(2)
      r_a=-1.e12
      r_p=1.e12

      do iv=1,Nvdat
         do ir=1,Nrdat
            xlibT(ir,iv)=0.
            D3T(ir,iv)=0.
            do ivel=1,Nvel
               v1libT(ivel,ir,iv)=0.
            enddo
         end do
      end do

      do ira=1,Nrani
         do ita=1,Ntani
            v0liba(ira,ita)=0.
            vr2liba(ira,ita)=0.
            vt2liba(ira,ita)=0.
            vrtliba(ira,ita)=0.
            vp1liba(ira,ita)=0.
            vp2liba(ira,ita)=0.
            vrliba(ira,ita)=0.
            vtliba(ira,ita)=0.
            vrtpliba(ira,ita)=0.
            vrpliba(ira,ita)=0.
            vtpliba(ira,ita)=0.
	    v0_last(ira,ita)=0.
         enddo
      enddo

      call derivs(pos0,dposdt)

C---- TIME LOOP:

999     continue
        nstep=nstep+1

        call derivs(pos0,dposdt)
        call step(pos0,xLz,dt,Nvar,dposdt,epsilon)

        if(abs(energy(pos0(1),pos0(2),pos0(3),pos0(4))/Etot-1.).gt.
     &       deltaE_max) then
           Esos=666.
           goto 333
        end if

        r = pos0(1)
        r_a = max(r_a,r)
        r_p = min(r_p,r)
        th = pos0(2)
        if(abs(th).ge.pi/2.) then
           Esos=666.
           goto 333
        endif
        v = sin(th)
        cv = cos(th)
        vr = pos0(3)
        vth = pos0(4)
        vphi = xLz/r/cv

        zmax = max(zmax,r*v)
        themax = max(themax,th)

        t = t + dt
        ts = ts + dt

c --- SOS stuff
        if (thold.le.0.0.and.th.ge.0.0) then

           ! check if we can finish the integration
           if(isos.ge.Nsosnew) then
              maxdiff=-1.e10
              do irc=1,Nrlib
                 do ivc=1,Nvlib
                    if(v0_last(irc,ivc).ne.0.0) then
                       maxdiff = max(maxdiff,abs(v0liba(irc,ivc)/t-
     &                      v0_last(irc,ivc))**2/v0_last(irc,ivc)**2)
                    end if
                    v0_last(irc,ivc)=v0liba(irc,ivc)/t
                 end do
              end do

              if(maxdiff.gt.0.01.or.maxdiff.eq.-1.e10) then
                 Nsosnew = Nsosnew + 30
              else
                 goto 333
              end if

           end if

           if(istart.ne.1) then
              isos=isos+1
              rsos(isos) = r
              Vrsos(isos) = abs(vr)
              tsos(isos)=ts
              ts=0.
              Eakt = energy(pos0(1),pos0(2),pos0(3),pos0(4))
              if(abs(Eakt/Etot-1.).gt.abs(Esos/Etot-1.)) then
                 Esos=Eakt
              end if
c              if(abs(Eakt/Etot-1.).gt.0.075) then
c                 goto 333
c              end if
           else
              istart=0
              ts=0.
              Esos=Etot
           end if
        endif
c --- end

        ir3d=IRneeR(r)
        iv3d=IVneeV(v)
        isign=1
        if(iv3d.lt.0) then
           iv3d=-iv3d
           isign=-1
        endif
        if (ir3d.ge.1.and.ir3d.le.irmax)
     &    D3T(ir3d,iv3d) = D3T(ir3d,iv3d) + dt

C get the internal moments

        ira=IRCneeR(r)
        if(ira.le.Nrani.and.ira.ge.1) then
           ita=(iv3d-1)/ivrat+1
           if(ita.le.Ntani.and.ita.ge.1) then
              v0liba(ira,ita)=v0liba(ira,ita)+dt
              vr2liba(ira,ita)=vr2liba(ira,ita)+dt*vr*vr
              vt2liba(ira,ita)=vt2liba(ira,ita)+dt*vth*vth
              vrtliba(ira,ita)=vrtliba(ira,ita)+dt*vr*vth*isign
              vp1liba(ira,ita)=vp1liba(ira,ita)+dt*vphi
              vp2liba(ira,ita)=vp2liba(ira,ita)+dt*vphi*vphi
              vrliba(ira,ita)=vrliba(ira,ita)+dt*vr
              vtliba(ira,ita)=vtliba(ira,ita)+dt*vth*isign
              vrtpliba(ira,ita)=vrtpliba(ira,ita)+dt*vr*vphi*vth*isign
              vrpliba(ira,ita)=vrpliba(ira,ita)+dt*vr*vphi
              vtpliba(ira,ita)=vtpliba(ira,ita)+dt*vphi*vth*isign
           endif
        endif

        rorb = rorb + r*dt

c        time1=time1+dtime(tarr)

        tmp1=vr*cv-vth*v
        tmp2=vphi
        z=r*v
        yt=r*cv

        do ip = 1,Npdat
           phi = ran1(iseed)*pi
           sp = sin(phi)
           cp = cos(phi)
           y = yt*sp
           rsky = sqrt(y*y+z*z)
           vsky = z/rsky
           irsky = IRneeR(rsky)
           ivsky = IVneeV(abs(vsky))
 112       if (irsky.le.irmax.and.irsky.gt.0.and.
     &          ivsky.le.ivmax.and.ivsky.gt.0) then
              xlibT(irsky,ivsky) = xlibT(irsky,ivsky)+dt
              tVc = tmp1*cp-tmp2*sp
              ivel=nint(tVc*vmult+vadd)
              if(ivel.ge.1.and.ivel.le.Nvel) then
                 v1libT(ivel,irsky,ivsky)=v1libT(ivel,irsky,ivsky)+dt
              endif
           endif
        enddo

        thold = pos0(2)

        call rk4(Nvar,pos0,dposdt,dt,pos,derivs)
        
        do ii=1,Nvar
          pos0(ii)=pos(ii)
        enddo

        if(t.lt.1.e10) goto 212
        iesc=1
        goto 333
 212    continue

        if(r.gt.2.2) then
           iesc=1
           goto 333
        endif
        if(nstep.gt.5.e5) then
           iesc=-1
           goto 333
        endif

       if (isos.lt.Nsos) goto 999 
c        if (isos.lt.Nsos.or.vr.gt.0.0) goto 999
C---- END TIME LOOP

 333  continue

C---- scale maps

      tt = t*float(Npdat)
      do ir=1,Nrlib
         sumv2(ir)=sumv2(ir)/tt*fac1*fac1
      enddo

      do iv=1,ivmax
         do ir=1,irmax
            xlibT(ir,iv) = xlibT(ir,iv) / tt
            D3T(ir,iv) = D3T(ir,iv) / t
            do ivel=1,Nvel
               v1libT(ivel,ir,iv)=v1libT(ivel,ir,iv)/tt
            enddo
         end do
      end do
      
      do ira=1,Nrani
         do ita=1,Ntani
            v0liba(ira,ita)=v0liba(ira,ita)/t
            vr2liba(ira,ita)=vr2liba(ira,ita)/t
            vt2liba(ira,ita)=vt2liba(ira,ita)/t
            vrtliba(ira,ita)=vrtliba(ira,ita)/t
            vp1liba(ira,ita)=vp1liba(ira,ita)/t
            vp2liba(ira,ita)=vp2liba(ira,ita)/t
            vrliba(ira,ita)=vrliba(ira,ita)/t
            vtliba(ira,ita)=vtliba(ira,ita)/t
            vrtpliba(ira,ita)=vrtpliba(ira,ita)/t
            vrpliba(ira,ita)=vrpliba(ira,ita)/t
            vtpliba(ira,ita)=vtpliba(ira,ita)/t
         enddo
      enddo

      rorb = rorb / t

666   RETURN
      END
