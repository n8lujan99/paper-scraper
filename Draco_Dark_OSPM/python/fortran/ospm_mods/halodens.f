      subroutine halodens(r,theta,xmass)
      include 'libdefs.h'
      real r,theta,xmass,rakt,w,xrho,xR,xZ,makt,mscale,xd
      real xh,xrhocrit

      if(ihalo.eq.3) then !non-singular isothermal spheroid (halo)
         rakt=r
         w=theta
         xrho=0.78722918/dis/dis
         xrho=xrho*v0*v0/qdm/qdm
         xR=rakt*cos(w)
         xZ=rakt*sin(w)
         xrho=xrho*((2.*qdm*qdm+1.)*rc*rc+xR*xR+2.*
     &        (1.-0.5/qdm/qdm)*xZ*xZ)
         xrho=xrho/(rc*rc+xR*xR+xZ*xZ/qdm/qdm)
         xrho=xrho/(rc*rc+xR*xR+xZ*xZ/qdm/qdm) ![xrho]=Msun/pc**3
      end if

      if(ihalo.eq.1) then !gamma-profile
         if(gamma.ne.0.0) then !Dehnen-Profil
            rakt=r
            w=theta
            mscale=rsgamma
            makt=r*sqrt(cos(w)*cos(w)+sin(w)*sin(w)/qdm/qdm)
            mscale=mscale*dis/206.265*1.e3 ![mscale]=pc
            makt=makt*dis/206.265*1.e3 ![makt]=pc
            xrho=xmgamma/4./pi*(3.-gamma)*mscale
            xrho=xrho/makt**(gamma)
            xrho=xrho/(mscale+makt)**(4.-gamma) ![xrho]=Msun/pc**3
         else ! Plummer-Profil
            rakt=r
            w=theta
            mscale=rsgamma
            makt=r*sqrt(cos(w)*cos(w)+sin(w)*sin(w)/qdm/qdm)
            mscale=mscale*dis/206.265*1.e3 ![mscale]=pc
            makt=makt*dis/206.265*1.e3 ![makt]=pc
            xrho=3.*xmgamma/4./pi/mscale/mscale/mscale
            xrho=xrho*(1.+(makt/mscale)*(makt/mscale))**(-5./2.) ![xrho]=Msun/pc**3
         end if
      end if

      if(ihalo.eq.2) then !nfw-profile
         xh=70.                 !hubble
         xhparam=xh/100.
         xrhocrit=2.7754996776e-7*xhparam**2 ![xrhocrit]=Msun/pc**3
         rakt=r
         w=theta
         mscale=rsnfw
         makt=r*sqrt(cos(w)*cos(w)+sin(w)*sin(w)/qdm/qdm)
         mscale=mscale*dis/206.265*1.e3 ![mscale]=pc
         makt=makt*dis/206.265*1.e3 ![makt]=pc
         xd=200./3.*cnfw*cnfw*cnfw/(log(1.+cnfw)-cnfw/(1.+cnfw))
         xrho=xrhocrit*xd/(makt/mscale)
         xrho=xrho/(1.+makt/mscale)/(1.+makt/mscale) ![xrho]=Msun/pc**3
      end if

      if(ihalo.eq.4) then !no DM-halo
         xrho=0.
      end if

      xmass=xrho
      
      return
      end


