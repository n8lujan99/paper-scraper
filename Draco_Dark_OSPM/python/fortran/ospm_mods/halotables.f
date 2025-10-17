      SUBROUTINE halotables()
      INCLUDE 'libdefs.h'
      parameter(limit=1000,lenw=limit*4)
      DIMENSION htab1(nlegdim,npot),htab2(nlegdim,npot)
      real work(lenw)
      integer iwork(limit)
      common /chfunc/ n
      external hfunc1,hfunc2

      data epsabs,epsrel /1.e-10,1.e-4/

C---- first calculate temporary tables:

      do n=0,nlegup,2
         nleg = n/2+1
         do irtab = 1,npot
            call qags(hfunc1,0.,RneeIRhalo(irtab),epsabs,epsrel,ss,
     $           abserr,neval,ier,limit,lenw,last,iwork,work)
c            if(ier.ne.0) print *,'                          first ',ier
            htab1(nleg,irtab) = ss*2./float(Nvdat)
            call qags(hfunc2,RneeIRhalo(irtab),5.1*RneeIRhalo(npot),
     $           epsabs,epsrel,ss,abserr,neval,ier,limit,
     $           lenw,last,iwork,work)
c            if(ier.ne.0) print *,'                          second ',ier
            htab2(nleg,irtab) = ss*2./float(Nvdat)
c            write(*,'("n, irtab = ",i4,i4,1x,f11.4,a1,$)') 
c     $           n,irtab,ss*2.,char(13)
c            call flush(6)
         enddo
      enddo
      

C---- now calculate real look-up tables

      do n=0,nlegup,2
      nleg = n/2+1

        do ir=1,npot

          r = RneeIRhalo(ir)

          htabv(nleg,ir) = r**(-n-1)*htab1(nleg,ir) + 
     &         r**n*htab2(nleg,ir)

          htabfr(nleg,ir) = float(n+1)*r**(-n-2)*htab1(nleg,ir) -
     &      float(n)*r**(n-1)*htab2(nleg,ir)

          htabfv(nleg,ir) = r**(-n-2)*htab1(nleg,ir) +
     &      r**(n-1)*htab2(nleg,ir)

        enddo

      enddo

      RETURN
      END

      function hfunc1(r)
      INCLUDE 'libdefs.h'
      common /chfunc/ n
      t1=0.
      do iv=1,Nvdat
         t1=t1+dalo(r,iv)*Pl(n,VneeIV(iv))
      enddo
      hfunc1=r**(n+2)*t1
      return
      end

      function hfunc2(r)
      INCLUDE 'libdefs.h'
      common /chfunc/ n
      t1=0.
      do iv=1,Nvdat
         t1=t1+dalo(r,iv)*Pl(n,VneeIV(iv))
      enddo
      hfunc2=r**(1-n)*t1
      return
      end

c      FUNCTION dalo(r,iv)
c      INCLUDE 'libdefs.h'
c      rir = log(a*r/b+1.)/a
c      irlo = int(rir)
c      irup = int(rir+1)
c 
c      rlo = log10(RneeIR(irlo))
c      rup = log10(RneeIR(irup))
c      dr = log10(r) - rlo
c      if(irup.gt.irmax+npot) then
c         dalo=hrho(irlo,iv)/(rlo-rup)*dr+hrho(irlo,iv)
c      else
c         dalo = 10**(
c     $        (log10(hrho(irlo,iv))-log10(hrho(irup,iv)))/(rlo-rup)*dr+
c     $        log10(hrho(irlo,iv)))
c      end if
c      RETURN
c      END

      FUNCTION dalo(r,iv)
      INCLUDE 'libdefs.h'
      xfactor = totlight/distance/distance/distance*arcsec*arcsec*
     &     arcsec/angrad/angrad/angrad 
      xr=r*angrad
      xtheta=asin(VneeIV(iv))
      call halodens(xr,xtheta,xhalodens)
      xhalodens=xhalodens/xfactor/gdennorm
      dalo=xhalodens
      RETURN
      END









