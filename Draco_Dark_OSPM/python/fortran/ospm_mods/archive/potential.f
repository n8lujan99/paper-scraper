C==============================================================================D
C     subroutine POTENTIAL calculates the Legengre polynomial-expanded
C       potential from the tables Tabr1, Tabr2, and Tabth
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE potential(r,v,VV)
      INCLUDE 'libdefs.h'
      real VVhalo,VVhalolo,VVhaloup

      rir=(log10(r)-rlgpotmin)*float(npot-1)/(rlgpotmax-
     &     rlgpotmin)+1.
c      rir=(log10(r)-log10(rpotmin))*float(npot-1)/(log10(rpotmax)-
c     &     log10(rpotmin))+1.
      irlo=int(rir)
      irup=int(rir+1.)
      
      rlo = RneeIRhalo(irlo)
      rup = RneeIRhalo(irup)
      
      rdiff = r - rlo
      VVlo = 0.
      VVup = 0.

C------linearly interpolate for the radial sums, for Legendre orders 0,2,...,nlegup
      
      do n=0,nlegup,2
         nleg = n/2+1
         
         if (irlo.eq.0) then
            term = 0.
         else
            term = Pl(n,v)*tabv(nleg,irlo)
         endif
         
         VVlo = VVlo + term
         VVup = VVup + Pl(n,v)*tabv(nleg,irup)
      enddo
      
      VV = VVlo + (VVup-VVlo)/(rup-rlo) * rdiff
      VV = -2.*pi*VV - hole/totlight/gdennorm/r

      RETURN
      END

