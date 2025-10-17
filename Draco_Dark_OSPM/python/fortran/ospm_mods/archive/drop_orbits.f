C=============================================================================D
C     subroutine ENTROPY really calculates the profit function (i.e. the
C       function to be maximized),  which is the "entropy"  minus a chi-
C       square term due to slit errors minus a chi-square term due to the
C       error in the target value of (GM/Ro)^-1.  It also calculates the
C       first and second derivatives of the profit function.
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE drop_orbits()
      INCLUDE 'moddefs.h'
      PARAMETER (TINY=1.e-30)
      COMMON/drop/iorbiter(Norbitm),Norbito

      iorb = 0

 101  continue
      iorb = iorb + 1
      if(w(iorb).le.TINY) then

 102     continue
         if(w(Norbit).le.TINY) then
            Norbit = Norbit - 1
            goto 102
         end if

         w(iorb) = w(Norbit)
         w(Norbit) = 0.

         xtemp = wphase(iorb)
         wphase(iorb) = wphase(Norbit)
         wphase(Norbit) = xtemp
         iorbiter(iorb) = Norbit
         iorbiter(Norbit) = iorb
         call newcmatrix(iorb)
         print *,'dropped orbit ',iorb,Norbit,Norbito

      end if

      if(iorb.lt.Norbit) goto 101

      RETURN
      END




