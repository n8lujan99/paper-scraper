C------------------------------------------------------------------------------D
C     subroutine PHASEWRITE writes out the orbital phase space densities
C
C     USED BY LIBRARY
C
C------------------------------------------------------------------------------D
      SUBROUTINE phasewrite()
      INCLUDE 'libdefs.h'

      open (unit=72,file='phase.out',status='unknown')
      write (72,*)(wphase(i),i=1,Norbtot)
      close(72)

      RETURN
      END
