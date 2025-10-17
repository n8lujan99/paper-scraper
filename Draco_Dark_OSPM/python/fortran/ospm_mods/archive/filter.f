C=============================================================================D
C     subroutine FILTER weeds out negative weights and sets them to TINY
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE filter()
      INCLUDE 'moddefs.h'
      PARAMETER (TINY=1.e-37)
      do iorb=1,Norbit
         w(iorb)=max(TINY,w(iorb))
      enddo
      RETURN
      END
