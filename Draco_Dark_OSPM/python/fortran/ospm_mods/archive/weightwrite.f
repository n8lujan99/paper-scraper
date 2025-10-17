C-----------------------------------------------------------------------------D
C     subroutine WEIGHTWRITE writes out the orbital weights of the maximum
C       entropy configuration
C
C     USED BY MODEL
C
C-----------------------------------------------------------------------------D
      SUBROUTINE weightwrite()
      INCLUDE 'moddefs.h'

      open (unit=55,file='weights.out',status='unknown')

      do i=1,Norbit+Nvel*Nvelb
        write (55,*)i,w(i)
      enddo

      close(55)
      RETURN
      END
