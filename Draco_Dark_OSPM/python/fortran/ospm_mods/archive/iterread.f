C==============================================================================D
C     subroutine ITERREAD reads in iteration parameters
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE iterread()
      INCLUDE 'moddefs.h'
      CHARACTER*8 dummy

      open (unit=77,file='iter.params',status='old')

      read (77,*)dummy,Niter
      read (77,*)dummy,apfac
      read (77,*)dummy,apfacmu
      read (77,*)dummy,ifit
      read (77,*)dummy,alphainc
      read (77,*)dummy,alphastop
      read (77,*)dummy,fracrv

      close(77)

      RETURN
      END
