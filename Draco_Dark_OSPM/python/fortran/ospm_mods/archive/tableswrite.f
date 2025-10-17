C==============================================================================D
C     subroutine TABLESWRITE writes out the tables for use by oneorbit.f
C       and/or count.f (so that the lengthy table calculation won't have
C       to be reperformed for each singular orbit run)
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE tableswrite()
      INCLUDE 'libdefs.h'

      open(unit=74,file='tables.out',status='unknown')
      do n=0,nlegup,2
         nleg = n/2+1
         do irtab = 1,npot
            write(74,*) nleg,irtab,tabv(nleg,irtab)
            write(74,*) nleg,irtab,tabfr(nleg,irtab)
            write(74,*) nleg,irtab,tabfv(nleg,irtab)
         end do
      end do
      close(74)

c      open(unit=74,file='tables.out',status='unknown')
c      write(74,*) tabv,tabfr,tabfv
c      close(74)
      RETURN
      END
