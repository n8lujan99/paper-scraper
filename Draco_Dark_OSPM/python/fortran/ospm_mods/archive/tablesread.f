C==============================================================================D
C     subroutine TABLESREAD reads in the tables for use by oneorbit.f
C       and/or count.f (so that the lengthy table calculation won't have
C       to be re-performed for each singular orbit run)
C
C==============================================================================D
      SUBROUTINE tablesread()
      INCLUDE 'libdefs.h'
      open(unit=74,file='tables.out',status='old')
      do n=0,nlegup,2
         nleg = n/2+1
         do irtab = 1,npot
            read(74,*) i1,i2,x1
            tabv(i1,i2)=x1
            read(74,*) i1,i2,x1
            tabfr(i1,i2)=x1
            read(74,*) i1,i2,x1
            tabfv(i1,i2)=x1
         end do
      end do
      close(74)
c      open(unit=74,file='tables.out',status='old')
c      read(74,*) tabv,tabfr,tabfv
c      close(74)
      RETURN
      END
