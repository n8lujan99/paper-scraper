C------------------------------------------------------------------------------D
C     subroutine WEIGHTREAD read the orbital weights from weights.out
C
C     USED BY MODEL
C
C------------------------------------------------------------------------------D
      SUBROUTINE weightread()
      INCLUDE 'moddefs.h'
      real wt(Norbitm+Nvel*Nvelbm)

      ilast=Norbit+Nvel*Nvelb

      open (unit=55,file='weights.out',status='old')

      do i=1,ilast
         read(55,*,end=666) ii,wt(ii)
      enddo
      do i=1,ilast
         w(i)=wt(i)
      enddo
      goto 676
 666  continue

      do i=1,Norb
         w(2*i-1)=wt(i)/2.
         w(2*i)=wt(i)/2.
      enddo
      do i=1,Nvel*Nvelb+1
         w(Norbit+i)=wt(Norb+i)
      enddo

 676  continue
      
      close(55)
      RETURN
      END
