C=============================================================================D
C     subroutine CMATRIX assembles the matrix Cm which is used to assemble
C       the matrix Am in spear()
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE cmatrix(Cm)
      INCLUDE 'moddefs.h'
      parameter(lds=Norbitm+Nvel*Nvelbm,lda=Nbin+Nvel*Nvelbm)
      DOUBLE PRECISION Cm(lda,lds)

      Nvtot=Nvel*Nvelb
      ilast = Norbit+Nvtot

C -- left two pieces
      do i=1,Norbit

C -- upper left piece
        do j=1,Nbin
          Cm(j,i) = dble(xlib(j,i))
        enddo

C -- lower left piece
        do irc=1,Nvelb
           do ivel=1,Nvel
              j=Nbin+(irc-1)*Nvel+ivel
              Cm(j,i)=dble(v1lib(ivel,irc,i))
           enddo
        enddo
      enddo

C -- right two pieces sans "1" entries in identity matrix (lower right)
      do i=Norbit+1,Norbit+Nvtot
        do j=1,Nbin+Nvel*Nvelb
          Cm(j,i) = 0.d0
        enddo
      enddo

C -- "1" entries in identity matrix except for last (lower right)
      do i=1,Nvtot
        Cm(Nbin+i,Norbit+i) = 1.d0
      enddo

      RETURN
      END
