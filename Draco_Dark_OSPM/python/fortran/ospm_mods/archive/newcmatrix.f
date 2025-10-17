C==============================================================================D
C     subroutine CMATRIX assembles the matrix Cm which is used to assemble
C       the matrix Am in spear()
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE newcmatrix(inew)
      INCLUDE 'moddefs.h'
      real xlibtemp(Nbin),v1libtemp(Nvel,Nvelbm)

      Nvtot=Nvel*Nvelb
      ilast = Norbit+Nvtot

      do i=1,Nbin
         xlibtemp(i)=xlib(i,inew)
      end do
      
      do i=1,Nvelb
         do ivel=1,Nvel
            v1libtemp(ivel,i) = v1lib(ivel,i,inew)
         end do
      end do

      do i=1,Nbin
         xlib(i,inew)=xlib(i,Norbit)
         xlib(i,Norbit)=xlibtemp(i)
      end do
      do i=1,Nvelb
         do ivel=1,Nvel
            v1lib(ivel,i,inew) = v1lib(ivel,i,Norbit)
            v1lib(ivel,i,Norbit) = v1libtemp(ivel,i)
         end do
      end do
         

C -- upper left piece
      do j=1,Nbin
         Cm(j,inew) = xlib(j,inew)
      enddo
      
C -- lower left piece
      do irc=1,Nvelb
         do ivel=1,Nvel
            j=Nbin+(irc-1)*Nvel+ivel
            Cm(j,inew)=v1lib(ivel,irc,inew)
         enddo
      enddo

      Norbit = Norbit-1
C -- right two pieces sans "1" entries in identity matrix (lower right)
      do i=Norbit+1,Norbit+Nvtot
        do j=1,Nbin+Nvel*Nvelb
          Cm(j,i) = 0.
        enddo
      enddo

C -- "1" entries in identity matrix except for last (lower right)
      do i=1,Nvtot
        Cm(Nbin+i,Norbit+i) = 1.
      enddo

      RETURN
      END
