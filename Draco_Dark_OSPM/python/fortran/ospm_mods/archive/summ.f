C=============================================================================D
C     subroutine SUMM dots w() into xlib() for Nbin spatial bins and
C       leaves the result in the vector summed()
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE SUMM()
      INCLUDE 'moddefs.h'
      DOUBLE PRECISION sum

      do ibin=1,Nbin
        sum=0.d0
        do iorb=1,Norbit
          sum=sum+dble(w(iorb)*xlib(ibin,iorb))
        enddo
        summed(ibin)=sngl(sum)
      enddo

      do irc=1,Nvelb
         do ivel=1,Nvel
            i = Nbin+(irc-1)*Nvel+ivel
            sum=0.d0
            do iorb=1,Norbit
               sum=sum+dble(w(iorb)*v1lib(ivel,irc,iorb))
            enddo
            summed(i)=sngl(sum)
         enddo
      enddo

      RETURN
      END
