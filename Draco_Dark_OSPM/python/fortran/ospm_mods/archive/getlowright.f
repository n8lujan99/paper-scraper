C=============================================================================D
C     subroutine GETLOWRIGHT assembles the lower left portion of the 
C       matrix Cm which is used to assemble the matrix Am in SPEAR
C
C=============================================================================D
      SUBROUTINE getlowright(Cm)
      INCLUDE 'moddefs.h'
      parameter(lds=Norbitm+Nvel*Nvelbm,lda=Nbin+Nvel*Nvelbm)
      DOUBLE PRECISION Cm(lda,lds)

      data big /1.e20/

      Nvtot=Nvel*ircmax
      ilast = Norbit+Nvtot+1
      xmu=w(ilast)

      do irc=1,ircmax
         do ivel=1,Nvel
            i=Nbin+(irc-1)*Nvel+ivel

c            diff=big
c            do k=1,nad
c               if(abs(veld(k)/xmu-velm(ivel)).lt.diff) then
c                  diff=abs(veld(k)/xmu-velm(ivel))
c                  idiff=k
c               endif
c            enddo

c            Cm(i,ilast)=-ad(idiff,irc)
            Cm(i,ilast)=dble(-sumad(ivel,irc))

         enddo
      enddo

      RETURN
      END
