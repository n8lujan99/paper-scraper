C==============================================================================D
C     subroutine PHASEREAD reads in phase-space weights from phase.out
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE phaseread()
      INCLUDE 'moddefs.h'
      real slush(Norbm)

      open(unit=68,file='norbit',status='old')
      read(68,*) Norb
      close(68)
      print *,'Norbit = ',Norb
      if(Norb.gt.Norbm) print *,'make Norbm bigger'
      if(Norbm.ne.Norbitm) then
         Norbit=2*Norb
      else
         Norbit=Norb
      endif

      open (unit=68,file='phase.out',status='old')
      read (68,*) (slush(i),i=1,Norb)
      close(68)

      do i=1,Norb
         if(Norbm.ne.Norbitm) then
            wphase(2*i-1)=slush(i)
            wphase(2*i)=slush(i)
         else
            wphase(i)=slush(i)
         endif
      enddo

      wmin=1.e20
      wmax=-1.e20
      do i=1,Norbit
         if(wphase(i).lt.wmin) then
            wmin=wphase(i)
            imin=i
         endif
         if(wphase(i).gt.wmax) then
            wmax=wphase(i)
            imax=i
         endif
      enddo
c      print *,'min and max phase vol : ',imin,wmin,imax,wmax

      RETURN
      END
