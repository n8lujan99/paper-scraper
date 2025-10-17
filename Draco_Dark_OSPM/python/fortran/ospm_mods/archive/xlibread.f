C==============================================================================D
C     subroutine XLIBREAD reads in the orbit library from xlib.out and
C       smears out lowest-r bin information over all angles
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE xlibread()
      INCLUDE 'moddefs.h'
      DIMENSION slush(Nvlib*Nrlib)

      open (unit=71,file='xlib.out',status='old')

C---- smear out first coarse bin in theta (i.e. sum over iv)

      do iorb=1,Norb
         if(Norb.ne.Norbit) then
            smlib(1,2*iorb-1)=0.
            smlib(1,2*iorb)=0.
         else
            smlib(1,iorb)=0.
         endif
         read(71,*) slush
         do iv=1,Nvlib
            if(Norb.ne.Norbit) then
              smlib(1,2*iorb-1)=smlib(1,2*iorb-1)+slush(iv)
              smlib(1,2*iorb)=smlib(1,2*iorb)+slush(iv)
           else
              smlib(1,iorb)=smlib(1,iorb)+slush(iv)
           endif
        enddo
        do ibin=2,Nbin
           if(Norb.ne.Norbit) then
              smlib(ibin,2*iorb-1)=slush(ibin+Nvlib-1)
              smlib(ibin,2*iorb)=slush(ibin+Nvlib-1)
           else
              smlib(ibin,iorb)=slush(ibin+Nvlib-1)
           endif
        enddo
      enddo

      close(71)

      RETURN
      END
