C==============================================================================D
C     function EKIN calculates the kinetic energy of a particle given the
C       dynamical variables r and v=sin(th) (Etot in common block)
C
C     USED BY LIBRARY
C
C==============================================================================D
      FUNCTION ekin(r,v)
      INCLUDE 'libdefs.h'

      call potential(r,v,VV)
      ekin = Etot - VV - xLz*xLz/r/r/(1.-v*v)/2.

      RETURN
      END
