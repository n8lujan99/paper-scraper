C==============================================================================D
C     function ENERGY calculates the energy of an orbit given r, th,
C       vr, vth, and xLz
C
C     USED BY LIBRARY
C
C==============================================================================D
      FUNCTION energy(r,th,vr,vth)
      INCLUDE 'libdefs.h'

      v = sin(th)
      vco = cos(th)
      vph = xLz/r/vco

      call potential(r,v,VV)

      energy = 0.5*(vr*vr+vth*vth+vph*vph)+VV

      RETURN
      END
