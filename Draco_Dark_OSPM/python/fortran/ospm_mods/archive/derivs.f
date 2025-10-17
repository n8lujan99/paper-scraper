C==============================================================================D
C     subroutine DERIVS calculates the derivatives used by RK4 using the
C       gradient of the Legendre polynomial-expanded potential
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE derivs(pos0,dposdt)
      INCLUDE 'libdefs.h'
      DIMENSION pos0(Nvar),dposdt(Nvar)

      r = pos0(1)
      th = pos0(2)
      v = sin(th)
      vr = pos0(3)
      vth = pos0(4)
      vph = xLz/r/cos(th)

      call force(r,v,frL,fvL)

      dposdt(1) = vr
      dposdt(2) = vth/r
      dposdt(3) = (vth*vth+vph*vph)/r + frL
      dposdt(4) = -tan(th)*vph*vph/r - vth*vr/r + fvL

      RETURN
      END
