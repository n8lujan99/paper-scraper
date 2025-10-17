C==============================================================================D
C     functions which convert between bin indices and physical quantities
C       IRneeR -- converts physical r (on [0,1]) into radial bin index ir
C       RneeIR -- converts radial bin index ir into physical r
C       IVneeV -- converts physical v=sin(th) into angular bin index iv
C       VneeIV -- converts angular bin index iv into physical v = sin(th)
C
C     Note: functions returning physical quantities from bin indices return
C       values corresponding to the BIN CENTER
C
C     USED BY LIBRARY AND MODEL
C
C==============================================================================D

      FUNCTION IRneeR(r)
      include 'bothdefs.h'
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IRneeR = int( log(a*r/b+cval)/a + 0.5 )
      RETURN
      END

      FUNCTION RneeIR(ir)
      include 'bothdefs.h'
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      RneeIR = b/a * (exp(a*float(ir))-cval)
      RETURN
      END

      FUNCTION IVneeV(v)
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IVneeV = int(sign(1.,v))*min(nint(float(ivmax)*abs(v)+0.5),ivmax)
      RETURN
      END

      FUNCTION VneeIV(iv)
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      VneeIV = (float(iv-1)+0.5)/float(ivmax)
      RETURN
      END

      FUNCTION IRCneeR(r)
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IRCneeR = (IRneeR(r)-1)/irrat + 1
      RETURN
      END

      FUNCTION IVCneeV(v)
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IVCneeV = (IVneeV(v)-1)/ivrat + 1
      RETURN
      END

      function RneeIRC(irc)
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      RneeIRC = .5*(RneeIR((irc - 1)*irrat + 1)+
     &     RneeIR((irc-1)*irrat + irrat))
      return
      end

      function VneeIVC(ivc)
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      VneeIVC = .5*(VneeIV((ivc - 1)*ivrat + 1)+
     &     VneeIV((ivc-1)*ivrat + ivrat))
      return
      end

      function RneeIRhalo(ir)
      include 'libdefs.h'
      RneeIRhalo=10**(rlgpotmin+(rlgpotmax-rlgpotmin)
     &     /float(npot-1)*(ir-1))
      return
      end

      function IRneeRhalo(r)
      include 'libdefs.h'
      IRneeRhalo=nint((log10(r)-rlgpotmin)/
     &     (rlgpotmax-rlgpotmin)*float(npot-1)+1.)
      return
      end

      FUNCTION RneeIRapo(ir)
      include 'libdefs.h'
      RneeIRapo = belz/aelz * (exp(aelz*float(ir))-celz)
      RETURN
      END

      FUNCTION RneeIRperi(ir)
      include 'libdefs.h'
      RneeIRperi = belz/aelz * (exp(aelz*float(ir))-celz)
      RETURN
      END
      
      FUNCTION IRneeRapo(r)
      include 'libdefs.h'
      IRneeRapo = int( log(aelz*r/belz+celz)/aelz + 0.5 )
      RETURN
      END

      FUNCTION IRneeRperi(r)
      include 'libdefs.h'
      IRneeRperi = int( log(aelz*r/belz+celz)/aelz + 0.5 )
      RETURN
      END
