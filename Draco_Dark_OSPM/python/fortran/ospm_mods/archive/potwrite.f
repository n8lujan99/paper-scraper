C==============================================================================D
C     subroutine POTWRITE samples the potential
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE potwrite()
      INCLUDE 'libdefs.h'

      open (unit=82,file='potential.out',status='unknown')

89    format (4(e12.6,2x))

      ddd = log(3.6/b)/999.
      ccc = b/3./exp(ddd)
      do irr=1,1000
        r = ccc*exp(ddd*float(irr))
        call potential(r,0.,VV0)
        call potential(r,.5,VV5)
        call potential(r,.9,VV9)
        write (82,89)r,VV0,VV5,VV9
      enddo

      r=0.5

      do ith=0,899
        th = float(ith)/10.*pi/180.
        v = sin(th)
        call potential(r,v,VV)
        write (82,89)180.*th/pi,-VV
      enddo

      close(82)

      open(unit=82,file='vcirc.out',status='unknown')
      xfac = sqrt(GG*totlight*arcsec/distance/angrad)
      xmin = log10(rpotmin)
      xmax = log10(rpotmax)
      ipnts = 200
      do ir=1,ipnts
         xr = 10**(xmin + (xmax-xmin)/float(ipnts-1)*float(ir-1))
         call force(xr,0.,xrforce,xvforce)
         write(82,*) ir,xr,sqrt(abs(xrforce)*xr)*xfac
      end do
      close(82)
      

      RETURN
      END
