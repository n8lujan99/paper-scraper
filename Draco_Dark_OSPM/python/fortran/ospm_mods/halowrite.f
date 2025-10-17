      subroutine halowrite()
      include 'libdefs.h'

      open(unit=72,file='halo.dat',status='unknown')
      write(72,*) ihalo
      write(72,*) dis
      write(72,*) v0,rc,qdm
      write(72,*) xmgamma,rsgamma,gamma
      write(72,*) cnfw,rsnfw
      close(72)

      return
      end
