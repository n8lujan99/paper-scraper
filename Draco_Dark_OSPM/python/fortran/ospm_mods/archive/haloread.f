      subroutine haloread()
      include 'libdefs.h'

      open(unit=72,file='halo.dat',status='old')
      read(72,*) ihalo
      read(72,*) dis
      read(72,*) v0,rc,qdm
      read(72,*) xmgamma,rsgamma,gamma
      read(72,*) cnfw,rsnfw
      close(72)

      open(unit=82,file='gden.norm',status='old')
      read(82,*) gdennorm
      close(82)

      return
      end
