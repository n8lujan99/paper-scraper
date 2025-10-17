      subroutine mlread()
      INCLUDE 'moddefs.h'
      open(unit=11,file='ml.out',status='old')
      read(11,*) xmu
      read(11,*) alp,ent,chi
      read(11,*) x1,x2,ent,chi
      xmu=1./sqrt(xmu)
      close(11)
      return
      end
