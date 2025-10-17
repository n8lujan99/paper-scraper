C=============================================================================D
C     subroutine VLIBREAD reads in the velocity libraries v0lib.out,
C       v1lib.out, and v2lib.out
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE vlibread()
      INCLUDE 'moddefs.h'
      real slush(Nvel*Nrlib*Nvlib),slush2(Nvel*Nrlib*Nvlib)
      real slusha(Nvel*Nrlib*Nvlib,Nstot)
      real nivbin(Nvelbm)
      character fileo(9)*40

      fileo(1)='v1libsc1.out'
      fileo(2)='v1libsc2.out'
      fileo(3)='v1libsc3.out'
      fileo(4)='v1libsc4.out'
      fileo(5)='v1libsc5.out'
      fileo(6)='v1libsc6.out'
      fileo(7)='v1libsc7.out'
      fileo(8)='v1libsc8.out'
      fileo(9)='v1libsc9.out'

      open(unit=2,file='see.dat',status='old')
      read(2,*)
      do i=1,Nstot
         read(2,*,end=667) i1,x2
         if(i1.ne.i) print *,'Order see.dat or else....'
         seeb(i1)=x2/angrad
      enddo
 667  continue
      close(2)

      do iorb=1,Norbit
         do ivbin=1,Nvelb
            do ivel=1,Nvel
               v1lib(ivel,ivbin,iorb)=0.
c               v1lib2(ivel,ivbin,iorb)=0.
            enddo
         enddo
      enddo

      open (unit=73,file='v1lib.out',status='old')
      izerofwhm=0
      ifile=0
      do i=1,Nstot
         if(seeb(i).gt.0) then
            ifile=ifile+1
            open (unit=73+ifile,file=fileo(ifile),status='old')
         else
            izerofwhm=1
         endif
      enddo

      do iorb=1,Norb
         read(73,*) slush
         do is=1,ifile
            read(73+is,*) slush2
            do i=1,Nvel*Nrlib*Nvlib
               slusha(i,is)=slush2(i)
            enddo
         enddo
         do ivbin=1,Nvelb
            do ir=1,Nrlib
               do iv=1,Nvlib
                  fracm=areakin(ivbin,itobin(ir,iv))
                  do ivel=1,Nvel
                     ic=(ir-1)*Nvlib*Nvel+(iv-1)*Nvel+ivel
                     if(seeb(isee(ivbin)).eq.0) then
                        value=slush(ic)
                     else
                        value=slusha(ic,isee(ivbin)-izerofwhm)
                     endif
                     if(Norb.ne.Norbit) then
                        v1lib(ivel,ivbin,2*iorb-1)=
     $                       v1lib(ivel,ivbin,2*iorb-1)+value*fracm
                        v1lib(Nvel-ivel+1,ivbin,2*iorb)=
     $                       v1lib(Nvel-ivel+1,ivbin,2*iorb)+value*fracm
c                        v1lib2(ivel,ivbin,2*iorb-1)=
c     $                       v1lib2(ivel,ivbin,2*iorb-1)+slush(ic)*fracm
c                        v1lib2(Nvel-ivel+1,ivbin,2*iorb)=
c     $                       v1lib2(Nvel-ivel+1,ivbin,2*iorb)+
c     $                       slush(ic)*fracm
                     else
                        v1lib(ivel,ivbin,iorb)=
     $                       v1lib(ivel,ivbin,iorb)+value*fracm
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(73)
      do i=1,ifile
         close(73+i)
      enddo

c - symmetrize losvd's if on minor axis

c      fracd=0.50
      fracd=1.
      do ivbin=1,Nvelb
         if(iminor(ivbin).eq.1) then
            if(isee(ivbin).eq.1) then
               fract=fracd
            else
               fract=1.
            endif
            do iorb=1,Norb
               do ivel=1,Nvel
                  v1=v1lib(ivel,ivbin,2*iorb-1)
                  v2=v1lib(ivel,ivbin,2*iorb)
c                  v1b=v1lib2(ivel,ivbin,2*iorb-1)
c                  v2b=v1lib2(ivel,ivbin,2*iorb)
                  v1lib(ivel,ivbin,2*iorb-1)=
     $                 (v1+v2*fract)/(1.+fract)
                  v1lib(ivel,ivbin,2*iorb)=
     $                 (v1*fract+v2)/(1.+fract)
c                  v1lib2(ivel,ivbin,2*iorb-1)=
c     $                 (v1b+v2b*fract)/(1.+fract)
c                  v1lib2(ivel,ivbin,2*iorb)=
c     $                 (v1b*fract+v2b)/(1.+fract)
               enddo
            enddo
         endif
      enddo

c      open(unit=13,file='junk',status='unknown')
c      do ivbin=1,Nvelb
c         do ivel=1,Nvel
c            write(13,*) ivbin,ivel,v1lib(ivel,ivbin,3367)
c         enddo
c      enddo
c      close(13)
c      print *,'done'

      RETURN
      END
