
      INCLUDE 'moddefs.h'
      real sumb(Nrlib,Nvlib,Nrlib,Nvlib),sumbn(Nrlib,Nvlib,Nrlib,Nvlib)
      real slush(Nvel*Nrlib*Nvlib)
      real v1libt(Nvel,Nrlib,Nvlib),v1libt2(Nvel,Nrlib,Nvlib)
      character filen(Nveld)*40,galname*40,fileo(9)*40

      fileo(1)='v1libsc1.out'
      fileo(2)='v1libsc2.out'
      fileo(3)='v1libsc3.out'
      fileo(4)='v1libsc4.out'
      fileo(5)='v1libsc5.out'
      fileo(6)='v1libsc6.out'
      fileo(7)='v1libsc7.out'
      fileo(8)='v1libsc8.out'
      fileo(9)='v1libsc9.out'

      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call galreadm(galname)
      call velbin(filen)

      open(unit=68,file='norbit',status='old')
      read(68,*) Norb
      close(68)
      print *,'Norbit = ',Norb
      if(Norb.gt.Norbm) print *,'make Norbm bigger'
      if(Norbm.ne.Norbitm) then
         Norbit=2*Norb
      else
         Norbit=Norb
      endif

      open(unit=2,file='see.dat',status='old')
      read(2,*)
      do i=1,Nstot
         read(2,*,end=667) i1,x2
         if(i1.ne.i) print *,'Order see.dat or else....'
         seeb(i1)=x2/angrad
      enddo
 667  continue
      close(2)

      ifile=0
      do is=1,Nstot
         if(seeb(is).gt.0.) then
            open (unit=1,file='v1lib.out',status='old')
            ifile=ifile+1
            open (unit=2,file=fileo(ifile),status='unknown')
            call getsum(is,seeb(is),sumb,sumbn)
            do iorb=1,Norb
               read(1,*) slush
               do ir=1,Nrlib
                  do iv=1,Nvlib
                     do ivel=1,Nvel
                        ic=(ir-1)*Nvlib*Nvel+(iv-1)*Nvel+ivel
                        v1libt(ivel,ir,iv)=slush(ic)
                     enddo
                  enddo
               enddo
               call convolve(sumb,sumbn,v1libt,v1libt2)
               do ir=1,Nrlib
                  do iv=1,Nvlib
                     sum1=0.
                     sum2=0.
                     do ivel=1,Nvel
                        ic=(ir-1)*Nvlib*Nvel+(iv-1)*Nvel+ivel
                        slush(ic)=v1libt2(ivel,ir,iv)
                        sum1=sum1+v1libt(ivel,ir,iv)
                        sum2=sum2+v1libt2(ivel,ir,iv)
                     enddo
                     if(float(iorb)/1000.0.eq.iorb/1000) then
                     write(*,'("Orbit number ",i5,a1,$)') iorb,char(13)
                     call flush(6)
                     endif
                  enddo
               enddo
               write(2,*) slush
            enddo
            close(1)
            close(2)
         endif
      enddo
      write(*,*)

      end
