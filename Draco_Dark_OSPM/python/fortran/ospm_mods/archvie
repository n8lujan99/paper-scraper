C-----------------------------------------------------------------------------D
C     subroutine VANIREAD reads the velocity stuff
C
C     USED BY VLOOK
C
C-----------------------------------------------------------------------------D
      SUBROUTINE highervaniread()
      INCLUDE 'moddefs.h'
      real vrpliba(Nrani,Ntani),vrliba(Nrani,Ntani),vtliba(Nrani,Ntani)
      real vrtpliba(Nrani,Ntani),vtpliba(Nrani,Ntani)

      open(unit=11,file='highervani.out',status='old')
      do iorb=1,Norb
         read(11,*) vrliba,vtliba,vrtpliba,vrpliba,vtpliba
         do ir=1,Nrani
            do it=1,Ntani
               if(Norb.ne.Norbit) then
                  vrpa(ir,it,2*iorb-1)=vrpliba(ir,it)
                  vrpa(ir,it,2*iorb)=vrpliba(ir,it)
                  vtpa(ir,it,2*iorb-1)=vtpliba(ir,it)
                  vtpa(ir,it,2*iorb)=vtpliba(ir,it)
                  vrtpa(ir,it,2*iorb-1)=vrtpliba(ir,it)
                  vrtpa(ir,it,2*iorb)=vrtpliba(ir,it)
                  vta(ir,it,2*iorb-1)=vtliba(ir,it)
                  vta(ir,it,2*iorb)=vtliba(ir,it)
                  vra(ir,it,2*iorb-1)=vrliba(ir,it)
                  vra(ir,it,2*iorb)=vrliba(ir,it)
               else
                  vrpa(ir,it,iorb)=vrpliba(ir,it)
                  vtpa(ir,it,iorb)=vtpliba(ir,it)
                  vrtpa(ir,it,iorb)=vrtpliba(ir,it)
                  vta(ir,it,iorb)=vtliba(ir,it)
                  vra(ir,it,iorb)=vrliba(ir,it)
               endif
            enddo
         enddo
      enddo
      close(11)


      RETURN
      END
