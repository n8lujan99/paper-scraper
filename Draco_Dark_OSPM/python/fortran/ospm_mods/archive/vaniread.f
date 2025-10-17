C-----------------------------------------------------------------------------D
C     subroutine VANIREAD reads the velocity stuff
C
C     USED BY VLOOK
C
C-----------------------------------------------------------------------------D
      SUBROUTINE vaniread(v0a,vr2a,vt2a,vp1a,vp2a)
      INCLUDE 'moddefs.h'
      real v0liba(Nrani,Ntani),vr2liba(Nrani,Ntani),vt2liba(Nrani,Ntani)
      real vp1liba(Nrani,Ntani),vp2liba(Nrani,Ntani)
      real vrtliba(Nrani,Ntani)
      real v0a(Nrani,Ntani,Norbitm),vr2a(Nrani,Ntani,Norbitm)
      real vt2a(Nrani,Ntani,Norbitm),vp1a(Nrani,Ntani,Norbitm)
      real vp2a(Nrani,Ntani,Norbitm)

      open(unit=11,file='vani.out',status='old')
      do iorb=1,Norb
         read(11,*) v0liba,vr2liba,vt2liba,vrtliba,vp1liba,vp2liba
c         read(11,*) v0liba,vr2liba,vt2liba,vp1liba,vp2liba
         do ir=1,Nrani
            do it=1,Ntani
               if(Norb.ne.Norbit) then
                  v0a(ir,it,2*iorb-1)=v0liba(ir,it)
                  v0a(ir,it,2*iorb)=v0liba(ir,it)
                  vr2a(ir,it,2*iorb-1)=vr2liba(ir,it)
                  vr2a(ir,it,2*iorb)=vr2liba(ir,it)
                  vt2a(ir,it,2*iorb-1)=vt2liba(ir,it)
                  vt2a(ir,it,2*iorb)=vt2liba(ir,it)
c                  vrta(ir,it,2*iorb-1)=vrtliba(ir,it)
c                  vrta(ir,it,2*iorb)=vrtliba(ir,it)
                  vp1a(ir,it,2*iorb-1)=vp1liba(ir,it)
                  vp1a(ir,it,2*iorb)=-vp1liba(ir,it)
                  vp2a(ir,it,2*iorb-1)=vp2liba(ir,it)
                  vp2a(ir,it,2*iorb)=vp2liba(ir,it)
               else
                  v0a(ir,it,iorb)=v0liba(ir,it)
                  vr2a(ir,it,iorb)=vr2liba(ir,it)
                  vt2a(ir,it,iorb)=vt2liba(ir,it)
c                  vrta(ir,it,iorb)=vrtliba(ir,it)
                  vp1a(ir,it,iorb)=vp1liba(ir,it)
                  vp2a(ir,it,iorb)=vp2liba(ir,it)
               endif
            enddo
         enddo
      enddo
      close(11)


      RETURN
      END
