C=============================================================================D
C     program MODEL calculates the weights of the orbits in xlib.dat such
C       that their assembly best matches the 2-D profile while maximizing
C       the entropy of the configuration
C=============================================================================D
      PROGRAM projmass
      INCLUDE 'moddefs.h'
      parameter(ifrac=2,nadt=1000)
      DOUBLE PRECISION S,dS(Norbitm+Nvel*Nvelbm)
      DOUBLE PRECISION ddS(Norbitm+Nvel*Nvelbm),rcond
      REAL tarray(2),times(1000),losvd(1000),losvdd(1000)
      real frac(ifrac),fwda(Nvelbm,ifrac),velin(nadt)
      character filen(Nveld)*40,galname*40
      character(len=4) iternam(500)

      integer itercount,i2

      data big,idone,igodef /1.e20,0,0/
      data frac /0.25,0.50/
      data fracrv,fracr1v /0.02,0.3/

      real rotfrac,wnorm,alphat,xmutarget
      integer ians

      itercount=0

      do i=1,199
         if(i.lt.100) then
            if(i.lt.10) then
               write(iternam(i),'(".00",i1)') i
            else
               write(iternam(i),'(".0",i2)') i
            end if
         else
            write(iternam(i),'(".",i3)') i
         end if
      end do

      time = etime(tarray)

      alphat = 0.

      ratioold=big
      ratio1old=big

C-- entropy types (type 0 and 1 have never worked well?), where
c   w is orbit weight, and p is 1/(phase volume) :
C   0 => S=w*w
C   1 => S=log(w*p)/p
C   2 => S=-w*log(w*p)

      iS = 2

C-- read in all input files
 
 10   format (' -> reading in orbit library and parameter files...')
      write(6,10)

      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call totmlread()
      call galreadm(galname)
      if(alphat.eq.0) then
         open(unit=63,file='gden.norm',status='old')
         read(63,*) x1
         read(63,*) xmutarget
         xmutarget = xmutarget*x1
         close(63)
         open(unit=94,file='gden.norm',status='old')
         read(94,*) xgdennorm !normalization
         read(94,*) xsml
         close(94)
         xsml=xsml*xgdennorm
         xmutarget=xmutarget*xgdennorm/xsml
         xmutarget=1./sqrt(xmutarget)
      endif
      call phaseread()
      call velbin(filen)
      call m3coarse()
      call smcoarse()
      call iterread()
      call xlibread()
      call d3read()
      call spacel()
      call vlibread()
      call vdataread(filen)

c-- fit either surface brightness or 3d density

      if(ifit.eq.1) then
         do ibin=1,Nbin
            do iorb=1,Norbit
               xlib(ibin,iorb)=smlib(ibin,iorb)
            enddo
            SMc(ibin)=SMcorig(ibin)
         enddo
      else
         do ibin=1,Nbin
            do iorb=1,Norbit
               xlib(ibin,iorb)=d3lib(ibin,iorb)
            enddo
            SMc(ibin)=d3c(ibin)
         enddo
      endif         

c-- if alpha=0 then use new weights, otherwise use old weights
      
      inweights=0
      if(alphat.ne.0) inweights=1

 11   format(50x,' .. done.')
      write(6,11)
      write(6,*)''

C-- initialize weights

      if (inweights.eq.0) then

c --- original version
c         do i=1,Norbit
c            w(i) = 1./float(Norbit)
c         end do
c --- end

c --- new version
         ! to check the orbit sampling ...
         wnorm=0.
         rotfrac = 0.5
         do i=1,Norb
            w(2*i-1) = 1./wphase(2*i-1) * rotfrac
            w(2*i) = 1./wphase(2*i) * (1.-rotfrac)
            wnorm = wnorm + w(2*i-1) + w(2*i)
         enddo
         do i=1,Norbit
            w(i) = w(i)/wnorm
         end do
         open(unit=92,file='weights.000',status='unknown')
         do i12=1,Norbit
            write(92,*) i12,w(i12)
         end do
         close(92)
         ! ... now set the initial weights
         wnorm=0.
         rotfrac = 0.75
         do i=1,Norb
            w(2*i-1) = 1./wphase(2*i-1) * rotfrac
            w(2*i) = 1./wphase(2*i) * (1.-rotfrac)
            wnorm = wnorm + w(2*i-1) + w(2*i)
         enddo
         do i=1,Norbit
            w(i) = w(i)/wnorm
         end do
c --- end


c-- set the xmu to the initial guess and set the data losvd's

        xmu=xmutarget
        call velmtod()

        sxmu=0.
        do irc=1,Nvelb
           do ivel=1,Nvel
              i=Norbit+(irc-1)*Nvel+ivel
              sum = 0.
              do j=1,Norbit
                 sum=sum+v1lib(ivel,irc,j)*w(j)
              enddo
              losvd(ivel)=sum            
              w(i)=sumad(ivel,irc)-sum
           enddo

           do ivel=1,nad
              losvdd(ivel)=ad(ivel,irc)
              velin(ivel)=veld(ivel,irc)
           enddo
           do i=1,ifrac
              call getfwhm(Nvel,velm,losvd,frac(i),fwm,xd1,xd2)
              call getfwhm(nad,velin,losvdd,frac(i),fwd,xd1,xd2)
              sxmu=sxmu+fwm/fwd
              fwda(irc,i)=fwd
           enddo
        enddo
      else
         call weightread()
         open(unit=11,file='ml.out',status='old')
         read(11,*) xmu
         xmu=1./sqrt(xmu)
         close(11)
         do irc=1,Nvelb
            do ivel=1,nad
               losvdd(ivel)=ad(ivel,irc)
               velin(ivel)=veld(ivel,irc)
            enddo
            do i=1,ifrac
               call getfwhm(nad,velin,losvdd,frac(i),fwd,xd1,xd2)
               fwda(irc,i)=fwd
            enddo
         enddo
      endif

      write(6,11)
      write(6,*)''

 12   format(' -> calculating Cm()...')
      write(6,12)

C-- calculate expanded "library" matrix Cm() used in constructing Am

      call cmatrix()

 13   format(' -> iteration:')
      write(6,13)

C-- write out zeroth iteration stuff

      iter=0
      call velmtod()
      call entropy(S,dS,ddS,iS)
      call comparewrite(iS, 6,0,iter,rcond)
      if(alphat.gt.0) print *,xmu,sxmu,pml,chi/alphat,alphat

      times(1) = dtime(tarray)

c-- open a file for the iteration data
      
      open(unit=13,file='iteration.out',status='unknown')

C-- iteration loop

 766  continue

c --- write out weights of iteration
      itercount=itercount+1
      open(unit=92,file='niteration.out',status='unknown')
      write(92,*) itercount
      close(92)
      open(unit=92,file='weights'//iternam(itercount),status='unknown')
      do i12=1,Norbit
         write(92,*) i12,w(i12)
      end do
      close(92)
c --- end


      do iter=1,Niter

        call velmtod()
        call entropy(S,dS,ddS,iS)

c-- spear is where all the cpu is being spent

        call spear(S,dS,ddS,rcond)

        write(6,11)
        write(6,*)''

C-- write out iteration information

        call comparewrite(iS, 6,1,iter,rcond)

C-- get the new xmu

        sxmu=0.
        istart=1
        do irc=istart,Nvelb
           do ivel=1,Nvel
              sum=0.
              do j=1,Norbit
                 sum=sum+v1lib(ivel,irc,j)*w(j)
              enddo
              losvd(ivel)=sum
              losvdd(ivel)=sumad(ivel,irc)
           enddo
           do i=1,ifrac
              call getfwhm(Nvel,velm,losvd,frac(i),fwm,xd1,xd2)
              call getfwhm(Nvel,velm,losvdd,frac(i),fwd,xd1,xd2)
              sxmu=sxmu+fwm/fwd
           enddo
        enddo
        sxmu=sxmu/float((Nvelb-istart+1)*ifrac)
        pml=xmu+xmu*(sxmu-1.)
        xmu=xmu+apfacmu*xmu*(sxmu-1.)
        if(alphat.eq.0) xmu=xmutarget

c-- check if everything is ok for this iteration

        if(alphat.gt.0) then
           chival=chi/alphat/float(Nvelb*Nvel)
           if(pml.gt.0) then
              pml=1./pml**2
           else
              pml=666
           endif
           print *,xmu,sxmu,pml,chi/alphat,alphat
           if(chival.lt.0.0005) then
              write(*,*) 'Hit the minimum chival'
              idone=1
              goto 1666
           endif
           if(chival.gt.50.and.alphat.ge.1.e-3) goto 666
           if(xlbad.gt.0.01.and.iter.gt.(Niter-1)
     $          .and.alphat.gt.3.e-4) goto 666
        endif

        if(1.d0+rcond.eq.1.d0) print *,'Possible Singular Matrix !!!'

        times(iter) = dtime(tarray)

c-- check if the chi^2 is changing, if not then get out

        if(alphat.gt.0.and.xlbad.lt.0.01) then
           ratio1=chi/alphat
           fracr1=abs(ratio1old-ratio1)
           ratio1old=ratio1
           if(fracr1.le.fracr1v) goto 1666
        endif

      enddo

C-- end iteration loop

C-- write out files

 1666 continue
 14   format (' -> writing out files..')
      write (6,14)

c-- write out the final iteration information

      open(unit=11,file='ml.out',status='unknown')
      write(11,*) 1./xmu**2
      write(11,*) alphat,ent,chi
      ratio=0.
      if(alphat.ne.0) ratio=chi/alphat
      write(11,*) hole,xinclin*180./pi,ent+chi,ratio,sxmu,xlbad
      close(11)
      open(unit=11,file='chi.out',status='unknown')
      do ivel=1,Nvelb
         write(11,*) rbin(ivel),chi1(ivel)
      enddo
      close(11)

c-- write out the current iteration information

      write(13,*) alphat,ratio,ent,ent+chi

c-- check is chi^2 is changing, if not then quit

      if(ratio.gt.0.and.alphat.gt..01) then
         fracr=abs(ratioold-ratio)
         ratioold=ratio
         if(fracr.lt.1.5) fracr1v=fracrv
         if(fracr.lt.fracrv) idone=1
      endif

c-- write out the weights

      call weightwrite()
      open (unit=73,file='compare.out',status='unknown')
      call comparewrite(istype,73,1,iter-1,rcond)
      close(73)

      avt=0.
      do i=1,iter
        avt=avt+times(i)
      enddo
      write(6,*)' '
      avt=avt/float(iter)
      time=etime(tarray)

      call iterwrite(iter,S,time,avt,Nbin)

      write(6,11)
      write(6,*)''
      write(6,*) 'Alpha used : ',alphat
      write(6,*)

      if(idone.eq.0) then

c-- calculate the surface brightness from the model projected light

         if(alphat.le.1.1e-10) then
            do ibin=1,Nbin
               sum=0.
               sum2=0.
               do iorb=1,Norbit
                  sum=sum+w(iorb)*smlib(ibin,iorb)
                  sum2=sum2+w(iorb)*d3lib(ibin,iorb)
               enddo
               print *,ibin,sum/SMcorig(ibin),sum2/d3c(ibin)
               SMcorig(ibin)=sum
            enddo
            ians = 1

c-- write out the new surface brightness values if desired

            if(ians.eq.1) then
               open(unit=11,file='SM.dat',status='unknown')
               do ibin=1,Nbin
                  write(11,*) itor(ibin),itov(ibin),
     $                 SMcorig(ibin)*totmass,RneeIRC(itor(ibin))*angrad
               enddo
               close(11)
               goto 666
            endif
         endif

c-- iterate with new alpha

         if(igodef.eq.0) then
            write(*,"('Continue? (1-yes) : '$)")
            read *,ians
            if(ians.eq.1) then
               write(*,"('Input alpha : '$)")
               read *,ans

c-- if alpha=666, then do an incremental change from now on

               if(ans.eq.666) then
                  alphat=alphat*alphainc
                  igodef=1
               else
                  alphat=ans
               endif
               goto 766
            endif
         else
            alphat=alphat*alphainc
            goto 766
         endif
      endif

 666  continue
      close(13)

      END
