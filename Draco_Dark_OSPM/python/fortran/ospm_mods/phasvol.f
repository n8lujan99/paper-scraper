C==============================================================================D
C     subroutine PHASVOL calculates the phase-space volume occupied
C       by an orbit by punching out the phase area internal to the
C       present orbit
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE phasvol(iorblo,iorbup,dLz,dE)
      INCLUDE 'libdefs.h'
      parameter(nmaximal=1000)
      DIMENSION wtemp(nmaximal),xitemp(nmaximal)
c      DIMENSION wtemp(Nvdat),xitemp(Nvdat)

333   format ('Vols: ',$)
335   format (1pe8.2,' ',$)

      if (iorbup.eq.iorblo) then
        if (wphase(iorblo).eq.7.77e18) then
          wphase(iorblo) = wphase(iorblo-1)
        else
          wphase(iorblo) = abs(wphase(iorblo)*dE*dLz)
        endif
        write (50,333)
        write (6,333)
        write (50,335)wphase(iorblo)
        write (6,335)wphase(iorblo)

      else

        ntemp = iorbup-iorblo+1

        write (50,333)
        write (6,333)


C------ load wphase() into wtemp() and store indices, for ntemp points

        do i=iorblo,iorbup
          ii = i-iorblo+1
          wtemp(ii) = abs(wphase(i))
          xitemp(ii) = float(ii)
        enddo


C------ fill up the unused portion of wtemp() with 1e19 for sort2() to chew on

        do ii=ntemp+1,nmaximal
          wtemp(ii) = 1.e19
          xitemp(ii) = float(ii)
        enddo
        

        call sort2(nmaximal,wtemp,xitemp)


c        do ii=ntemp+1,Nvdat
c          wtemp(ii) = 1.e19
c          xitemp(ii) = float(ii)
c        enddo

c        call sort2(Nvdat,wtemp,xitemp)


C------ check for bananas

        do ii=2,ntemp
          if (wtemp(ii).eq.7.77e18) wtemp(ii)=wtemp(ii-1)
        enddo


C------ difference adjacent annuli

        iflag=0
        do ii=ntemp,2,-1
          if (wtemp(ii).eq.wtemp(ii-1)) then
            iflag=ii
            goto 111
          else
            wtemp(ii) = abs(wtemp(ii)-wtemp(ii-1))
          endif
        enddo

111     if (iflag.ne.0) then
          do ii=iflag,ntemp
            wtemp(ii)=wtemp(ii-1)
          enddo
        endif

        do ix=1,iorbup-iorblo+1
          ii = ix-1+iorblo
          wphase(ii) = abs(wtemp(ix)*dE*dLz)
        enddo

        do ix = iorblo,iorbup
          write (50,335)wphase(ix)
          write (6,335)wphase(ix)
        enddo

      endif

      write (50,*)''
      write (6,*)''

      RETURN
      END
