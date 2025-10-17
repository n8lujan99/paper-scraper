C==============================================================================D
C     subroutine AREA calculates the total (i.e. non-annular) s.o.s. area
C       of orbit iorb and puts it in wphase(iorb) [later converted to an
C       annular area in VOLUME()]
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE area(iorb,isos,xm)
      INCLUDE 'libdefs.h'
      real xtemp(Nsos),ytemp(Nsos)


C---- if we have one point it's a failed orbit:

      if(isos.eq.1) then
         wphase(iorb) = 7.77e18
         goto 999
      endif

      if (isos.lt.Nsos) then
        do i=isos+1,Nsos
          rsos(i) = 1.e12
          Vrsos(i) = 1.e12
        enddo
      endif

      xrmax=-1.e20
      xvrmax=-1.e20
      do i=1,isos
         xrmax=max(xrmax,rsos(i))
         xvrmax=max(xvrmax,Vrsos(i))
      end do

      do i=1,Nsos
         xtemp(i)=rsos(i)/xrmax
         ytemp(i)=Vrsos(i)/xvrmax
      end do

C---- rearrange rsos() and Vrsos() in order of increasing r
      call sort2(Nsos,rsos,Vrsos)
      call sort2(Nsos,xtemp,ytemp)

      ize=0
      xsum=0.
      xsum2=0.
      do i=2,isos-1
         x1=(ytemp(i)-ytemp(i-1))/(xtemp(i)-xtemp(i-1))
         x2=(ytemp(i+1)-ytemp(i))/(xtemp(i+1)-xtemp(i))
c         x3=0.5*(xtemp(i+1)-xtemp(i-1))
         x3=1.
         xsum2=xsum2+x1*x1+x2*x2
         xsum=xsum+(x2-x1)/x3
      end do

c      print*
c      print*,'XSUM2',xsum2,xsum
c      print*

      xm=xsum/float(isos-2)

      if(xm.gt.-1.5) then 
         wphase(iorb) = 7.77e18
         goto 999
      end if
         
C---- calculate the fourth moment divided by the second moment = xm
c      call sosmom(Nsos,xm,rsos,Vrsos,isos)
      
C---- if xm < 0.75 it's a single-island s.o.s.

c      if (xm.lt.0.55) then
      do i=2,isos
         wphase(iorb) = wphase(iorb) + (min(Vrsos(i-1),Vrsos(i)) +
     &        abs(Vrsos(i)-Vrsos(i-1))/2.) * (rsos(i)-rsos(i-1))
      enddo

C---- if xm > 0.75 it's a weird orbit

C     (give it an area of 7.77e18 and later in phasvol() we'll convert
C     it to the previous orbit's phase volume)

c      else
c         wphase(iorb) = 7.77e18
c      endif

 999  continue

      RETURN
      END
