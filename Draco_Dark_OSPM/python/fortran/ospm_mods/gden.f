C This program calculates the mass in each binr for the library
C program using a surface brightness profile and specified flattening

      PROGRAM gden
      INCLUDE 'bothdefs.h'
      parameter(nmax=5000)
      real rd(nmax),td(nmax),d(nmax,nmax),s(nmax,nmax),mass(nmax,nmax)
      real den(nmax),radl(nmax),radu(nmax),thl(nmax),thu(nmax)
      real ratml(150,150),rin(nmax),tin(nmax)
      real w,sinw,rs,m200,rho0,r0
      real hparam,cparam,rhocrit,pizahl,r200,rhoakt,ratmin,ratmax
      real stml,cosw,rakt
      real rhoin(nmax,nmax),sbin(nmax,nmax)
      real r1,r2,r01,r21,th1,th2,th01,th21,d11,d21,d12,d22,f
      real f2,halomass(nmax,nmax),masscore,lightcore,halocore
      real normfac,stmass,inc
      integer ind1,ind2
      character filen*40
      external func,funci,func1,func2,func3,func4,funcelz
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat

      COMMON/relzbin/aelz,belz,relzmin,relzmax

      COMMON/halo/qdm,ihalo,cnfw,rsnfw,xmgamma,rsgamma,
     &     gamma,rc,v0,dis,gdennorm
      common/cfunci/axra,axrt
      common/cfunc1/thlo,thup
      common/cfunc1a/rin,tin,rhoin,s,rmax,mr,na,mass,halomass
      common/cfunc2a/sbin

      !defaults
      cnfw=0.
      rsnfw=1.0
      mgamma=0.
      rsgamma=1.0
      gamma=1.0
      v0=0.
      rc=1.0
c - Hubble Constant -
      hparam=70. 
c - end -

      rhocrit=2.7754e-11*hparam**2
      pizahl=3.141592654


      print*
      print*,'************'
      print*,'* NEW gden *'
      print*,'************'
      print*


c - get rmax
      open(unit=1,file='galaxy.params',status='old')
      read(1,*)
      read(1,*) filen,distance
      read(1,*) filen,rmax
      close(1)

      dis=distance

 1    write(*,"('Datafile : '$)")
      read *,filen
      open(unit=1,file=filen,status='old',err=1)
      print *,'Distance and Rmax is : ',distance,rmax
      write(*,"('Rmin : '$)")
      read *,rmin

      axra=1.0
      axrt=1.0
      write(*,"('Inclination : '$)")
      read *,inc
      inc=inc/180.0*pi
      write(*,"('BH mass : '$)")
      read *,bhmass

      stml=1.
      qdm=1.
c      call qr2('stellar mass-to-light and flatness of DM ',
c     &  'gden.def',stml,qdm)
      write(*,"('stellar mass-to-light and flatness of DM: '$)")
      read *, stml,qdm

c      call qi1('DM halo (1-Dehnen/Plummer 2-NFW 3-loghalo 4-none
c     &5-logtot) ','gden.def',ihalo)
      write(*,"('DM halo (1-Dehnen/Plummer 2-NFW 3-loghalo 4-none
     &5-logtot: '$)")
      read *, ihalo
c      ihalo=4

      if(ihalo.eq.1) then !gamma/Plummer
c         call qr3('M [Msun] R_s [kpc] and gamma ','gden.def',
c     &        xmgamma,rsgamma,gamma)
         write(*, *) 'gamma/Plummer'
         write(*,"('M [Msun] R_s [kpc] and gamma: '$)")
         read *,  xmgamma,rsgamma,gamma
         rsgamma = rsgamma / distance * 206.265 ![rsgamma] = arcsec
      end if
      if(ihalo.eq.2) then !NFW
c         call qr2('c and R_s [kpc] ','gden.def',cnfw,rsnfw)
         write(*,"('NFW c and R_s [kpc]: '$)")
         read *, cnfw,rsnfw
         rsnfw = rsnfw / distance * 206.265 ![rsnfw] = arcsec
      end if
      if(ihalo.eq.3.or.ihalo.eq.5) then !nonsingular isothermal sphere
c         call qr2('v_0 [km/s] and R_c [kpc] ','gden.def',v0,rc)
         print *, 'nonsingular isothermal sphere'
         write(*,"('v_0 [km/s] and R_c [kpc]: '$)")
         read *, v0, rc
         rc0 = rc   ! [rc0] = kpc
         rc = rc / distance * 206.265 ![rc] = arcsec
      end if
c      call savdef    !libquest routine

      irmax=Nrdat
      ivmax=Nvdat
      ircmax=Nrlib
      ivcmax=Nvlib

      print *,'inclination = ',inc*180./pi
      open(unit=11,file='gal.dat',status='unknown')
      write(11,*) 'BHmass',bhmass
      write(11,*) 'inclination',inc*180./pi
      write(11,*) 'Stellar_ML_ratio',stml
      if( ihalo.eq.3.or.ihalo.eq.5 ) then
        write(11,*) 'Halo_v_0[km/s]', v0
        write(11,*) 'Halo_Rc[kpc]', rc0
        write(11,*) 'Halo_Rc[arcsec]', rc
      end if
      write(11,*) 'DMhalo_flattening', qdm

      close(11)

      irrat=irmax/ircmax
      ivrat=ivmax/ivcmax
      rmin=rmin/rmax

      a=rtbiso(func,.0001,.9,1.e-5)
      b=a*rmin/(exp(a)-cval)
 
      open(unit=14,file='ab.dat',status='unknown')
      write(14,*) a,b
      close(14)

      irelzmax=Nrelz
      relzmin=xorbitmin/rmax
      relzmax=xorbitmax/rmax
      relzmax=max(relzmax,1.2)

      aelz=rtbiso(funcelz,.0001,.9,1.e-5)
      belz=aelz*relzmin/(exp(aelz)-celz)
 
      open(unit=14,file='abelz.dat',status='unknown')
      write(14,*) aelz,belz
      write(14,*) relzmin,relzmax
      close(14)
 
      print *,'rmax = ',rmax*b/a*(exp(a*irmax)-1.0)
      print *,'outer-to-inner ratio = ',
     $     (exp(a*(irmax-1))-exp(a*(irmax-5)))/
     $     (exp(a*6.0)-exp(a*2.0))

      step=0
      xold=-1.
      do i=1,nmax
         read(1,*) xnew
         if(i.gt.1.and.xnew.ne.xold) goto 2400
         xold=xnew
         step=step+1
      enddo
 2400 continue
      close(1)

      nstep=nint(step)
c      step=2
c      td(1)=0.
c      td(2)=pi/2.
      open(unit=1,file=filen,status='old',err=1)
      mr=0
      do i=1,nmax
         mr=mr+1
         na=0
c         read(1,*,end=666) x1,x3,x4
         do j=1,nstep
            na=na+1
            read(1,*,end=666) x1,x2,x3,x4
            td(na)=x2            !theta (rad)
            tin(na)=td(na)
c            d(mr,na)=log10(1.e1*exp(x3))   !(ln(1e10Lsun/kpc3) -> Lsun/pc3)
c            s(mr,na)=log10(1.e4*exp(x4))   !(ln(1e10Lsun/kpc2) -> Lsun/pc2)
            d(mr,na)=log10(x3)   !(ln(1e10Lsun/kpc3) -> Lsun/pc3)
            s(mr,na)=log10(x4)   !(ln(1e10Lsun/kpc2) -> Lsun/pc2)
            call halodens(x1,x2,xaddmass)
c            xaddmass=0.
            mass(mr,na)=log10(10**(d(mr,na))*stml+xaddmass)
            if(ihalo.ne.4) then
               halomass(mr,na)=log10(xaddmass)
            else
               halomass(mr,na)=0.
            end if
c            write(*,'("reading ",i3,i3,a1,$)') i,j,char(13)
         enddo
         rd(mr)=log10(x1)               !radius (arcsec)
         rin(mr)=log10(x1)
      enddo
      
 666  continue
      close(1)

      mr=mr-1
      na=step

      call getbinp(irmax,ivmax,rmin,rmax,radl,radu,thl,thu,Nrelz)

      do i=1,mr
         do j=1,na
            rhoin(i,j)=d(i,j)
            sbin(i,j)=s(i,j)
         end do
      end do

c --- interpolate density on the grid
      do ir=1,irmax
c         rlo=radl(ir)
c         rup=radu(ir)
         rlo=(RneeIR(ir)+RneeIR(ir-1))/2.
         rup=(RneeIR(ir)+RneeIR(ir+1))/2.
         rakt=log10(RneeIR(ir)*rmax)
         do iv=1,ivmax
           w=(thl(iv)+thu(iv))/2.
           do ind1=2,mr
              if (rakt.lt.rd(ind1)) goto 2401
           end do
 2401      continue
           if (ind1.gt.mr) ind1=mr
           do ind2=2,na
              if (w.lt.td(ind2)) goto 2402
           end do
 2402      continue
           if (ind2.gt.na) ind2=na
           r1=rd(ind1-1)           !radius
           r2=rd(ind1)
           th1=td(ind2-1)          !theta
           th2=td(ind2)
           r01=rakt-r1
           r21=r2-r1
           th21=th2-th1
           th01=w-th1
           
           d11=rhoin(ind1-1,ind2-1)       !luminosity
           d21=rhoin(ind1,ind2-1)
           d12=rhoin(ind1-1,ind2)
           d22=rhoin(ind1,ind2)
           
           f = d11+(d21-d11)/r21*r01+(d12-d11)/th21*th01+
     &          ((d22-d12-d21+d11)/(th21*r21))*(th01*r01)

           d(ir,iv)=f
           
        end do
      end do

      do ir=1,irmax
         rlo=(RneeIR(ir)+RneeIR(ir-1))/2.
         rup=(RneeIR(ir)+RneeIR(ir+1))/2.
         rd(ir)=log10(RneeIR(ir)*rmax)
      end do

      do iv=1,ivmax
         td(iv)=(thl(iv)+thu(iv))/2.
      end do
c --- end
      open(unit=11,file='dL.dat',status='unknown')
      open(unit=12,file='SB.dat',status='unknown')
      open(unit=13,file='ratML.dat',status='unknown')
      open(unit=24,file='dM.dat',status='unknown')
      open(unit=25,file='dMhalo.dat',status='unknown')
      
      fac1=4.*pi  ! galaxy has TWO halves -> 2*2*pi
      fac2=(distance*1.e6/arcsec*rmax)**3
      fac3=(distance*1.e6/arcsec*rmax)**2
      sum1=0.
      sum2=0.

      do ir=1,irmax
c         rlo=radl(ir)
c         rup=radu(ir)
         rlo=(RneeIR(ir)+RneeIR(ir-1))/2.
         rup=(RneeIR(ir)+RneeIR(ir+1))/2.
         do iv=1,ivmax
            thlo=thl(iv)
            thup=thu(iv)
            xvol=4./3.*pi *(rup*rup*rup-rlo*rlo*rlo)/float(ivmax)
c            sd=10**d(ir,iv)*xvol/fac1
            call qromb(func1,rlo,rup,sd)
            call qromb(func2,rlo,rup,ss)
            call qromb(func3,rlo,rup,smass)
            if(ihalo.ne.4) then
               if(ihalo.eq.3.and.v0.eq.0) then
                  shalo=0
               else
                  call qromb(func4,rlo,rup,shalo)
               endif
            else
               shalo=0.
            end if
            xmval=smass*fac1*fac2
            xhaloval=shalo*fac1*fac2
            dval=sd*fac1*fac2
            sval=ss*4.*fac3
            sum1=sum1+dval
            sum2=sum2+sval
            if(iv.eq.ivmax) den(ir)=sum1
            write(11,*) ir,iv,dval
            write(12,*) ir,iv,sval
            write(24,*) ir,iv,xmval
            write(25,*) ir,iv,xhaloval
            ratml(ir,iv)=xmval/dval
         enddo
      enddo
      print *,sum1,sum2,sum1/sum2
      stmass = sum1*stml

      open(unit=31,file='coredat.out',status='unknown')
      rlo = 1.e-6
      rup=(RneeIR(1)+RneeIR(0))/2.
      masscore = 0.
      halocore = 0.
      lightcore= 0.
      print*,rlo,rup
      do iv=1,ivmax
         thlo=thl(iv)
         thup=thu(iv)
         xvol=4./3.*pi *(rup*rup*rup-rlo*rlo*rlo)/float(ivmax)
         call qromb(func1,rlo,rup,sd)
         call qromb(func2,rlo,rup,ss)
         call qromb(func3,rlo,rup,smass)
         if(ihalo.ne.4) then
            if(ihalo.eq.3.and.v0.eq.0) then
               shalo=0
            else
               call qromb(func4,rlo,rup,shalo)
            endif
         else
            shalo=0.
         end if
         masscore=masscore+smass*fac1*fac2
         halocore=halocore+shalo*fac1*fac2
         lightcore=lightcore+sd*fac1*fac2
      enddo
      write(31,*) masscore,halocore,lightcore

      ratmax=-1.e20
      ratmin=1.e20
      
      do ir=1,irmax
         do iv=1,ivmax
            ratmax=max(ratmax,ratml(ir,iv))
            ratmin=min(ratmin,ratml(ir,iv))
         end do
      end do

c      if(bhmass.ne.0.0) then
c         normfac=bhmass/100.
c      else
c         normfac=(ratmax+ratmin)/2.
c      end if
      normfac = 1.0

      do ir=1,irmax
         do iv=1,ivmax
            write(13,*) ir,iv,stml/normfac
         enddo
c         write(*,'("r=",f5.1," arcsec : major M/L=",f8.1,
c     &        " minor M/L=",f8.1)') RneeIR(ir)*rmax,ratml(ir,1),
c     &        ratml(ir,ivmax)
      enddo

      close(11)
      close(12)
      close(13)
      close(24)
      close(25)

      open(unit=11,file='gden.norm',status='unknown')
      write(11,*) normfac
      write(11,*) stml/normfac
      close(11)

      call halowrite()
      

      end

      function test(theta)
      real theta
      test = 1.e7*cos(theta)
      return
      end

      function func1(x)
      external func1a
      common/cfunc1/thlo,thup
      common/cff/rp
      rp=x
      call qromb2(func1a,thlo,thup,s)
      func1=s*x*x
      return
      end

      function func3(x)
      external func3a
      common/cfunc1/thlo,thup
      common/cff/rp
      rp=x
      call qromb2(func3a,thlo,thup,s)
      func3=s*x*x
      return
      end

      function func4(x)
      external func4a
      common/cfunc1/thlo,thup
      common/cff/rp
      rp=x
      call qromb2(func4a,thlo,thup,s)
      func4=s*x*x
      return
      end

      function func2(x)
      external func2a
      common/cfunc1/thlo,thup
      common/cff/rp
      rp=x
      call qromb2(func2a,thlo,thup,s)
      func2=s*x
      return
      end

      function func1a(th)
      parameter(nmax=5000)
      real rin(nmax),tin(nmax),rhoin(nmax,nmax),s(nmax,nmax)
      real mass(nmax,nmax),halomass(nmax,nmax)
      real y1,y2,y3,y4,t,u
      common/cfunci/axra,axrt
      common/cfunc1a/rin,tin,rhoin,s,rmax,mr,na,mass,halomass
      common/cff/rp
      v=sin(th)
      cv=sqrt(1.-v*v)
      x=rp*cv
      y=rp*v

      rp2 = rp*rmax*(sqrt(cv**2+(v/axrt)**2))
      
      rmd=log10(rp2)

      DO i=2,mr
        IF (rmd.lt.rin(i)) goto 3000
      ENDDO
 3000 continue
      IF (i.gt.mr) i=mr

      DO j=2,na
        IF (th.lt.tin(j)) goto 3001
      END DO
 3001 continue
      IF (j.gt.na) j=na

      r1=rin(i-1) !radius
      r2=rin(i)
      th1=tin(j-1)       !theta
      th2=tin(j)

      r01=rmd-r1
      r21=r2-r1
      th21=th2-th1
      th01=th-th1

      d11=rhoin(i-1,j-1)        !luminosity
      d21=rhoin(i,j-1)
      d12=rhoin(i-1,j)
      d22=rhoin(i,j)

      f = d11+(d21-d11)/r21*r01+(d12-d11)/th21*th01+
     &  ((d22-d12-d21+d11)/(th21*r21))*(th01*r01)

      func1a=cv*10**f

      return
      end      

      function func2a(th)
      parameter(nmax=5000)
      real rin(nmax),tin(nmax),rhoin(nmax,nmax),s(nmax,nmax)
      real mass(nmax,nmax),halomass(nmax,nmax),sbin(nmax,nmax)
      common/cfunci/axra,axrt
      common/cfunc1a/rin,tin,rhoin,s,rmax,mr,na,mass,halomass
      common/cfunc2a/sbin
      common/cff/rp
      v=sin(th)
      x=rp*sqrt(1.-v*v)
      y=rp*v

      rmd=log10(rp*rmax)

      do i=2,mr
        if (rmd.lt.rin(i)) goto 3100
      end do
 3100 continue
      if(i.gt.mr) i=mr

      do j=2,na
        if(th.lt.tin(j)) goto 3101
      end do
 3101 continue
      if(j.gt.na) j=na

      r1=rin(i-1) !radius
      r2=rin(i)
      th1=tin(j-1)       !theta
      th2=tin(j)
      r01=rmd-r1
      r21=r2-r1
      th21=th2-th1
      th01=th-th1

      s11=sbin(i-1,j-1)        !"surface brightness"=luminosity
      s21=sbin(i,j-1)
      s12=sbin(i-1,j)
      s22=sbin(i,j)


      sp = s11+(s21-s11)/r21*r01+(s12-s11)/th21*th01+
     &  ((s22-s12-s21+s11)/(th21*r21))*(th01*r01)


      func2a=10**sp

      return
      end    


      function func3a(th)
      parameter(nmax=5000)
      real rin(nmax),tin(nmax),s(nmax,nmax),rhoin(nmax,nmax)
      real mass(nmax,nmax),halomass(nmax,nmax)
      real y1,y2,y3,y4,t,u
      common/cfunci/axra,axrt
      common/cfunc1a/rin,tin,rhoin,s,rmax,mr,na,mass,halomass
      common/cff/rp
      v=sin(th)
      cv=sqrt(1.-v*v)
      x=rp*cv
      y=rp*v

      rp2 = rp*rmax*(sqrt(cv**2+(v/axrt)**2))
      
      rmd=log10(rp2)

      DO i=2,mr
        IF (rmd.lt.rin(i)) goto 3000
      ENDDO
 3000 continue
      IF (i.gt.mr) i=mr

      DO j=2,na
        IF (th.lt.tin(j)) goto 3001
      END DO
 3001 continue
      IF (j.gt.na) j=na

      r1=rin(i-1) !radius
      r2=rin(i)
      th1=tin(j-1)       !theta
      th2=tin(j)

      r01=rmd-r1
      r21=r2-r1
      th21=th2-th1
      th01=th-th1

      d11=mass(i-1,j-1)        !luminosity
      d21=mass(i,j-1)
      d12=mass(i-1,j)
      d22=mass(i,j)

      f = d11+(d21-d11)/r21*r01+(d12-d11)/th21*th01+
     &  ((d22-d12-d21+d11)/(th21*r21))*(th01*r01)

      func3a=cv*10**f

      return
      end


      function func4a(th)
      parameter(nmax=5000)
      real rin(nmax),tin(nmax),s(nmax,nmax),rhoin(nmax,nmax)
      real halomass(nmax,nmax),mass(nmax,nmax)
      real y1,y2,y3,y4,t,u
      common/cfunci/axra,axrt
      common/cfunc1a/rin,tin,rhoin,s,rmax,mr,na,mass,halomass
      common/cff/rp
      v=sin(th)
      cv=sqrt(1.-v*v)
      x=rp*cv
      y=rp*v

      rp2 = rp*rmax*(sqrt(cv**2+(v/axrt)**2))
      
      rmd=log10(rp2)

      DO i=2,mr
        IF (rmd.lt.rin(i)) goto 3000
      ENDDO
 3000 continue
      IF (i.gt.mr) i=mr

      DO j=2,na
        IF (th.lt.tin(j)) goto 3001
      END DO
 3001 continue
      IF (j.gt.na) j=na

      r1=rin(i-1) !radius
      r2=rin(i)
      th1=tin(j-1)       !theta
      th2=tin(j)

      r01=rmd-r1
      r21=r2-r1
      th21=th2-th1
      th01=th-th1

      d11=halomass(i-1,j-1)        !luminosity
      d21=halomass(i,j-1)
      d12=halomass(i-1,j)
      d22=halomass(i,j)

      f = d11+(d21-d11)/r21*r01+(d12-d11)/th21*th01+
     &  ((d22-d12-d21+d11)/(th21*r21))*(th01*r01)

      func4a=cv*10**f

      return
      end


      FUNCTION IRneeR(r)
      include 'bothdefs.h'
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IRneeR = int( log(a*r/b+cval)/a + 0.5 )
      RETURN
      END
      FUNCTION RneeIR(ir)
      include 'bothdefs.h'
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      RneeIR = b/a * (exp(a*float(ir))-cval)
      RETURN
      END
      FUNCTION IVneeV(v)
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IVneeV = int(sign(1.,v))*min(nint(float(ivmax)*abs(v)+0.5),ivmax)
      RETURN
      END
      FUNCTION VneeIV(iv)
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      VneeIV = (float(iv-1)+0.5)/float(ivmax)
      RETURN
      END
      FUNCTION IRCneeR(r)
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IRCneeR = (IRneeR(r)-1)/irrat + 1
      RETURN
      END
      FUNCTION IVCneeV(v)
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      IVCneeV = (IVneeV(v)-1)/ivrat + 1
      RETURN
      END

      FUNCTION RneeIRapo(ir)
      include 'libdefs.h'
      RneeIRapo = belz/aelz * (exp(aelz*float(ir))-celz)
      RETURN
      END

      FUNCTION IRneeRapo(r)
      include 'libdefs.h'
      IRneeRapo = int( log(aelz*r/belz+celz)/aelz + 0.5 )
      RETURN
      END

      function func(x)
      include 'bothdefs.h'
      COMMON/rbin/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat
      func=rmin*exp(x*irmax)-rmin*cval-exp(x)+cval
      return
      end

      function funcelz(x)
      include 'libdefs.h'
      funcelz=relzmin*exp(x*Nrelz)-relzmin*celz-
     &     relzmax*exp(x)+relzmax*celz
      return
      end

      function funci(x)
      common/cfunci/axra,axrt
      funci=cos(x)**2+axrt*axrt*sin(x)**2-axra*axra
      return
      end

      subroutine getbinp(irmax,ivmax,rmin,angrad,radl,radu,thl,thu,ielz)
      real radl(irmax),radu(irmax),thl(ivmax),thu(ivmax)
      integer ielz
      PARAMETER (ndiv=200000,nsample=1000)
      parameter (pitod=57.29578)

      open (unit=20,file='bindemo_r.out',status='unknown')
      open (unit=21,file='bindemo_v.out',status='unknown')
      open (unit=22,file='bin_r.out',status='unknown')
      open (unit=23,file='bin_v.out',status='unknown')
 
      write (20,*)'              bin                  bin'
      write (21,*)'              bin                  bin'
      write (20,*)'             lower      bin       upper'
      write (21,*)'             lower      bin       upper'
      write (20,*)' ir   irc     edge     center      edge'
      write (21,*)' iv   ivc     edge     center      edge'
      write (22,*)' irc    rlow      rmid      rup    bin_size'
      write (23,*)' irv    vlow      vmid      vup    bin_size'
      write (20,*)'---   ---   --------   --------   -------'
      write (21,*)'---   ---   --------   --------   -------'
      write (22,*)'----   ------    ------    ------    ------'
      write (23,*)'----   ------    ------    ------    ------'
 100  format (i3,3x,i3,3(3x,f8.6),2(3x,f8.4))
 101  format (1x,i3,4(2x,f10.3))
      rledge=0.
      vledge=0.
      do i=0,ndiv+ndiv/10
        rlo = float(i)/float(ndiv)
        vlo = float(i)/float(ndiv)
        rup = float(i+1)/float(ndiv)
        vup = float(i+1)/float(ndiv)
        irlo = IRneeR(rlo)
        ivlo = IVneeV(vlo)
        irup = IRneeR(rup)
        ivup = IVneeV(vup)
        if (irlo.ne.irup.and.irlo.ge.0.and.irup.ge.0) then
          rmid = RneeIR(irlo)
          redge = .5 * (rlo+rup)
          write (20,100)irlo,IRCneeR(rmid),rledge,rmid,redge,
     $         rmid*angrad,redge*angrad
          if(irlo.gt.0.and.irlo.le.irmax) then
             radl(irlo)=rledge
             radu(irlo)=redge
          endif
          rledge=redge
        endif
        if (ivlo.ne.ivup) then
          vmid = VneeIV(ivlo)
          vedge = .5 * (vlo+vup)
          write (21,100)ivlo,IVCneeV(vmid),vledge,vmid,vedge
          if(ivlo.gt.0.and.ivlo.le.ivmax) then
             thl(ivlo)=asin(vledge)
             thu(ivlo)=asin(vedge)
          endif
          vledge=vedge
        endif
      enddo
      vlo=vedge
      vedge=1.0
      vmid=(vedge+vlo)/2.
      write (21,100)ivup,IVCneeV(vmid),vledge,vmid,vedge
      thl(ivup)=asin(vledge)
      thu(ivup)=asin(vedge)

      rloold=rmin/2.
      vloold=0.
      do i=0,ndiv+ndiv/10
        rlo = float(i)/float(ndiv)
        vlo = float(i)/float(ndiv)
        rup = float(i+1)/float(ndiv)
        vup = float(i+1)/float(ndiv)
        irlo = IRCneeR(rlo)
        ivlo = IVCneeV(vlo)
        irup = IRCneeR(rup)
        ivup = IVCneeV(vup)
        if (irlo.ne.irup.and.irlo.ge.0.and.irup.ge.0) then
          redge = .5 * (rlo+rup)
          rmid=(redge+rloold)/2.
          write (22,101)IRCneeR(rmid),rloold*angrad,rmid*angrad,
     $         redge*angrad,(redge-rloold)*angrad
          rloold=redge
        endif
        if (ivlo.ne.ivup) then
          vedge = .5 * (vlo+vup)
          vmid=(vedge+vloold)/2.
          thlo=pitod*asin(vloold)
          thup=pitod*asin(vedge)
          thmid=(thlo+thup)/2.
          write (23,101)IVCneeV(vmid),thlo,thmid,thup,thup-thlo
          vloold=vedge
        endif
      enddo
      vedge=1.0
      vmid=(vedge+vloold)/2.
      thlo=pitod*asin(vloold)
      thup=pitod*asin(vedge)
      thmid=(thlo+thup)/2.
      write (23,101)IVCneeV(vmid),thlo,thmid,thup,thup-thlo

      close(20)
      close(21)
      close(22)
      close(23)

      open(unit=24,file='binelz.out',status='unknown')
      write(24,'("  i  r/angrad   r [arcsec]")')
      do i=1,ielz
 104     format(i3,2x,f8.5,3x,f8.3)
         write(24,104) i,RneeIRapo(i),RneeIRapo(i)*angrad
      end do
      close(24)

      open (unit=30,file='binsample.out',status='unknown')
      do i = 1,nsample
        r = float(i)/float(nsample)
        v = float(i)/float(nsample)
        ir = IRneeR(r)
        iv = IVneeV(v)
	if(ir.ge.0) then
        write (30,*)i,ir,iv,r,v
	end if
      enddo
      close(30)

      return
      END
 
      FUNCTION rtbiso(func,x1,x2,xacc)
      INTEGER JMAX
      REAL rtbiso,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) print *,'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbiso=x1
        dx=x2-x1
      else
        rtbiso=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbiso+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbiso=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      print *,'too many bisections in rtbis'
      END

      SUBROUTINE qromb2(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd2(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint2(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      print *,'too many steps in qromb2'
      END
      SUBROUTINE polint2(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)print *,'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
      SUBROUTINE trapzd2(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      print *,'too many steps in qromb'
      END
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)print *,'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
