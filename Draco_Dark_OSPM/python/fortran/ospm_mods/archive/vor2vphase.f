c ********************************************************************
c * CELLS2VPHASE: transforms the output of the voronoi-tesselation   * 
c *               into phasevols of the orbits                       *
c ********************************************************************
      program vor2vphase
      INCLUDE 'libdefs.h'
      parameter(nmax=Norbit,nlim=2*Nrdat)
      real xmin,xmax,ymin,ymax
      real sos_sum(nmax),dEdLz(nmax),E(nmax),Lz(nmax)
      real phasevol(nmax)
      real area,I_3(nmax),sum
      real E_fac,Lz_fac,t_sos,E_last,Lz_last
      real diff(nmax),xedge(5),yedge(5)
      real E_arr(nlim,nlim),Lz_arr(nlim,nlim)
      real d_E,d_Lz
      real ra,rp,phi_ra,phi_rp
      integer iedge(5),iLz(nmax),iE(nmax),nseq
      integer iE_max(nlim),iE_min(nlim)
      integer iEstep_up(nlim,nlim),iEstep_lo(nlim,nlim)
      integer iLzstep_up(nlim,nlim),iLzstep_lo(nlim,nlim)
      character lf*4

      real dEdLz_seq(nmax),E_seq(nmax),Lz_seq(nmax)
      integer iLz_seq(nmax),iE_seq(nmax),bad_edges,n_seq

      integer i,j,norb,iorb,i1,i2,ix,iEmax,iEmin
      integer iLzmax,iLzmin,iE_last,iLz_last,i_e,i_lz
      integer iE_step_min,iLz_step_min,iLz_step,iE_step
      integer iorblast,iseq,il_act,ie_act,imod
      real x1,x2,x3,x4,xarea,xjunk,xnorm,ynorm

      integer irperi

      real hits

c --- get library
      call totmlread()
      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call dataread()
      call galaxyread()
      call haloread()
      call spacel()
      cden = dL(1,1)*ratML(1,1) / vol2d(1)
      core = 4.*pi*cden*rmin**3/(3.+cslope)
      call density()
      call tablesread()
      if(iquad.ne.0) call forcefix()
      fac1 = sqrt(GG*totlight*arcsec/distance/angrad)
      distance = distance/1.e3
c --- end

c --- initialization
      do i=1,nlim
         do j=1,nlim
            E_arr(i,j)=1.e15
            Lz_arr(i,j)=1.e15
         end do
      end do

      open(unit=8,file='nseq.out',status='old')
      read(8,*) nseq
      close(8)

      open(unit=68,file='norbit',status='old')
      read(68,*) Norb
      close(68)

      open(unit=11,file='integrals.out',status='old')
      do i=1,Norb
         read(11,*) x1,x2,x3,x4
         E(i)=log(-x2*fac1*fac1)
         Lz(i)=log(x1*fac1*angrad*distance/arcsec)
      end do
      close(11)

      do i=1,Norb
         sos_sum(i)=0.
      end do
c --- end



c ================================
c =                              =
c =       calculate d(E,Lz)      =
c =                              =
c ================================


c --- read orbital E- and Lz-bins
      xmax=-1.e10
      ymax=-1.e10
      xmin=1.e10
      ymin=1.e10
      open(unit=8,file='sos2.out',status='old')

 800  continue
      read(8,*,end=778) iorb,i1,i2,x1,x2
      iLz(iorb)=i2
      iE(iorb)=i1
      xmax=max(xmax,float(i1))
      xmin=min(xmin,float(i1))
      ymax=max(ymax,float(i2))
      ymin=min(ymin,float(i2))
      goto 800

 778  continue
      close(8)
      iEmax=nint(xmax)
      iEmin=nint(xmin)
      iLzmax=nint(ymax)
      iLzmin=nint(ymin)
c --- end 


c --- convert orbit energies to seequence energies
      iLz_last=-1
      iE_last=-1
      n_seq=0
      do i=1,Norb
         if(iLz(i).ne.iLz_last.or.iE(i).ne.iE_last) then
            n_seq=n_seq+1
            iLz_seq(n_seq)=iLz(i)
            iE_seq(n_seq)=iE(i)
            E_seq(n_seq)=E(i)
            Lz_seq(n_seq)=Lz(i)
            E_arr(iE_seq(n_seq),iLz_seq(n_seq))=E(i)
            Lz_arr(iE_seq(n_seq),iLz_seq(n_seq))=Lz(i)
         end if
         iLz_last=iLz(i)
         iE_last=iE(i)
      end do
c --- end


c --- determine iE_step and iLz_step
c      iE_step=10000
c      iLz_step=10000
c      do i=1,n_seq
c         if(iE_seq(i).ne.iEmin) then
c            if(iE_seq(i)-iEmin.lt.iE_step) then
c               iE_step=iE_seq(i)-iEmin
c            end if
c         end if
c         if(iLz_seq(i).ne.iLzmin) then
c            if(iLz_seq(i)-iLzmin.lt.iLz_step) then
c               iLz_step=iLz_seq(i)-iLzmin
c            end if
c         end if
c      end do

      iE_step = iEstep
      iLz_step = iLzstep

c --- end


c --- determine local iE_max and iE_min
      do i=iLzmin,iLzmax
         xmax=-1000
         xmin=1000
         do j=1,n_seq
            if(iLz_seq(j).eq.i) then
               xmax=max(float(iE_seq(j)),xmax)
               xmin=min(float(iE_seq(j)),xmin)
            end if
         end do
         iE_max(i)=nint(xmax)
         iE_min(i)=nint(xmin)
      end do
c --- end


c      call pgbegin(0,'/cps',1,1)
c      call pgsch(1.2)
c      call pgscf(2)
c      call gminmax(n_seq,E_seq,xmin,xmax)
c      call gminmax(n_seq,Lz_seq,ymin,ymax)
c      call pgsch(1.5)
c      call pgenv(xmin,xmax,ymin,ymax,0,0)
c      call pglabel('ln |E| [km\\U2\\D/s\\U2\\D]',
c     &     'ln L\\Dz\\U [km/s kpc]','')
c      call pgsch(1.2)
c      call pgpoint(n_seq,E_seq,Lz_seq,17)


c --- calculate d(E,Lz) for sequences from displaced grid
      do i=1,n_seq
         xedge(5)=E_arr(iE_seq(i),iLz_seq(i))
         yedge(5)=Lz_arr(iE_seq(i),iLz_seq(i))
         ! upper left edge
         i_lz=iLz_seq(i)+iLz_step
         i_e=iE_seq(i)+iE_step
         rp=(RneeIRperi(iLz_seq(i))+RneeIRperi(i_lz))/2.
         ra=(RneeIRapo(iE_seq(i))+RneeIRapo(i_e))/2.
         if(abs(rp-ra)/ra.lt.0.0001) then
            irperi = IRneeRperi(rp)
            rp=RneeIRperi(irperi)
            ra=.5*(RneeIRperi(irperi)+RneeIRperi(irperi+1))
         end if

         if(rp.eq.0.0) then
            print*,'rp = 0',rp,ra
            rp = ra/10.
         end if

         call potential(rp,0.,phi_rp)
         call potential(ra,0.,phi_ra)
         x1=sqrt(2.*abs(phi_ra-phi_rp)/abs(ra**2-rp**2))*ra*rp
         x2=phi_rp+0.5*x1*x1/rp/rp
         yedge(4)=log(x1*fac1*angrad*distance/arcsec)
         xedge(4)=log(abs(x2*fac1*fac1))
         iedge(4)=1

         ! lower left edge
         i_lz=iLz_seq(i)-iLz_step
         i_e=iE_seq(i)+iE_step
         if(i_lz.lt.1) then
            rp=(RneeIRperi(iLz_seq(i))
     &           +RneeIRperi(iLz_seq(i)+iLzstep))/2.
            rp=2*log(RneeIRperi(iLz_seq(i)))-log(rp)
            rp=exp(rp)
            if(rp.lt.0.0) then
               rp=.75*RneeIRperi(iLz_seq(i))
            end if
         else
            rp=(RneeIRperi(iLz_seq(i))+RneeIRperi(i_lz))/2.
         end if
         ra=(RneeIRapo(iE_seq(i))+RneeIRapo(i_e))/2.
         if(abs(rp-ra)/ra.lt.0.0001) then
            irperi = IRneeRperi(rp)
            rp=RneeIRperi(irperi)
            ra=.5*(RneeIRperi(irperi)+RneeIRperi(irperi+1))
         end if

         if(rp.eq.0.0) then
            print*,'rp = 0',rp,ra
            rp = ra/10.
         end if

         call potential(rp,0.,phi_rp)
         call potential(ra,0.,phi_ra)
         x1=sqrt(2.*abs(phi_ra-phi_rp)/abs(ra**2-rp**2))*ra*rp
         x2=phi_rp+0.5*x1*x1/rp/rp
         yedge(3)=log(x1*fac1*angrad*distance/arcsec)
         xedge(3)=log(abs(x2*fac1*fac1))
         iedge(3)=1

         ! lower right edge
         i_lz=iLz_seq(i)-iLz_step
         i_e=iE_seq(i)-iE_step
         if(i_lz.lt.1) then
            rp=(RneeIRperi(iLz_seq(i))
     &           +RneeIRperi(iLz_seq(i)+iLzstep))/2.
            rp=2*log(RneeIRperi(iLz_seq(i)))-log(rp)
            rp=exp(rp)
            if(rp.lt.0.0) then
               rp=.75*RneeIRperi(iLz_seq(i))
            end if
         else
            rp=(RneeIRperi(iLz_seq(i))+RneeIRperi(i_lz))/2.
         end if
         if(i_e.lt.1) then
            ra=(RneeIRapo(iE_seq(i))
     &           +RneeIRapo(iE_seq(i)+iEstep))/2.
            ra=2*log(RneeIRapo(iE_seq(i)))-log(ra)
            ra=exp(ra)
            if(rp.lt.0.0) then
               ra=.75*RneeIRapo(iE_seq(i))
            end if
         else
            ra=(RneeIRapo(iE_seq(i))+RneeIRapo(i_e))/2.  
         end if
         if(abs(rp-ra)/ra.lt.0.0001) then
            irperi = IRneeRperi(rp)
            rp=RneeIRperi(irperi)
            ra=.5*(RneeIRperi(irperi)+RneeIRperi(irperi+1))
         end if

         if(rp.eq.0.0) then
            print*,'rp = 0',rp,ra
            rp = ra/10.
         end if

         call potential(rp,0.,phi_rp)
         call potential(ra,0.,phi_ra)
         x1=sqrt(2.*abs(phi_ra-phi_rp)/abs(ra**2-rp**2))*ra*rp
         x2=phi_rp+0.5*x1*x1/rp/rp
         yedge(2)=log(x1*fac1*angrad*distance/arcsec)
         xedge(2)=log(abs(x2*fac1*fac1))
         iedge(2)=1

         ! upper right edge
         if(iLz_seq(i).eq.iE_seq(i)) then
            xedge(1)=xedge(5)
            yedge(1)=yedge(5)
            iedge(1)=1
         else
            i_lz=iLz_seq(i)+iLz_step
            i_e=iE_seq(i)-iE_step
            rp=(RneeIRperi(iLz_seq(i))+RneeIRperi(i_lz))/2.
            if(i_e.lt.1) then
               ra=(RneeIRapo(iE_seq(i))
     &              +RneeIRapo(iE_seq(i)+iEstep))/2.
               ra=2*log(RneeIRapo(iE_seq(i)))-log(ra)
               ra=exp(ra)
               if(rp.lt.0.0) then
                  ra=.75*RneeIRapo(iE_seq(i))
               end if
            else
               ra=(RneeIRapo(iE_seq(i))+RneeIRapo(i_e))/2.  
            end if
            if(abs(rp-ra)/ra.lt.0.0001) then
               irperi = IRneeRperi(rp)
               rp=RneeIRperi(irperi)
               ra=.5*(RneeIRperi(irperi)+RneeIRperi(irperi+1))
            end if

            if(rp.eq.0.0) then
               print*,'rp = 0',rp,ra
               rp = ra/10.
            end if

            call potential(rp,0.,phi_rp)
            call potential(ra,0.,phi_ra)
            x1=sqrt(2.*abs(phi_ra-phi_rp)/abs(ra**2-rp**2))*ra*rp
            x2=phi_rp+0.5*x1*x1/rp/rp
            yedge(1)=log(x1*fac1*angrad*distance/arcsec)
            xedge(1)=log(abs(x2*fac1*fac1))
            iedge(1)=1
         end if

c --- old version
c         call ELzarea(xedge,yedge,iedge,xarea,bad_edges)
c         dEdLz_seq(i)=xarea*exp(E_seq(i)+Lz_seq(i))
c --- end of old version

c --- new version
         call ELzarea2(xedge,yedge,iedge,xarea,bad_edges,hits)
         dEdLz_seq(i)=xarea
c --- end of new version

         write(*,'("nseq=",i4," iLz=",i4," iE=",i4," frac=",f9.6a1,$)') 
     &        n_seq-i,iLz_seq(i),iE_seq(i),hits*100.,char(13)
         if(dEdLz_seq(i).gt.0.0.and.dEdLz_seq(i).lt.1.e10) then
            ! OK
         else
            dEdLz_seq(i)=1.e-5
         end if
         call plotELzarea2(xedge,yedge,E_seq,Lz_seq,n_seq)
c         call plotELzarea3(xedge,yedge,E_seq,Lz_seq,n_seq)
      end do

c      call pgend()

c --- end      

      ! convert sequential d(E,Lz) into orbital d(E,Lz)
      do j=1,n_seq
         do i=1,Norb
            if(iLz(i).eq.iLz_seq(j).and.iE(i).eq.iE_seq(j)) then
              dEdLz(i)=dEdLz_seq(j)
            end if
         end do
      end do


c ================================
c =                              =
c =       calculate dI3          =
c =                              =
c ================================
         

c --- read and sum up the orbital Voronoi-areas weighted with the SOS-periods
      iorblast=0
      do iseq=1,nseq
         if(iseq.lt.10) write(lf,'("000",i1)') iseq
         if(iseq.ge.10.and.iseq.lt.100) write(lf,'("00",i2)') iseq
         if(iseq.ge.100.and.iseq.lt.1000) write(lf,'("0",i3)') iseq
         if(iseq.ge.1000) write(lf,'(i4)') iseq
         open(unit=25,file=lf//'.vor.cells',status='old')
         open(unit=26,file=lf//'.vor.sos',status='old')

 700     continue
         read(25,*,err=101,end=777) il_act,ie_act,i_mod,area
         goto 102
 101     i_mod=1
 102     continue
         read(26,*) x1,x2,xjunk,iorb,t_sos,xnorm,ynorm
         if(i_mod.ne.1.and.iorb.ne.0) then
            sos_sum(iorb)=sos_sum(iorb)+area*t_sos*xnorm*ynorm
         end if
         goto 700

 777     continue
         close(26)
         close(25)
      end do


c --- define and write out the I_3 integral
      open(unit=45,file='periapo.out',status='old')
      do i=1,Norb
         read(45,*) i1,i2,xperi,xapo,xmin,xmax
         I_3(i) = (xmax-xmin)/(xapo-xperi)
      end do
      close(45)
      Lz_last=10**(Lz(Norb))
      E_last=10**(E(Norb))
      do i=Norb-1,1,-1
         Lz_last=10**(Lz(i))
         E_last=10**(E(i))
      end do
      open(unit=8,file='integrals2.out',status='unknown')
      do i=1,Norb
         write(8,*) exp(Lz(i))/fac1/angrad/distance*arcsec,
     &        -exp(E(i))/fac1/fac1,I_3(i)
      end do
      close(8)
      open(unit=8,file='integrals3.out',status='unknown')
      do i=1,Norb
         write(8,*) 2*i-1,iLz(i),iE(i),I_3(i)
         write(8,*) 2*i,iLz(i),iE(i),I_3(i)
      end do
      close(8)
c --- end


c --- calculate the final phasevols
      open(unit=21,file='vor.failure',status='unknown')
      open(unit=22,file='dvphase.out',status='unknown')
      do i=1,Norb
         write(22,*) i,iLz(i),iE(i),sos_sum(i),dEdLz(i)
         if(sos_sum(i)*dEdLz(i).eq.0.0) then 
            if(sos_sum(i).eq.0.0) then
               print*,'no I3 for orbit:',i
               write(21,*) 'no I3 for orbit:',i
            end if
            if(dEdLz(i).eq.0.0) then
               print*,'no dELz for orbit:',i
               write(21,*) 'no dELz for orbit:',i
            end if
            phasevol(i)=1.e10    !set vol=tiny, if no value exists
         else
            phasevol(i)=abs(1./(sos_sum(i)*dEdLz(i)))
         end if
      end do

      close(21)
      close(22)

      ! normalization to avoid numerical problems in model.x
      sum=0.0
      do i=1,Norb
         sum=sum+1./phasevol(i)
      end do
      do i=1,Norb
         phasevol(i)=phasevol(i)*sum
      end do
c --- end


c --- write out volumes
      open(unit=9,file='phase.vor',status='unknown')
      write(9,*)(phasevol(i),i=1,Norb)
      close(9)
c --- end   

      end


      subroutine ELzarea2(xedge,yedge,iedge,xarea,bad_edges,hits)
      parameter(nx=750,ny=750)
c      parameter(nx=5,ny=5)
      real xedge(5),yedge(5),dia1(2),dia2(2)
      real w1,w2,w3,w4,w,xp,xa,ya,da,db
      integer iedge(5),bad_edges
      real dotprod,hits
      
      xp=3.1415926539

      bad_edges=4
      itemp=0
      itemp2=0
      do i=1,4
         bad_edges=bad_edges-iedge(i)
         if(iedge(i).ne.0.and.itemp.eq.0) itemp=i
         if(iedge(i).ne.0.and.itemp.ne.0) itemp2=i
      end do

      if(bad_edges.eq.3) then
         xarea=4.*(abs(xedge(5)-xedge(itemp))*
     &        abs(yedge(5)-yedge(itemp)))
      end if
      if(bad_edges.eq.2) then
         if(iedge(4).eq.0.and.iedge(3).eq.0) then
            xedge(3)=xedge(5)-abs(xedge(2)-xedge(5))
            yedge(3)=yedge(2)
            xedge(4)=xedge(5)-abs(xedge(1)-xedge(5))
            yedge(4)=yedge(1)
            bad_edges=0
         end if
         if(iedge(2).eq.0.and.iedge(3).eq.0) then
            xedge(2)=xedge(1)
            yedge(2)=yedge(5)-abs(yedge(1)-yedge(5))
            xedge(3)=xedge(4)
            yedge(3)=yedge(5)-abs(yedge(4)-yedge(5))
            bad_edges=0
         end if
         if(bad_edges.ne.0) then
            xarea=0.
            do i=1,4
               if(iedge(i).eq.1) then
                  xarea=xarea+(abs(xedge(5)-xedge(i))*
     &                 abs(yedge(5)-yedge(i)))
               end if
            end do
         end if
      end if
      if(bad_edges.eq.1) then
         do i=1,4
            if(iedge(i).eq.0) then
               xedge(i)=xedge(5)
               yedge(i)=yedge(5)
            end if
         end do
         bad_edges=0
      end if
      xmin=1.e10
      xmax=-1.e10
      ymin=1.e10
      ymax=-1.e10
      do i=1,4
         xmin=min(xmin,xedge(i))
         xmax=max(xmax,xedge(i))
         ymin=min(ymin,yedge(i))
         ymax=max(ymax,yedge(i))
      end do
      xmin=min(1.1*xmin,0.9*xmin)
      ymin=min(1.1*ymin,0.9*ymin)
      xmax=max(1.1*xmax,0.9*xmax)
      ymax=max(1.1*ymax,0.9*ymax)
      if(bad_edges.eq.0) then
         xarea=0.
         da=(xmax-xmin)/float(nx-1)
         db=(ymax-ymin)/float(ny-1)
         hits=0.
         do ix=1,nx
            xa=xmin+(xmax-xmin)/float(nx-1)*float(ix-1)
            do iy=1,ny
               ya=ymin+(ymax-ymin)/float(ny-1)*float(iy-1)
               x1=xedge(1)-xa
               y1=yedge(1)-ya
               x2=xedge(2)-xa
               y2=yedge(2)-ya
               xprod=x1*x2+y1*y2
               xabs1=sqrt(x1*x1+y1*y1)
               xabs2=sqrt(x2*x2+y2*y2)
               w1=acos(xprod/xabs1/xabs2)
               x1=xedge(2)-xa
               y1=yedge(2)-ya
               x2=xedge(3)-xa
               y2=yedge(3)-ya
               xprod=x1*x2+y1*y2
               xabs1=sqrt(x1*x1+y1*y1)
               xabs2=sqrt(x2*x2+y2*y2)
               w2=acos(xprod/xabs1/xabs2)
               x1=xedge(3)-xa
               y1=yedge(3)-ya
               x2=xedge(4)-xa
               y2=yedge(4)-ya
               xprod=x1*x2+y1*y2
               xabs1=sqrt(x1*x1+y1*y1)
               xabs2=sqrt(x2*x2+y2*y2)
               w3=acos(xprod/xabs1/xabs2)
               x1=xedge(4)-xa
               y1=yedge(4)-ya
               x2=xedge(1)-xa
               y2=yedge(1)-ya
               xprod=x1*x2+y1*y2
               xabs1=sqrt(x1*x1+y1*y1)
               xabs2=sqrt(x2*x2+y2*y2)
               w4=acos(xprod/xabs1/xabs2)
               w=w1+w2+w3+w4
               if(abs(w-2.*xp)/abs(w).lt.0.01) then
                  xarea=xarea+exp(xa+ya)*da*db
                  hits=hits+1.
               end if
            end do
         end do
      end if
      hits=hits/float(nx)/float(ny)
      return
      end


      subroutine ELzarea(xedge,yedge,iedge,xarea,bad_edges)
      real xedge(5),yedge(5),dia1(2),dia2(2)
      integer iedge(5),bad_edges
      real dotprod
      bad_edges=4
      itemp=0
      itemp2=0
      do i=1,4
         bad_edges=bad_edges-iedge(i)
         if(iedge(i).ne.0.and.itemp.eq.0) itemp=i
         if(iedge(i).ne.0.and.itemp.ne.0) itemp2=i
      end do

      if(bad_edges.eq.3) then
         xarea=4.*(abs(xedge(5)-xedge(itemp))*
     &        abs(yedge(5)-yedge(itemp)))
      end if
      if(bad_edges.eq.2) then
         if(iedge(4).eq.0.and.iedge(3).eq.0) then
            xedge(3)=xedge(5)-abs(xedge(2)-xedge(5))
            yedge(3)=yedge(2)
            xedge(4)=xedge(5)-abs(xedge(1)-xedge(5))
            yedge(4)=yedge(1)
            bad_edges=0
         end if
         if(iedge(2).eq.0.and.iedge(3).eq.0) then
            xedge(2)=xedge(1)
            yedge(2)=yedge(5)-abs(yedge(1)-yedge(5))
            xedge(3)=xedge(4)
            yedge(3)=yedge(5)-abs(yedge(4)-yedge(5))
            bad_edges=0
         end if
         if(bad_edges.ne.0) then
            xarea=0.
            do i=1,4
               if(iedge(i).eq.1) then
                  xarea=xarea+(abs(xedge(5)-xedge(i))*
     &                 abs(yedge(5)-yedge(i)))
               end if
            end do
         end if
      end if
      if(bad_edges.eq.1) then
         do i=1,4
            if(iedge(i).eq.0) then
               xedge(i)=xedge(5)
               yedge(i)=yedge(5)
            end if
         end do
         bad_edges=0
      end if
      if(bad_edges.eq.0) then
         dia1(1)=xedge(1)-xedge(3)
         dia1(2)=yedge(1)-yedge(3)
         dia2(1)=xedge(2)-xedge(4)
         dia2(2)=yedge(2)-yedge(4)
         xdia1=sqrt(dia1(1)**2+dia1(2)**2)
         xdia2=sqrt(dia2(1)**2+dia2(2)**2)
         dotprod=0.
         do i=1,2
            dotprod=dotprod+dia1(i)*dia2(i)
         end do
         xsinalpha=sqrt(1.-(dotprod/xdia1/xdia2)**2)
         xarea=.5*xdia1*xdia2*xsinalpha
      end if
      return
      end


      subroutine plotELzarea(xedge,yedge)
      real xedge(5),yedge(5),dia1(2),dia2(2)
      integer iedge(5),icol
c      call pgpoint(1,xedge(5),yedge(5),17)
      xedge(5)=xedge(1)
      yedge(5)=yedge(1)
c      call pgline(5,xedge,yedge)
      return
      end

      subroutine plotELzarea2(xedge,yedge,E_seq,Lz_seq,n_seq)
      parameter(nmax=10000)
      real xedge(5),yedge(5),dia1(2),dia2(2)
      real E_seq(nmax),Lz_seq(nmax)
      integer iedge(5),icol,n_seq
c      call gminmax(n_seq,E_seq,xmin,xmax)
c      call gminmax(n_seq,Lz_seq,ymin,ymax)
      xedge(5)=xedge(1)
      yedge(5)=yedge(1)
c      call pgsci(2)
c      call pgline(5,xedge,yedge)
c      call pgsci(1)
      return
      end

      subroutine plotELzarea3(xedge,yedge,E_seq,Lz_seq,n_seq)
      parameter(nmax=10000)
      real xedge(5),yedge(5),dia1(2),dia2(2)
      real E_seq(nmax),Lz_seq(nmax)
      integer iedge(5),icol,n_seq
c      call gminmax(n_seq,E_seq,xmin,xmax)
c      call gminmax(n_seq,Lz_seq,ymin,ymax)
c      call pgenv(xmin,xmax,ymin,ymax,0,0)
c      call pgpoint(1,xedge(5),yedge(5),17)
      xedge(5)=xedge(1)
      yedge(5)=yedge(1)
c      call pgsci(2)
c      call pgline(5,xedge,yedge)
c      call pgsci(1)
      return
      end
