      program plotlibsos
      INCLUDE 'libdefs.h'
c      parameter(nmax=100000)
      parameter(nmax=norbsmax*Nsos)
      real lz(2*Norbit),erg(2*Norbit)
      real xmin,xmax,ymin,ymax,fac1
      real x(nmax),y(nmax),d_E,d_Lz
      real rs(Norbit,Nsos),vs(Norbit,Nsos)
      integer ebin(Norbit),lzbin(Norbit),nsosorb(Norbit)
      integer e_seq(2*Norbit),lz_seq(2*Norbit),nseq
      character cinc*40,clabel(200)*3
      integer nc,iml,Norb2

      do i=1,200
         write(clabel(i),'(i3)') i
      end do

      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call totmlread()
      call galaxyread()

      fac1 = sqrt(GG*(totlight+hole)*arcsec/distance/angrad)
      distance = distance/1.e3

      open(unit=68,file='norbit',status='old')
      read(68,*) Norb
      close(68)
      Norb2=2*Norb

      open(unit=11,file='integrals.out',status='old')
      do i=1,Norb
         read(11,*) x1,x2,x3,x4
         erg(2*i-1)=log10(-x2*fac1*fac1)
         erg(2*i)=log10(-x2*fac1*fac1)
         lz(2*i-1)=log10(x1*fac1*angrad*distance/arcsec)
         lz(2*i)=log10(x1*fac1*angrad*distance/arcsec)
      end do
      close(11)

      call pgbegin(0,'/cps',3,3)
      call pgsch(1.2)
      call pgscf(2)

c - plot orbitsampling in (E-Lz)-plane
      do i=1,Norb
         lz(2*i-1)=-10**lz(2*i-1)
         lz(2*i)=10**lz(2*i)
      end do
      call gminmax(Norb2,erg,xmin,xmax)
      call gminmax(Norb2,lz,ymin,ymax)
      call pgenv(xmin,xmax,ymin,ymax,0,10)
      call pgpoint(Norb2,erg,lz,17)
      call pglabel('log(-E [km\\U2\\D/s\\U2\\D])',
     &     'L\\Dz\\U [km/s kpc]','')
c - end

c - read SOSs
      open(unit=43,file='sos1.out',status='old')
      open(unit=44,file='sos2.out',status='old')
      isos=0
      nseq=0
      iorblast=1
      do i=1,90000
         do j=1,90000
            read(44,*,end=777) iorb,iEtemp,iLztemp,d_Lz,d_E
            read(43,*) x1,x2,x3
            if(iorb.ne.iorblast) then
               isos=0
            end if
            isos=isos+1
            rs(iorb,isos)=log10(x1*angrad*distance/arcsec)
            vs(iorb,isos)=x2*fac1
            nsosorb(iorb)=isos
            ebin(iorb)=iEtemp
            lzbin(iorb)=iLztemp
            iorblast=iorb
         end do
      end do
 777  continue
      close(43)
      close(44)
c --- end

c --- plot sos
      open(unit=44,file='sos2.out',status='old')
      read(44,*) iorb,iEtemp,iLztemp,d_Lz,d_E,i3
      ilzlast=iLztemp
      ielast=iEtemp
      close(44)
      open(unit=40,file='integrals.out',status='old')
      isos=0
      do i=1,Norb
         if(ebin(i).eq.ielast.and.lzbin(i).eq.ilzlast) then
            read(40,*) xlz,xe,x3,x4
            do j=1,nsosorb(i)
               isos=isos+1
               x(isos)=rs(i,j)
               y(isos)=vs(i,j)
            end do
         else
            nseq=nseq+1
            e_seq(nseq)=ielast
            lz_seq(nseq)=ilzlast
            call gminmax(isos,x,xmin,xmax)
            call gminmax(isos,y,ymin,ymax)
            call pgenv(xmin,xmax,ymin,ymax,0,10)
            iml=nint(abs(xe*100.))
            call pgnumb(iml,-2,0,cinc,nc)
            call puttext(0.95,0.95,1.,'E = -'//cinc(1:nc),
     &           xmin,xmax,ymin,ymax)
            iml=nint(abs(xlz*10000.))
            call pgnumb(iml,-4,0,cinc,nc)
            call puttext(0.95,0.875,1.,'L\\Dz\\U = '//cinc(1:nc),
     &           xmin,xmax,ymin,ymax)
            icol=0
            do i2=1,Norb
               if(ebin(i2).eq.ielast.and.lzbin(i2).eq.ilzlast) then
                  do j=1,nsosorb(i2)
                     x(j)=rs(i2,j)
                     y(j)=vs(i2,j)
                  end do
                  if(icol.eq.15) icol=0
                  icol=icol+1
                  call pgsci(icol)
                  call pgpoint(nsosorb(i2),x,y,17)
                  call pgsci(1)
               end if
            end do
            call pglabel('r [kpc] -- iLz='//clabel(lz_seq(nseq))//
     &           ' iE='//clabel(e_seq(nseq)),'v\\Dr\\U [km/s]','')
            isos=0   
            read(40,*) xlz,xe
            ilzlast=lzbin(i)
            ielast=ebin(i)
            do j=1,nsosorb(i)
               isos=isos+1
               x(isos)=rs(i,j)
               y(isos)=vs(i,j)
            end do
         end if
      end do
      nseq=nseq+1
      e_seq(nseq)=ielast
      lz_seq(nseq)=ilzlast
      call gminmax(isos,x,xmin,xmax)
      call gminmax(isos,y,ymin,ymax)
      call pgenv(xmin,xmax,ymin,ymax,0,10)
      iml=nint(abs(xe*100.))
      call pgnumb(iml,-2,0,cinc,nc)
      call puttext(0.95,0.95,1.,'E = -'//cinc(1:nc),
     &     xmin,xmax,ymin,ymax)
      iml=nint(abs(xlz*10000.))
      call pgnumb(iml,-4,0,cinc,nc)
      call puttext(0.95,0.875,1.,'L\\Dz\\U = '//cinc(1:nc),
     &     xmin,xmax,ymin,ymax)
      call pgpoint(isos,x,y,17)
      call pglabel('r [kpc] -- iLz='//clabel(lz_seq(nseq))//
     &     ' iE='//clabel(e_seq(nseq)),'v\\Dr\\U [km/s]','')
      isos=0   
      do j=1,nsosorb(i)
         isos=isos+1
         x(isos)=rs(i,j)
         y(isos)=vs(i,j)
      end do
      close(40)
c --- end

      call pgend

      end
