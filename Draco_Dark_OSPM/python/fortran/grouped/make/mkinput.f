
      parameter(nmax=100000)
      real ra(nmax),dec(nmax)
      character cdate*8,cshot*3,file1*120,cifu*3,a3*13,a1*3
      character a11*11,a12*12

      open(unit=1,file="input",status='old')
      read(1,*) cdate,cshot
      read(1,*) nsim
      read(1,*) cifu
      read(1,*) idum
      read(1,*) wavel,waveh
      read(1,*) fluxl,fluxh
      read(1,*) xdiff
      read(1,*) wsig
      close(1)

      xdiff=xdiff/3600.

      call getseed(idum)

c      file1="/data/00115/gebhardt/detect/"//cdate//"v"//cshot//
c      file1="/scratch/00115/gebhardt/detect/"//cdate//"v"//cshot//
c     $     "/dithall.use"
      file1="info"

      open(unit=1,file=file1,status='old')
      n=0
      rmin=1e10
      rmax=-rmin
      dmin=rmin
      dmax=rmax
      do i=1,nmax
         read(1,*,end=666) i1,x1,x2,a3
         if(a3(11:13).eq.cifu) then
            n=n+1
            ra(n)=x1
            dec(n)=x2
            rmin=min(rmin,ra(n))
            rmax=max(rmax,ra(n))
            dmin=min(dmin,dec(n))
            dmax=max(dmax,dec(n))
         endif
      enddo
 666  continue      
      close(1)

      cosd=cos((dmax+dmin)/2./57.29578)
      rmax=rmax+2./3600/cosd
      rmin=rmin-2./3600/cosd
      dmax=dmax+2./3600.
      dmin=dmin-2./3600.
      rr=rmax-rmin
      dr=dmax-dmin
      diffmin=1.5

      wdiff=waveh-wavel
      fdiff=fluxh-fluxl
      cdec=cos(dec(1)/57.3)
      open(unit=11,file='out',status='unknown')
      do i=1,nsim
c         iran=nint(ran2(idum)*(n-1))+1
c         ra1=ra(iran)
c         dec1=dec(iran)
c         roff=-xdiff+2.*ran2(idum)*xdiff
c         ra1=ra1+roff/cdec
c         dec1=dec1-xdiff+2.*ran2(idum)*xdiff

         ntry=0
 887     continue
         rtry=rmin+ran2(idum)*rr
         dtry=dmin+ran2(idum)*dr
         ntry=ntry+1
         if(ntry.gt.100) then
            print *,"What???"
            goto 888
         endif
         do j=1,n
            diff=sqrt(cosd*cosd*(rtry-ra(j))**2+(dtry-dec(j))**2)
            diff=diff*3600.
            if(diff.lt.diffmin) goto 888
         enddo
         goto 887
 888     continue
         ra1=rtry
         dec1=dtry

         w1=wavel+ran2(idum)*wdiff
         f1=fluxl+ran2(idum)*fdiff
         f1=f1*1e-17
         write(11,1101) ra1,dec1,w1,f1,wsig
c         write(11,*) ra1,dec1,w1,f1
      enddo
      close(11)

      goto 667
      open(unit=1,file="rdetall",status='old')
      do i=1,10000
         read(1,*,end=667) a1,x2,x3,i4,i5,i6,a11,a12
         if(a11(5:7).eq.cifu) then
            write(*,1102) "rs1",x2,x3,35,i5,i6,a11,a12,
     $           "1.7 3 3.5 0.15 3 110"
         endif
      enddo
 667  continue
      close(1)

 1101 format(2(f10.6,1x),f7.2,1x,1pe12.5,1x,1pe11.3)
 1102 format(a3,2(1x,f10.6),1x,i2,1x,i4,1x,i2,1x,a11,1x,a12,1x,a20)

      end

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
 11             continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

      subroutine getseed(iseed)
      real xran(10000)
      iseed=-1
 668  continue
      open(unit=1,file='rannum.dat',status='old',err=11)
      nran=0
      do i=1,1000
         read(1,*,end=666) x1
         nran=nran+1
         xran(nran)=x1
      enddo
 666  close(1)
      open(unit=1,file='rannum.dat',status='unknown')
      do i=2,nran
         write(1,*) xran(i)
      enddo
      close(1)
      iseed=-nint(xran(1)*100000)
      if(nran.le.3) then
         call mkrannum(iseed)
         goto 668
      endif
      goto 667
      close(1)
 11   continue
      close(1)
      iseed=-nint(ran2(iseed)*100000.)
      call mkrannum(iseed)
 667  continue
      return
      end
      subroutine mkrannum(iseed)
      open(unit=1,file='rannum.dat',status='unknown')
      do i=1,1000
         write(1,*) ran2(iseed)
      enddo
      close(1)
      end
