
      parameter(nmax=100000)
      real ra(nmax),dec(nmax),dfw(nmax),dff(nmax),wv(5),wh(5),xl(nmax)
      character cdate*8,cshot*3,file1*120,cifu*3,a3*13,a1*3
      character a11*11,a12*12,file2*120

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

      open(unit=1,file='n1',status='old')
      read(1,*) n1
      close(1)
      call getseed(idum)
c      idum=n1

      open(unit=1,file="ndfsim",status='old')
      read(1,*) ndfsim,ndfsim0
      close(1)
      print *,ndfsim,idum

      wv(1)=2.0
      wv(2)=3.0
      wv(3)=4.0
      wv(4)=5.0
      wv(5)=6.0
      wh(1)=0.10
      wh(2)=0.15
      wh(3)=0.32
      wh(4)=0.30
      wh(5)=0.12

      file2="/work/00115/gebhardt/sim/simscript/dfX.txt"
c      write(file2(36:36),3001) ndfsim0
      write(file2(38:38),3001) ndfsim0
 3001 format(i1)
      open(unit=1,file=file2,status='old')
      ndfsimb=0
      do i=1,ndfsim
         read(1,*)
         ndfsimb=ndfsimb+1
      enddo
      do i=1,nsim
         read(1,*,end=880) i1,x2,x3,x4
         dfw(i)=1215.666*(1.+x2)
         dff(i)=x4
         xl(i)=x3
         ndfsimb=ndfsimb+1
      enddo
 880  continue
      close(1)
      open(unit=11,file="ndfsim",status='unknown')
      write(11,*) ndfsimb,ndfsim0
      close(11)

      xdiff=xdiff/3600.

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
      wr=wv(5)-wv(1)
      diffmin=1.5

      wdiff=waveh-wavel
      fdiff=fluxh-fluxl
      cdec=cos(dec(1)/57.3)
      open(unit=11,file='out',status='unknown')
      do i=1,nsim
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

c- now get wsig
 889     wtry=wv(1)+ran2(idum)*wr
         call xlinint(wtry,5,wv,wh,w0)
         xv=ran2(idum)
         if(xv.gt.w0) goto 889
         wsig=wtry

c         w1=wavel+ran2(idum)*wdiff
         f1=fluxl+ran2(idum)*fdiff
         w1=dfw(i)
c         f1=dff(i)
         f1=f1*1e-17
         write(11,1101) ra1,dec1,w1,f1,wsig,xl(i)
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

 1101 format(2(f10.6,1x),f7.2,1x,1pe12.5,1x,1pe11.3,1x,1pe11.3)
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
      open(unit=1,file='../rannum.dat',status='old',err=11)
      nran=0
      do i=1,1000
         read(1,*,end=666) x1
         nran=nran+1
         xran(nran)=x1
      enddo
 666  close(1)
      open(unit=1,file='../rannum.dat',status='unknown')
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
      open(unit=1,file='../rannum.dat',status='unknown')
      do i=1,1000
         write(1,*) ran2(iseed)
      enddo
      close(1)
      end
      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
