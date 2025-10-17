

      read *,rms
      nsum=3

      call getseed(iseed)

      open(unit=1,file='in',status='old')
      open(unit=11,file='out',status='unknown')
      do j=1,100000
         sum=0
         do k=1,3
            read(1,*,end=666) x1,x2
            x2=x2+rms*gasdev(iseed)
            sum=sum+x2
            if(k.eq.2) wave=x1
         enddo
         write(11,*) wave,sum/3.
      enddo
 666  continue
      close(1)
      close(11)

      end

	subroutine getseed(iseed)
	real xran(10000)
	iseed=-1
 668	continue
	open(unit=1,file='rannum.dat',status='old',err=11)
	nran=0
	do i=1,1000
	   read(1,*,end=666) x1
	   nran=nran+1
	   xran(nran)=x1
	enddo
 666	close(1)
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
 11	continue
	close(1)
	iseed=-nint(ran2(iseed)*100000.)
 	call mkrannum(iseed)
 667	continue
	return
	end
	subroutine mkrannum(iseed)
	open(unit=1,file='rannum.dat',status='unknown')
	do i=1,1000
	   write(1,*) ran2(iseed)
	enddo
	close(1)
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
11      continue
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
