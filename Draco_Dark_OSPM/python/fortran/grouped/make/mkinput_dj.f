c     - 236.7796630859   48.9740753174  3536.69   64.2067856278   101
c-  236.7737415  48.9707814   ifu013    11.440   -13.220   161.402  -463.161  multi_412_013_043_LU_001.ixy  20240731T045328.8     exp01 

      parameter(nmax=100000)
      real ra(nmax),dec(nmax),wv(10),wh(10)
      character adum*6,a8*28,amp(nmax)*14,amp0*14

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
      wr=wv(5)-wv(1)
      
      open(unit=1,file='pos.txt',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,adum,x4,x5,x6,x7,a8
         n=n+1
         ra(n)=x1
         dec(n)=x2
         amp(n)=a8(7:20)
      enddo
 666  continue
      close(1)

      open(unit=1,file="indj",status='old')
      open(unit=11,file="out",status='unknown')
      do i=1,1000000
         read(1,*,end=667) x1,x2,x3,x4,i5
         dmin=1e10
         do j=1,n
            dist=sqrt((x1-ra(j))**2+(x2-dec(j))**2)
            if(dist.lt.dmin) then
               dmin=dist
               amp0=amp(j)
            endif
         enddo

 889     wtry=wv(1)+ran2(idum)*wr
         call xlinint(wtry,5,wv,wh,w0)
         xv=ran2(idum)
         if(xv.gt.w0) goto 889
         wsig=wtry

         write(11,*) x1,x2,x3,x4,wsig,i5,amp0(1:11)
      enddo
 667  continue
      close(1)
      close(11)

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
