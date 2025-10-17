C	THIS PROGRAM LOOKS FOR ROTATION IN A GLOBULAR CLUSTER
C	BY FITTING THE FUNCTION V=VROT*COS(PA-PA0) TO THE VELOCITIES
C	OF THE STARS.  THE STARS ARE CHOSEN TO FALL IN A SPECIFIED
C	RANGE OF RADII
C
C	2-JAN-92: CONVERTED TO SUN, INCREASED ARRAYS TO 1000, AND MADE ANGLE
C	  FROM NORTH THROUGH EAST (ASSUMING THAT X IS EAST AND Y IS NORTH).
C

        parameter (nmax=200000)
	real th(nmax),v(nmax),r(nmax),xplt(nmax),yplt(nmax)
        real xl(2),yl(2),ve(nmax),ypltl(nmax),ypltu(nmax)
        real vsim(10000)
	character toplabel*48,lab1*40,filein*50
	parameter (radtodeg=57.29578,pi=3.141593,twopi=6.283185)

        ierr=1
        ilab=1
C
C
C	-- READ IN THE DATA
	N=0
	VPMAX=0.0
 111    call qc1('Velocity file ','rotate.def',filein)
        call savdef
        open(unit=1,file=filein,status='old',err=111)
C	-- GET THE MEAN VELOCITY
	call qr1('MEAN VELOCITY ','rotate.def',vmean)
1	CONTINUE
	N=N+1
2	READ (1,*,END=10) X,Y,V(N),ve(n)
C100	FORMAT(7X,F6.2,F6.2,1X,F7.2)
	R(N)=SQRT(X*X+Y*Y)
        th(n)=0.
        if(y.lt.0.) th(n)=pi
        if(x.gt.0.) then
           if(y.eq.0.) th(n)=pi/2.
           if(y.gt.0.) th(n)=atan(x/y)
           if(y.lt.0.) th(n)=pi-atan(-x/y)
        elseif(x.lt.0.) then
           if(y.eq.0.) th(n)=3*pi/2.
           if(y.gt.0.) th(n)=twopi-atan(-x/y)
           if(y.lt.0.) th(n)=pi+atan(x/y)
        endif
	V(N)=V(N)-VMEAN
	VPMAX=MAX(VPMAX,abs(V(N)))
	GOTO 1
C
10	CONTINUE
	CLOSE (UNIT=1)
	N=N-1
C
C	-- FIND THE ROTATION AMPLITUDE AND PHASE FOR STARS IN THE
C	   SPECIFIED RANGE OF RADIUS
	call qr2('RANGE OF RADII (ARCMINS) ','rotate.def',rl,rh)
        call savdef
	S1=0.0
	S2=0.0
	S3=0.0
	S4=0.0
	S5=0.0
	S6=0.0
	NB=0
        sumr=0.
	DO I=1,N
	  IF ((R(I).GT.RL).AND.(R(I).LT.RH)) THEN
            sumr=sumr+r(i)
	    STH=SIN(TH(I))
	    CTH=COS(TH(I))
	    NB=NB+1
	    S1=S1+CTH*CTH
	    S2=S2+STH*CTH
	    S3=S3+STH*STH
	    S4=S4+V(I)*CTH
	    S5=S5+V(I)*STH
	    S6=S6+V(I)*V(I)
	  END IF
	END DO
        print *,'Average Radius is ',sumr/float(nb)
C
	STH=S3*S1-S2*S2
	A=(S3*S4-S2*S5)/STH
	B=(S1*S5-S2*S4)/STH
	SIGVSQ=(S6+2.00*(A*(B*S2-S4)-B*S5)+A*A*S1+B*B*S3)/FLOAT(NB-1)
	SIGASQ=S3*SIGVSQ/STH
	SIGBSQ=S1*SIGVSQ/STH
C
	STH=A*A
	CTH=B*B
	VROT=SQRT(STH+CTH)
	TH0=ATAN(B/ABS(A))
	IF (A.LT.0.0) TH0=PI-TH0
	IF (TH0.LT.0.0) TH0=TWOPI+TH0
	SIGVR=SQRT(STH*SIGASQ+CTH*SIGBSQ)/VROT
	SIGTH0=SQRT(CTH*SIGASQ+STH*SIGBSQ)/(VROT*VROT)
	WRITE (6,101) NB,RL,RH,VROT,SIGVR,RADTODEG*TH0,RADTODEG*SIGTH0,
     $       SQRT(SIGVSQ)
101	FORMAT(1X,I5,' POINTS WITH',F6.1,'<R<',F6.1,' :'/1X,'VROT=',
     $        F6.2,'+-',F4.1,3X,'PA0=',F5.1,'+-',F3.0,3X,'DISP=',F6.2)
C
C	-- PLOT THE POINTS AND THE FIT
	CALL PGBEGIN(0,'?',1,1)
        call pgpap(0.,1.)
        call pgscf(2)
        call pgsch(1.2)
        vpmax=vpmax+vpmax/10.
c	CALL PGWINDOW(0.0,360.0,-VPMAX,VPMAX)
c	CALL PGBOX('BCNSTA',90.0,9,'BCNST',10.0,2)
        call pgenv(-10.,370.,-vpmax,vpmax,0,0)
        xl(1)=0.
        xl(2)=360.
        yl(1)=0.
        yl(2)=0.
        call pgline(2,xl,yl)
	TOPLABEL='VELOCITY VS. PA FOR STARS WITH        <R<       '
c	WRITE(TOPLABEL(32:38),102) RL
c	WRITE(TOPLABEL(42:48),102) RH
102	FORMAT(F7.2)
c	CALL PGLABEL('PA (DEG)','V-<V> (KM/S)','')
	CALL PGLABEL('PA (DEG)','fw/<fw>','')
c	if(ilab.eq.1) CALL PGLABEL('','',TOPLABEL)
C
	LAB1='V\Drot\U=     '//'\(2233)'
        print *,vrot,sigvr
	WRITE(LAB1(10:14),103) VROT
103	FORMAT(F5.1)
	WRITE(LAB1(24:28),104) SIGVR
104	FORMAT(F5.1)
	if(ilab.eq.1) CALL PGMTEXT('R',1.5,0.05,0.0,LAB1)
	LAB1='PA\Do\U=    '//'\(2233)'
	WRITE(LAB1(9:11),106) nint(RADTODEG*TH0)
106	FORMAT(i3)
	WRITE(LAB1(19:21),105) nint(RADTODEG*SIGTH0)
105	FORMAT(i3)
	if(ilab.eq.1) CALL PGMTEXT('R',1.5,0.60,0.0,LAB1)
C
C	-- THE POINTS
	NB=0
	DO I=1,N
	  IF ((R(I).GT.RL).AND.(R(I).LT.RH)) THEN
	    NB=NB+1
	    XPLT(NB)=RADTODEG*TH(I)
	    YPLT(NB)=V(I)
            ypltl(nb)=v(i)+ve(i)
            ypltu(nb)=v(i)-ve(i)
	  END IF
	END DO
	CALL PGPOINT(NB,XPLT,YPLT,17)
        if(ierr.eq.1) call pgerry(nb,xplt,ypltl,ypltu,1.)
        open(unit=11,file='rotate.out',status='unknown')
        do i=1,nb
           vel=vrot*cos(xplt(i)/radtodeg-th0)
           write(11,*) yplt(i)-vel,ve(i)
        enddo
        close(11)

        nsim=10000
        call getsim(nb,xplt,yplt,nsim,vsim)
        call sort(nsim,vsim)
        i90=nint(.90*nsim)
        i95=nint(.95*nsim)
        i99=nint(.99*nsim)
        print *,"Get 90,95,99 = ",vsim(i90),vsim(i95),vsim(i99)

C
C	-- THE FIT
        open(unit=11,file='rotate2.out',status='unknown')
	DO I=1,91
	  XPLT(I)=4.00*FLOAT(I-1)
	  YPLT(I)=VROT*COS(6.981317E-2*FLOAT(I-1)-TH0)
          write(11,*) xplt(i),yplt(i)
	END DO
        close(11)
        call pgsci(2)
	if(ilab.eq.1) CALL PGLINE(91,XPLT,YPLT)
        call pgsci(1)
        vrot=8.7
        th0=90./radtodeg
	DO I=1,91
	  yplt(I)=vrot*cos(6.981317E-2*float(i-1)-th0)
	END DO
        call pgsci(4)
        call pgsls(4)
c	CALL PGLINE(91,xplt,yplt)
        call pgsci(1)

	CALL PGEND
C
	STOP
C
C	-- FILE OPEN ERROR
98	CONTINUE
	WRITE (6,*) 'ERROR OPENING VELOCITY FILE'
	STOP
	END

        subroutine getsim(n,x,y,ns,ys)
        real x(n),y(n),ys(ns),yn(10000)
	parameter (radtodeg=57.29578,pi=3.141593,twopi=6.283185)
        idum=-1
        call biwgt(y,n,xb,xs)
        disp=xs
        print *
        print *,"With a dispersion = ",disp
        do is=1,ns
           do i=1,n
              yn(i)=0.+disp*gasdev(idum)
           enddo
           S1=0.0
           S2=0.0
           S3=0.0
           S4=0.0
           S5=0.0
           S6=0.0
           DO I=1,N
              radian=x(i)/radtodeg
              STH=SIN(radian)
              CTH=COS(radian)
              S1=S1+CTH*CTH
              S2=S2+STH*CTH
              S3=S3+STH*STH
              S4=S4+yn(i)*CTH
              S5=S5+yn(i)*STH
              S6=S6+yn(i)*yn(i)
           END DO
           STH=S3*S1-S2*S2
           A=(S3*S4-S2*S5)/STH
           B=(S1*S5-S2*S4)/STH
           SIGVSQ=(S6+2.00*(A*(B*S2-S4)-B*S5)+A*A*S1+B*B*S3)/FLOAT(N-1)
           SIGASQ=S3*SIGVSQ/STH
           SIGBSQ=S1*SIGVSQ/STH
           STH=A*A
           CTH=B*B
           ys(is)=SQRT(STH+CTH)
        enddo
        return
        end
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
