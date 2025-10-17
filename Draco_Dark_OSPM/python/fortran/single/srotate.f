	parameter(nmax=30000)
	real th(nmax),r(nmax),v(nmax),ev(nmax),ev2(nmax),s(nmax)
	real vin(nmax),ev2in(nmax),thin(nmax)
	real p(3),xi(3,3),ai(3,3),var(3,3),rtmp(nmax),vos(nmax)
	real siga(nmax),sigau(nmax),vrota(nmax),vrotau(nmax)
	real ra(nmax),th0a(nmax),th0au(nmax),xt(nmax),xb(nmax)
	real wk(nmax)
	integer indx(3),iwk(nmax)
        character file1*80
	parameter(RADTODEG=57.29578,PI=3.141593,
	1    TWOPI=6.283185,DEGTORAD=1.745329E-2)
c
        common /vels/ nin,ifit,vin,ev2in,thin,th0i

	ilog=1

c
c	-- Get the first guess of pa0
	call qr1('Initial guess for the PA ','srotate.def',th0i)
	call qi1('Fix (0) or Fit (1) PA ','srotate.def',ifit)

c -- ifit=0 forces the PA; ifit=1 fits for it; 
c    ifit=2 fits 2nd moment only;
c    ifit<0 sets all errors to be the same (-ifit)

 101    call qc1('Input file ','srotate.def',file1)
        open(unit=1,file=file1,status='old',err=101)
	call savdef
	th0i=degtorad*th0i
C
C	-- READ IN THE DATA
	N=0
	ravg=0.0
	rmin=1.0e10
	rmax=-1.0e10
1	CONTINUE
	N=N+1
2	READ (1,*,END=10) X,Y,V(N),ev(n)
	if(ifit.lt.0) ev(n)=-ifit
	R(n)=SQRT(X*X+Y*Y)
	IF (Y.EQ.0.) THEN
	  IF (X.EQ.0.) THEN
	    GOTO 2
	  ELSE IF (X.GT.0.) THEN
	    TH(N)=0.50*PI
	  ELSE
	    TH(N)=1.50*PI
	  END IF
	ELSE
	  TH(N)=ATAN(X/ABS(Y))
	END IF
	IF (Y.LT.0.0) TH(N)=PI-TH(N)
	IF (TH(N).LT.0.0) TH(N)=TWOPI+TH(N)
c	-- Find max and min r for scaling
	rmin=min(rmin,r(n))
	rmax=max(rmax,r(n))
	ravg=ravg+r(n)
	GOTO 1
C
10	CONTINUE
	CLOSE (UNIT=1)
	N=N-1
	ravg=ravg/float(n)
	if(ifit.lt.0) ifit=2
	write (6,*) 'range of radii: ',rmin,rmax,' avg radius = ',ravg
C
c	-- Start by finding the max likelihood mean and dispersion for the
c	   entire sample
	call dispmaxl(n,v,ev,vm0,vm0unc,sig0,sig0unc,covar0,ier)
	write (6,108) n,vm0,vm0unc,sig0,sig0unc,covar0,ier
 108	format(1x,i6,' stars: <v> = ',f7.2,'+-',f5.2,'  sig0 = ',
	1    f6.2,'+-',f5.2/13x,'covar = ',f7.5,' ier = ',i2)
c
c	-- Get the mean velocity to use
	call qr1('mean velocity to use ','srotate.def',vm0)
	call savdef
c
c	-- Scale the velocities by this mean and dispersion
	vpmax=0.0
	do i=1,n
	  v(i)=(v(i)-vm0)/sig0
	  ev2(i)=(ev(i)/sig0)**2
	  if (abs(v(i)).gt.vpmax) vpmax=abs(v(i))
c	  write (6,*) i,v(i),ev(i),ev2(i)
	end do
	vpmax=vpmax*sig0
c
c	-- Sort the data into radial order
	do i=1,n
	   rtmp(i)=r(i)
	enddo
	call sort3(n,r,th,v,wk,iwk)
	call sort2(n,rtmp,ev2)
c
c	-- Calculate the rotation in radial bins
50	continue
	call qi2('Number of stars per bin and min number',
	1    'srotate.def',npb,nmin)
	call savdef
	if (npb.lt.0) goto 90

	open(unit=11,file='srotate.out',status='unknown')
	open(unit=12,file='srotate2.out',status='unknown')

	ivela=1
	if(ivela.eq.1) then
	   call pgbegin(0,'?',1,1)
	else
	   call pgbegin(0,'?',1,1)
	endif
	call pgpap(0.,1.)
	call pgscf(2)
	call pgsch(1.5)
	call pgslw(2)

	imin=1-npb
	imax=n+npb
	n2=0
	do iall=imin,imax
	   jmin=max(1,iall-npb/2)
	   jmax=min(n,iall+npb/2)
	   nin=0
	   sum=0.
	   do j=jmin,jmax
	      nin=nin+1
	      vin(nin)=v(j)
	      ev2in(nin)=ev2(j)
	      thin(nin)=th(j)
	      sum=sum+r(j)
c	      print *,nin,vin(nin),ev2in(nin),thin(nin)
	   enddo
	   if(nin.lt.nmin) goto 766
	   ravg=sum/float(nin)
c	-- Set up to find the maximum likelihood
c	dispersion
	   p(1)=1.0
c	rotation amplitude
	   p(2)=0.1
c	position angle
	   p(3)=th0i
c	
	   do i=1,3
	      do j=1,3
		 xi(j,i)=0.0
	      end do
	      xi(i,i)=1.0
	   end do
	   ftol=1.0e-4
c
c	-- Find the maximum likelihood solution and unscale the parameters
	   if(ifit.eq.0) np=2
	   if(ifit.eq.1) np=3
	   if(ifit.eq.2) np=1
c	   print *,nin,p
	   call powell2(p,xi,np,3,ftol,iter,fret)
c	   print *,iter
c	write (6,109) iter,fret
c109	format(1x/'Found minimum in ',i3,' iterations, fret =',f7.2/1x)
c
	   vm=vm0
c	-- If the dispersion is negative, flip its sign (it only enters **2)
	   if (p(1).lt.0.) p(1)=-p(1)
	   sigma=sig0*p(1)
	   if (p(2).lt.0.) then
c	  -- force vrot to be positive by changing th0 by 180 degrees
	      p(2)=-p(2)
	      if(np.eq.3) then
	      if (p(3).gt.pi) then
		 p(3)=p(3)-pi
	      else
		 p(3)=p(3)+pi
	      end if
	      endif
	   end if
	   vrot=sig0*p(2)
	   th0=radtodeg*p(3)
C
c	-- Calculate s=sin(th-th0) with max like solution.
	   do i=1,nin
	      s(i)=sin(th(i)-p(3))
	   end do
c
c	-- Find the uncertainties in the parameters and the covariances.
c	   These are correct in the limit of large sample sizes.  The basic
c	   method is calculating the information matrix.  Its inverse is
c	   the covariance matrix.
c
c	-- Zero the I matrix elements and set up an identity matrix for
c	   later inversion.
	   do i=1,np
	      do j=1,np
		 ai(i,j)=0.
		 var(i,j)=0.
	      end do
	      var(i,i)=1.0
	   end do
c
c	-- Sum df(v_k)/dtheta_i/dtheta_j over the scaled velocities
c	   with the max likelihood parameters to get I_i_j.
	   p1sq=p(1)*p(1)
	   do i=1,nin
	      tmp=(v(i)-p(2)*s(i))/(p1sq+ev2(i))
	      dlfdp1=p(1)*(tmp*(v(i)-p(2)*s(i))-1.00)/(p1sq+ev2(i))
	      dlfdp2=tmp*s(i)
	      dlfdp3=-tmp*p(2)*cos(th(i)-p(3))
	      ai(1,1)=ai(1,1)+dlfdp1*dlfdp1
	      ai(1,2)=ai(1,2)+dlfdp1*dlfdp2
	      ai(1,3)=ai(1,3)+dlfdp1*dlfdp3
	      ai(2,2)=ai(2,2)+dlfdp2*dlfdp2
	      ai(2,3)=ai(2,3)+dlfdp2*dlfdp3
	      ai(3,3)=ai(3,3)+dlfdp3*dlfdp3
	   end do
	   ai(2,1)=ai(1,2)
	   ai(3,1)=ai(1,3)
	   ai(3,2)=ai(2,3)
c       
c	-- Find the covariance matrix V by inverting I (using numerical
c	   Recipes routines with the previously prepared identity matrix V)
	   call ludcmp(ai,np,3,indx,parity)
	   do i=1,np
	      call lubksb(ai,np,3,indx,var(1,i))
	   end do
c
c	-- Put the scales back in the covariance matrix
	   do i=1,np
	      do j=1,np
		 if (i.eq.3) then
		    var(i,j)=var(i,j)*radtodeg
		 else
		    var(i,j)=var(i,j)*sig0
		 end if
		 if (j.eq.3) then
		    var(i,j)=var(i,j)*radtodeg
		 else
		    var(i,j)=var(i,j)*sig0
		 end if
	      end do
	   end do
c
c	-- write out the results
	   vrotunc=sqrt(var(2,2))
	   sigmaunc=sqrt(var(1,1))
	   if(ifit.eq.1) th0unc=sqrt(var(3,3))
	   if(ifit.eq.0.or.ifit.eq.2) th0unc=0.
	   n2=n2+1
	   if(ilog.eq.1) then
	      ra(n2)=log10(ravg)
	   else
	      ra(n2)=ravg
	   endif
	   siga(n2)=sigma
	   sigau(n2)=sigmaunc
	   vrota(n2)=vrot
	   vrotau(n2)=vrotunc
	   th0a(n2)=th0
	   th0au(n2)=th0unc
	   vos(n2)=vrota(n2)/siga(n2)
	   write(11,1101) ra(n2),vrota(n2),vrotau(n2),siga(n2),
	1	sigau(n2),th0a(n2),th0au(n2),nin
	   write(12,*) 10**ra(n2),siga(n2),siga(n2)-sigau(n2),
	1	siga(n2)+sigau(n2)
 766	   continue
	enddo
 1101	format(f9.3,6(1x,f8.3),1x,i4)
	close(11)
	close(12)

	rmaxp=ra(n2)
c	rmaxp=log10(3.)
	if(ivela.eq.0) then
	if(ilog.eq.1) then
        call pgenv(ra(1),rmaxp,0.,15.,0,10)
c        call pgenv(ra(1),ra(n2),200.,440.,0,10)
	else
	   call pgenv(ra(1),ra(n2),0.,20.,0,0)
	endif
	call pgline(n2,ra,siga)
	endif
	do i=1,n2
	   xt(i)=siga(i)+sigau(i)
	   xb(i)=siga(i)-sigau(i)
	enddo
	if(ivela.eq.0) then
	call pgsls(4)
	call pgline(n2,ra,xt)
	call pgline(n2,ra,xb)
	call pgsls(1)
	call pglabel('R','\gs','')
	endif
	
	if(ilog.eq.1) then
	   call pgenv(log10(0.9),rmaxp,0.,14.,0,10)
c	   call pgenv(ra(1),rmaxp,0.,10.,0,10)
	else
	   call pgenv(ra(1),rmaxp,0.,10.,0,0)
	endif
	call pgslw(5)
	call pgline(n2,ra,vrota)
	do i=1,n2
	   xt(i)=vrota(i)+vrotau(i)
	   xb(i)=vrota(i)-vrotau(i)
	enddo
	call pgslw(3)
	call pgsls(4)
	call pgline(n2,ra,xt)
	call pgline(n2,ra,xb)
	call pgsls(1)
	call pgslw(2)
	call pglabel('Radius (arcmin)',
	1    'Rotation Amplitude (km s\U-1\D)','')

	if(ivela.eq.0) then
	if(ilog.eq.1) then
	   call pgenv(ra(1),rmaxp,-100.,320.,0,10)
	else
	   call pgenv(ra(1),rmaxp,-100.,320.,0,0)
	endif
	call pgline(n2,ra,th0a)
	do i=1,n2
	   xt(i)=th0a(i)+th0au(i)
	   xb(i)=th0a(i)-th0au(i)
	enddo
	call pgsls(4)
	call pgline(n2,ra,xt)
	call pgline(n2,ra,xb)
	call pgsls(1)
	call pglabel('R','PA','')
	if(ilog.eq.1) then
	   call pgenv(ra(1),rmaxp,0.,1.0,0,10)
	else
	   call pgenv(ra(1),rmaxp,0.,1.,0,0)
	endif
	call pgline(n2,ra,vos)
	call pglabel('R','V/\gs','')
	endif

	call pgend
c
90	continue
	STOP
C
C	-- FILE OPEN ERROR
98	CONTINUE
	WRITE (6,*) 'ERROR OPENING VELOCITY FILE'
	STOP
	END
c
c
	real function func(p)
	parameter(nmax=30000)
	real p(3),v(nmax),ev2(nmax),th(nmax)
	PARAMETER (TWOPI=6.283185)
c
	common /vels/ nv,ifit,v,ev2,th,th0i
c
c	-- Calculate -log(likelihood)
	func=0.
	sig2=p(1)*p(1)
	if(ifit.eq.0) p(3)=th0i
	if(ifit.eq.2) then
	   p(2)=0.
	   p(3)=th0i
	endif
	if(abs(p(2)/p(1)).lt.1.e-6) p(3)=th0i
c - get the weights
        wd=(float(nv)/2.)**3
        w0=float(nv)/2.
	do i=1,nv
          w=1.-(abs(float(i)-w0))**3/wd
c	  w=1.
	  func=func+w*(((v(i)-(p(2)*sin(th(i)-p(3))))**2)/
	1      (2.0*(sig2+ev2(i))) + 0.50*alog(twopi*(sig2+ev2(i))))
	enddo
c
c	write (6,*) 'func',p(1),p(2),p(3),func
	return
	end

	SUBROUTINE powell2(p,xi,n,np,ftol,iter,fret)
	INTEGER iter,n,np,NMAX,ITMAX
	REAL fret,ftol,p(np),xi(np,np),func
	EXTERNAL func
	PARAMETER (NMAX=20,ITMAX=1000)
C       U    USES func,linmin
	INTEGER i,ibig,j
	REAL del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
	fret=func(p)
	do 11 j=1,n
	   pt(j)=p(j)
 11	continue
	iter=0
 1	iter=iter+1
	fp=fret
	ibig=0
	del=0.
	do 13 i=1,n
	   do 12 j=1,n
	      xit(j)=xi(j,i)
 12	   continue
	   fptt=fret
	   call linmin(p,xit,n,fret)
	   if(abs(fptt-fret).gt.del)then
	      del=abs(fptt-fret)
	      ibig=i
	   endif
 13	continue
	if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
	if(abs(p(2)/p(1)).lt.1.e-6) then
	   p(2)=0.
	   p(3)=th0i
	   return
	endif
	if(iter.eq.ITMAX) print *,'powell exceeding maximum iterations'
	do 14 j=1,n
	   ptt(j)=2.*p(j)-pt(j)
	   xit(j)=p(j)-pt(j)
	   pt(j)=p(j)
 14	continue
	fptt=func(ptt)
	if(fptt.ge.fp)goto 1
	t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
	if(t.ge.0.)goto 1
	call linmin(p,xit,n,fret)
	do 15 j=1,n
	   xi(j,ibig)=xi(j,n)
	   xi(j,n)=xit(j)
 15	continue
	goto 1
	END

	subroutine dispmaxl(n,v,ev,vm,vmunc,sigma,sigmaunc,covar,ier)
c
c	Subroutine to find the dispersion and mean of a set of velocities
c	using maximum likelihood.  This approach automatically removes the
c	effect of the velocity meansurement uncertainties from the estimated
c	dispersion.  Note that the estimate of the dispersion is biased
c	low by at least sqrt((N-1)/N), but that this bias is always smaller
c	than the uncertainty in the dispersion.
c
c	The iterative scheme used to solve the maximum likelihood
c	equations if very similar to that used by Gunn and Griffin (1979,
c	AJ 84, 752) to fit the velocity scale parameter of King models.
c
c	The uncertainties for the mean velocity and the dispersion are
c	calculated from the diagonal elements of the covariance matrix and
c	are correct in the limit of large sample size, N (they are as
c	accurate as sigma/sqrt(2N) is for the simple dispersion estimator).
c	These uncertainties assume that the two parameters are independent
c	and so are reasonable as long as the normalized covariance is much
c	less than 1.  The covariance seems to always be small if the sample
c	size is reasonable and the velocity uncertainties are smaller than
c	the dispersion.  Simulations show that the uncertainties are
c	correct for samples at least as small as 20.  (The uncertainties are
c	correct on average; the uncertainty in the estimated uncertainty of
c	course increases as the sample size gets smaller.)  The simulations
c	also show that the bias in the estimated dispersion is always
c	larger than sqrt((N-1)/N).  For example, it is about twice this
c	when the velocity uncertainties are half as large as the dispersion
c	and N = 20 - 50.
c
c	For more discussion, see Pryor and Meylan 1993, in Structure and
c	Dynamics of Globular Clusters, edited by Meylan and Djorgovski,
c	ASP Conference Series Vol. 50, p. 357.  If you have questions or
c	comments, send email to Tad Pryor: pryor@physics.rutgers.edu.
c
c	Parameters:
c	  n: (input) sample size
c	  v: (input) array of radial velocities
c	  ev: (input) array of radial velocity uncertainties
c	  vm: (output) mean radial velocity
c	  vmunc: (output) uncertainty in mean velocity
c	  sigma: (output) velocity dispersion
c	  sigmaunc: (output) uncertainty in the velocity dispersion
c	  covar: (output) normalized covariance between vm and sigma
c	  ier: (output) error flag. 0 - no error; 1 - iteration failed
c	       to converge, no parameters calculated; 2 - singular
c	       information matrix, no uncertainties calculated.
c
	real v(n),ev(n),i11,i22,i12
	data icntm/1000/
c
c	-- Find the max. likelihood mean velocity and velocity scale.
	ier=0
	sum2=0.0
	sum3=0.0
c
c	-- For the first iteration ignore the velocity uncertainties.
	do 506 i=1,n
	  sum2=sum2+v(i)
	  sum3=sum3+v(i)*v(i)
506	continue
	vm=sum2/float(n)
	v02=(sum3-vm*vm*float(n))/float(n)
c
c	-- Now iterate with the uncertainties
	dvm=1.0
	dv02=1.0
	icnt=1
507	if (((dvm.gt.0.01).or.(dv02.gt.0.01)).and.(icnt.lt.icntm)) then
	  vmt=vm
	  v0t2=v02
	  sum1=0.0
	  sum2=0.0
	  sum3=0.0
	  sum4=0.0
	  sum5=0.0
	  sum6=0.0
	  do 508 i=1,n
	    tmp=1.00/(1.00+ev(i)*ev(i)/v0t2)
	    sum1=sum1+tmp
	    tmp2=tmp*tmp
	    sum2=sum2+tmp
	    sum3=sum3+v(i)*tmp
	    sum4=sum4+tmp2
	    sum5=sum5+v(i)*tmp2
	    sum6=sum6+v(i)*v(i)*tmp2
508	  continue
	  vm=sum3/sum2
	  v02=(sum6+vm*(vm*sum4-2.0*sum5))/sum1
	  v02=max(1.e-9,v02)
	  dvm=abs(vm-vmt)/vm
	  dv02=abs(v02-v0t2)/v02
	  icnt=icnt+1
	  goto 507
	end if
c
c	-- Return if convergence failed.
	if (icnt.eq.icntm) then
	  print *,'here'
	  ier=1
	  vm=0.
	  sigma=0.
	  vmunc=0.
	  sigmaunc=0.
	  covar=0.
	  return
	end if
c
	sigma=sqrt(v02)
c
c	-- Now estimate the uncertainities in the parameters
c	   Start by accumulating some sums one last time to estimate the
c	   information matrix, I.  I11 = d^2L/dvm^2,  I22 = d^2L/dsigma^2,
c	   I12 = I21 = d^2L/(dvm*dsigma).
	sum1=0.0
	sum2=0.0
	sum3=0.0
	sum4=0.0
	do 509 i=1,n
	  tmp=1.00/(1.00+ev(i)*ev(i)/v02)
	  sum1=sum1+tmp
	  tmp2=tmp*tmp
	  tmp3=(v(i)-vm)/sigma
	  sum2=sum2+tmp2*tmp3
	  tmp3=tmp3*tmp3
	  sum3=sum3+(tmp3+2.0)*tmp2
	  sum4=sum4+tmp3*tmp2*tmp
509	continue
	i11=sum1/v02
	i12=2.0*sum2/v02
	i22=(sum1-sum3+4.0*sum4)/v02
c
c	-- The next step is to invert the information matrix to get the
c	   covariance matrix V.  The uncertainty in vm is sqrt(V11) and the
c	   uncertainty in sigma is sqrt(V22).  The covariance is V12.
	det=i11*i22-i12*i12
	if (det.eq.0.) then
c	  -- Singular information matrix.
	  ier=2
	  vmunc=0.
	  sigmaunc=0.
	  covar=0.
	else
	  v11=i22/det
	  v22=i11/det
	  v12=-i12/det
	  vmunc=sqrt(v11)
	  sigmaunc=sqrt(v22)
	  covar=v12/sqrt(v11*v22)
	end if
c
c	-- Done
	return
	end

