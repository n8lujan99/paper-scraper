	PROGRAM ROSTAT

c*****************************************************************************
c
c    VERSION 1.2  February.1991
c
c--- this routine contains the current versions of statistical
c    routines being tested by T. Beers, K. Flynn, and K. Gebhardt for robust
c    estimation of simple statistics
c 
c*****************************************************************************
c
c-- CHANGES SINCE V1.0 (see update document for more information)
c
c-- 10/10/90 -- include the original estimate of the biweight when
c               error corrections are desired
c
c-- 10/10/90 -- change dip test so it doesn't loop forever
c
c-- 10/10/90 -- add check so KSREJ does not clip to N=5
c
c-- 10/11/90 -- changed calculation of acceleration constant
c
c--  1/ 1/91 -- changed b1 & b2 test and added probability calculation
c
c-- CHANGES SINCE V1.1 (see update document for more information)
c
c--  1/28/91 -- included the term of 1/1.349 in the calculation of the 
c               confidence interval for the median and fourth spread
c
c--  1/28/91 -- the calculation for the confidence interval for the 3-sigma
c               mean and stand. dev. was changed in confidence and PRINTER 
c               was updated
c
c--  1/28/91 -- PRINTER was changed so it now writes out the proper T-VALUES
c               which are calculated in the program
c
c--  1/28/91 -- Changed calculation of B1-test so it gives the left side 
c               probability when the skewness is less than 0
c
c--  1/30/91 -- Fixed default number of decimal places for output
c
c--  1/31/91 -- Added check so that SIGMA does not clip to N = 5
c
c--  2/ 1/91 -- Added check so that KURTSKEW is skipped when N < 9
c
c--  2/ 1/91 -- Changed calculation of B2-test so it gives the left side
c	        probability when the kurtosis is less than 3
c
c--  2/ 1/91 -- The subroutine KSREJ now calls EDFGOF directly, rather than
c               going through the entire NORMALITY subroutine
c
c--  2/ 1/91 -- Updated PRINTER so that exit diagnostics ( values = -999 )
c               are properly output
c
c-- CHANGES SINCE V1.2 (see update document for more information)
c
c--  4/18/91 -- Updated main program so that it will not call TAIL,KSREJ,TRANS
c               if n is less than 7, but so that it will still call SIGMA to
c               n = 5 (check is inside SIGMA)
c
c--  4/18/91 -- Updated CONFIDENCE so that it will return intervals = 0 if
c               the estimate of XSDEV = -999, indicating that it was not
c               calculated earlier in the program
c   	
c--  4/24/91 -- Corrected typo in data statement of WEXT.  Apparently W test
c               is unaffected.
c
c--  4/26/91 -- Added a check in the W probability calculation of WEXT to
c		guard against underflows in call to ALNORM
c
c--  4/29/91 -- Changed KURTSKEW so it can handle larger datasets
c		without blowing up.
c
c*****************************************************************************

        integer iunit2,iunit3,iunit4,ounit3,ounit4
        real mletter(8),mcconf(32,4),conf(4)
        dimension xdata(50000),odata(50000),xletter(15,2),depths(15)
	dimension xerr(50000),oerr(50000)
        dimension xmsc(4),xmdc(4),xtsc(4),xjkc(4),xjlc(4)
        dimension xxscl(4),xxscu(4),xjbc(4),xjblc(4),xjgc(4),xjglc(4)
        dimension clvalues(4),cuvalues(4),tvalues(12),test(26)
        dimension sletter(15),rletter(8),cletter(8)
        dimension r(50000*2)
        dimension zstar(50000),zbig(50000),jbig(50000)

	character filename*32

        data zero,cee /0.0,2.99792458E+05/
        common iunit2,iunit3,iunit4,ounit3,ounit4
	common /print/ jtnor,jtgap,jtcon,jtchi,jtasy,jtcorr,jtot
c        data iunit2,iunit3,iunit4,ounit3,ounit4 /2,3,4,9,10/
	iunit2=2
	iunit3=3
        iunit4=4
	ounit3=9
	ounit4=10

c--- open up the necessary files

        open(unit=2,file='/work/00115/gebhardt/progs/ttable.dat',
     $       status='old')
        open(unit=3,file='/work/00115/gebhardt/progs/cmin.dat',
     $       status='old')
        open(unit=4,file='/work/00115/gebhardt/progs/cplu.dat',
     $       status='old')
        open(unit=8,file='/work/00115/gebhardt/progs/corr.dat',
     $       status='old')
        open(unit=9,file='rostat.out',status='unknown')
        open(unit=10,file='rostat.res',status='unknown')

c--- get the data

 1      write(*,"('INPUT DATA FILE : '$)")
	read(*,'(a)') filename
	open(unit=20,file=filename,status='old',err=1)
	write(*,"('Including Velocity Errors?  1-yes, 0-no : '$)")
	read *,ierr
	nin=0
	j=0

c--- 	typical inputs

	write(*,"('What is number of simulations in Boot? : '$)")
	read *,nsims

	conf(1)=68	
	conf(2)=90
	conf(3)=95
	conf(4)=99

	do 331 i=1,10000
	if (ierr) 50,50,51
51	read(20,*,end=15) odata(i),oerr(i)
	goto 52
50	read(20,*,end=15) odata(i)
52	nin=nin+1
331	continue

15	xlbiwt = zero	

c-----  If cosmological errors are requested we assume by default we
c       are using ROSTAT for clusters of galaxies, and approximate errors
c	are calculated taking the biweight estimates of location and scale,
c       then reanalyzing the data set after scaling by CEE

	if (ierr.eq.1)  then		

	j=0
	do 392 i=1,nin
	if (odata(i).ne.-32368.) then
	j=j+1
	xerr(j)=oerr(i)
	xdata(j)=odata(i)
	endif
392	continue

        CALL XBIWT (XDATA,J,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1)

	vt2=xlbiwt

c--  make the cosmological correction here

	do 393 i=1,j
393	xdata(i)=(xdata(i)-xlbiwt)/(1.+xlbiwt/cee)

	else

	do 394 i=1,nin
	if (odata(i).ne.-32368.) then
	j=j+1
	xdata(j)=odata(i)
	endif
394	continue

	endif	
	n=j
395 	continue

c--- here we initialize the estimators

        xave      = zero               
        xsdev     = zero               
        xadev     = zero               
        xskew     = zero               
        xcurt     = zero               
        xvar      = zero             
        xdf       = zero               
	xdfs	  = zero	     
        xsave     = zero             
        xssdev    = zero             
	xkave 	  = zero	     
	xkdev     = zero	     
        xtrim5    = zero             
        xtrim10   = zero             
        xtrim20   = zero             
        xtrim30   = zero             
        xtrim40   = zero             
        xmed      = zero             
        xbmed     = zero             
        xtrimn    = zero             
        xmid      = zero             
        xmadm     = zero             
        xadm      = zero             

        xlbiwt    = zero             
        xsbiwt    = zero             
        xlbiwt1   = zero             
        xsbiwt1   = zero           
        xljack1   = zero           
        xsjack1   = zero           
        xlljack1  = zero           
        xlsjack1  = zero           
        xljack2   = zero           
        xsjack2   = zero           
        xlljack2  = zero           
        xlsjack2  = zero           
  	xljack6a  = zero	   
	xsjack6a  = zero	   
	xlljack6a = zero	   
	xlsjack6a = zero	   
        xljack6   = zero           
        xsjack6   = zero           
        xlljack6  = zero           
        xlsjack6  = zero           
        xljack7   = zero           

        xsjack7   = zero           
        xlljack7  = zero           
        xlsjack7  = zero           
        xsgap     = zero           
	tindex1   = zero	   
	tindex2   = zero	   
	pow       = zero	   
	tee       = zero	   


c--- obtain some initial information about the data set

c--- sort the data and find the set of letter values

        CALL LETTERS (XDATA,N,DEPTHS,XLETTER,MLETTER,
     +                SLETTER,RLETTER,CLETTER,NLETTER)
	
        xmed = xletter(1,1)        
	xfl  = xletter(2,1)		
	xfu  = xletter(2,2)		
	xdf  = abs(xfu-xfl)	        
	xdfs = xdf/1.349		

c--- perform calculations and report back

        CALL MOMENT (XDATA,N,XAVE,XADEV,XSDEV,XVAR,XSKEW,XCURT)

        CALL xTRIM (XDATA,N,.05,XTRIM5)

        CALL xTRIM (XDATA,N,.1,XTRIM10)

        CALL xTRIM (XDATA,N,.2,XTRIM20)

        CALL XMIDMEAN (XDATA,N,XMID)

        CALL xTRIM (XDATA,N,.3,XTRIM30)

        CALL xTRIM (XDATA,N,.4,XTRIM40)

        CALL BMED (XDATA,N,XBMED)

        CALL TRIMEAN (XDATA,N,XMED,XFL,XFU,XTRIMN)

        CALL XMAD (XDATA,N,XMED,XMADM)

        CALL XAD (XDATA,N,XMED,XADM)

        CALL XBIWT (XDATA,N,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1)

        CALL GAPPER (XDATA,N,XSGAP,ZSTAR,NBIG,JBIG,ZBIG)

	if (ierr.eq.1) then 		

	d1sum = 0
        do 55 I=1,n
55      d1sum = d1sum+1./(xerr(i)*xerr(i))        

        davg=(n/d1sum)**0.5

C----- obtain corrections for clusters of galaxies

	cicorr=davg*davg	
        scorr=davg*davg/(1.+vt2/cee)**2		
	sicorr = scorr*(1.+scorr/(2.*xsbiwt*xsbiwt))/n	

	endif
	

        CALL JACKNIFE  (XDATA,N,1,XLJACK1,XSJACK1,	
     +                 XLLJACK1,XLSJACK1,XLJACK2,       
     +		       XSJACK2,XLLJACK2,XLSJACK2)

        CALL JACKNIFE  (XDATA,N,6,XLJACK6a,XSJACK6a,	
     +                 XLLJACK6a,XLSJACK6a,XLJACK6,     
     +		       XSJACK6,XLLJACK6,XLSJACK6)

        CALL JACKNIFE  (XDATA,N,7,XLJACK7,XSJACK7,	
     +                 XLLJACK7,XLSJACK7,XJ71,XJ72,
     +                 XJ73,XJ74)

	if (n .ge. 7) then 
	        CALL TAIL (XDATA,N,TINDEX1,TINDEX2)
		CALL KSREJ (XDATA,N,XKAVE,XKDEV,NKTOT)
		CALL TRANS (N,NLETTER,XLETTER,POW)
	endif

        CALL SIGMA (XDATA,N,XSAVE,XSSDEV,NSTOT)	

c--	If isin is set equal to 1 then the canonical confidence interval
c       used will be calculated using the 3-sigma mean and standard
c       deviation.  If isin = 0, the mean and standard deviation will be used .

	isin=1

	if (isin.eq.0) then
	nsin=n
	xsin=xsdev
	else if (isin.eq.1) then
	nsin=nstot
	xsin=xssdev
	endif

769 	    CALL CONFIDENCE (N,NSIN,XSIN,XDFS,XSBIWT,XSJACK2,
     +           XLSJACK2,XSJACK6,XLSJACK6,XSJACK7,XLSJACK7,
     +           XMSC,XMDC,XTSC,XJKC,XJLC,XJBC,XJBLC,
     +           XJGC,XJGLC,XXSCL,XXSCU,
     +           CLVALUES,CUVALUES,TVALUES)

	if(nsims.eq.0) goto 779    

        CALL MCCON (XDATA,N,CONF,XAVE,XSDEV,XLBIWT,XSBIWT,
     +              MCCONF,NSIMS) 

c---  for mcconf(i,j) : j=1,2,3,4 --> interval percentages
c		        i= 1-8 are mean limits
c			i= 9-16 are biwt loc limits
c			i= 17-24 are the stand dev limits
c			i= 25-32 are biwt spread limits

779     nt = 26                    

        CALL NORMALITY (XDATA,N,XAVE,XSDEV,XSKEW,XCURT,
     +                  XSBIWT,XMED,TEST)

	CALL SYMMETRY (XDATA,N,TEE)

c--- fill up a array of results to be saved and printed   

        r(1) = n                   
        j=1

        do 100 i=1,n               
        j=j+1
        r(j) = xdata(i)
100     continue

        j=j+1
        r(j) = nletter             
        do 101 i=1,nletter         
	j=j+1                      
	r(j) = depths(i)	   
        j=j+1                         
        r(j)=xletter(i,1)          
        j=j+1
        r(j)=xletter(i,2)          
	j=j+1
	r(j) =mletter(i)	   
	j=j+1
	r(j)=sletter(i)            
	j=j+1
	r(j) = rletter(i)          
	j=j+1
	r(j) = cletter(i)	   
101     continue

        j=j+1                      
        r(j) = xave
        j=j+1
        r(j) = xadev           
        j=j+1
        r(j) = xsdev           
        j=j+1
        r(j) = xvar
        j=j+1
        r(j) = xskew          
        j=j+1
        r(j) = xcurt          


        j=j+1                 
        r(j) = xave           
        j=j+1
        r(j) = nstot          
        j=j+1
        r(j) = xsave          
	j=j+1
	r(j) = nktot		
	j=j+1
	r(j) = xkave		
        j=j+1
        r(j) = xtrim5           
        j=j+1
        r(j) = xtrim10          
        j=j+1
        r(j) = xtrim20          
        j=j+1
        r(j) = xmid             
        j=j+1
        r(j) = xtrim30          
        j=j+1
        r(j) = xtrim40          
        j=j+1
	r(j) = xmed             
        j=j+1
        r(j) = xbmed            
        j=j+1
        r(j) = xtrimn           
        j=j+1
        r(j) = xlbiwt           
        j=j+1
        r(j) = XLBIWT1          
        j=j+1
        r(j) = xljack1          
        j=j+1
        r(j) = xsjack1          
        j=j+1
        r(j) = xlljack1         
        j=j+1
        r(j) = xlsjack1         
        j=j+1
        r(j) = xlboot           


        j=j+1                   
        r(j) = xsdev            
        j=j+1
        r(j) = nstot            
        j=j+1
        r(j) = xssdev           
        j=j+1
	r(j) = nktot		
	j=j+1
	r(j) = xkdev		
	j=j+1
        r(j) = xadev            
        j=j+1
        r(j) = xadm             
        j=j+1
        r(j) = xdf              
        j=j+1
        r(j) = xdf/1.349        
        j=j+1
        r(j) = xmadm            
        j=j+1
        r(j) = xmadm/0.6745     
        j=j+1
        r(j) = xsbiwt           
        j=j+1
        r(j) = xsbiwt1          
        j=j+1
        r(j) = xljack2          
        j=j+1
        r(j) = xsjack2          
        j=j+1
        r(j) = xlljack2         
        j=j+1
        r(j) = xlsjack2         
        j = j+1
	r(j) = xljack6		
	j=j+1
	r(j) = xsjack6          
	j=j+1
	r(j) = xlljack6		
	j = j+1
	r(j) = xlsjack6         
	j = j+1
	r(j) = xljack7		
	j=j+1
	r(j) = xsjack7          
	j=j+1
	r(j) = xlljack7		
	j = j+1
	r(j) = xlsjack7         
	j = j+1
        r(j) = xsboot           

        j=j+1                   
        r(j) = tindex1
        j=j+1
        r(j) = tindex2          

	j=j+1			
	r(j) = xsgap

	j=j+1
	r(j) = vt2		
	j=j+1
	r(j) = ierr             

	jtnor = n+7*nletter+60

        do 103 i=1,nt           
        j=j+1                   
        r(j) = test(i)
103     continue


	jtgap = jtnor+nt
	ngap = n-1

	j=j+1
	r(j) = ngap		      

        do 112 i=1,ngap               
        j=j+1
	r(j)= zstar(i)
112     continue

	j=j+1
	r(j) = nbig		      
	do 113 i=1,nbig		      
	j=j+1			      
	r(j) = jbig(i)
	j=j+1
	r(j) = zbig(i)
113	continue

	jtcon = 2+jtgap+ngap+2*nbig

        do 104 i= 1,4                 
        j=j+1
        r(j) = xmsc(i)                
104     continue		      

        do 105,i=1,4
        j=j+1
        r(j) = xmdc(i)                
105     continue                      

        do 106 i=1,4
        j=j+1
        r(j) = xtsc(i)                
106     continue                      

        do 107 i=1,4
        j=j+1
        r(j) = xjkc(i)                
107     continue                      

        do 108 i=1,4
        j=j+1
        r(j) = xjlc(i)                
108     continue                      

	do 121 i = 1,4
	j = j+1			      
	r(j) = xjbc(i)		      
121	continue

	do 122 i = 1,4
	j=j+1			      
	r(j) = xjblc(i)		      
122	continue

	do 123 i = 1,4
	j=j+1			      
	r(j) = xjgc(i)		      
123	continue

	do 124 i = 1,4	
	j=j+1			      
	r(j) = xjglc(i)		      
124	continue

	do 125 i = 1,4
	j=j+1			      
	r(j) = mcconf(1,i)	      
	j=j+1
	r(j) = mcconf(2,i)
125	continue

	do 126 i = 1,4
	j=j+1			      
	r(j) = mcconf(3,i)	      
	j=j+1
	r(j) = mcconf(4,i)
126	continue

	do 127 i = 1,4
	j=j+1			      
	r(j) = mcconf(5,i)	      
	j=j+1
	r(j) = mcconf(6,i)
127	continue

	do 128 i = 1,4
	j=j+1			      
	r(j) = mcconf(9,i)	      
	j=j+1
	r(j) = mcconf(10,i)
128	continue

	do 129 i = 1,4
	j=j+1			      
	r(j) = mcconf(11,i)	      
	j=j+1
	r(j) = mcconf(12,i)
129	continue

	do 130 i = 1,4
	j=j+1			      
	r(j) = mcconf(13,i)	      
	j=j+1
	r(j) = mcconf(14,i)
130	continue

	do 131 i = 1,4
	j=j+1			      
	r(j) = mcconf(17,i)	      
	j=j+1
	r(j) = mcconf(18,i)
131	continue

	do 132 i = 1,4
	j=j+1			      
	r(j) = mcconf(19,i)	      
	j=j+1
	r(j) = mcconf(20,i)
132	continue

	do 133 i = 1,4
	j=j+1			      
	r(j) = mcconf(21,i)	      
	j=j+1
	r(j) = mcconf(22,i)
133	continue

	do 134 i = 1,4
	j=j+1			      
	r(j) = mcconf(25,i)	      
	j=j+1
	r(j) = mcconf(26,i)
134	continue

	do 135 i = 1,4
	j=j+1			      
	r(j) = mcconf(27,i)	      
	j=j+1
	r(j) = mcconf(28,i)
135	continue

	do 136 i = 1,4
	j=j+1			      
	r(j) = mcconf(29,i)	      
	j=j+1
	r(j) = mcconf(30,i)
136	continue

	do 137 i = 1,4
	j=j+1			      
	r(j) = mcconf(7,i)	      
	j=j+1
	r(j) = mcconf(8,i)
137	continue

	do 138 i = 1,4
	j=j+1			      
	r(j) = mcconf(15,i)	      
	j=j+1
	r(j) = mcconf(16,i)
138	continue

	do 139 i = 1,4
	j=j+1			      
	r(j) = mcconf(23,i)	      
	j=j+1
	r(j) = mcconf(24,i)
139	continue

	do 140 i = 1,4
	j=j+1			      
	r(j) = mcconf(31,i)	      
	j=j+1
	r(j) = mcconf(32,i)
140	continue

	j = j+1
	r(j)=isin		      
	jtchi = jtcon+165	      
				      
				      
        do 109 i=1,4
        j=j+1
        r(j) = xxscl(i)               
        j=j+1                         
        r(j) = xxscu(i)               
109     continue

        do 110 i=1,4                  
        j=j+1
        r(j) = clvalues(i)           
        j=j+1                        
        r(j) = cuvalues(i)           
110     continue

        do 111 i=1,12                
        j=j+1                        
        r(j) = tvalues(i)            
111	continue		     
			             
	
	jtasy = jtchi + 28

	j=j+1
	r(j) = tee		     
	j=j+1
	r(j) = pow		     

	jtcorr = jtasy + 2

	j=j+1
	r(j) = cicorr			
	j=j+1				
	r(j)= scorr			
	j=j+1				
	r(j) = sicorr			
					
        jtot = j

c--- print out the result array in a readible format

	iin=ounit3
        CALL PRINTER (iin,FILENAME,R)

        write (ounit4,*) jtnor,jtgap,jtcon,jtchi,jtasy,jtcorr,jtot
        do 116 j=1,jtot
	write (ounit4,*) r(j)
116     continue

999     stop
        end


c------------------------------------------------------------------------------
        SUBROUTINE BIASAD (EST,NSIMS,ACC,SIMVAL,Z,LIMIT1,LIMIT2,
     +                     LIMITA1,LIMITA2,IFAULT)
c------------------------------------------------------------------------------

c--- Finds adjustment required to implement bias-corrected percentile method
c
c    uses functions ALNORM and PPND Algorithms AS 66 and AS 111
c
        dimension simval(nsims)
        data half,two/0.5E0,2.0E0/
        j = 0
        k = 0
c
c--- m = (number of values less than est) + (half the
c         number equal to est), rounded up if not integral
c
        do 30 m = 1,nsims
        if (simval(m) - est) 10,20,30
 10     j = j+1
 20     k = k+1
 30     continue
        m = (j + k + 1)/2
        if (m .eq. 0) ifault = 3
        if (m .eq. nsims) ifault = 4
        if (ifault .gt. 0) return
        zed = two*ppnd(float(m)/float(nsims),ifault)
c
c--- ifail cannot exceed 0 simce m .ge. 1 and m .le. nsims - 1
c
        fn1 = float(nsims + 1)
        limit1 = fn1*alnorm(zed - z,.false.) + half
        limit2 = fn1*alnorm(zed + z,.false.) + half

c
c--   calculate limits from the acceleration constant
c
	xzed1 = zed/two + ((zed/two - z)/(1-acc*(zed/two-z)))
	xzed2 = zed/two + ((zed/two + z)/(1-acc*(zed/two+z)))
	limita1 = fn1*alnorm(xzed1,.false.) + half
	limita2 = fn1*alnorm(xzed2,.false.) + half

        return
        end


c------------------------------------------------------------------------------
        SUBROUTINE BMED (XDATA,N,XBMED)
c------------------------------------------------------------------------------

c--- This subroutine calculates the broadened median for a set of N
c    ORDERED statistics in XDATA. The value of the broadened median
c    is returned as XBMED.  XBMED is defined differently for different
c    values of N, and is not defined for N<5.  XBMED is calculated by
c    taking a weighted average of the inner depths as specified 
c    page 313 of UREDA. The weights are as follows:
c
c    5 <= n <= 12   n odd  -->   middle 3 ordered stats are averaged
c                   n even -->   middle 4 ordered stats are averaged
c                                with weights 1/6 1/3 1/3 1/6
c
c    n > 12         n odd  -->   middle 5 ordered stats are averaged
c                   n even -->   middle 6 ordered stats are averaged
c                                with weights 1/10,1/5,1/5,1/5,1/5,1/10
c
c    The broadened median can also be easily found by taking an alpha
c    trimmed mean with alpha defined as:
c
c    5 <= n <= 12   alpha = (.5 - 1.5/n)
c
c    n > 12         alpha = (.5 - 2.5/n)
        
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        dimension xdata(n)
        data n2,d6,d5,d3,d10,d1 /2,6.0,5.0,3.0,10.0,1.0/
        data n1,n3 /1,3/
 
        do3 = d1/d3
        do5 = d1/d5
        do6 = d1/d6
        do10 = d1/d10
 
c--- bmed is not defined for n<5
 
        if (n .le. 12) then
           if (float(n)/float(n2) - int(n/n2) .eq. 0) then
                i1 = n/n2 - n1
                i2 = n/n2
                i3 = n/n2 + n1
                i4 = n/n2 + n2
                xbmed = do6*(xdata(i1)+xdata(i4))+do3*(xdata(i2)
     +                  +xdata(i3))
           else
                i1 = int(n/n2)
                i2 = int(n/n2) + n1
                i3 = int(n/n2) + n2
                xbmed = (xdata(i1) + xdata(i2) + xdata(i3))*do3
           endif
        else
           if (float(n)/float(n2) - int(n/n2) .eq. 0) then
                i1 = n/n2 - n2
                i2 = n/n2 - n1
                i3 = n/n2
                i4 = n/n2 + n1
                i5 = n/n2 + n2
                i6 = n/n2 + n3
                xbmed = do10*(xdata(i1)+xdata(i6)) + do5*(xdata(i2)
     +                  +xdata(i3)+xdata(i4)+xdata(i5))
           else
                i1 = int(n/n2) - n1
                i2 = int(n/n2)
                i3 = int(n/n2) + n1
                i4 = int(n/n2) + n2
                i5 = int(n/n2) + n3
                xbmed = do5*(xdata(i1)+xdata(i2)+xdata(i3)+xdata(i4)
     +                  +xdata(i5))
           endif
 
        endif
 
        return
        end


c------------------------------------------------------------------------------
        SUBROUTINE BOOTSTRAP (X,N,IP,NSIMS,X1SIM,X2SIM)
c------------------------------------------------------------------------------
 
c--- This routine returns a BOOTSTRAP array (X1SIM,X2SIM) of the 
c    simulated values.  The initial seed is set by a call to the system clock.
c
c******************************************************************************
 
        real*4 x(n),xdata(50000),x1sim(nsims),x2sim(nsims)
        integer iran1(50000)

c--- fill up an array of random numbers of length n (between 1 and n)

c	xx=SECNDS(0.0)			
c	k1=2*jifix(10.*xx)+1
c	k2=2*jifix(100.*xx)+1

c	do 6 i=1,100
c	xx=RAN(K1)
c6	xx=RAN(K2)

	idum = -1

        do 20 j=1,nsims

        do 5 i=1,n
        iran1(i) = 1 + int(n*RAN1(IDUM))
5       continue
 
        do 10 i=1,n
        xdata(i) = x(iran1(i))
10      continue

	if (ip.eq.1 .or. ip.eq.2)then
	x1sim(j) = PARAM(XDATA,N,1)
	x2sim(j) = PARAM(XDATA,N,2)

	else if (ip.eq.6) then
	CALL XBIWT(XDATA,N,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1)
	x1sim(j) = xlbiwt
	x2sim(j) = xsbiwt
	endif

20      continue

        return
        end


c------------------------------------------------------------------------------
        SUBROUTINE CONFIDENCE (N,NSTOT,XSDEV,XDF,XSBIWT,XSJACK2,
     +                XLSJACK2,XSJACK6,XLSJACK6,XSJACK7,XLSJACK7,
     +                XMSC,XMDC,XTSC,XJKC,XJLC,XJBC,XJBLC,
     +                XJGC,XJGLC,XXSCL,XXSCU,
     +                CLVALUES,CUVALUES,TVALUES)
c------------------------------------------------------------------------------

c--- CONFIDENCE calculates the confidence intervals for a set of location
c    parameters, in this case the straight mean XAVE, the median XMED,
c    and the bi-weight estimator XLBIWT.  The confidence intervals are
c    returned in the variables:

c		LOCATIONS and their CONFIDENCE INTERVALS

c               XMSC:   uses mean and standard deviation
c               XMDC:   uses median and fourth spread
c               XTSC:   uses bi-weight location and spread
c
c
c		SCALES and their CONFIDENCE INTERVALS
c
c               XJKC:   uses jacknifed estimates 
c               XJLC:   uses log - jacknifed estimates
c               XJBC:   uses jacknifed biwt estimates
c               XJBLC:  uses log - jacknifed estimates
c               XJGC:   uses jacknifed gapped estimates
c               XJGLC:  uses log -jacknifed gapped estimates
c               XXSC:   uses standard deviation (chi-square)
c                       (returns lower and upper values (XXSCL,XXSCU))
c
c               NOTE:   the above variables are arrays with return the
c                       68% (1), 90% (2), 95% (3), and 99% (4)
c                       confidence values
c
c            TVALUES:   tvalues of the confidence estimator
c
c                       (1) 0.68t for XMSC, XJKC, and XJLC
c                       (2) 0.90t for XMSC, XJKC, and XJLC
c                       (3) 0.95t for XMSC, XJKC, and XJLC
c                       (4) 0.99t for XMSC, XJKC, and XJLC
c
c                       (5) 0.68t for XMDC
c                       (6) 0.90t for XMDC
c                       (7) 0.95t for XMDC
c                       (8) 0.99t for XMDC
c
c                       (9) 0.68t for XTSC
c                       (10) 0.90t for XTSC
c                       (11) 0.95t for XTSC
c                       (12) 0.99t for XTSC
c
c
c            CVALUES:   cvalues of the confidence estimator
c
c                       (1) 0.68c for XXSC (u and l)
c                       (2) 0.90c fpr XXSC
c                       (3) 0.95c for XXSC
c                       (4) 0.99c for XXSC
c
c
c--- NOTE that the appropriate t-values are also calculated and returned
c    these values are also useful for standard hypothesis testing
c
c
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        integer n
        dimension tvalues(12),clvalues(4),cuvalues(4)
        dimension xmsc(4),xmdc(4),xtsc(4),xjkc(4),xjlc(4),
     +       xxscl(4),xxscu(4),xjbc(4),xjblc(4),xjgc(4),xjglc(4)

c--- obtain the appropriate t-values and c-values for each test 
 
        ndof = n-1
	nsdof = nstot-1
        tprob = .68
        t68 = TVALUE (NDOF,TPROB)
	ts68 = TVALUE (NSDOF,TPROB)
        c68l = CVALUE (0,NSDOF,TPROB)
        c68u = CVALUE (1,NSDOF,TPROB)
        tprob = .90
        t90 = TVALUE (NDOF,TPROB)
        ts90 = TVALUE (NSDOF,TPROB)
        c90l = CVALUE (0,NSDOF,TPROB)
        c90u = CVALUE (1,NSDOF,TPROB)
        tprob = .95
        t95 = TVALUE (NDOF,TPROB)
        ts95 = TVALUE (NSDOF,TPROB)
        c95l = CVALUE (0,NSDOF,TPROB)
        c95u = CVALUE (1,NSDOF,TPROB)
        tprob = .99
        t99 = TVALUE (NDOF,TPROB)
        ts99 = TVALUE (NSDOF,TPROB)
        c99l = CVALUE (0,NSDOF,TPROB)
        c99u = CVALUE (1,NSDOF,TPROB)

c--- XMSC, XJKC, and XJLC use the straight percentage points of the
c    t-distribution with N-1 degrees of freedom
 
        tvalues(1) = t68 
        tvalues(2) = t90
        tvalues(3) = t95
        tvalues(4) = t99 

c--- XMDC uses the approximation t* = t(n-1 DOF)/1.075
 
        tvalues(5) = t68/1.075
        tvalues(6) = t90/1.075
        tvalues(7) = t95/1.075
        tvalues(8) = t99/1.075
 
c--- XTSC uses the approximation t* = t with int(0.7(n-1)) DOF
 
        idof = int (0.7*ndof)
        tprob = .68
        tvalues(9) = TVALUE (IDOF,TPROB)
        tprob = .90
        tvalues(10) = TVALUE (IDOF,TPROB)
        tprob = .95
        tvalues(11) = TVALUE (IDOF,TPROB)
        tprob = .99
        tvalues(12) = TVALUE (IDOF,TPROB)
 
c--- XXSC uses the percentage points of the chi-square distribution

        clvalues(1) = c68l
        clvalues(2) = c90l
        clvalues(3) = c95l
        clvalues(4) = c99l

        cuvalues(1) = c68u
        cuvalues(2) = c90u
        cuvalues(3) = c95u
        cuvalues(4) = c99u

c--- obtain interval estimators for central locations

	if(xsdev .ne. -999) then 
	        xmsc(1) = ts68*xsdev/nstot**0.5
       		xmsc(2) = ts90*xsdev/nstot**0.5
       		xmsc(3) = ts95*xsdev/nstot**0.5
        	xmsc(4) = ts99*xsdev/nstot**0.5
 	endif

        xmdc(1) = tvalues(5)*xdf/n**0.5
        xmdc(2) = tvalues(6)*xdf/n**0.5
        xmdc(3) = tvalues(7)*xdf/n**0.5
        xmdc(4) = tvalues(8)*xdf/n**0.5
 
        xtsc(1) = tvalues(9)*xsbiwt/n**0.5
        xtsc(2) = tvalues(10)*xsbiwt/n**0.5
        xtsc(3) = tvalues(11)*xsbiwt/n**0.5
        xtsc(4) = tvalues(12)*xsbiwt/n**0.5

c--- obtain interval estimators for scale

        xjkc(1) = tvalues(1)*xsjack2
        xjkc(2) = tvalues(2)*xsjack2
        xjkc(3) = tvalues(3)*xsjack2
        xjkc(4) = tvalues(4)*xsjack2

        xjlc(1) = tvalues(1)*xlsjack2
        xjlc(2) = tvalues(2)*xlsjack2
        xjlc(3) = tvalues(3)*xlsjack2
        xjlc(4) = tvalues(4)*xlsjack2
 
        xjbc(1) = tvalues(1)*xsjack6
        xjbc(2) = tvalues(2)*xsjack6
        xjbc(3) = tvalues(3)*xsjack6
        xjbc(4) = tvalues(4)*xsjack6

        xjblc(1) = tvalues(1)*xlsjack6
        xjblc(2) = tvalues(2)*xlsjack6
        xjblc(3) = tvalues(3)*xlsjack6
        xjblc(4) = tvalues(4)*xlsjack6
 
        xjgc(1) = tvalues(1)*xsjack7
        xjgc(2) = tvalues(2)*xsjack7
        xjgc(3) = tvalues(3)*xsjack7
        xjgc(4) = tvalues(4)*xsjack7

        xjglc(1) = tvalues(1)*xlsjack7
        xjglc(2) = tvalues(2)*xlsjack7
        xjglc(3) = tvalues(3)*xlsjack7
        xjglc(4) = tvalues(4)*xlsjack7

	if(xsdev .ne. -999) then
        	xxscu(1) = sqrt(xsdev**2/(clvalues(1)/float(nsdof)))
        	xxscu(2) = sqrt(xsdev**2/(clvalues(2)/float(nsdof)))
       	        xxscu(3) = sqrt(xsdev**2/(clvalues(3)/float(nsdof)))
	        xxscu(4) = sqrt(xsdev**2/(clvalues(4)/float(nsdof)))
       
	        xxscl(1) = sqrt(xsdev**2/(cuvalues(1)/float(nsdof)))
       		xxscl(2) = sqrt(xsdev**2/(cuvalues(2)/float(nsdof)))
        	xxscl(3) = sqrt(xsdev**2/(cuvalues(3)/float(nsdof)))
        	xxscl(4) = sqrt(xsdev**2/(cuvalues(4)/float(nsdof)))
	endif

        return
        end


c------------------------------------------------------------------------------
        SUBROUTINE DIPTST (X,N,DIP,XL,XU,IFAULT,GCM,LCM,MN,MJ)
c------------------------------------------------------------------------------

c--- Algorithm AS 217 Appl. Stat. (1985) Vol.34, No.3
c    Hartigan, P.M.
c
c--- This routine provides an estimator which guages the probability that
c    the data is drawn from a parent population having a single mode.  As
c    such, it can be used as a conservative test for multi-modality of
c     a distribution.

c--- Does the dip calculation for an ordered vector X
c    using the greatest convex minorant and the least concave
c    majorant, skepping through the data using the change
c    points of these distributions.  It returns the dip
c    statistic 'DIP' and the modal interval (XL,XU)
c
c--- NOTE that the array X is expected to already be sorted low
c    to high on input
c
c******************************************************************************
 
        dimension x(n), mn(n),mj(n),lcm(n)
        integer gcm(n), high
c
c--- check that n is positive
c
        ifault = 1
        if (n .le. 0) return
        ifault = 0
c
c--- check if n is one
c
        if (n .eq. 1) goto 4
c
c--- check that x is sorted
c
        ifault = 2
        do 3 k = 2,n
        if (x(k) .lt. x(k-1)) return
3       continue
        ifault = 0
c
c--- check for all values of x identical and for 1 .lt. n .lt. 4
c
        if (x(n) .gt. x(1) .and. n .ge. 4) goto 5
4       xl = x(1)
        xu = x(n)
        dip = 0.0
        return

c
c--- low contains the index of the current estimate of the 
c    lower end of the modal interval, high contains the
c    index for the current upper end.
c

5       fn = float(n)
        low = 1
        high = n
        dip = 1.0/fn
        xl = x(low)
        xu = x(high)
c
c--- establish the indices over which combination is
c    necessary for the convex minorant fit
c
        mn(1) = 1
        do 28 j=2,n
        mn(j) = j-1
25      mnj = mn(j)
        mnmnj = mn(mnj)
        a = float(mnj - mnmnj)
        b = float(j - mnj)
        if (mnj .eq. 1 .or. (x(j) - x(mnj))*a .lt. (x(mnj)-x(mnmnj))
     +                                             *b) goto 28
        mn(j) = mnmnj
        goto 25
28      continue
c
c--- establish the indices over which combination is
c    necessary for the concave majorant fit.
c
        mj(n) = n
        na = n-1
        do 34 jk = 1,na
        k = n - jk
        mj(k) = k + 1
32      mjk = mj(k)
        mjmjk = mj(mjk)
        a = float(mjk - mjmjk)
        b = float(k - mjk)
        if (mjk .eq. n .or. (x(k) - x(mjk))*a .lt. (x(mjk)-x(mjmjk))
     +                                                *b) goto 34
        mj(k) = mjmjk
        goto 32
34      continue
c
c--- start the cycling
c    collect the change points for the gcm from high to low
c
	icount=0
40      ic = 1
        gcm(1) = high
42      igcm1 = gcm(ic)
        ic = ic + 1
        gcm(ic) = mn(igcm1)
        if (gcm(ic) .gt. low) goto 42
        icx = ic
c
c--- collect the change points for the lcm from low to high
c
        ic = 1
        lcm(1) = low
44      lcm1 = lcm(ic)
        ic = ic + 1
        lcm(ic) = mj(lcm1)
        if(lcm(ic) .lt. high) goto 44
        icv = ic
c
c--- icx,ix,ig are counters for the convex minorant
c    icv,iv,ih are counters for the concave majorant
c
        ig = icx
        ih = icv
c
c--- find the largest distance greater the 'DIP'
c    between the gcm and the lcm from low th high
c
        ix = icx - 1
        iv = 2
        d = 0.0
        if(icx .ne. 2 .or. icv .ne. 2) goto 50
        d=1.0/fn
        goto 65
50      igcmx = gcm(ix)
        lcmiv = lcm(iv)
        if (igcmx .gt. lcmiv) goto 55
c
c--- if the next point of either the gcm or lcm is
c    from the lcm then calculate distance here
c
        lcmiv1 = lcm(iv - 1)
        a = float(lcmiv - lcmiv1)
        b = float(igcmx - lcmiv1 - 1)
        dx = (x(igcmx) - x(lcmiv1)*a)/(fn*(x(lcmiv) - x(lcmiv1)))
     +                                              - b/fn
        ix = ix - 1
        if (dx .lt. d) goto 60
        d = dx
        ig = ix +1
        ih = iv
        goto 60
c
c--- if the next point of either the gcm or lcm is
c    from the gcm then calculate distance here
c
55      lcmiv = lcm(iv)
        igcm = gcm(ix)
        igcm1 = gcm(ix + 1)
        a = float(lcmiv - igcm1 + 1)
        b = float(igcm - igcm1)
        dx = a/fn - ((x(lcmiv) - x(igcm1))*b)/(fn*(x(igcm)
     +                                  - x(igcm1)))
        iv = iv + 1
        if (dx .lt. d) goto 60
        d = dx
        ig = ix + 1
        ih = iv - 1
60      if (ix .lt. 1) ix = 1
        if (iv .gt. icv) iv = icv
        if (gcm(ix) .ne. lcm(iv)) goto 50
65      if (d .lt. dip) goto 100
c
c--- caclulate the dips for the current low and high
c
c    the dip for the convex minorant
c
        dl = 0.0
        if (ig .eq. icx) goto 80
        icxa = icx - 1
        do 76 j = ig,icxa
        temp = 1.0/ fn
        jb = gcm(j+1)
        je = gcm(j)
        if (je - jb .le. 1) goto 74
        if (x(je) .eq. x(jb)) goto 74
        a = float(je - jb)
        const = a/(fn*(x(je) - x(jb)))
        do 72 jr = jb,je
        b = float(jr - jb +1)
        t = b/ fn - (x(jr) - x(jb)) * const
        if (t .gt. temp) temp = t
72      continue
74      if (dl .lt. temp) dl = temp
76      continue
c
c--- the dip for the concave majorant
c
80      du = 0.0
        if (ih .eq. icv) goto 90
        icva = icv - 1
        do 88 k = ih,icva
        temp = 1.0/fn
        kb = lcm(k)
        ke = lcm(k + 1)
        if (ke - kb .le. 1) goto 86
        if (x(ke) .eq. x(kb)) goto 86
        a = float(ke - kb)
        const = a/ (fn * (x(ke) - x(kb)))
        do 84 kr = kb,ke
        b = float(kr - kb - 1)
        t = (x(kr) - x(kb)) * const - b/fn
        if (t .gt. temp) temp = t
84      continue
86      if (du .lt. temp) du = temp
88      continue
c
c---    determine the current maximum
c
90      dipnew = dl
        if (du .gt. dl) dipnew = du 
        if (dip .lt. dipnew) dip = dipnew
        low = gcm(ig)
        high = lcm(ih)
c
c--- recycle
c
	icount=icount+1
	if(icount.gt.100) then
 	 dip=-999
	 xl=-999
	 xu=-999

	return
	endif

        goto 40
c
100     dip = .5 * dip
        xl = x(low)
        xu = x(high)

        return
        end


c------------------------------------------------------------------------------
	SUBROUTINE EDFGOF (N, X, ITEST, XBAR, S2, D, V, W2, W2P, U2,
     +                     U2P, A2, A2P, IFAULT)
c------------------------------------------------------------------------------

c---	Algorithm AS 248.1  Applied Statistics (91989), Vol. 38. No. 3
c
c	Tests a sample for uniformity, normality, or exponentiality
c	using goodness-of-fit statistics based on the empirical 
c	distribution function

	integer i, ifault, itest, n
c	real x(n), z(50000), xbar, s2, d, v, w2, w2p, u2, u2p,
 	real x(50000), z(50000), xbar, s2, d, v, w2, w2p, u2, u2p, 
     +      a2, a2p, rn,
     +      rootn, rn2, zero, pt01, pt05, pt1, pt11, pt12, pt155,
     +      pt16, pt2, pt24, pt26, pt3, pt35, pt4, pt5, pt6, pt75, pt8,
     +      pt82, pt85, one, twop25, wcut2(3), wcoef2(4,3), wcut3(3),
     +      wcoef3(4,3), ucut2(3), ucoef2(4,3), ucut3(3), ucoef3(4,3),
     +      acut2(3), acoef2(4,3), acut3(3), acoef3(4,3)
	data zero, pt01, pt1, pt05, pt11, pt12, pt155, pt16, pt2,
     +      pt24 / 0.0, 0.01, 0.05, 0.1, 0.11, 0.12, 0.155, 0.16, 0.2,
     +      0.24/
	data pt26, pt3, pt35, pt4, pt5, pt6, pt75, pt8, pt82, pt85,
     +      one / 0.26, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75, 0.8, 0.82, 0.85,
     +      1.0 /
	data twop25 / 2.25 /
	data wcut2, wcut3, ucut2, ucut3, acut2, acut3 / 0.0275, 0.051,
     +      0.092, 0.035, 0.074, 0.160, 0.0262, 0.048, 0.094, 0.029,
     +      0.062, 0.120, 0.200, 0.340, 0.600, 0.260, 0.510, 0.950 /
	data ((wcoef2(i,j), j=1,3), i=1,4) / -13.953, 775.5, -12542.61,
     +      -5.903, 179.546, -1515.29, 0.886, -31.62, 10.897, 1.111,
     +      -34.242, 12.832 /
	data ((wcoef3(i,j), j=1,3), i=1,4) / -11.334, 459.098, -5652.1,
     +      -5.779, 132.89, -866.58, 0.586, -17.87, 7.417, 0.447,
     +      -16.592, 4.849 /
	data ((ucoef2(i,j), j=1,3), i=1,4) / -13.642, 766.31,
     +      -12432.74, -6.3328, 214.57, -2022.28, 0.8510, -32.006, 
     +      -3.45, 1.325, -38.918, 16.45 /
	data ((ucoef3(i,j), j=1,3), i=1,4) / -11.703, 542.5, -7574.59,
     +      -6.3288, 178.1, -1399.49, 0.8071, -25.166, 8.44, 0.7663,
     +      -24.359, 4.539 /
	data ((acoef2(i,j), j=1,3), i=1,4) / -13.436, 101.14, -223.73,
     +      -8.318, 42.796, -59.938, 0.9177, -4.279, -1.38, 1.2937,
     +      -5.709, 0.0186 /
	data ((acoef3(i,j), j=1,3), i=1,4) / -12.2204, 67.459, -110.3,
     +      -6.1327, 20.218, -18.663, 0.9209, -3.353, 0.300, 0.731,
     +      -3.009, 0.15 /

	ifault = 1
	if (itest .lt. 1 .or. itest .gt. 3) return
	ifault = 2
	if ( n .lt. 2 .or. (itest .eq. 2 .and. n .eq. 2)) return
	ifault = 3
	do 10 i= 2, n
	   if (x(i) .lt. x(i-1)) return
10	continue
	rn = n
	rootn = sqrt(rn)
	rn2 = rn * rn
	if(itest .eq. 1) then
c
C	   Test for uniformity (F(X) is completely specified)
c
	   ifault = 0
	   CALL STATS (N, X, D, V, W2, U2, A2, IFAULT)
	   if (ifault .ne. 0) return
C
C	   Modifications when F(X) is completely specified
C
	   d = d * (rootn + pt12 + pt11 / rootn)
	   v = v * (rootn + pt155 + pt24 / rootn)
	   w2 = (w2 - pt4 / rn + pt6 / rn2) * (one + one / rn)
	   u2 = (u2 - pt1 / rn + pt1 / rn2) * (one + pt8 / rn)
 	   return
	else
C
C 	   Estimate the mean by XBAR
C
 	   xbar = zero
	   do 20 i = 1, n
	      xbar = xbar + x(i)
20	   continue
	   xbar = xbar / rn
	   if (itest .eq. 2) then
C
C	   Test for normality (MU and SIGMA**2 unspecified)
C	   First estimate the variance by S**2
C
	   	ifault = 5
		s2 = zero
		do 30 i = 1, n
		   s2 = s2 + (x(i) - xbar) ** 2
30		continue
		s2 = s2 / (rn - one)
		if (s2 .le. zero) return
C
C	   Compute Z(I) = F((XI)-XBAR)/S)
C
		s = sqrt(s2)
		do 40 i = 1, n
		   z(i) = ALNORM((X(I) - XBAR) / S, .FALSE.)
		if (z(i).le.0) z(i)=1.E-07
		if (z(i).ge.1) z(i)=1.-1.E-07
40	  	continue
		CALL STATS (N, Z, D, V, W2, U2, A2, IFAULT)
		if (ifaulT .NE. 0) return
C
C	   Modifications when F(X) is the normal distribution
C
		d = d * (rootn - pt01 + pt85 / rootn)
		v = v * (rootn + pt05 + pt82 / rootn)
		w2 = w2 * (one + pt5 / rn)
		w2p=0
		if (w2.gt.2) goto 41  		
		w2p = PVALUE (W2, WCUT2, WCOEF2)
41		u2 = u2 * (one + pt5 / rn)
		u2p=0
		if (u2.gt.2) goto 42		
		u2p = PVALUE (U2, UCUT2, UCOEF2)
42		a2 = a2 * (one + pt75 / rn + twop25 / rn2)
		a2p=0
		if (a2.gt.2) return
		a2p = PVALUE (A2, ACUT2, ACOEF2)
		return
	   else
C
C	Test for exponentiality (scale parameter unspecified)
C
		ifault = 6
		if (xbar .le. zero) return
		do 50 i = 1, n
		   z(i) = one - exp(-x(i) / xbar)
50		continue
		CALL STATS (N, Z, D, V, W2, U2, A2, IFAULT)
		if (ifault .ne. 0) return
C
C	Modifications when F(X) is the exponential distribution
C
		d = (d - pt2 / rn) * (rootn + pt26 + pt5 / rootn)
		v = (v - pt2 / rn) * (rootn + pt24 + pt35 / rootn)
		w2 = w2 * (one + pt16 / rn)
		w2p = PVALUE (W2, WCUT3, WCOEF3)
		u2 = u2 * (one + pt16 / rn)
		u2p = PVALUE (U2, UCUT3, UCOEF3)
		a2 = a2 * (one + pt3 / rn)
		a2p = PVALUE (A2, ACUT3, ACOEF3)
		return
	   endif
	endif
	end

	
c------------------------------------------------------------------------------
        SUBROUTINE GAPPER (XDATA,N,XSGAP,ZSTAR,NBIG,JBIG,ZBIG)
c------------------------------------------------------------------------------

c---    This routine returns a list of gaps in a data set, weighted by
c       their location with respect to the middle of the data.  It also
c       obtains values for these gaps relative to their average size.  This
c       normalized gap provides a measure of how large a given gap is
c       compared to the rest of the gaps.  Wainer and Schacht (1978;
c       Psychometrika 43, 203) assess the liklihood that a gap of given
c       size could occur in a normal distribution.

c                ZSTAR -- weighted and normalized gaps
c                NBIG  -- number of "big" gaps (> 2.25)
c                JBIG  -- indices of lower data value giving rise
c                         to a "big" gap
c                ZBIG  -- the size of the "big" gap

c       Thr routine also reports an estimate of scale , XSGAP, based on
c       the sum of the weighted gaps (see Wainer and Thissen 1976;
c       Psychometrika 41, 9.)

c***************************************************************************

        data pi /3.1415927/
	real xdata(n),zstar(n),ygap(50000),zbig(50000)
	real big(50000)
	integer jbig(50000)

c--- determine gap, weight, weighted gap vectors, and sum of gaps

        ysum = 0.
        do 10 i=1,n-1
        gap = xdata(i+1) - xdata(i)
        weight = i*(n-i)
        ygap(i) = (gap*weight)**0.5
	big(i) = i
        ysum = ysum + gap*weight
10      continue

c--- determine mid-mean of the weighted gaps

        CALL SORT2 (N-1,YGAP,BIG)

        CALL XMIDMEAN (YGAP,N-1,GMID)

c--- obtain standardized weighted gaps

	if(gmid.eq.0) gmid = 1

        do 20 i=1,n-1
        zstar(i) = ygap(i)/gmid
20      continue

c--- save gaps larger than 2.25 as possibly significant

	j=0
        do 30 i=1,n-1
        if(zstar(i).ge.2.25) then
        j=j+1
        zbig(j)=zstar(i)               
	jbig(j)=nint(big(i))	       
        endif
30	continue     

        nbig = j                       

c--- obtain scale estimator using sum of gaps

        xsgap = sqrt(pi)/(n*(n-1))*ysum

        return
        end

c------------------------------------------------------------------------------
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
c------------------------------------------------------------------------------

c---  taken from Numerical Recipes

c******************************************************************************

      parameter (itmax=100,eps=3.e-7)
      GLN=GAMMLN(A)
      gold=0.
      a0=1.
      a1=x
      b0=0.
      b1=1.
      fac=1.
      do 11 n=1,itmax
        an=float(n)
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if(a1.ne.0.)then
          fac=1./a1
          g=b1*fac
          if(abs((g-gold)/g).lt.eps)go to 1
          gold=g
        endif
11    continue

      print *,'A too large, ITMAX too small'
1	exp1=-x+a*alog(x)-gln
	if (exp1.lt.-50.) exp1=-50.
      gammcf=exp(exp1)*g

      return
      end

c------------------------------------------------------------------------------
      SUBROUTINE GSER(GAMSER,A,X,GLN)
c------------------------------------------------------------------------------

c---  taken from Numerical Recipes

c******************************************************************************

      parameter (itmax=100,eps=3.e-7)
      GLN=GAMMLN(A)
      if(x.le.0.) then
        if(x.lt.0.) print *,'gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps)go to 1
11    continue
      print *,'a TOO LARGE, itmax TOO SMALL'
1     gamser=sum*exp(-x+a*log(x)-gln)

      return
      end

c------------------------------------------------------------------------------
        SUBROUTINE JACKNIFE  (XDATA,N,IP,X1LJACK,X1SJACK,
     +                       X1LLJACK,X1LSJACK,X2LJACK,
     +			     X2SJACK,X2LLJACK,X2LSJACK)
c------------------------------------------------------------------------------

c--- given an array XDATA of length N, and external function PARAM,
c    this routine calculates the expected error in the estimate
c    of the supplied function via the statistical jacknife.  In other
c    words, it obtains the spread in the parameter when calculated by
c    dropping one data point at a time.  For reference,
c    see the article by Diaconis and Efron in Scientific American.
c
c****************************************************************************

        implicit real*4 (a-h,o-z)
        real xdata(n),xout(50000),x1jparam(50000),x1p(50000),x1lp(50000)
        real x2jparam(50000),x2p(50000),x2lp(50000)

c--- find the parameter considering all values

	if (ip.eq.1 .or. ip.eq.2) then		
	x1all = PARAM(XDATA,N,1)
	x2all = PARAM(XDATA,N,2)

	else if (ip.eq.6) then
	call XBIWT(XDATA,N,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1)
        x1all = xlbiwt
	x2all = xsbiwt

	else if (ip.eq.7) then		
	x1all = PARAM(XDATA,N,IP)
	x2all = x1all
	endif

	if(x1all.gt.0)	then
	x1lall = alog10(x1all)
	else
	x1lall = 0
	endif
        
	if(x2all.gt.0) then
 	x2lall = alog10(x2all)
	else
	x2lall = 0
	endif

c--- create the pseudovalues data array with the missing value identified
c    by setting it to -32368., then calling PARAM, and repeating process

        do 10 i=1,n
          temp=xdata(i)                         
          xdata(i)=-32368.                      

	if (ip.eq.1 .or. ip.eq.2) then
	x1jparam(i)=PARAM(XDATA,N,1)
	x2jparam(i)=PARAM(XDATA,N,2)

	else if (ip.eq.6) then
	  CALL UPDATE(N,XDATA,M,XOUT)
	  CALL XBIWT(XOUT,M,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1)
          x1jparam(i) = xlbiwt			
	  x2jparam(i) = xsbiwt

	else if (ip.eq.7) then
	x1jparam(i)=PARAM(XDATA,N,IP)
	x2jparam(i)=x1jparam(i)
	endif

	if(x1jparam(i).gt.0) x1 = alog10(x1jparam(i))
	if(x2jparam(i).gt.0) x2 = alog10(x2jparam(i))

          x1p(i) = n*x1all - (n-1)*x1jparam(i)  
	  x2p(i) = n*x2all - (n-1)*x2jparam(i)
          x1lp(i) = n*x1lall - (n-1)*x1  
	  x2lp(i) = n*x2lall - (n-1)*x2
          xdata(i)=temp                  

10      continue

c--- obtain estimate of central location and spread

        x1sum = 0
	x2sum = 0
        x1xsum1 = 0
	x2xsum1 = 0
	x1lsum = 0
        x2lsum = 0
	x1xlsum1 = 0
        x2xlsum1 = 0
	x1xsum2 = 0
        x2xsum2 = 0
	x1xlsum2 = 0
        x2xlsum2 = 0

        do 20 i=1,n
  	x1sum = x1sum + x1p(i)
        x2sum = x2sum + x2p(i)
	x1xsum1 = x1xsum1 + x1p(i)*x1p(i)
        x2xsum1 = x2xsum1 + x2p(i)*x2p(i)
	x1lsum = x1lsum + x1lp(i)
        x2lsum = x2lsum + x2lp(i)
 	x1xlsum1 = x1xlsum1 + x1lp(i)*x1lp(i)
        x2xlsum1 = x2xlsum1 + x2lp(i)*x2lp(i)
20      continue

        x1xsum2 = x1sum*x1sum/n
	x2xsum2 = x2sum*x2sum/n
 	x1xlsum2 = x1lsum*x1lsum/n
        x2xlsum2 = x2lsum*x2lsum/n

        x1ljack = x1sum/n
	x2ljack = x2sum/n
	x1sjack = (abs(x1xsum1-x1xsum2)/(n*(n-1)))**0.5
        x2sjack = (abs(x2xsum1-x2xsum2)/(n*(n-1)))**0.5
        x1lljack = x1lsum/n
	x2lljack = x2lsum/n
	x1lsjack = (abs(x1xlsum1-x1xlsum2)/(n*(n-1)))**0.5
        x2lsjack = (abs(x2xlsum1-x2xlsum2)/(n*(n-1)))**0.5
	
        return
        end


c------------------------------------------------------------------------------
	SUBROUTINE KSREJ (DATA,N,KSAVE,KSDEV,NTOT)
c------------------------------------------------------------------------------

c--- This routine performs a "K-S" clipping of the data, similar to the
c    cleaning technique described in Yahil and Vidal (1977).  The procedure
c    first checks (with the K-S test) whether the sample is consistent with
c    a Gaussian distribution.  If not the procedure then deletes 
c    the data point most deviant from the sample mean, and redetermines the
c    sample mean and standard deviation.  Each point is tested in turn
c    in order of their apparent deviation, until no further rejection occurs.
c    The list is now considered "clean," and the final estimates of mean and
c    standard deviation are accepted.
c
c--- This routine will exit back to the main program if the initial number of
c    data points  N <= 5, or if the sample is clipped to NTOT=5 
c****************************************************************************

        implicit real*4 (a-h,o-z)
        dimension data(n),tdata(50000),sdata(50000),ssdata(50000)
        real ksave,ksdev

	if(n.le.5) then		
	   ksave=-999
	   ksdev=-999
	   ntot = n
	return
	endif

c--- copy the original data set

        do 5 i=1,n
        sdata(i) = data(i)
5       continue

c--- find intial estimate of sample mean and deviation

        ave  = PARAM(SDATA,N,1)              
        sdev = PARAM(SDATA,N,2)              
        nout=0

c--- test once to see if all is OK

	ITEST=2  

c	CALL EDFGOF(N,XDATA,ITEST,XBAR,S2,D,V,W2,W2P,U2,U2P,
	CALL EDFGOF(N,DATA,ITEST,XBAR,S2,D,V,W2,W2P,U2,U2P,
     +  A2,A2P,IFAULT)

	if(d.le..775) test=.25
	if(d.gt..775.and.d.le..819)  test=.15
	if(d.gt..819.and.d.le..895)  test=.10
	if(d.gt..895.and.d.le..995)  test=.05
	if(d.gt..995.and.d.le.1.035) test=.025
	if(d.gt.1.035) test=.01

	if (test.gt.0.20) then		
	ksave = ave
	ksdev = sdev
	ntot = n

	return
	endif

c--- create a list of data ordered by their distance from the sample mean

20      continue

        do 10 i=1,n-nout
          tdata(i) = abs(sdata(i)-ave)
10      continue

        CALL SORT2 (N-NOUT,TDATA,SDATA)              

c--- drop first deviant point

          temp = sdata(n-nout)
          sdata(n-nout) = -32368.
          ave  = PARAM(SDATA,N-NOUT,1)             
          sdev = PARAM(SDATA,N-NOUT,2)             

c    test if new mean and standard deviation give an acceptable Gaussian fit
c    if so, then delete offending point from consideration

	CALL UPDATE (N-NOUT,SDATA,M,SSDATA)
	CALL SORT (M,SSDATA)

	if(m.eq.5) then                 
	  ksave=-999
	  ksdev=-999
	  ntot=5
	  return
	endif

c--- test to see if consistent with a Gaussian
	
	ITEST=2 

	CALL EDFGOF(M,SSDATA,ITEST,XBAR,S2,D,V,W2,W2P,U2,U2P,
     +  A2,A2P,IFAULT)

	if(d.le..775) test=.25
	if(d.gt..775.and.d.le..819)  test=.15
	if(d.gt..819.and.d.le..895)  test=.10
	if(d.gt..895.and.d.le..995)  test=.05
	if(d.gt..995.and.d.le.1.035) test=.025
	if(d.gt.1.035) test=.01

        if(test.le.0.20) then 
          nout = nout + 1
        else
          ksave=ave
          ksdev=sdev
          ntot=n-nout-1
          return
        endif
        goto 20

        end

c------------------------------------------------------------------------------
	SUBROUTINE KURTSKEW (N,SKEW,KURT,Z1,Z2,K2)
c------------------------------------------------------------------------------

c---	This subroutine obtains the statistics z1 and z2, which are
c	normally distributed under the hypothesis of normality of
c	the underlying population.  It also returns the statistic k2, 
c	which is distributed like Chi^2 when the underlying population
c	is normally distributed.
c
c	The user inputs the sample skewness, SKEW, and the sample kurtosis,
c	KURT.  For details of the calculation see  D'Agostino, Belanger,
c	and D'Agostino	(1990), American Statistician 44, 316.
c
c	For values of n<9 this should not be used
c
c*****************************************************************************

	implicit real*4 (a-h,o-z)
	real*4 k2,kurt,nn
	integer n

	nn = n

	if (nn.lt.9) return	
	
c---	This is the skewness test

	y = skew*((nn+1)*(nn+3)/(6*(nn-2)))**0.5
	beta2 = 3*(nn*nn+27*nn-70)*(nn+1)*(nn+3)/((nn-2)*(nn+5)
     +          *(nn+7)*(nn+9))
	w = (-1.+(2.*(beta2-1.))**0.5)**0.5
	d = 1./(log(w))**0.5
	a = (2./(w**2-1.))**0.5
	z1 = d*log(y/a+((y/a)**2+1.)**0.5)

c---	This is the kurtosis test

	e = 3*(nn-1)/(nn+1)
	var =24*nn*(nn-2)*(nn-3)/((nn+1)**2*(nn+3)*(nn+5))
	x = (kurt-e)/(var)**0.5
	beta1 = (6*(nn**2-5*nn+2)/((nn+7)*(nn+9))*(6*
     +          (nn+3)*(nn+5)/(nn*(nn-2)*(nn-3)))**0.5)**2
	A = 6.+8./(beta1**0.5)*(2./(beta1**0.5)+(1.+4./beta1)**0.5)
	exp1=(1.-(2./A))/(1.+x*(2./(A-4.))**0.5)
	exp2=abs(exp1)**0.333
	if (exp1.lt.0.) exp2=-exp2
	z2 = (1.-(2./9./A)-exp2)/
     +       (2./(9.*A))**0.5

c---	This is the omnibus test

	k2 = (z1)**2+(z2)**2

	return
	end

c------------------------------------------------------------------------------
        SUBROUTINE LETTERS (XDATA,N,DEPTHS,XLV,MLV,SLV,RLV,CLV,NLV)
c------------------------------------------------------------------------------

c--- NOTE -- this routine also sorts the data low to high 
c--- For the batch of values in XDATA, find the selected quantiles known
c    as the letter values.  Upon exit, XLV contains the letter values,
c    D contains the corresponding depths, and NLV is the number of pairs
c    of letter values.  Specifically, XLV(1,1) and XLV(1,2) are both
c    set equal to the median, whose depth DEPTHS(1), is (N+1)/2.  The rest of
c    the letter values come in pairs and are stored in XLV in order from
c    the hinges out to the extremes.  Thus XLV(2,1) and XLV(2,2) are the
c    lower hinge and upper hinge, respectively, and XLV(NLV,1) and XLV(NLV,2)
c    are the lower extreme (minimum) and upper extremem (maximum), respectively.

c
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        integer n,nlv,i,j,k,pt1,pt2
        dimension xdata(n),depths(15),xlv(15,2),slv(15),glv(8),rlv(8)
	real mlv(8),clv(8),cclv(50000,8)
        data glv/1.,1.349,2.301,3.068,3.726,4.308,4.836,5.320/
 
c--- read in correction factors

	if (n.gt.100) then
	xn=float(n)
	cclv(n,1)=1.-(2./xn)*(1.-.86*((xn/2.)-(n/2)))
	cclv(n,2)=1.71-(3./xn)*(1.-.8*((xn/4.)-(n/4)))
	cclv(n,3)=2.27-(5./xn)*(1.-.77*((xn/8.)-(n/8)))
	cclv(n,4)=2.76-(10./xn)*(1.-.76*((xn/16.)-(n/16)))
	cclv(n,5)=3.19-(15./xn)*(1.-.76*((xn/32.)-(n/32)))
	cclv(n,6)=3.58-(30./xn)*(1.-.75*((xn/64.)-(n/64)))
	cclv(n,7)=3.94-(55./xn)*(1.-.75*((xn/128.)-(n/128)))
	do 10 i=1,7
	cclv(n,i)=1.349*cclv(n,i)
10	continue
	cclv(n,8)=1.0
	else

	do 11 i=1,100
	read (8,*) id,(cclv(i,j),j=1,7)	
	cclv(i,8) = 1.0	
11	continue

	rewind (unit=8)
	endif

c--- sort the data low to high
 
        CALL SORT (N,XDATA)
 
c--- handle median separately because it is not a pair of lv

        depths(1) = float(n+1)/2.0
        j = (n/2) + 1
        pt2 = n + 1 - j
        xlv(1,1) = (xdata(j) + xdata(pt2))/2.0
        xlv(1,2) = xlv(1,1)
 
        k=n
        i=2
 
20      k= (k+1)/2
        j = (k/2) + 1
        depths(i) = float(k+1)/2.0
        pt2 = k + 1 - j
        xlv(i,1) = (xdata(j) + xdata(pt2))/2.0
        pt1 = n-k+j
        pt2 = n+1-j
        xlv(i,2) = (xdata(pt1) + xdata(pt2))/2.0
 
        i = i + 1
        if(depths(i-1) .gt. 2.0) goto 20

        nlv = i
	if(nlv.gt.8) nlv=8
        depths(i) = 1.0
        xlv(i,1) = xdata(1)
        xlv(i,2) = xdata(n)
 
c--- obtain list of letter mids,spreads, and ratios to Gaussian values
 
        slv(1)=0
        rlv(1)=0        
        do 30 i=2,min(8,nlv)
          slv(i)= abs(xlv(i,1)-xlv(i,2))
          rlv(i) = slv(i)/glv(i)
	  mlv(i) = (xlv(i,1)+xlv(i,2))/2.
30      continue

c--- obtain corrections to spreads as a function of n

	do 40 i=2,min(8,nlv)
	  clv(i) = slv(i)/cclv(n,i-1)		
40	continue

        return
        end     


c------------------------------------------------------------------------------
        SUBROUTINE MCCON (X,N,CONF,XAVE,XSDEV,XLBIWT,XSBIWT,
     +                    MCCONF,NSIMS)
c------------------------------------------------------------------------------
 
c--- This routine obtains the BOOTSTRAP confidence
c    intervals for an input parameter
 
c*****************************************************************************
 
        real msim(50000),ssim(50000),bsim(50000),sbsim(50000)
        real conf(4)
        real x(n),xout(50000),mcconf(32,4)
	real infav(50000),infsd(50000),infbt(50000),infsbt(50000)
        logical lsb,lbb
 
c--- obtain sampling distribution
 
        CALL BOOTSTRAP (X,N,1,NSIMS,MSIM,SSIM)

	CALL BOOTSTRAP (X,N,6,NSIMS,BSIM,SBSIM) 
c
c--    calculate the acceleration constants for the data set
c      using the Jacknife Estimation of Frangos and Schucany
c      Comp. Stat. & Data Anal. 9 (1990) 271-281
c
        do 15 i=1,n
          temp=x(i)                         
          x(i)=-32368.  

	  CALL UPDATE(N,X,M,XOUT)
	  xnave = PARAM(XOUT,M,1)
          xnsdev = PARAM(XOUT,M,2)
	  CALL XBIWT(XOUT,M,XNLBIWT,XNSBIWT,XNLBIWT1,XNSBIWT1)

	infav(i) = (n-1)*(xave-xnave)
	infsd(i) = (n-1)*(xsdev-xnsdev)
	infbt(i) = (n-1)*(xlbiwt-xnlbiwt)
	infsbt(i) = (n-1)*(xsbiwt-xnsbiwt)
	x(i) = temp			    
15	continue

	snum1=0
	sden1=0
	snum2=0
	sden2=0
	snum3=0
	sden3=0
	snum4=0
	sden4=0

	do 16 i=1,n
	snum1= snum1 + infav(i)**3
	sden1= sden1 + infav(i)**2
	snum2= snum2 + infsd(i)**3
	sden2= sden2 + infsd(i)**2
	snum3= snum3 + infbt(i)**3
	sden3= sden3 + infbt(i)**2
	snum4= snum4 + infsbt(i)**3
	sden4= sden4 + infsbt(i)**2
16	continue
	aave=snum1/(6*(sden1**1.5))
	asdev=snum2/(6*(sden2**1.5))
	abwt=snum3/(6*(sden3**1.5))
	asbwt=snum4/(6*(sden4**1.5))


c--- obtain monte-carlo estimators of intervals

c--- get confidence iunterval of the mean from MSIM

	do 10 i=1,4
 
        lsb=.true.
        sbound=0.0
        lbb=.true.
        bbound=0.0
        ifault = 0
 
        CALL MONTE (XAVE,CONF(i),NSIMS,AAVE,MSIM,LSB,
     +  SBOUND,LBB,BBOUND,MCCONF(1,i),MCCONF(2,i),MCCONF(3,i),
     +  MCCONF(4,i),MCCONF(5,i),MCCONF(6,i),MCCONF(7,i),
     +  MCCONF(8,i),IFAULT)

C--  get biweight conf int on location from BSIM

        lsb=.true.
        sbound=0.0
        lbb=.true.
        bbound=0.0
        ifault = 0
 
        CALL MONTE (XLBIWT,CONF(i),NSIMS,ABWT,BSIM,LSB,
     +  SBOUND,LBB,BBOUND,MCCONF(9,i),MCCONF(10,i),MCCONF(11,i),
     +  MCCONF(12,i),MCCONF(13,i),MCCONF(14,i),MCCONF(15,i),
     +  MCCONF(16,i),IFAULT)

c--- get confidence interval of the stand dev from SSIM
 
        lsb=.true.
        sbound=0.0
        lbb=.true.
        bbound=0.0
        ifault = 0

        CALL MONTE (XSDEV,CONF(i),NSIMS,ASDEV,SSIM,LSB,
     +  SBOUND,LBB,BBOUND,MCCONF(17,i),MCCONF(18,i),MCCONF(19,i),
     +  MCCONF(20,i),MCCONF(21,i),MCCONF(22,i),MCCONF(23,i),
     +  MCCONF(24,i),IFAULT)

c--  get biweight spread conf limits from SBSIM

        lsb=.true.
        sbound=0.0
        lbb=.true.
        bbound=0.0
        ifault = 0

        CALL MONTE (XSBIWT,CONF(i),NSIMS,ASBWT,SBSIM,LSB,
     +  SBOUND,LBB,BBOUND,MCCONF(25,i),MCCONF(26,i),MCCONF(27,i),
     +  MCCONF(28,i),MCCONF(29,i),MCCONF(30,i),MCCONF(31,i),
     +  MCCONF(32,i),IFAULT)

10	continue

        return
        end
 

c------------------------------------------------------------------------------
	SUBROUTINE MDIAN1 (X,N,XMED)
c------------------------------------------------------------------------------

c---  Taken from Numerical Recipes, page 460.
c     Given an array X of N numbers, returns their median value XMED, and 
c     the array is sorted low to high
c
c*******************************************************************************

	real x(n)
	call sort(n,x)
	n2=n/2
	if (2*n2.eq.n) then
		xmed=0.5*(x(n2)+x(n2+1))
	else
		xmed=x(n2+1)
	endif
	return
	end


c------------------------------------------------------------------------------
        SUBROUTINE MOMENT(DATA,N,AVE,ADEV,SDEV,VAR,SKEW,CURT)
c------------------------------------------------------------------------------

c--- given an array DATA of length N, this routine returns its
c    mean AVE, average deviation ADEV, standard deviation SDEV,
c    variance VAR, skewness SKEW, and kurtosis CURT
c 
c--- stolen (unabashedly) from NUMERICAL RECIPES
c
c**************************************************************************** 
 
      implicit real*4 (a-h,o-z)
      dimension data(n)

      if(n.le.1) print *,'n MUST BE AT LEAST 2'
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      do 12 j=1,n
        s=data(j)-ave
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=var/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/float(n)/(var*(n-1)/float(n))**1.5
        curt=curt/float(n)/(var*(n-1)/float(n))**2  
      else
	skew = 0.0
	curt = 0.0
      endif
      return
      end

c------------------------------------------------------------------------------
        SUBROUTINE MONTE(EST,CONFI,NSIMS,ACC,SIMVAL,LSB,
     +                   SBOUND,LBB,BBOUND,ALOW,AHI,BLOW,BHI,
     +                   CLOW,CHI,DLOW,DHI,IFAULT)
c------------------------------------------------------------------------------

c
c--- Algorithm AS 214 Appl. Stat. (1985) Vol. 34, No.3
c   
c--- Sets up Monte Carlo confidence intervals
c    uses function PPND -  AS 111
c    uses function ALNORN -- AS 66
c
c--- The user supplies the estimate EST of a parameter, and MONTE
c    obtains the appropriate confidence intervals
c
c               ALOW: lower conf. limit, symmetric
c               AHI:  upper conf. limit, symmetric
c               BLOW: lower conf. limit, percentile 
c               BHI:  upper conf. limit, percentile
c               CLOW: lower conf. limit, bias corrected
c               CHI:  upper conf. limit, bias corrected
c               DLOW: lower conf. limit, bias corrected accelerated
c               DHI:  upper conf. limit, bias corrected accelerated
c
c***************************************************************************
 
        logical lsb,lbb
        dimension simval(nsims)
        data least/5/
        data zero,half,one,two,hun,big
     +       /0.0E0,0.5E0,1.0E0,2.0E0,1.0E2,1.e20/
        ifault = 0
        sbnd = sbound
        bbnd = bbound
        if (confi .ge. hun) ifault = 6
        if (confi .le. zero) ifault = 7
        if (ifault .gt. 0) return
c
c--- symmetric mci
c
c--- find mean and variance
c
        v1 = zero
        v2 = zero
        do 10 j = 1,nsims
        fj = float(j)
        v2 = v2 + (fj - one)*(simval(j) - v1)**2/fj
        v1 = (simval(j) + (fj - one)*v1)/fj
 10     continue
 
        stderr = zero
        if (v2 .gt. zero) stderr = sqrt(v2/float(nsims - 1))
        alpha = half*(hun - confi)/hun
        z = -ppnd(alpha,ifault)
        if (ifault .eq. 0) goto 20
        ifault = 6
        return
 20     ahi = z*stderr
        alow = est - ahi
        ahi = est + ahi
c
c--- calculate bias adjustment, so that the number
c    of values that must be ordered is known
c
        call biasad(est,nsims,acc,simval,z,limit1,limit2,
     +              limita1,limita2,ifault)
	
        limitl = alpha*float(nsims + 1) + half
        limitu = (one - alpha)*float(nsims + 1) + half
	if (limitl .lt. least) ifault = 5
        if (limitl .eq. limitu) ifault = 7
        if (ifault .gt. 0) return
        l1 = max0(limit1, 2*limitl)
        l2 = min0(limit2, 2*limitu - nsims - 1)
	if (lsb) sbnd = -big
        if (lbb) bbnd = big
c
c--- select and order the l1 smallest and the nsims+1-l2 largest
c    values
c

	call sort(nsims,simval)
        if (sbnd .le. simval(1)) goto 40
        if (lsb) goto 30
        ifault = 1
        return
 30     sbnd = simval(1)*two
        if (sbnd .gt. zero) sbnd = -sbnd
 40     if (bbnd .ge. simval(nsims)) goto 60
        if (lbb) goto 50
        ifault = 2
        return
 50     bbnd = simval(nsims)*two
        if (bbnd .lt. zero) bbnd = -bbnd
c
c--- equal tails mci (percentile method)
c
 60	if (limitl .gt. 0) blow = simval(limitl)
	if (limitl .le. 0) blow = simval(1)
	if (limitu .le. nsims) bhi = simval(limitu)
	if (limitu .gt. nsims) bhi = simval(nsims)
c
c--- bias-corrected percentile method
c
        clow = sbnd
        if (limit1 .gt. 0) clow = simval(limit1)
	if (limit1.le.0) clow = simval(1)
        chi = bbnd
        if (limit2 .le. nsims) chi = simval(limit2)
	if (limit2.gt.nsims) chi = simval(nsims)

c--  accelerated bias corrected

        dlow = sbnd
        if (limita1 .gt. 0) dlow = simval(limita1)
	if (limita1.le.0) dlow = simval(1)
        dhi = bbnd
        if (limita2 .le. nsims) dhi = simval(limita2)
	if (limita2.gt.nsims) dhi = simval(nsims)
	
        return
        end


c------------------------------------------------------------------------------
        SUBROUTINE NORMALITY (XDATA,N,XAVE,XSDEV,XSKEW,XCURT,
     +                        XSBIWT,XMED,TEST)
c------------------------------------------------------------------------------

c--- This routine obtains the tests for non-normality of a data sample. 
c    A number of tests are performed.  In order, these are:
c
c    (1)  The a- , u-, w- tests referred to by Yahil
c         and Vidal (1977) (Ap.J. 214, 347) with an improved W-test
c         calculation due to Royston (A.S. 181) -- critical values for
c         a and u in REFS, critical value for w calculated and returned
c
c    (2)  The b1 and b2 tests (essentially a renormalized set of third and
c         fourth moments) -- probability given in D'Agostino, Belanger,
c         D'Agostino (1990), The American Statistician, 44, 316-321
c
c    (3)  The i-test described by Iglewicz (1983) (UREDA, p.426) -- critical
c         values are calculated and returned
c
c    (4)  The DIP statistic test for unimodality (Hartigan A.S. 217) --
c         critical values given in Hartigan and Hartigan (1985) Ann. Stat.
c         13, 70.
c
c    (5)  A modified K-S statistic -- see Goodness of Fit (D-Agistino and
c         Stephens, Chp. 4) 
c
c    (6)  A modified Kuiper statistic -- see Goodness of Fit (D-Agistino and
c         Stephens, Chp. 4) 
c
c    (7) A modified Cramer Von-Mises -- see Goodness of Fit (D-Agistino and
c         Stephens, Chp. 4) 
c
c    (8) A modified Watson test -- see Goodness of Fit (D-Agistino and
c         Stephens, Chp. 4) 
c
c    (9) A modified Anderson-Darling test --see Goodness of Fit (D-Agistino and
c         Stephens, Chp. 4) 
c
c****************************************************************************

        implicit real*4 (a-h,o-z)
        dimension xdata(n),a(50000)
        real test(26),i90,i95,istat,isum,ntilda,k2
        integer mn(50000),mj(50000),gcm(50000),lcm(50000)

	if (n.lt.5) return

	do 9 i=1,26
9	test(i) = 0.0

c----------------------  TEST 1 -- A, U, W ---------------------

        a1sum=0.0
        do 10 i=1,n                             
          a1sum = a1sum + abs(xdata(i)-xave)
10      continue

        a1 = a1sum/(n*xsdev)
        test(1) = a1

        u1 = (xdata(n)-xdata (1))/xsdev         
        test(2) = u1

c--- a slightly different version of the W-test

        n2=n/2
        ifault = 0
	ssq = 0.

c--- obtain SSQ

        do 19 i=1,n
19      ssq = ssq+(xdata(i)-xave)**2

        CALL WCOEF (A,N,N2,EPS,IFAULT)
        CALL WEXT (XDATA,N,SSQ,A,N2,EPS,W,PW,IFAULT)

	test(3) = w
        test(4) = pw

c----------------------  TEST 2 -- b1,b2 TESTS -----------------------

	if (n.lt.9) then		
	test(23)=-999
	test(24)=-999
	test(25)=-999
	test(26)=-999
	goto 323
	endif

	CALL KURTSKEW (N,XSKEW,XCURT,ZB1,ZB2,K2)

	test(23)=ALNORM(ZB1,.true.)	
	if (xskew.lt.0.) test(23)=1.-test(23)	
	test(24)=ALNORM(ZB2,.true.)		
	if (xcurt.lt.3.) test(24)=1.-test(24)	
	test(25)=K2
	test(26)=GAMMQ(1.,0.5*K2)		
323     test(5) = xskew
        test(6) = xcurt

c----------------------  TEST 3 -- I-TEST ----------------------------

        isum = 0.0
        do 100 i=1,n
          isum = isum + (xdata(i)-xmed)**2
100     continue

        istat = isum/((n-1)*xsbiwt*xsbiwt)
        test(7) = istat

C--- use an approximate formula for the percentage points

        ntilda = alog10(float(n-1))
        i90 = 0.6376-1.1535*ntilda+0.1266*ntilda*ntilda
        if (n.lt.50) then
          i95 = 1.9065-2.5465*ntilda+0.5652*ntilda*ntilda
        else
          i95 = 0.7824-1.1021*ntilda+0.1021*ntilda*ntilda
        endif

        test(8) = 10**i90+0.982
        test(9) = 10**i95+0.982

c----------------------  TEST 4 -- DIP TEST ----------------------------

 	if(n.ge.7) CALL DIPTST (XDATA,N,DIP,XL,XU,IFAULT,GCM,LCM,MN,MJ)
        test(10) = dip
        test(11) = xl
        test(12) = xu

c----------------------  TEST 5 -- MODIFIED KS-TEST--------------------

c
c-- Use the algorithm of Davis and Stephens
c

	ITEST=2  

	CALL EDFGOF(N,XDATA,ITEST,XBAR,S2,D,V,W2,W2P,U2,U2P,
     +  A2,A2P,IFAULT)
	test(13) = d

	if(d.le..775) test(14)=.25
	if(d.gt..775.and.d.le..819)  test(14)=.15
	if(d.gt..819.and.d.le..895)  test(14)=.10
	if(d.gt..895.and.d.le..995)  test(14)=.05
	if(d.gt..995.and.d.le.1.035) test(14)=.025
	if(d.gt.1.035) test(14)=.01

c---------------------- TEST 6 -- MODIFIED KUIPER STATISTIC ---------

	test(15) = v

	if(v.le.1.32) test(16)=.25
	if(v.gt.1.32.and.v.le.1.386)  test(16)=.15
	if(v.gt.1.386.and.v.le.1.489)  test(16)=.10
	if(v.gt.1.489.and.v.le.1.585)  test(16)=.05
	if(v.gt.1.585.and.v.le.1.693) test(16)=.025
	if(v.gt.1.693) test(16)=.01

c-------------- TEST 7 -- MODIFIED CRAMER-VON-MISES STATISTIC -------

	test(17) = W2
	test(18) = W2P 

c----------------------- TEST 8 -- MODIFIED WATSON STATISTIC ---------

	test(19) = U2
	test(20) = U2P 

c--------------- TEST 9 -- MODIFIED ANDERSON-DARLING STATISTIC--------

	test(21) = A2
	test(22) = A2P 

        return
        end


c----------------------------------------------------------------------------- 
        SUBROUTINE NSCOR2(S,N,N2,IFAULT)
c----------------------------------------------------------------------------- 

c--- Algorithm AS 177.3  Appl. Stat. (1982) Vol. 31, No.2
c
c--- Approximation for rankits.
c
        real s(n2),eps(4),dl1(4),dl2(4),gam(4),lam(4),bb,d,
     +       b1,an,ai,e1,e2,l1,correc,ppnd
        data eps(1),eps(2),eps(3),eps(4)
     +       /.419885e0,.450536e0,.456936e0,.468488e0/,
     +       dl1(1),dl1(2),dl1(3),dl1(4)
     +       /.112063e0,.121770e0,.239299e0,.215159e0/,
     +       dl2(1),dl2(2),dl2(3),dl2(4)
     +       /.080122e0,.111348e0,-.211867e0,-.115049e0/,
     +       gam(1),gam(2),gam(3),gam(4)
     +       /.474798e0,.469051e0,.208597e0,.259784e0/,
     +       lam(1),lam(2),lam(3),lam(4)
     +       /.282765e0,.304856e0,.407708e0,.414093e0/,
     +       bb/-.283833e0/,d/-.106136e0/,b1/.5641896e0/
        ifault = 3
        if (n2 .ne. n/2) return
        ifault = 1
        if (n .le. 1) return
        ifault = 0
        if (n .gt. 2000) ifault = 2
        s(1) = b1
        if  (n .eq. 2) return
c
c--- calculate normal areas for 3 largest rankits
c
        an = n
        k = 3
        if (n2 .lt. k) k = n2
        do 5 i = 1,k
        ai = i
        e1 = (ai - eps(i))/(an + gam(i))
        e2 = e1**lam(i)
        s(i) = e1 + e2*(dl1(i) + e2*dl2(i))/ an - correc(i,n)
 5      continue
        if (n2 .eq. k) goto 20
c
c--- calculate normal areas for remaining rankits
c
        do 10 i = 4,n2
        ai = i
        l1 = lam(4) + bb/(ai + d)
        e1 = (ai - eps(4))/ (an + gam(4))
        e2 = e1**l1
        s(i) = e1 + e2*(dl1(4) + e2*dl2(4))/ an - correc(i,n)
 10     continue
c
c--- convert normal tail areas to normal deviates
c
 20     do 30 i = 1,n2
 30     s(i) = -ppnd(s(i),ifault)
        return
        end


c----------------------------------------------------------------------------- 
	SUBROUTINE PRINTER (OUNIT,FILENAME,R)
c----------------------------------------------------------------------------- 

c--- This routine writes out the results of the statistical routine
c    ROSTAT in a format which is comprehensible by the average human

c******************************************************************************

	integer ounit
	real*4 r(50000*2)
	character*32 filename
	common/print/ jtnor,jtgap,jtcon,jtchi,jtasy,jtcorr,jtot

 3      format(/)
 2	format(//)
 1	format('1')

c---	Declare the number of decimals you want to output

	if(abs(r(2)).gt.50.) id=1
	if(abs(r(2)).le.50.) id=3

c--- write out sorted dat

	write (ounit,1)
	write (ounit,5) filename
 5	format(t25,'FILE  =  ',a)

	write (ounit,2)
	
	n = ifix(r(1))
	write (ounit,6) n
 6	format(t25,'SORTED DATA, N = ',i3)
	write (ounit,7)
 7	format (t25,'====================')
	write (ounit,3)
	write (ounit,10) (r(j),j=2,n+1)
10	format(5(1x,f10.2))

c--- letter value matrix

	write (ounit,2)
	write (ounit,2)
	write (ounit,11) 
11	format(t25,'LETTER VALUE MATRIX')
	write (ounit,111)
111	format(t25,'===================')
	write (ounit,2)

	write (ounit,12)
12	format(t10,'DEPTH',T20,'LOWER',T30,'UPPER',T40,'MID',T50,'SPREAD',
     +         T60,'RATIO',t70,'PSEUDO')
	write (ounit,3)

	nl = ifix(r(n+2))
	do 14 i=1,nl
	k = n+3+(i-1)*7
	write (ounit,13) i,r(k),r(k+1),r(k+2),r(k+3),r(k+4),r(k+5),r(k+6)
14	continue

13	format(1x,i2,t5,f10.2,t15,f10.2,t25,f10.2,t35,f10.2,
     +         t45,f10.2,t55,f10.2,t65,f10.2)


c--- classical estimators

	write (ounit,1)
	write (ounit,5) filename
	write (ounit,2)
	write (ounit,15) 
15	format(t25,'CLASSICAL ESTIMATORS')
	write (ounit,115)
115	format(t25,'====================')
	write (ounit,2)

	k=n+7*nl
	write (ounit,16) 
16	format(t10,'MEAN',t20,'SDEV',T30,'ADEV',T40,'SKEW',T50,
     +         'CURT',T60,'VAR')
	WRITE (OUNIT,3)
	write (ounit,17) r(k+3),r(k+5),r(k+4),r(k+7),r(k+8),r(k+6)
17	format(t5,f10.2,t15,f10.2,t25,f10.2,t35,f10.2,t45,f10.2,
     +         t57,E10.4)

c--- location estimates

	write (ounit,2)
	write (ounit,20)
20	format(t25,'ESTIMATES OF CENTRAL LOCATION')
	write (ounit,201)
201	format(t25,'=============================')
	write (ounit,2)

	write (ounit,21)
21	format (t10,'MEAN',t20,'3S-MEAN',T30,'N-3SIG',T40,'KS-MEAN',T50,
     +          'N-KS  ',T60,'MEDIAN')
	write (ounit,3)
	write(ounit,22) r(k+9),r(k+11),ifix(r(k+10)),r(k+13),
     +  ifix(r(k+12)),r(k+20)

22	format (t5,f10.2,t15,f10.2,t25,i10,t35,f10.2,
     +  t45,i10,t55,f10.2)

	write (ounit,2)
	write (ounit,211)
211	format(T8,'TRIM-5',t18,'TRIM-10',T28,'TRIM-20',T38,'MID-MEAN',
     +           T48,'TRIM-30',T58,'TRIM-40')
	write (ounit,3) 
	write(ounit,221) r(k+14),r(k+15),r(k+16),r(k+17),r(k+18),r(k+19)
221	format (t5,f10.2,t15,f10.2,t25,f10.2,t35,f10.2,
     +          t45,f10.2,t55,f10.2)

	write (ounit,2)

	write (ounit,23) 
23	format (t10,'BMED',T20,'TRIMEAN',T30,'BIWT',T40,'BIWT-1')
	write (ounit,3)
	write (ounit,24) r(k+21),r(k+22),r(k+23),r(k+24)

24 	format (t5,f10.2,t15,f10.2,t25,f10.2,t35,f10.2)

	if (r(k+60).eq.0.) goto 2000
	write (ounit,2)
	write (ounit,*) '  The above central location estimates are '
	write (ounit,*) '  with respect to the biweight estimate = ',
     +                  r(k+59)
2000	continue

c--- scale estimates

	write (ounit,1)
	write (ounit,5) filename
	write (ounit,2)
	write (ounit,25)
25	format(t25,'ESTIMATES OF DATA SCALE')
	write (ounit,251)
251	format(t25,'=======================')
	write (ounit,2)

	write (ounit,26)
26	format(t10,'SDEV',T20,'3S-DEV',T30, 'N-3SIG',T40,'KS-DEV',
     +         T50,'N-KS')
	write (ounit,3)

	write (ounit,27) r(k+30),r(k+32),ifix(r(k+31)),
     +                   r(k+34),ifix(r(k+33))
27	format(t5,f10.2,t15,f10.2,t25,i10,t35,f10.2,t45,i10)
	write (ounit,2)

	write (ounit,28)
28	format(T10,'SIG-F',T20,'SIG-MAD',t30,'SIG-GAP',t40,'S-BIWT',
     +  T50,'S-BIWT1')
	write (ounit,3)
	write (ounit,29) r(k+38),r(k+40),r(k+58),r(k+41),r(k+42)
29	format (t5,f10.2,t15,f10.2,t25,f10.2,t35,f10.2,t45,f10.2)

c--- confidence intervals

	write (ounit,1)
	write (ounit,5) filename
	write (ounit,3)

	if (r(jtcon+165).eq.1) write (ounit,50)
	if (r(jtcon+165).eq.0) write (ounit,502)
50	format(t5,'ERRORS IN LOCATION: FORMULAE(XMSC IS 3-SIGMA CLIPPED)')
502	format(t5,'ERRORS IN LOCATION: FORMULAE(  XMSC IS NOT CLIPPED  )')

	write (ounit,501)
501	format(t5,'====================================================')

	write (ounit,3)
	write (ounit,51)
51	format(t20,'LOCATION',T30,'PLU/MIN',T40,'LOWER',T50,'UPPER')

	write(ounit,52)
52	format(t5,'XMSC')	

	if (r(jtcon+165).eq.1.) k=n+7*nl+11	
	if (r(jtcon+165).eq.0.) k=n+7*nl+9	

	do 55 i=1,4
	l=jtcon+i
	write(ounit,53) r(k),r(l),r(k)-r(l),r(k)+r(l)
55	continue
53	format(t15,f10.2,t25,f10.2,t35,f10.2,t45,f10.2)

	write(ounit,54)
54	format(t5,'XMDC')	

	k=n+7*nl+20
	do 56 i=5,8
	l=jtcon+i
	write(ounit,53) r(k),r(l),r(k)-r(l),r(k)+r(l)
56	continue

	write(ounit,57)
57	format(t5,'XTSC')	

	k=n+7*nl+23
	do 59 i=9,12
	l=jtcon+i
	write(ounit,53) r(k),r(l),r(k)-r(l),r(k)+r(l)
59	continue

	write(ounit,3)
	write(ounit,80)
80	format(t20,'ERRORS IN LOCATION:  BOOTSTRAP')
	write(ounit,801)
801	format(t20,'==============================')

	write (ounit,3)
	write(ounit,81)
81	format(t20,'LOCATION',T30,'LENGTH',T40,'LOWER',T50,'UPPER')

	write(ounit,82)
82	format(t5,'MEAN-1')

	k=n+7*nl+3
	do 85 i=37,44,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
85	continue

	write(ounit,86)
86	format(t5,'MEAN-2')

	k=n+7*nl+3
	do 87 i=45,52,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
87	continue

	write(ounit,88)
88	format(t5,'MEAN-3')

	k=n+7*nl+3
	do 89 i=53,60,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
89	continue

	write(ounit,90)
90	format(t5,'MEAN-4')

	k=n+7*nl+3
	do 91 i=133,140,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
91	continue

	write(ounit,810)
810	format(t5,'BIWT-1')	

	k=n+7*nl+23
	do 811 i=61,68,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
811	continue

	write(ounit,812)
812	format(t5,'BIWT-2')	

	k=n+7*nl+23
	do 813 i=69,76,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
813	continue

	write(ounit,814)
814	format(t5,'BIWT-3')	

	k=n+7*nl+23
	do 815 i=77,84,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
815	continue

	write(ounit,92)
92	format(t5,'BIWT-4')	

	k=n+7*nl+23
	do 93 i=141,148,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
93	continue

	write (ounit,1)
	write (ounit,5) filename
	write (ounit,3)

	if(r(jtcon+165).eq.1 ) write (ounit,65) 
	if(r(jtcon+165).eq.0 ) write (ounit,655)
65	format(t10,'ERRORS IN SCALE:  CLASSICAL (3-SIGMA CLIPPED)')
655	format(t10,'ERRORS IN SCALE:  CLASSICAL (  NO CLIPPING  )')
 
	write (ounit,651)
651	format(t10,'=============================================')
	write (ounit,3)
	write (ounit,66)
66	format(t15,'SCALE',T25,'LO-ERROR',T35,'HI-ERROR',T45,'LOWER',
     +         T55,'UPPER')

	write (ounit,3)
	write (ounit,642)
642	format(t5,'XXSC')

	if (r(jtcon+165).eq.1.) k=n+7*nl+32	
	if (r(jtcon+165).eq.0.) k=n+7*nl+30	
	do 67 i=1,4
	l=jtchi+1+(i-1)*2
	write (ounit,68) r(k),r(k)-r(l),r(l+1)-r(k),r(l),r(l+1)
67	continue

68	format(t10,f10.2,t20,f10.2,t30,f10.2,t40,f10.2,t50,f10.2)

c	write (ounit,641)
c641	format(t5,'XJLC')

c	k=n+7*nl+45
c	do 64 i=17,20
c	l=jtcon+i
c	rv = 10**r(k)
c	rlv = 10**(r(k)-r(l))
c	ruv = 10**(r(k)+r(l))
c	rle = rv-rlv
c	rue = ruv-rv
c	write(ounit,68) rv,rle,rue,rlv,ruv
c64	continue

c	write (ounit,643)
c643	format(t5,'XJBLC')

c	k=n+7*nl+49
c	do 644 i=25,28
c	l=jtcon+i
c	rv = 10**r(k)
c	rlv = 10**(r(k)-r(l))
c	ruv = 10**(r(k)+r(l))
c	rle = rv-rlv
c	rue = ruv-rv
c	write(ounit,68) rv,rle,rue,rlv,ruv
c644	continue

c	write (ounit,645)
c645	format(t5,'XJGLC')

c	k=n+7*nl+53
c	do 646 i=33,36
c	l=jtcon+i
c	rv = 10**r(k)
c	rlv = 10**(r(k)-r(l))
c	ruv = 10**(r(k)+r(l))
c	rle = rv-rlv
c	rue = ruv-rv
c	write(ounit,68) rv,rle,rue,rlv,ruv
c646	continue

	write (ounit,3)
	write (ounit,667) 
667	format(t20,'ERRORS IN SCALE:  JACKNIFE')
	write (ounit,668)
668	format(t20,'==========================')
	write(ounit,3)
	write (ounit,511)
511	format(t20,'SCALE',T30,'PLU/MIN',T40,'LOWER',T50,'UPPER')
	write (ounit,3)

	write(ounit,60)
60	format(t5,'XJKC')	

	k=n+7*nl+30
	do 62 i=13,16
	l=jtcon+i
	write(ounit,53) r(k),r(l),r(k)-r(l),r(k)+r(l)
62	continue

	write(ounit,601)
601	format(t5,'XJBC')	

	k=n+7*nl+41
	do 621 i=21,24
	l=jtcon+i
	write(ounit,53) r(k),r(l),r(k)-r(l),r(k)+r(l)
621	continue

	write(ounit,602)
602	format(t5,'XJGC')	

	k=n+7*nl+58
	do 622 i=29,32
	l=jtcon+i
	write(ounit,53) r(k),r(l),r(k)-r(l),r(k)+r(l)
622	continue

c--- chi-square and t-values and optional corrections

	write (ounit,1)
	write (ounit,5) filename
	write (ounit,2)
	write (ounit,70) n 
70	format(t5,'CHI-SQUARE AND T-VALUES FOR N = ',i3)
	write (ounit,701)
701	format(t5,'===================================')
	write (ounit,2)

	write (ounit,71)
71	format(t10,'CHI-LOW',T20,'CHI-HIGH')
	write (ounit,3)

	do 72 i=1,4
	l=jtchi+9+(i-1)*2
	write (ounit,73) r(l),r(l+1)
72	continue
73	format(t5,f10.2,t15,f10.2)

	write (ounit,2)
	write (ounit,2)
	write (ounit,74)
74	format(t5,'T-VALUES')
	write (ounit,741)
741	format(t5,'========')
	write (ounit,3)

	write (ounit,75) r(jtchi+17),r(jtchi+18),r(jtchi+19),r(jtchi+20)
75	format (t5,'t(n-1)',t17,f10.2,t27,f10.2,t37,f10.2,t47,
     +         f10.2)
	write (ounit,3)
	write (ounit,76) r(jtchi+21),r(jtchi+22),r(jtchi+23),r(jtchi+24)
76	format (t5,'t(n-1)/1.075',t17,f10.2,t27,f10.2,t37,f10.2,
     +         t47,f10.2)
	write (ounit,3)
	write (ounit,77) r(jtchi+25),r(jtchi+26),r(jtchi+27),r(jtchi+28)
77	format (t5,'t(0.7*(n-1))',t17,f10.2,t27,f10.2,t37,f10.2,
     +         t47,f10.2)

	write (ounit,2)
	write (ounit,2)
	write (ounit,774)
774	format(t5,'OPTIONAL COSMOLOGICAL CORRECTIONS')
	write (ounit,7741)
7741	format(t5,'=================================')
	write (ounit,3)

	write (ounit,78) r(jtcorr+1),r(jtcorr+2),r(jtcorr+3)
78	format(' LOCATION CONFIDENCE CORRECTION - add in quadrature: ',
     +       f10.2,
     +       // ' SCALE CORRECTION - subtract in quadrature: ',f10.2,
     +       // ' SCALE CONFIDENCE CORRECTION - add in quadrature: ',
     +       f10.2)
	 
	write(ounit,1)
	write(ounit,5) filename
	write(ounit,3)

	write(ounit,816)
816	format(t20,'ERRORS IN SCALE:  BOOTSTRAP')
	write(ounit,817)
817	format(t20,'===========================')
	write(ounit,3)

	write(ounit,818)
818	format(t20,'SCALE',T30,'LENGTH',T40,'LOWER',T50,'UPPER')

	write(ounit,819)
819	format(t5,'SDEV-1')

	k=n+7*nl+30
	do 820 i=85,92,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
820	continue

	write(ounit,821)
821	format(t5,'SDEV-2')

	k=n+7*nl+30
	do 822 i=93,100,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
822	continue

	write(ounit,823)
823	format(t5,'SDEV-3')

	k=n+7*nl+30
	do 824 i=101,108,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
824	continue

	write(ounit,831)
831	format(t5,'SDEV-4')

	k=n+7*nl+30
	do 832 i=149,156,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
832	continue

	write(ounit,825)
825	format(t5,'S-BWT-1')

	k=n+7*nl+41
	do 826 i=109,116,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
826	continue

	write(ounit,827)
827	format(t5,'S-BWT-2')

	do 828 i=117,124,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
828	continue

	write(ounit,829)
829	format(t5,'S-BWT-3')

	do 830 i=125,132,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
830	continue

	write(ounit,833)
833	format(t5,'S-BWT-4')

	do 834 i=157,164,2
	l=jtcon+i
	write(ounit,53) r(k),r(l+1)-r(l),r(l),r(l+1)
834	continue

c--- tails

	write (ounit,1)
	write (ounit,5) filename
	write (ounit,2)
	write (ounit,30)
30	format(t25,'TAILS OF THE DATA SET')
	write (ounit,301)
301	format(t25,'=====================')
	write (ounit,2)

	write (ounit,311) r(jtasy+1),r(jtasy+2)
311	format(t5,'ASYMMETRY INDEX FOR DATA SET:  ',f6.3,//
     +         t5,'TRANSFORMATION POWER:  ',f6.3,//)
	k = n+7*nl
	write (ounit,31) r(k+56),r(k+57)
31	format (t5,'RAW TAIL INDEX FOR DATA SET:  ',f5.3,//
     +          t5,'SCALED TAIL INDEX FOR DATA SET:  ',f5.3,//)
	write (ounit,32)
32	format(t5,'COMPARISON INDICES FOR OTHER DISTRIBUTIONS',//
     +         ,t5,'UNIFORM      =     0.84'/
     +         ,t5,'TRIANGULAR   =     0.91'/
     +         ,t5,'GAUSSIAN     =     1.00'/
     +         ,T5,'CN(.05,3)    =     1.02'/
     +         ,T5,'CN(.05,10)   =     1.03'/
     +         ,T5,'LOGISTIC     =     1.05'/
     +         ,T5,'CN(.20,3)    =     1.08'/
     +         ,T5,'DBL EXP      =     1.23'/
     +         ,T5,'CN(.20,10)   =     1.26'/
     +         ,T5,'SLASH        =     1.43'/
     +         ,T5,'CAUCHY       =     1.62'////)

c--- gaps

	write (ounit,33)
33	format(t25,'WEIGHTED GAPS IN THE DATA')
	write (ounit,331)
331	format(t25,'=========================')
	write (ounit,2)

	ngap = ifix(r(jtgap+1))
	write (ounit,34) (r(k),k=jtgap+2,jtgap+1+ngap)
34	format(5(1x,f7.3))

	nbig = ifix(r(jtgap+2+ngap))
	write (ounit,2)
	write (ounit,35) nbig
35	format(t5,'THERE ARE  ',i3,' BIG GAPS IN THE DATA SET')
	write (ounit,2)

	write (ounit,36)
36	format(t5,'IBIG',T15,'ZBIG')
	write (ounit,3)

	do 37 i=1,nbig	
	k=jtgap+3+ngap+(i-1)*2
	write (ounit,38) ifix(r(k)),r(k+1)
37	continue
38	format (t5,i3,t15,f7.3)

c--- normality test results

	write (ounit,1)
	write (ounit,5) filename
	write (ounit,2)

	write (ounit,40) 
40	format (t25,'RESULTS OF THE NORMALITY TESTS')
	write (ounit,401)
401	format (t25,'==============================')
	write (ounit,2)

	write (ounit,41) r(jtnor+1)
	write (ounit,42) r(jtnor+2)
	write (ounit,43) r(jtnor+3),r(jtnor+4)
	write (ounit,44) r(jtnor+5),r(jtnor+23)
	write (ounit,45) r(jtnor+6),r(jtnor+24)
	write (ounit,1053) r(jtnor+25),r(jtnor+26)
	write (ounit,46) r(jtnor+7),r(jtnor+8),r(jtnor+9)
	write (ounit,47) r(jtnor+10),r(jtnor+11),r(jtnor+12)
	write (ounit,48) r(jtnor+13),r(jtnor+14)
	write (ounit,49) r(jtnor+15),r(jtnor+16)
	write (ounit,1050) r(jtnor+17),r(jtnor+18)
	write (ounit,1051) r(jtnor+19),r(jtnor+20)
	write (ounit,1052) r(jtnor+21),r(jtnor+22)

41	format(T3,'A-TEST   STATISTIC:      ',F8.3)
42	FORMAT(T3,'U-TEST   STATISTIC:      ',F8.3)
43	FORMAT(T3,'W-TEST   STATISTIC:   ',F8.3,T37,'PROBABILITY: ',F8.3)
44	FORMAT(T3,'B1-TEST  STATISTIC:   ',F8.3,T37,'PROBABILITY: ',F8.3)
45	FORMAT(T3,'B2-TEST  STATISTIC:   ',F8.3,T37,'PROBABILITY: ',F8.3)
1053    FORMAT(T3,'B1 & B2 OMNI TEST :      ',F8.3,T37,'PROBABILITY: ',
     +         F8.3)
46	FORMAT(T3,'I-TEST   STATISTIC:   ',F8.3,T37,'90% POINT:   ',F8.3,
     +         T59,'95% POINT:  ',F8.3)
47	FORMAT(T3,'DIP      STATISTIC:      ',F8.3,T37,'LOWER MODAL: ',
     +         F8.2,T59,'UPPER MODAL:',F8.2)
48	FORMAT(T3,'KS       STATISTIC:   ',F8.3,T37,'PROBABILITY: ',F8.3)
49	FORMAT(T3,'KUIPER V-STATISTIC:   ',F8.3,T37,'PROBABILITY: ',F8.3)
1050    FORMAT(T3,'CRAMER VON-MISES W2-STAT:',F8.3,T37,'PROBABILITY: ',
     +         F8.3)
1051    FORMAT(T3,'WATSON U2-STATISTIC:     ',F8.3,T37,'PROBABILITY: ',
     +         F8.3)
1052    FORMAT(T3,'ANDERSON-DARLING A2-STAT:',F8.3,T37,'PROBABILITY: ',
     +         F8.3)

	return
	end


c------------------------------------------------------------------------------
        SUBROUTINE SIGMA (DATA,N,SXAVE,SXSDEV,NTOT)
c------------------------------------------------------------------------------

c--- This routine performs a "3-sigma" clipping of the data using the
c    cleaning technique described in Yahil and Vidal (1977).  The procedure
c    calls for deleting the data point most deviant from the sample mean,
c    redetermining the sample mean and standard deviation.  If the removed
c    point deviates by more than 3 standard deviations from the new mean
c    it is rejected as a contamination.  Each point in the list is tested
c    in order of their apparent deviation, until no further rejection occurs.
c    The list is now considered "clean," and the final estimates of mean and
c    standard deviation are accepted.
c
c--- Note that this routine will exit to the main program if the input number
c    of data point  N <= 5, or if clipped to NTOT=5
c****************************************************************************

        implicit real*4 (a-h,o-z)
        dimension data(n),tdata(50000),sdata(50000)

	if(n.le.5) then
	  sxave=-999
	  sxsdev=-999
	  ntot=n
	return
	endif

c--- copy the original data set

        do 5 i=1,n
        sdata(i) = data(i)
5       continue

c--- find intial estimate of sample mean and deviation

        ave  = PARAM(SDATA,N,1)              
        sdev = PARAM(SDATA,N,2)              
        nout=0

c--- create a list of data ordered by their distance from the sample mean

20      continue

	if(n-nout.eq.5) then		     
	  sxave=-999
	  sxsdev=-999
	  ntot=5
	return
	endif

        do 10 i=1,n-nout
          tdata(i) = abs(sdata(i)-ave)
10      continue

        CALL SORT2 (N-NOUT,TDATA,SDATA)      

c--- drop first deviant point

          temp = sdata(n-nout)
          sdata(n-nout) = -32368.
          ave  = PARAM(SDATA,N-NOUT,1)       
          sdev = PARAM(SDATA,N-NOUT,2)       

c--- test if most deviant point is greater than 3-sigma from sample mean
c    calculated without including it in the estimate,
c    if so, then delete it from consideration

        if(abs(temp-ave) .gt. 3*sdev) then   
          nout = nout + 1
        else
          sdata(n-nout)=temp
          ave =  PARAM(SDATA,N-NOUT,1)  
          sdev = PARAM(SDATA,N-NOUT,2)  
          sxave=ave
          sxsdev=sdev
          ntot=n-nout
          return
        endif
        goto 20

        end


c------------------------------------------------------------------------------
        SUBROUTINE SORT(N,RA)
c------------------------------------------------------------------------------
 
c--- Routine to do a heapsort of a data array RA
c    Stolen (unabashedly) from NUMERICAL RECIPES
c
c**************************************************************************** 
 
      implicit real*4 (a-h,o-z) 
      dimension ra(n)
 
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end


c------------------------------------------------------------------------------
      SUBROUTINE SORT2(N,RA,RB)
c------------------------------------------------------------------------------
 
c--- Routine to do a heapsort of a data array RA and RB
c    Stolen (unabashedly) from NUMERICAL RECIPES
c
c**************************************************************************** 
 
      implicit real*4 (a-h,o-z)
      dimension ra(n),rb(n)
 
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
          rrb=rb(l)
        else
          rra=ra(ir)
          rrb=rb(ir)
          ra(ir)=ra(1)
          rb(ir)=rb(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            rb(1)=rrb
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            rb(i)=rb(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
        rb(i)=rrb
      go to 10
      end


c------------------------------------------------------------------------------
	SUBROUTINE STATS (N, Z, D, V, W2, U2, A2, IFAULT)
c------------------------------------------------------------------------------

C	   ALGORITHM AS 248.2 APPL.STATIST. (1989), VOL. 38, NO. 3)
C
C	   Computes the goodness-of-fit statistics D, V, W**2, U**2
C	   and A**2 from the transformed Z values
	
	
	integer i, ifault, n, ni
	real*4 z(n), d, v, w2, u2, a2, ri, rn, dplus, dminus, d1, d2, 
     +  sumz, twon, zm1, a2sum, zero, small, half, one, two, twelve
C
C	   Initialize constants
C
	data zero / 0.0 / , small / 1.0e-37 / , half / 0.5 /,
     +      one / 1.0 /, two / 2.0 /, twelve / 12.0 /

	rn = n
	twon = two*rn
C
C	   Calculating the Kologorov statistics DPLUS, DMINUS, D
C	   and the Kuiper statistic V.
C
	dplus = zero
	dminus = zero
	ifault = 4
	do 10 i = 1, n
	   if (z(i) .le. zero .or. z(i) .ge. one) return
	   ri = i
	   d1 = ri / rn - z(i)
	   if (d1 .gt. dplus) dplus = d1
	   d2 = z(i) - (ri - one) / rn
	   if (d2 .gt. dminus) dminus = d2
10	continue
	ifault = 0
	d = dplus
	if (dminus .gt. dplus) d = dminus
	v = dplus + dminus
C
C	   Calculating the Cramer-Von Mises statistic W2 and the
C	   Watson statistic U2.
C
	w2 = one / (twelve * rn)
	sumz = zero
	do 20 i = 1,n
	w2 = w2 + (z(i) - (two * i - one) / twon) ** 2
	   sumz = sumz + z(i)
20	continue
	u2 = w2 - rn * (sumz / rn - half) ** 2
c
C          Calculating the Anderson-Darling statistic A2
C
	a2sum = zero
	ni = n
	do 30 i = 1, n
	   if (z(i) .lt. small) z(i) = small
	   zm1 = one - z(ni)
	   if (zm1 .lt. small) zm1 = small
	   a2sum = a2sum + (two * i - one) * (log(z(i)) + log(zm1))
	   ni = ni - 1
30	continue
	a2 = (-a2sum) / rn - rn
	return
	end


c------------------------------------------------------------------------------
	SUBROUTINE SYMMETRY (XDATA,N,TEE)
c------------------------------------------------------------------------------

c-- for a reference see Finch,S."Robust Univariate Test of Symmetry",J. of
c   American Stat. Assoc., June,1977,Vol.72 No.358,p.387.
c   This compares the gaps in the ordered data set from the right side to 
c   the left side of the median as a test of the symmetry.  
c
c******************************************************************************

	real num,den,inv,yi,yn,vi,prod,tee,wii
	real xdata(n),wi(50000)

	num=0
	den=0

c-- calculate the weights for each data point

	do 10 i=1,n/2
	wii=0.
	do 20 j=i,n-i
	inv=1./float(j)
	wii=wii+inv
20	continue
	wi(i)=wii
10	continue

c-- get the gaps in the data and calculte the statistic

	do 30 i=1,n/2
	yi = xdata(i+1) - xdata(i)
	yn = xdata(n-i+1) - xdata(n-i)
	if (yi+yn.ne.0) vi = (yn-yi)/(yi+yn)
	prod = wi(i)*vi
	num = num + prod

c-- have to normalize by the sum of the sqaure of the weights 

	den = den + (wi(i)*wi(i))
30	continue

	tee = num/(den**0.5)

	return
	end


c------------------------------------------------------------------------------
        SUBROUTINE TAIL (XDATA,N,TINDEX1,TINDEX2)
c------------------------------------------------------------------------------

c--- TAIL provides a robust estimate of the weight of the tails
c    of a symmetric distribution  (UREDA 322).  Essentially,
c    this means we calculate:
c
c                               Q(.90) - Q(.10)
c                   TINDEX1  =   ---------------
c                               Q(.75) - Q(.25)
c
c    where Q is simply the rank statistic closest to the ideal
c    depth DEPTHS.
c
c    TINDEX1 can then be compared to similar values as calculated
c    for a variety of interesting distributions.
c
c    TINDEX2 is the scaled version of TINDEX1 relative to a gaussian
c
c--- NOTE -- TAIL assumes that XDATA is sorted low to high on input
c
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        dimension xdata(n)      
        data dt/0.33333/
 
c--- obtain estimate of Q parameter
 
        d90 = (n+dt)*.90+ dt
        d10 = (n+dt)*.10+ dt
        d75 = (n+dt)*.75+ dt
        d25 = (n+dt)*.25+ dt
 
        ig90 = int(d90)
        ig10 = int(d10)
        ig75 = int(d75)
        ig25 = int(d25)
 
        r90 = d90 - ig90
        r10 = d10 - ig10
        r75 = d75 - ig75
        r25 = d25 - ig25

	if (n .le. 7) then
		tindex1 = 0.0
		tindex2 = 0.0
		goto 10
	endif
 
        q90 = (1.-r90)*xdata(ig90)+r90*xdata(ig90+1)
        q10 = (1.-r10)*xdata(ig10)+r10*xdata(ig10+1)
        q75 = (1.-r75)*xdata(ig75)+r75*xdata(ig75+1)
        q25 = (1.-r25)*xdata(ig25)+r25*xdata(ig25+1)
 
c--- obtain estimate of TINDEX
 
        tindex1 = (q90-q10)/(q75-q25)
        tindex2 = tindex1/1.90         

10      return
        end


c------------------------------------------------------------------------------
	SUBROUTINE TRANS(N,NLV,XLV,POW)
c------------------------------------------------------------------------------

c-- Trans gives the power for the data set to make a more symmetric
c   distribution (for a discussion see UREDA pp.105-111).
c
c--- CODE UNDER REVIEW  USE WITH CAUTION.
c**********************************************************************

	real*4 pp(8),pow
	real*4 xlv(15,2)

c-- from the midsummaries a value pp will be found for each and
c   the power trans will be the median of those values

	xm = xlv(1,1)
	do 10 i=2,nlv
	x=((xlv(i,1)-xm)**2 + (xlv(i,2)-xm)**2)/(4.*xm)
	y=((xlv(i,1)+xlv(i,2))/2.)-xm
	pp(i-1)=1.-(y/x)
10	continue

	nlvv=nlv-1

c-- need to get the median of the pp(i)

	CALL MDIAN1(PP,NLVV,XMED)
	pow=xmed
	return
	end


c------------------------------------------------------------------------------
        SUBROUTINE xTRIM (XDATA,N,ALPHA,TRIMMED)
c------------------------------------------------------------------------------
 
c--- The Trimmed mean finds the mean of the set of N ordered statistics
c    after trimming the percentage ALPHA from both ends of the set.
c    The value TRIMMED is returned as the mean of the trimmed set
c    of ordered statistics.  The weights given to the ordered stats
c    and the formula for computation of the trimmed mean can be
c    found on page 311 of UREDA.
c
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        dimension xdata(n)
        data d1,d2,n1,n2,zero/1.0,2.0,1,2,0.0/
 
        ig = int(alpha * n)
        r = (alpha*n) - float(ig)
        sum1 = zero
        do 11 i = ig+n2,n-ig-n1
        sum1 = sum1 + xdata(i)
11      continue
        sum3 = (d1-r)*(xdata(ig+n1) + xdata(n-ig))
        trimmed = (d1/(float(n)*(d1 -(d2*alpha))))*(sum3+sum1)
 
        return  
        end

 
c----------------------------------------------------------------------------- 
        SUBROUTINE TRIMEAN (XDATA,N,XMED,XFL,XFU,TRIMN)
c----------------------------------------------------------------------------- 
 
c--- The TRIMEAN is calculated for the set of N ORDERED statistics
c    passed in the array XDATA. The value of the TRIMEAN is returned
c    as TRIMN. The TRIMEAN is given by the formula
c
c                     TRIMEAN = .25(FL + 2M + FU)
c
c    where FL,M,FU are the lower fourth, median, and upper fourth
c    respectively and are passed in the array XLETTER. For more
c    information see pages 313,314 in UREDA.
c 
c***************************************************************************** 
 
 
        implicit real*4 (a-h,o-z)
        dimension xdata(n)
        data do4,d2/0.25,2.0/

        trimn = do4*(xfu + (d2*xmed) + xfl)
        return
        end


c----------------------------------------------------------------------------- 
	SUBROUTINE UPDATE (N,XDATA,M,XOUT)
c----------------------------------------------------------------------------- 

	real*4 xdata(n),xout(n)

c--- only copy data from XDATA to XOUT if it is not = -32368

	j=0
	miss=0

	do 10 i=1,n
  	  if (xdata(i).ne.-32368.) then
	    j=j+1
	    xout(j)=xdata(i)
	  else
	    miss = miss + 1
	  endif
10	continue

	m = j

	return
        end


c----------------------------------------------------------------------------- 
        SUBROUTINE XAD (XDATA,N,XMED,XADM)
c----------------------------------------------------------------------------- 
 
c--- The XAD subroutine calculates the Mean Absolute Deviation from
c    the sample median. The median, M , is subtracted from each 
c    ORDERED statistic and then the absolute value is taken.
c    The AD is then defined to be the mean of these differences.
c    new set of statistics and is returned as XADM. The AD can
c    be defined:
c
c                  ADM = mean{ abs(x(i) - M) }
c 
c    where the x(i) are the values passed in the array XDATA, and
c    the median, M, is passed in the array XLETTER. 
c    For more information see page 408 in UREDA.
c
c 
c****************************************************************************
 
        implicit real*4 (a-h,o-z)
        dimension xdata(n)
        data zero/0.0/
 
        xsum= zero
        do 11 i = 1,n
        xsum = abs(xdata(i) - xmed) + xsum
11      continue
        
        xadm = xsum/n

        return
        end


c----------------------------------------------------------------------------- 
        SUBROUTINE XBIWT (XDATA,N,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1)
c----------------------------------------------------------------------------- 

c--- The subroutine XBIWT provides an estimator of the location and
c    scale of the data set XDATA.  The scale uses the Biweight function
c    in the general formula of "A-estimators." This formula is given
c    on page of 416 in UREDA (formula 4). The BIWEIGHT scale estimate
c    is returned as the value XSBIWT. The BIWEIGHT function is given
c    by:
c
c                                  u((1-u*u)**2)     abs(u) <= 1
c                         f(u) =
c                                  0                 abs(u) >  1
c
c    where u is defined by
c
c                         u = (XDATA(I) - M) / c*MAD  .
c
c    M, MAD, and c are the median, the median absolute deviation from
c    the median, and the tuning constant respectively. The tuning
c    constant is a parameter which is chosen depending on the sample
c    size and the specific function being used for the scale estimate.
c    (See page 417 in UREDA).  Here we take c = 9.0.
c
c--- The biweght location is found using the formula:
c
c                         T = M + (sums)
c
c                         where M is the sample median and sums are
c                         as given on page 421 in UREDA
c
c                         the tuning constant c is set to 6.0 for calculation
c                         of the location as reccommended by Tukey ()
c
c--- NOTE that the biweight is meant to be an iterated estimator, but one
c    commonly only takes the first step of the iteration.  Here we report
c    both the one-step estimators (XLBIWT1, XSBIWT1) and the preferred
c    fully iterated versions (XLBIWT, XSBIWT).
c
c -- 03/09/92 - a check to see if the iteration converged was put in,
c               to make it quicker
c
c****************************************************************************

        implicit real*4 (a-h,o-z)
        dimension xdata(n),u1(50000),u2(50000),
     +  xlb(11),xsb(11)
        data zero,d6,d1,d5,d9/0.0,6.0,1.0,5.0,9.0/
	
c---	sort the data and find the median
	
	CALL MDIAN1(XDATA,N,XM)

c---	call xmad to find the median absolute deviation

	CALL XMAD(XDATA,N,XM,XMADM)

c---    must choose value of the tuning constant "c"
c       here c = 6.0 for the location estimator and
c       9.0 for the scale estimator

        c1 = d6
        c2 = d9

	if (xmadm.le..0001) then
	xlbiwt=xm
	xlbiwt1=xm
	xsbiwt=xmadm
	xsbiwt1=xmadm
	goto 20
	endif

        do 11 i = 1,n
	u1(i) = (xdata(i) - xm)/(c1*xmadm)
        u2(i) = (xdata(i) - xm)/(c2*xmadm)
11      continue

        s1 = zero
        s2 = zero
        s3 = zero
        s4 = zero

        do 12 i = 1,n
        if (abs(u2(i)) .lt. d1) then
            s1 = s1+(((xdata(i)-xm)**2)*(d1-(u2(i)*u2(i)))**4)
            s2 = s2+((d1-u2(i)*u2(i))*(d1-(d5*u2(i)*u2(i))))
        endif
        if (abs(u1(i)) .lt. d1) then
            s3 = s3+(xdata(i)-xm)*(d1-u1(i)*u1(i))**2
            s4 = s4+(d1-u1(i)*u1(i))**2
        endif
12      continue

c--- here are the one-step estimators

        xlbiwt1 = xm+s3/s4
        xsbiwt1 = float(n)/(float(n-1))**0.5*s1**0.5/abs(s2)

c--- now obtain the fully-iterated versions

c--- solve for new estimates of u1 and u2

        xlb(1) = xlbiwt1
        xsb(1) = xsbiwt1

        do 15 j = 2,11   

        xmm = xlb(j-1)
        do 13 i = 1,n
        u1(i) = (xdata(i) - xmm)/(c1*xmadm)
        u2(i) = (xdata(i) - xmm)/(c2*xmadm)
13      continue

        s1 = zero
        s2 = zero
        s3 = zero
        s4 = zero

        do 14 i = 1,n
        if (abs(u2(i)) .lt. d1) then
            s1 = s1+(((xdata(i)-xmm)**2)*(d1-(u2(i)*u2(i)))**4)
            s2 = s2+((d1-u2(i)*u2(i))*(d1-(d5*u2(i)*u2(i))))
        endif
        if (abs(u1(i)) .lt. d1) then
            s3 = s3+(xdata(i)-xmm)*(d1-u1(i)*u1(i))**2
            s4 = s4+(d1-u1(i)*u1(i))**2
        endif
14      continue

        xlb(j) = xlb(j-1)+s3/s4
        xsb(j) = float(n)/(float(n-1))**0.5*s1**0.5/abs(s2)

        if(abs((xlb(j)-xlb(j-1))/xlb(j)).lt.1.e-5) then
           xlbiwt=xlb(j)
           xsbiwt=xsb(j)
           return
        endif
	
15      continue

        xlbiwt = xlb(11)
        xsbiwt = xsb(11)

20	continue
        return
        end


c----------------------------------------------------------------------------- 
        SUBROUTINE XMAD (XDATA,N,XMED,XMADM)
c----------------------------------------------------------------------------- 
 
c--- The XMAD subroutine calculates the Median Absolute Deviation from
c    the sample median. The median, M , is subtracted from each 
c    ORDERED statistic and then the absolute value is taken. This new
c    set of of statistics is then resorted so that they are ORDERED
c    statistics. The MAD is then defined to be the median of this
c    new set of statistics and is returned as XMADM. The MAD can
c    be defined:
c
c                   XMADM = median{ abs(x(i) - M) }
c 
c    where the x(i) are the values passed in the array XDATA, and
c    the median, M, is passed in the array XLETTER. The set of stats
c    in the brackets is assumed to be resorted. For more information
c    see page 408 in UREDA.
c
c 
c****************************************************************************
 
        implicit real*4 (a-h,o-z)
        dimension xdata2(50000),xdata(n)
        data dhalf,n1,n2/0.5,1,2/
 
        do 11 i = 1,n
        xdata2(i) = abs(xdata(i) - xmed)
11      continue
        
        CALL SORT(N,XDATA2)
 
        if (float(n)/float(n2) - int(n/n2) .eq. 0) then
                i1 = n/n2
                i2 = n/n2 + n1
                xmadm = dhalf*(xdata2(i1) + xdata2(i2))
        else    
                i1 = int(n/n2) + n1
                xmadm = xdata2(i1)
        endif

        return
        end


c----------------------------------------------------------------------------- 
        SUBROUTINE XMIDMEAN(XDATA,N,XMID)
c----------------------------------------------------------------------------- 
 
c--- This subroutine calculates the MIDMEAN for a set of N ORDERED
c    statistics stored in XDATA. The MIDMEAN is defined to be the
c    mean of the central 50% of the data set. This corresponds to 
c    a 25% Trimmed Mean. The value of the MIDMEAN is returned as 
c    XMID and is defined:
c
c                   XMID = TRIM(.25) 
c
c    where TRIM(.25) is the 25% trimmed mean as defined above. For
c    more information on the MIDMEAN or TRIMMED MEAN see pages 312,
c    313 and pages 307, 308 in UREDA.
 
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        dimension xdata(n)
        data n1,n2,d1,d2,df,zero/1,2,1.0,2.0,0.25,0.0/

        ig = int(df*n)
        r = (df*n) - float(ig)
        sum1 = zero
        do 11 i=ig + n2,n-ig-n1
        sum1 = sum1 + xdata(i)
11      continue
        sum3 = (d1-r)*(xdata(ig+n1) + xdata(n-ig))
        xmid = (d1/(float(n)*(d1-(d2*df))))*(sum3+sum1)

        return
        end


c----------------------------------------------------------------------------- 
        SUBROUTINE WCOEF(A,N,N2,EPS,IFAULT)
c----------------------------------------------------------------------------- 

c--- Algorithm AS 181.1 Appl. Stat. (1982) Vol.31, No.2
c
c--- Obtain array A of weights for calculating W
c
        real*4 a(n2),c4(2),c5(2),c6(3)
        data c4(1),c4(2)/.6869,.1678/,c5(1),c5(2)/.6647,.2412/,
     +       c6(1),c6(2),c6(3)/.6431,.2806,.0875/
        data rsqrt2/.70710678/,zero/0.0/,half/.5/,one/1.0/,
     +       two/2.0/,six/6.0/,seven/7.0/,eight/8.0/,thirt/13./
        ifault = 1
        if (n .le. 2) return
        ifault = 3
        if (n/2 .ne.n2) return
        ifault = 2
        if (n .gt. 2000) return
        ifault = 0
        if (n .le. 6) goto 30
c
c--- n .gt. 6 calculate rankits using approximate routine NSCOR2
c    (as177)
c
        CALL NSCOR2(A,N,N2,IFAULT)
        sastar = zero
        do 10 j = 2,n2
 10     sastar = sastar + a(j)*a(j)
        sastar = sastar*eight
        nn = n
        if (n .le. 20) nn = nn-1
        an = nn
        a1sq = exp(alog(six*an + seven) - alog(six*an + thirt)
     +         + half*(one + (an - two)*alog(an + one) - (an - one)
     +         *alog(an + two)))
        a1star = sastar/(one/a1sq - two)
        sastar = sqrt(sastar+two*a1star)
        a(1) = sqrt(a1star)/sastar
        do 20 j = 2,n2
 20     a(j) = two*a(j)/sastar
        goto 70
c
c--- n .le. 6 use exact values for weights
c
 30     a(1) = rsqrt2
        if (n .eq. 3) goto 70
        n3 = n - 3
        goto (40,50,60), n3
 40     do 45 j = 1,2
 45     a(j) = c4(j)
        goto 70
 50     do 55 j = 1,2
 55     a(j) = c5(j)
        goto 70
 60     do 65 j = 1,3
 65     a(j) = c6(j)

c
c--- calculate the minimum possible value of w
c

 70     eps = a(1)*a(1)/(one - one/float(n))
        return
        end


c----------------------------------------------------------------------------- 
        SUBROUTINE WEXT(X,N,SSQ,A,N2,EPS,W,PW,IFAULT)
c----------------------------------------------------------------------------- 

c--- Algorithm AS 181 Appl. Stat. (1982) Vol.31, No. 2
c
c--- Calculates Shapiro and Wilk W statistic and its sig. level
c
c--- NOTE that WEXT assumes array X is sorted low to high on entry
c
c    W is evaluated for N in the range 3 < N < 2000
c
c******************************************************************************

        real*4 x(n),a(n2),lamda,wa(3),wb(4),wc(4),wd(6),we(6),
     +       wf(7),c1(5,3),c2(5,3),c(5),unl(3),unh(3)
        integer nc1(3),nc2(3)
        logical upper
        data wa(1),wa(2),wa(3)/.118898,.133414,.327907/,
     +  wb(1),wb(2),wb(3),wb(4)/-.37542,-.492145,-1.124332,-.199422/,
     +  wc(1),wc(2),wc(3),wc(4)/-3.15805,.729399,3.01855,1.558776/,
     +  wd(1),wd(2),wd(3),wd(4),wd(5),wd(6)/
     +  .480385,.318828,0.0,-.0241665,.00879701,.002989646/,
     +  we(1),we(2),we(3),we(4),we(5),we(6)/
     +  -1.91487,-1.37888,-.04183209,.1066339,-.03513666,
     +  -.01504614/,
     +  wf(1),wf(2),wf(3),wf(4),wf(5),wf(6),wf(7)/
     +  -3.73538,-1.015807,-.331885,.1773538,-.01638782,
     +  -.03215018,.003852646/
        data c1(1,1),c1(2,1),c1(3,1),c1(4,1),c1(5,1),
     +       c1(1,2),c1(2,2),c1(3,2),c1(4,2),c1(5,2),
     +       c1(1,3),c1(2,3),c1(3,3),c1(4,3),c1(5,3)/
     +       -1.26233,1.87969,.0649583,-.0475604,-.0139682,
     +       -2.28135,2.26186,0.0,0.0,-.00865763,
     +       -3.30623,2.76287,-.83484,1.20857,-.507590/
        data c2(1,1),c2(2,1),c2(3,1),c2(4,1),c2(5,1),
     +       c2(1,2),c2(2,2),c2(3,2),c2(4,2),c2(5,2),
     +       c2(1,3),c2(2,3),c2(3,3),c2(4,3),c2(5,3)/
     +       -.287696,1.78953,-.180114,0.0,0.0,
     +       -1.63638,5.60924,-3.63738,1.08439,0.0,
     +       -5.991908,21.04575,-24.58061,13.78661,-2.835295/
        data unl(1),unl(2),unl(3)/-3.8,-3.0,-1.0/,
     +       unh(1),unh(2),unh(3)/8.6,  5.8, 5.4/
        data nc1(1),nc1(2),nc1(3)/5,5,5/,
     +       nc2(1),nc2(2),nc2(3)/3,4,5/
        data pi6/1.90985932/, stqr/1.04710755/, upper/.true./,
     +       zero/0.0/,tqr/.75/,one/1.0/,onept4/1.4/,three/3.0/,
     +       five/5.0/
        ifault = 1
        pw = one
        w = one
        if (n .le. 2) return
        ifault = 3
        if (n/2 .ne. n2) return
        ifault = 2
        if (n .gt. 2000) return
c
c--- calculate w
c
        ifault = 0
        w = zero
        an = n
        i = n
        do 10 j = 1,n2
        w = w + a(j)*(x(i) - x(j))
        i = i - 1
 10     continue
        w = w*w/ssq
        if (w .lt. one) goto 20
        w = one
        return
c
c--- get significance level of w
c
 20     if (n .le. 6) goto 100
c
c--- n between 7 and 2000... transform w to y, get mean and sd,
c    standardize and get significance level
c
        if (n .gt. 20) goto 30
        al = alog(an) - three
        lamda = poly(wa,3,al)
        ybar = exp(POLY(WB,4,AL))
        sdy = exp(POLY(WC,4,AL))
        goto 40
 30     al = alog(an) - five
        lamda = POLY(WD,6,AL)
        ybar = exp(POLY(WE,6,AL))
        sdy = exp(POLY(WF,7,AL))
 40     y = (one - w)**lamda
        z = (y - ybar)/sdy
        
c--- if z is very large ( > 10), reset to 10
	
	if(z.gt.10) z = 10
	pw = ALNORM(Z,UPPER)
        return
c
c--- deal with n less than 7 (exact significance level for n=3)
c
 100    if (w .le. eps) goto 160
        ww = w
        if (n .eq. 3) goto 150
        un = alog((w-eps)/(one - w))
        n3 = n - 3
        if (un .lt. unl(n3)) goto 160
        if (un .gt. onept4)  goto 120
        nc = nc1(n3)
        do 110 i = 1,nc
 110    c(i) = c1(i,n3)
        eu3 = exp(POLY(C,NC,UN))
        goto 140
 120    if (un .gt. unh(n3)) return
        nc = nc2(n3)
        do 130 i = 1,nc
 130    c(i) = c2(i,n3)
        un = alog(un)
        eu3 = exp(exp(POLY(C,NC,UN)))
 140    ww = (eu3 + tqr)/(one + eu3)
 150    pw = pi6*(atan(sqrt(ww/(1.0 - ww))) - stqr)
        return
 160    pw = zero
        return
        end

 
c**************************AUXILLIARY FUNCTIONS*******************************
 

c----------------------------------------------------------------------------- 
        REAL FUNCTION ALNFAC(J)
c----------------------------------------------------------------------------- 
c
c--- Algorithm AS 177.2 Appl. Stat. (1982) Vol.31,No.2
c
c--- Natural logarithm of factorial for non-negative argument
c
        real*4 r(7),one,half,ao,three,four,fourtn,fortty,fivfty,w,z
        data r(1),r(2),r(3),r(4),r(5),r(6),r(7)/0.0e0,0.0e0,
     +       .69314718056e0,1.79175946923e0,3.17805383035e0,
     +       4.78749174178e0,6.57925121101e0/
        data one,half,ao,three,four,fourtn,fortty,fivfty/
     +       1.0e0,.5e0,.918938533205e0,3.0e0,4.0e0,14.0e0,420.0e0
     +       ,5040.0e0/
        if (j .ge. 0) goto 10
        alnfac = one
        return
 10     if (j .ge. 7) goto 20
        alnfac = r(j+1)
        return
 20     w = j + 1
        z = one / (w*w)
        alnfac = (w - half)*alog(w) - w + ao + (((four - three*z)
     +           *z - fourtn)*z + fortty)/(fivfty*w)
        return
        end


c----------------------------------------------------------------------------- 
        FUNCTION ALNORM(X,UPPER)
c----------------------------------------------------------------------------- 

c--- Algorithm AS 66  Appl. Stat. (1973) Vol. 22, No.3
c
c--- Evaluates the tail area of the standardised normal curve
c    from x to infinity if upper is .true. or from minus
c    infinity to x if upper is .false.
c
        real*4 ltone,utzero,zero,half,one,con,z,y,x
        logical upper,up
c
c--- ltone and utzero must be set to suit the particular
c    computer (see introductory text)
c
        data ltone,utzero/7.0,18.66/
        data zero,half,one,con/0.0,.5,1.0,1.28/
        up = upper
        z = x
        if (z .ge. zero) goto 10
        up = .not.up
        z = -z
 10     if ((z .le. ltone) .or. (up .and. (z .le. utzero))) goto 20
        alnorm = zero
        goto 40
 20     y = half*z*z
        if (z .gt. con) goto 30
c
c
        alnorm = half - z*(.398942280444 - .399903438504*y/
     +           (y + 5.75885480458 - 29.8213557808 /
     +           (y + 2.62433121679 + 48.6959930692 /
     +           (y + 5.92885724438))))
        goto 40
c
c
 30     alnorm = .398942280385*exp(-y)/
     +           (z + 3.805e-8 + 1.00000615302/
     +           (z + 3.98064794e-4 + 1.98615381364/
     +           (z - .151679116635 + 5.29330324926/
     +           (z + 4.8385912808 - 15.1508972451/
     +           (z + .742380924027 + 30.789933034/
     +           (z + 3.99019417011))))))
c
c
 40     if (.not.up) alnorm = one - alnorm
        return
        end


c----------------------------------------------------------------------------- 
        REAL FUNCTION CORREC(I,N)
c----------------------------------------------------------------------------- 

c--- Algorithm AS 177.4 Appl. Stat. (1982) Vol.31,No.2
c
c--- Calculates correction for tail area of normal distribution 
c    corresponding to ith largest rankit in sample size n.
c
        real*4 c1(7),c2(7),c3(7),an,mic,c14
        data c1(1),c1(2),c1(3),c1(4),c1(5),c1(6),c1(7)
     +       /9.5e0,28.7e0,1.9e0,0.0e0,-7.0e0,-6.2e0,-1.6e0/,
     +       c2(1),c2(2),c2(3),c2(4),c2(5),c2(6),c2(7)
     +       /-6.195e3,-9.569e3,-6.728e3,-17.614e3,-8.278e3,-3.570e3,
     +       1.075e3/,
     +       c3(1),c3(2),c3(3),c3(4),c3(5),c3(6),c3(7)
     +       /9.338e4,1.7516e5,4.1040e5,2.157e6,2.376e6,2.065e6,
     +       2.065e6/,
     +       mic/1.0e-6/,c14/1.9e-5/
        correc = c14
        if (i*n .eq. 4) return
        correc = 0.0
        if (i .lt. 1 .or. i .gt. 7) return
        if (i .ne. 4 .and. n .gt. 20) return
        if (i .eq. 4 .and. n .gt. 40) return
        an = n
        an = 1.0/(an*an)
        correc = (c1(i) + an*(c2(i) + an*c3(i))) * mic
        return
        end


c----------------------------------------------------------------------------- 
	FUNCTION CVALUE (HL,NDOF,PROB)
c----------------------------------------------------------------------------- 
	 
c--- CVALUE returns the value of chi-square corresponding to
c    the area of the Gaussian PROB (expressed in decimal notation)
c    and NDOF degrees of freedom via a simple look-up table
c
c--- NOTE that only the values 0.68, 0.90, 0.95, and 0.99 are allowed for
c    the value of PROB
c
c--- NOTE that this routine returns the hi or low values according to whether
c    the parameter HL = 0 (low) or 1 (high)
 
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        integer iunit2,iunit3,iunit4,ounit3,ounit4
        common iunit2,iunit3,iunit4,ounit3,ounit4
        integer ndof,hl
 
c--- check if PROB is one of the allowed inputs
 
        if(prob.ne.0.68.and.prob.ne.0.90.and.prob.ne.0.95
     +     .and.prob.ne.0.99) then
        print *, ' IMPROPER INPUT VALUE OF PROB, PLEASE CHECK '
        return
        endif
 
        if(hl.eq.0) then
        iunit=iunit3
        else
        iunit=iunit4
        endif
 
        do 10 i=1,200
        read (iunit,*,end=20) dof,c68,c90,c95,c99
        if(ifix(dof).eq.ndof) go to 20
10      continue
 
20      if(prob.eq.0.68) then
        cvalue = c68
        else if(prob.eq.0.90) then
        cvalue = c90
        else if(prob.eq.0.95) then
        cvalue = c95
        else
        cvalue = c99
        endif
 
        rewind(iunit)
        return
 
        end

c----------------------------------------------------------------------------- 
      FUNCTION GAMMLN(XX)
c----------------------------------------------------------------------------- 

c---  taken from Numerical Recipes

c*****************************************************************************

      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)

      return
      end

c----------------------------------------------------------------------------- 
      FUNCTION GAMMQ(A,X)
c----------------------------------------------------------------------------- 

c---  taken from Numerical Recipes, used to calculate the chi-square prob.

c*****************************************************************************

      if(x.lt.0..or.a.le.0.) print *,'gammq'
      if(x.lt.a+1.)then
        CALL GSER(GAMSER,A,X,GLN)
        gammq=1.-gamser
      else
        CALL GCF(GAMMQ,A,X,GLN)

      endif
      return
      end

c----------------------------------------------------------------------------- 
        FUNCTION PARAM (XDATA,N,IP)
c----------------------------------------------------------------------------- 

c--- This routine is a user-supplied function which calculates the
c    value of a statistic PARAM, chosen according to the value of IP

c---            IP = 1;         sample mean
c               IP = 2;         sample standard deviation
c               IP = 3;         sample geometric mean (sum of 1/DATA(I))
c               IP = 4;         sample harmonic mean (sum of 1/DATA(I)**2
c               IP = 5;         sample pairwise harmonic mean
c		IP = 6;		BIWEIGHT scale estimator
c		IP = 7;		GAPPED scale estimator
c
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z) 
        integer n,ip,miss
	integer jbig(50000)
        real xdata(n),xout(50000)
	real zstar(50000),zbig(50000)
	
c--- MEAN

        if(ip.eq.1) then 
          sum = 0.0
          miss = 0
          do 10 i=1,n
            if(xdata(i).ne.-32368.) then
              sum = sum + xdata(i)
            else
              miss = miss + 1            
            endif
10	  continue
 
          ave = sum/(n-miss)
          param = ave
 
c--- STANDARD DEVIATION
 
        else if (ip.eq.2) then
 
          sum = 0.0
          miss = 0
          do 20 i=1,n
            if(xdata(i).ne.-32368.) then
              sum = sum + xdata(i)
            else
              miss = miss + 1            
            endif
20	  continue
 
          ave = sum/(n-miss)
 
          ssum= 0.0
          do 21 i=1,n
            if(xdata(i).ne.-32368.) then
              ssum=ssum+(xdata(i)-ave)**2
            endif
21	   continue
 
          sdev = (ssum/(n-miss-1))**0.5
          param = sdev
 
c--- GEOMETRIC MEAN

        else if (ip.eq.3) then

          sum = 0.0
          miss = 0
          do 30 i=1,n
            if(xdata(i).ne.-32368.) then
              sum = sum + 1/xdata(i)
            else
              miss = miss + 1
            endif
30	  continue
 
          gave = 1./(sum/(n-miss))
          param = gave
 
 
c--- HARMONIC MEAN
 
         else if (ip.eq.4) then
 
          sum = 0.0
          miss = 0
          do 40 i=1,n
            if(xdata(i).ne.-32368.) then
              sum = sum + 1/xdata(i)**2
            else
              miss = miss + 1
            endif
40	   continue
 
          have = sqrt(1./(sum/(n-miss)))
          param = have
 
c--- PAIRWISE HARMONIC MEAN                     
 
        else if (ip.eq.5) then
 
          sum = 0.0
          miss = 0
          do 50 i=1,n
            if(xdata(i).ne.-32368.) then
              sum = sum + 1/xdata(i)
            else
              miss = miss + 1
            endif
50 	   continue
 
          pave = 1./(sum/(n-miss))
          param = pave

c--- BIWEIGHT SCALE ESTIMATOR
        
	else if (ip.eq.6) then

c---	Rewrite array dropping the value = -32368.

	CALL UPDATE (N,XDATA,M,XOUT)

c---	now obtain biweight estimator for the new array

	CALL XBIWT (XOUT,M,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1)
	param = xsbiwt

c--- GAPPER SCALE ESTIMATOR

	else if (ip.eq.7) then

c---	Rewrite array dropping valeu = -32368.

	CALL UPDATE (N,XDATA,M,XOUT)

c---	now obtain gapped estimator for the new array

	CALL GAPPER (XOUT,M,XSGAP,ZSTAR,NBIG,JBIG,ZBIG)
	param = xsgap

	endif
        return

	end


c----------------------------------------------------------------------------- 
        FUNCTION POLY(C,NORD,X)
c----------------------------------------------------------------------------- 

c--- Algorithm AS 181.2 Appl. Stat. (1982) Vol.31, No.2
c
c--- Calculates the algebraic polynomial of order nord-1 with
c    array of coefficients c.  Zero order coefficient is c(1).
c
 
        real*4 c(nord)
        poly = c(1)
        if (nord .eq. 1) return
        p = x*c(nord)
        if (nord .eq. 2) goto 20
        n2 = nord - 2
        j = n2 + 1
        do 10 i = 1,n2
        p = (p + c(j))*x
        j = j-1
 10     continue
 20     poly = poly + p
        return
        end


c----------------------------------------------------------------------------- 
        REAL FUNCTION PPND(P,IFAULT)
c----------------------------------------------------------------------------- 

c--- Algorithm AS 111  Appl. Stat. (1977), Vol. 26, No.1
c
c--- Produces normal deviate corresponding to lower tail area of P
c    real version for eps = 2**(-31)
c    The hash sums are the sums of the moduli of the coefficients
c    they have no inherent meanings but are included for use in
c    checking transcriptions.
c    Standard functions  abs,alog and sqrt are used.
c
        real*4 zero,split,half,one
        real*4 a0,a1,a2,a3,b1,b2,b3,b4,c0,c1,c2,c3,d1,d2
        real*4 p,q,r
        data zero/0.0e0/,half/.5e0/,one/1.0e0/
        data split /.42e0/
        data a0/   2.50662823884e0/
        data a1/ -18.61500062529e0/
        data a2/  41.39119773534e0/
        data a3/ -25.44106049637e0/
        data b1/  -8.47351093090e0/
        data b2/  23.08336743743e0/
        data b3/ -21.06224101826e0/
        data b4/   3.13082909833e0/
c
c--- hash sum ab 143.70383558076
c
        data c0/   -2.78718931138e0/
        data c1/   -2.29796479134e0/
        data c2/    4.85914127135e0/
        data c3/    2.32121276858e0/
        data d1/    3.54388924762e0/
        data d2/    1.63706781897e0/
c
c--- hash sum cd 17.43746520924
c
        ifault = 0
        q = p - half
        if (abs(q) .gt. split) goto 1
        r = q*q
        ppnd = q*(((a3*r + a2)*r + a1)*r + a0)/
     +         ((((b4*r + b3)*r + b2)*r + b1)*r + one)
        return
 1      r = p
        if (q .gt. zero) r = one - p
        if (r .le. zero) goto 2
        r = sqrt(-alog(r))
        ppnd = (((c3*r + c2)*r + c1)*r + c0)/
     +         ((d2*r + d1)*r + one)
        if (q .lt. zero) ppnd = -ppnd
        return
 2      ifault = 1
        ppnd = zero
        return
        end


c----------------------------------------------------------------------------- 
	REAL FUNCTION PVALUE (T, CUT, COEF)
c----------------------------------------------------------------------------- 

C	ALGORITHM AS 248.3 APPL.STATIST. (1989), VOL. 38, NO. 3)
C
C	Computes the approximate significance level for W**2, U**2,
C	or A**2 when testing for normality or exponentiality
C
	integer i, j
	real*4 t, cut(3), coef(4,3), one
	data one / 1.0 /

	i = 4
	do 10 j = 1, 3
	   if (t .lt. cut(j)) i = i - 1
10	continue
	pvalue = exp (coef(i, 1) + t * (coef( i, 2) + t * coef(i, 3)))
	if (i . lt. 3) pvalue = one - pvalue
	return
	end


c----------------------------------------------------------------------------- 
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

      FUNCTION RAN1old(IDUM)
c----------------------------------------------------------------------------- 

c--- returns a random deviate from a uniform distribution between 0.0 and 1.0
c
c    STOLEN (unabashedly) from Numerical Recipes
c
c******************************************************************************

      dimension r(97)
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
      parameter (m3=243000,ia3=4561,ic3=51349)
      data iff /0/
      if (idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,97
          ix1=mod(ia1*ix1+ic1,m1)
          ix2=mod(ia2*ix2+ic2,m2)
          r(j)=(float(ix1)+float(ix2)*rm2)*rm1
11      continue
        idum=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j.gt.97.or.j.lt.1) print *,'ran1'
      ran1=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end


c----------------------------------------------------------------------------- 
        FUNCTION TVALUE (NDOF,PROB)
c----------------------------------------------------------------------------- 
 
c--- TVALUE returns the value of Student's t corresponding to
c    the area of the Gaussian PROB (expressed in decimal notation)
c    and NDOF degrees of freedom via a simple look-up table
c
c--- NOTE that only the values 0.68, 0.90, 0.95, and 0.99 are allowed for
c    the value of PROB
 
c**************************************************************************** 
 
        implicit real*4 (a-h,o-z)
        integer iunit2,iunit3,iunit4,ounit3,ounit4
        common iunit2,iunit3,iunit4,ounit3,ounit4

c--- check if PROB is one of the allowed inputs
 
        if(prob.ne.0.68.and.prob.ne.0.90.and.prob.ne.0.95
     +     .and.prob.ne.0.99) then
        print *, ' IMPROPER INPUT VALUE OF PROB, PLEASE CHECK '
        return
        endif
 
        do 10 i=1,200
        read (iunit2,*,end=20) dof,t68,t90,t95,t99
        if(ifix(dof).eq.ndof) go to 20
10      continue
 
20      if(prob.eq.0.68) then
        tvalue = t68
        else if(prob.eq.0.90) then
        tvalue = t90
        else if(prob.eq.0.95) then
        tvalue = t95
        else
        tvalue = t99
        endif
 
        rewind(iunit2)
        return
 
        end
 
