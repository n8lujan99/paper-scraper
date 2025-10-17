C-----------------------------------------------------------------------------D
C     subroutine COMPAREWRITE writes out the light values of each bin used
C       in the chi-square minimization
C
C     USED BY MODEL
C
C-----------------------------------------------------------------------------D
      SUBROUTINE comparewrite(istype,io,ifilter,iter,rcond)
      INCLUDE 'moddefs.h'
      DIMENSION errA(Nbin),errR(Nbin)
      DOUBLE PRECISION S,dS(Norbitm+Nvel*Nvelbm)
      DOUBLE PRECISION ddS(Norbitm+Nvel*Nvelbm),rcond

111   format ('____________________________________________________',
     &  '___________________________')
735   format ('worst rel.:  ir,iv = ',i4,',',i3,'   error = ',f14.8)
738   format ('       mean relative error = ',f14.8)
888   format ('iter=',i3,'   m/l=',f11.5,'    S=',d12.6,'   S-type=',
     &  i1,'    rcond=',d12.6)

      write (io,111)
      write (io,*)''
      write (io,*)'iteration number ',iter
      write (io,*)''

      if (ifilter.eq.1) call filter()
      call summ()
      call comparer(errA,errR)
      call minmax(errA,errR,errAmin,errRmin,errAmax,errRmax,
     &     ibAmin,ibAmax,ibRmin,ibRmax)
      call qual(errA,errR,aveA,aveR,adevA,adevR,sdevA,sdevR)

      xlbad=errRmax

      write (io,735)itor(ibRmax),itov(ibRmax),errRmax
      write (io,*)''
      write (io,738)aveR
      write (io,*)''

c 801  format(<Nvelb>(1x,f7.2))
 801  format(40(1x,f13.2))
      write(io,801) (chi1(i)/float(Nvel),i=1,Nvelb)

      write (io,*)''

      call entropy(S,dS,ddS,istype)
      write (io,888)iter,1./xmu**2,S,istype,rcond
      write (io,*)''
      ent=sngl(S)

      RETURN
      END
