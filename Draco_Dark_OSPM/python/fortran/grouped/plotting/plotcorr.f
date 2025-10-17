
      parameter(nmax=10000)
      real s(nmax),dd(nmax),dr(nmax),rr(nmax),corr(nmax)
      real corr2(nmax)

      f=0.206719
      off=-1.92
      read *,f

      f2=0.21
      off2=2.82
      xfp=0.2
      
      open(unit=1,file='in',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4
         n=n+1
         s(n)=x1
         dd(n)=x2
         dr(n)=x3
         rr(n)=x4
      enddo
 666  continue
      close(1)

      do i=1,n
         corr(i)=(dd(i)-2.*f*dr(i)+f*f*rr(i))/(f*f*rr(i))
c         corr(i)=(dd(i)-2.*f*dr(i)+f*f*rr(i))/(f*f*rr(i))+off
c         corr(i)=corr(i)*s(i)*s(i)
         corr2(i)=dd(i)+xfp*xfp*f2*f2*rr(i)-2.*dr(i)*xfp*f2-
     $        (1.-xfp)*2.*f2*dr(i)+(1.-xfp)*2.*xfp*f2*f2*rr(i)+
     $        (1.-xfp)**2*f2*f2*rr(i)
         corr2(i)=corr2(i)/((1.-xfp)**2*f2*f2*rr(i))+2.598
         corr2(i)=corr2(i)*s(i)
         print *,s(i),corr(i),corr2(i)
      enddo

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      xmin=70.
      xmax=180.
      ymin=-5.
      ymax=10.
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("s(Mpc)","s\U2\DX","")
      call pgslw(6)
      call pgsci(2)
      call pgline(n,s,corr)
      call pgsci(4)
c      call pgline(n,s,corr2)
      call pgslw(2)
      call pgsci(1)

      call pgend

      end
