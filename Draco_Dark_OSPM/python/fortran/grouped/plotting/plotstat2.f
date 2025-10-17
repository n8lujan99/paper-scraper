c 4.9, 5.2, 5.8, 6.5, 11
      parameter(nmax=100)
      real xo1(nmax),xi1(nmax),xs1(nmax)
      real xo2(nmax),xi2(nmax),xs2(nmax)
      real xo3(nmax),xi3(nmax),xs3(nmax)
      real xo4(nmax),xi4(nmax),xs4(nmax)

      n1=0
      n2=0
      n3=0
      n4=0
      open(unit=1,file='stat.use',status='old')
      do i=1,1000
         read(1,*,end=666) i1,i2,x3,x4,x5,x6
         if(x6.eq.4.9) then
            n1=n1+1
            xo1(n1)=x5
            xi1(n1)=x3/x5
            xs1(n1)=x4/x5
         endif
         if(x6.eq.5.2) then
            n2=n2+1
            xo2(n2)=x5
            xi2(n2)=x3/x5
            xs2(n2)=x4/x5
         endif
         if(x6.eq.5.8) then
            n3=n3+1
            xo3(n3)=x5
            xi3(n3)=x3/x5
            xs3(n3)=x4/x5
         endif
         if(x6.eq.6.5) then
            n4=n4+1
            xo4(n4)=x5
            xi4(n4)=x3/x5
            xs4(n4)=x4/x5
         endif
      enddo
 666  continue
      close(1)

      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.9)
      call pgslw(3)

      xmin=6.0
      xmax=33.
      ymin=0.6
      ymax=1.2
      
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pgsci(1)
      call pgline(n1,xo1,xi1)
      call pgsci(2)
      call pgline(n2,xo2,xi2)
      call pgsci(3)
      call pgline(n3,xo3,xi3)
      call pgsci(4)
      call pgline(n4,xo4,xi4)
      call pgsci(1)
      call pglabel("f\Doutput\U(10\U-17\Dergs/s/cm\U2\D)"
     $     ,"f\Dinput\U/f\Doutput","")

      call pgsci(1)
      call pgptxt(22.,1.12,0.,0.,"S/N=4.9")
      call pgsci(2)
      call pgptxt(22.,1.07,0.,0.,"S/N=5.2")
      call pgsci(3)
      call pgptxt(22.,1.02,0.,0.,"S/N=5.8")
      call pgsci(4)
      call pgptxt(22.,0.97,0.,0.,"S/N=6.5")
      call pgsci(1)

      xmin=6.0
      xmax=33.
      ymin=0.2
      ymax=0.35
      
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pgsci(1)
      call pgline(n1,xo1,xs1)
      call pgsci(2)
      call pgline(n2,xo2,xs2)
      call pgsci(3)
      call pgline(n3,xo3,xs3)
      call pgsci(4)
      call pgline(n4,xo4,xs4)
      call pgsci(1)
      call pglabel("f\Doutput\U(10\U-17\Dergs/s/cm\U2\D)",
     $     "sigma/f\Doutput","")

      call pgend

      end
