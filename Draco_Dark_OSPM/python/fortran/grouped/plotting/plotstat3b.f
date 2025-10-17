
      parameter(nmax=100)
      real xh0(nmax),xh1(nmax),xh2(nmax),xl(2),yl(2)
      real xin(nmax),yin(nmax),xin2(nmax),yin2(nmax)

      fac1o=1.28
      fac2o=1.85

c      fac1=1.35
c      fac2=2.05

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(3)

      xmin=0.0
      xmax=23.
      ymin=0.
      ymax=0.38
      
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("f\Dinput\U(10\U-17\Dergs/s/cm\U2\D)"
     $     ,"Relative Number of Objects","")

      call pgslw(6)
      open(unit=1,file='out',status='old')
      do j=1,100
c      do j=1,1
         read(1,*,end=666) i1,x2,x3,x4,i5,x6

         fac1=fac1o
         fac2=fac2o
         fsn=min(2.,1.+0.2*x6/4.9)
         ffl=1.-0.035*(x4/4.9)**2
         print *,fsn,ffl,fsn*ffl
         fac1=fac1*fsn*ffl
         fac2=fac2*fsn*ffl

         xp=x2
         xpo=x4
         nh=0
         do i=1,i5
            read(1,*) x1,x2,x3
            nh=nh+1
            xh0(nh)=x1
            xh1(nh)=x2
            xh2(nh)=x3
         enddo
         call pgsci(j)
         call pgline(nh,xh0,xh2)
         if(j.eq.1) then
            do i=1,nh
               xin(i)=xh0(i)
               yin(i)=xh2(i)
            enddo
         endif
      enddo
 666  continue
      close(1)
      
      call pgsls(4)
      do i=1,nh
         xin2(i)=xin(i)*fac1
         yin2(i)=yin(i)/fac1
      enddo
      call pgsci(2)
      call pgline(nh,xin2,yin2)
      do i=1,nh
         xin2(i)=xin(i)*fac2
         yin2(i)=yin(i)/fac2
      enddo
      call pgsci(3)
      call pgline(nh,xin2,yin2)
      call pgsls(1)

      call pgsci(1)
c      call pgptxt(0.7,0.32,0.,0.,"Flux\Dout\U")
c      call pgptxt(0.7,0.29,0.,0.,"4.8-5.0")
c      call pgsci(2)
c      call pgptxt(0.7,0.27,0.,0.,"6.0-7.0")
c      call pgsci(3)
c      call pgptxt(0.7,0.25,0.,0.,"7.0-9.0")
c      call pgsci(4)
c      call pgptxt(0.7,0.23,0.,0.,"11-15")

      call pgend

      end
