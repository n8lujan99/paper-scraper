
      parameter(nmax=100)
      real xh0(nmax),xh1(nmax),xh2(nmax),xl(2),yl(2)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(3)

      xmin=0.0
      xmax=16.
      ymin=0.
      ymax=0.34
      
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("f\Dinput\U(10\U-17\Dergs/s/cm\U2\D)"
     $     ,"Relative Number of Objects","")

      call pgslw(6)
      open(unit=1,file='out',status='old')
      do j=1,100
c      do j=1,1
         read(1,*,end=666) i1,x2,x3,x4,i5,x6
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
         xl(1)=xp
         xl(2)=xp
         yl(1)=ymin
         yl(2)=ymax
c         call pgline(2,xl,yl)
         xl(1)=xpo
         xl(2)=xpo
         call pgsls(4)
         call pgline(2,xl,yl)
         call pgsls(1)
      enddo
 666  continue
      close(1)

      call pgsci(1)
      call pgptxt(0.7,0.32,0.,0.,"Flux\Dout\U")
      call pgptxt(0.7,0.29,0.,0.,"4.8-5.0")
      call pgsci(2)
      call pgptxt(0.7,0.27,0.,0.,"6.0-7.0")
      call pgsci(3)
      call pgptxt(0.7,0.25,0.,0.,"7.0-9.0")
      call pgsci(4)
      call pgptxt(0.7,0.23,0.,0.,"11-15")

      call pgend

      end
