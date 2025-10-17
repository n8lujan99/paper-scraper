
      parameter(nmax=100)
      real xh0(nmax),xh1(nmax),xh2(nmax),xl(2),yl(2)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(3)

      xmin=0.0
      xmax=2.0
      ymin=0.
      ymax=1.2
      
      call pgenv(xmin,xmax,ymin,ymax,0,0)
c      call pglabel("f\Dinput\U(10\U-17\Dergs/s/cm\U2\D)"
c     $     ,"Relative Number of Objects","")
      call pglabel("f\Dinput\U/f\Doutput"
     $     ,"Relative Number of Objects","")

      xl(1)=1.
      xl(2)=1.
      yl(1)=ymin
      yl(2)=ymax
      call pgsls(4)
      call pgline(2,xl,yl)
      call pgsls(1)

      call pgslw(6)
      open(unit=1,file='out',status='old')
      do j=1,100
c      do j=1,1
         read(1,*,end=666) i1,x2,x3,x4,i5,x6
         xp=x2
         xpo=x4
         nh=0
         sum=0.
         do i=1,i5
            read(1,*) x1,x2,x3
            nh=nh+1
c            xh0(nh)=x1
            xh0(nh)=x1/xpo
            xh1(nh)=x2
            xh2(nh)=x3
            sum=sum+xh2(nh)
         enddo
         print *,sum,xh0(2)-xh0(1)
         sum=sum*(xh0(2)-xh0(1))*1.5
         do i=1,nh
            xh2(i)=xh2(i)/sum
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
c         call pgline(2,xl,yl)
         call pgsls(1)
      enddo
 666  continue
      close(1)

      call pgsci(1)
      call pgptxt(1.5,1.1,0.,0.,"Flux\Dout\U")
      call pgptxt(1.5,1.00,0.,0.,"4.8-5.0")
      call pgsci(2)
      call pgptxt(1.5,0.92,0.,0.,"6.0-7.0")
      call pgsci(3)
      call pgptxt(1.5,0.84,0.,0.,"7.0-9.0")
      call pgsci(4)
      call pgptxt(1.5,0.76,0.,0.,"11-15")

      call pgend

      end
