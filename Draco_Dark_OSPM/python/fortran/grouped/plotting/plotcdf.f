
      parameter(nmax=1000000)
      real r(nmax),ri(nmax),di(nmax),ro(nmax),do(nmax),sn(nmax)
      real xnum(nmax),xrad(nmax),snl(10),snh(10)

      iseed=-1

      open(unit=1,file='j1s',status='old')

      n=0
      sum=0.
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,x5,x6
         n=n+1
         r(n)=x1
         ri(n)=x2
         di(n)=x3
         ro(n)=x4
         do(n)=x5
         sn(n)=x6
         sum=sum+di(n)
      enddo
 666  continue
      close(1)
      d0=sum/float(n)
      cosd=cos(d0/57.29578)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      call pgenv(0.,2.0,0.,1.,0,0)
      call pglabel("Radial Offset","CDF","")
      call pgslw(5)

      snl(1)=4.5
      snh(1)=5.0
      snl(2)=5.0
      snh(2)=6.0
      snl(3)=6.0
      snh(3)=9.0

      aerr0=0.35

      do j=1,3
         n2=0
         do i=1,n
            if(sn(i).ge.snl(j).and.sn(i).lt.snh(j)) then
               n2=n2+1
               xrad(n2)=r(i)
               aerr=aerr0*gasdev(iseed)
               xrad(n2)=sqrt(r(i)*r(i)+aerr*aerr)
            endif
         enddo
         call sort(n2,xrad)
         do i=1,n2
            xnum(i)=float(i)/float(n2)
         enddo
         call pgsci(j)
         call pgline(n2,xrad,xnum)
      enddo

      end
         
