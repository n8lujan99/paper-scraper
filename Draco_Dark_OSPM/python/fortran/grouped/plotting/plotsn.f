
      parameter(nmax=300000)
      real sn(nmax),x(nmax),y(nmax)
      integer isn(nmax)

      xmin=4.5
      xmax=6.5
      ymin=0.
      ymax=0.35
      
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=1,file='sn.dat',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,i2
         n=n+1
         sn(n)=x1
         isn(n)=i2
      enddo
 666  continue
      close(1)

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("S/N","False positive rate","")

      ns=10
      do i=1,ns
         snp=xmin+(xmax-xmin)*float(i-1)/float(ns-1)
         x(i)=snp
         sum=0.
         sumt=0.
         do j=1,n
            if(sn(j).ge.x(i)) then
               sumt=sumt+1
               if(isn(j).eq.0) sum=sum+1.
            endif
         enddo
         y(i)=sum/sumt
         print *,x(i),y(i),sumt-sum
      enddo
      call pgslw(4)
      call pgline(ns,x,y)

      call pgend
      end
