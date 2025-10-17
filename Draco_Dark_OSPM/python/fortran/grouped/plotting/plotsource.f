
      parameter(nmax=300000)
      real sni(nmax),sns(nmax),x(nmax),y(nmax),y2(nmax)

      xmin=5.
      xmax=8.
      ymin=0.
      ymax=6.
      
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=1,file='inverse.dat',status='old')
      ni=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,x5,x6,x7,x8,x9
         ni=ni+1
         sni(ni)=x9
      enddo
 666  continue
      close(1)
      open(unit=1,file='sources.dat',status='old')
      ns=0
      do i=1,nmax
         read(1,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9
         ns=ns+1
         sns(ns)=x9
      enddo
 667  continue
      close(1)

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      print *,ni,ns
      n=50
      xifu=11306.
      do i=1,n
         sn=xmin+(xmax-xmin)*float(i-1)/float(n-1)
         x(i)=sn
         sum=0.
         do j=1,ns
            if(sns(j).gt.x(i)) sum=sum+1
         enddo
         sum=sum*0.9
         y(i)=sum/xifu
         sumi=0.
         do j=1,ni
            if(sni(j).gt.x(i)) sumi=sumi+1
         enddo
         sumi=sumi*0.8
         y2(i)=sumi/xifu
         print *,x(i),y(i),sumi/sum*4.1/2.5,sumi/sum
      enddo
      call pgline(n,x,y)
      call pgsci(4)
      call pgline(n,x,y2)

      x(1)=xmin
      x(2)=xmax
      y(1)=4.1
      y(2)=4.1
      call pgsci(2)
      call pgline(2,x,y)

      call pgend
      end
