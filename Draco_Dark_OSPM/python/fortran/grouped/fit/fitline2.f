
      parameter(nmax=8000)
      real x(nmax),y(nmax),sig(nmax),xsig(nmax),ysig(nmax)
      real xin(nmax),yin(nmax),diff(nmax)
      integer icol(nmax)
      character c1*40,c2*40

      ifit=1
      c2='in'
      open(unit=1,file=c2,status='old')
      
      ilog=0
      xmin=0.5
      xmax=5.
      ymin=4.
      ymax=35.
      
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      xoff=8.5
      xoff=0.
      yoff=0.
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,i3
         n=n+1
         x(n)=x1
         xin(n)=x(n)-xoff
         y(n)=x2
         yin(n)=y(n)
         icol(n)=i3
         xsig(n)=0.1
         ysig(n)=0.5
      enddo
 666  continue
      close(1)

      if(ilog.eq.0) then
         call pgenv(xmin,xmax,ymin,ymax,0,0)
      else
         xmin=log10(xmin)
         xmax=log10(xmax)
         call pgenv(xmin,xmax,ymin,ymax,0,10)
      endif
      call pgsls(1)
      call pgsch(0.6)
      do i=1,n
         call pgsci(icol(i))
         xp=x(i)
         call pgpt1(xp,y(i),17)
      enddo
      call pgsci(1)

      if(ifit.eq.1) then
         na=50
         as=2.0
         ae=3.0
         nb=50
         bs=4.
         be=6.
         nc=7
         cs=1.8
         ce=2.5
         fmin=1.e10
         do ia=1,na
            at=as+(ae-as)*float(ia-1)/float(na-1)
            do ib=1,nb
               bt=bs+(be-bs)*float(ib-1)/float(nb-1)
               do ic=1,nc
                  ct=cs+(ce-cs)*float(ic-1)/float(nc-1)
                  do i=1,n
c                     yv=at+bt*x(i)+ct/25.*x(i)*x(i)
                     yv=at+bt*x(i)+ct/5./5./5.*x(i)*x(i)*x(i)
                     diff(i)=abs(yv-y(i))
                  enddo
                  call biwgt(diff,n,xb,xs)
                  if(xb.lt.fmin) then
                     fmin=xb
                     abest=at
                     bbest=bt
                     cbest=ct
                  endif
               enddo
            enddo
         enddo
         a=abest
         b=bbest
         c=cbest

         open(unit=11,file='fitout',status='unknown')
         write(11,*) a,b,c
         close(11)

         do i=1,100
            x(i)=xmin+(xmax-xmin)*float(i-1)/float(100-1)
c            y(i)=a+b*x(i)+c/25.*x(i)*x(i)
            y(i)=a+b*x(i)+c/5./5./5.*x(i)*x(i)*x(i)
         enddo
         call pgslw(7)
         call pgsci(2)
         call pgline(100,x,y)
      endif

 865  continue

      call pgslw(3)
      call pgsci(1)
      call pgsch(1.5)
      open(unit=11,file='labels.dat',status='old',err=866)
      read(11,*,err=866,end=866) c1
      call pgmtxt('B',2.5,0.5,0.5,c1)
      read(11,*,err=866,end=866) c1
      call pgmtxt('L',2.0,0.5,0.5,c1)
      read(11,*,err=866,end=866) c1
      call pgmtxt('T',1.5,0.5,0.5,c1)
 866  close(11)

      call pgend
      end
