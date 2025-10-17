
      parameter(nmax=18000000)
      real x(nmax),sn(nmax),xp(1000),yp(1000),xn(nmax)
      real snl(10),snh(10),yp1(1000)
      integer ic(10)
      character file1*80

 1    continue
c 1    write(*,"('Data file : '$)")
c     read *,file1
      file1="zhist"
      open(unit=1,file=file1,status='old',err=1)

c      read *,s0
      s0=1.5
      
      n=0
      xmin=1.9
      xmax=3.4
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4
c     s1=sqrt(x2*x2-s0*s0)
         s1=x2
c         s1=x2/1.005
c         s1=x2/1.1
         if(x3.le.3950.and.s1.le.4.8.or.x4.gt.8.) then
         else
            n=n+1
            x(n)=x1
            sn(n)=x2
         endif
      enddo
 666  continue
      close(1)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      call pgenv(xmin,xmax,0.,1.2,0,0)
      call pgslw(5)
      
      snl(1)=4.79
      snh(1)=5.2
      snl(2)=5.2
      snh(2)=5.7
      snl(3)=5.7
      snh(3)=10.
      ic(1)=1
      ic(2)=2
      ic(3)=4
      
      do iall=1,3
         nn=0
         do i=1,nmax
            if(sn(i).gt.snl(iall).and.sn(i).le.snh(iall)) then
               nn=nn+1
               xn(nn)=x(i)
            endif
         enddo

         nbin=100
         xbin=(xmax-xmin)/float(nbin-1)
         ymax=0.
         xtot=0.
         do i=1,nbin
            xlo=xmin+xbin*float(i-1)
            xhi=xlo+xbin
            nb=0
            do j=1,nn
               if(xn(j).ge.xlo.and.xn(j).lt.xhi) nb=nb+1
            enddo
c            print *,(xhi+xlo)/2.,nb
            xp(i)=(xhi+xlo)/2.
            yp(i)=float(nb)
            xtot=xtot+yp(i)
            ymax=max(ymax,yp(i))
         enddo
         yrat=0.
         do i=1,nbin
            yp(i)=yp(i)/ymax
            if(iall.eq.1) then
               yp1(i)=yp(i)
            else
               yrat=yrat+yp(i)/yp1(i)
            endif
         enddo
         yrat=yrat/float(nbin)
         print *,yrat,xtot
         if(iall.gt.1) then
            do i=1,nbin
               yp(i)=yp(i)/yrat
            enddo
         endif

         call pgsci(ic(iall))
         call pgline(nbin,xp,yp)
      enddo

      call pgend

      end
