
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax)
      real xb(nmax),yb(nmax)
      character file1*80,file2*80,c1*20,ctitle*15

      idumb=1 ! 0 is normal, 1 is per A
      wavec=40.
c      wavec=130.
      nf=18
      ibin=3
      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      ctitle="            "
      open(unit=1,file='title',status='old',err=111)
      read(1,*) ctitle
 111  close(1)


      call pgbegin(0,'?',4,4)
c      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=3500.
      xmax=5500.

      open(unit=1,file='splist',status='old')

      nl=0
      ic=0
      do il=1,20000
         read(1,*,end=666) file1,wave
c         read(1,*,end=666) file1,xnorm
         call getname(file1,iname)
         xmin=wave-wavec
         xmax=wave+wavec
         open(unit=2,file=file1,status='old',err=888)
         goto 889
 888     continue
         print *,"Does not exist: ",file1
         goto 866
 889     continue
         ymin=1e10
         ymax=-ymin
         n=0
         x10=0.
         ymin=1e10
         ymax=-ymin
         do i=1,nmax
            read(2,*,end=667) x1,x2
c            read(2,*,end=667) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
            if(x2.ne.0.and.x1.gt.xmin.and.x1.lt.xmax) then
               n=n+1
               x(n)=x1
               y(n)=x2
c               y(n)=x2*xnorm
               if(idumb.eq.1) y(n)=y(n)/2.
c               ymin=min(ymin,y(n))
c               ymax=max(ymax,y(n))
            endif
         enddo
 667     continue
         close(2)
         if(x10.gt.5.) goto 866
         nbb=0
         do j=1,n,ibin
            nbb=nbb+1
            istart=max(0,j-ib1)
            iend=istart+ibin-1
            if(iend.gt.n) then
               iend=n
               istart=n-ibin+1
            endif
            sum=0.
            nb=0
            do is=istart,iend
               sum=sum+y(is)
               nb=nb+1
               yb(nb)=y(is)
               xb(nb)=x(is)
            enddo
            call biwgt(yb,nb,xbb,xsb)
            yn(nbb)=xbb
            call biwgt(xb,nb,xbb,xsb)
            xn(nbb)=xbb
            ymin=min(ymin,yn(nbb))
            ymax=max(ymax,yn(nbb))
         enddo
         c1=file1(1:19)
         icp=1
         call pgsls(1)
         call pgslw(2)
         call pgsci(1)
c         call pgsch(1.8)
         call pgsch(1.2)
         ybit=(ymax-ymin)/15.
         ymin=ymin-ybit
         ymax=ymax+ybit
c         ymax=max(ymax,1.5)
c         ymin=-0.3
c         ymax=1.0

c         if(il.eq.1) then
            call pgenv(xmin,xmax,ymin,ymax,0,0)
            if(idumb.eq.0) then
               call pglabel('Wavelength (\(2078))',
     $              '10\U-17\D ergs/cm\U2\D/s','')
            else
               call pglabel('Wavelength (\(2078))',
     $              '10\U-17\D ergs/cm\U2\D/s/\(2078)',ctitle)
            endif
c         endif
         call pgslw(5)
c         call pgenv(xmin,xmax,ymin,ymax,0,0)
c         if(il.eq.1) call pgenv(xmin,xmax,0.,ymax,0,0)
c         call pgsci(il)
         ic=ic+1
         call pgsci(ic)
         if(ic.eq.13) ic=0
         call pgsci(1)
c         call pgline(n,x,y)
         call pgline(nbb,xn,yn)
         call pgsch(1.8)
c         call pgsci(2)

         do j=1,80
            if(file1(j:j).eq.".") then
               nf=j-1
               goto 555
            endif
         enddo
 555     continue
         call pgsch(2.3)
         call pgslw(1)
c         call pgmtxt('T',-1.15,0.5,0.5,file1(5:nf-8))
         call pgmtxt('T',-1.15,0.5,0.5,file1(21:iname))
         call pgsch(1.5)
         call pgsci(1)
 866     continue
      enddo
 666  continue
      close(1)

      call pgend

      end

      subroutine getname(f,ic)
      character f*60,c*15

      do i=1,60
         if(f(i:i).eq.".") then
            ic=i-1
            goto 666
         endif
      enddo
 666  continue

      return
      end
