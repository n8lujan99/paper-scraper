
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),suma(200),xsum(200)
      real xb(nmax),yb(nmax),xin(20),xine(20),yin(nmax,200)
      real xisum(200),xin1(nmax),xin2(nmax),yin2(nmax,200)
      real siga(nmax),xisig(nmax),xin3(nmax),xin2e(nmax)
      character cdone(nmax)*10
      character file1*80,file2*80,c1*18,afib(nmax)*12,cname*10

      idumb=1 ! 0 is normal, 1 if per A

      wavec=50.
      wc=5.
      nf=18
      ibin=7
c      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=3500.
      xmax=5500.

      open(unit=1,file='splista',status='old')

      nl=0
      ic=0
      nfibt=0
      ymin=1e10
      ymax=-ymin
      do il=1,20000
         read(1,*,end=666) file1,wave
         xmin=wave-wavec
         xmax=wave+wavec
         open(unit=2,file=file1,status='old',err=888)
         goto 889
 888     continue
         print *,"Does not exist: ",file1
         goto 667
 889     continue
         read(2,*) ifib
         nfibt=nfibt+ifib
         n=0
         x10=0.
         do i=1,nmax
            read(2,*,end=667) x1,(xin(j),xine(j),j=1,ifib)
            if(x1.gt.xmin.and.x1.lt.xmax) then
               n=n+1
               x(n)=x1
               do j=1,ifib
                  ntmp=nfibt-ifib+j
                  yin(n,ntmp)=xin(j)
                  afib(ntmp)=file1(8:19)
                  ymin=min(ymin,xin(j))
                  ymax=max(ymax,xin(j))
               enddo
            endif
         enddo
 667     continue
         close(2)
 966     continue
      enddo
 666  continue
      close(1)

      ic=0
      xcmin=wave-wc
      xcmax=wave+wc
      do j=1,nfibt
         sum=0.
         sumf=0.
         nsum=0
         nsumf=0
         do k=1,n
            y(k)=yin(k,j)
            if(x(k).gt.xcmin.and.x(k).lt.xcmax) then
               weight=1.-abs(x(k)-wave)/wc
c               sum=sum+y(k)
               sum=sum+y(k)*weight
               nsum=nsum+1
               xin1(nsum)=y(k)
            else
               sumf=sumf+y(k)
               nsumf=nsumf+1
               xin2(nsumf)=y(k)
            endif
         enddo
         call biwgt(xin1,nsum,xb1,xs1)
         call biwgt(xin2,nsumf,xb2,xs2)
         suma(j)=sum
c         suma(j)=sum/float(nsum)-sumf/float(nsumf)
c         suma(j)=(xb1-xb2)/xs2
c         suma(j)=(sum-xb2)/xs2
         xsum(j)=suma(j)
         xisum(j)=float(j)
         siga(j)=xs2
         xin3(j)=xs2
         xisig(j)=float(j)
      enddo

      call biwgt(xin3,nfibt,xbsig,xssig)
c      xscut=xbsig+4.*xssig
      xscut=xbsig+30.*xssig

      call sort2(nfibt,xsum,xisum)

      np=7
      np=13

      call biwgt(xsum,nfibt,xbs,xss)

      call pgsls(1)
      call pgslw(2)
      call pgsci(1)
      ybit=(ymax-ymin)/30.
      ymin=ymin-ybit
      ymax=ymax+ybit
      ymax=min(2.7,ymax)
      if(idumb.eq.1) ymax=min(0.8,ymax)
      ymin=max(-0.7,ymin)
c      ymax=min(6.7,ymax)
c      ymin=max(-6.7,ymin)
c      ymax=1.

      call pgsch(1.5)
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      if(idumb.eq.0) then
         call pglabel('Wavelength (\(2078))',
     $        '10\U-17\D ergs/cm\U2\D/s','')
      else
         call pglabel('Wavelength (\(2078))',
     $        '10\U-17\D ergs/cm\U2\D/s/\(2078)','')
      endif
      call pgsch(1.8)
      call pgmtxt('T',1.,0.5,0.5,file1(22:31))
      call pgsch(1.5)
c      call pgslw(4)
      call pgslw(2)
      ic=0
      do j=1,np
         nft=nint(xisum(nfibt-j+1))
         do k=1,n
            y(k)=yin(k,nft)
            if(idumb.eq.1) y(k)=y(k)/2.
            yin2(j,k)=y(k)
         enddo
         ic=ic+1
         if(ic.eq.7) ic=ic+1
         if(ic.gt.13) ic=1
         call pgsci(ic)
         call pgsls(1)
         if(siga(nft).gt.xscut) then
            call pgsls(4)
            call pgslw(1)
         endif
         call pgline(n,x,y)
         call pgslw(2)
         call pgsls(1)
      enddo
      do k=1,n
         sum=0.
         nsum1=0
         do j=1,np
            nft=nint(xisum(nfibt-j+1))
            xin1(j)=yin2(j,k)
            if(siga(nft).lt.xscut) then
               sum=sum+yin2(j,k)
               nsum1=nsum1+1
            endif
         enddo
         call biwgt(xin1,np,xb1,xs1)
         xin2(k)=sum/float(nsum1)
         xin2e(k)=xs1/sqrt(float(nsum1))
      enddo
      call pgslw(9)
      call pgsci(1)
      call pgline(n,x,xin2)
      call getstat(n,x,xin2,wave,sig1,xs1)
      call pgslw(2)
      call pgsci(1)

      open(unit=11,file='outspec',status='unknown')
      do i=1,n
         write(11,*) x(i),xin2(i),xin2e(i)
      enddo
      close(11)

      call pgend
 1001 format(i5,1x,a10,3(1x,f8.2),1x,i1,2(1x,f8.3))
 1002 format(a10,1x,f7.2,1x,f7.2,1x,i2,4(1x,f8.3))

      end

      subroutine getstat(n,x,y,w,sig,xs)
      real x(n),y(n),xin(10000)

      wc1=4.5
      wc2=8.0
      w1l=w-wc1
      w1h=w+wc1
      w2l=w-wc2
      w2h=w+wc2

      sum=0.
      ns=0
      nsum=0
      do i=1,n
         if(x(i).ge.w1l.and.x(i).le.w1h) then
            sum=sum+y(i)
            nsum=nsum+1
         endif
         if(x(i).le.w2l.or.x(i).ge.w2h) then
            ns=ns+1
            xin(ns)=y(i)
         endif
      enddo
      call biwgt(xin,ns,xb,xs)
      if(nsum.gt.0) then
         sum=sum/float(nsum)
         xs=xs/sqrt(float(nsum))
         sig=(sum-xb)/xs
      else
         sum=0.
         sig=0.
      endif
      
      return
      end
