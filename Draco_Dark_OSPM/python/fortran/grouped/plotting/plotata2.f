
      real x(1000),y(1000),xb(1000),yb(1000),xn(1000),yn(1000)
      character filei(4)*10,cfile*40,clab(4)*7

      ibin=11
      ib1=(ibin-1)/2
      xib=float(ibin)

      filei(1)='l1'
      filei(2)='l2'
      filei(3)='l3'
      filei(4)='l4'
      clab(1)='000-205'
      clab(2)='301-318'
      clab(3)='319-330'
      clab(4)='403-427'

      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)
      xmin=3500.
      xmax=5500.
      ymin=0.65
      ymax=1.45

      do ifile=1,4
         open(unit=1,file=filei(ifile),status='old')
         call pgsci(1)
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pglabel("Wavelength (\(2078))","Relative Throughput","")
         call pgmtxt('T',-1.5,0.2,0.5,clab(ifile))
         ic=1
         do iall=1,1000
            read(1,*,end=666) cfile
            if(cfile(6:7).eq.'LL') then
               open(unit=2,file=cfile,status='old')
               n=0
               do i=1,1000
                  read(2,*,end=667) x1,x2
                  if(x2.gt.0.6) then
                     n=n+1
                     x(n)=x1
                     y(n)=x2
                  endif
               enddo
 667           continue
               close(2)

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
               enddo
               ic=ic+1
               if(ic.eq.14) ic=1
               call pgsci(ic)
c               call pgline(n,x,y)
               call pgline(nbb,xn,yn)
            endif
         enddo
 666     continue

         close(1)
      enddo

      call pgend
      end
