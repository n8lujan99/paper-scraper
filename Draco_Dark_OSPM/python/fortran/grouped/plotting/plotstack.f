
      parameter(nmax=35000)
      real xa1(nmax),xa2(nmax),xa3(nmax),w(10)
      real xw1(nmax),xv1(nmax),xw2(nmax),xv2(nmax)
      character file1*80

      w(1)=1215.55
      w(2)=1215.60
      w(3)=1215.65
      w(4)=1215.70

      open(unit=1,file="j1",status='old')
      
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         n=n+1
         xa1(n)=x1
         xa2(n)=x2
      enddo
 666  continue
      close(1)

      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)

      xmin=10.
      xmax=25.
      ymin=-0.06
      ymax=0.07

      do iw=1,4
         wave0=w(iw)
         n1=0
         n2=0
         do i=1,n
            if(xa1(i).gt.wave0) then
               n1=n1+1
               xw1(n1)=xa1(i)-wave0
               xv1(n1)=xa2(i)
            endif
            if(xa1(i).le.wave0) then
               n2=n2+1
               xw2(n2)=wave0-xa1(i)
               xv2(n2)=xa2(i)
            endif
         enddo
         call pgsci(1)
         call pgslw(5)
         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pgsci(2)
         call pgline(n1,xw1,xv1)
         call pgsci(4)
         call pgline(n2,xw2,xv2)
      enddo

      call pgend

      end
