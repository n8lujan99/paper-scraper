      parameter(nmax=10000)
      real xa(100,nmax),x(100),yin(100),xin(nmax)
      real yl(100),yh(100)
      character file1*40

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      xmin=3500.
      xmax=5500.
      ymin=0.
      ymax=0.2

      xlo=4600.
      xhi=5100.
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("Wavelength (\(2078))","Relative Throughput","")
      
      call pgsch(0.7)

      open(unit=1,file='list',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) file1
         n=n+1
         open(unit=2,file=file1,status='unknown')
         sum=0.
         nsum=0
         do j=1,nmax
            read(2,*,end=667) x1,x2
            if(x1.gt.xlo.and.x1.lt.xhi) then
               sum=sum+x2
               nsum=nsum+1
            endif
         enddo
 667     rewind(2)
         sum=sum/float(nsum)
         nw=0
         do j=1,nmax
            read(2,*,end=668) x1,x2
            nw=nw+1
            x2=x2*0.15/sum
            call pgpt1(x1,x2,17)
            x(j)=x1
            xa(j,n)=x2
         enddo
 668     continue
         close(2)
      enddo
 666  continue
      close(1)

      do j=1,nw
         do i=1,n
            xin(i)=xa(j,i)
         enddo
         call biwgt(xin,n,xb,xs)
         yin(j)=xb
         yl(j)=xb-xs
         yh(j)=xb+xs
         print *,x(j),yin(j),yl(j),yh(j)
      enddo

      call pgsci(2)
      call pgsch(1.5)
      call pgslw(5)
      call pgline(nw,x,yin)
      call pgsls(4)
      call pgline(nw,x,yl)
      call pgline(nw,x,yh)

      call pgend
      end
