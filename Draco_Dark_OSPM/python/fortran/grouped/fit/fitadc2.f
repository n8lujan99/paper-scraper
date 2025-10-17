Complete
      real wave(10),xadc(10000,1000),xin(10000),adcf(10)
      real xadcy(10000,1000),adcfy(10),xly(10000),xhy(10000)
      real wadc(5),adc(5),xl(10000),xh(10000)
      character file1*40

      open(unit=1,file='listin',status='old')

      nt=0
      do j=1,10000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         read(2,*)
         n=0
         nt=nt+1
         do i=1,10
            read(2,*,end=667) x1,x2,x3,x4,x5
            n=n+1
            wave(n)=x1
            xadc(nt,n)=x2
            xadcy(nt,n)=x4
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

      do i=1,10
         do j=1,nt
            xin(j)=xadc(j,i)
         enddo
         call biwgt(xin,nt,xb,xs)
         adcf(i)=xb
         print *,wave(i),xb,xs
         xl(i)=xb-xs
         xh(i)=xb+xs

         do j=1,nt
            xin(j)=xadcy(j,i)
         enddo
         call biwgt(xin,nt,xb,xs)
         adcfy(i)=xb
         print *,wave(i),xb,xs
         xly(i)=xb-xs
         xhy(i)=xb+xs
      enddo

      nadc=5
      wadc(1)=3500.
      wadc(2)=4000.
      wadc(3)=4500.
      wadc(4)=5000.
      wadc(5)=5500.
      adc(1)=-0.71
      adc(2)=-0.34
      adc(3)=-0.085
      adc(4)=0.08
      adc(5)=0.20

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      call pgenv(3400.,5600.,-0.9,0.5,0,0)
      call pglabel('Wavelength (\(2078))','ADC (arcsec)','2020-2021')
      call pgslw(5)
      call pgline(5,wadc,adc)
      call pgsci(2)
      call pgline(10,wave,adcf)
      call pgsls(4)
      call pgslw(2)
      call pgline(10,wave,xl)
      call pgline(10,wave,xh)

      call pgsls(1)
      call pgsci(4)
      call pgline(10,wave,adcfy)
      call pgsls(4)
      call pgslw(2)
      call pgline(10,wave,xly)
      call pgline(10,wave,xhy)

      call pgend



      end
