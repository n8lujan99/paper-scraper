Marked
      real wave(100),fall(100,10000),xin(10000),tpn(10000)
      real rbc(10),xin2(10000)
      character file1*40

      open(unit=1,file='list',status='old')
      n=0
      do i=1,10000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         read(2,*,end=667) x1,x2,x3
         if(x3.le.0) goto 667
         rewind(2)
         n=n+1
         nw=0
         do j=1,100
            read(2,*,end=667) x1,x2,x3
            nw=nw+1
            fall(nw,n)=x2
            wave(nw)=x1
         enddo
 667     continue
         close(2)
      enddo
 666  continue
      close(1)

      wmin=4240.
      wmax=5140.

      do i=1,n
         nin=0
         do j=1,nw
            if(wave(j).ge.wmin.and.wave(j).le.wmax) then
               nin=nin+1
               xin(nin)=fall(j,i)
            endif
         enddo
         call biwgt(xin,nin,xb,xs)
         tpn(i)=xb
      enddo

c - Robin's coefficients
      rbc(1)=-0.01972573
      rbc(2)=+0.43904144
      rbc(3)=-3.86654830
      rbc(4)=+16.7814045
      rbc(5)=-35.7178268
      rbc(6)=+29.7167950

      open(unit=11,file='out',status='unknown')

      do j=1,nw
         do i=1,n
            xin(i)=fall(j,i)*0.14/tpn(i)
            xin2(i)=fall(j,i)
         enddo
         call biwgt(xin,n,xb,xs)
         call biwgt(xin2,n,xb2,xs2)
         xrbc=wave(j)/1000.
         fit=0.
         do ifit=1,6
            fit=(fit*xrbc)+rbc(ifit)
         enddo
         fit=fit*0.135/0.10

         write(11,1101) wave(j),xb,xs,xb/fit,n
c         write(*,1101) wave(j),xb,xs,xb/fit,n
         if(j.eq.1) xb1=xb
         if(j.eq.15) print *,xb/xb1,xb2
      enddo
      close(11)

 1101 format(1x,f6.1,3(2x,f8.5),2x,i5)

      end
