
      parameter(nmax=10000)
      real xr(nmax),xd(nmax),dl(nmax),du(nmax)
      real rn(nmax),dn(nmax),pn(nmax)
      real*8 d1,d2
      integer ip(nmax)
      character a1*12,an(nmax)*12,a3*6,a8*28,a9*17,a10*5
      parameter(radtodeg=57.29578)

      open(unit=1,file='/scratch/00115/gebhardt/detect/coordfix.dat',
     $     status='old')
      nf=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,i5
         nf=nf+1
         xr(nf)=x1
         xd(nf)=x2
         dl(nf)=x3
         du(nf)=x4
         ip(nf)=i5
      enddo
 666  continue
      close(1)

      open(unit=1,file='/scratch/00115/gebhardt/detect/radec.dat',
     $     status='old')
      n=0
      do i=1,nmax
         read(1,*,end=667) a1,x2,x3,x4
         n=n+1
         an(n)=a1
         rn(n)=x2
         dn(n)=x3
         pn(n)=x4
      enddo
 667  continue
      close(1)

      read *,a1

      do i=1,n
         if(a1.eq.an(i)) then
            rp=rn(i)
            dp=dn(i)
            pp=pn(i)
            goto 777
         endif
      enddo
      print *,"Not Here! ",a1
 777  continue
      
      do i=1,nf
         if(dp.gt.dl(i).and.dp.lt.du(i)) then
            if(pp.gt.180.and.ip(i).eq.1) then
               xcor=xr(i)
               dcor=xd(i)
               goto 888
            endif
            if(pp.le.180.and.ip(i).eq.0) then
               xcor=xr(i)
               dcor=xd(i)
               goto 888
            endif
         endif
      enddo
      xcor=0.
      dcor=0.
 888  continue

      xcor=xcor/cos(dp/radtodeg)

      open(unit=11,file='out',status='unknown')
      nt=0
      open(unit=1,file=a1//"/dithall.use",status='old')
      do i=1,1000000
         read(1,*,end=669) d1,d2,a3,x4,x5,x6,x7,a8,a9,a10
         d1=d1-dble(xcor)/3600.d0
         d2=d2-dble(dcor)/3600.d0
         write(11,1101) sngl(d1),sngl(d2),a3,x4,x5,x6,x7,a8,a9,a10
         nt=nt+1
      enddo
 669  continue
      close(1)
      close(11)
      print *,a1," ",nt

c      write(*,1001) a1,xnew,ynew,x1,x2
 1101 format(1x,f11.7,1x,f11.7,3x,a6,4(2x,f8.3),2x,a28,2x,a17,5x,a5)
 1001 format(a12,4(1x,f8.4))
      end
