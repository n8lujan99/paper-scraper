
      parameter (narrm1=1032,narrm2=1032)
      real xd(narrm1,narrm2),xp(1032),yp(1032)
      real*8 xin(narrm1*2)
      integer naxes(2)
      character file1*180,a1*18,a2*14
      logical simple,extend,anyf

c      read *,file1

      open(unit=1,file='list',status='old')
      open(unit=11,file='out',status='unknown')
      iext=16

      do iall=1,10000
         read(1,*,end=666) file1
         a1="d"//file1(38:45)//"s"//file1(62:64)//file1(66:70)
         a2=file1(84:97)
         ier=0
         iread=0

         call ftopen(51,file1,iread,iblock,ier)
         call ftmahd(51,iext,ihd,ier)
         call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         if(naxis.eq.1) naxes(2)=1
         ncol=naxes(1)
         nrow=naxes(2)
         call ftg2de(51,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
         call ftclos(51,ier)

c      call pgbegin(0,'?',1,1)
c      call pgpap(0.,1.)
c      call pgscf(2)
c      call pgsch(1.5)
c      call pgslw(2)

c      call pgenv(1.,1032.,0.,1e4)

         xhigh=2000
         istart=10
         iend=1020
c         istart=1
c         iend=35
         
         nhigh=0
         nc=0
         do j=1,nrow
            do i=1,ncol
               xin(i)=dble(xd(i,j))
               xin(i+ncol)=0.d0
            enddo
            call drealft(xin,ncol,1)
            do i=istart,iend
               k=2*i-1
               xp(i)=float(i)
               yp(i)=sngl(
     $              sqrt(xin(k)*xin(k)+xin(k+1)*xin(k+1)))
               if(yp(i).gt.xhigh) nhigh=nhigh+1
c               print *,i,yp(i)
            enddo

         enddo
         write(11,1101) a1,a2,nhigh
      enddo
 666  continue
      close(1)
      close(11)

 1101 format(a18,1x,a14,1x,i6)
      end
