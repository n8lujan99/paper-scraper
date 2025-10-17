
      parameter (narrm1=2000,narrm2=2000)
      real xd(narrm1,narrm2),w(1036),wc(1036)
      real xd2(narrm1,narrm2),rnew(1036),xin(10000)
      real rin(10000),win(10000)
      integer naxes(2)
      character file1*180,cifu*3,camp*2,cspec*3
      logical simple,extend,anyf


      read *,cifu,cspec

      camp="LL"
      file1="/scratch/03946/hetdex/lib_calib/202106/i"//cifu//
     $     "a"//camp//"wave.dat"
      open(unit=1,file=file1,status='old')
      n=0
      do i=1,1036
         read(1,*,end=667) x1,x2
         n=n+1
         w(n)=x1
         wc(n)=x1+x2
      enddo
 667  continue
      close(1)

      file1="/scratch/03946/hetdex/lib_calib/202106/i"//cifu//
     $     "a"//camp//"cbwt.fits"
      iext=1

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         call ftclos(im1,ier)
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      ier=0
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      print *,naxes(1),naxes(2)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      file1="rres"//cspec//camp//".fits"
      iext=1

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         call ftclos(im1,ier)
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      ier=0
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      print *,naxes(1),naxes(2)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,narrm1,ncol,nrow,xd2,anyf,ier)
      call ftclos(im1,ier)

      do iw=1,1036
         nin=0
         do j=1,nrow
            do i=1,ncol
               win(i)=xd(i,j)+wc(iw)
               rin(i)=xd2(i,j)
c               print *,iw,i,j,xd(i,j),xd2(i,j),wc(iw)
            enddo
            call xlinint(w(iw),ncol,win,rin,r0)
            nin=nin+1
            xin(nin)=r0
         enddo
         call biwgt(xin,nin,xb,xs)
         rnew(iw)=xb
         print *,w(iw),rnew(iw),xs,nin
      enddo
            

 706  continue
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
