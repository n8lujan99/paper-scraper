
      parameter(nmax=10000)
      real x(nmax),y(nmax),xn(nmax),yn(nmax),yw2(nmax)
      real xb(nmax),yb(nmax),yb2(nmax),yn2(nmax),y2(nmax)
      real xw(nmax),yw(nmax),ywa(nmax,1000),xin(nmax)
      character file1*12,file2*80,c1*18,cname*25

      xlya=1215.666

      idumb=1 ! 0 is normal, 1 if per A
      nf=18
      ibin=7
      ibin=1
      ib1=(ibin-1)/2
      xib=float(ibin)

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xwh=85.
      xmin=xlya-xwh
      xmax=xlya+xwh

      xminl=xlya-2.5
      xmaxl=xlya+2.5

      xwb=0.5
c      xwb=1.5
      nw=nint((xwh*2.)/xwb)
      do i=1,nw
         xw(i)=xmin+xwb*float(i-1)
      enddo

      open(unit=1,file='list',status='old')

      nl=0
      do il=1,20000
         read(1,*,end=666) file1,wave
         file2="spec/"//file1//".fits"
         z=wave/xlya-1.
         call getfits(file2,n,x,y)
         suml=0.
         do j=1,n
            x(j)=x(j)/(1.+z)
            if(x(j).gt.xminl.and.x(j).lt.xmaxl) suml=suml+y(j)
         enddo
c         if(suml.gt.4000.) goto 866
         if(suml.gt.3000.) goto 866
c         if(suml.gt.2000.) goto 866
         nl=nl+1
         print *,file1
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

         do j=1,nw
            call xlinint(xw(j),nbb,xn,yn,y0)
            ywa(j,nl)=y0
         enddo

 866     continue
      enddo
 666  continue
      close(1)

      print *,nl
      do i=1,nw
         nin=0.
         sum=0.
         do j=1,nl
            if(ywa(i,j).ne.0) then
               nin=nin+1
               xin(nin)=ywa(i,j)
               sum=sum+ywa(i,j)
            endif
         enddo
         sum=sum/float(nin)
         call biwgt(xin,nin,xb0,xs0)
         yw(i)=xb0
c         yw(i)=sum
      enddo

      xmin=1175.
c      xmax=1300.
      xmax=1260.
      call pgenv(xmin,xmax,-15.,20.,0,0)
      call pgline(nw,xw,yw)

      do i=1,nw
         yw2(i)=5.
         if(xw(i).lt.1203.) yw2(i)=10.
         if(xw(i).gt.1232.) yw2(i)=10.
      enddo
      call pgsci(4)
      call pgslw(5)
c      call pgline(nw,xw,yw2)

      ibin=6
      ib1=(ibin-1)/2
      xib=float(ibin)
      nbb=0
      n=nw
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
            sum=sum+yw(is)
            nb=nb+1
            yb(nb)=yw(is)
            xb(nb)=xw(is)
         enddo
         call biwgt(yb,nb,xbb,xsb)
         yn(nbb)=sum/float(nb)
c         yn(nbb)=xbb
         call biwgt(xb,nb,xbb,xsb)
         xn(nbb)=xbb
      enddo
      call pgsci(2)
      call pgslw(9)
      call pgline(nbb,xn,yn)

      call pgend

      end

      subroutine getfits(file1,n,x,y)
      parameter(nmax=10000)
      real x(nmax),y(nmax)
      integer naxes(2)
      character file1*80
      logical simple,extend,anyf

      im1=0
      ier=0
      iext=2
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)

      n=3681
      do i=1,n
         call ftgcve(im1,1,i,1,1,0.,wave,anyf,ier)
         call ftgcve(im1,3,i,1,1,0.,flux,anyf,ier)
         x(i)=wave
         y(i)=flux
      enddo

 706  continue
      call ftclos(im1,ier)

      return
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
c      if(xp.le.x(1)) yp=y(1)
c      if(xp.ge.x(n)) yp=y(n)
      if(xp.le.x(1)) yp=0.
      if(xp.ge.x(n)) yp=0.
      return
      end
