
      parameter (nmax=2000,narrm=2000)
      real xw(nmax),xs(nmax),xd(nmax,nmax),xin(nmax)
      real xday(10000),xmjd(10000)
      integer naxes(2)
      character file1*120,filea(10000)*120,comm*120
      character cspec*25,cdate*8,cshot*3
      logical simple,extend,anyf

      cspec="multi_312_045_065_LU.fits"
      open(unit=1,file="listin",status='old',err=706)

 1001 format(a8,1x,a3)
      n=0
      do i=1,10000
c         read(1,1001,end=666) cdate,cshot,x3
         read(1,*,end=666) cdate,cshot,x3
         file1="/data/00115/gebhardt/red1/reductions/"//cdate//
     $        "/virus/virus0000"//cshot//"/exp01/virus/"//cspec
         open(unit=2,file=file1,status='old',err=707)
         n=n+1
         filea(n)=file1
         if(n.eq.1) then
            mjd0=x3
            print *,cdate," ",cshot
         endif
         xday(n)=(x3-mjd0)/365.
         xmjd(n)=x3
c         print *,filea(n)
 707     continue
         close(2)
      enddo
 666  continue

      print *,n,xday(n),mjd0

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

c      xmin=5300.
c      xmax=5500.
c      ymin=0.5
c      ymax=3.
      xmin=0.
      xmax=xday(n)
      ymin=2.
      ymax=12.
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("Years since 201808","5461 Relative to Sky continuum"
     $     ,"")

      w1=5300.
      w2=5450.
      w3=5458.
      w4=5464.
      ier=0
      im1=50

      open(unit=11,file='out',status='unknown')

      do i=1,n
         ier=0
c         call ftgiou(im1,ier)
         call ftopen(im1,filea(i),0,iblock,ier)
         call ftgkye(im1,"TRAJCDEC",dec,file1,ier2)
         call ftgkye(im1,"STRUCTAZ",az,file1,ier2)
c- this is wavelength
         call ftmahd(im1,12,ihd,ier)
         call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=1032
         nrow=112
         if(ier.ne.0) goto 888
         call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd,anyf,ier)
         do j=1,ncol
            xw(j)=xd(j,66)
         enddo
c- this is sky
         call ftmahd(im1,17,ihd,ier)
         call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd,anyf,ier)
         do j=1,ncol
            xs(j)=xd(j,66)
         enddo
         call ftclos(im1,ier)

         nw=0
         sum=0.
         do j=1,ncol
            if(xw(j).gt.w1.and.xw(j).lt.w2) then
               nw=nw+1
               xin(nw)=xs(j)
            endif
            if(xw(j).gt.w3.and.xw(j).lt.w4) then
               sum=sum+xs(j)
            endif
         enddo
         call biwgt(xin,nw,xb,xs)
         if(xb.gt.100.) then
            do j=1,ncol
               xs(j)=xs(j)/xb
            enddo
            sum=sum/xb

c            call pgline(ncol,xw,xs)
            call pgpt1(xday(i),sum,17)
            if(dec.lt.-100) then
               dec=999.
               az=999.
            endif
            write(11,1101) xday(i),sum,xmjd(i),dec,az
         endif
c         read *
 888     continue
         call ftclos(im1,ier)
      enddo
      close(11)

 1101 format(f7.4,1x,f7.3,1x,f10.2,1x,f10.6,1x,f10.2)
 706  continue
      end
