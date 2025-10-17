Skipped
      real arr1(1036,112*12)
      real arr2(1036,112*12)
      real arr3(1036,112*12)
      character a1*12,file1*80,a2*6,a8*24,file2*120,a9*16,a10*5

      rad1=600.

      read *,ra,dec,rad2,wave0,wc

      icol=nint((wave0-3470.)/2.+1)
      icolh=nint((wave0+wc-3470.)/2.+1)
      icoll=nint((wave0-wc-3470.)/2.+1)

      open(unit=1,file="/data/00115/gebhardt/detect/radec.dat",
     $     status='old')

      cosd=cos(dec/57.3)

      open(unit=11,file='out',status='unknown')

      do i=1,100000
         read(1,*,end=666) a1,x2,x3,x4
         rad=sqrt( (cosd*(ra-x2))**2 + (dec-x3)**2 )
         rad=rad*3600.
         iopen=0
         if(rad.lt.rad1) then
            file1="/data/00115/gebhardt/detect/"//a1//
     $           "/dithall.use"
            open(unit=2,file=file1,status='old')
            do j=1,100000
               read(2,*,end=667) x1,x2,a2,x4,x5,x6,x7,a8,a9,a10
               dx=cosd*(ra-x1)*3600.
               dy=(dec-x2)*3600.
               rad=sqrt( (cosd*(ra-x1))**2 + (dec-x2)**2 )
               rad=rad*3600.
               if(rad.lt.rad2) then
                  if(iopen.eq.0) then
                     file2="/data/00115/gebhardt/calfits/"//a1(1:6)//
     $                    "/"//a1(1:12)//"_"//a8(7:17)//"_cal.fits"
                     print *,file2
                     call getfits(file2,1,arr1,ier)
                     call getfits(file2,2,arr2,ier)
                     call getfits(file2,5,arr3,ier)
                     iopen=1
                  endif
                  if(a8(19:20).eq."LL") iamp=1
                  if(a8(19:20).eq."LU") iamp=2
                  if(a8(19:20).eq."RL") iamp=3
                  if(a8(19:20).eq."RU") iamp=4
                  if(a10.eq."exp01") iexp=1
                  if(a10.eq."exp02") iexp=2
                  if(a10.eq."exp03") iexp=3
                  read(a8(22:24),2001) irowa
                  irow=(iamp-1)*112+(iexp-1)*448+irowa
c                  print *,iamp,iexp,irowa,icol,irow
                  call getsum(arr1,arr2,arr3,icol,icolh,icoll,
     $                 irow,sum1,sum2)
                  write(11,1101) dx,dy,sum1,sum2,ra,dec,
     $                 a1,a8(1:18),a8(19:20),a8(22:24),a10,icol,irow
               endif
            enddo
 667        continue
            close(2)
         endif
      enddo
 666  continue
 1101 format(6(1x,f10.5),1x,a12,1x,a17,1x,a2,1x,a3,1x,a5,2(1x,i4))
 2001 format(i3)
      close(1)
      close(11)

      end

      subroutine getsum(arr1,arr2,arr3,icol,icolh,icoll,irow,sum1,sum2)
      real arr1(1036,112*12)
      real arr2(1036,112*12)
      real arr3(1036,112*12)

      j=irow
      sum1=0.
      sum2=0.
      sume1=0.
      do i=icoll,icolh
         if(arr1(i,j).ne.0.and.arr2(i,j).ne.0.) then
            sum1=sum1+arr1(i,j)/arr2(i,j)/arr2(i,j)
            sum2=sum2+arr3(i,j)/arr2(i,j)/arr2(i,j)
            sume1=sume1+1./arr2(i,j)/arr2(i,j)
         endif
      enddo

      if(sume1.gt.0.) then
         sum1=sum1/sume1
         sum2=sum2/sume1
      else
         sum1=0.
         sum2=0.
      endif

      return
      end

      subroutine getfits(file1,iext,arr,ier)
      real arr(1036,112*12)
      integer naxes(2)
      character file1*120
      logical simple,extend,anyf

      ier=0

      im1=50
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)

      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1036,ncol,nrow,arr,anyf,ier)

      call ftclos(im1,ier)
 706  continue
      if(ier.ne.0) print *,"No file for: ",file1
      return
      end
