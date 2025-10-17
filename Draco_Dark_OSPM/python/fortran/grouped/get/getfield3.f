Marked
      real arr1(1036,112*12)
      real arr2(1036,112*12)
      real arr3(1036,112*12)
      real sum1a(1036),sum2a(1036),sume1a(1036),sumg(1036),sumge(1036)
      character a1*12,file1*80,a2*6,a8*24,file2*120,a9*16,a10*5

      rad1=600.

      read *,ra,dec,rad2,wave0,wc,sumcut

      icol=nint((wave0-3470.)/2.+1)
      icolh=nint((wave0+wc-3470.)/2.+1)
      icoll=nint((wave0-wc-3470.)/2.+1)

      do i=1,1036
         sum1a(i)=0.
         sum2a(i)=0.
         sume1a(i)=0.
         sumg(i)=0.
         sumge(i)=0.
      enddo

      open(unit=1,file="/data/00115/gebhardt/detect/radec.dat",
     $     status='old')

      cosd=cos(dec/57.3)

      open(unit=11,file='out1',status='unknown')
      open(unit=12,file='out2',status='unknown')
      open(unit=13,file='out3',status='unknown')
      open(unit=14,file='out4',status='unknown')

      nsumg=0
      nsumt=0
      nid=0
      do iall=1,100000
         read(1,*,end=666) a1,x2,x3,x4
         rad=sqrt( (cosd*(ra-x2))**2 + (dec-x3)**2 )
         rad=rad*3600.
         iopen=0
         if(rad.lt.rad1) then
            file1="/data/00115/gebhardt/detect/"//a1//
     $           "/dithall.use"
            open(unit=2,file=file1,status='old')
            do jall=1,100000
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
                  nid=nid+1
c                  print *,iamp,iexp,irowa,icol,irow
                  call getsum(arr1,arr2,arr3,icol,icolh,icoll,
     $                 irow,sum1,sum2)
                  nsumt=nsumt+1
                  if(sum1.gt.sumcut) nsumg=nsumg+1
                  write(11,1101) dx,dy,sum1,sum2,x1,x2,
     $                 a1,a8(1:18),a8(19:20),a8(22:24),a10,icol,irow,
     $                 rad,nid
                  j=irow
                  do i=1,1036
                     wave=3470.+float(i-1)*2.
                     icol=nint((wave0-3470.)/2.+1)
                     if(arr1(i,j).ne.0.and.arr2(i,j).ne.0.) then
                        sum1a(i)=sum1a(i)+arr1(i,j)/arr2(i,j)/arr2(i,j)
                        sum2a(i)=sum2a(i)+arr3(i,j)/arr2(i,j)/arr2(i,j)
                        sume1a(i)=sume1a(i)+1./arr2(i,j)/arr2(i,j)
                        if(sum1.gt.sumcut) then
                           sumg(i)=sumg(i)+arr1(i,j)/arr2(i,j)/arr2(i,j)
                           sumge(i)=sumge(i)+1./arr2(i,j)/arr2(i,j)
                        endif
                     endif
                     write(12,1201) wave,arr1(i,j),
     $                    arr2(i,j),arr3(i,j),nid
                  enddo
               endif
            enddo
 667        continue
            close(2)
         endif
      enddo
 666  continue
      close(1)
      close(11)
      close(12)

      print *,nsumg,nsumt
      do i=1,1036
         wave=3470.+float(i-1)*2.
         if(sume1a(i).gt.0.) then
            sum1=sum1a(i)/sume1a(i)
            sum2=sum2a(i)/sume1a(i)
            sume1=sqrt(1./sume1a(i))
         else
            sum1=0.
            sum2=0.
            sume1=0.
         endif
         if(sumge(i).gt.0.) then
            sumg1=sumg(i)/sumge(i)
            sumge1=sqrt(1./sumge(i))
         else
            sumg1=0.
            sumge1=0.
         endif
c         write(13,1301) wave,sum1,sume1,sum2
         write(13,1301) wave,sumg1,sumge1,sum2
      enddo
      close(13)

 1101 format(6(1x,f10.5),1x,a12,1x,a17,1x,a2,1x,a3,1x,a5,2(1x,i4),
     $     1x,f10.3,1x,i5)
 1201 format(f7.2,1x,3(f10.3,1x),i5)
 1301 format(f7.2,3(1x,f10.4))
 2001 format(i3)
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
