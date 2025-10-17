skipped
      parameter(nmax=1000000)
      real rall(10000),dall(10000),rout(nmax),dout(nmax)
      real rl(nmax),dl(nmax)
      real sum1a(1036),sum2a(1036),sume1a(1036),sumg(1036),sumge(1036)
      integer idout(nmax),idl(nmax)
      character a1*12,file1*80,a2*6,a8*24,file2*120,a9*16,a10*5,id2*12
      character aname(10000)*12,aout(nmax)*12,id2l(nmax)*12

      rad1=600.
      rad2=2.

      open(unit=1,file="/scratch/00115/gebhardt/detect/radec.dat",
     $     status='old')

      nall=0
      do iall=1,10000
         read(1,*,end=666) a1,x2,x3,x4
         nall=nall+1
         aname(nall)=a1
         rall(nall)=x2
         dall(nall)=x3
      enddo
 666  continue
      close(1)

      ncheck=100
      ntot=0
      open(unit=1,file='listin',status='old')
      nl=0
      do i=1,nmax
         read(1,*,end=888) ra,dec,id,id2
         nl=nl+1
         rl(nl)=ra
         dl(nl)=dec
         idl(nl)=id
         id2l(nl)=id2
      enddo
 888  continue
      close(1)

      print *,nall,nl

      open(unit=11,file='out',status='unknown')
      do il=1,nl
         ra=rl(il)
         dec=dl(il)
         id=idl(il)
         id2=id2l(il)

         if((float(il)/float(ncheck)-int(il/ncheck)).eq.0) then
            print *,il,ntot
         endif

         cosd=cos(dec/57.3)

         do iall=1,nall
            a1=aname(iall)
            if(a1.eq.id2) then
               ntot=ntot+1
               idout(ntot)=id
               rout(ntot)=ra
               dout(ntot)=dec
               aout(ntot)=a1
               goto 777
            endif
            x2=rall(iall)
            x3=dall(iall)
            rad=sqrt( (cosd*(ra-x2))**2 + (dec-x3)**2 )
            rad=rad*3600.
            if(rad.lt.rad1) then
               file1="/scratch/00115/gebhardt/detect/"//a1//
     $              "/dithall.use"
               open(unit=2,file=file1,status='old')
               do jall=1,100000
                  read(2,*,end=667) x1,x2,a2,x4,x5,x6,x7,a8,a9,a10
c                  dx=cosd*(ra-x1)*3600.
c                  dy=(dec-x2)*3600.
                  rad=sqrt( (cosd*(ra-x1))**2 + (dec-x2)**2 )
                  rad=rad*3600.
                  if(rad.lt.rad2) then
                     ntot=ntot+1
                     idout(ntot)=id
                     rout(ntot)=ra
                     dout(ntot)=dec
                     aout(ntot)=a1
                     write(11,*) id,ra,dec,a1
                     goto 667
                  endif
               enddo
 667           continue
               close(2)
            endif
 777        continue
         enddo
      enddo

c      do i=1,ntot
c      enddo
      close(11)

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
