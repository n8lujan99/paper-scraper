
      parameter (narrm=100,narrm3=2060)
      real xd(narrm,narrm,narrm3),xr(narrm3),xp(narrm3)
      real xda(5,narrm,narrm,narrm3),spec(10,narrm3),xin(narrm3)
      real xba(10),xgrid(1000),ygrid(1000)
      integer naxes(3)
      integer iarr1(5,1000),iarr2(5,1000),iarr3(5,1000),iarr4(5,1000)
      character file1*80,a1*9,aname(5)*9,aout*13
      logical simple,extend,anyf

      open(unit=1,file='cen.dat',status='old')
      nf=0
      do i=1,5
         read(1,*) a1,x2,x3
         nf=nf+1
         aname(nf)=a1
         xoff=30.-x2-1.
         yoff=30.-x3-1.
         n=0
         do ix=1,20
            ixs=nint(float(ix)*3.-2.-xoff)
            ixe=ixs+2
            do iy=1,20
               n=n+1
               iys=nint(float(iy)*2.-2.-yoff)
               iye=iys+2
               iarr1(nf,n)=ixs
               iarr2(nf,n)=ixe
               iarr3(nf,n)=iye
               iarr4(nf,n)=iye
               if(i.eq.1) then
                  x1=float(ix)*3.-2.
                  x2=x1+2.
                  y1=float(iy)*3.-2.
                  y2=y1+2.
                  xgrid(n)=(x1+x2)/2.
                  ygrid(n)=(y1+y2)/2.
               endif
            enddo
         enddo
c         print *,n,iarr1(i,66),iarr2(i,66)
      enddo
 666  continue
      close(1)
      
      do i=1,nf
         file1=aname(i)
         iext=2

         im1=0
         ier=0
         call ftgiou(im1,ier)
         iread=0
         call ftopen(im1,file1,iread,iblock,ier)
         call ftmahd(im1,iext,ihd,ier)
         call ftghpr(im1,3,simple,ibit,naxis,naxes,ipc,igc,extend,ier)

         ncol=naxes(1)
         nrow=naxes(2)
         num=naxes(3)
         call ftg3de(im1,igc,0.,narrm,narrm,ncol,nrow,num,xd,anyf,ier)
         call ftclos(im1,ier)
         do i1=1,ncol
            do i2=1,nrow
               do i3=1,num
                  xda(i,i1,i2,i3)=xd(i1,i2,i3)
               enddo
            enddo
         enddo
      enddo

      naxis=1
      naxes(1)=num
      im1=0
      ier=0
      iblock=1
      igc=0
      aout="spec    .fits"
      open(unit=11,file='out',status='unknown')
      do i=1,n
         iout=1000+i
         write(aout(5:8),2001) iout
         nt=0
         do ia=1,nf
            ixs=iarr1(ia,i)
            ixe=iarr2(ia,i)
            iys=iarr3(ia,i)
            iye=iarr4(ia,i)
            if(ixs.ge.2.and.ixe.le.59.and.
     $           iys.ge.2.and.iye.le.59) then
               nt=nt+1
               do i3=1,num
                  sum=0.
                  do i1=ixs,ixe
                     do i2=iys,iye
                        sum=sum+xda(ia,i1,i2,i3)
                     enddo
                  enddo
                  spec(nt,i3)=sum
               enddo
            endif
         enddo
c- now normalize and then take biweight
         do it=1,nt
            nin=0
            do j=300,1000
               nin=nin+1
               xin(nin)=spec(it,j)
            enddo
            call biwgt(xin,nin,xb,xs)
            xba(it)=xb
            print *,i,it,xb,xs
         enddo
         do j=1,num
            do it=1,nt
               xin(it)=spec(it,j)/xba(it)
            enddo
            call biwgt(xin,nt,xb,xs)
            xp(j)=xb
         enddo

         call ftinit(50,aout,iblock,ier)
         call ftphps(50,-32,naxis,naxes,ier)
         call ftp2de(50,igc,narrm3,naxes(1),1,xp,ier)
         call ftclos(50,ier)

         write(11,1101) aout,xgrid(i),ygrid(i),nt

      enddo
      close(11)

 706  continue
 1101 format(a13,1x,f7.2,1x,f7.2,1x,i4)
 2001 format(i4)

      end
