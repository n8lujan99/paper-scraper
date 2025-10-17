Skipped
      real arr1(1032,1700),xprof(1032,112,15),xa(1032*17000)
      real arr2(1032,1700),xprof2(1032,112,15)
      real xin(100000),yin(100000),xiv(1000)
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      read *,file1
      iext=1

      call getfits(file1,arr1,ier1)
      write(file1(7:7),1001) "U"
      call getfits(file1,arr2,ier2)

      do i=1,1032
         do j=1,112
            do k=1,15
               jp=(j-1)*15+k
               xprof(i,j,k)=arr1(i,jp)
            enddo
         enddo
      enddo

      do i=1,1032
         do j=1,112
            do k=1,15
               jp=(j-1)*15+k
               xprof2(i,j,k)=arr2(i,jp)
            enddo
         enddo
      enddo

      np=0
      nt=40
      xmax=1010.
      xmin=20.
      do i=1,nt
         xp=xmin+(xmax-xmin)*float(i-1)/float(nt-1)
         ixp=nint(xp)
         do j=1,112
            if(ier1.eq.0) then
               x1=xprof(ixp,j,7)
               x2=xprof(ixp,j,15)
               x3=xprof(ixp,j,1)
               x2=(x2+x3)/2.
               xnum=x1-2.*x2
               xden=x1+2.*x2
               if(xnum.ge.0.and.xnum.le.1) then
                  if(xden.gt.0.and.xden.le.1) then
                     xcon=xnum/xden
                     np=np+1
                     xa(np)=xcon
                  endif
               endif
            endif
            if(ier2.eq.0) then
               x1=xprof2(ixp,j,7)
               x2=xprof2(ixp,j,15)
               x3=xprof2(ixp,j,1)
               x2=(x2+x3)/2.
               xnum=x1-2.*x2
               xden=x1+2.*x2
               if(xnum.ge.0.and.xnum.le.1) then
                  if(xden.gt.0.and.xden.le.1) then
                     xcon=xnum/xden
                     np=np+1
                     xa(np)=xcon
                  endif
               endif
            endif
         enddo
      enddo
      print *,np

      call sort(np,xa)
      do i=1,np
         xin(i)=float(i)/float(np)
         yin(i)=xa(i)
      enddo

      ni=10
      xiv(1)=0.01
      xiv(2)=0.15
      xiv(3)=0.25
      xiv(4)=0.35
      xiv(5)=0.45
      xiv(6)=0.55
      xiv(7)=0.65
      xiv(8)=0.75
      xiv(9)=0.85
      xiv(10)=0.99

      open(unit=11,file='out',status='unknown')
      do i=1,ni
         call xlinint(xiv(i),np,xin,yin,yv)
         write(11,*) yv,xiv(i)
      enddo      
      close(11)

 706  continue
 1001 format(a1)
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end

      subroutine getfits(file1,arr,ier)
      real arr(1032,1700)
      integer naxes(2)
      character file1*100
      character camp*2,cspecid*3,cifu*3,cifupos*3
      logical simple,extend,anyf

      im1=0
      ier=0
      iext=1
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1032,ncol,nrow,arr,anyf,ier)

      call ftclos(im1,ier)
 706  continue
      if(ier.ne.0) print *,"No file for: ",file1
      return
      end
