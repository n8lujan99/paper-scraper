
      parameter (nmax=2000,narrm=2000)
      real val(nmax,narrm),arr(narrm,narrm),xb(nmax),vala(narrm)
      real tarray(2),vala2(narrm),arrs(narrm,narrm),vala3(narrm)
      real vala4(narrm)
      real val1(nmax,narrm),val2(nmax,narrm),val3(nmax,narrm)
      real xb1(nmax),xb2(nmax),xb3(nmax),arr1(narrm,narrm)
      real xin(nmax*narrm),sky(1036)
      integer naxes(2)
      character file1*120,filea(nmax)*120,comm*120
      logical simple,extend,anyf

      open(unit=1,file="listin",status='old',err=706)
      read *,iamp

      tim=etime(tarray)
      n=0
      do i=1,2000
c         read(1,*,end=666) file1,i2
         read(1,*,end=666) file1
c         i2=8
         n=n+1
         filea(n)=file1
         if(i.eq.1) print *,file1
      enddo
 666  continue

      nuse=nint(0.2*float(n))
      print *,nuse
      print *

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,filea(1),iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',filea(1)
         goto 706
      endif
      call ftmahd(im1,5,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      call ftclos(im1,ier)
      if(naxis.eq.1) naxes(2)=1
      if(naxes(1).gt.narrm.or.naxes(2).gt.narrm) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=naxes(2)

      norig=n
      nexp=3
      if(iamp.eq.1) istart=1
      if(iamp.eq.2) istart=113
      if(iamp.eq.3) istart=225
      if(iamp.eq.4) istart=337
      iexp=1
      ir1=(iexp-1)*448+istart
      ir2=ir1+111
c      print *,ir1,ir2,iexp
      do irow=ir1,ir2
c     print *,"irow: ",irow
         irowa=irow-istart+1
         do iexp=1,nexp
            irowe=irow+(iexp-1)*448
            n=norig
            do iall=1,n
               i=(iexp-1)*n+iall
c               call getsky(filea(iall),iexp,sky)
               call ftfiou(-1,ier)
               call ftgiou(im2,ier)
               call ftopen(im2,filea(iall),0,iblock,ier)
c- this is ff residual
               call ftmahd(im2,5,ihd,ier)
c               ip=(irow-1)*ncol+1
               ip=(irowe-1)*ncol+1
               call ftgpve(im2,ig,ip,ncol,0.,vala,anyf,ier)
c- this is local res
               call ftmahd(im2,1,ihd,ier)
               call ftgpve(im2,ig,ip,ncol,0.,vala2,anyf,ier)
c- this is error
               call ftmahd(im2,2,ihd,ier)
               call ftgpve(im2,ig,ip,ncol,0.,vala3,anyf,ier)
c- this is count
               call ftmahd(im2,3,ihd,ier)
               call ftgpve(im2,ig,ip,ncol,0.,vala4,anyf,ier)
               call ftclos(im2,ier)
c- convert sky to flux and write into vala2
c               do icol=1,ncol
c                  if(vala4(icol).ne.0) then
c                     sky(icol)=sky(icol)*vala2(icol)/vala4(icol)
c                  else
c                     sky(icol)=0.
c                  endif                    
c               enddo
c               do icol=1,ncol
c                  vala2(icol)=sky(icol)
c               enddo
               do icol=1,ncol
                  if(vala(icol).ne.0.and.vala2(icol).gt.0) then
c                     val(i,icol)=vala(icol)/vala2(icol)
                     val(i,icol)=vala(icol)
                     val1(i,icol)=vala(icol)
                     val2(i,icol)=vala2(icol)
                     val3(i,icol)=vala3(icol)
                  else
                     val(i,icol)=0.
                     val1(i,icol)=0.
                     val2(i,icol)=0.
                     val3(i,icol)=0.
                  endif
               enddo
            enddo
         enddo
         n=norig*nexp

         do icol=1,ncol
            nin=0
            nin1=0
            nin2=0
            nin3=0
            do i=1,n
               if(val(i,icol).ne.0.) then
                  nin=nin+1
                  xb(nin)=val(i,icol)
               endif
            enddo
            do i=1,n
               if(val1(i,icol).ne.0.) then
                  nin1=nin1+1
                  xb1(nin1)=val1(i,icol)
               endif
            enddo
            do i=1,n
               if(val2(i,icol).ne.0.) then
                  nin2=nin2+1
                  xb2(nin2)=val2(i,icol)
               endif
            enddo
            do i=1,n
               if(val3(i,icol).ne.0.) then
                  nin3=nin3+1
                  xb3(nin3)=val3(i,icol)
               endif
            enddo
            call biwgt(xb,nin,xl,xs)
            call biwgt(xb1,nin1,xl1,xs1)
            call biwgt(xb2,nin2,xl2,xs2)
            call biwgt(xb3,nin3,xl3,xs3)
            if(nin.gt.nuse) then
               arr(icol,irowa)=xl
               arrs(icol,irowa)=xs
            else
               arr(icol,irowa)=0.
               arrs(icol,irowa)=0.
            endif
            if(xl3.gt.0.and.nin3.gt.nuse) then
               arr1(icol,irowa)=xs*xl2/xl3
            else
               arr1(icol,irowa)=1.
            endif
         enddo
      enddo

      nin=0
      do icol=10,1021
         do irow=2,111
            if(arr1(icol,irow).ne.1) then
               nin=nin+1
               xin(nin)=arr1(icol,irow)
            endif
         enddo
      enddo
      if(nin.gt.1000) then
         call biwgt(xin,nin,xl,xs)
         do icol=1,1032
            do irow=1,112
               if(arr1(icol,irow).ne.1) 
     $              arr1(icol,irow)=arr1(icol,irow)/xl
            enddo
         enddo
      endif

      naxes(2)=112
      call ftfiou(-1,ier)
      call ftinit(52,'imcmb.fits',iblock,ier)
      call ftphps(52,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftpkyj(52,"N",n,comm,ier)
      call ftp2de(52,igc,narrm,ncol,nrow,arr,ier)
      call ftclos(52,ier)

      call ftfiou(-1,ier)
      call ftinit(52,'imcmbs.fits',iblock,ier)
      call ftphps(52,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftpkyj(52,"N",n,comm,ier)
      call ftp2de(52,igc,narrm,ncol,nrow,arrs,ier)
      call ftclos(52,ier)

      call ftfiou(-1,ier)
      call ftinit(52,'imcmb1.fits',iblock,ier)
      call ftphps(52,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftpkyj(52,"N",n,"Number of files",ier)
      call ftpkye(52,"norm",xl,5,"Normalization",ier)
      call ftp2de(52,igc,narrm,ncol,nrow,arr1,ier)
      call ftclos(52,ier)

 706  continue
      end

      subroutine getsky(file1,iexp,sky)
      real sky(1036),w(100000),s(100000)
      character file1*120,date*8,shot*3,cexp*1,file2*26

      date=file1(1:8)
      shot=file1(10:12)
      write(cexp,1001) iexp
      file2="d"//date//"s"//shot//"exp0"//cexp//"sky.dat"

      open(unit=2,file=file2,status='old')
      n=0
      do i=1,100000
         read(2,*,end=666) x1,x2
         n=n+1
         w(n)=x1
         s(n)=x2
      enddo
 666  continue
      close(2)

      do i=1,1036
         ws=3470.+float(i-1)*2.
         call xlinint(ws,n,w,s,wo)
         sky(i)=wo
      enddo
      
 1001 format(i1)

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
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
