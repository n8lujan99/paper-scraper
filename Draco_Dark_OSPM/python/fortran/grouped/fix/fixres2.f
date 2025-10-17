
      parameter (narrm0=36000,narrm=1040)
      real xd(narrm,narrm0),xp(narrm),wp(narrm)
      real ws(narrm0),xs(narrm0),ws2(narrm0),xs2(narrm0)
      real xres(narrm0),xd2(narrm,narrm0),xsuba(narrm)
      real xr(narrm,112),xin(narrm0*narrm)
      integer naxes(2),naxesr(2)
      character file1*100,c1*14,cs1(1000)*3,cs2(1000)*2
      logical simple,extend,anyf

      file1='in3'
      open(unit=1,file=file1,status='old')
      read(1,*)
      ns=0
      do i=1,narrm0
         read(1,*,end=666) c1
         ns=ns+1
         cs1(ns)=c1(1:3)
         cs2(ns)=c1(13:14)
      enddo
 666  continue
      close(1)

      file1='in.fits'
      iext=1

      im1=51
      ier=0
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      if(naxes(1).gt.narrm.or.naxes(2).gt.narrm0) then
         write(*,"('Arrays too small - make narrm bigger')")
         goto 706
      endif
      ncol=naxes(1)
      nrow=max(1,naxes(2))
      call ftg2de(im1,igc,0.,narrm,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      iall=0
      do i=1,ns
         file1="/scratch/00115/gebhardt/res/res"//cs1(i)//cs2(i)
     $        //".fits"
         im1=51
         call ftopen(im1,file1,iread,iblock,ier)
         if(ier.ne.0) then
            write(*,*) 'Error opening image : ',file1
         endif
         call ftmahd(im1,iext,ihd,ier)
         call ftghpr(im1,2,simple,ibit,naxis,naxesr,ipc,igc,extend,ier)
         ncolr=naxesr(1)
         nrowr=max(1,naxesr(2))
         call ftg2de(im1,igc,0.,narrm,ncolr,nrowr,xr,anyf,ier)
         call ftclos(im1,ier)
         print *,i,ier
         do irow=1,nrowr
            iall=iall+1
            do icol=1,ncolr
c               xd2(icol,iall)=xd(icol,iall)-xr(icol,irow)
               xd2(icol,iall)=xr(icol,irow)
            enddo
         enddo
      enddo

      fac1=8.
      fac2=20.
      nfac=10
      smin=1.e10
      do ifac=1,nfac
         fact=fac1+float(ifac-1)*(fac2-fac1)/float(nfac-1)
         nin=0
         do irow=1,nrow
            do icol=50,500
               nin=nin+1
               xin(nin)=xd(icol,irow)-fact*xd2(icol,irow)
            enddo
         enddo
         call biwgt(xin,nin,xbf,xsf)
         print *,ifac,fact,xbf,xsf
         if(xsf.lt.smin) then
            smin=xsf
            facmin=fact
         endif
      enddo

      print *,facmin
      fac=facmin
      nin=0
      do irow=1,nrow
         do icol=1,ncol
            xd2(icol,irow)=xd(icol,irow)-fac*xd2(icol,irow)
         enddo
      enddo

      call ftgiou(im1,ier)
      call ftinit(im1,'out2.fits',iblock,ier)
      call ftphps(im1,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
         goto 706
      endif
      call ftp2de(im1,igc,narrm,naxes(1),naxes(2),xd2,ier)
      call ftclos(im1,ier)
 
 706  continue
      end

      subroutine smcont(ns,ws,xs,ws2,xs2)
      parameter(nmax=36000)
      real ws(nmax),xs(nmax),ws2(nmax),xs2(nmax),win(nmax),xin(nmax)
      real xws(nmax),xbs(nmax)

      nstep=500
      n2=200
      xfrac=0.2
      nb=0
      do i=1,ns,nstep
         nb=nb+1
         do j=1,nstep
            win(j)=ws(i)
            xin(j)=xs(i)
         enddo
         call biwgt(xin,nstep,xb,xs)
         call biwgt(xin,n2,xb,xs)
         call biwgt(win,nstep,xbw,xs)
         xws(nb)=xbw
         xbs(nb)=xb
      enddo

      do i=1,ns
         call xlinint(ws(i),nb,xws,xbs,xval)
         xs2(i)=xs(i)-xval
         xs(i)=xs(i)/20.
         xs2(i)=xs2(i)/20.
         xs2(i)=max(xs2(i),0.)
      enddo

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

