
      parameter (narrm0=36000,narrm=1040)
      real xd(narrm,narrm0),xr(narrm),xp(narrm),wp(narrm)
      real ws(narrm0),xs(narrm0),ws2(narrm0),xs2(narrm0)
      real xres(narrm0),xin(narrm0),xd2(narrm,narrm0),xsuba(narrm)
      integer naxes(2)
      character file1*160
      logical simple,extend,anyf

      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      file1='in2'
      open(unit=1,file=file1,status='old')
      ns=0
      do i=1,narrm0
         read(1,*,end=666) x1,x2
         ns=ns+1
         ws(ns)=x1
         xs(ns)=x2
      enddo
 666  continue
      close(1)
      call smcont(ns,ws,xs,ws2,xs2)

      file1='in.fits'
      iext=1

      im1=51
      ier=0
      iread=0
c      call ftinit(im1,file1,iblock,ier)
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

      nrall=0
      irmax=nrow-111
      do ir=1,irmax,112
         ir1=ir
         
      ispec=1
c      ir1=8510
      ir2=ir1+111
      iflip=0
      isum=0

      suma=0.
      nsum=0
      do i=1,ncol
         nr=0
         sum=0.
         sum1=0.
         sum2=0.
         do j=ir1,ir2
            if(nint(xd(i,j)).ne.0) then
               nr=nr+1
               xr(nr)=xd(i,j)
               sum=sum+xr(nr)
               sum1=sum1+xd(i,j)*float(j)
               sum2=sum2+xd(i,j)
            endif
         enddo
         call biwgt(xr,nr,xb,xs)
         if(nr.gt.0) then
            sum=sum/float(nr)
            if(sum2.gt.0) then
               suma=suma+sum1/sum2
               nsum=nsum+1
            endif
         else
            sum=-666.
         endif   
         xp(i)=xb
         if(isum.eq.1) xp(i)=sum
         wp(i)=3470.+2.*float(i-1)
      enddo

      xblo=1.e10
      nuse=500
      do ifac=1,50
         xfac=0.+float(ifac)/40.
         do i=1,nuse
            call xlinint(wp(i),ns,ws,xs2,xval)
            xsub=max(0.,xfac*xval)
            xres(i)=xp(i)-xsub
            xin(i)=xres(i)*xres(i)
         enddo
         call biwgt(xin,nuse,xb,xs)
         if(xb.lt.xblo) then
            xblo=xb
            xfacl=xfac
         endif
      enddo
      xfac=xfacl
      do i=1,ncol
         call xlinint(wp(i),ns,ws,xs2,xval)
         xsub=max(0.,xfac*xval)
         xres(i)=xp(i)-xsub
         xsuba(i)=xsub
      enddo
      print *,xfac,ir1,ir2,nrall
      do j=1,112
         nrall=nrall+1
         do i=1,ncol
            xd2(i,nrall)=xd(i,nrall)-xsuba(i)
         enddo
      enddo

      call pgenv(3500.,5500.,-20.,20.,0,0)
      call pgsci(2)
      call pgline(ns,ws,xs)
      call pgsci(3)
c      call pgline(ns,ws,xs2)
      call pgsci(1)
      call pgline(ncol,wp,xp)
      call pgsci(3)
      call pgline(ncol,wp,xres)
      call pgsci(1)
      enddo
      call pgend

      ier=0
      call ftclos(im1,ier)
      call ftgiou(im1,ier)
      call ftinit(im1,'out.fits',iblock,ier)
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

