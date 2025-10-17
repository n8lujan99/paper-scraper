
      parameter(nmax=70)
      real xarr(1344),yarr(1344)
      real arr(nmax,nmax),xin(1344*1036)
      real*8 da,db
      integer naxes(3),iarr(1344)
      character file1*40,comm*132
      parameter(pi=3.14159,radtodeg=57.29578)

      open(unit=1,file='radec.in',status='old')
      read(1,*) ra0,dec0
      close(1)
      open(unit=1,file='radec2.dat',status='old')
      read(1,*) x1,x2,pa0
      open(unit=1,file='fit.in',status='old')
      read(1,*) da,db
      close(1)

      cdec=cos(dec0/radtodeg)
      cosp=cos(pa0/radtodeg)
      sinp=sin(pa0/radtodeg)

      do i=1,nmax
         do j=1,nmax
            arr(i,j)=1.
         enddo
      enddo

      dx=2.0
      nx=31
      ny=31
      nw=1036
      xcen=2.*float(nx)/dx
      xcenwcs=float(nx)/2.

      open(unit=1,file='in',status='old')

      nin=0
      nfib=0
      do j=1,1344*1036
         read(1,*,end=666) x1,x2,x3,x4,x5
         xd=3600.*cdec*(x1-ra0)
         yd=3600.*(x2-dec0)
         xdn=cosp*xd-sinp*yd
         ydn=sinp*xd+cosp*yd
         xarr(j)=xdn
         yarr(j)=ydn
         call getflag(x1,x2,da,db,iflag)
         iarr(j)=iflag
         nfib=nfib+1
      enddo
 666  continue
      close(1)
      print *,"read in",nfib

      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      iblock=1
      igc=0
      ier=0

      xoff=xcen
      yoff=xcen
      do i=1,nx
         xp=float(i)*dx-xoff
         do j=1,ny
            yp=float(j)*dx-yoff
            diff=3.0
            radmin=1e10
            do k=1,nfib
               x=xarr(k)
               y=yarr(k)
               rad=sqrt((xp-x)**2+(yp-y)**2)
               if(rad.lt.radmin) then
                  radmin=rad
                  kmin=k
               endif
            enddo
            arr(i,j)=float(iarr(kmin))
         enddo
      enddo

      pa0=-pa0
      radfac=1./(3600./dx)
      cd11=radfac*cos(pa0/radtodeg)
      cd12=-radfac*sin(pa0/radtodeg)
      cd21=radfac*sin(pa0/radtodeg)
      cd22=radfac*cos(pa0/radtodeg)
       
      call ftinit(50,'image3d.fits',iblock,ier)
      call ftphps(50,-32,naxis,naxes,ier)
      call ftpkye(50,'CRVAL1',ra0,5,"RA",ier)
      call ftpkye(50,'CRVAL2',dec0,5,"DEC",ier)
      call ftpkye(50,'CROTA2',pa0,5,"PA",ier)
      call ftpkye(50,'CD1_1',cd11,5,"",ier)
      call ftpkye(50,'CD1_2',cd12,5,"",ier)
      call ftpkye(50,'CD2_2',cd22,5,"",ier)
      call ftpkye(50,'CD2_1',cd21,5,"",ier)
      call ftpkye(50,'CRPIX1',xcenwcs,5,"X0",ier)
      call ftpkye(50,'CRPIX2',xcenwcs,5,"Y0",ier)
      call ftpkys(50,'CUNIT1',"deg","",ier)
      call ftpkys(50,'CUNIT2',"deg","",ier)
      call ftpkys(50,'CTYPE1',"RA---TAN","",ier)
      call ftpkys(50,'CTYPE2',"DEC--TAN","",ier)
      
      print *,naxes(1),naxes(2)
      call ftp2de(50,igc,nmax,naxes(1),naxes(2),arr,ier)
      call ftclos(50,ier)

c      call ftinit(50,'image3df.fits',iblock,ier)
c      call ftphps(50,-32,naxis,naxes,ier)
c      if(ier.ne.0) then
c         print *,'Error in output file ',ier
c      endif
c      print *,naxes(1),naxes(2)
c      call ftp3de(50,igc,100,100,naxes(1),naxes(2),ntot,arr2,ier)
c      call ftclos(50,ier)

      end

      subroutine getflag(x1,x2,da,db,iflag)
      real*8 da,db,d1,d2

      dmin0=8.
      rtot=0.2
      nstep=10000
      rmin=x1-rtot/2.
      rmax=rmin+rtot
      dmin=1e10
      cosd=cos(x2/57.3)
      do j=1,nstep
         r=rmin+float(j-1)/float(nstep-1)*(rmax-rmin)
         d1=da+db*dble(r)
         d2=3600.d0*dabs(d1-dble(x2))
         ds=sngl(d2)
         rs=3600.*abs(r-x1)*cosd
         dist=sqrt(ds*ds+rs*rs)
         dmin=min(dmin,dist)
      enddo
      iflag=1
      if(dmin.lt.dmin0) iflag=0

      return
      end
