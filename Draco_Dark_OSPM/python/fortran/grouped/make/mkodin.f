
      parameter(nmax=10000)
      real arr(1000,1000),f1(nmax)
      integer naxes(2)
      parameter(pi=3.14159)

      fn=0.25
      fp=0.2
      fp=0.

      see=1.4
      rms=1.8 
      rms=1.0
      flux0=8.0
      ng=700.
      
      size=10.2
      dx=0.27
      
      nx=nint(size/dx)
      ny=nx
      icen=nint(float(nx)/2.)
      xcen0=float(icen)
      ycen0=xcen0

      sigs=see/2.35

      open(unit=1,file='in.dat',status='old')
      read(1,*)
      read(1,*)
      nf=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         nf=nf+1
         f1(nf)=x2
      enddo
 666  continue
      close(1)
      print *,nf
      nfp=nint(float(nf)*fp)
      print *,nfp
      ng=nf      

      idum=-1
      
      do i=1,1000
         do j=1,1000
            arr(i,j)=0.
         enddo
      enddo

      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      iblock=1
      igc=0
      ier=0

      gn=1./sqrt(2.*pi*sigs*sigs)*dx*dx
      sum=0.
      do ig=1,ng
         do i=1,nx
            do j=1,ny
               arr(i,j)=arr(i,j)+fn*gasdev(idum)
            enddo
         enddo

         if(ig.gt.nfp) then
c            flux=flux0
            flux=f1(ig)
            xcen=xcen0+rms*gasdev(idum)/dx
            ycen=ycen0+rms*gasdev(idum)/dx
            do i=1,nx
               x=float(i)
               do j=1,ny
                  y=float(j)
                  r=dx*sqrt((xcen-x)**2+(ycen-y)**2)
                  w=(r)/sigs
                  gaus=flux*gn*exp(-w*w/2.)/sqrt(2.*pi*sigs**2)
                  sum=sum+gaus
                  arr(i,j)=arr(i,j)+gaus
               enddo
            enddo
         endif
      enddo

      sum=0.
      do i=1,nx
         do j=1,ny
            arr(i,j)=arr(i,j)/float(ng)
            sum=sum+arr(i,j)
         enddo
      enddo
      print *,sum,flux0
           
      call ftinit(50,'image.fits',iblock,ier)
      call ftphps(50,-32,naxis,naxes,ier)
      if(ier.ne.0) then
         print *,'Error in output file ',ier
      endif
      call ftp2de(50,igc,1000,naxes(1),naxes(2),arr,ier)
      call ftclos(50,ier)

      end
