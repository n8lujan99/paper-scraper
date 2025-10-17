
      parameter(nmax=10000)
      real x(nmax),y1(nmax),y2(nmax)
      character cname*50
      
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)
      xmin=log10(0.1)
      xmax=log10(200.)
      ymin=log10(2.e-8)
      ymax=log10(1.e-3)
      call pgenv(xmin,xmax,ymin,ymax,0,30)
      call pglabel('Radius (arcsec)','\gn (L\D\(2281)\U/pc\U3\D)','')
      call pgslw(3)

      open(unit=1,file='list',status='old')

      ic=1
      do iall=1,100
         read(1,*,end=667) cname,rs
      

      open(unit=2,file=cname)
      n=0
      do i=1,nmax
         read(2,*,end=666) x1,x2
         n=n+1
         x(n)=log10(x1/rs)
         y1(n)=log10(x2)
      enddo
 666  continue
      close(2)

      ic=ic+1
      if(iall.eq.9) then
         ic=1
         call pgslw(7)
      endif
      call pgsci(ic)
      call pgline(n,x,y1)

      enddo
 667  continue
      close(1)
      call pgend

      end
      
