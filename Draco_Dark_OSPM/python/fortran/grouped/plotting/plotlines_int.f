
      parameter(nmax=10000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(100)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      character file1*80,file2*80,c1*3

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=0.
      xmax=11.
      ymin=0.1
      ymax=3.

      open(unit=1,file='linelist',status='old')

      do ia=1,700
         read(1,*,end=666) file1,ic
         open(unit=2,file=file1,status='old')
         if(ia.eq.1) then
            call pgsci(1)
            call pgenv(xmin,xmax,ymin,ymax,0,0)
            call pglabel('Radius(arcmin)',
     $              '\gs\Dr\U/\gs\Dt','')
         endif
         n=0
         do i=1,8000
            read(2,*,end=667) x1,x2
            n=n+1
            x(n)=x1/60.
            y(n)=x2
         enddo
 667     continue
         close(2)
         call pgslw(8)
         call pgsci(ic)
         call pgline(n,x,y)
         call pgslw(1)
      enddo
 666  continue
      close(1)

      call pgend

      end
