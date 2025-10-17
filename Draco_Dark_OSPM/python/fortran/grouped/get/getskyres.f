Skipped
      parameter(nmax=50000)
      real xd(1036,34944),xin(34944),yall(1036,2000)
      real ws(nmax),sky(nmax),xp(1036),yp(1036)
      integer naxes(2)
      character file1*80,date*8,shot*3,exp*5
      logical simple,extend,anyf

      xntry=1.0

      irat=1

      call pgbegin(0,'?',1,1)
      call pgask(.false.)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

      xmin=3500.
      xmax=5500.
      if(irat.eq.0) then
         ymin=-1.
         ymax=4.
      else
         ymin=-0.02
         ymax=0.07
      endif

c      call pgenv(xmin,xmax,ymin,ymax,0,0)
c      call pglabel("Wavelength","Counts","")

      open(unit=2,file='listin',status='old')
      nall=0
      do iall=1,10000
      read(2,*,end=888) date,shot,exp
      nall=nall+1

      file1="/scratch/00115/gebhardt/alldet/output/d"//date//"s"//
     $     shot//exp//"sky.dat"
      open(unit=1,file=file1,status='unknown')
      ns=0
      do i=1,nmax
         read(1,*,end=666) x1,x2
         ns=ns+1
         ws(ns)=x1
         sky(ns)=x2
      enddo
 666  continue
      close(1)

      file1="/scratch/00115/gebhardt/alldet/output/d"//date//"s"//
     $     shot//exp//"sub.fits"
      iext=1

      im1=0
      ier=0
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      naxis=2
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1036,ncol,nrow,xd,anyf,ier)
      call ftclos(im1,ier)

      do i=1,ncol
         nin=0
         do j=1,nrow
            if(xd(i,j).ne.0) then
               nin=nin+1
               xin(nin)=xd(i,j)
            endif
         enddo
         call biwgt(xin,nin,xb,xs)
         ntry=nint(xntry*float(nin))
         call biwgt(xin,ntry,xb,xs)
         iwave=3470+2*(i-1)
         call xlinint(float(iwave),ns,ws,sky,sky0)
         if(nin.lt.100) xb=0.
         if(nin.lt.100) xs=0.
c         write(11,*) iwave,xb,xs,nin,sky0

         xp(i)=float(iwave)
         yp(i)=xb/sky0
         yall(i,nall)=yp(i)
      enddo

      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel("Wavelength","Counts","")

      call pgline(ncol,xp,yp)

      open(11,file='out',status='unknown')
      do i=1,ncol
         do j=1,nall
            xin(j)=yall(i,j)
         enddo
         call biwgt(xin,nall,xb,xs)
         yp(i)=xb
         write(11,*) xp(i),yp(i),xs,nall
      enddo
      close(11)
      call pgsci(2)
      call pgslw(4)
      call pgline(ncol,xp,yp)
      call pgsci(1)
      call pgslw(2)

      enddo
 888  close(2)

      call pgend

 706  continue
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
