Marked
      parameter(nmax=2000)
      real w(nmax),xd1(nmax,nmax),xd2(nmax,nmax)
      real w0(nmax),f0(nmax),s(nmax,nmax),xin(nmax),s0(nmax,nmax)
      integer naxes(2)
      character file1*60
      logical simple,extend,anyf

      open(unit=1,file='list',status='old')

      do i=1,1036
         w(i)=3470.+float(i-1)*2.
      enddo

      ier=0
      im1=50
      nall=0
      do iall=1,10000
         read(1,*,end=666) file1
         nall=nall+1
         call ftopen(im1,file1,0,iblock,ier)
         call ftmahd(im1,12,ihd,ier)
         call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=1032
         nrow=112
         call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd1,anyf,ier)
         call ftmahd(im1,16,ihd,ier)
         call ftg2de(im1,igc,0.,nmax,ncol,nrow,xd2,anyf,ier)
         call ftclos(im1,ier)
         do j=1,112
            do i=1,1032
               w0(i)=xd1(i,j)
               f0(i)=xd2(i,j)
            enddo
            do i=1,1036
               call xlinint(w(i),1032,w0,f0,f1)
               s(i,j)=f1
            enddo
         enddo

         do i=1,1036
            nin=0
            do j=1,112
               if(s(i,j).ne.0) then
                  nin=nin+1
                  xin(nin)=s(i,j)
               endif
            enddo
            call biwgt(xin,nin,xb,xs)
            n2=nint(0.941*float(nin))
            call biwgt(xin,n2,xb,xs)
            s0(i,nall)=xb
c            print *,i,w(i),s0(i,nall)
         enddo
         print *,nall,ier

      enddo
 666  continue
      close(1)

      open(unit=11,file='sky.out',status='unknown')
      do i=1,1036
         nin=0
         do j=1,nall
            if(s0(i,j).ne.0) then
               nin=nin+1
               xin(nin)=s0(i,j)
            endif
         enddo
         call biwgt(xin,nin,xb,xs)
         write(11,*) w(i),xb,xs
      enddo
      close(11)

      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.le.x(1)) yp=0.
      if(xp.ge.x(n)) yp=0.
      return
      end
      subroutine xlinint2(xp,n,x,y,yp,jin,jout)
      real x(n),y(n)
      do j=jin,n-1
         if(xp.ge.x(j).and.xp.le.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            jout=j
            return
         endif
      enddo
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      jout=1
      return
      end
