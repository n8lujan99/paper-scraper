
      real xn(7),xf(7,100),xc(7,100),xin(100),yin(100)
      real xini(100),yini(100),xci(7,100),yinis(100)
      character ctitle*40

      nw=7
      nw=3
      open(unit=1,file='in',status='old')
      read(1,*)
      if(nw.eq.7) read(1,*) x0,xn(1),xn(2),xn(3),xn(4),xn(5),xn(6),xn(7)
      if(nw.eq.3) read(1,*) x0,xn(1),xn(2),xn(3)
      xs=xn(2)
      if(xn(2).gt.90.) xs=xn(3)
      n=0
      do i=1,100
         if(nw.eq.7) read(1,*) x0,x1,x2,x3,x4,x5,x6,x7
         if(nw.eq.3) read(1,*) x0,x1,x2,x3
         do j=1,nw
           xf(j,i)=x0/xn(j)*xs
c           xf(j,i)=x0
         enddo
         xc(1,i)=x1
         xc(2,i)=x2
         xc(3,i)=x3
         xc(4,i)=x4
         xc(5,i)=x5
         xc(6,i)=x6
         xc(7,i)=x7
      enddo
 666  continue
      close(1)

      ctitle="          "
      open(unit=1,file='title',status='old',err=456)
      read(1,*) ctitle
 456  close(1)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=0.0
      xmax=15.0
      xmax=45.0
c      xmax=100.0
      ymin=0.0
c      ymax=1.0
      ymax=1.05
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Flux (1e-17 ergs/cm\U2\D/s)','Completeness',ctitle)

      call pgslw(8)
      do iw=1,nw
         do j=1,100
            xin(j)=xf(iw,j)
            yin(j)=xc(iw,j)
         enddo
         ic=9-iw
         if(ic.ge.6) ic=ic+2
         call pgsci(ic)
         call pgline(100,xin,yin)
      enddo

      do j=1,100
         xini(j)=xf(5,j)
      enddo
      do j=1,100
         xini(j)=1.0+float(j-1)*0.5
c         xini(j)=1.0+float(j-1)*1.0
      enddo
      do iw=1,nw
         do j=1,100
            xin(j)=xf(iw,j)
            yin(j)=xc(iw,j)
         enddo
         do j=1,100
            call xlinint(xini(j),100,xin,yin,y0)
            xci(iw,j)=y0
         enddo
      enddo

      do j=1,100
         do iw=1,nw
            xin(iw)=xci(iw,j)
         enddo
         call biwgt(xin,nw,xb,xs0)
         yini(j)=xb
         yinis(j)=xs0
      enddo      

      call pgsci(1)
      call pgsls(4)
      call pgline(100,xini,yini)
      do j=1,99
         if(yini(j).gt.0.95) yini(j+1)=min(1.,yini(j)+0.01)
      enddo
      do j=100,2,-1
         if(yini(j).lt.yini(j-1)) yini(j-1)=yini(j)
      enddo
      x1=xs*0.9
      y1=0.2
      do j=1,100
         if(xini(j).lt.x1) then
            xlow=y1+(1.-y1)*xini(j)/x1
            yini(j)=yini(j)*xlow
         endif
      enddo
      call pgslw(10)
      call pgsls(1)
      call pgline(100,xini,yini)

      open(unit=11,file='out4',status='unknown')
      write(11,1102) 0.5,xs
      do j=1,100
         write(11,1101) xini(j),yini(j),yinis(j)
      enddo
      close(11)

      call pgend

 1101 format(3(f10.7,1x))
 1102 format(2(f10.7,1x))
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
