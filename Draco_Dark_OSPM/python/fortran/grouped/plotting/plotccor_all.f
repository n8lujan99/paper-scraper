
      real f(1000),cc(1000),xin(10000),fn(1000)
      real fa(1000),ccall(1000,10000),fc0(1000),cc0(1000)
      character c1*19

      inorm=1
      
c      call pgbegin(0,'?',1,1)
      call pgbegin(0,'?',3,2)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(3)

      xmin=0.0
      xmax=40.0
      ymin=0.0
      ymax=1.0

      na=100
      do i=1,na
         fa(i)=xmin+(xmax-xmin)*float(i-1)/float(na-1)
      enddo

      open(unit=1,file='cc_norm.txt',status='old')
      ncc=0
      do i=1,1000
         read(1,*,end=777) x1,x2
         ncc=ncc+1
         fc0(ncc)=x1
         cc0(ncc)=x2
      enddo
 777  continue
      close(1)
      
      open(unit=1,file='l1',status='old')
      open(unit=11,file='out',status='unknown')

      ic=0
      nc=0
      nt=0
      do i=1,10000
         read(1,*,end=666) c1
         f50=-1.
         rms=-1.
         open(unit=2,file=c1,status='old',err=667)
         read(2,*,err=668) x1,x2
         f50=x2
         goto 669
 668     continue
         f50=100.
         rms=100.
 669     continue
         if(f50.gt.45.) then
            goto 667
         endif
         nt=nt+1
         xin(nt)=f50
         if(nc.eq.0) then
            call pgsci(1)
            call pgenv(xmin,xmax,ymin,ymax,0,0)
            call pglabel("Flux","Completeness","")
         endif
         nc=nc+1
         if(nc.eq.900) nc=0
         do j=1,100
            read(2,*) x1,x2
            f(j)=x1
            cc(j)=x2
         enddo
         do j=1,100
            fn(j)=f(j)/f50*11.
         enddo

         rms=0.
         do j=1,na
            call xlinint(fa(j),100,fn,cc,cct)
            ccall(j,nt)=cct
            if(fa(j).lt.18.) rms=rms+(cct-cc0(j))**2
         enddo

         ic=ic+1
         if(ic.eq.15) ic=1
         call pgsci(ic)
         if(inorm.eq.0) then
            call pgline(100,f,cc)
         else
            if(rms.lt.0.2) call pgline(100,fn,cc)
         endif
 667     close(2)
         write(11,*) c1," ",f50,rms
      enddo
 666  continue
      close(1)
      close(11)
      call pgend

      call biwgt(xin,nt,xb,xs)
      print *,xb,xs

      open(unit=11,file='out2',status='unknown')
      do i=1,na
         do j=1,nt
            xin(j)=ccall(i,j)
         enddo
         call biwgt(xin,nt,xb,xs)
         write(11,*) fa(i),xb,xs
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
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
