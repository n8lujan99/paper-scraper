
      parameter(nmax=100000)

      real v(nmax),y(nmax),yl(nmax),yh(nmax),a(6),yherm(nmax),covar(6,6)
      real xp(nmax),yp(nmax)
      integer ia(6)
      character file1*50
      data big/1.e30/

      parameter(pi=3.1415926539)

      imc=1
      if(imc.eq.0) file1="blist"
      if(imc.eq.1) file1="blist2"
      yminc=0.003

      open(unit=1,file=file1,status='old')

 1001 format(7x,i3)

c      call pgbegin(0,'?',3,3)
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.2)

      vmin=-1500.
      vmax=1500.
      ymin=-0.005
      ymax=0.125

      do if=1,nmax
         read(1,*,end=666) file1,xp1,xp2
         do i=1,40
            if(file1(i:i).eq.' ') then
               nfile=i-1
               goto 966
            endif
         enddo
 966     continue
         vref=0
         open(unit=2,file=file1(1:nfile),status='old')
         vref=0.
         if(imc.eq.0) read(2,*) vref,x2,x3,ntot
         nl=0
         do i=1,nmax
            if(imc.eq.0) then
               read(2,*,end=667) i1,x2,x3
            else
               read(2,*,end=667) x2,x3,x4,x5
            endif
            nl=nl+1
            v(nl)=x2+vref
            y(nl)=x3
            if(imc.eq.1) then
               yd1=x5-x3
               yd1=max(yd1,yminc)
               x5=x3+yd1

               yd2=x3-x4
               yd2=max(yd2,yminc)
               x4=x3-yd2

               if(yd1.gt.yd2) x4=x3-yd1
               yl(nl)=max(0.,x4)
               yh(nl)=x5
            endif
         enddo
 667     continue
         close(2)

         xpix=11.1
         xhalf=1./xpix-0.008
c         xhalf=1./xpix
         x1=0.5+xp1/xpix-xhalf
         x2=x1+2.*xhalf
         y1=0.5+xp2/xpix-xhalf
         y2=y1+2.*xhalf

         call pgsci(14)
         call pgslw(1)
         call pgvport(x1,x2,y1,y2)
         call pgwindow(vmin,vmax,ymin,ymax)
         call pgbox('bc',0.,0,'bc',0.,0)
         call pgsci(2)
         call pgslw(1)
         if(imc.eq.1) then
            np=0
            do i=1,nl
               np=np+1
               xp(np)=v(i)
               yp(np)=yh(i)
            enddo
            do i=nl,1,-1
               np=np+1
               xp(np)=v(i)
               yp(np)=yl(i)
            enddo
            call pgpoly(np,xp,yp)
         endif
         call pgsci(4)
         call pgslw(4)
         call pgline(nl,v,y)
      enddo
 666  continue
      close(1)

      call pgsci(1)
      call pgvport(0.05,0.95,0.05,0.95)
      call pgwindow(-0.25,0.25,-0.25,0.25)
      call pgbox('bcstn',0.,0,'bcstn',0.,0)

      call pgend

      end
