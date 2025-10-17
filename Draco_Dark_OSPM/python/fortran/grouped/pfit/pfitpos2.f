
      parameter(nmax=100000)

      real v(nmax),y(nmax),yl(nmax),yh(nmax),a(6),yherm(nmax),covar(6,6)
      real xp(nmax),yp(nmax),v0(nmax),y0(nmax)
      integer ia(6)
      character file1*50,file2*5,file3*6
      data big/1.e30/

      parameter(pi=3.1415926539)

      vref0=-75.

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

      rms=0.
      nrms=0
      do if=1,nmax
         read(1,*,end=666) file1,xp1,xp2
         do i=1,40
            if(file1(i:i).eq.' ') then
               nfile=i-1
               goto 966
            endif
         enddo
 966     continue
         open(unit=2,file=file1(1:nfile),status='old')
         file2="mod"//file1(2:3)
c         file3="mod1"//file1(2:3)
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
            v(nl)=x2+vref+vref0
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
         sum0=0.
         do i=1,nl
            sum0=sum0+y(i)
            v0(i)=v(i)
            y0(i)=y(i)
         enddo
         sum0=sum0*(v(2)-v(1))
         open(unit=2,file=file2,status='old',err=556)
c         open(unit=2,file=file3,status='old',err=556)
         nm=0
         sum1=0.
         do i=1,nmax
            read(2,*,end=555) x1,x2
            nm=nm+1
            v(nm)=x1
            y(nm)=x2
            sum1=sum1+y(nm)
         enddo
 555     continue
         sum1=sum1*(v(2)-v(1))
         do i=1,nm
            y(i)=y(i)/sum1*sum0
            call xlinint(v(i),nl,v0,yh,yh0)
            call xlinint(v(i),nl,v0,yl,yl0)
            xp0=(yh0+yl0)/2.
            err=(yh0-yl0)/2.
            rms=rms+((y(i)-xp0)/err)**2
            nrms=nrms+1
         enddo
         call pgsci(4)
         call pgslw(6)
         call pgline(nm,v,y)
 556     continue
         close(2)
      enddo
 666  continue
      close(1)
      rms=sqrt(rms/float(nrms))
      open(unit=11,file='rms.out',status='old')
      write(11,*) rms
      close(11)
      print *,rms

      call pgsci(1)
      call pgvport(0.05,0.95,0.05,0.95)
      call pgwindow(-0.25,0.25,-0.25,0.25)
      call pgbox('bcstn',0.,0,'bcstn',0.,0)

      call pgend

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
