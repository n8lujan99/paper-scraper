
      parameter(nmax=10000)
      real x(nmax),y1(nmax),y2(nmax),xl(2),yl(2)
      
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.2)
      call pgscf(2)
      call pgslw(2)

      z1=0.
      z2=4.0
      nz=500
      do i=1,nz
         x(i)=z1+float(i-1)/float(nz-1)*(z2-z1)
         y1(i)=66.*sqrt(0.30*(1+x(i))**3+0.70)
         y2(i)=y1(i)/(1.+x(i))
         print *,x(i),y1(i)
      enddo

      xmin=0.
      xmax=3.8
      ymin=61.
      ymax=430.
      call pgvport(0.15,0.85,0.5,0.85)
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bcst',0.,0,'bcnst',0.,0)
      call pgsch(1.2)
      call pgmtxt('L',2.1,0.5,0.5,"H (km/s/Mpc)")

      xl(1)=1.88
      xl(2)=1.88
      yl(1)=ymin
      yl(2)=ymax
      call pgslw(3)
      call pgsls(4)
      call pgline(2,xl,yl)
      xl(1)=3.52
      xl(2)=3.52
      call pgline(2,xl,yl)
      call pgsls(1)
      call pgsci(2)
      call pgslw(6)
      call pgline(nz,x,y1)
      call pgslw(2)
      call pgsci(1)

      ymin=58.
      ymax=84.
      call pgvport(0.15,0.85,0.15,0.5)
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bcnst',0.,0,'bcnst',0.,0)
      call pgmtxt('L',2.1,0.5,0.5,"H/(1+z) (km/s/Mpc)")
      call pgmtxt('B',2.1,0.5,0.5,"z")
      xl(1)=1.88
      xl(2)=1.88
      yl(1)=ymin
      yl(2)=ymax
      call pgslw(3)
      call pgsls(4)
      call pgline(2,xl,yl)
      xl(1)=3.52
      xl(2)=3.52
      call pgline(2,xl,yl)
      call pgsls(1)
      call pgsci(2)
      call pgslw(6)
      call pgline(nz,x,y2)

      call pgend

      end
      
