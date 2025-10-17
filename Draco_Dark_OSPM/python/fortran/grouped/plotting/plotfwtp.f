
      real y(10000),f(10000)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=1,file='fwhm.all',status='old')
      n=0
      do i=1,10000
         read(1,*,end=666) x1,x2
         n=n+1
         y(n)=x1+2017.
         f(n)=x2
      enddo
 666  continue
      close(1)

      xmin=0.+2017.
      xmax=5.05+2017.
      xmax=2023.8
      ymin=1.2
      ymax=2.7
      call pgsch(1.0)
      call pgvport(0.15,0.85,0.73,0.95)
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bcmst',0.,0,'bcnst',0.,0)
      call pgpt(n,y,f,17)
      call pgmtxt('L',2.1,0.5,0.5,"FWHM")

      open(unit=1,file='tp.all',status='old')
      n=0
      do i=1,10000
         read(1,*,end=667) x1,x2,x3,x4
         n=n+1
         y(n)=x1+2017.
         f(n)=x4
      enddo
 667  continue
      close(1)

      xmin=0.+2017.
      xmax=5.05+2017.
      xmax=2023.8
      ymin=0.08
      ymax=0.18
      call pgsch(1.0)
      call pgvport(0.15,0.85,0.51,0.73)
      call pgwindow(xmin,xmax,ymin,ymax)
c      call pgbox('bcst',0.,0,'bcnst',0.,0)
      call pgbox('bncst',0.,0,'bcnst',0.,0)
      call pgpt(n,y,f,17)
      call pgmtxt('L',2.1,0.5,0.5,"tp at 4940")
      call pgmtxt('B',2.3,0.5,0.5,"Year")

      goto 888

      open(unit=1,file='tp.all',status='old')
      n=0
      do i=1,10000
         read(1,*,end=668) x1,x2,x3,x4
         n=n+1
         y(n)=x1+2017.
         f(n)=x2
      enddo
 668  continue
      close(1)

      xmin=0.+2017.
      xmax=5.05+2017.
      xmax=6.40+2017.
      ymin=3.6
      ymax=4.7
      call pgsch(1.0)
      call pgvport(0.15,0.85,0.29,0.51)
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bcst',0.,0,'bcnst',0.,0)
      call pgpt(n,y,f,17)
      call pgmtxt('L',2.1,0.5,0.5,"tp 4940/3540")

      open(unit=1,file='nlae.dat',status='old')
      n=0
      do i=1,10000
         read(1,*,end=669) x1,x2
         n=n+1
         y(n)=float(i)/12.+2017.
         f(n)=x2
      enddo
 669  continue
      close(1)

      xmin=0.+2017.
      xmax=5.05+2017.
      xmax=6.40+2017.
      ymin=1.0
      ymax=3.4
      call pgsch(1.0)
      call pgvport(0.15,0.85,0.07,0.29)
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bcnst',0.,0,'bcnst',0.,0)
      call pgpt(n,y,f,17)
      call pgmtxt('L',2.1,0.5,0.5,"N\DLAE\U/N\DIFU")
      call pgmtxt('B',2.1,0.5,0.5,"Year")
 888  continue

      call pgend

      end



