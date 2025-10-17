
      parameter(nmax=10000)
      real w1(nmax,10),f1(nmax,10),t1(nmax,10)
      real w2(nmax),f2(nmax),t2(nmax)
      real e1(nmax,10),e2(nmax)
      character filea(10)*40,nfilea(10)*40
      data big/1.e10/

c      file1='bin02.ascii'
c      file2='bin03.ascii'
c      file3='region039.ascii'
c      nfile1='R=0.13"'
c      nfile2='R=0.24"'
c      nfile3='R=0.6"'
      filea(1)='f02fit.spec'
      filea(2)='f03fit.spec'
      filea(3)='f04fit.spec'
      filea(4)='f06fit.spec'
      nfilea(1)='R=0.075"'
      nfilea(2)='R=0.15"'
      nfilea(3)='R=0.25"'
      nfilea(4)='R=0.45"'

      nfile=4
      do j=1,nfile
         open(unit=1,file=filea(j),status='old')
         n1=0
         do i=1,nmax
            read(1,*,end=666) x1,x2,x3,x4,i5
            n1=n1+1
            w1(n1,j)=x1
            f1(n1,j)=x2
            t1(n1,j)=x3
            err=0.
            if(i5.eq.1) err=big
            e1(n1,j)=err
         enddo
 666     continue
         close(1)
      enddo

      call pgbegin(0,'?',1,1)
      call pgscf(2)
      call pgsch(1.0)
      call pgslw(2)

      xmin=2.25e4
      xmax=2.40e4
      ymin=0.65
      ymax=1.15

      do ip=1,4
         if(ip.eq.1) call pgvport(0.20,0.80,0.75,0.95)
         if(ip.eq.2) call pgvport(0.20,0.80,0.55,0.75)
         if(ip.eq.3) call pgvport(0.20,0.80,0.35,0.55)
         if(ip.eq.4) call pgvport(0.20,0.80,0.15,0.35)
         call pgwindow(xmin,xmax,ymin,ymax)
         if(ip.eq.1) call pgbox('bcst',0.,0,'bcnst',0.,0)
         if(ip.eq.2) call pgbox('bcst',0.,0,'bcnst',0.,0)
         if(ip.eq.3) call pgbox('bcst',0.,0,'bcnst',0.,0)
         if(ip.eq.4) call pgbox('bcnst',0.,0,'bcnst',0.,0)
         call pgsci(1)
         do i=1,n1
            w2(i)=w1(i,ip)
            f2(i)=f1(i,ip)
            e2(i)=e1(i,ip)
            t2(i)=t1(i,ip)
         enddo
         call plotline(n1,w2,f2,e2)
         call pgsci(2)
         call pgline(n1,w2,t2)
         call pgsci(1)
         call pgptxt(2.364e4,1.08,0.,0.,nfilea(ip))
      enddo

      call pgsch(1.5)
c      call pglabel('Wavelength (\(2078))','','')
      call pgmtxt('L',1.7,2.0,0.5,'Normalized Flux')
      call pgmtxt('B',1.9,0.5,0.5,'Wavelength (\(2078))')

      call pgend
      end

      subroutine plotline(n,x,y,ye)
      parameter(nmax=50000)
      real x(n),y(n),ye(n),xp(nmax),yp(nmax)
      data big/1.e10/
      ystart=ye(1)
      istart=1
      idone=0
      do j=1,nmax
         ibad=1
         if(ystart.lt.big) ibad=0
         np=0
         do i=istart,n
            if(ye(i).eq.ystart) then
               np=np+1
               xp(np)=x(i)
               yp(np)=y(i)
            else
               ystart=ye(i)
               istart=i
               goto 666
            endif
         enddo
         idone=1
 666     if(ibad.eq.1) then
            call pgsls(4)
            call pgline(np,xp,yp)
            call pgsls(1)
         else
            call pgline(np,xp,yp)
         endif
         if(idone.eq.1) goto 766
      enddo
 766    continue

      return
      end
