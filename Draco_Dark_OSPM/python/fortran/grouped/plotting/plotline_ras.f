
      parameter(nmax=10000)
      real x(nmax),y(nmax),y1(nmax),xa(100,nmax)
      character cname*50

      nsm=11
      nsmh=5
c      nsm=23
c      nsmh=11
      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)
      xmin=1195.
      xmax=1235.
      xmin=1180.
      xmax=1250.
      ymin=-0.2
      ymax=2.1
      call pgenv(xmin,xmax,ymin,ymax,0,00)
      call pglabel('Wavelength','Flux','')
      call pgslw(1)

      open(unit=1,file='list',status='old')

      ic=1
      ia=0
      do iall=1,100
         read(1,*,end=667) cname
      
         open(unit=2,file=cname,status='old')
         n=0
         nb=0
         do i=1,nmax
            read(2,*,end=666) x1,x2
            n=n+1
            x(n)=x1
            y(n)=x2
c            if(x1.le.1205..or.x1.ge.1222.) then
            if(x1.le.1190..or.x1.ge.1240.) then
               nb=nb+1
               y1(nb)=y(n)
            endif
         enddo
 666     continue
         close(2)
         call biwgt(y1,nb,xb,xs)

         do i=1,n
c            y(i)=y(i)-y1(i)
c            y(i)=y1(i)
c     y(i)=y1(i)-y(i)
            y(i)=y(i)/xb
         enddo
         do i=1,n
            sum=0.
            jstart=max(i-nsmh,1)
            jend=min(jstart+nsm,n)
            ns=0
            do j=i,i+nsm
               sum=sum+y(j)
               ns=ns+1
            enddo
            y1(i)=sum/float(ns)
         enddo
         ic=ic+1
         if(ic.eq.14) ic=2
         call pgsci(ic)
c         call pgline(n,x,y1)
         ia=ia+1
         do i=1,n
            xa(ia,i)=y1(i)
         enddo

      enddo
 667  continue
      close(1)

      do i=1,n
         do j=1,ia
            y(j)=xa(j,i)
         enddo
         call biwgt(y,ia,xb,xs)
         y1(i)=xb
      enddo
      call pgsci(1)
      call pgslw(6)
      call pgline(n,x,y1)

      cname='prof.lw'
      open(unit=2,file=cname,status='old')
      n=0
      nb=0
      do i=1,nmax
         read(2,*,end=668) x1,x2
         n=n+1
         x(n)=x1
         y(n)=x2
         if(x1.le.1205..or.x1.ge.1222.) then
            nb=nb+1
            y1(nb)=y(n)
         endif
      enddo
 668  continue
      close(2)
      call biwgt(y1,nb,xb,xs)
      do i=1,n
         y(i)=y(i)/xb
      enddo
      call pgsci(2)
      call pgline(n,x,y)
      
      call pgend

      end
      
