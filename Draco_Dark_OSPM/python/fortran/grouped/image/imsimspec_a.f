
      parameter (narrm1=1000,narrm2=440000)
      real xd(narrm1,narrm2),x(narrm1),y(narrm1)
      real xl(10000),yl(10000),ylf(10000)
      real*8 dx2
      integer naxes(2)
      character file1*180
      logical simple,extend,anyf

      jplot=24134
      iplot=0
      open(unit=1,file="lw_in",status="old")
      read(1,*)
      nl=0
      do i=1,10000
         read(1,*,end=669) x1,dx2
         if(x1.lt.1211..or.x1.gt.1219.) then
            nl=nl+1
            xl(nl)=x1
            yl(nl)=sngl(dx2/2.5d40)
            yl(nl)=(yl(nl)+5.)/5.5
         endif
      enddo
 669  continue
      close(1)

      file1='b0_chunk0.fits'

      ncol=1000
      nrow=100000
      open(unit=1,file="impact0.txt",status="old")
      do i=1,nrow
         read(1,*) (x(j),j=1,1000)
         do j=1,1000
            xd(j,i)=x(j)
         enddo
      enddo
      close(1)

      open(unit=1,file="wavelength.txt",status="old")
      n=0
      do i=1,1000
         read(1,*) x1
         n=n+1
         x(n)=x1
         ylf(n)=1.
         if(x1.gt.1201..and.x1.lt.1231.) ylf(n)=0.
      enddo
      close(1)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.4)
      call pgslw(2)

      rms0=1e10
      do j=1,nrow
         ymax=0.
         ymin=1e10
         rms=0.
         do i=1,ncol
            y(i)=exp(-xd(i,j))
            ymax=max(ymax,y(i))
            ymin=min(ymin,y(i))
            if(x(i).gt.1190.and.x(i).lt.1240.) rms=rms+(y(i)-ylf(i))**2
         enddo
         rms=sqrt(rms)
         if(rms.lt.rms0) then
            rms0=rms
            jmin=j
         endif
         if(j.eq.jplot) iplot=1
         if(iplot.eq.1) then
            call pgenv(1162.,1265.,0.,1.05,0,0)
            call pgline(n,x,y)
            call pgsci(2)
            call pgline(n,x,ylf)
            call pgsci(1)
            iplot=0
         endif
      enddo
      call pgend
      print *,jmin,rms0

 706  continue
      end
