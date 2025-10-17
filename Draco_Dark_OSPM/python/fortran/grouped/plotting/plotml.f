
      parameter(nmax=1000)
      real suma(nmax),suma2(nmax),suml2(nmax),sumh2(nmax)
      real sumall(nmax,nmax),sumalt(nmax,nmax),suml(nmax),sumh(nmax)
      real xlight(nmax),xlights(nmax),xmll(nmax),xmlt(nmax),xrad(nmax)
      integer iang(10)

      character file1*50,gname*50,file2*80

      na=5
      na=2
      open(unit=1,file='galint.dat',status='old')
      read(1,*,end=667) gname,file1,rmin,rmax,vmin,vmax,rhalf,icore
c         read(1,*) file2
      close(1)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgslw(2)
      call pgsch(1.4)

      open(unit=1,file='bindemo_r.out',status='old')
      do i=1,5
         read(1,*)
      enddo
      do i=1,80
         read(1,*,end=664) i1,i2,x3,x4,x5,x6
         xrad(i)=log10(x6)
      enddo
 664  continue
      close(1)

      open(unit=1,file='dL_mod.dat',status='old')
      sum2=0.
      do i=1,80
         sum=0.
         do j=1,20
            read(1,*,end=665) i1,i2,x3
            sum=sum+x3
         enddo
         sum2=sum2+sum
         xlight(i)=sum
         xlights(i)=sum2
      enddo
 665  continue
      close(1)

      open(unit=1,file='mlist',status='old')
      nt=0
      do iall=1,1000
         read(1,*,end=667) file2
         read(file2(10:14),*) xml
         read(file2(16:24),*) xbh
         nt=nt+1
         open(unit=2,file=file2,status='old')
         sum2=0.
         sum2=xbh
         do i=1,80
            sum=0
            do j=1,20
               read(2,*,end=666) i1,i2,x3
               sum=sum+x3
            enddo
            sum2=sum2+sum
            xmll(i)=sum/xlight(i)/xml
            xmlt(i)=sum2/xlights(i)/xml
c            print *,xrad(i),xmll(i),xmlt(i)
            sumall(i,nt)=xmll(i)
            sumalt(i,nt)=xmlt(i)
         enddo
 666     continue
         close(2)
      enddo
 667  continue

      open(11,file='ml.out',status='unknown')
      do i=1,80
         avg=0.
         rmax=-1e10
         rmin=1e10
         avg2=0.
         rmax2=-1e10
         rmin2=1e10
         do j=1,nt
            avg=avg+sumall(i,j)
            rmax=max(rmax,sumall(i,j))
            rmin=min(rmin,sumall(i,j))
            avg2=avg2+sumalt(i,j)
            rmax2=max(rmax2,sumalt(i,j))
            rmin2=min(rmin2,sumalt(i,j))
         enddo
         avg=avg/float(nt)
         suma(i)=avg
         suml(i)=rmin
         sumh(i)=rmax
         avg2=avg2/float(nt)
         suma2(i)=avg2
         suml2(i)=rmin2
         sumh2(i)=rmax2
         write(11,1001) 10**xrad(i),suma(i),suml(i),sumh(i),suma2(i),
     $        suml2(i),sumh2(i)
      enddo
      close(11)

      xmin=log10(0.05)
      xmax=log10(20.)
      ymin=0.8
      ymax=6.
      call pgslw(2)
      call pgenv(xmin,xmax,ymin,ymax,0,10)
      call pglabel('Radius (")','M/L (local and enclosed)','')
      call pgsci(1)
      call pgslw(6)
      call pgsls(1)
      call pgline(80,xrad,suma)
      call pgsls(4)
      call pgline(80,xrad,suml)
      call pgline(80,xrad,sumh)

      call pgsci(4)
      call pgsls(1)
      call pgline(80,xrad,suma2)
      call pgsls(4)
      call pgline(80,xrad,suml2)
      call pgline(80,xrad,sumh2)

      call pgend

 1001 format(7(1x,f8.3))
      end

      subroutine getlim(file2,iang,rm)
      parameter(nb=10)
      integer iang(nb)
      character file2*50

      do i=1,nb
         iang(i)=0
      enddo
      rm=0.
      open(unit=12,file=file2,status='old')
      do i=1,1000
         read(12,*,end=666) i1,x2
         iang(i1)=1
         rm=max(rm,x2)
      enddo
 666  continue
      close(12)
      return
      end
