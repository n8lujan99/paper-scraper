Marked
      parameter(nmax=10000)
      real f50(nmax),zl(nmax),xl(nmax),z(nmax)
      real zv(nmax),vol(nmax),xlum(nmax),sn(nmax)
      real xall(10,13)
      real xwc(7),xhc(7),xcc(199,7),xic(199)
      real*8 dx2
      integer*8 i1

      open(unit=1,file='f50b.dat',status='old')
      ns=0
      do i=1,nmax
         read(1,*,end=666) i1,x2
         ns=ns+1
         f50(ns)=x2
      enddo
 666  continue
      close(1)

      open(unit=1,file='zvol.dat',status='old')
      nv=0
      do i=1,nmax
         read(1,*,end=667) x1,x2
         nv=nv+1
         zv(nv)=x1
         vol(nv)=x2
      enddo
 667  continue
      close(1)

c- get the flux to luminosity as a function of z
      open(unit=1,file='FtoLz2.dat',status='old')
      nz=0
      do i=1,10000
         read(1,*,end=668) x1,dx2
         nz=nz+1
         zl(nz)=x1
         xl(nz)=sngl(dx2/1.d40)
c         xl(nz)=(1.+x1)*sngl(dx2/1.d40)
      enddo
 668    continue
      close(1)

      zmin=1.88
      zmax=3.55
      nz2=100
      do i=1,nz2
         z(i)=zmin+(zmax-zmin)*float(i-1)/float(nz2-1)
      enddo

      nsn=6
      sn(1)=4.8
      sn(2)=5.0
      sn(3)=5.5
      sn(4)=6.0
      sn(5)=6.5
      sn(6)=7.0

      nlum=13
      xlum(1)=42.3
      xlum(2)=42.4
      xlum(3)=42.5
      xlum(4)=42.6
      xlum(5)=42.7
      xlum(6)=42.8
      xlum(7)=42.9
      xlum(8)=43.0
      xlum(9)=43.1
      xlum(10)=43.2
      xlum(11)=43.3
      xlum(12)=43.4
      xlum(13)=43.5

      do isn=1,nsn
      sn0=sn(isn)
      call readsncc(sn0,xwc,xhc,xic,xcc)
      do ilum=1,nlum
      x0=xlum(ilum)
      nt0a=0
      nt1a=0
      v0=0.
      v1=0.

      xtot=0.
      do iz=1,nz2
         z0=z(iz)
         call xlinint(z0,nz,zl,xl,xconv)
         call xlinint(z0,nv,zv,vol,vol0)
         f0=10**(x0-40.)/xconv*1.e17
         wave=1215.67*(1.+z0)
         xmin=1e10
         do iw=1,7
            if(abs(wave-xwc(iw)).lt.xmin) then
               xmin=abs(wave-xwc(iw))
               imin=iw
            endif
         enddo
c         print *,wave,xwc(imin),xhc(imin),f0
         do i=1,ns
            f50u=f50(i)*sn0
            fmult=f50u/xhc(imin)
            call getint(f0,fmult,imin,xic,xcc,xcc0)
c            if(i.eq.1) print *,i,f50u,xhc(imin),fmult,xcc0,f0
c            read *
            xtot=xtot+xcc0*vol0
         enddo
      enddo
      xall(isn,ilum)=xtot/1.e7
      enddo
      enddo

      do ilum=1,nlum
         print *,xlum(ilum),(xall(i,ilum),i=1,nsn)
      enddo

      end

      subroutine getint(f0,fmult,imin,xic,xcc,xcc0)
      real xcc(199,7),xic(199)
      f0t=f0/fmult
      if(f0t.le.xic(1)) then
         xcc0=0.
         return
      endif
      xcc0=1.
      do i=1,198
         xlo=xic(i)
         xhi=xic(i+1)
         if(f0t.ge.xlo.and.f0t.lt.xhi) then
            xcc0=xcc(i,imin)
            return
         endif
      enddo      

      return
      end

      subroutine readsncc(sn0,xwc,xhc,xic,xcc)
      real xwc(7),xhc(7),xcc(199,7),xic(199)
      character file1*9
      file1="sn4.8.use"
      write(file1(3:5),1001) sn0
 1001 format(f3.1)
      print *,file1
      open(unit=1,file=file1,status='old')
      read(1,*) x1,xwc(1),xwc(2),xwc(3),xwc(4),xwc(5),xwc(6),xwc(7)
      read(1,*) x1,xhc(1),xhc(2),xhc(3),xhc(4),xhc(5),xhc(6),xhc(7)
      do i=1,199
         read(1,*,end=555) x1,x2,x3,x4,x5,x6,x7,x8
         xic(i)=x1
         xcc(i,1)=x2
         xcc(i,2)=x3
         xcc(i,3)=x4
         xcc(i,4)=x5
         xcc(i,5)=x6
         xcc(i,6)=x7
         xcc(i,7)=x8
      enddo
 555    continue
      close(1)

      return
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
