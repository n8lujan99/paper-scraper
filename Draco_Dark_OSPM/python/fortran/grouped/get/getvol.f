Marked
      parameter(nmax=10000)
      real f50(nmax),zl(nmax),xl(nmax),z(nmax)
      real zv(nmax),vol(nmax),xlum(nmax),sn(nmax)
      real xall(10,13),afsn(nmax),asn(nmax)
      real anew(nmax),snew(nmax)
      real*8 dx2
      integer*8 i1

      open(unit=1,file='f50.dat',status='old')
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
      enddo
 668  continue
      close(1)

      open(unit=1,file='sn_all.sav',status='old')
      read(1,*) x1,f50s
      nsa=0
      do i=1,100
         read(1,*,end=669) x1,x2
         nsa=nsa+1
         afsn(nsa)=x1
         asn(nsa)=x2
      enddo
 669  continue
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
         do ilum=1,nlum
            x0=xlum(ilum)
            nt0a=0
            nt1a=0
            v0=0.
            v1=0.
            
            do iz=1,nz2
               z0=z(iz)
               call xlinint(z0,nz,zl,xl,xconv)
               call xlinint(z0,nv,zv,vol,vol0)
               f0=10**(x0-40.)/xconv*1.e17
c               f0=f0*0.6
               f0=f0*0.2
               f1=f0*4.8/sn0
               nt0=0
               sum0=0.
               nt1=0
               sum1=0.
               do i=1,ns
                  ratf=f50(i)/f50s
                  do ia=1,nsa
                     anew(ia)=ratf*afsn(ia)
                     snew(ia)=asn(ia)
                  enddo
                  call xlinint(f0,nsa,anew,snew,s0)
                  call xlinint(f1,nsa,anew,snew,s1)
                  if(f0.gt.f50(i)) nt0=nt0+1
                  if(f1.gt.f50(i)) nt1=nt1+1
                  sum0=sum0+s0
                  sum1=sum1+s1
               enddo
               rat=0.
               if(nt1.gt.0) rat=float(nt0)/float(nt1)
               nt0a=nt0a+nt0
               nt1a=nt1a+nt1
               v0=v0+float(nt0)
               v1=v1+float(nt1)
               v0=v0+vol0*float(nt0)
               v1=v1+vol0*float(nt1)
c               v0=v0+vol0*sum0
c               v1=v1+vol0*sum1
            enddo
            rat=0.
            if(v1.gt.0) rat=v0/v1
            xall(isn,ilum)=rat
         enddo
      enddo

      do ilum=1,nlum
         print *,xlum(ilum),(xall(i,ilum),i=1,nsn)
      enddo

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
