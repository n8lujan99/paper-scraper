
      parameter(nmax=1000)
      real z(nmax),he(nmax),he1(nmax),he2(nmax),heb(nmax)
      integer itd(nmax)
      character lab*10

      open(unit=1,file='herr.dat',status='old')
      read(1,*)
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,x5
         n=n+1
         z(n)=x1
         he(n)=x4
         heb(n)=x5
      enddo
 666  continue
      close(1)

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.4)
      call pgslw(3)
      call pgscf(2)

c      call pgenv(0.,3.5,0.,2.4,0,0)
c      call pgenv(0.,3.,0.,1.3,0,0)
      call pgenv(0.,3.,0.,2.7,0,0)
      call pgslw(4)
c      call pgline(n,z,he)
      call pgsls(4)
c      call pgline(n,z,heb)
      call pgsls(1)

      call pglabel('z','Percent Accuracy on Combined Distance','')
      call pgsci(4)

      ilc=0
      ilc=ilc+1
c      call pgsch(1.1)
c      call pgptxt(2.65,yp,0.,0.,"DEX k<0.2")
      call pgsch(2.5)
      ilc=ilc+1
      yp=2.3-0.1*ilc
      call pgsch(1.1)
c      call pgptxt(2.65,yp,0.,0.,"DEX k<0.3")
      call pgsch(2.5)
      ilc=ilc+1
      yp=2.3-0.1*ilc
      call pgsch(1.1)
c      call pgptxt(2.65,yp,0.,0.,"DEX +3pt")
      call pgsch(2.5)
      ilc=ilc+1
      yp0=0.4
      xp0=2.45
      yp=yp0-0.1*ilc
      call pgsch(1.2)
c      call pgptxt(2.65,yp,0.,0.,"DEX +IM")
      call pgptxt(xp0,0.35,0.,0.,"HETDEX")
      call pgsch(2.5)

      call pgsci(3)
      open(unit=1,file='euclid.h',status='old')
      read(1,*)
      n=0
      do i=1,1000
         read(1,*,end=778) x1,x2,x3,ic,it,ilab,lab
         n=n+1
         z(n)=x1
         he(n)=x2
      enddo
 778  continue
      close(1)

      call pgslw(5)
      call pgsci(ic)
c      call pgline(n,z,he)
      call pgslw(4)
      ilc=ilc+1
      yp=yp0-0.1*ilc
      call pgsch(1.2)
c      call pgptxt(xp0,0.27,0.,0.,lab)
      call pgsch(2.5)

      call pgsci(2)
      open(unit=1,file='desi.h',status='old')
      read(1,*)
      n=0
      do i=1,1000
         read(1,*,end=777) x1,x2,x3,ic,it,ilab,lab
         n=n+1
         z(n)=x1
         he(n)=x2
         itd(n)=it
         itype=it
      enddo
 777  continue
      close(1)

      call pgsci(2)
      call pgslw(4)
      call pgsch(2.0)
      do i=1,n
         call pgpt1(z(i),he(i),itd(i))
      enddo
c      call pgpt(n,z,he,itype)
      call pgsls(4)
      call pgslw(5)
c      call pgline(n,z,he)
      call pgsls(1)
c      call pgline(5,z,he)
      call pgslw(4)
      ilc=ilc+1
      yp=yp0-0.1*ilc
      call pgsch(1.2)
      call pgptxt(xp0,0.24,0.,0.,lab)
      call pgsch(2.5)

      fac=sqrt(1./5.)
      open(unit=1,file='cosmic_variance_distance.dat',status='old')
      read(1,*)
      n=0
      do i=1,nmax
         read(1,*,end=668) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,
     $        x13,x14,x15,x16,x17
         n=n+1
         z(n)=(x1+x2)/2.
         he(n)=x14*fac
      enddo
 668  continue
      close(1)    
      call pgsci(1)
      call pgsls(2)
      call pgslw(5)
      call pgline(n,z,he)
      call pgslw(4)

      call pgsci(4)
      call pgsls(1)
      call pgsch(1.2)
      call pgsch(4.0)
      call pgpt1(2.4,0.66,17)
c      y=0.66*1.3
c      call pgpt1(2.4,y,21)
      call pgsch(1.2)
      open(unit=1,file='herr2.dat',status='old')
      do i=1,1000
         read(1,*,end=555) x1,lab
c         call pgpt1(2.4,x1,8)
         xch=1.0+(0.2)*float(i-1)
         call pgsch(xch)
c         call pgarro(2.4,x1,2.4,x1-0.16)
c         call pgerr1(2,2.4,x1,0.,3.0)
      enddo
 555  continue
      close(1)

      call pgend
      end
