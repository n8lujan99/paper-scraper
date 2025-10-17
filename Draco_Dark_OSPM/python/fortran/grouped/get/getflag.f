Skipped
      parameter(nmax=100000)
      real ra(nmax),dec(nmax),wave(nmax),sn(nmax),flux(nmax)
      integer iraw(nmax),jraw(nmax),iuse(nmax)
      integer ixs(nmax),ixe(nmax),iys(nmax),iye(nmax)
      character a5*3,a6*8,a7*3,a9*12,a10*24,c20*20
      character csim(nmax)*3,cdate(nmax)*8,cshot(nmax)*3
      character cfield(nmax)*12,cfib(nmax)*24,aflag(nmax)*24
      character aflag1(nmax)*20

c     flag1 is 5200, which is to remove 5194-5197 and 5200-5205
c     flag2 is remove 5456-5466      
      open(unit=1,file='j1',status='old')
      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3,x4,a5,a6,a7,x8,a9,a10,i11,i12
         n=n+1
         ra(n)=x1
         dec(n)=x2
         wave(n)=x3
         sn(n)=x4
         csim(n)=a5
         cdate(n)=a6
         cshot(n)=a7
         flux(n)=x8
         cfield(n)=a9
         cfib(n)=a10
         iraw(n)=i11
         jraw(n)=i12
         iuse(n)=1
         if(wave(n).gt.3534..and.wave(n).lt.3556.) iuse(n)=0
      enddo
 666  continue
      close(1)

c this is bad fiber
      open(unit=1,file='flag0',status='old')
      nf=0
      do i=1,nmax
         read(1,*,end=667) a10
         nf=nf+1
         aflag(nf)=a10
      enddo
 667  continue
      close(1)
      do i=1,n
         do j=1,nf
            if(cfib(i).eq.aflag(j)) iuse(i)=0
         enddo
      enddo

c     flag1 is 5200, which is to remove 5194-5197 and 5200-5205
      open(unit=1,file='flag1',status='old')
      nf=0
      do i=1,nmax
         read(1,*,end=668) c20
         nf=nf+1
         aflag1(nf)=c20
      enddo
 668  continue
      close(1)
      do i=1,n
         do j=1,nf
            c20=cfib(i)
            if(wave(i).gt.5194..and.wave(i).lt.5197.) then
               if(c20.eq.aflag1(j)) iuse(i)=0
            endif
            if(wave(i).gt.5200..and.wave(i).lt.5205.) then
               if(c20.eq.aflag1(j)) iuse(i)=0
            endif
         enddo
      enddo

c     flag2 is remove 5456-5466      
      open(unit=1,file='flag2',status='old')
      nf=0
      do i=1,nmax
         read(1,*,end=669) c20
         nf=nf+1
         aflag1(nf)=c20
      enddo
 669  continue
      close(1)
      do i=1,n
         do j=1,nf
            c20=cfib(i)
            if(wave(i).gt.5456..and.wave(i).lt.5466.) then
               if(c20.eq.aflag1(j)) iuse(i)=0
            endif
         enddo
      enddo

c     this is bad pixel in raw
      open(unit=1,file='badpix.list',status='old')
      nf=0
      do i=1,nmax
         read(1,*,end=670) c20,i2,i3,i4,i5
         nf=nf+1
         aflag1(nf)=c20
         ixs(nf)=i2
         ixe(nf)=i3
         iys(nf)=i4
         iye(nf)=i5
      enddo
 670  continue
      close(1)
      do i=1,n
         c20=cfib(i)
         ix=iraw(i)
         iy=jraw(i)
         do j=1,nf
            if(ix.gt.ixs(j).and.ix.le.ixe(j).and.
     $         iy.gt.iys(j).and.iy.le.iye(j)) then
               if(c20.eq.aflag1(j)) iuse(i)=0
            endif
         enddo
      enddo

      open(unit=11,file='out',status='unknown')
      do i=1,n
         if(iuse(i).eq.1) write(11,1101) ra(i),dec(i),wave(i),sn(i),
     $        csim(i),cdate(i),cshot(i),flux(i),cfib(i)
      enddo
      close(11)

 1101 format(2(f9.5,1x),f8.3,1x,f6.2,1x,a3,1x,a8,1x,a3,1x,f8.3,1x,a24)
      end
