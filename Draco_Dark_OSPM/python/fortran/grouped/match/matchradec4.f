
      parameter(nmax=1000000)
      real w(nmax),xf(nmax),sn(nmax),sb(nmax),sr(nmax)
      real*8 dra(nmax),ddec(nmax),drad,dx1,dx2
      integer iok(nmax),iok2(nmax),id(nmax)
      character aname(nmax)*12,a6*11

      parameter(radtodeg=57.29578)

      rad=2.5

      open(unit=1,file='d4',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) dx1,dx2,i5,a6
         n=n+1
         dra(n)=dx1
         ddec(n)=dx2
         id(n)=i5
         aname(n)=a6(1:8)//"v"//a6(9:11)
         iok(n)=1
      enddo
 666  continue
      close(1)

      cosd=cos(sngl(ddec(1))/radtodeg)

      do i=1,n-1
         if(iok(i).eq.1) then
            do j=i+1,n
               drad=(dble(cosd)*(dra(i)-dra(j)))**2+(ddec(i)-ddec(j))**2
               drad=3600.d0*dsqrt(drad)
               radc=sngl(drad)
               if(radc.lt.rad) then
                  iok(j)=0
               endif
            enddo
         endif
      enddo

      open(unit=11,file='radec.out',status='unknown')
      do i=1,n
         if(iok(i).eq.1) then
            write(11,1011) dra(i),ddec(i),id(i),aname(i)
         endif
      enddo
      close(11)
 888  continue

 1011 format(2(1x,f11.7),1x,i10,1x,a12)
      end
