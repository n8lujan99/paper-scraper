
      parameter(nmax=1000000)
      real xi(nmax),xo(nmax),xl(10000),xu(10000)
      real xin(nmax)
      character file1*5

      ximin=9.75
      ximax=30.25
      xstep=0.5
      nstep=nint((ximax-ximin)/xstep)+1
      do i=1,nstep-1
         xl(i)=ximin+(ximax-ximin)*float(i-1)/float(nstep-1)
         xu(i)=ximin+(ximax-ximin)*float(i)/float(nstep-1)
      enddo

      open(unit=1,file='listts',status='old')

      do iall=1,1000
         read(1,*,end=666) file1
         open(unit=2,file=file1,status='old')
         n=0
         do j=1,nmax
            read(2,*,end=667) x1,x2
            n=n+1
            xi(n)=x1
            xo(n)=x2
         enddo
 667     continue
         close(2)

         do j=1,nstep-1
            nin=0
            do i=1,n
               if(xi(i).ge.xl(j).and.xi(i).lt.xu(j)) then
                  nin=nin+1
                  xin(nin)=xo(i)
               endif
            enddo
            call biwgt(xin,nin,xb,xs)
            xcen=(xl(j)+xu(j))/2.
c            print *,xcen,xb,xs,nin," ",file1
            print *,xcen,xcen-xb,xs/xb,nin," ",file1
         enddo
      enddo
 666  continue
      close(1)

      end

