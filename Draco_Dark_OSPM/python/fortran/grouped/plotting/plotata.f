
      real xsp(100),xif(100),xin(100),yvl(100),yvr(100)
      integer ispec(100),ifua(100)
      character file1*14,c2*3,c3*3,cname(100)*7,lab*7

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)
      call pgenv(0.,61.,0.7,1.45,0,0)
      call pglabel("Index","Relative Throughput","")

      open(unit=1,file='list',status='old')

      do i=1,100
         read(1,*,end=666) file1,c2,c3
         read(file1(2:4),1001) ifu
         cname(i)=c2//"_"//c3
         read(c2,1001) ispec(i)
         read(c3,1001) ifua(i)
         xsp(i)=float(ispec(i))
         xif(i)=float(ifua(i))
         xin(i)=float(i)
         call getw(file1,y1,y2,y3)
         yvl(i)=y2
         write(file1(6:7),1002) "RL"
         call getw(file1,y1,y2,y3)
         yvr(i)=y2
         nt=i
      enddo
 666  continue
      close(1)

      call sort2(nt,xsp,xin)
c      call sort2(nt,xif,xin)

      itop=0
      do i=1,nt
         ip=nint(xin(i))
         call pgsci(1)
         call pgpt1(float(i),yvl(ip),17)
         call pgsci(2)
         call pgpt1(float(i),yvr(ip),17)
         call pgsci(1)
         xpos=float(i)+0.5
         ypos=0.75
         itop=itop+1
         if(itop.eq.2) then
            ypos=1.39
            itop=0
         endif
         lab=cname(ip)
         call pgsch(0.7)
         call pgptxt(xpos,ypos,90.,0.5,lab)
         call pgsch(1.5)
      enddo


      call pgend
 1001 format(i3)
 1002 format(a2)

      end

      subroutine getw(file1,y1,y2,y3)
      real xin1(2000),xin2(2000),xin3(2000)
      character file1*14

      w1=3700.
      w2=4100.
      w3=4400.
      w4=4800.
      w5=5100.
      w6=5400.

      open(unit=2,file=file1,status='old')
      n1=0
      n2=0
      n3=0
      do j=1,2000
         read(2,*,end=667) x1,x2
         if(x1.gt.w1.and.x1.lt.w2) then
            n1=n1+1
            xin1(n1)=x2
         endif
         if(x1.gt.w3.and.x1.lt.w4) then
            n2=n2+1
            xin2(n2)=x2
         endif
         if(x1.gt.w5.and.x1.lt.w6) then
            n3=n3+1
            xin3(n3)=x2
         endif
      enddo
 667  continue
      close(2)
      call biwgt(xin1,n1,y1,xs)
      call biwgt(xin2,n2,y2,xs)
      call biwgt(xin3,n3,y3,xs)

      return
      end
