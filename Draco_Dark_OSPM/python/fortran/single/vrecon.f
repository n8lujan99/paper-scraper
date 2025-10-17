
      parameter(nmax=1000)
      real arr1(1032,112),arr2(1032,112),arr3(1032,112),arr4(1032,112)
      real xifu(50,50)
      character file1*100,file2*100
      character camp*2,cspecid*3,cifu*3,cifupos*3
      character cfib*50,corig*80
      
      nx=50
      ny=50
      dx=1

      read(*,*) file1
      open(unit=1,file=file1,status='old')
c      open(unit=1,file='list',status='old')

      do iall=1,nmax
         call setzero(arr1)
         call setzero(arr2)
         call setzero(arr3)
         call setzero(arr4)

         read(1,*,end=666) file1
         file2=file1
         call getfits(file1,arr1,cspecid,camp,cifu,cifupos,corig)
         write(file2(24:25),1001) "LU"
         file1=file2
         call getfits(file1,arr2,cspecid,camp,cifu,cifupos,corig)
         write(file2(24:25),1001) "RL"
         file1=file2
         call getfits(file1,arr3,cspecid,camp,cifu,cifupos,corig)
         write(file2(24:25),1001) "RU"
         file1=file2
         call getfits(file1,arr4,cspecid,camp,cifu,cifupos,corig)
         
         call getifucen(arr1,arr2,arr3,arr4,cspecid,cifupos,cifu,xifu,
     $        nx,ny,dx)

         call writefits(xifu,nx,ny,cifupos,corig)

      enddo
 666  continue
      close(1)
 1001 format(a2)

      end

      subroutine setzero(arr)
      real arr(1032,112)
      do i=1,1032
         do j=1,112
            arr(i,j)=0.
         enddo
      enddo
      return
      end

      subroutine getfits(file1,arr,cspecid,camp,cifu,cifupos,corig)
      real arr(1032,112)
      integer naxes(2)
      character file1*100
      character camp*2,cspecid*3,cifu*3,cifupos*3,corig*80
      logical simple,extend,anyf

      im1=0
      ier=0
      iext=1
      call ftgiou(im1,ier)
      iread=0
      call ftopen(im1,file1,iread,iblock,ier)
      if(ier.ne.0) then
c         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftmahd(im1,iext,ihd,ier)
      call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(im1,igc,0.,1032,ncol,nrow,arr,anyf,ier)

      call ftgkys(im1,"SPECID",cspecid,file1,ier)
      call ftgkys(im1,"AMP",camp,file1,ier)
      call ftgkys(im1,"IFUID",cifu,file1,ier)
      call ftgkys(im1,"IFUSLOT",cifupos,file1,ier)
      call ftgkys(im1,"NAME0",corig,file1,ier)
      if(ier.ne.0) then
         ier=0
         corig=""
      endif

      call ftclos(im1,ier)
 706  continue
c      if(ier.ne.0) print *,"No file for: ",file1
      return
      end

      subroutine getifucen(arr1,arr2,arr3,arr4,cspecid,cifupos,
     $        cifu,xifu,nx,ny,dx)

      real arr1(1032,112),arr2(1032,112),arr3(1032,112),arr4(1032,112)
      real xifu(50,50),xfib(112),xfiba(448),xpos(448),ypos(448)
      real xposa(448),yposa(448),xt(448),yt(448),xin(448),yin(44)
      integer ifib(112)
      character camp*2,cspecid*3,cifu*3,cifupos*3
      character file1*100

      wave1=400.
      wave2=800.

      if(cifu.eq."004") then
         file1="/work2/00115/gebhardt/maverick/scripts/back/ifu1.txt"
      else
         file1="/work2/00115/gebhardt/maverick/scripts/back/ifu0.txt"
      endif
      open(unit=2,file=file1,status='old')
      do i=1,448
         read(2,*) i1,x2,x3
         xpos(i)=x2+24.5
         ypos(i)=x3+24.5
      enddo
      close(2)

      if(cifu.eq.'003'.or.cifu.eq.'004'
     $     .or.cifu.eq.'005'.or.cifu.eq.'008') then
         do i=1,224
            xt(i)=xpos(449-i)
            yt(i)=ypos(449-i)
         enddo
         do i=1,224
            xpos(224+i)=xt(i)
            ypos(224+i)=yt(i)
         enddo
      endif
      if(cifu.eq.'007') then
         xp=xpos(38)
         yp=ypos(38)
         xpos(38)=xpos(39)
         ypos(38)=ypos(39)
         xpos(39)=xp
         ypos(39)=yp
      endif
      if(cifu.eq.'025') then
         x209=xpos(209)
         y209=ypos(209)
         x210=xpos(210)
         y210=ypos(210)
         x211=xpos(211)
         y211=ypos(211)
         x212=xpos(212)
         y212=ypos(212)
         x213=xpos(213)
         y213=ypos(213)
         x214=xpos(214)
         y214=ypos(214)
         xpos(209)=x210
         ypos(209)=y210
         xpos(210)=x211
         ypos(210)=y211
         xpos(211)=x212
         ypos(211)=y212
         xpos(212)=x214
         ypos(212)=y214
         xpos(213)=x209
         ypos(213)=y209
         xpos(214)=x213
         ypos(214)=y213
c         xpos(209)=x213
c         ypos(209)=y213
c         xpos(210)=x209
c         ypos(210)=y209
c         xpos(211)=x210
c         ypos(211)=y210
c         xpos(212)=x211
c         ypos(212)=y211
c         xpos(213)=x214
c         ypos(213)=y214
c         xpos(214)=x212
c         ypos(214)=y212
      endif
      if(cifu.eq.'030') then
         xp=xpos(446)
         yp=ypos(446)
         xpos(446)=xpos(447)
         ypos(446)=ypos(447)
         xpos(447)=xp
         ypos(447)=yp
      endif
      if(cifu.eq.'038') then
         xp=xpos(303)
         yp=ypos(303)
         xpos(303)=xpos(304)
         ypos(303)=ypos(304)
         xpos(304)=xp
         ypos(304)=yp
      endif
      if(cifu.eq.'041') then
         xp=xpos(252)
         yp=ypos(252)
         xpos(252)=xpos(253)
         ypos(252)=ypos(253)
         xpos(253)=xp
         ypos(253)=yp
      endif

      nfib=0
      call getfibb(arr1,xfib,wave1,wave2)
      do j=1,112
         nfib=nfib+1
         xfiba(nfib)=xfib(j)
         xposa(nfib)=xpos(225-j)
         yposa(nfib)=ypos(225-j)
      enddo

      call getfibb(arr2,xfib,wave1,wave2)
      do j=1,112
         nfib=nfib+1
         xfiba(nfib)=xfib(j)
         xposa(nfib)=xpos(113-j)
         yposa(nfib)=ypos(113-j)
      enddo

      call getfibb(arr3,xfib,wave1,wave2)
      do j=1,112
         nfib=nfib+1
         xfiba(nfib)=xfib(j)
         xposa(nfib)=xpos(337-j)
         yposa(nfib)=ypos(337-j)
      enddo

      call getfibb(arr4,xfib,wave1,wave2)
      do j=1,112
         nfib=nfib+1
         xfiba(nfib)=xfib(j)
         xposa(nfib)=xpos(449-j)
         yposa(nfib)=ypos(449-j)
      enddo

      do i=1,nfib
         xin(i)=xfiba(i)
      enddo
      call biwgt(xin,nfib,xb,xs)
      do i=1,nfib
         if(xfiba(i).lt.-4.*xs) xfiba(i)=0.
      enddo

      radm=6.0
      rfw=4.0
c      rfw=2.0
      rsig=rfw/2.35

      do i=1,nx
         do j=1,ny
            xifu(i,j)=0.
         enddo
      enddo

      itype=0
      if(itype.eq.1) then
      do i=1,nx
         xp=float(i-1)*dx
         do j=1,ny
            yp=float(j-1)*dx
            nin=0
            do k=1,nfib
               rad=sqrt((xp-xposa(k))**2+(yp-yposa(k))**2)
               if(rad.lt.radm) then
                  nin=nin+1
                  xin(nin)=rad
                  yin(nin)=xfiba(k)
               endif
            enddo
            if(nin.ge.2) then
               call gfit(nin,xin,yin,rsig,gaus)
            else
               gaus=0.
            endif
            xifu(i,j)=gaus
         enddo
      enddo

      else
      sumg=0.
      do i=1,nx
         xp=float(i-1)*dx
         do j=1,ny
            yp=float(j-1)*dx
            rad=sqrt((xp-24.5)**2+(yp-24.5)**2)
            if(rad.lt.radm) then
               w=rad/rsig
               gaus=(exp(-w*w/2.))
               sumg=sumg+gaus
            endif
         enddo
      enddo
      fac=float(nx*ny)/float(nfib)
      sumg=sumg/fac

      sum1=0.
      do k=1,nfib
         flux=xfiba(k)
         sum1=sum1+flux
         if(flux.ne.0) then
            do i=1,nx
               xp=float(i-1)*dx
               do j=1,ny
                  yp=float(j-1)*dx
                  rad=sqrt((xp-xposa(k))**2+(yp-yposa(k))**2)
                  if(rad.lt.radm) then
                     w=rad/rsig
                     gaus=(exp(-w*w/2.))/sumg
                     xifu(i,j)=xifu(i,j)+gaus*flux
                  endif
               enddo
            enddo
         endif
      enddo
      endif

      return
      end

      subroutine writefits(xifu,nx,ny,cifupos,corig)
      real xifu(50,50)
      integer naxes(2)
      character cname*40,cifupos*3,corig*80

c      cname="Co_"//cifupos//".fits"
c            CoFeS20190121T040040.1_085_sci.fits
      if(corig(1:1).ne."2") then
         cname="Coaaa20000000T000000.0_"//cifupos//"_sci.fits"
      else
         cname="CoFeS"//corig(1:18)//cifupos//"_sci.fits"
      endif

      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      iblock=1
      igc=0
      ier=0
      call ftinit(50,cname,iblock,ier)
      call ftphps(50,-32,naxis,naxes,ier)
      call ftp2de(50,igc,nx,naxes(1),naxes(2),xifu,ier)
      call ftclos(50,ier)

      return
      end

      subroutine getfibb(arr,xfib,w1,w2)
      real arr(1032,112),xfib(112),xin(1032)

      do j=1,112
         ilo=nint(w1)
         ihi=nint(w2)
         nin=0
         do i=ilo,ihi
            nin=nin+1
            xin(nin)=arr(i,j)
         enddo
         call biwgt(xin,nin,xb,xs)
         xfib(j)=xb
      enddo
      return
      end

      subroutine gfit(nin,xin,yin,rsig,gaus)
      real xin(nin),yin(nin)

      call sort2(nin,xin,yin)
      sum=0.
      ntot=min(3,nin)
      do i=1,ntot
         w=xin(i)/rsig
         gp=(exp(-w*w/2.))
         sum=sum+yin(i)/gp
      enddo
      gaus=sum/float(ntot)

      return
      end
