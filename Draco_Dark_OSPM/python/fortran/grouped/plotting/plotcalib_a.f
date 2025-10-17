
      parameter(nmax=10000)
      real xin(nmax),yin(nmax),yin2(nmax)
      real xv1(nmax),xv2(nmax),yc(nmax),xl(2),yl(2)
      integer ic(nmax)
      character file1*30,c3*30,cdate*6,af(50)*6
      character a3*6,a4*3,a5*3,a6*2,aold*3
      character adate(nmax)*6,aspec(nmax)*3,acon(nmax)*3,aamp(nmax)*2

      af(1)='201701'
      af(2)='201702'
      af(3)='201703'
      af(4)='201704'
      af(5)='201705'
      af(6)='201706'
      af(7)='201707'
      af(8)='201708'
      af(9)='201709'
      af(10)='201710'
      af(11)='201711'
      af(12)='201712'
      af(13)='201801'
      af(14)='201802'
      af(15)='201803'
      af(16)='201804'
      af(17)='201805'
      af(18)='201806'
      af(19)='201807'
      af(20)='201808'
      af(21)='201809'
      af(22)='201810'
      af(23)='201811'
      af(24)='201812'
      af(25)='201901'
      af(26)='201902'
      af(27)='201903'
      af(28)='201904'
      af(29)='201905'
      af(30)='201906'
      af(31)='201907'
      af(32)='201908'
      af(33)='201909'
      af(34)='201910'
      af(35)='201911'
      af(36)='201912'
      af(37)='202001'
      nmth=37

      call pgbegin(0,'?',2,2)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

c      call pgenv(0.,37.,980.,1001.,0,0)
c      call pgenv(0.,37.,988.,996.,0,0)
c      call pglabel("months since 201701","delta pixel","")

      open(unit=1,file='listin',status='old')
      ntot=0
      do i=1,10000
         read(1,*,end=777) x1,x2,a3,a4,a5,a6
         ntot=ntot+1
         xv1(ntot)=x1
         xv2(ntot)=x2
         adate(ntot)=a3
         aspec(ntot)=a4
         acon(ntot)=a5
         aamp(ntot)=a6
      enddo
 777  continue
      rewind(1)

      open(unit=1,file='listin2',status='old')
      do i=1,nmax
         read(1,*,end=666) a3
         nin=0
         nc=0
         ymin=1e10
         ymax=-ymin
         do j=1,ntot
            if(a3.eq.aspec(j)) then
               nin=nin+1
               yin(nin)=xv1(j)
               if(aamp(j).eq.'LL') ic(nin)=1
               if(aamp(j).eq.'LU') ic(nin)=2
               if(aamp(j).eq.'RL') ic(nin)=3
               if(aamp(j).eq.'RU') ic(nin)=4
               ymin=min(ymin,yin(nin))
               ymax=max(ymax,yin(nin))
               do k=1,nmth
                  if(af(k).eq.adate(j)) then
                     xin(nin)=float(k)
                  endif
               enddo
               if(nin.eq.1) aold=acon(j)
               if(acon(j).ne.aold) then
                  aold=acon(j)
                  nc=nc+1
                  yc(nc)=xin(nin)
               endif
            endif
         enddo
         ybit=(ymax-ymin)/10.
         ymin=ymin-ybit
         ymax=ymax+ybit
         call pgenv(0.,38.,ymin,ymax,0,0)
         call pglabel("months since 201701",
     $        "Relative Normalization at 3700","")
         call pgmtxt('T',1.0,0.5,0.5,"IFUSLOT= "//a3)
         do k=1,nin
            call pgsci(ic(k))
            call pgpt1(xin(k),yin(k),17)
         enddo
         call pgsci(1)
         do k=1,nc
            xl(1)=yc(k)
            xl(2)=yc(k)
            yl(1)=ymin
            yl(2)=ymax
            call pgline(2,xl,yl)
         enddo
      enddo
 666  continue

      call pgend

      end
