c     Karl's plotter program modified by dor
c     reads new demographics file
c     makes a four-panel plot with evocative symbols.
c
      call readem(n) 
      pause 
c      
      print *, ' plotting M-M relation' 
      call pgbegin(0,'?',1,1)
c     set width as big as possible,aspect=1.
      call pgpap(0.,1.)
c     top panel 
      call pgpanl (1, 1)
c
      call pgvport(0.15,0.55,0.60, 1.00)
c     plot the L correlation 
      call  lplot (1, n )
c
c - now plot the sigma correlation
c
      print *, ' plotting M-sig relation' 
      call pgvport(0.55,0.95,0.60, 1.00)
      call sigplot (1, n) 
c
c     bottom panel 
c
      call pgpanl (1, 2)
c
c     plot L correlation again

      call pgvport(0.15,0.55,0.10, 0.50)
      call  lplot (2, n )
c
c - now plot the sigma correlation again 
c
      call pgvport(0.55,0.95,0.10, 0.50)
      call sigplot (2, n) 
      call pgend
      stop 
      end
c
      function icolor(method)
      dimension icols(20)
      data icols/10*0, 2, 11, 1, 9, 8, 5*0/
c     11 - stars, 12 - gas, 13 - masers, 14 - reverb, 15-2i, 16-bertola
c        red,       blue,      black,     goosedung,  orange 
c
      icolor = icols(method)
c      print *, ' method: ', method, ' color: ', icolor
      return
      end
c
      integer function marker(htype)
      character*4  htype
c
      if (htype(1:1).eq.'E') then 
         marker = -8 
         return
      endif
c
c     S0's are tricky: S0 or SB0
      if(htype(1:2).eq.'S0') then 
         marker = 2269
         return 
      endif
      if(htype(1:3).eq.'SB0') then 
         marker = 2269
         return 
      endif
c
c     all other spirals
      if (htype(1:1).eq.'S') then 
         marker = 2268
         return
      endif
c     any remaining cases
      marker = 18
      print *, ' in marker, failed to find common hubble type'
      return
      end
c
c
      subroutine lplot (itime, n)
c
      parameter(nmax=1000)
      real m, xl(2), yl(2)
      character cname*9, htype*4, type*4 
      common/plotdata/ xlbulge(nmax), xltot(nmax), ttype(nmax), 
     1     bh(nmax), bhlo(nmax), bhhi(nmax), 
     2     sigma(nmax), errsigma(nmax), 
     3     distance(nmax), re(nmax), xml(nmax), 
     4     bherr(nmax), meth(nmax), idetect(nmax) 
      common/names/ cname(nmax), htype(nmax)
c
      xmin=log10(1.0e8)
      xmax=log10(3.e11)
      ymin=log10(8.e5)
      ymax=log10(6.e9)
c
c     font, height, width 
c      call pgscf(2)
      call pgsch(1.3)
      call pgslw(2.5)
      call pgsci(1)
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bclnst',0.,0,'bclnst',0.,0)
c      call pglabel('','M\\D\\(0850)\\U (M\\D\\(2281)\\U)','')
      call pgmtxt('L', 2.3,0.5,0.5, 
     1     'M\\D\\(0850)\\U (M\\D\\(2281)\\U)')
      if (itime.eq.2) then 
         call pgmtxt('B',2.3,0.5,0.5,'L\\DB\\U(bulge)/L\\D\\(2281)\\U')
      endif

      slope=0.1
      xl(1)=xmin
      xl(2)=xmax
c Magorrian
      yl(1)=log10(5.e-3)-1.11+1.2*(xl(1))
      yl(2)=log10(5.e-3)-1.11+1.2*(xl(2))
c      call pgline(2,xl,yl)
c KR95
      yl(1)=log10(2.5e-3)-1.11+1.2*(xl(1))
      yl(2)=log10(2.5e-3)-1.11+1.2*(xl(2))
c      call pgline(2,xl,yl)
c QSO light
      yl(1)=-4.349+1.2*xl(1)
      yl(2)=-4.349+1.2*xl(2)
      call pgsls(4)
      call pgline(2,xl,yl)
      call pgsls(1)

      ns=0
      bbit=0.3
      call pgsah(1,50.,0.5)
c     
      if (itime.eq.1) then 
         do 200 i=1,n
c     methods to skip: 
            if(meth(i).eq.15) go to 199
            type = htype(i)
            ipt1 = marker(type)
            call pgpt1(xlbulge(i),bh(i),ipt1)
c
         if(idetect(i).eq.1) then
            call pgsch(1.1)
            print *, ' -----------------------------------------'
            print *,  cname(i)
            print *, cname(i), xlbulge(i), sigma(i), bh(i), ipt1
            call pgslw(2.5)
            if(bhlo(i).le.0) then
               call pgsch(0.6)
               print *,  ' plot nondetect: ', cname(i), xlbulge(i),bh(i)
               xlb = xlbulge(i)
               call pgarro(xlb,bh(i),xlb,bh(i)-bbit)
            endif
         endif
 199     continue
 200  continue
      endif
      if (itime.eq.2) then 
         do 300 i=1,n
c     methods to skip: 
            if(meth(i).eq.15) go to 299
c
         if(idetect(i).eq.1) then
            call pgsch(1.1)
            print *, ' -----------------------------------------'
            print *,  cname(i)
            print *, cname(i), xlbulge(i), sigma(i), bh(i), ipt1
            call pgslw(2.5)
            call pgpt1(xlbulge(i), bh(i), meth(i))
            if(bhlo(i).gt.0) then
               call pgerry(1,xlbulge(i),bhhi(i),bhlo(i),1.)
            else
               call pgsch(0.6)
               print *,' plot nondetect: ',cname(i),xlbulge(i),bh(i)
               xlb = xlbulge(i)
               call pgarro(xlb,bh(i),xlb,bh(i)-bbit)
            endif
         endif
 299     continue
 300  continue
      endif
      return
      end
c
c
      subroutine sigplot (itime, n) 
c
      parameter(nmax=1000)
      real m, xl(2), yl(2)
      character cname*9, htype*4, type*4 
      common/plotdata/ xlbulge(nmax), xltot(nmax), ttype(nmax), 
     1     bh(nmax), bhlo(nmax), bhhi(nmax), 
     2     sigma(nmax), errsigma(nmax), 
     3     distance(nmax), re(nmax), xml(nmax), 
     4     bherr(nmax), meth(nmax), idetect(nmax) 
      common/names/ cname(nmax), htype(nmax)
c
      xmin=log10(50.)
      xmax=log10(440.)
      ymin=log10(8.e5)
      ymax=log10(6.e9)
      call pgsci(1)
      call pgsch(1.3)
c
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bclnst',0.,0,'bclst',0.,0)
      if (itime.eq.2) then 
         call pgmtxt('B',2.3,0.5,0.5,'\\gs\\De\\U  (km/s)')
      endif

      bbit=0.3
      call pgsah(1,50.,0.5)
c
      if (itime.eq.1) then 
         do 200 i =1, n
c     methods to skip: 
            if(meth(i).eq.15) go to 199
c
            if(idetect(i).eq.1) then
               type = htype(i)
               print *, ' -----------------------------------------'
               print *,  cname(i)
               ipt1 = marker(type)
               call pgsch(1.1)
               print *, cname(i), xlbulge(i), sigma(i), bh(i),type,ipt1
               slog = alog10(sigma(i))
               call pgpt1(slog,bh(i),ipt1)
               call pgslw(2.5)
               if(bhlo(i).le.0) then
                  call pgsch(0.6)
                  print *,  ' plot nondetect: ',cname(i),sigma(i),bh(i)
                  call pgsci(icolor(meth(i)))
                  call pgarro(slog,bh(i),slog,bh(i)-bbit)
               endif
            endif
c
 199        continue
 200     continue
c     
      endif
c
      if (itime.eq.2) then 
c
         do 300 i =1, n
c     methods to skip: 
            if(meth(i).eq.15) go to 299
c
            if(idetect(i).eq.1) then
               print *, ' -----------------------------------------'
               print *,  cname(i)
               type = htype(i)
               call pgsch(1.1)
               print *, cname(i), xlbulge(i),sigma(i),bh(i), type, ipt1
               slog = alog10(sigma(i))
               call pgpt1(slog,bh(i),meth(i))
               call pgslw(2.5)
               if(bhlo(i).gt.0) then
                  call pgerry(1,slog,bhhi(i),bhlo(i),1.)
               endif
            else
               call pgsch(0.6)
               print *,  ' plot nondetect: ', cname(i), sigma(i), bh(i)
               call pgsci(icolor(meth(i)))
               slog = alog10(sigma(i))
               call pgarro(slog,bh(i),slog,bh(i)-bbit)
            endif
c
 299        continue
 300     continue
      endif
c
      if (itime.eq.1) call legend1 (9.6, 1.75, .25, .05)
      if (itime.eq.2) call legend2 (9.6, 1.75, .25, .05)
      return
      end 
c
      subroutine legend2 (top, xleft, delta, hspace)
c     detectable bh mass, etc
      print *, ' in legend' 
      call pgsch(1.30)
      call pgsci(1)
      xp = xleft
      xl = xleft+hspace
      spaceleft = xp -2.*hspace 
c     traditional bulges, circle 
      y = top - delta 
      yl = y - delta/4.
      call pgpt (1, xp, y, 11)
      call pgtext (xl, yl,  'Stars')
c
      y = y - delta 
      yl = y - delta/4.
      call pgpt (1, xp, y, 12)
      call pgtext (xl, yl,  'Gas')
c
      y = y - delta 
      yl = y - delta/4.
      call pgpt (1, xp, y,  13)
      call pgtext (xl, yl,  'Masers')
c
      return
      end
c
      subroutine legend1 (top, xleft, delta, hspace)
c     detectable bh mass, etc
      print *, ' in legend' 
      call pgsch(1.30)
      call pgsci(1)
      xp = xleft
      xl = xleft+hspace
      spaceleft = xp -2.*hspace 
c     traditional bulges, circle 
      y = top - delta 
      yl = y - delta/4.
      call pgpt (1, xp, y, -8)
      call pgtext (xl, yl,  'Elliptical')
c
      y = y - delta 
      yl = y - delta/4.
      call pgpt (1, xp, y, 2269)
      call pgtext (xl, yl,  'S0')
c
      y = y - delta 
      yl = y - delta/4.
      call pgpt (1, xp, y,  2268)
      call pgtext (xl, yl,  'Spirals')
c
      return
      end
c
      subroutine readem(n)
c
      parameter(nmax=1000)
      character a1*120,a2*120,a3*120,a4*120
      character cname*9, htype*4, type*4 
      common/plotdata/ xlbulge(nmax), xltot(nmax), ttype(nmax), 
     1     bh(nmax), bhlo(nmax), bhhi(nmax), 
     2     sigma(nmax), errsigma(nmax), 
     3     distance(nmax), re(nmax), xml(nmax), 
     4     bherr(nmax), meth(nmax), idetect(nmax) 
c
      common/names/ cname(nmax), htype(nmax)
c
      print *, ' input from file dem.dat' 
      open(unit=1,file='dem.dat',status='old')

      n = 0
      do 100 i=1,nmax
         read(1,1001,end=101) a1
         print *, a1 
         if(a1(1:1).ne.'#') then
c
            n = n + 1
c
            read(a1,1002) a2,a3,a4
            read(a4,*) x1, x2,x3,x4,xi1,x5,x6,x7, x8, x9, x10, x11 
            cname(n)=a2
            htype(n) = a3
            print *, cname(n), htype(n)
            print *, '        ', x1, x2, x3, x4, xi1, x5, x6, 
     1           x7, x8, x9, x10, x11 
c
c           x1 is ttype 
            ttype(n) = x1
c
c           x2 is probably blue bulge magnitude
c           xlbulge is the log of the blue bulge luminosity
c
            if(x2.lt.0.) then
               xlbulge(n)=-0.4*(x2-5.48)
            else
               xlbulge(n)=alog10(x2)
            endif
c
c                x3 is TOTAL magnitude, probably blue 
c                ltot is the log of the total  luminosity
c
            if(x3.lt.0.) then
               xltot(n)=-0.4*(x3-5.48)
            else
               xltot(n)=alog10(x3)
            endif
c
            if(x4.gt.0) then
c                 x4 is BH mass, gt.0 is successful detection 
               idetect(n)= 1
               bh(n)=alog10(x4)
            else
c                 nondetection 
               x4=-x4
               idetect(n)= 0
               bh(n)=alog10(x4)
            endif
c
c     Method: 
            meth(n)=xi1
c
            if(x5.gt.0) then
c                  1 sigma upper and lower errbars on bh mass 
               bhlo(n)=alog10(x5)
               bhhi(n)=alog10(x6)
               bherr(n)=alog10(x6-x5)/2.
            else
               bhlo(n)=x5
               bhhi(n)=x6
               bherr(n)=0.
            endif
c
            sigma(n)=x7
            errsigma(n) = x8 
            distance(n) = x9
            re(n) = x10
            xml(n) = x11
         endif
 100  continue
c
 101  continue
      print *, '   ' 
      print *, ' *** conclusion of input' 
c
      close(1)
      print *, ' stored quantities ' 
      write (6, 1004) 
      do 200 i = 1, n 
c         print *, i, cname(i), htype(i), xlbulge(i), xltot(i), ttype(i), 
c     1        bh(i), bhlo(i), bhhi(i), sigma(i), errsigma(i), bherr(i), 
c     2        meth(i), idetect(i), distance(i), re(i), xml(i)
         write (6, 1003) i, cname(i), htype(i), xlbulge(i), 
     1        xltot(i), ttype(i), 
     2        bh(i), bhlo(i), bhhi(i), bherr(i), sigma(i), errsigma(i), 
     3        meth(i), idetect(i), distance(i), re(i), xml(i)
 200  continue
      write (6, 1004) 
c
      return
 1001 format(a120)
 1002 format(a9,a3,a120)
 1003 format(1x, i4, 1x, a9, 1x, a4, 1x, f7.3, 
     1     1x, f7.3, 1x, f6.1, 
     2     6(1x, f8.3), 
     3     2(1x, i4) , 3(1x, f8.3)) 
 1004 format (2x, 'index name type     Lbulge  Ltotl    ttype',  
     1  '     bh', 5x, ' bhlo    bhhi   bherr  sigma   errsig  ', 
     2  ' method  detct distnc     Re     M/L' )
      end 

