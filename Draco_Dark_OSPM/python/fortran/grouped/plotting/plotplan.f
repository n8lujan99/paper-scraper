
      parameter(nmax=100000)
      real x(nmax),y(nmax),ye(nmax),ya(nmax,100),xin(100)
      real yel(nmax),yeu(nmax),ydiff(nmax),ydiffn(nmax)
      real xl(2),yl(2),xtril(100),y2p(nmax),x2p(nmax)
      real xpred(nmax),ypred(nmax)
      character file1*80,file2*80,c1*3,ctri(100)*4

      xtri=121./365.
      xtril(1)=20./365.+2018.
      xtril(2)=xtril(1)+xtri
      xtril(3)=xtril(2)+xtri
      xtril(4)=xtril(3)+xtri
      xtril(5)=xtril(4)+xtri
      xtril(6)=xtril(5)+xtri
      xtril(7)=xtril(6)+xtri
      xtril(8)=xtril(7)+xtri
      xtril(9)=xtril(8)+xtri
      xtril(10)=xtril(9)+xtri
      xtril(11)=xtril(10)+xtri
      ctri(1)="18-1"
      ctri(2)="18-2"
      ctri(3)="18-3"
      ctri(4)="19-1"
      ctri(5)="19-2"
      ctri(6)="19-3"
      ctri(7)="20-1"
      ctri(8)="20-2"
      ctri(9)="20-3"
      ctri(10)="21-1"
      ctri(11)="21-2"

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)
      call pgslw(2)

c - diff from Jan 1 to Sept 1 is 122 days
c   let's center on Feb 1, 2018, which is mjd=58150.3
      x0=58130.3
      xt1=58211.3-x0
      xt2=58331.2-x0
      xt3=58452.2-x0
      xt4=xt3+121.
      xt5=xt4+121.
      xt6=xt5+121.
      xt7=xt6+121.
      xt8=xt7+121.
      xt9=xt8+121.
      xt10=xt9+121.
      xt11=xt10+121.
      
      xt1=xt1/365.+2018
      xt2=xt2/365.+2018
      xt3=xt3/365.+2018
      xt4=xt4/365.+2018
      xt5=xt5/365.+2018
      xt6=xt6/365.+2018
      xt7=xt7/365.+2018
      xt8=xt8/365.+2018
      xt9=xt9/365.+2018
      xt10=xt10/365.+2018
      xt11=xt11/365.+2018

      xmin=0.+2018
      xmax=2022.0
      ymin=0.
      ymax=4200.
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      call pglabel('Year','Cumulative Number of Observations','')

      open(unit=1,file='inall',status='old')
      n=0
      n2=0
      do i=1,100000
         read(1,*,end=667) x1
         x1=x1-7
         if(x1.gt.0) then
            n=n+1
            x(n)=x1
            x(n)=x(n)/365.+2018.
            xcount=float(n)/9.
            fac=0.9
            y(n)=xcount*fac
            xpred(n)=x(n)
            ypred(n)=y(n)
            if(x(n).gt.2018.9) then
               n2=n2+1
               x2p(n2)=x(n)
               y2p(n2)=y(n)-250.
            endif
         endif
      enddo
 667  continue
      close(1)
      npred=n
      call pgslw(4)
      call pgsls(1)
      call pgline(n,x,y)
      call pgsls(1)
c      call pgline(n2,x2p,y2p)
      open(unit=11,file='out',status='unknown')
      do i=1,n
         write(11,*) x(i),y(i)
      enddo
      close(11)

      open(unit=1,file='dex.dat',status='old')
      n=0
      fac=0.95
      do i=1,100000
         read(1,*,end=666) x1
         n=n+1
         x(n)=x1-x0
         if(x1.le.x0) xn0=float(n)
         y(n)=float(n)*fac
      enddo
 666  continue
      close(1)
      do i=1,n
         y(i)=y(i)-xn0
      enddo
      call pgsci(2)
      do i=1,n
         x(i)=x(i)/365.+2018
c         y(i)=y(i)+90.
c         y(i)=y(i)+230.
      enddo
      call pgline(n,x,y)
      do i=1050,n
         y(i)=y(i)+90.
      enddo
      call pgsls(2)
c      call pgline(n,x,y)
      call getct(0.,xt1,npred,xpred,ypred,n,x,y,0)
      call getct(xt1,xt2,npred,xpred,ypred,n,x,y,1)
      call getct(xt2,xt3,npred,xpred,ypred,n,x,y,1)
      call getct(xt3,xt4,npred,xpred,ypred,n,x,y,1)
      call getct(xt4,xt5,npred,xpred,ypred,n,x,y,1)
      call getct(xt5,xt6,npred,xpred,ypred,n,x,y,1)
      call getct(xt6,xt7,npred,xpred,ypred,n,x,y,1)
      call getct(xt7,xt8,npred,xpred,ypred,n,x,y,1)
      call getct(xt8,xt9,npred,xpred,ypred,n,x,y,1)
      call getct(xt9,xt10,npred,xpred,ypred,n,x,y,1)
      call getct(xt10,xt11,npred,xpred,ypred,n,x,y,0)

      call pgsci(1)
      call pgsls(4)
      yl(1)=ymin
      yl(2)=ymax
      xl(1)=xt1
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt2
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt3
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt4
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt5
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt6
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt7
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt8
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt9
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt10
      xl(2)=xl(1)
      call pgline(2,xl,yl)
      xl(1)=xt11
      xl(2)=xl(1)
      call pgline(2,xl,yl)

      call pgsls(1)
      call pgsci(1)
      call pgsch(0.8)
      ic=0
      do i=2,11
         if(ic.eq.0) xp=100.
         if(ic.eq.1) xp=250.
         ic=ic+1
         call pgptxt(xtril(i),xp,0.0,0.5,ctri(i))
         if(ic.eq.2) ic=0
      enddo

      call pgend

      end

      subroutine getct(x1,x2,npred,xpred,ypred,n,x,y,iplot)
      real xpred(npred),ypred(npred),x(n),y(n)
      character clab*5

      call xlinint(x1,npred,xpred,ypred,y1p)
      call xlinint(x2,npred,xpred,ypred,y2p)
      call xlinint(x1,n,x,y,y1)
      call xlinint(x2,n,x,y,y2)
      if(x1.eq.0) y1=0
      if(x1.eq.0) y1p=0
      
      ysp=y2p-y1p
      ys=y2-y1
      r1=(ys-ysp)/ysp*100.
      r2=(y2-y2p)/y2p*100.
      xpos=(x2+x1)/2.
      ypos=4000.
      call pgsch(0.7)
      write(clab,1001) r1
      if(iplot.eq.1) call pgptxt(xpos,ypos,0.,0.5,clab)
      write(clab,1001) r2
      ypos=ypos-150.
      if(iplot.eq.1) call pgptxt(xpos,ypos,0.,0.5,clab)
      call pgsch(1.0)
 1001 format(f5.1)
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
      if(xp.lt.x(1)) yp=y(1)
      if(xp.gt.x(n)) yp=y(n)
      return
      end
