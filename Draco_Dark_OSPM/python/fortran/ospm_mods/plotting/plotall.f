
      parameter(nmax=10000)
      real hole(nmax),xmu(nmax),alp(nmax),ent(nmax)
      real chi(nmax),xl(2),yl(2),xlbad(nmax),holen(nmax),chin(nmax)
      integer nname(nmax),idi(nmax),idh(nmax)
      character fdum*40,filen(nmax)*40,lab(nmax)*14,gname*40
      
      parameter(pi=3.1415926539)

      data big/1.e20/

      open(unit=1,file='plotall.param',status='old')
      read(1,*) gname
      read(1,*) cmin,cmax,xmin,xmax,bh0,ilog
      close(1)
c      read(1,*) fdum
c      inset=0
c      if(fdum(1:1).eq.'0') then
c         rewind(1)
c         read(1,*)
c         read(1,*) 
c         read(1,*) x1,bhmin2,bhmax2,cmin2,cmax2
c         print *,bhmin2,bhmax2,cmin2,cmax2
c         inset=1
c      else
c         rewind(1)
c         read(1,*)
c         read(1,*)      
c      endif
      open(unit=11,file='plotall.out',status='unknown')
      open(unit=1,file='plotall.in',status='old')

      n=0
      do i=1,nmax
         read(1,*,end=666) fdum
         n=n+1
         nn=0
         do ii=1,70
            if(fdum(ii:ii).eq.' ') then
               goto 867
            else
               nn=nn+1
            endif
         enddo
 867     continue
         filen(n)=fdum
         nname(n)=nn
      enddo
 666  continue
      close(1)

      bhmax=-big
      bhmin=big
      i=0
      do iall=1,n
         fdum=filen(iall)
         open(unit=1,file=fdum(1:nname(iall)),status='old',
     $        err=767)
         i=i+1
         filen(i)=fdum
         read(1,*) x1b
         read(1,*) x2,x3,x4
         read(1,*) xbh,xincp,xd,chival,xd,xd6
         close(1)
         xmu(i)=x1b
         if(xbh.eq.0) then
            hole(i)=log10(bh0)
         else
            hole(i)=log10(xmu(i)*xbh)
         endif
         if(ilog.eq.0) hole(i)=xmu(i)*xbh
         alp(i)=x2
         if(alp(i).lt.1) print *,fdum(1:nname(i)),alp(i),chival,
     $        '  LOW ALPHA'
         ent(i)=x3+x4
         chi(i)=chival
         xlbad(i)=xd6
         write(11,*) hole(i),xmu(i),chi(i)," "//fdum(1:nname(iall))
c         print *,fdum(1:nname(iall)),hole(i),chi(i),xlbad(i)
         bhmin=min(bhmin,hole(i))
         bhmax=max(bhmax,hole(i))
 767     close(1)
      enddo
      close(11)
      n=i

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgsch(1.5)
      call pgscf(2)

      if(ilog.eq.1) then
         call pgenv(bhmin-.3,bhmax+.4,cmin,cmax,0,10)
      else
         call pgenv(xmin,xmax,cmin,cmax,0,0)
         call pgsfs(2)
         if(inset.eq.1) call pgrect(bhmin2,bhmax2,cmin2,cmax2)
      endif
      ltype=0
      lthick=1
      np=0
      do i=1,n
         np=np+1
         holen(np)=hole(i)
         chin(np)=chi(i)
         fdum=filen(i)
         print *,fdum(1:nname(i)),hole(i),chi(i),xlbad(i)
      enddo
      call sort2(np,holen,chin)
      ltype=ltype+1
      if(ltype.eq.4) then
         ltype=0
         lthick=lthick+2
      endif
      call pgslw(lthick)
      call pgline(np,holen,chin)
      if(ilog.eq.1) then
         call pgpoint(np-1,holen(2),chin(2),17)
         call pgsch(2.2)
         call pgslw(3)
         call pgpoint(1,holen(1),chin(1),28)
      else
         call pgpoint(np,holen,chin,17)
      endif
      call pgsls(1)
      call pgslw(1)
      call pgsch(1.5)
      call pgmtxt('B',2.5,0.5,0.5,
     $     'M\\D\\(0850)\\U (M\\D\\(2281)\\U)')
      call pglabel('','\\Gx\\U2\\D','')
      call pgmtxt('T',-1.7,0.5,0.5,gname)

      if(inset.eq.1) then
         call pgsch(0.9)
         call pgvport(0.60,0.82,0.60,0.72)
         call pgswin(bhmin2,bhmax2,cmin2,cmax2)
         call pgbox('bcnst',0.,0,'bcs',0.,0)
         ltype=0
         lthick=1
         np=0
         do i=1,n
            np=np+1
            holen(np)=10**hole(i)
            chin(np)=chi(i)
         enddo
         call sort2(np,holen,chin)
         ltype=ltype+1
         if(ltype.eq.4) then
            ltype=0
            lthick=lthick+2
         endif
         call pgslw(lthick)
         call pgline(np,holen,chin)
         call pgpoint(np,holen,chin,17)
         call pgsls(1)
         call pgslw(1)
      endif

      call pgend

      end
      FUNCTION gammq(a,x)
      REAL a,gammq,x
CU    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
      SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      REAL arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
