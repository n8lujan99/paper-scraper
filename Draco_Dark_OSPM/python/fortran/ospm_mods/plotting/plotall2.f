
      parameter(nmax=10000)
      real hole(nmax),xmu(nmax),alp(nmax),ent(nmax)
      real chi(nmax),xl(2),yl(2),xlbad(nmax),holen(nmax),chin(nmax)
      real as(nmax),cs(nmax)
      integer nname(nmax),idi(nmax),idh(nmax)
      character fdum*40,filen(nmax)*40,lab(nmax)*14,gname*40
      
      parameter(pi=3.1415926539)

      data big/1.e20/

      open(unit=1,file='plotall.param',status='old')
      read(1,*) gname
      read(1,*) cmin,cmax,xmin,xmax,bh0,ilog
      close(1)
      open(unit=11,file='plotall.out',status='unknown')
      open(unit=1,file='plotall.in',status='old')

      open(unit=2,file='plotiter.out',status='old')
      ns=0
      do i=1,nmax
         read(2,*,end=868) x1,x2
         ns=ns+1
         as(ns)=x1
         cs(ns)=x2
      enddo
 868  continue
      close(2)

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
         read(1,*) x1b
         read(1,*) x2,x3,x4
         read(1,*) xbh,xincp,xd,chival,xd,xd6
         close(1)
         xmu(i)=x1b
         if(xbh.eq.0) then
c            hole(i)=log10(1.e6)
            hole(i)=log10(bh0)
         else
            hole(i)=log10(xmu(i)*xbh)
         endif
         if(ilog.eq.0) hole(i)=xmu(i)*xbh
         alp(i)=x2
         if(alp(i).lt.1) print *,fdum(1:nname(i)),alp(i),chival,
     $        '  LOW ALPHA'
         call xlinint(alp(i),ns,as,cs,cdiff)
         chival=chival-cdiff
         ent(i)=x3+x4
c         chi(i)=x5
         chi(i)=chival
         xlbad(i)=xd6
         write(11,*) hole(i),xmu(i),chi(i),cdiff,
     $        " "//fdum(1:nname(iall))
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
         print *,fdum(1:nname(i)),hole(i),chi(i)
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
