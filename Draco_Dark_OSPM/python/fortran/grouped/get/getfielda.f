Marked
      parameter(nmax=100000000,nmax2=5000)
      real ra(nmax),dec(nmax),rfcen(nmax2),dfcen(nmax2)
      integer istart(nmax2),isa(nmax2),iend(nmax2),iea(nmax2)
      character a3*6,a8*28,a9*17,a10*5,a11*8,a12*3
      character cspec(nmax)*18,cfield(nmax)*12,cfcen(nmax2)*12
      character cold*12,cnew*12,cdetect*12,crep(100)*12

      rad1=600.
      radd=3.

      open(unit=1,file="/work/00115/gebhardt/hdr/data/radec.dat",
     $     status='old')
      nf=0
      do i=1,nmax2
         read(1,*,end=668) cold,x2,x3
         nf=nf+1
         cfcen(nf)=cold
         rfcen(nf)=x2
         dfcen(nf)=x3
      enddo
 668  continue
      close(1)

      open(unit=1,file="in2",status='old')
 
      cold='666'
      n=0
      do i=1,nmax
c      do i=1,5000000
         read(1,*,end=666) x1,x2,a3,x4,x5,x6,x7,a8,a9,
     $        a10,a11,a12
         n=n+1
         ra(n)=x1
         dec(n)=x2
         cspec(n)=a8(7:24)
         cfield(n)=a11//"v"//a12
         if(cold.ne.cfield(n)) then
            do j=1,nf
               if(cfield(n).eq.cfcen(j)) then
                  istart(j)=n
                  goto 888
               endif
            enddo
            print *,"Not here ",cfield(n)
 888        continue
            cold=cfield(n)
         else
            iend(j)=n
         endif
         if(i/1000000.eq.float(i)/1000000.) then
            write(6,"('Reading ',i10,a1,$)") i,char(13)
            call flush(6)
         endif
      enddo
 666  continue
      close(1)

      print *

      open(unit=1,file="in",status='old')
      open(unit=11,file="out",status='unknown')
      do j=1,1000000
         read(1,*,end=667) x1,x2,x3,i4,a8,cdetect,x7
         if(j/100.eq.float(j)/100.) then
            write(6,"('Doing ',i7,a1,$)") j,char(13)
            call flush(6)
         endif
         cosd=cos(x2/57.3)
         nt=0
         do i=1,nf
            rad=sqrt( (cosd*(rfcen(i)-x1))**2 + (dfcen(i)-x2)**2 )
            rad=rad*3600.
            if(rad.lt.rad1) then
               nt=nt+1
               isa(nt)=istart(i)
               iea(nt)=iend(i)
            endif
         enddo
         cold="new"
         nrep=0
         do k=1,nt
            i1=isa(k)
            i2=iea(k)
            do i=i1,i2
c         do i=1,n
c               if(cfield(i).ne.cold) then
               rad=sqrt( (cosd*(ra(i)-x1))**2 + (dec(i)-x2)**2 )
               rad=rad*3600.
               if(rad.lt.radd) then
                  cold=cfield(i)
                  nrep=nrep+1
                  crep(nrep)=cfield(i)
                  goto 777
               endif
c              endif
            enddo
 777        continue
         enddo
         if(nrep.gt.1) then
            do i=1,nrep
               write(11,1101) j,nrep,x1,x2,x3,i4,a8,cdetect,crep(i),x7
            enddo
         endif
 999     continue
      enddo
 667  continue
      close(1)
      close(11)

 1101 format(i7,1x,i4,2(1x,f10.5),1x,f7.2,1x,i8,1x,a25,1x,a12,1x,a12,
     $     1x,f7.2)

      end
