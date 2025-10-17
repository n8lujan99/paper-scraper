Skipped
      parameter(nmax=100000000,nmax2=5000)
      real ra(nmax),dec(nmax),rfcen(nmax2),dfcen(nmax2)
      integer istart(nmax2),isa(nmax2),iend(nmax2),iea(nmax2)
      character a3*6,a8*28,a9*17,a10*5,a11*8,a12*3
      character cspec(nmax)*18,cfield(nmax)*12,cfcen(nmax2)*12
      character cold*12,cnew*12,cdetect*12,crep(100)*12,a1*20

      rad1=600.
      radd=3.

      open(unit=1,file="/work/00115/gebhardt/hdr/data/radec.all",
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

c      open(unit=1,file='sob.txt',status='old')
      open(unit=1,file='odin1.txt',status='old')
      read(1,*)
      do iall=1,10000
c         read(1,*,end=667) a1,x2,x3,x4
         read(1,*,end=667) i1,x2,x3
         cosd=cos(x3/57.3)
         nt=0
         do i=1,nf
            rad=sqrt( (cosd*(rfcen(i)-x2))**2 + (dfcen(i)-x3)**2 )
            rad=rad*3600.
            if(rad.lt.rad1) then
               nt=nt+1
c               print *,a1,cfcen(i),x2,x3,x4
               print *,i1,cfcen(i),x2,x3
            endif
         enddo
      enddo
 667  continue
      close(1)
      close(11)

 1101 format(i7,1x,i4,2(1x,f10.5),1x,f7.2,1x,i8,1x,a25,1x,a12,1x,a12,
     $     1x,f7.2)

      end
