
      parameter(nmax=10000)
      real d(nmax),rd(nmax),vel(nmax*10),vlo(nmax),vhi(nmax)
      real ybin(nmax),yg(nmax)
      character c0*10
      parameter(gee=4.3e-3)

      read *,xln,xhn,yln,yhn,xnormu

      open(unit=1,file='pfit.in',status='old')
      read(1,*) c0,xcen
      read(1,*) c0,ycen
      read(1,*) c0,angd
      read(1,*) c0,xinc
      read(1,*) c0,vsig
      read(1,*) c0,vsig_s
      read(1,*) c0,d(1)
      read(1,*) c0,d(2)
      read(1,*) c0,rd(1)
      read(1,*) c0,rd(2)
      close(1)

c      xcen=-0.03
c      ycen=0.03

      xln=xln-xcen
      xhn=xhn-xcen
      yln=yln-ycen
      yhn=yhn-ycen

      astopc=4.84*16.
      bh=6.5e9

c      angd=-20.
c      xinc=25.
      vcor=sin(xinc/57.29)
c      vsig=200.
c      vsig_s=500.

c      d(1)=1.
c      d(2)=0.1
c      rd(1)=0.02
c      rd(2)=0.15
      nd=2

      nvbin=29
      vmin=-1750.
      vmax=1750.
      vdiff=(vmax-vmin)/float(nvbin-1)
      do i=1,nvbin
         v=vmin+(vmax-vmin)*float(i-1)/float(nvbin-1)
         vlo(i)=v-vdiff/2.
         vhi(i)=v+vdiff/2.
         ybin(i)=0.
      enddo

      nstepx=100
      nstepy=100
      nstepr=100
      nv=0
      iseed=-1
c- integrate over the spatial element
      do ix=1,nstepx
         x=xln+(xhn-xln)*float(ix-1)/float(nstepx-1)
      do iy=1,nstepy
         y=yln+(yhn-yln)*float(iy-1)/float(nstepy-1)

         ang=57.29*atan(y/x)
         dang=abs(angd-ang)
         rlo=sqrt(x*x+y*y)
         vcor2=cos(dang/57.29)
         rhi=rd(nd)
         sum=0.
         sum2=0.
c- integrate along the line of sight
         do j=1,nstepr
            r0=rlo+(rhi-rlo)*float(j-1)/float(nstepr-1)
            call xlinint(r0,nd,rd,d,d0)
            if(r0.lt.rd(1)) d0=0.
            dran=ran2(iseed)
            if(dran.lt.d0) then
               rpc=r0*astopc
               v=vcor*vcor2*sqrt(gee*bh/rpc)
               v1=v+vsig*gasdev(iseed)
               v2=-(v+vsig*gasdev(iseed))
               sum=sum+x/r0*v*d0
               sum2=sum2+d0
               nv=nv+1
               vel(nv)=v1
               nv=nv+1
               vel(nv)=v2
               do iv=1,nvbin
                  if(v1.ge.vlo(iv).and.v1.lt.vhi(iv)) then
                     ybin(iv)=ybin(iv)+1
                     goto 666
                  endif
               enddo
 666           continue
               do iv=1,nvbin
                  if(v2.ge.vlo(iv).and.v2.lt.vhi(iv)) then
                     ybin(iv)=ybin(iv)+1
                     goto 667
                  endif
               enddo
 667           continue
            endif
         enddo
      enddo
      enddo

      sum=0
      sumg=0.
      do i=1,nvbin
         sum=sum+ybin(i)
         w=(vlo(i)+vhi(i))/2./vsig_s
         yg(i)=exp(-w*w/2.)
         sumg=sumg+yg(i)
      enddo
      open(unit=11,file='out',status='unknown')
      xnorm=sum/float(nstepx*nstepy*nstepr)
      xnorm=xnormu
      do i=1,nvbin
         y1=ybin(i)/sum
         y2=yg(i)/sumg
         yout=y1*xnorm+y2*(1.-xnorm)
         write(11,*) (vlo(i)+vhi(i))/2.,yout
      enddo
      close(11)

      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end
