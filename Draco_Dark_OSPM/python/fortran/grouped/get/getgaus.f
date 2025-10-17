Marked

      xfw=1.8
      sig=xfw/2.35
      xs=0.
      xe=sig*10.
      n=2000
      sum=0.
c      area=sqrt(2.*3.14159*sig**2)
      area=2.*3.14159*sig**2
      xdel=(xe-xs)/float(n)
c      area=area/xdel
      area=area/xdel/xdel
      xold=0.
      sum2=0.
      sum3=0.

      do j=1,n
         y=xs+float(j-1)*(xe-xs)/float(n-1)
         do i=1,n
            x=xs+float(i-1)*(xe-xs)/float(n-1)
            r=sqrt(x*x+y*y)
            x=r
            xp=x
            x=x/sig
            g=exp(-x*x/2.)/area
            sum=sum+g
            sum2=sum2+g*xp
            sum3=sum3+g*xdel*xdel
         enddo
      enddo
c      print *,2.*sum
c      print *,2.*sum2/sum
      print *,4.*sum
      avgr=sum2/sum
      print *,avgr,3.14159*avgr*avgr,4.*sum3
      end
