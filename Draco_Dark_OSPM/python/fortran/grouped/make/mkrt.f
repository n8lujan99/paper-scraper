
c      print *,"Enter N: "
      read *,ntot
      istart=1
      do i=1,999
         i1=1+(i-1)*10000
         i2=i*10000
         iout=100+i
         print *,"rgetcorr",i1,i2,iout
         if(i2.gt.ntot) goto 666
      enddo
 666  continue

      end
