
      parameter (narrm1=2000,narrm2=100000)
      character file1*180,ttype(3)*10,class*6
c      character nullstr*1,name*8,cname(narrm2)*18
      character nullstr*1,name*8,cname(narrm2)*24,cspec*15
      logical simple,extend,anyf

      file1="specObj-dr12.fits"

      im1=0
      ier=0
      iread=0
      call ftgiou(im1,ier)
      call ftopen(im1,file1,iread,iblock,ier)

      call ftmahd(im1,2,ihd,ier)
c      call ftgkns(im1,'TTYPE',1,3,ttype,nfound,ier)
      open(unit=11,file='out',status='unknown')
      do i=1,100000000
         call ftgcvj(im1,57,i,1,1,0.,iplate,anyf,ier)
         call ftgcvj(im1,58,i,1,1,0.,itile,anyf,ier)
         call ftgcvj(im1,59,i,1,1,0.,imjd,anyf,ier)
         call ftgcvj(im1,60,i,1,1,0.,ifiberid,anyf,ier)
         call ftgcvj(im1,61,i,1,1,0.,ibjid,anyf,ier)
         call ftgcve(im1,62,i,1,1,0.,ra,anyf,ier)
         call ftgcve(im1,63,i,1,1,0.,dec,anyf,ier)
         call ftgcve(im1,66,i,1,1,0.,z,anyf,ier)
c         print *,i,ra,dec,ier
         call ftgcvs(im1,64,i,1,1,nullstr,class,anyf,ier)
         cspec="0000 00000 0000"
         if(iplate.lt.10) write(cspec(4:4),1001) iplate
         if(iplate.ge.10.and.iplate.lt.100) 
     $        write(cspec(3:4),1002) iplate
         if(iplate.ge.100.and.iplate.lt.1000) 
     $        write(cspec(2:4),1003) iplate
         if(iplate.ge.1000.and.iplate.lt.10000) 
     $        write(cspec(1:4),1004) iplate
         if(imjd.lt.10) write(cspec(10:10),1001) imjd
         if(imjd.ge.10.and.imjd.lt.100) 
     $        write(cspec(9:10),1002) imjd
         if(imjd.ge.100.and.imjd.lt.1000) 
     $        write(cspec(8:10),1003) imjd
         if(imjd.ge.1000.and.imjd.lt.10000) 
     $        write(cspec(7:10),1004) imjd
         if(imjd.ge.10000.and.imjd.lt.100000) 
     $        write(cspec(6:10),1005) imjd
         if(ifiberid.lt.10) write(cspec(15:15),1001) ifiberid
         if(ifiberid.ge.10.and.ifiberid.lt.100) 
     $        write(cspec(14:15),1002) ifiberid
         if(ifiberid.ge.100.and.ifiberid.lt.1000) 
     $        write(cspec(13:15),1003) ifiberid
         if(ifiberid.ge.1000.and.ifiberid.lt.10000) 
     $        write(cspec(12:15),1004) ifiberid
         if(ifiberid.ge.10000.and.ifiberid.lt.100000) 
     $        write(cspec(11:15),1005) ifiberid
         if(class(1:4).eq.'STAR') write(11,1101) i,ra,dec,z,class,
     $        cspec(1:4),cspec
      enddo
      close(11)

      call ftclos(im1,ier)

 706  continue
 1001 format(i1)
 1002 format(i2)
 1003 format(i3)
 1004 format(i4)
 1005 format(i5)
 1101 format(i8,3(1x,f10.5),1x,a5,1x,a4,1x,a15)
      end
