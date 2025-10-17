
      real xd(1032,112)
      integer naxes(2)
      character file1*120,cdate*8,cshot*3,cexp*5,cspec*24
      logical simple,extend,anyf

      open(unit=1,file='list',status='old')
      open(unit=11,file='out',status='unknown')

      do iall=1,1000000
         read(1,*,end=666) cdate,cshot,cexp,cspec,icol,id
         read(cspec(22:24),1001) irow
         file1="/data/00115/gebhardt/red1/reductions/"//cdate//
     $        "/virus/virus0000"//cshot//"/"//cexp//"/virus/"//
     $        cspec(1:20)//".fits"

         ier=0
         iread=0
         call ftopen(51,file1,iread,iblock,ier)
         iext=16
         call ftmahd(51,iext,ihd,ier)
         call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=naxes(1)
         nrow=naxes(2)
         call ftg2de(51,igc,0.,1032,ncol,nrow,xd,anyf,ier)
         call ftclos(51,ier)
         izero=0
         if(xd(icol,irow).eq.0) izero=1
         if(xd(icol-1,irow).eq.0) izero=1
         if(xd(icol+1,irow).eq.0) izero=1
         if(izero.eq.1) then 
            print *,cspec(7:20)," ",icol,irow,id
            if(id.lt.10) then
               write(11,1101) cdate,cshot,"_",cspec(7:17),"_",id
            elseif(id.ge.10.and.id.lt.100) then
               write(11,1102) cdate,cshot,"_",cspec(7:17),"_",id
            elseif(id.ge.100.and.id.lt.1000) then
               write(11,1103) cdate,cshot,"_",cspec(7:17),"_",id
            elseif(id.ge.1000.and.id.lt.10000) then
               write(11,1104) cdate,cshot,"_",cspec(7:17),"_",id
            endif
         endif
      enddo
 666  continue
      close(1)

 1001 format(i3)
 1101 format(a8,a3,a1,a11,a1,i1)
 1102 format(a8,a3,a1,a11,a1,i2)
 1103 format(a8,a3,a1,a11,a1,i3)
 1104 format(a8,a3,a1,a11,a1,i4)
      end
