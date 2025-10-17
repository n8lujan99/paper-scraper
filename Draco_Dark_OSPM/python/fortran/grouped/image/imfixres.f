
      parameter (narrm1=5000,narrm2=5000)
      real xd(narrm1,narrm2),xderr(narrm1,narrm2),resmap(1032,112)
      integer naxes(2),ix(narrm1),iyl(narrm1),iyh(narrm1)
      integer iy2dl(narrm1),iy2dh(narrm1)
      character file1*120,kname*120,kname0*160
      character cpos*1,chalf*1,cspecid*3,cname*5,a6*5
      logical simple,extend,anyf

      read *,file1

      im1=0
      ier=0
c      call ftgiou(51,ier)
      iread=0
      call ftopen(51,file1,iread,iblock,ier)
      call ftinit(50,'out.fits',iblock,ier)
      if(ier.ne.0) then
         call ftclos(51,ier)
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      call ftgkye(51,'SPECID',specid,file1,ier)
      ispec=nint(specid)
      cspecid="000"
      if(ispec.lt.10) write(cspecid(3:3),1001) ispec
      if(ispec.ge.10.and.ispec.lt.100) write(cspecid(2:3),1002) ispec
      if(ispec.ge.100) write(cspecid(1:3),1003) ispec
      call ftgkys(51,'CCDPOS',cpos,file1,ier)
      call ftgkys(51,'CCDHALF',chalf,file1,ier)
      cname=cspecid//cpos//chalf

c- get the res file
      file1="/scratch/03946/hetdex/lib_calib/reschi/res"//cname//".fits"
      call ftopen(52,file1,iread,iblock,ier)
      if(ier.eq.0) then
         iext=1
         call ftmahd(52,iext,ihd,ier)
         call ftghpr(52,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=naxes(1)
         nrow=max(1,naxes(2))
         call ftg2de(52,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
         call ftclos(52,ier)
         do i=1,ncol
            do j=1,nrow
               resmap(i,j)=xd(i,j)
            enddo
         enddo
      else
         ier=0
         do j=1,112
            do i=1,1032
               resmap(i,j)=0.
            enddo
         enddo
      endif
      do j=1,112
         do i=1,1032
            if(resmap(i,j).lt.-0.3) resmap(i,j)=0.
         enddo
      enddo

      call ftmahd(51,17,ihd,ier)
      call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(51,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
      do i=1,ncol
         do j=1,nrow
            resmap(i,j)=resmap(i,j)*xd(i,j)
         enddo
      enddo

c - these are the arrays to mask:
c    2 the 2d error 
c    3 the 2d sky subtracted
c   11 extracted spectrum
c   15 extracted error 
c   16 extracted sky subtracted

c - now replace      
      do iext=1,17
         ier=0
         call ftmahd(51,iext,ihd,ier)
         call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=naxes(1)
         nrow=naxes(2)
         call ftg2de(51,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
         call ftcopy(51,50,0,ier)
         if(iext.eq.16) then
            do i=1,ncol
               do j=1,nrow
                  xd0=xd(i,j)
                  xd(i,j)=xd(i,j)-resmap(i,j)
               enddo
            enddo
         endif
         call ftp2de(50,igc,narrm1,naxes(1),naxes(2),xd,ier)
      enddo
      call ftclos(51,ier)
      call ftclos(50,ier)

 706  continue
 1001 format(i1)
 1002 format(i2)
 1003 format(i3)
      end
