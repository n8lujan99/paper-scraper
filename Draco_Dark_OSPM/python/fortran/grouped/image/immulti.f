
      parameter (narrm1=5000,narrm2=5000)
      real xd(narrm1,narrm2),xderr(narrm1,narrm2)
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
      file1="ctrap.all"
      open(unit=1,file=file1,status='old')
      nm=0
      do i=1,10000
         read(1,*,end=666) i1,i2,i3,x4,i5,a6
         if(a6.eq.cname) then
            nm=nm+1
            ix(nm)=i1
            iyl(nm)=i2
            iyh(nm)=i3
         endif
      enddo
 666  continue
      close(1)
c - these are the arrays to mask:
c    2 the 2d error 
c    3 the 2d sky subtracted
c   11 extracted spectrum
c   15 extracted error 
c   16 extracted sky subtracted

c - first get 2d locations from trace
      call ftmahd(51,13,ihd,ier)
      call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
      ncol=naxes(1)
      nrow=naxes(2)
      call ftg2de(51,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
      do i=1,nm
         iy2dl(i)=nint(xd(ix(i),iyl(i)))
         iy2dh(i)=nint(xd(ix(i),iyh(i)))
         if(iyl(i).lt.3) iy2dl(i)=1
         if(iyh(i).gt.110) iy2dh(i)=1032
      enddo

c - now replace      
      do iext=1,17
         ier=0
         call ftmahd(51,iext,ihd,ier)
         call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=naxes(1)
         nrow=naxes(2)
         call ftg2de(51,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
         call ftcopy(51,50,0,ier)
         if(nm.gt.0) then
            if(iext.eq.2.or.iext.eq.3) then
               do i=1,nm
                  do j=iy2dl(i),iy2dh(i)
                     xd(ix(i),j)=0.
                  enddo
               enddo
            endif
            if(iext.eq.11.or.iext.eq.15.or.iext.eq.16) then
               do i=1,nm
                  do j=iyl(i),iyh(i)
                     xd(ix(i),j)=0.
                  enddo
               enddo
            endif
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
