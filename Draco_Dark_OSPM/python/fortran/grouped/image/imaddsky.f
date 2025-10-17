
      parameter (narrm1=1036,narrm2=1344,nfmax=312)
      real xd(narrm1,narrm2),vala(1036),val(1036,1344)
      real specs(1036,nfmax*112),xd1(narrm1,narrm2),xd2(narrm1,narrm2)
      integer naxes(2),iamp(4)
      character n1*11,n2*12,file1*120,file2*120,a1*14
      character cexp(3)*5
      logical simple,extend,anyf

      ncol=1036
      nrow=112
      iread=0

      do ir1=1,1344
         do icol=1,1036
            val(icol,ir1)=0.
         enddo
      enddo
      
      read *,n1,n2

      file1="in.fits"
c      file1="/data/00115/gebhardt/calfits/"//n2(1:6)//"/"//
c     $     n2//"_"//n1//"_cal.fits"

      cexp(1)="exp01"
      cexp(2)="exp02"
      cexp(3)="exp03"
      ir1=0
      iamp(1)=0
      iamp(2)=0
      iamp(3)=0
      iamp(4)=0
      do j=1,3
c         file2="/scratch/00115/gebhardt/alldet/output/d"
         file2="/scratch/03261/polonius/science_reductions"//
     $        "/alldet/output/d"
     $        //n2(1:8)//
     $        "s"//n2(10:12)//cexp(j)//"amp.dat"
         open(unit=1,file=file2,status='old',err=666)
         read(1,*)
         do i=1,10000
            read(1,*,end=666) a1
            if(a1(1:11).eq.n1) then
               if(a1(13:14).eq."LL") iamp(1)=(i-1)*112+1
               if(a1(13:14).eq."LU") iamp(2)=(i-1)*112+1
               if(a1(13:14).eq."RL") iamp(3)=(i-1)*112+1
               if(a1(13:14).eq."RU") iamp(4)=(i-1)*112+1
c               istart=i
c               goto 666
            endif
         enddo
 666     continue
         close(1)
c         istart=(istart-1)*112+1

c         file2="/scratch/00115/gebhardt/alldet/output/d"//n2(1:8)
         file2="/scratch/03261/polonius/science_reductions/"//
     $        "/alldet/output/d"//n2(1:8)
     $        //"s"//n2(10:12)//cexp(j)//"sub.fits"
         im1=0
         ier=0
         call ftgiou(im1,ier)
         call ftopen(im1,file2,iread,iblock,ier)
         call ftghpr(im1,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         ncol=naxes(1)
         nrow=naxes(2)
         call ftg2de(im1,igc,0.,1036,ncol,nrow,specs,anyf,ier)
         call ftclos(im1,ier)
         do ia=1,4
            if(iamp(ia).eq.0) then
               ir1=ir1+112
            else
               istart=iamp(ia)
               iend=istart+111
               do irow=istart,iend
c               do irow=istart,istart+447
                  ir1=ir1+1
                  do icol=1,1036
                     val(icol,ir1)=specs(icol,irow)
                  enddo
               enddo
            endif
         enddo
      enddo

      ier=0
      call ftopen(51,file1,iread,iblock,ier)
      call ftinit(50,'out.fits',iblock,ier)
      if(ier.ne.0) then
         call ftclos(51,ier)
         write(*,*) 'Error opening image : ',file1
         goto 706
      endif
      do iext=1,4
         call ftmahd(51,iext,ihd,ier)
         call ftghpr(51,2,simple,ibit,naxis,naxes,ipc,igc,extend,ier)
         if(naxis.eq.1) naxes(2)=1
         if(naxes(1).gt.narrm1.or.naxes(2).gt.narrm2) then
            print *,iext,naxes(1),naxes(2)
            write(*,"('Arrays too small - make narrm bigger')")
            goto 706
         endif
         ncol=naxes(1)
         nrow=naxes(2)
         if(iext.eq.1) nrow1=nrow
         call ftg2de(51,igc,0.,narrm1,ncol,nrow,xd,anyf,ier)
         if(iext.eq.1) then
            do j=1,nrow
               do i=1,ncol
                  xd1(i,j)=xd(i,j)
               enddo
            enddo
         endif
         if(iext.eq.3) then
            do j=1,nrow
               do i=1,ncol
                  xd2(i,j)=xd(i,j)
               enddo
            enddo
         endif
         call ftcopy(51,50,0,ier)
         call ftp2de(50,igc,narrm1,naxes(1),naxes(2),xd,ier)
      enddo
      ncol=1036
      nrow=1344
      naxes(1)=ncol
      naxes(2)=nrow
      if(nrow1.ne.nrow) then
         print *,"Missing Amps in: ",n2
      endif

      nzero=0
      do j=1,nrow
         do i=1,ncol
            if(xd1(i,j).ne.0.and.xd2(i,j).ne.0) then
               factor=xd1(i,j)/xd2(i,j)
            else
               factor=0.
            endif
            xtmp=factor*val(i,j)
            val(i,j)=xtmp
            if(val(i,j).eq.0) nzero=nzero+1
         enddo
      enddo

      print *,"Nzero= ",nzero,n2," ",n1
      
      call ftiimg(50,-32,2,naxes,ier)
      call ftp2de(50,igc,1036,naxes(1),naxes(2),val,ier)
      call ftpkls(50,'EXTNAME','fullsky',"Label",ier)

      if(ier.ne.0) call ftdelt(50,ier)

      call ftclos(51,ier)
      call ftclos(50,ier)

 706  continue
      end
