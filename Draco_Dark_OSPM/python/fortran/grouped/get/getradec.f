Skipped
      character cfile*6,file1*60,file2*60

      read *,ra1,ra2,dec1,dec2

      file1="/work/00115/gebhardt/findsource/imaging/"
      file2=file1(1:40)//"poslist"

      open(unit=1,file=file2,status='old')
      open(unit=11,file='out',status='unknown')
      do i=1,1000
         read(1,*,end=666) cfile,x2,x3
         if((ra1.ge.x2.and.ra1.le.x3).or.
     $      (ra2.ge.x2.and.ra2.le.x3).or.
     $      (x2.ge.ra1.and.x3.le.ra2)) then
            file2=file1(1:40)//"output/"//cfile
            open(unit=2,file=file2,status='old')
            print *,cfile
            do j=1,1000000
               read(2,*,end=667) ra,dec
               if(ra.ge.ra1.and.ra.le.ra2.and.
     $              dec.ge.dec1.and.dec.le.dec2) then
                  write(11,*) ra,dec
               endif
            enddo
 667        continue
            close(2)
         endif
      enddo
 666  continue
      close(1)
      close(11)

      end
