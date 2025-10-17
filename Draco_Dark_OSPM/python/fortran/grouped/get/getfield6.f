Skipped
      parameter(nmax=100000000,nmax2=10000)
      integer*8 i6,i7
      real rfcen(nmax2),dfcen(nmax2)
      character cdate*8,cshot*3,cname*12,cfcen(nmax2)*12
      character file1*60,a12*11,cf2*12

      rad1=660.
      rad2=2.

      file1="/scratch/projects/hetdex/detect/dithall/"

      open(unit=1,file="justdex",status='old')
      nf=0
      do i=1,nmax2
         read(1,*,end=668) cdate,cshot,cname,x4,x5,x6,x7,x8
         if(x5>0) then
            nf=nf+1
            cfcen(nf)=cdate//"v"//cshot
            rfcen(nf)=x5*15.
            dfcen(nf)=x6
         endif
      enddo
 668  continue
      close(1)

      open(unit=1,file='in',status='old')
      open(unit=11,file='out1',status='unknown')
      open(unit=12,file='out2',status='unknown')
      
      do iall=1,10000000
         read(1,*,end=667) x1,x2,x3,x4,x5,i6,i7,x8,x9,x10,x11,a12
         ra=x9
         dec=x10
         cosd=cos(dec/57.3)
         nt=0
         do i=1,nf
            rad=sqrt( (cosd*(rfcen(i)-ra))**2 + (dfcen(i)-dec)**2 )
            rad=rad*3600.
            if(rad.lt.rad1) then
               cf2=a12(1:8)//"v"//a12(9:11)
               if(cf2.ne.cfcen(i)) then
                  file1="/scratch/projects/hetdex/detect/dithall/"
     $                 //cfcen(i)//".dithall"
                  open(unit=2,file=file1,status='old',err=445)
                  do j=1,110000
                     read(2,*,end=445) y1,y2
                     rad=sqrt((cosd*(ra-y1))**2+(dec-y2)**2)
                     rad=rad*3600.
                     if(rad.lt.rad2) goto 444
                  enddo
                  goto 445
 444              continue
                  nt=nt+1
                  write(11,*) i7,ra,dec,cf2," ",cfcen(i)
 445              continue
                  close(2)
               endif
            endif
         enddo
         write(12,*) i7,nt+1
      enddo
 667  continue
      close(1)
      close(11)
      close(12)

      end
