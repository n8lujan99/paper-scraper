
      parameter(nfmax=1000)
      real f(1000),c(nfmax,1000),c1(1000)
      character f1*100
      
      open(unit=1,file="l1",status="old")

      call pgbegin(0,'?',1,1)
      call pgpap(0.,1.)
      call pgscf(2)
      call pgsch(1.5)
      call pgslw(2)

      open(unit=11,file="out",status="unknown")
      do iflux=1,7
         call pgslw(2)
         call pgenv(0.,30.,0.,0.25,0,0)
         call pgslw(1)
         ic=1

         n=0
         do i=1,nfmax
            read(1,*,end=666) f1
            ibad=0
            n=n+1
            open(unit=2,file=f1,status="old")
            if(iflux.gt.1) then
               do j=1,iflux-1
                  read(2,*) i1,x2,x3,x4,i5,x6
                  do k=1,i5
                     read(2,*)
                  enddo
               enddo
            endif               
            do j=1,1
               read(2,*,end=667) i1,x2,x3,x4,i5,x6
               if(i5.ne.79) then
c                  print *,f1
                  n=n-1
                  ibad=1
                  goto 668
               endif
               nf=i5
               i1o=i1
               x2o=x2
               x3o=x3
               x4o=x4
               i5o=i5
               x6o=x6
               do k=1,i5
                  read(2,*) x1,x2,x3
                  f(k)=x1
                  c(n,k)=x3
                  c1(k)=x3
               enddo
            enddo
 667        continue
c            print *,nf
            ic=ic+1
            if(ic.eq.14) ic=1
            call pgsci(ic)
            call pgline(nf,f,c1)
 668        close(2)
         enddo
 666     continue
         rewind(1)

c         print *,nf
         sumt=0.
         do i=1,nf
            sum=0.
            do j=1,n
               sum=sum+c(j,i)
            enddo
            sum=sum/float(n)
            c1(i)=sum
            sumt=sumt+c1(i)
         enddo
c         print *,sumt
         write(11,*) 666,666.,666.,x4o,i5o,x6o
         do i=1,nf
            c1(i)=c1(i)/sumt
            write(11,*) f(i),c1(i),c1(i)
         enddo

         call pgsci(1)
         call pgslw(13)
         call pgline(nf,f,c1)
         call pgsci(1)
      enddo
      close(1)
      close(11)

      call pgend

      end
