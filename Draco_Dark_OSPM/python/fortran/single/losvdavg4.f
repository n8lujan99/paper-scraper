
      parameter(nmax=1000,nfm=100)
      real v(nmax),ya(nfm,nmax),yin(nfm),yl(nfm,nmax),yh(nfm,nmax)
      character file1*30,file2*30

      open(unit=1,file='listbin',status='old')
      do nf=1,nfm
         read(1,*,end=666) file1
         nft=nf
         open(unit=2,file=file1,status='old')
         n=0
         do i=1,nmax
            read(2,*,end=667) x1,x2,x3,x4
            n=n+1
            v(n)=x1
            ya(nf,n)=x2
            yl(nf,n)=x3
            yh(nf,n)=x4
         enddo
 667     continue
         close(2)
         nt=n
      enddo
 666  continue
      close(1)
      
      open(unit=11,file='losvd.out',status='unknown')
      do i=1,nt
         do j=1,nft
            yin(j)=ya(j,i)
         enddo
         call biwgt(yin,nft,xb,xs)
         do j=1,nft
            yin(j)=yl(j,i)
         enddo
         call biwgt(yin,nft,xbl,xsl)
         do j=1,nft
            yin(j)=yh(j,i)
         enddo
         call biwgt(yin,nft,xbh,xsh)
         write(11,*) v(i),xb,xbl,xbh
      enddo
      close(11)

      end


         
