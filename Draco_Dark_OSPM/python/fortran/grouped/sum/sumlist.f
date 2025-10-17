
	parameter(nmax=100)
	real b(nmax)
	real*8 x2,x3,x4
	real*8 s(nmax),s2(nmax),s3(nmax)
	integer*8 i2,is(nmax)
	character c1*20

	do i=1,nmax
	   is(i)=0
	   s(i)=0.d0
	   s2(i)=0.d0
	   s3(i)=0.d0
	enddo
	
	open(unit=1,file='list',status='old')
	do i=1,10000
	   read(1,*,end=666) c1
	   open(unit=2,file=c1,status='old')
	   nt=0
	   do j=1,nmax
c	      read(2,*,end=667) x1,x2
	      read(2,*,end=667) x1,x2,x3,x4
	      nt=nt+1
c	      is(j)=is(j)+i2
	      s(j)=s(j)+x2
	      s2(j)=s2(j)+x3
	      s3(j)=s3(j)+x4
	      b(j)=x1
	   enddo
 667	   continue
	   close(2)
	enddo
 666	continue
	close(1)

	do i=1,nt
	   print *,b(i),sngl(s(i)),sngl(s2(i)),sngl(s3(i))
c	   print *,b(i),sngl(s(i))
	enddo

	end
