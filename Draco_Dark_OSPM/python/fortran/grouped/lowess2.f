      subroutine lowess2(x, y, yerr, n, f, nsteps, delta, ys, rw, res)
      integer n
      integer nsteps
      real x(n), y(n), f, delta, ys(n), rw(n)
      real res(n), yerr(n)
      integer nright, min0, max0, i, j, ifix
      integer iter, last, m1, m2, ns, nleft
      real abs, cut, cmad, r, d1, d2
      real c1, c9, alpha, denom, float
      logical ok
      if (n .ge. 2) goto 1
         ys(1) = y(1)
         return
c at least two, at most n points
 1    continue
      if(f.lt.0) then
         ns=nint(-f)
      else
      ns = max0(min0(ifix(f*float(n)), n), 2)
      endif
      iter = 1
         goto  3
   2     iter = iter+1
   3     if (iter .gt. nsteps+1) goto  22
c robustness iterations
         nleft = 1
         nright = ns
c index of prev estimated point
         last = 0
c index of current point
         i = 1
   4        if (nright .ge. n) goto  5
c move nleft, nright to right if radius decreases
               d1 = x(i)-x(nleft)
c if d1<=d2 with x(nright+1)==x(nright), lowest fixes
               d2 = x(nright+1)-x(i)
               if (d1 .le. d2) goto  5
c radius will not decrease by move right
               nleft = nleft+1
               nright = nright+1
               goto  4
c fitted value at x(i)
   5        call lowest(x, y, n, x(i), ys(i), nleft, nright, res, iter
     +     .gt. 1, rw, ok)
c            print *,'out of lowest',i,x(i),y(i),ys(i),res(i)
            if (.not. ok) ys(i) = y(i)
c all weights zero - copy over value (all rw==0)
            if (last .ge. i-1) goto 9
               denom = x(i)-x(last)
c skipped points -- interpolate
c non-zero - proof?
               j = last+1
                  goto  7
   6              j = j+1
   7              if (j .ge. i) goto  8
                  alpha = (x(j)-x(last))/denom
                  ys(j) = alpha*ys(i)+(1.0-alpha)*ys(last)
                  goto  6
   8           continue
c last point actually estimated
   9        last = i
c x coord of close points
            cut = x(last)+delta
            i = last+1
               goto  11
  10           i = i+1
  11           if (i .gt. n) goto  13
c find close points
               if (x(i) .gt. cut) goto  13
c i one beyond last pt within cut
               if (x(i) .ne. x(last)) goto 12
                  ys(i) = ys(last)
c exact match in x
                  last = i
  12           continue
               goto  10
c back 1 point so interpolation within delta, but always go forward
  13        i = max0(last+1, i-1)
  14        if (last .lt. n) goto  4
c residuals
         do  15 i = 1, n
            res(i) = y(i)-ys(i)
  15        continue
         if (iter .gt. nsteps) goto  22
c compute robustness weights except last time
         do  16 i = 1, n
c            rw(i) = abs(res(i))
            rw(i) = yerr(i)
  16        continue
         call sort(n,rw)
         m1 = n/2+1
         m2 = n-m1+1
c 6 median abs resid
c         cmad = 3.0*(rw(m1)+rw(m2))
         cmad = 1.5*rw(n)
         c9 = .999*cmad
         c1 = .001*cmad
         do  21 i = 1, n
c            r = abs(res(i))
            r = yerr(i)
            if (r .gt. c1) goto 17
               rw(i) = 1.
c near 0, avoid underflow
               goto  20
  17           if (r .le. c9) goto 18
                  rw(i) = 0.
c near 1, avoid underflow
                  goto  19
  18              rw(i) = (1.0-(r/cmad)**2)**2
c                  print *,i,cmad,yerr(i),rw(i)
  19        continue
  20        continue
  21        continue
         goto  2
  22  return
      end
