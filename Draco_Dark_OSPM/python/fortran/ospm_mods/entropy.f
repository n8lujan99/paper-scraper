C=============================================================================D
C     subroutine ENTROPY really calculates the profit function (i.e. the
C       function to be maximized),  which is the "entropy"  minus a chi-
C       square term due to slit errors minus a chi-square term due to the
C       error in the target value of (GM/Ro)^-1.  It also calculates the
C       first and second derivatives of the profit function.
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE entropy(S,dS,ddS,istype)
      INCLUDE 'moddefs.h'
      DOUBLE PRECISION S,dS(Norbitm+Nvel*Nvelbm)
      DOUBLE PRECISION ddS(Norbitm+Nvel*Nvelbm)
      DOUBLE PRECISION ww,wwphase,sum,dNorbit
      DOUBLE PRECISION sumt,sumt2,sumdata

      data big /1.e20/
      common/cfracnew/fracnew

      dNorbit = dble(Norbit)

C---- calculate the "entropy" and its first 2 derivatives

      S=0.d0
      chi=0.

      if (istype.eq.0) then

        do i=1,Norbit
          ww = dble(w(i))
          wwphase = dble(wphase(i))
          S =  S + 1.d0/wwphase - ww*ww*wwphase
          dS(i) = -2.d0*ww*wwphase
          ddS(i) = -2.d0*wwphase
        enddo

      else if (istype.eq.1) then

        do i=1,Norbit
          ww = dble(w(i))
          wwphase = dble(wphase(i))
          S = S + 1.d0/wwphase * dlog(ww*wwphase*dNorbit)
          dS(i) = 1.d0/ww/wwphase
          ddS(i) = -1.d0/ww/ww/wwphase
        enddo
        
      else

        do i=1,Norbit
          ww = dble(w(i))
          wwphase = dble(wphase(i))
          S = S - ww*dlog(ww*wwphase)
          dS(i) = -1.d0 - dlog(ww*wwphase)
          ddS(i) = -1.d0/ww

        enddo
      endif

C---- now calculate contribution of the slit errors

      do ivbin=1,Nvelb
         chi1(ivbin)=0.
         sumt=0.d0
         sumt2=0.d0
         sumdata=0.d0
         do ivel=1,Nvel
            do j=1,Norbit 
               sumt=sumt+dble(v1lib(ivel,ivbin,j)*w(j))
               if(sadfer(ivel,ivbin).ne.-666.) then
                  sumt2=sumt2+dble(v1lib(ivel,ivbin,j)*w(j))
               end if
            enddo
            sumdata=sumdata+dble(sumad(ivel,ivbin))
         enddo
c         fracnew=sumt2/sumdata
         fracnew=sumt2/sumt

         do ivel=1,Nvel

            i=Norbit+(ivbin-1)*Nvel+ivel

            sum=0.d0
            do j=1,Norbit
               sum=sum+dble(v1lib(ivel,ivbin,j)*w(j))
            enddo

            ww=dble(sumad(ivel,ivbin)*fracnew)-sum
            if(sadfer(ivel,ivbin).eq.-666.) then
               den=(1.e6)**2
            else
               den=(sadfer(ivel,ivbin)*fracnew)**2
            end if

c - weight the HST bins more by some factor
c            if(isee(ivbin).eq.1) den=den/2.

            S=S-ww*ww*dble(alphat/den)
            chi=chi+sngl(ww*ww)*alphat/den
            dS(i)=-2.d0*ww*dble(alphat/den)
            ddS(i)=-2.d0*dble(alphat/den)
            if(ddS(i).eq.0.) ddS(i)=1.
            chi1(ivbin)=chi1(ivbin)+sngl(ww*ww)/den
         enddo
      enddo

C---- now calculate contribution of the target (GM/Ro)^-1 error

c      S= S - (dxmu-dxmut)**2/(dbeta*dxmut)**2
c      dS(ilast) = -2.d0*(dxmu-dxmut)/(dbeta*dxmut)**2
c      ddS(ilast) = -2.d0/(dbeta*dxmut)**2

      RETURN
      END
