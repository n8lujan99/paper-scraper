C=============================================================================D
C     subroutine VELMTOD finds which velocity bins of the data correspond
C       to the velocity bins edges of the model
C
C=============================================================================D
      SUBROUTINE velmtod()
      INCLUDE 'moddefs.h'
      real v1(nad),v2(nad)

      vbin=velm(2)-velm(1)

      do irc=1,Nvelb
c-- first find velocity edges for data
         v1(1)=veld(1,irc)*xmu-(veld(2,irc)-veld(1,irc))*xmu/2.
         v2(1)=veld(1,irc)*xmu+(veld(2,irc)-veld(1,irc))*xmu/2.
         v1(nad)=veld(nad,irc)*xmu-
     $        (veld(nad,irc)-veld(nad-1,irc))*xmu/2.
         v2(nad)=veld(nad,irc)*xmu+
     $        (veld(nad,irc)-veld(nad-1,irc))*xmu/2.
         do k=2,nad-1
            v1(k)=veld(k,irc)*xmu-(veld(k,irc)-veld(k-1,irc))*xmu/2.
            v2(k)=veld(k,irc)*xmu+(veld(k+1,irc)-veld(k,irc))*xmu/2.
         enddo
         sum1=0.
         do i=1,nad
            sum1=sum1+ad(i,irc)
         enddo
         sum2=0.
c-- now find it for the model and then sum the data
         do i=1,Nvel
            vpos1=velm(i)-vbin/2.
            vpos2=velm(i)+vbin/2.
            sum=0.
            sume=0.
            do k=1,nad
               frac=0.
               if(vpos1.lt.v1(k).and.vpos2.gt.v2(k)) frac=1.
               if(vpos1.gt.v1(k).and.vpos1.lt.v2(k)) 
     $              frac=(v2(k)-vpos1)/(v2(k)-v1(k))
               if(vpos2.gt.v1(k).and.vpos2.lt.v2(k)) 
     $              frac=(vpos2-v1(k))/(v2(k)-v1(k))
               sum=sum+frac*ad(k,irc)
               sume=sume+frac*adfer(k,irc)
            enddo
            if(velm(i).lt.v1(1).or.velm(i).gt.v2(nad)) then
               sum=0.
               sume=-666.
            endif
            sumad(i,irc)=sum
            sadfer(i,irc)=sume
            sum2=sum2+sumad(i,irc)
         enddo
       
         do i=1,Nvel
c            sumad(i,irc)=sumad(i,irc)*sum1/sum2
c            sadfer(i,irc)=sadfer(i,irc)*sum1/sum2
         enddo

      enddo

      return
      end
