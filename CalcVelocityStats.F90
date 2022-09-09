!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcVelocityStats.F90                          !
!    CONTAINS: subroutine CalcVelocityStats               !
!                                                         ! 
!    PURPOSE: Calculates the maximum and rms velocity     !
!     and outputs in max_vel.out and rms_vel.out          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcVelocityStats
      use param
      use local_arrays, only: qr,qz,qt
      use decomp_2d, only: xstart,xend,nrank
      use mpih
      implicit none
      real    :: my_vmax2,my_vmax3,my_vmax1,volt
      real    :: qtrmsinst,qrrmsinst,qzrmsinst
      integer :: jc,kc,ic,ip,jp,kp

      vmax=-100.d0

      qtrmsinst = 0.0d0
      qrrmsinst = 0.0d0
      qzrmsinst = 0.0d0

      volt = pi*(rext**2-rint**2)

      do ic=xstart(3),xend(3)
       ip=ic+1
       do jc=xstart(2),xend(2)
        jp=jc+1
        do kc=1,nrm
         kp=kc+1
         vmax(1) = max(vmax(1),abs(qt(kc,jc,ic)))
         vmax(2) = max(vmax(2),abs(qr(kc,jc,ic)))
         vmax(3) = max(vmax(3),abs(qz(kc,jc,ic)))
         qtrmsinst = qtrmsinst +  &
        (qt(kc,jc,ic)+qt(kc,jc,ip))**2*pvol(kc)*pi*0.25d0 
         qrrmsinst = qrrmsinst +  &
       (qr(kc,jc,ic)*usrc(kc)+qr(kp,jc,ic)*usrc(kp))**2 &
         *pvol(kc)*pi*0.25d0
         qzrmsinst = qzrmsinst +  &
        (qz(kc,jc,ic)+qz(kc,jp,ic))**2*pvol(kc)*pi*0.25d0

        enddo
       enddo
      enddo

      call MpiMaxRealScalar(vmax(1))
      call MpiMaxRealScalar(vmax(2))
      call MpiMaxRealScalar(vmax(3))

      call MpiSumRealScalar(qtrmsinst)
      call MpiSumRealScalar(qrrmsinst)
      call MpiSumRealScalar(qzrmsinst)

      qtrmsinst = qtrmsinst / (float(nthm)*float(nzm)*volt)
      qrrmsinst = qrrmsinst / (float(nthm)*float(nzm)*volt)
      qzrmsinst = qzrmsinst / (float(nthm)*float(nzm)*volt)
       
 546   format(4(1x,e14.6))

      if (ismaster) write(95,546) time, vmax(1), vmax(2), vmax(3)
      if (ismaster) write(98,546) time,qtrmsinst,qrrmsinst,qzrmsinst

      return   
      end     
