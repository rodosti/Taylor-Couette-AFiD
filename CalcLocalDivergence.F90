!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcLocalDivergence.F90                        !
!    CONTAINS: subroutine CalcLocalDivergence             !
!                                                         ! 
!    PURPOSE: Calculates the local divergence of the      !
!     intermediate velocity field for the correction step !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcLocalDivergence
      use param
      use local_arrays, only: qt,qr,qz,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real    :: usdt,dqcap   

      usdt = 1.d0/(al*dt)

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,qt,qr,qz,dxt,dxz,usdt) &
!$OMP& SHARED(dph,nrm,xend,usrm,udx2m) &
!$OMP& PRIVATE(ic,jc,kc,ip,jp,kp) &
!$OMP& PRIVATE(dqcap)
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,nrm
             kp=kc+1
             dqcap= ((qt(kc,jc,ip)-qt(kc,jc,ic))*dxt  &
                    +(qr(kp,jc,ic)-qr(kc,jc,ic))*udx2m(kc))*usrm(kc)  &
                    +(qz(kc,jp,ic)-qz(kc,jc,ic))*dxz

              dph(kc,jc,ic)=dqcap*usdt
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
