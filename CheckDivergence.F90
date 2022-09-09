!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CheckDivergence.F90                            !
!    CONTAINS: subroutine CheckDivergence                 !
!                                                         ! 
!    PURPOSE: Check the maximum divergence of velocity    !
!     in the domain                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CheckDivergence(qmax)
      use param
      use local_arrays, only: qr,qz,qt
      use mpih
      use decomp_2d, only: xstart,xend,DECOMP_2D_COMM_CART_X
      implicit none
      real,intent(out) :: qmax
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap
        
      qmax=0.d0                                                     
!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,qt,qr,qz,dxt,dxz,usrm,udx2m) &
!$OMP& PRIVATE(ic,jc,kc,ip,jp,kp) &
!$OMP& PRIVATE(dqcap) &
!$OMP& REDUCTION(max:qmax)
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,nrm
            kp=kc+1
             dqcap=((qt(kc,jc,ip)-qt(kc,jc,ic))*dxt &
                    +(qr(kp,jc,ic)-qr(kc,jc,ic))*udx2m(kc))*usrm(kc) &
                    +(qz(kc,jp,ic)-qz(kc,jc,ic))*dxz
              qmax = max(abs(dqcap),qmax)          
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
     
      call MpiAllMaxRealScalar(qmax)
      
      return     
      end         
