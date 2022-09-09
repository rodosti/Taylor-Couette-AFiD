!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_QZ.F90                       !
!    CONTAINS: subroutine SolveImpEqnUpdate_QZ            !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in the axial direction, and updates it to           !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SolveImpEqnUpdate_QZ
      use param
      use local_arrays, only : rhs,qz
      use decomp_2d, only: xstart,xend
      implicit none
      real, dimension(nr) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,info,ipkv(nrm),ic
      real :: betadx,ackl_b,actot
      real :: amkT(nrm-1),ackT(nrm),apkT(nrm-1),appk(nr-3)

      betadx=beta*al


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,rhs,ap3j,am3j) &
!$OMP& SHARED(ac3j,amkT,ackT,apkT,appk) &
!$OMP& SHARED(ipkv,qz,betadx,icylsl,ocylsl) &
!$OMP& SHARED(ac3j_fs,am3j_fs,ap3j_fs) &
!$OMP& PRIVATE(ic,jc,kc,ackl_b,fkl,info,actot,ackl) &
!$OMP& PRIVATE(apkl,amkl)
      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

          actot = (ac3j(1)*icylsl(jc,ic)+ac3j_fs(1)*(1-icylsl(jc,ic)))
          ackl_b=1.0d0/(1.0d0-(actot*betadx))
          ackl(1)=1.0d0
          apkl(1)=-ap3j(1)*betadx*ackl_b
          amkl(1)=0.0d0

          do kc=2,nrm-1
            ackl_b=1.0d0/(1.0d0-ac3j(kc)*betadx)
            amkl(kc)=-am3j(kc)*betadx*ackl_b
            ackl(kc)=1.0d0
            apkl(kc)=-ap3j(kc)*betadx*ackl_b
          enddo

          actot = (ac3j(nrm)*ocylsl(jc,ic)+ac3j_fs(2)*(1-ocylsl(jc,ic)))
          ackl_b=1.0d0/(1.0d0-(actot*betadx))
          amkl(nrm)=-am3j(nrm)*betadx*ackl_b
          ackl(nrm)=1.0d0
          apkl(nrm)=0.0d0

          amkT=amkl(2:nrm)
          apkT=apkl(1:(nrm-1))
          ackT=ackl(1:nrm)

          call dgttrf(nrm,amkT,ackT,apkT,appk,ipkv,info)
          do kc=1,nrm
            ackl_b=1.0d0/(1.-ac3j(kc)*betadx)
            fkl(kc)=rhs(kc,jc,ic)*ackl_b
          end do

          call dgttrs('N',nrm,1,amkT,ackT,apkT,appk,ipkv,fkl,nrm,info)

          do kc=1,nrm
            qz(kc,jc,ic)=qz(kc,jc,ic) + fkl(kc)
          end do
        end do
      end do
!$OMP END PARALLEL DO

      return
      end
