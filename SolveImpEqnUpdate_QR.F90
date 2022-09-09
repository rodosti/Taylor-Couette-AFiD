!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_QR.F90                       !
!    CONTAINS: subroutine SolveImpEqnUpdate_QR            !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in the radial direction, and updates it to          !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpEqnUpdate_QR
      use param
      use local_arrays, only : qr,rhs
      use decomp_2d, only: xstart,xend
      implicit none
      real, dimension(nr) :: amkl,apkl,ackl, fkl
      real :: amkT(nr-1),apkT(nr-1)
      real :: appk(nr-2)
      real :: ackT(nr)
      integer :: jc,kc,info,ic
      integer :: ipkv(nr)
      real :: betadx,ackl_b

      betadx=beta*al

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,nrm
        ackl_b=1.0d0/(1.0d0-ac2j(kc)*betadx)
        amkl(kc)=-am2j(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap2j(kc)*betadx*ackl_b
      enddo
      amkl(nr)=0.d0
      apkl(nr)=0.d0
      ackl(nr)=1.d0

      amkT=amkl(2:nr)
      apkT=apkl(1:(nr-1))
      ackT=ackl(1:nr)

      call dgttrf(nr,amkT,ackT,apkT,appk,ipkv,info)

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,rhs) &
!$OMP& SHARED(ac2j,amkT,ackT,apkT,appk) &
!$OMP& SHARED(ipkv,qr,betadx,nr) &
!$OMP& PRIVATE(ic,jc,kc,ackl_b,fkl,info)
      do ic=xstart(3),xend(3)
          do jc=xstart(2),xend(2)
            fkl(1)= 0.d0
          do kc=2,nrm
            ackl_b=1.0d0/(1.0d0-ac2j(kc)*betadx)
            fkl(kc) = rhs(kc,jc,ic)*ackl_b
          enddo
            fkl(nr)= 0.d0
          
          call dgttrs('N',nr,1,amkT,ackT,apkT,appk,ipkv,fkl,nr,info)

          do kc=2,nrm
            qr(kc,jc,ic)=qr(kc,jc,ic) + fkl(kc)
          enddo
          enddo
      end do
!$OMP END PARALLEL DO


      return
      end
