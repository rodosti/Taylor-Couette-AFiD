!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_QT.F90                       !
!    CONTAINS: subroutine SolveImpEqnUpdate_QT            !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in the azimuthal direction, and updates it to       !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpEqnUpdate_QT
      use param
      use local_arrays, only : rhs, qt
      use decomp_2d, only: xstart,xend
      implicit none
      real, dimension(nrm) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,info,ipkv(nrm),ic
      real :: betadx,ackl_b,actot
      real :: amkT(nrm-1),ackT(nrm),apkT(nrm-1),appk(nr-3)

      betadx=beta*al


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,rhs) &
!$OMP& SHARED(ac1j,am1j,ap1j,amkT,ackT,apkT,appk) &
!$OMP& SHARED(ipkv,qt,betadx,dt) &
!$OMP& SHARED(icylsl,ocylsl,al,dvinner_t) &
!$OMP& SHARED(ac1j_fs,am1j_fs,ap1j_fs) &
!$OMP& PRIVATE(ic,jc,kc,ackl_b,fkl,info) &
!$OMP& PRIVATE(apkl,amkl,ackl,actot)
      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

          actot = (ac1j(1)*float(icylsl(jc,ic)) &
                  +ac1j_fs(1)*(1-float(icylsl(jc,ic))))
          ackl_b=1.0d0/(1.0d0-actot*betadx)
          ackl(1)=1.0d0
          apkl(1)=-ap1j(1)*betadx*ackl_b
          amkl(1)=(am1j(1)*float(icylsl(jc,ic)) + &
           am1j_fs(1)*(1-float(icylsl(jc,ic))))

          do kc=2,nrm-1
            ackl_b=1.0d0/(1.0d0-ac1j(kc)*betadx)
            amkl(kc)=-am1j(kc)*betadx*ackl_b
            ackl(kc)=1.0d0
            apkl(kc)=-ap1j(kc)*betadx*ackl_b
          enddo

          actot = (ac1j(nrm)*float(ocylsl(jc,ic)) &
                  +ac1j_fs(2)*(1-float(ocylsl(jc,ic))))
          ackl_b=1.0d0/(1.0d0-actot*betadx)
          amkl(nrm)=-am1j(nrm)*betadx*ackl_b
          ackl(nrm)=1.0d0
          apkl(nrm)=0.0d0

          amkT=amkl(2:nrm)
          apkT=apkl(1:(nrm-1))
          ackT=ackl(1:nrm)

          call dgttrf(nrm,amkT,ackT,apkT,appk,ipkv,info)

          ackl_b=1.0d0/(1.-ac1j(1)*betadx)
!          dvinner_t=0.0
          fkl(1)=(rhs(1,jc,ic)+0.5*al*dt*dvinner_t*am1j(1)/ren)*ackl_b

          do kc=2,nrm
            ackl_b=1.0d0/(1.-ac1j(kc)*betadx)
            fkl(kc)=rhs(kc,jc,ic)*ackl_b
          end do

          call dgttrs('N',nrm,1,amkT,ackT,apkT,appk,ipkv,fkl,nrm,info)

          do kc=1,nrm
            qt(kc,jc,ic)=qt(kc,jc,ic) + fkl(kc)
          end do
        end do
      end do
!$OMP END PARALLEL DO

      return
      end
