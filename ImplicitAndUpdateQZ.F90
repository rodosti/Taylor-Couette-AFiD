!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateQZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateQZ             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (axial) direction and         !
!     call the implicit solver. After this routine, the   !
!     axial velocity has been updated to the new          !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateQZ
      use param
      use local_arrays, only: qz,rhs,ru3,qcap,pr
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc,jm
      integer :: km,kp,ic
      real    :: alre,udxz
      real    :: amm,acc,app,dpx33,dqz3

      alre=al/ren
      udxz = dxz*al

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,qz,pr) &
!$OMP& SHARED(kmv,kpv,am3j_fs,ac3j_fs,ap3j_fs) &
!$OMP& SHARED(al,ga,ro,alre,dt,qcap) &
!$OMP& SHARED(rhs,ru3,am3j,ac3j,ap3j,icylsl) &
!$OMP& PRIVATE(ic,jc,kc,kp,jm,km) &
!$OMP& PRIVATE(amm,acc,app,udxz,ocylsl) &
!$OMP& PRIVATE(dqz3,dpx33)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      jm=jc-1

!     Inner cylinder

      kc = 1
      km=kmv(kc)
      kp=kpv(kc)
      app=ap3j(kc)
      acc=(ac3j(kc)*icylsl(jc,ic) + & 
           ac3j_fs(1)*(1-icylsl(jc,ic)))
      amm=(am3j(kc)*icylsl(jc,ic) + & 
           am3j_fs(1)*(1-icylsl(jc,ic)))

!   second radial derivatives of qz
      dqz3=qz(kp,jc,ic)*app & 
          +qz(kc,jc,ic)*acc & 
          +qz(km,jc,ic)*amm

!   pressure gradient

      dpx33=(pr(kc,jc,ic)-pr(kc,jm,ic))*udxz

      rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*ru3(kc,jc,ic) & 
                          +alre*dqz3-dpx33)*dt 

      ru3(kc,jc,ic)=qcap(kc,jc,ic)

!     Inner points

      do kc=2,nrm
      km=kmv(kc)
      kp=kpv(kc)
      amm=am3j(kc)
      acc=ac3j(kc)
      app=ap3j(kc)

!   second radial derivatives of qz
            dqz3=qz(kp,jc,ic)*app & 
                +qz(kc,jc,ic)*acc & 
                +qz(km,jc,ic)*amm

!  component of grad(pr) along axial direction

            dpx33=(pr(kc,jc,ic)-pr(kc,jm,ic))*udxz

            rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*ru3(kc,jc,ic) & 
                          +alre*dqz3-dpx33)*dt 

            ru3(kc,jc,ic)=qcap(kc,jc,ic)
      enddo


!     Outer cylinder

      kc = nrm
      km=kmv(kc)
      kp=kpv(kc)
      app=ap3j(kc)
      app=(ap3j(kc)*float(ocylsl(jc,ic)) + & 
           ap3j_fs(2)*(1-float(ocylsl(jc,ic))))
      acc=(ac3j(kc)*float(ocylsl(jc,ic)) + & 
           ac3j_fs(2)*(1-float(ocylsl(jc,ic))))


!   second radial derivatives of qz
      dqz3=qz(kp,jc,ic)*app &
          +qz(kc,jc,ic)*acc &
          +qz(km,jc,ic)*amm

!    pressure gradient
      dpx33=(pr(kc,jc,ic)-pr(kc,jm,ic))*udxz

      rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*ru3(kc,jc,ic) &
             +alre*dqz3-dpx33)*dt 

      ru3(kc,jc,ic)=qcap(kc,jc,ic)

      enddo
      enddo
!$OMP END PARALLEL DO

      call SolveImpEqnUpdate_QZ


      return
      end
