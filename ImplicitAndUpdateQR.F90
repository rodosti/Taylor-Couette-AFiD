!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateQR.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateQR             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the R (radial) direction and call   !
!     the implicit solver. After this routine, the        !
!     radial velocity has been updated to the new         !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine ImplicitAndUpdateQR
      use param
      use local_arrays, only: qr,ru2,pr,rhs,dph
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: kc,jc,ic
      integer :: kpp,kmm
      real    :: alre,udx2
      real    :: amm,acc,app
      real    :: d33qr,dpx22


      alre=al/ren

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,qr,pr) &
!$OMP& SHARED(kmv,kpv,am2je,ac2je,ap2je) &
!$OMP& SHARED(dxr,al,ga,ro,alre,dt,dph) &
!$OMP& SHARED(rhs,ru2,rc,udx2c) &
!$OMP& PRIVATE(ic,jc,kc,kmm,kpp) &
!$OMP& PRIVATE(amm,acc,app) &
!$OMP& PRIVATE(d33qr,dpx22,udx2)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,nrm
       kmm=kc-1
       kpp=kc+1
       amm=am2je(kc)
       acc=ac2je(kc)
       app=ap2je(kc)
       udx2 = rc(kc)*al*udx2c(kc)


!   second derivative of qr radially

            d33qr=qr(kpp,jc,ic)*app & 
                 +qr(kc,jc,ic)*acc & 
                 +qr(kmm,jc,ic)*amm


!   component of grad(pr) along radial direction
            dpx22=(pr(kc,jc,ic)-pr(kmm,jc,ic))*udx2

            rhs(kc,jc,ic)=(ga*dph(kc,jc,ic)+ro*ru2(kc,jc,ic) & 
                          +alre*d33qr-dpx22)*dt

            ru2(kc,jc,ic)=dph(kc,jc,ic)
      enddo
      enddo
      enddo

      call SolveImpEqnUpdate_QR

      qr(1,:,:)=0.0d0
      qr(nr,:,:)=0.0d0

      
      return
      end
