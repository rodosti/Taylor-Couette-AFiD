!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectPressure.F90                            !
!    CONTAINS: subroutine CorrectPressure                 !
!                                                         ! 
!    PURPOSE: Apply the pressure correction to the        !
!     pressure                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CorrectPressure
      use param
      use local_arrays, only: pr,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real    :: be,amm,acc,app
!
!    the pressure is evaluated at the center of the box.
!
!     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
!
      be=beta
!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(pr,dph,be,amph,acph,apph,udx1qm) &
!$OMP& SHARED(xstart,xend,nrm,kmv,kpv,dxzq) &
!$OMP& PRIVATE(ic,jc,kc) &
!$OMP& PRIVATE(im,jm,km,ip,jp,kp) &
!$OMP& PRIVATE(amm,acc,app)
      do ic=xstart(3),xend(3)
       im=ic-1
       ip=ic+1
       do jc=xstart(2),xend(2)
        jm=jc-1
        jp=jc+1
        do kc=1,nrm
         kp=kpv(kc)
         km=kmv(kc)
         amm=apph(kc)
         acc=acph(kc)
         app=apph(kc)
         pr(kc,jc,ic)=pr(kc,jc,ic)+dph(kc,jc,ic)-be*( &
          (dph(kc,jc,ip)-2.0*dph(kc,jc,ic)+dph(kc,jc,im))*udx1qm(kc)+ &
          (dph(kc,jp,ic)-2.0*dph(kc,jc,ic)+dph(kc,jm,ic))*dxzq+ &
          (dph(kp,jc,ic)*app+dph(kc,jc,ic)*acc+dph(km,jc,ic)*amm))
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO
      return
      end
