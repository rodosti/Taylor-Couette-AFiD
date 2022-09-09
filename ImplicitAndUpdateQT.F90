!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateQT.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateQT             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the t (azimuthal) direction and     !
!     call the implicit solver. After this routine, the   !
!     azimuthal velocity has been updated to the new      !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateQT
      use param
      use local_arrays, only: qt,dq,ru1,rhs,pr
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,jc,ic,imm
      integer :: kmm,kpp
      real    :: alre,amm,acc,app
      real    :: d33qt,dpx11

      real :: mck1

      alre=al/ren


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,qt,pr) &
!$OMP& SHARED(kmv,kpv,usrm) &
!$OMP& SHARED(dxt,al,ga,ro,alre,dt,dq) &
!$OMP& SHARED(rhs,ru1,virotext,virotint) &
!$OMP& SHARED(ap1j,ac1j,am1j,vinner_t) &
!$OMP& SHARED(icylsl,ocylsl,ap1j_fs,ac1j_fs,am1j_fs) &
!$OMP& PRIVATE(ic,jc,kc,imm,kmm,kpp) &
!$OMP& PRIVATE(amm,acc,app) &
!$OMP& PRIVATE(d33qt,dpx11)
      do ic=xstart(3),xend(3)
      imm=ic-1
      do jc=xstart(2),xend(2)

!RO   Inner cylinder
      kc=1
      kpp=2
      app=ap1j(kc)
      acc=(ac1j(kc)*float(icylsl(jc,ic)) +  &
           ac1j_fs(1)*(1-float(icylsl(jc,ic))))
      amm=(am1j(kc)*float(icylsl(jc,ic)) +  &
           am1j_fs(1)*(1-float(icylsl(jc,ic))))

!   radial derivative of qt

       d33qt=qt(kpp,jc,ic)*app &
            +qt(kc,jc,ic)*acc &
            +(virotint+vinner_t)*amm


!   grad(pr) along theta direction

        dpx11=(pr(kc,jc,ic)-pr(kc,jc,imm))*dxt*al*usrm(kc)


        rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ru1(kc,jc,ic)  &
                          +alre*d33qt-dpx11)*dt

        ru1(kc,jc,ic)=dq(kc,jc,ic)

!RO   inner points
      do kc=2,nrm-1
      kmm=kc-1
      kpp=kc+1
      amm=am1j(kc)
      acc=ac1j(kc)
      app=ap1j(kc)


!   radial derivative of qt
            d33qt=qt(kpp,jc,ic)*app+qt(kc,jc,ic)*acc+qt(kmm,jc,ic)*amm
      
!   grad(pr) along theta direction
            dpx11=(pr(kc,jc,ic)-pr(kc,jc,imm))*dxt*al*usrm(kc)


            rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ru1(kc,jc,ic) &
                          +alre*d33qt-dpx11)*dt

            ru1(kc,jc,ic)=dq(kc,jc,ic)
      enddo


!RO   Outer cylinder

        kc=nrm
        kmm=nrm-1
        amm=am1j(kc)
        app=(ap1j(kc)*float(ocylsl(jc,ic)) + & 
             ap1j_fs(2)*(1-float(ocylsl(jc,ic))))
        acc=(ac1j(kc)*float(ocylsl(jc,ic)) + & 
             ac1j_fs(2)*(1-float(ocylsl(jc,ic))))

!   radial derivative of qt

        d33qt=virotext*app+qt(kc,jc,ic)*acc+qt(kmm,jc,ic)*amm

!   grad(pr) along theta direction

        dpx11=(pr(kc,jc,ic)-pr(kc,jc,imm))*dxt*al*usrm(kc)


        rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ru1(kc,jc,ic) &
                           +alre*d33qt-dpx11)*dt

        ru1(kc,jc,ic)=dq(kc,jc,ic)

       enddo
      enddo
!$OMP END PARALLEL DO

      call SolveImpEqnUpdate_QT

      return
      end
