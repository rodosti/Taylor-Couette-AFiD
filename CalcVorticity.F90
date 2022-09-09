!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcVorticity.F90                              !
!    CONTAINS: subroutine CalcVorticity                   !
!                                                         ! 
!    PURPOSE: Calculates the three components of          !
!     vorticity using the instantaneous velocity          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcVorticity
      use param
      use local_arrays, only: vortt,vortr,vortz,qr,qz,qt
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      real    :: durdz,duzdr,dutdr,durdt,duzdt,dutdz


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(qt,qr,qz,usrc,usrm,dxz,udx2mm) &
!$OMP& SHARED(dxt,rm,nrm,kmv,kpv) &
!$OMP& SHARED(vortt,vortr,vortz,xstart,xend) &
!$OMP& PRIVATE(im,ip,ic,jm,jp,jc,kc,kmm,kp,kpp) &
!$OMP& PRIVATE(duzdr,duzdt,dutdz,dutdr,durdt,durdz)
      do ic=xstart(3),xend(3)
       im=ic-1
       ip=ic+1
       do jc=xstart(2),xend(2)
        jm=jc-1
        jp=jc+1
        do kc=1,nrm
         kmm=kmv(kc)
         kpp=kpv(kc)
         kp=kc+1

!     Azimuthal Vorticity
!     d u_r    d u_z
!     ----  - -----
!      d z      d r


      durdz  = ((qr(kp,jp,ic)*usrc(kp)+qr(kc,jp,ic)*usrc(kc)) & 
               -(qr(kp,jm,ic)*usrc(kp)+qr(kc,jm,ic)*usrc(kc)))*dxz*0.25

      duzdr  = ((qz(kpp,jc,ic)+qz(kpp,jp,ic)) & 
               -(qz(kmm,jc,ic)+qz(kmm,jp,ic)))*0.5*udx2mm(kc)

      vortt(kc,jc,ic) = durdz - duzdr

!     Radial Vorticity
!     1 d uz    d ut
!     ------  - ------
!     r d t     d z

      duzdt  = ((qz(kc,jc,ip)+qz(kc,jp,ip)) & 
               -(qz(kc,jc,im)+qz(kc,jp,im)))*0.25*dxt*usrm(kc)

      dutdz  = ((qt(kc,jp,ic)+qt(kc,jp,ip)) & 
               -(qt(kc,jm,ic)+qt(kc,jm,ip)))*0.25*dxz 

      vortr(kc,jc,ic) = duzdt - dutdz

!     Axial Vorticity
!  1  d (r u_t)      1 d u_r
!  -  ---------   -  - -----
!  r     d r         r  dt


       dutdr = ((qt(kpp,jc,ic)+qt(kpp,jc,ip))*rm(kpp) & 
               -(qt(kmm,jc,ic)+qt(kmm,jc,ip))*rm(kmm)) & 
               *0.5*udx2mm(kc)*usrm(kc)

       durdt = ((qr(kc,jc,ip)*usrc(kc)+qr(kp,jc,ip)*usrc(kp)) & 
               -(qr(kc,jc,im)*usrc(kc)+qr(kp,jc,im)*usrc(kp)))*dxt*0.25 

       vortz(kc,jc,ic) = dutdr - durdt


      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
