!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsQZ.F90                            !
!    CONTAINS: subroutine ExplicitTermsQZ                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the z (axial) dimension             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ExplicitTermsQZ
      use param
      use local_arrays, only: qt,qr,qz,qcap
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc,kmm,kpp
      integer :: kp,jm,jp,ic,im,ip
      real    :: h32,h33,h31
      real    :: udxt,udxz
      real    :: udxtq,udxzq
      real    :: dqz1,dqz2

      udxz=dxz*0.25
      udxzq=dxzq/ren

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,qt,qr,qz,dxt) &
!$OMP& SHARED(kmv,kpv,usrm,dxtq,ren,udx2rm) &
!$OMP& SHARED(udxz,udxzq,qcap,rm) &
!$OMP& PRIVATE(ic,jc,kc,im,ip,kmm,kpp,kp) &
!$OMP& PRIVATE(jm,jp,udxt,udxtq) &
!$OMP& PRIVATE(h31,h32,h33,dqz1,dqz2)
      do ic=xstart(3),xend(3)
       im=ic-1
       ip=ic+1
       do jc=xstart(2),xend(2)
        jm=jc-1
        jp=jc+1
        do kc=1,nrm
         kp=kc+1
         kmm=kmv(kc)
         kpp=kpv(kc)
         udxt=dxt*0.25*usrm(kc)
         udxtq=dxtq/(ren*rm(kc)**2)

!
!  right hand side of the momentum equation (qt qz)
!
!
!            1    d  q_x q_t 
!           ---  -----------
!            r    d   t      
!

      h31=((qt(kc,jc,ip)+qt(kc,jm,ip)) &
          *(qz(kc,jc,ip)+qz(kc,jc,ic)) &
          -(qt(kc,jc,ic)+qt(kc,jm,ic)) &
          *(qz(kc,jc,ic)+qz(kc,jc,im)))*udxt


!
!    qz qr term
!
!
!            1    d  q_x q_r 
!           ---  -----------
!            r    d   r      
!
         h32=(((qr(kp,jc,ic)+qr(kp,jm,ic)) &
              *(qz(kpp,jc,ic)+qz(kc,jc,ic))) &
             -((qr(kc,jc,ic)+qr(kc,jm,ic)) &
              *(qz(kc,jc,ic)+qz(kmm,jc,ic))))*udx2rm(kc)*0.25

!
!    qz qz term
!
!
!                 d  q_x q_x 
!                -----------
!                 d   x      
!
         h33=((qz(kc,jp,ic)+qz(kc,jc,ic))*(qz(kc,jp,ic)+qz(kc,jc,ic)) &
             -(qz(kc,jc,ic)+qz(kc,jm,ic))*(qz(kc,jc,ic)+qz(kc,jm,ic)) &
             )*udxz


!
!   second derivatives of qz azimuthally
!
            dqz1=(qz(kc,jc,im) &
             -2.0*qz(kc,jc,ic) &
                 +qz(kc,jc,ip))*udxtq
!
!   second derivatives of qz axially
!
            dqz2=(qz(kc,jm,ic)  &
             -2.0*qz(kc,jc,ic)  &
                 +qz(kc,jp,ic))*udxzq


          qcap(kc,jc,ic) =-(h31+h32+h33)+dqz1+dqz2
            
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
 
