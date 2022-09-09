!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsQR.F90                            !
!    CONTAINS: subroutine ExplicitTermsQR                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the r (radial) dimension            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine ExplicitTermsQR
      use param
      use local_arrays, only: qr,qz,qt,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kpp,kmm
      real    :: udxtq,udxzq,qte,qtw, coriol
      real    :: h11n,h22,h23,udxt,udxz,h21
      real    :: d11qr,d22qr,d11qre


      udxzq=dxzq/ren
      udxz=dxz*0.25
!
!     h term for the qr momentum equation at i+1/2,j,k+1/2
!

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,qt,qr,qz,dxt) &
!$OMP& SHARED(kmv,kpv,usrc,ren,dxtq,udx1vis2) &
!$OMP& SHARED(udxz,udxzq,dph,nrm,udx2c,rc,romeg) &
!$OMP& PRIVATE(ic,jc,kc,im,ip,kmm,kp,kpp) &
!$OMP& PRIVATE(jm,jp,udxt,udxtq) &
!$OMP& PRIVATE(h21,h22,h23,d11qr,d22qr) &
!$OMP& PRIVATE(h11n,qte,qtw,d11qre,coriol)
      do ic=xstart(3),xend(3)
       im=ic-1
       ip=ic+1
       do jc=xstart(2),xend(2)
        jm=jc-1
        jp=jc+1
        do kc=2,nrm
         udxt=dxt*0.25*usrc(kc)
         udxtq=dxtq/(ren*rc(kc)**2)
         kmm=kmv(kc)
         kpp=kpv(kc)
         kp=kc+1

!
!   qt qt term  not differentiated
!
!             ( q_t )^2
!
      h11n=((qt(kc,jc,ip)+qt(kmm,jc,ip)  &
            +qt(kc,jc,ic)+qt(kmm,jc,ic))*0.25)**2


!
!
!
!             1   d  q_t q_r 
!            --- -----------
!             r   d   t      
!


      h21=( (qt(kc,jc,ip)+qt(kmm,jc,ip))  &
           *(qr(kc,jc,ip)+qr(kc,jc,ic))  &
           -(qt(kc,jc,ic)+qt(kmm,jc,ic))  &
           *(qr(kc,jc,ic)+qr(kc,jc,im)))*udxt


!
!     qr qr term
!
!
!                 d  q_r q_r / r
!                ---------------
!                 d   r      
!
      h22=( (qr(kp,jc,ic)*usrc(kp)+qr(kc,jc,ic)*usrc(kc)) &
           *(qr(kp,jc,ic)+qr(kc,jc,ic)) &
           -(qr(kc,jc,ic)*usrc(kc)+qr(kmm,jc,ic)*usrc(kmm)) &
           *(qr(kc,jc,ic)+qr(kmm,jc,ic)) &
          )*udx2c(kc)*0.25


!
!     qr qz term
!
!
!                 d  q_x q_r 
!                -----------
!                 d   x      
!
      h23=((qz(kc,jp,ic)+qz(kmm,jp,ic))*(qr(kc,jp,ic)+qr(kc,jc,ic))  &
          -(qz(kc,jc,ic)+qz(kmm,jc,ic))*(qr(kc,jc,ic)+qr(kc,jm,ic))  &
          )*udxz



!
!   second derivative of qt with respect to x1
!
!
!              2      d  q_t    
!          - ------  --------
!            Re r     d   t      
!
      qte=(qt(kc,jc,ip)+qt(kmm,jc,ip))
      qtw=(qt(kc,jc,ic) +qt(kmm,jc,ic))
      d11qre=-(qte-qtw)*udx1vis2(kc)


!
!   second derivative of qr azimuthally
!
            d11qr=(qr(kc,jc,ip)  &
              -2.0*qr(kc,jc,ic)  &
                  +qr(kc,jc,im))*udxtq

!
!   second derivative of qr axially
!
            d22qr=(qr(kc,jp,ic)  &
              -2.0*qr(kc,jc,ic)  &
                  +qr(kc,jm,ic))*udxzq

!
!
!   coriolis term
!
!
!              1            
!             -----   r q_t
!              Ro             
!
      coriol= (qt(kc,jc,ic)+qt(kc,jc,ip)+qt(kmm,jc,ic)+  &
               qt(kmm,jc,ip))*0.25*rc(kc)*romeg


            dph(kc,jc,ic)=-(h21+h22+h23)+h11n+d11qr+d22qr+d11qre+coriol
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
 
