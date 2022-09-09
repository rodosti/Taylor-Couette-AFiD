!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsQT.F90                            !
!    CONTAINS: subroutine ExplicitTermsQT                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the t (theta) dimension             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ExplicitTermsQT
      use param
      use local_arrays, only: qr,qz,dph,qt,dq
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      real    :: h11,h12,h13,udxt,udxz
      real    :: udxtq,udxzq, qre,qrw, coriol
      real    :: d11qt,d22qt,d11qre,h12d, h12n

      udxzq=dxzq/ren
      udxz=dxz*0.25

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,qt,qr,qz,dxt,dxz) &
!$OMP& SHARED(kmv,kpv,ren,dxtq,udx1vis1,rm,romeg) &
!$OMP& SHARED(udxz,udxzq,dq,nrm,usrc,usrm,udx2rm) &
!$OMP& PRIVATE(ic,jc,kc,im,ip,kmm,kp,kpp) &
!$OMP& PRIVATE(jm,jp,udxt,udxtq,coriol) &
!$OMP& PRIVATE(h11,h12,h13,d11qt,d22qt) &
!$OMP& PRIVATE(h12d,h12n,qre,qrw,d11qre)
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
         udxt=dxt*0.25*usrm(kc)
         udxtq=dxtq/(ren*rm(kc)**2) 


!
!             1   d  q_t q_t 
!            --- -----------  azimuthal non-linear term
!             r   d   t      
!
        h11=((qt(kc,jc,ip)+qt(kc,jc,ic))*(qt(kc,jc,ip)+qt(kc,jc,ic)) & 
            -(qt(kc,jc,ic)+qt(kc,jc,im))*(qt(kc,jc,ic)+qt(kc,jc,im)) & 
          )*udxt



!
!    qt qr term
!
!             1   d  q_t q_r     q_t q_r
!            --- -----------  +  --------    kc > 2
!             r   d   r            r^2
!
!

        h12d=((qr(kp,jc,ic)+qr(kp,jc,im))*(qt(kpp,jc,ic)+qt(kc,jc,ic)) & 
             -(qr(kc,jc,ic)+qr(kc,jc,im))*(qt(kc,jc,ic)+qt(kmm,jc,ic)) & 
              )*0.25*udx2rm(kc)

         h12n=  qt(kc,jc,ic)*usrm(kc) * ( & 
                           (qr(kp,jc,ic)+qr(kp,jc,im))*usrc(kp) & 
                          +(qr(kc,jc,ic)+qr(kc,jc,im))*usrc(kc)  )*0.25
         h12=h12d+h12n

      

!
!   qt qz term
!
!
!                 d  q_t q_x 
!                -----------
!                 d   x      
!
         h13=((qz(kc,jp,ic)+qz(kc,jp,im))*(qt(kc,jp,ic)+qt(kc,jc,ic)) &
             -(qz(kc,jc,ic)+qz(kc,jc,im))*(qt(kc,jc,ic)+qt(kc,jm,ic)) &
             )*udxz


!
!   second derivative of qt azimuthally
!
            d11qt=(qt(kc,jc,ip)   &
                  -2.0*qt(kc,jc,ic)   &
                  +qt(kc,jc,im))*udxtq   


!
!   second derivative of qt axially
!
            d22qt=(qt(kc,jp,ic)  & 
                  -2.0*qt(kc,jc,ic)  & 
                  +qt(kc,jm,ic))*udxzq

!
!   first derivative of qr with respect to x1
!
!
!              2      d  q_r / r
!            ------  -----------
!            Re r^2   d   t      
!
        qre=qr(kp,jc,ic)*usrc(kp)+qr(kc,jc,ic)*usrc(kc)
        qrw=qr(kp,jc,im)*usrc(kp)+qr(kc,jc,im)*usrc(kc)
        d11qre=(qre-qrw)*udx1vis1(kc)*0.5


      coriol=   -((qr(kc,jc,ic)+qr(kc,jc,im))*usrc(kc) +  & 
                 (qr(kp,jc,ic)+qr(kp,jc,im))*usrc(kp)  )*0.25*romeg 


 
        dq(kc,jc,ic)=-(h11+h12+h13)+d11qre+d22qt+d11qt+coriol
 
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
 
