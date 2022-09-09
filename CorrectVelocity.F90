!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine CorrectVelocity
      use param
      use local_arrays, only: qr,qz,dph,qt
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,jm,kc,km,ic,im
      real    :: locdph


!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(qt,qr,qz,dph,dxt,dxz) &
!$OMP& SHARED(xstart,xend,nrm,kmv,dt,al) &
!$OMP& SHARED(usrm,rudx2c) &
!$OMP& PRIVATE(ic,jc,kc) &
!$OMP& PRIVATE(im,jm,km,locdph)
      do ic=xstart(3),xend(3)
        im=ic-1
        do jc=xstart(2),xend(2)
          jm=jc-1
          do kc=1,nrm
          km=kmv(kc)
          locdph=dph(kc,jc,ic)
          qt(kc,jc,ic)=qt(kc,jc,ic)- &
            (locdph-dph(kc,jc,im))*dt*al*dxt*usrm(kc)
          qr(kc,jc,ic)=qr(kc,jc,ic)- &
            (locdph-dph(km,jc,ic))*dt*al*rudx2c(kc)
          qz(kc,jc,ic)=qz(kc,jc,ic)- &
            (locdph-dph(kc,jm,ic))*al*dt*dxz
        enddo 
       enddo
      enddo
!$OMP END PARALLEL DO
      
      qr(1,:,:)=0.0d0
      qr(nr,:,:)=0.0d0

      return
      end

