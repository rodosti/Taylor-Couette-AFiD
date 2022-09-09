!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InnerCylVelocity.F90                           !
!    CONTAINS: subroutine InnerCylVelocity                !
!                                                         ! 
!    PURPOSE: Update time-dependent inner cylinder        !
!       velocity save change for implicit terms           ! 
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InnerCylVelocity
      use param
      implicit none

      vinner_t = sizepert*sin(2.*pi*omegapert*time) 
      dvinner_t= 2.*pi*omegapert*sizepert*cos(2.*pi*omegapert*time)*dt
       
      return
      end
