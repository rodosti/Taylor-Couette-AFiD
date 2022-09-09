!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         ! 
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity in time            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TimeMarcher
      use param
      use local_arrays
      use mpih
      use decomp_2d
      implicit none
      integer :: ns
      
      beta=dt/ren*0.5d0

      call InnerCylVelocity
      call SlipBCs

      do ns=1,nsst                                                 

        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        call ExplicitTermsQT
        call ExplicitTermsQR
        call ExplicitTermsQZ

        call ImplicitAndUpdateQT
        call ImplicitAndUpdateQR
        call ImplicitAndUpdateQZ

        call update_halo(qt,lvlhalo,.TRUE.)
        call update_halo(qr,lvlhalo,.TRUE.)
        call update_halo(qz,lvlhalo,.TRUE.)

        call CalcLocalDivergence
        call SolvePressureCorrection

        call update_halo(dph,lvlhalo,.TRUE.)
        
        call CorrectVelocity 
        call CorrectPressure 

        call update_halo(qt,lvlhalo,.TRUE.)
        call update_halo(qr,lvlhalo,.TRUE.)
        call update_halo(qz,lvlhalo,.TRUE.)
        call update_halo(pr,lvlhalo,.TRUE.)
 

      enddo

        if(mod(time,tout).lt.dt) then
         if(ismaster) then
          write(6,*) ' ---------------------------------------- '
          write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
         endif
        endif

      return                                                            
      end                                                               
