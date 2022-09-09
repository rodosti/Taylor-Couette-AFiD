!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: OpenLogs.F90                                   !
!    CONTAINS: subroutine OpenLogs, CloseLogs             !
!                                                         ! 
!    PURPOSE: (1) Initialization routine. Open all log    !
!     files and reset if necessary. (2) Finalization      !
!     routine. Close all log files                        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine OpenLogs
      use param
      implicit none

!EP    vmax.out  in vmaxv.F
      open(95,file='vmax.out',status='unknown',access='sequential', &
       position='append')

!EP    torquenu.out  in densmc.F
      open(97,file="torquenu.out",status='unknown',access='sequential', &
       position='append')

!EP    rms_vel.out  in vmaxv.F
      open(98,file="rms_vel.out",status='unknown',access='sequential', &
       position='append')

!EP   dissipnu.out in balance.F
      open(92,file='dissipnu.out',status='unknown',access='sequential', &
       position='append')

      if(resetlogstime) then    
       rewind(92)
       rewind(95)
       rewind(97)
       rewind(98)
      endif

      return
      end   
      


      subroutine CloseLogs
      implicit none

      close(91)
      close(95)
      close(97)
      close(98)

      return 
      end
