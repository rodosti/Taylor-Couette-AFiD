!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatReadReduceWrite.F90                        !
!    CONTAINS: subroutine StatReadReduceWrite             !
!                                                         ! 
!    PURPOSE: Reduce statistics across processors, read   !
!     old statistics if required and write them out       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine StatReadReduceWrite(var,filename,dsetname)
      use param
      use decomp_2d, only: xstart,xend,nrank
      use mpih
      implicit none
      integer :: jc,kc
      character*30,intent(in) :: filename, dsetname
      real :: var(1:nrm,xstart(2):xend(2))
      real :: var_old(1:nrm,1:nzm)
      real :: var_new(1:nrm,1:nzm)

      var_new=0.0d0
      var_new(1:nrm,xstart(2):xend(2))=var

      call MpiSumReal1D(var_new,nrm*nzm)

      if (ismaster) then

       if(readstats) then
        call HdfReadSerialReal2D(filename,dsetname,nrm,nzm,var_old)
        var_new = var_new + var_old 
       endif 
 
       call HdfWriteSerialReal2D(filename,dsetname,nrm,nzm,var_new)

      end if

      return
      end

