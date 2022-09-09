!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcStats
      use param
      use local_arrays, only: qt,qr,qz,pr,vortt,vortr,vortz
      use decomp_2d, only: xstart,xend,nrank
      use stat_arrays
      use mpih
      implicit none
      real :: usnthm
      integer :: i,ip,j,k,kp

      nstatsamples = nstatsamples + 1

      usnthm = 1.0/nthm

      call CalcVorticity

      do i=xstart(3),xend(3)
        ip=i+1
        do j=xstart(2),xend(2)
          do k=1,nrm
           kp = kpv(k)
           vt_me(k,j) = vt_me(k,j) + qt(k,j,i)*usnthm
           vr_me(k,j) = vr_me(k,j) + qr(k,j,i)*usnthm*usrc(k)
           vz_me(k,j) = vz_me(k,j) + qz(k,j,i)*usnthm

           wt_me(k,j) = wt_me(k,j) + vortt(k,j,i)*usnthm
           wr_me(k,j) = wr_me(k,j) + vortr(k,j,i)*usnthm
           wz_me(k,j) = wz_me(k,j) + vortz(k,j,i)*usnthm

           vt_rms(k,j) = vt_rms(k,j) + qt(k,j,i)**2*usnthm
           vr_rms(k,j) = vr_rms(k,j) + (qr(k,j,i)*usrc(k))**2*usnthm
           vz_rms(k,j) = vz_rms(k,j) + qz(k,j,i)**2*usnthm

           wt_rms(k,j) = wt_rms(k,j) + vortt(k,j,i)**2*usnthm
           wr_rms(k,j) = wr_rms(k,j) + vortr(k,j,i)**2*usnthm
           wz_rms(k,j) = wz_rms(k,j) + vortz(k,j,i)**2*usnthm

           pr_me(k,j)  = pr_me(k,j)  + pr(k,j,i)*usnthm
           pr_rms(k,j) = pr_rms(k,j) + pr(k,j,i)**2*usnthm

           vtvr_me(k,j) = vtvr_me(k,j) +  &
               (qt(k,j,i)+qt(k,j,i+1)) &
              *(qr(k,j,i)*usrc(k)+qr(kp,j,i)*usrc(kp))*usnthm*0.25

           vtvz_me(k,j) = vtvz_me(k,j) +  &
               (qt(k,j,i)+qt(k,j,i+1)) &
              *(qz(k,j,i)+qr(k,j+1,i))*usnthm*0.25

           vrvz_me(k,j) = vrvz_me(k,j) +  &
               (qz(k,j,i)+qz(k,j+1,i)) &
              *(qr(k,j,i)*usrc(k)+qr(kp,j,i)*usrc(kp))*usnthm*0.25


            end do
         end do
      end do


      return  
      end
 
!***********************************************************************
      subroutine WriteStats
      use mpih
      use param
      use stat_arrays
      use hdf5
      use decomp_2d, only: xstart,xend,nrank

      implicit none

      integer hdf_error, jc, kc, iadd

      integer(HID_T) :: file_id

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_full

      integer :: ndims,nstatsamples_old

      character*30 filnamgrid,filnamdat,dsetnam

      filnamgrid = trim('./stats/stafield_master.h5')

!RO   Reduce, read and write 2D statistics

      dsetnam=trim('vt_mean')
      filnamdat = trim('./stats/vt_mean.h5')
      call StatReadReduceWrite(vt_me,filnamdat,dsetnam)
      dsetnam=trim('vr_mean')
      filnamdat = trim('./stats/vr_mean.h5')
      call StatReadReduceWrite(vr_me,filnamdat,dsetnam)
      dsetnam=trim('vz_mean')
      filnamdat = trim('./stats/vz_mean.h5')
      call StatReadReduceWrite(vz_me,filnamdat,dsetnam)

      dsetnam=trim('wt_mean')
      filnamdat = trim('./stats/wt_mean.h5')
      call StatReadReduceWrite(wt_me,filnamdat,dsetnam)
      dsetnam=trim('wr_mean')
      filnamdat = trim('./stats/wr_mean.h5')
      call StatReadReduceWrite(wr_me,filnamdat,dsetnam)
      dsetnam=trim('wz_mean')
      filnamdat = trim('./stats/wz_mean.h5')
      call StatReadReduceWrite(wz_me,filnamdat,dsetnam)

      dsetnam=trim('vt_rms')
      filnamdat = trim('./stats/vt_rms.h5')
      call StatReadReduceWrite(vt_rms,filnamdat,dsetnam)
      dsetnam=trim('vr_rms')
      filnamdat = trim('./stats/vr_rms.h5')
      call StatReadReduceWrite(vr_rms,filnamdat,dsetnam)
      dsetnam=trim('vz_rms')
      filnamdat = trim('./stats/vz_rms.h5')
      call StatReadReduceWrite(vz_rms,filnamdat,dsetnam)

      dsetnam=trim('wt_rms')
      filnamdat = trim('./stats/wt_rms.h5')
      call StatReadReduceWrite(wt_rms,filnamdat,dsetnam)
      dsetnam=trim('wr_rms')
      filnamdat = trim('./stats/wr_rms.h5')
      call StatReadReduceWrite(wr_rms,filnamdat,dsetnam)
      dsetnam=trim('wz_rms')
      filnamdat = trim('./stats/wz_rms.h5')
      call StatReadReduceWrite(wz_rms,filnamdat,dsetnam)

      dsetnam=trim('vtvr_mean')
      filnamdat = trim('./stats/vtvr_mean.h5')
      call StatReadReduceWrite(vtvr_me,filnamdat,dsetnam)
      dsetnam=trim('vtvz_mean')
      filnamdat = trim('./stats/vtvz_mean.h5')
      call StatReadReduceWrite(vtvz_me,filnamdat,dsetnam)
      dsetnam=trim('vrvz_mean')
      filnamdat = trim('./stats/vrvz_mean.h5')
      call StatReadReduceWrite(vrvz_me,filnamdat,dsetnam)

      dsetnam=trim('pr_mean')
      filnamdat = trim('./stats/pr_mean.h5')
      call StatReadReduceWrite(pr_me,filnamdat,  dsetnam)
      dsetnam=trim('pr_rms')
      filnamdat = trim('./stats/pr_rms.h5')
      call StatReadReduceWrite(pr_rms,filnamdat,dsetnam)
      dsetnam=trim('disste')
      filnamdat = trim('./stats/disste.h5')
      call StatReadReduceWrite(disste,filnamdat,dsetnam)

      dsetnam = trim('averaging_time')
      if (ismaster) then
       if(readstats) then
        call HdfSerialReadIntScalar(dsetnam,filnamgrid,nstatsamples_old)
        nstatsamples = nstatsamples + nstatsamples_old
       else
        call HdfCreateBlankFile(filnamgrid)
       endif

       call HdfSerialWriteIntScalar(dsetnam,filnamgrid,nstatsamples)

       dsetnam = trim('R_cordin')
       call HdfSerialWriteReal1D(dsetnam,filnamgrid,rm,nrm)

       dsetnam = trim('Th_cordin')
       call HdfSerialWriteReal1D(dsetnam,filnamgrid,tm,nthm)

       dsetnam = trim('Z_cordin')
       call HdfSerialWriteReal1D(dsetnam,filnamgrid,zm,nzm)

       dsetnam = trim('Reynolds Number')
       call HdfSerialWriteRealScalar(dsetnam,filnamgrid,ren)

       dsetnam = trim('R_Omega')
       call HdfSerialWriteRealScalar(dsetnam,filnamgrid,romeg)


      end if


      return  
      end
  
!***********************************************************************

      subroutine InitStats
      use param
      use stat_arrays
      implicit none

!EP   Read or initialize stat arrays

      nstatsamples = 0

!EP   Initialize to 0

       vt_me    =0.0d0
       vr_me    =0.0d0
       vz_me    =0.0d0
       wt_me    =0.0d0
       wr_me    =0.0d0
       wz_me    =0.0d0
       pr_me    =0.0d0

       vt_rms   =0.0d0
       vr_rms   =0.0d0
       vz_rms   =0.0d0
       wt_rms   =0.0d0
       wr_rms   =0.0d0
       wz_rms   =0.0d0
       pr_rms   =0.0d0

       vtvr_me  =0.0d0
       vtvz_me  =0.0d0
       vrvz_me  =0.0d0
       disste = 0.0d0

      return
      end
