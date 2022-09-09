!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         ! 
!    PURPOSE: Initialization routine. Allocates and       !
!       zeroes all variables used in the code             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine InitVariables
      use param
      use mpih
      use mpi_param
      use local_arrays
      use decomp_2d
      use stat_arrays
      use AuxiliaryRoutines
      implicit none
      

      !-------------------------------------------------
      ! Arrays with ghost cells
      !-------------------------------------------------
      call AllocateReal3DArray(qt,1,nr,xstart(2)-lvlhalo, & 
       xend(2)+lvlhalo, xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(qr,1,nr,xstart(2)-lvlhalo, & 
       xend(2)+lvlhalo, xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(qz,1,nr,xstart(2)-lvlhalo, & 
       xend(2)+lvlhalo, xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(pr,1,nr,xstart(2)-lvlhalo, & 
       xend(2)+lvlhalo, xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(dph,1,nr,xstart(2)-lvlhalo, & 
       xend(2)+lvlhalo, xstart(3)-lvlhalo,xend(3)+lvlhalo)

      !-----------------------------------------------
      ! Arrays without ghost cells
      !-----------------------------------------------
      call AllocateReal3DArray(rhs,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(dq,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(qcap,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(ru1,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(ru2,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(ru3,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(vortt,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(vortr,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))
      call AllocateReal3DArray(vortz,1,nr,xstart(2),xend(2), &
       xstart(3),xend(3))

      call AllocateReal1DArray(tc,1,nth)
      call AllocateReal1DArray(tm,1,nth)
      call AllocateReal1DArray(rc,1,nr)
      call AllocateReal1DArray(rm,1,nr)
      call AllocateReal1DArray(zc,1,nz)
      call AllocateReal1DArray(zm,1,nz)

      call AllocateReal1DArray(g2rm,1,nr)
      call AllocateReal1DArray(g2rc,1,nr)
      call AllocateReal1DArray(g2rmm,1,nr)
      call AllocateReal1DArray(g2rcc,1,nr)
      call AllocateReal1DArray(pvol,1,nr)
      call AllocateReal1DArray(usrc,1,nr)
      call AllocateReal1DArray(usrm,1,nr)

      call AllocateReal1DArray(udx1vis1,1,nr)
      call AllocateReal1DArray(udx1vis2,1,nr)
      call AllocateReal1DArray(udx1qm,1,nr)
      call AllocateReal1DArray(udx1qc,1,nr)
      call AllocateReal1DArray(udx2c,1,nr)
      call AllocateReal1DArray(udx2m,1,nr)
      call AllocateReal1DArray(udx2cc,1,nr)
      call AllocateReal1DArray(udx2mm,1,nr)
      call AllocateReal1DArray(rudx2c,1,nr)
      call AllocateReal1DArray(udx2rm,1,nr)


      call AllocateInt1DArray(kmc,1,nr)
      call AllocateInt1DArray(kpc,1,nr)
      call AllocateInt1DArray(kmv,1,nr)
      call AllocateInt1DArray(kpv,1,nr)
      call AllocateInt1DArray(kup,1,nr)
      call AllocateInt1DArray(kum,1,nr)

      call AllocateReal1DArray(ap1j,1,nr)
      call AllocateReal1DArray(am1j,1,nr)
      call AllocateReal1DArray(ac1j,1,nr)
      call AllocateReal1DArray(ap2j,1,nr)
      call AllocateReal1DArray(am2j,1,nr)
      call AllocateReal1DArray(ac2j,1,nr)
      call AllocateReal1DArray(ap2je,1,nr)
      call AllocateReal1DArray(am2je,1,nr)
      call AllocateReal1DArray(ac2je,1,nr)
      call AllocateReal1DArray(ap3j,1,nr)
      call AllocateReal1DArray(am3j,1,nr)
      call AllocateReal1DArray(ac3j,1,nr)
      call AllocateReal1DArray(apph,1,nr)
      call AllocateReal1DArray(amph,1,nr)
      call AllocateReal1DArray(acph,1,nr)
      call AllocateReal1DArray(ak1,1,nth)
      call AllocateReal1DArray(ap,1,nth)
      call AllocateReal1DArray(ak3,1,nz)
      call AllocateReal1DArray(ao,1,nz)

      call AllocateInt2DArray(icylsl,1,nzm,1,nthm)
      call AllocateInt2DArray(ocylsl,1,nzm,1,nthm)

      call AllocateReal2DArray(vt_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vr_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vz_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(wt_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(wr_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(wz_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(pr_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vt_rms,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vr_rms,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vz_rms,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(wt_rms,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(wr_rms,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(wz_rms,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(pr_rms,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vtvr_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vtvz_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(vrvz_me,1,nrm,xstart(2),xend(2))
      call AllocateReal2DArray(disste,1,nrm,xstart(2),xend(2))
          
      end subroutine InitVariables

