!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateVariables.F90                        !
!    CONTAINS: subroutine DeallocateVariables             !
!                                                         ! 
!    PURPOSE: Finalizing routine. Deallocates all         !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DeallocateVariables
      use param
      use mpi_param
      use local_arrays
      use decomp_2d
      use stat_arrays
      use AuxiliaryRoutines
      implicit none

      call DestroyReal3DArray(qt)
      call DestroyReal3DArray(qr)
      call DestroyReal3DArray(qz)
      call DestroyReal3DArray(pr)
      call DestroyReal3DArray(dph)

      call DestroyReal3DArray(rhs)
      call DestroyReal3DArray(dq)
      call DestroyReal3DArray(qcap)
      call DestroyReal3DArray(ru1)
      call DestroyReal3DArray(ru2)
      call DestroyReal3DArray(ru3)
      call DestroyReal3DArray(vortt)
      call DestroyReal3DArray(vortr)
      call DestroyReal3DArray(vortz)

      call DestroyReal1DArray(tc)
      call DestroyReal1DArray(tm)
      call DestroyReal1DArray(rc)
      call DestroyReal1DArray(rm)
      call DestroyReal1DArray(zc)
      call DestroyReal1DArray(zm)

      call DestroyReal1DArray(g2rm)
      call DestroyReal1DArray(g2rc)
      call DestroyReal1DArray(pvol)
      call DestroyReal1DArray(usrc)
      call DestroyReal1DArray(usrm)
      call DestroyReal1DArray(g2rmm)
      call DestroyReal1DArray(g2rcc)

      call DestroyReal1DArray(udx1vis1)
      call DestroyReal1DArray(udx1vis2)
      call DestroyReal1DArray(udx1qm)
      call DestroyReal1DArray(udx1qc)
      call DestroyReal1DArray(udx2c)
      call DestroyReal1DArray(udx2m)
      call DestroyReal1DArray(udx2cc)
      call DestroyReal1DArray(udx2mm)
      call DestroyReal1DArray(rudx2c)
      call DestroyReal1DArray(udx2rm)


      call DestroyInt1DArray(kmc)
      call DestroyInt1DArray(kpc)
      call DestroyInt1DArray(kmv)
      call DestroyInt1DArray(kpv)
      call DestroyInt1DArray(kup)
      call DestroyInt1DArray(kum)

      call DestroyReal1DArray(ap1j)
      call DestroyReal1DArray(am1j)
      call DestroyReal1DArray(ac1j)
      call DestroyReal1DArray(ap2j)
      call DestroyReal1DArray(am2j)
      call DestroyReal1DArray(ac2j)
      call DestroyReal1DArray(ap2je)
      call DestroyReal1DArray(am2je)
      call DestroyReal1DArray(ac2je)
      call DestroyReal1DArray(ap3j)
      call DestroyReal1DArray(am3j)
      call DestroyReal1DArray(ac3j)
      call DestroyReal1DArray(apph)
      call DestroyReal1DArray(amph)
      call DestroyReal1DArray(acph)
      call DestroyReal1DArray(ak1)
      call DestroyReal1DArray(ap)
      call DestroyReal1DArray(ak3)
      call DestroyReal1DArray(ao)

      call DestroyInt2DArray(icylsl)
      call DestroyInt2DArray(ocylsl)

      call DestroyReal2DArray(vt_me)
      call DestroyReal2DArray(vr_me)
      call DestroyReal2DArray(vz_me)
      call DestroyReal2DArray(wt_me)
      call DestroyReal2DArray(wr_me)
      call DestroyReal2DArray(wz_me)
      call DestroyReal2DArray(pr_me)
      call DestroyReal2DArray(vt_rms)
      call DestroyReal2DArray(vr_rms)
      call DestroyReal2DArray(vz_rms)
      call DestroyReal2DArray(wt_rms)
      call DestroyReal2DArray(wr_rms)
      call DestroyReal2DArray(wz_rms)
      call DestroyReal2DArray(pr_rms)
      call DestroyReal2DArray(vtvr_me)
      call DestroyReal2DArray(vtvz_me)
      call DestroyReal2DArray(vrvz_me)
      call DestroyReal2DArray(disste)
          
      end subroutine DeallocateVariables

