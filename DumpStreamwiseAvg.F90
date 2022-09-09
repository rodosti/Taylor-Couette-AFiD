!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DumpStreamwiseAvg.F90                          !
!    CONTAINS: subroutine DumpStreamwiseAvg               !
!                                                         ! 
!    PURPOSE: Write down streamwise averages of the       !
!     velocity fields.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DumpStreamwiseAvg
      use mpih
      use param
      use local_arrays, only: qt,qr,qz
      use decomp_2d, only: xstart,xend,nrank
      implicit none
      integer :: k,j,itime,i
      real :: vtcc(nrm,nzm)
      real :: vrcc(nrm,nzm)
      real :: vzcc(nrm,nzm)
      real :: usnthm
      character*30 :: filnamvt
      character*30 :: filnamvr
      character*30 :: filnamvz
      character*30 :: dsetnam
      character*8 ipfi
      character*70 :: filnamxmf

      vtcc = 0.0d0
      vrcc = 0.0d0
      vzcc = 0.0d0
      usnthm = 1.0 / float(nthm)

       do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
         do k=1,nrm
           vtcc(k,j) = vtcc(k,j)+qt(k,j,i)*usnthm
           vrcc(k,j) = vrcc(k,j)+qr(k,j,i)*usrc(k)*usnthm
           vzcc(k,j) = vzcc(k,j)+qz(k,j,i)*usnthm
          enddo
         enddo
        end do

      call MpiSumReal1D(vtcc,nrm*nzm)
      call MpiSumReal1D(vrcc,nrm*nzm)
      call MpiSumReal1D(vzcc,nrm*nzm)

      if (ismaster) then

      itime=nint(1000.*time+0.005)
      write(ipfi,99) itime
   99 format(i8.8)

      filnamvt='stave/slabvt_ax_'//ipfi//'.h5'
      filnamvr='stave/slabvr_ax_'//ipfi//'.h5'
      filnamvz='stave/slabvz_ax_'//ipfi//'.h5'

      dsetnam='vt'
      call HdfWriteSerialReal2D(filnamvt,dsetnam,nrm,nzm,vtcc)
      dsetnam='vr'
      call HdfWriteSerialReal2D(filnamvr,dsetnam,nrm,nzm,vrcc)
      dsetnam='vz'
      call HdfWriteSerialReal2D(filnamvz,dsetnam,nrm,nzm,vzcc)

      endif


      return
      end subroutine DumpStreamwiseAvg

      
