!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DumpThCuts.F90                                 !
!    CONTAINS: subroutine DumpThCuts                      !
!                                                         ! 
!    PURPOSE: Write down a constant azimuth cut of the    !
!     velocity fields.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DumpThCuts
      use mpih
      use param
      use local_arrays, only: qt,qr,qz
      use decomp_2d, only: xstart,xend,nrank
      implicit none
      integer :: k,j,itime,i
      real :: vtcc(nrm,nzm)
      real :: vrcc(nrm,nzm)
      real :: vzcc(nrm,nzm)
      character*30 :: filnamvt
      character*30 :: filnamvr
      character*30 :: filnamvz
      character*30 :: dsetnam
      character*8 ipfi
      character*70 :: filnamful

      vtcc = 0.0d0
      vrcc = 0.0d0
      vzcc = 0.0d0

       do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
         do k=1,nrm
           if(i.eq.1) then
            vtcc(k,j) = vtcc(k,j)+qt(k,j,i)
            vrcc(k,j) = vrcc(k,j)+qr(k,j,i)*usrc(k)
            vzcc(k,j) = vzcc(k,j)+qz(k,j,i)
           end if
          enddo
         enddo
        end do

      call MpiAllSumReal1D(vtcc,nrm*nzm)
      call MpiAllSumReal1D(vrcc,nrm*nzm)
      call MpiAllSumReal1D(vzcc,nrm*nzm)

      if (ismaster) then

      itime=nint(1000.*time+0.005)
      write(ipfi,99) itime
   99 format(i8.8)


       filnamvt='thcut/slabvt_'//ipfi//'.h5'
       dsetnam=trim('vt')
       call HdfWriteSerialReal2D(filnamvt,dsetnam,nrm,nzm,vtcc)

       filnamvr='thcut/slabvt_'//ipfi//'.h5'
       dsetnam=trim('vr')
       call HdfWriteSerialReal2D(filnamvr,dsetnam,nrm,nzm,vrcc)

       filnamvz='thcut/slabvz_'//ipfi//'.h5'
       dsetnam=trim('vz')
       call HdfWriteSerialReal2D(filnamvz,dsetnam,nrm,nzm,vzcc)

      endif

      return
      end subroutine DumpThCuts

