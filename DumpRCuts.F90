      subroutine DumpRCuts
      use param
      use local_arrays, only: qt,qr,qz
      use stat3_param
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: i,j,m,itime
      real :: vtcc(xstart(2):xend(2),xstart(3):xend(3))
      real :: vzcc(xstart(2):xend(2),xstart(3):xend(3))
      real :: vrcc(xstart(2):xend(2),xstart(3):xend(3))
      character*70 :: filnam, dsetnam
      character*1 :: charm
      character*8 :: ipfi

!RO   File writing part

      itime=nint(1000.*time+0.005)
      write(ipfi,99) itime
   99 format(i8.8)

!EP   Slabs
!EP   cell center only qr

      write(*,*) 'nslabr',nslabr
      do m=1,nslabr
!$OMP  PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
        do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           vzcc(j,i) = qz(kslab(m),j,i)
           vtcc(j,i) = qt(kslab(m),j,i)
           vrcc(j,i) = (qr(kslab(m),j,i)*usrc(kslab(m))+ &
                        qr(kslab(m)+1,j,i)*usrc(kslab(m)+1))*0.5
          enddo
         enddo
!$OMP  END PARALLEL DO
      write(charm,28) m
   28 format(i1.1)
      filnam=trim('./rcuts/slab'//charm//'vt_'//trim(ipfi)//'.h5')
      dsetnam=trim('vt')
      call HdfWriteReal2D(dsetnam,filnam,vtcc)

      dsetnam=trim('vr')
      filnam=trim('./rcuts/slab'//charm//'vr_'//trim(ipfi)//'.h5')
      call HdfWriteReal2D(dsetnam,filnam,vrcc)

      dsetnam=trim('vz')
      filnam=trim('./rcuts/slab'//charm//'vz_'//trim(ipfi)//'.h5')
      call HdfWriteReal2D(dsetnam,filnam,vzcc)

      enddo

      return
      end subroutine DumpRCuts

!==================================================

      subroutine InitDumpRCuts
      use param
      use decomp_2d, only: xstart,xend
      use stat3_param
      implicit none
      integer :: i,k,j
      character(len=4) :: dummy

!EP   Read from stst3.in
      
      if (nslabr.gt.9) nslabr=9
      open(unit=19,file='rcuts.in',status='old')
        read(19,301) dummy
        read(19,*) (rslab(i),i=1,nslabr)
301     format(a4)                
      close(19)


!EP   Compute which kslab corresponds to which rslab
       kslab=1

        do k=2,nrm
          do j=1,nslabr
            if(rm(k).gt.rslab(j).and.rm(k-1).lt.rslab(j)) then
             kslab(j) = k
            endif
          enddo
        end do
      

!EP   Write probe and slab locations
      
      open(unit=23,file='rcuts/rcutslocs.out',status='unknown')
        rewind(23)
        write(23,*) (kslab(i),i=1,nslabr)
      close(23)

      return
      end subroutine InitDumpRCuts

!==================================================
      
      subroutine WriteRCut(var,filnam)
      USE param
      use mpih
      USE hdf5
      use decomp_2d, only: xstart,xend,nrank

      IMPLICIT none

      real, intent(in) :: var(xstart(2):xend(2),xstart(3):xend(3))

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: timespace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset

      integer(HSIZE_T) :: dims(2)
      integer(HSIZE_T) :: dims2(1)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims,ndims2

      real :: tprfi
      integer :: itime

      character*70,intent(in) :: filnam
      character*70 :: namfile
      character*8 :: ipfi

!RO   File writing part

      itime=nint(1000.*time+0.005)
      write(ipfi,99) itime
   99 format(i8.8)

      namfile=trim('./rcuts/'//trim(filnam)//trim(ipfi)//'.h5')

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts
   
      ndims = 2

      dims(1)=nzm
      dims(2)=nthm

      ndims2 = 1
      dims2(1)=1

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

      data_count(1) = xend(2)-xstart(2)+1
      data_count(2) = xend(3)-xstart(3)+1

      data_offset(1) = xstart(2)-1
      data_offset(2) = xstart(3)-1


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)
 
      call h5pset_fapl_mpio_f(plist_id, comm, info, &
        hdf_error)

      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'var', H5T_NATIVE_DOUBLE, &
                      filespace, dset, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                           data_offset, data_count, hdf_error)

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                             hdf_error)

       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
        var(xstart(2):xend(2),xstart(3):xend(3)), dims,  &
        hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
        xfer_prp = plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   Write time
      if(nrank.eq.0) then
      call h5fopen_f(namfile, H5F_ACC_RDWR_F, file_id, hdf_error)
      call h5screate_simple_f(ndims2, dims2, timespace, hdf_error) 
      call h5dcreate_f(file_id, 'time', H5T_NATIVE_DOUBLE, &
                      timespace, dset, hdf_error)

      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, time, &
             dims2,hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      endif
      return                                                          
      end subroutine WriteRCut
