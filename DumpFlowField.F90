!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DumpFlowField.F90                              !
!    CONTAINS: subroutine DumpFlowField                   !
!                                                         ! 
!    PURPOSE: Write down the full flow snapshot           !
!     statistical use. Very similar to                    !
!      WriteFlowField.F90                                 !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DumpFlowField
      use param
      use local_arrays, only: qt,qr,qz
      implicit none
      integer :: itime
      character*30 :: filnam1,dsetname
      character*8 ipfi

      itime=nint(1000.*time+0.005)
      write(ipfi,99) itime
   99 format(i8.8)


      filnam1 = trim('fullff/qt_'//trim(ipfi)//'.h5')
      call HdfWriteRealHalo3D(filnam1,qt)
      filnam1 = trim('fullff/qr_'//trim(ipfi)//'.h5')
      call HdfWriteRealHalo3D(filnam1,qr)
      filnam1 = trim('fullff/qz_'//trim(ipfi)//'.h5')
      call HdfWriteRealHalo3D(filnam1,qz)

      if (ismaster) then !EP only write once
       filnam1 = trim('fullff/ms_'//trim(ipfi)//'.h5')
       call HdfCreateBlankFile(filnam1)
       dsetname = trim('time')
       call HdfSerialWriteRealScalar(dsetname,filnam1,time)
       dsetname = trim('rint')
       call HdfSerialWriteRealScalar(dsetname,filnam1,rint)
       dsetname = trim('rext')
       call HdfSerialWriteRealScalar(dsetname,filnam1,rext)
       dsetname = trim('istr')
       call HdfSerialWriteIntScalar(dsetname,filnam1,istr)
       dsetname = trim('str')
       call HdfSerialWriteRealScalar(dsetname,filnam1,str)
       dsetname = trim('alx3')
       call HdfSerialWriteRealScalar(dsetname,filnam1,alx3)
       dsetname = trim('lamb')
       call HdfSerialWriteRealScalar(dsetname,filnam1,lamb)

      endif

      end subroutine DumpFlowField
