!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: WriteFlowField.F90                             !
!    CONTAINS: subroutine WriteFlowField                  !
!                                                         ! 
!    PURPOSE: Write down the full flow snapshot for       !
!     restarting the simulation at a later date           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine WriteFlowField
      use param
      use local_arrays, only: qt,qr,qz
      implicit none
      character*30 :: filnam1,dsetname

      filnam1 = trim('continua_qt.h5')
      call HdfWriteRealHalo3D(filnam1,qt)
      filnam1 = trim('continua_qr.h5')
      call HdfWriteRealHalo3D(filnam1,qr)
      filnam1 = trim('continua_qz.h5')
      call HdfWriteRealHalo3D(filnam1,qz)

      if (ismaster) then !EP only write once
       filnam1 = trim('continua_master.h5')
       call HdfCreateBlankFile(filnam1)
 
       dsetname = trim('nr')
       call HdfSerialWriteIntScalar(dsetname,filnam1,nr)
       dsetname = trim('nth')
       call HdfSerialWriteIntScalar(dsetname,filnam1,nth)
       dsetname = trim('nz')
       call HdfSerialWriteIntScalar(dsetname,filnam1,nz)
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

      call MpiBarrier

      end subroutine WriteFlowField
