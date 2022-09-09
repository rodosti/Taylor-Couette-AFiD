!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcMaxCFL.F90                                 !
!    CONTAINS: subroutine CalcMaxCFL                      !
!                                                         ! 
!    PURPOSE: Calculates the maximum value of U/dx to     !
!     determine the new value of the time-step            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine CalcMaxCFL(cflm)
      use param
      use local_arrays, only: qr,qz,qt
      use decomp_2d
      use mpih
      implicit none
      real,intent(inout)    :: cflm
      integer :: j,k,jp,kp,i,ip
      real :: qcf
      
      cflm=0.00000001d0
                                                                        
!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,nrm,qt,qr,qz) &
!$OMP& SHARED(dxt,dxz,usrm,udx2rm) &
!$OMP& PRIVATE(i,j,k,ip,jp,kp,qcf) &
!$OMP& REDUCTION(max:cflm)
      do i=xstart(3),xend(3)
         ip=i+1
        do j=xstart(2),xend(2)
          jp=j+1
          do k=1,nrm
           kp=k+1
                qcf=(abs((qt(k,j,ip)+qt(k,j,i))*0.5*dxt*usrm(k)) &
                     +abs((qr(kp,j,i)+qr(k,j,i))*0.5*udx2rm(k)) &
                     +abs((qz(k,jp,i)+qz(k,j,i))*0.5*dxz))
                cflm = max(cflm,qcf)

      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
            
      call MpiAllMaxRealScalar(cflm)

      return  
      end                                                               
