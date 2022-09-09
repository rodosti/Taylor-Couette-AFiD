!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitPressureSolver.F90                         !
!    CONTAINS: subroutine InitPressureSolver              !
!                                                         ! 
!    PURPOSE: Initialization routines. Compute the metric !
!     terms and modified wavenumbers for the pressure     !
!     correction                                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine InitPressureSolver
      use param
      use decomp_2d_fft
      implicit none
      integer  :: kc,km,kp
      integer :: nzmh,nzmp,j,i,nthmh,nthmp
      real :: a22icc, a22icp, ac2, ugmmm
    
!RO   Initialize tridiag matrices
      do kc=1,nrm               
       km=kmv(kc)
       kp=kpv(kc)
       a22icc=rc(kc)*kmc(kc)*dxrq/g2rc(kc)
       a22icp=rc(kp)*kpc(kc)*dxrq/g2rc(kp)
       ac2=-(a22icc+a22icp)
       ugmmm=1./(rm(kc)*g2rm(kc))
       amph(kc)=a22icc*ugmmm
       apph(kc)=a22icp*ugmmm  
       acph(kc)=ac2*ugmmm      
      end do

!RO   Initialize FFTW
      nthmh=nthm/2+1
      nthmp=nthmh+1
      nzmh=nzm/2+1
      nzmp=nzmh+1


!RO    Modified wave number definition axial

      do j=1,nzmh                                         
       ao(j)=(j-1)*2.*pi                                     
      end do

      do j=nzmp,nzm                                     
       ao(j)=-(nzm-j+1)*2.*pi                              
      end do

      do j=1,nzm                                    
       ak3(j)=2.*(1.-cos(ao(j)/nzm))*(float(nzm)/alx3)**2 
      end do

!RO    Modified wave number definition azimuthal

      do i=1,nthmh  
       ap(i)=(i-1)*2.d0*pi
      end do
      do i=nthmp,nthm  
       ap(i)=-(nthm-i+1)*2.d0*pi 
      end do

      do i=1,nthm                       
       ak1(i)=2.*(1.-cos(ap(i)/nthm))*dxtq
      end do

!EP   Planning
      call decomp_2d_fft_init

      return
      end
      
