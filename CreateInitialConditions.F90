!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity                             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateInitialConditions
      use param
      use local_arrays, only: qr,qz,qt
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: j,k,i
      real :: eps,xxx,yyy
      real :: eta, Alam, Blam,lamsol

      eta = rint / rext
      Alam=(virotext-eta*eta*virotint)/(1-eta*eta)
      Blam=(virotint-virotext)*rint*rint/(1-eta*eta)


      eps=0.1d0
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,nrm
          lamsol=Alam*rm(k)+Blam/rm(k)
          qt(k,j,i)=0.0d0
          yyy=zm(j)/alx3
          qr(k,j,i)=eps*(sin(yyy*2*pi)+sin(yyy*4*pi) &
                        +sin(yyy*8*pi)+sin(yyy*6*pi))
          xxx=(rm(k)-rint)/(rext-rint)
          qz(k,j,i)=(sin(xxx*2*pi)*(1-abs(sin(yyy*2*pi))) &
                    *sin(yyy*2*pi)*eps) +  &
                    (sin(xxx*2*pi)*(1-abs(sin(yyy*4*pi))) &
                    *sin(yyy*4*pi)*eps) +  &
                    (sin(xxx*2*pi)*(1-abs(sin(yyy*6*pi))) &
                    *sin(yyy*6*pi)*eps) +  &
                    (sin(xxx*2*pi)*(1-abs(sin(yyy*8*pi))) &
                    *sin(yyy*8*pi)*eps) 
         enddo
        enddo
      enddo

      return                                                            
      end                                                               
