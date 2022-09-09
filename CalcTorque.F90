!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcTorque.F90                                 !
!    CONTAINS: subroutine CalcTorque                      !
!                                                         ! 
!    PURPOSE: Calculates the torque at both inner and     !
!     outer cylinders and writes it in torquenu.out       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       

      subroutine CalcTorque
      use param
      use local_arrays, only: qt
      use decomp_2d, only: xstart,xend, nrank, DECOMP_2D_COMM_CART_X
      use mpih
 
      IMPLICIT none
      integer i,j
 
      real reninv, my_surface
      real anusslow, grtlow
      real anussup, grtup
      real del, deln,jlam, jlami
 
      reninv = 1.0/ren
      anusslow = 0.d0
      del  = -1.0/(rc(2)-rc(1))
      anussup = 0.d0
      deln = -1.0/(rc(nr)-rc(nrm))

      jlam=2.0*rint*rint*rext*rext*(virotint/rint & 
            -virotext/rext)/(rext*rext-rint*rint)*reninv
      jlami=1.0/jlam

      do i=xstart(3),xend(3)
      do j = xstart(2),xend(2)
          grtlow = ((qt(2,j,i)*usrm(2)+qt(1,j,i)*usrm(1))*0.5d0 &
                    -virotint*usrc(1))*del
          anusslow = anusslow + grtlow*icylsl(j,i)
          grtup = (virotext*usrc(nr) & 
       -(qt(nrm-1,j,i)*usrm(nrm-1)+qt(nrm,j,i)*usrm(nrm))*0.5d0)*deln
          anussup = anussup + grtup*ocylsl(j,i)
       end do
      end do

      anusslow = anusslow * rc(2)**3 * reninv * jlami
      anussup = anussup * rc(nrm)**3 * reninv * jlami
 
      call MpiSumRealScalar(anusslow)
      call MpiSumRealScalar(anussup)

      my_surface = 1.0/(nthm*nzm)
      anusslow = anusslow * my_surface
      anussup = anussup * my_surface

      if(ismaster) write(97,546) time, anusslow, anussup

 546  format(3(1x,e14.6))

      return
      end


