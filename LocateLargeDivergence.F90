!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: LocateLargeDivergence.F90                      !
!    CONTAINS: subroutine LocateLargeDivergence           !
!                                                         ! 
!    PURPOSE: Debugging routine. Output the location(s)   !
!     of excessive divergence.                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine LocateLargeDivergence
      use param
      use local_arrays, only: qr,qz,qt
      use mpih
      use decomp_2d, only: xstart,xend,nrank
      implicit none
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap
        
      if(nrank.eq.0) write(*,*) "I   J   K   RANK"
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,nrm
            kp=kc+1
             dqcap=((qt(kc,jc,ip)-qt(kc,jc,ic))*dxt &
                    +(qr(kp,jc,ic)-qr(kc,jc,ic))*udx2m(kc))*usrm(kc) &
                    +(qz(kc,jp,ic)-qz(kc,jc,ic))*dxz
          if (abs(dqcap).gt.resid) then
          write(*,*) ic,jc,kc,nrank,abs(dqcap)
          write(*,*) "qt",(qt(kc,jc,ip)-qt(kc,jc,ic))*dxt
          write(*,*) "qr",(qr(kp,jc,ic)-qr(kc,jc,ic))*udx2m(kc)*usrm(kc)
          write(*,*) "qz",(qz(kc,jp,ic)-qz(kc,jc,ic))*dxz
              endif
      enddo
      enddo
      enddo
      
      return     
      end         
