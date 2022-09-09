!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcDissipation.F90                            !
!    CONTAINS: subroutine CalcDissipation                 !
!                                                         ! 
!    PURPOSE: Calculates the instantaneous kinetic        !
!     energy dissipation and writes it in dissipnu.out    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcDissipation
      use mpih
      use param
      use local_arrays,only: qt,qr,qz
      use decomp_2d, only: xstart,xend,nrank
      use stat_arrays

      implicit none
      integer :: i,j,k
      integer :: imm,ipp,jmm,jpp,kmm,kpp,kp
      real :: udx3_m,udx3_c
      real :: e11,e12,e13,e22,e23,e33
      real :: nute,volt, lamdissip
      real :: dissipte,my_dissipte

      
      my_dissipte = 0.0d0

!
!                      1  |         | 2
!      dissipation:  ---- | nabla  u|
!                     Re  |         |
!

!$OMP  PARALLEL DO &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(xstart,xend,qt,qr,qz) &
!$OMP& SHARED(usrc,udx2mm,g2rm) &
!$OMP& SHARED(nthm,nzm,nrm,ren,pi,pvol) &
!$OMP& SHARED(kpv,kmv,rm,dxz,udx2cc) &
!$OMP& SHARED(disste,usrm,dxt) &
!$OMP& PRIVATE(i,j,k,imm,ipp,jmm,jpp,kp,kpp,kmm) &
!$OMP& PRIVATE(udx3_m,udx3_c,dissipte) &
!$OMP& PRIVATE(e11,e12,e13,e22,e23,e33) &
!$OMP& REDUCTION(+:my_dissipte)
      do i=xstart(3),xend(3)
       imm= i-1
       ipp= i+1
       do j=xstart(2),xend(2)
        jmm=j-1
        jpp=j+1
        do k=1,nrm
         kp=k+1
         kpp=kpv(k)
         kmm=kmv(k)
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           1  d u_t     u_r
!           -  ----- +  ----                    e11
!           r   d_t       r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          e11=(qt(k,j,ipp)-qt(k,j,imm))*dxt*usrm(k)*0.5 & 
           + qr(k,j,i)*usrc(k)**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 d u_t / r     1  d u_r
!              r    -----    +  -  -----                e12
!                   d_r         r   dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          e12=(qt(kpp,j,i)*usrm(kpp)-qt(kmm,j,i)*usrm(kmm)) & 
                     *udx2mm(k)*rm(k) & 
             + (qr(k,j,ipp)-qr(k,j,imm))*usrc(k)*dxt*0.5  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               d u_t        1  d u_z
!                -----    +  -  -----                   e13
!                d_z         r   dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


           e13=(qt(k,jpp,i)-qt(k,jmm,i))*0.5*dxz &
             + (qz(k,j,ipp)-qz(k,j,imm))*0.5*dxt*usrm(k)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             d u_r 
!             -----                                     e22
!              d_r 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          e22=(qr(kp,j,i)*usrc(kp)-qr(kmm,j,i)*usrc(kmm))*udx2cc(k)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           d u_z         d u_r
!           -----    +    -----                          e23
!            d_r           d_z 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          e23=(qz(kpp,j,i)-qz(kmm,j,i))*udx2mm(k)  & 
             + (qr(k,jpp,i)-qr(k,jmm,i))*usrc(k)*dxz*0.5


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           d u_z  
!           -----                                       e33
!            d_z   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        e33=(qz(k,jpp,i)-qz(k,jmm,i))*dxz*0.5



       dissipte = 2.0*(e11**2+e22**2+e33**2)+  &
               e12**2+e13**2+e23**2


       my_dissipte = my_dissipte+dissipte*pvol(k)*pi
 
!$OMP CRITICAL
       disste(k,j) =  disste(k,j) + dissipte / ren / real(nthm) 
!$OMP END CRITICAL

       end do
       end do
       end do
!$OMP  END PARALLEL DO


      call MpiSumRealScalar(my_dissipte)

      volt = pi*(rext**2-rint**2)

      lamdissip =  4.0*rext**2*rint**2/((rext+rint)**2)* & 
               (((virotint/rint-virotext/rext)/ & 
                 (rext-rint))**2)

      nute = my_dissipte / (float(nthm)*float(nzm)*volt*lamdissip)

      if(ismaster) write(92,*) time,nute

      return   
      end
