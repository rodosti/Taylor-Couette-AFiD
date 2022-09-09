!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         ! 
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateGrid
      use param
      use mpi
      implicit none
      integer  :: j,k,nrmo,nclip,i,kc,km,kp
      real, dimension(1:nr) :: etaz
      real, dimension(1:nr+500) :: etazm
      real :: tstr, z2dp
      real :: x3,etain,delet
      real:: a22,a22m,a22p
      real :: apjjc, amjjc
      real :: urpp, urpc, ugmm

      pi=2.d0*asin(1.d0)

!RO   Metric terms for all directions 

      dxt=(float(nthm)*lamb)/(2.*pi)
      dxz=float(nzm)/alx3
      dxr=float(nrm)/1.0d0

      dxtq=dxt*dxt
      dxzq=dxz*dxz
      dxrq=dxr*dxr


!
!     prepare indices for direction normal to non-slip walls
!
      do kc=1,nrm
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.nrm) kpv(kc)=kc
      end do

      do kc=1,nrm
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
      enddo

!
!     AZIMUTHAL COORDINATE DEFINITION (uniform) thetac: value at grid, thetam:midpoint
!
      do i=1,nth
        tc(i)= float(i-1)/dxt
      end do
      do i=1,nthm
        tm(i)= (float(i-1)+0.5d0)/dxt
      end do

!
!     AXIAL COORDINATE DEFINITION (uniform)
!

      do j=1,nz
        zc(j)= float(j-1)/dxz
      end do
      do j=1,nzm
        zm(j)= (float(j-1)+0.5d0)/dxz
      end do

!
!     RADIAL COORDINATE DEFINITION
!

!
!     UNIFORM GRID
!

      tstr=tanh(str)

      if (istr.eq.0) then
        do k=1,nr
          x3=real(k-1)/real(nrm)
          rc(k)=(rext-rint)*x3+rint
        enddo
      endif


!       
!      CLUSTERING AT THE EXTERNAL RADIAL WALL 
!                       and  
!             CLUSTERING AT THE AXIS 
!      

        if (istr.eq.4) then
         rc(1)=rint
         do k=2,nr
          z2dp=float(2*k-nr-1)/float(nrm)
          rc(k)=(1+tanh(str*z2dp)/tstr)*0.5*(rext-rint)+rint
         end do
        end if



      if(istr.eq.6) then
      nclip = int(str)
      nrmo = nr+nclip+nclip
      do k=1,nrmo
        etazm(k)=+cos(pi*(float(k)-0.5)/float(nrmo))
      end do
      do k=1,nr
        etaz(k)=etazm(k+nclip)
      end do
      delet = etaz(1)-etaz(nr)
      etain = etaz(1)
      do k=1,nr
        etaz(k)=etaz(k)/(0.5*delet)
      end do
      rc(1) = rint
      do k=2,nrm
        rc(k) = (rext-rint)*(1.-etaz(k))*0.5+rint
      end do
      rc(nr) = rext
      endif
      
!-----------------------------------------
!
!     STAGGERED COORDINATES AND
!     METRIC QUANTITIES
!

      do k=1,nrm
        rm(k)=(rc(k)+rc(k+1))*0.5d0
        g2rm(k)=(rc(k+1)-rc(k))*dxr
      end do

      g2rc(1)=(rc(2)-rc(1))*dxr
      g2rmm(1) =(rm(2)-rm(1))*dxr
      g2rcc(1) =(rc(2)-rc(1))*dxr

      do k=2,nrm
        g2rc(k)=(rc(k+1)-rc(k-1))*dxr*0.5d0
        g2rmm(k) = (rm(k+1)-rm(k-1))*dxr
        g2rcc(k) = (rc(k+1)-rc(k-1))*dxr
      end do

      g2rc(nr)=(rc(nr)-rc(nrm))*dxr
      g2rmm(nr) =(rm(nr)-rm(nrm))*dxr
      g2rcc(nr) =(rc(nr)-rc(nrm))*dxr


      do k=1,nrm
        udx1vis1(k) = 2.d0*dxt/(rm(k)**2*ren)
        udx2rm(k) = dxr/(g2rm(k)*rm(k))
        udx2m(k) = dxr/g2rm(k)
        udx2mm(k) = dxr/g2rmm(k)
        udx2c(k) = dxr/g2rc(k)
        rudx2c(k) = udx2c(k)*rc(k)
        usrm(k) = 1.d0/rm(k)
        udx1qm(k) = dxtq/rm(k)**2
        pvol(k) = (rc(k+1)**2-rc(k)**2)
      end do
      udx2c(nr) = dxr/g2rc(nr)

       do k=1,nr
        udx1vis2(k) = dxt/(rc(k)*ren)
        usrc(k) = 1.d0/rc(k)
        udx1qc(k) = dxtq/rc(k)**2
        udx2cc(k) = dxr/g2rcc(k)
      end do


!
!     WRITE GRID INFORMATION
!
      if(ismaster) then
      open(unit=78,file='radcor.out',status='unknown')
      do k=1,nr
        write(78,345) k,rc(k),rm(k),g2rc(k),g2rm(k)
      end do
      close(78)
 345  format(i4,4(2x,e23.15))

      endif

!
!  ***********  coefficients for q1   inner points
!

      do kc=2,nrm-1
       kp=kc+1
       km=kc-1
       a22=dxrq/g2rm(kc)/rm(kc)
       a22p= +a22*rc(kp)/g2rc(kp)
       a22m= +a22*rc(kc)/g2rc(kc)
       ap1j(kc)=a22p
       am1j(kc)=a22m
       ac1j(kc)=-(a22p+a22m)-1./rm(kc)**2
      end do


!
!    external bound. conditions
!
      kc=nrm
      kp=kc+1
      a22=dxrq/(g2rm(kc)*rm(kc))
      apjjc=a22*rc(kp)/g2rc(kp)
      amjjc=a22*rc(kc)/g2rc(kc)
      ap1j(kc)=2.*apjjc
      ac1j(kc)=-2.*apjjc-amjjc-1./rm(kc)**2
      am1j(kc)=amjjc

      ap1j_fs(1)=0.
      ac1j_fs(1)=-amjjc
      am1j_fs(1)=amjjc


!
!    internal boundary conditions
!
      kc=1
      kp=kc+1
      urpp= 1./(g2rm(kc)*g2rc(kp))
      urpc= 1./(g2rm(kc)*g2rc(kc))
      a22=dxrq/g2rm(kc)/rm(kc)
      a22p= +a22*rc(kp)/g2rc(kp)
      a22m= +a22*rc(kc)/g2rc(kc)
      ap1j(kc)=a22p
      am1j(kc)=2.*a22m
      ac1j(kc)=-(a22p+2.*a22m)-1./rm(kc)**2
      ap1j_fs(1)= rc(kp)**3/(rm(kc)**2*rm(kp))*dxrq*urpp
      ac1j_fs(1)= - rc(kp)**3/rm(kc)**3*dxrq*urpp
      am1j_fs(1)=0.

!
!  ***********  coefficients for q2   
!

      am2j(1)=0.
      ap2j(1)=0.
      ac2j(1)=1.
      do kc=2,nrm
       km=kc-1
       kp=kc+1
       a22=rc(kc)*dxrq/g2rc(kc)
       a22p=1./(rm(kc)*g2rm(kc))
       a22m=1./(rm(km)*g2rm(km))
       ap2j(kc)=a22*a22p
       am2j(kc)=a22*a22m
       ac2j(kc)=-(a22*a22p+a22*a22m)
      end do
      am2j(nr)=0.
      ap2j(nr)=0.
      ac2j(nr)=1.

!
!  ***********  coefficients for q2  (explicit part)
!

!RO   Attention- this was kc=3,nrm in old code ??
      do kc=2,nrm
       km=kc-1
       kp=kc+1
       ugmm = dxrq/g2rc(kc)
       a22p=rm(kc)/(rc(kp)*g2rm(kc))
       a22m=rm(km)/(rc(km)*g2rm(km))
       a22= (rm(kc)/g2rm(kc) + rm(km)/g2rm(km))/rc(kc)
       ap2je(kc)=a22p*ugmm
       am2je(kc)=a22m*ugmm
       ac2je(kc)=-a22*ugmm  - 1./rc(kc)**2
      end do

!
!  ***********  coefficients for q3   inner points
!
      do kc=2,nrm-1
       kp=kc+1
       km=kc-1
       a22=dxrq/g2rm(kc)/rm(kc)
       a22p= +a22*rc(kp)/g2rc(kp)
       a22m= +a22*rc(kc)/g2rc(kc)
       ap3j(kc)=a22p
       am3j(kc)=a22m
       ac3j(kc)=-(a22p+a22m)
      end do

!     
!    internal bound. conditions  q3
!     

      kc=1
      kp=kc+1
      a22=dxrq/g2rm(kc)/rm(kc)
      a22p= +a22*rc(kp)/g2rc(kp)
      a22m= +a22*rc(kc)/g2rc(kc)
      am3j(kc)= 0.
      ap3j(kc)= a22p
      ac3j(kc)= -(a22p+2.*a22m)
      am3j_fs(1)= 0.
      ap3j_fs(1)= a22p
      ac3j_fs(1)= -(a22p)

!     
!    external bound. conditions  q3
!     
      kc=nrm
      kp=kc+1
      km=kc-1
      a22=dxrq/g2rm(kc)/rm(kc)
      a22m=+a22*rc(kc)/g2rc(kc)
      a22p=+a22*rc(kp)/g2rc(kp)
      am3j(kc)= a22m
      ap3j(kc)=0.
      ac3j(kc)=-(a22m +2.*a22p)
      am3j_fs(2)= a22m
      ap3j_fs(2)= 0.
      ac3j_fs(2)= -(a22m)

      call SlipBCs

      return                                                            
      end                                                               
