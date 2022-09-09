!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none
!==========================================================			
!       read from input file bou.in
!==========================================================
        integer   :: nth, nz, nr
        integer   :: nsst, nread, ntst, ireset
        real      :: tout,tmax
        real      :: alx3,str
        integer   :: istr
        real      :: rint,rext
        real      :: dt,resid,limitCFL
        integer   :: starea
        real      :: tsta, dtmax,cfllim 
        integer   :: nson,idtv
        integer   :: nsect
        integer   :: minstrt, minstrz, multmin
        integer   :: nslabr, dpthcut, dpstave, dpfullf
        real      :: virotext,virotint,lamb
        real      :: ren, romeg, dtmin,limitVel
        real      :: tframe2d,tframe3d,walltimemax
        real      :: omegapert, sizepert
!=================================================
!       end of input file
!=================================================
        real :: time
!******* Grid parameters**************************
        real :: dxz,dxr,dxt
        real :: dxzq,dxrq,dxtq
!        
        real, allocatable, dimension(:) :: tc,tm
        real, allocatable, dimension(:) :: zc,zm
        real, allocatable, dimension(:) :: rc,rm,g2rc,g2rm
        real, allocatable, dimension(:) :: g2rmm,g2rcc
        real, allocatable, dimension(:) :: pvol
!====================================================
!******* QUANTITIES FOR DERIVATIVES******************
        real, allocatable, dimension(:) :: usrm, usrc
        real, allocatable, dimension(:) :: udx1vis1, udx1vis2
        real, allocatable, dimension(:) :: udx1qm, udx1qc
        real, allocatable, dimension(:) :: udx2c, udx2m
        real, allocatable, dimension(:) :: udx2cc, udx2mm
        real, allocatable, dimension(:) :: rudx2c, udx2rm
!==========================================================
!******* Grid indices**************************************
        integer, allocatable, dimension(:) :: kmc,kpc,kmv,kpv,kup,kum
!===========================================================
!******* Metric coefficients *******************************
        real, allocatable, dimension(:) :: ap1j,am1j,ac1j
        real, allocatable, dimension(:) :: ap2j,am2j,ac2j
        real, allocatable, dimension(:) :: ap2je,am2je,ac2je
        real, allocatable, dimension(:) :: ap3j,am3j,ac3j
        real, allocatable, dimension(:) :: apph,amph,acph
        real, dimension(2) :: ap1j_fs,ac1j_fs,am1j_fs
        real, dimension(2) :: ap3j_fs,ac3j_fs,am3j_fs
        integer, allocatable, dimension(:,:)  :: icylsl, ocylsl
        real  :: vinner_t,dvinner_t
!============================================================
!******* Variables for FFTW and Poisson solver****************
        real, allocatable, dimension(:) :: ak1,ap
        real, allocatable, dimension(:) :: ak3,ao
        
!===========================================================
!******* Other variables ***********************************
        integer  :: nthm, nzm, nrm
        real :: pi
        real :: al,ga,ro,aldto
        real :: beta
        real :: qqmax,qqtot
        real :: re
        integer :: ntime,infig
        integer, parameter:: ndv=3
        real, dimension(1:ndv) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        logical :: ismaster = .false.
        logical :: dumpslabr=.false.
        logical :: dumpslabt=.false.
        logical :: dumpstave=.false.
        logical :: dumpfullf=.false.
        logical :: readflow=.false.
        logical :: readstats=.false.
        logical :: resetlogstime=.false.
        logical :: variabletstep=.true.

      end module param
      
!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******

      module local_arrays
      use param
        implicit none
        real, allocatable, dimension(:,:,:) :: qt,qr,qz
        real, allocatable, dimension(:,:,:) :: pr,rhs
        real, allocatable, dimension(:,:,:) :: ru1,ru2,ru3
        real, allocatable, dimension(:,:,:) :: dph,qcap,dq
        real, allocatable, dimension(:,:,:) :: vortt,vortr,vortz
      end module local_arrays

!===============================================================

      module stat_arrays
       implicit none
       real, allocatable, dimension(:,:) :: vt_me,vt_rms 
       real, allocatable, dimension(:,:) :: vr_me,vr_rms 
       real, allocatable, dimension(:,:) :: vz_me,vz_rms 
       real, allocatable, dimension(:,:) :: wt_me,wt_rms 
       real, allocatable, dimension(:,:) :: wr_me,wr_rms 
       real, allocatable, dimension(:,:) :: wz_me,wz_rms 
       real, allocatable, dimension(:,:) :: pr_me,pr_rms
       real, allocatable, dimension(:,:) :: vtvr_me,vtvz_me,vrvz_me
       real, allocatable, dimension(:,:) :: disste
       integer :: nstatsamples
      end module stat_arrays
!=====================================================       
      module stat3_param
        implicit none
        integer :: kslab(1:9)
        real    :: rslab(1:9)
      end module stat3_param

!=====================================================       
      module mpih
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer, parameter :: master=0
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer, parameter :: lvlhalo=1
      end module mpih

      
      module mpi_param
        implicit none
        integer :: istart,iend,jstart,jend, kstart,kend
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 
      end module mpi_param
!====================================================
      module fftw_params
!        use param, only: m2m,m2mh,m1m
        use iso_c_binding

        type, bind(C) :: fftw_iodim
           integer(C_INT) n, is, os
        end type fftw_iodim

        interface
          type(C_PTR) function fftw_plan_guru_dft(rank,dims, &
           howmany_rank,howmany_dims,in,out,sign,flags)  &
           bind(C, name='fftw_plan_guru_dft')
           import
           integer(C_INT), value :: rank
           type(fftw_iodim), dimension(*), intent(in) :: dims
           integer(C_INT), value :: howmany_rank
           type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
           complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
           complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
           integer(C_INT), value :: sign
           integer(C_INT), value :: flags
         end function fftw_plan_guru_dft

           type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims, &
             howmany_rank,howmany_dims,in,out,flags)  &
             bind(C, name='fftw_plan_guru_dft_r2c')
             import
             integer(C_INT), value :: rank
             type(fftw_iodim), dimension(*), intent(in) :: dims
             integer(C_INT), value :: howmany_rank
             type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
             real(C_DOUBLE), dimension(*), intent(out) :: in
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
             integer(C_INT), value :: flags
           end function fftw_plan_guru_dft_r2c
           
           type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims, &
             howmany_rank,howmany_dims,in,out,flags)  &
             bind(C, name='fftw_plan_guru_dft_c2r')
             import
             integer(C_INT), value :: rank
             type(fftw_iodim), dimension(*), intent(in) :: dims
             integer(C_INT), value :: howmany_rank
             type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
             real(C_DOUBLE), dimension(*), intent(out) :: out
             integer(C_INT), value :: flags
           end function fftw_plan_guru_dft_c2r


       end interface

        integer FFTW_PATIENT, FFTW_FORWARD, FFTW_BACKWARD
        parameter (FFTW_PATIENT=32)   
        parameter (FFTW_FORWARD=-1)   
        parameter (FFTW_BACKWARD=1)   
        type(C_PTR) :: fwd_guruplan_y,bwd_guruplan_y 
        type(C_PTR) :: fwd_guruplan_z,bwd_guruplan_z
        logical :: planned=.false.

      end module fftw_params
