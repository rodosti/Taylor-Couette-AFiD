!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: AuxiliaryRoutines.F90                          !
!    CONTAINS: subroutines Allocate*,Destroy*             !
!                                                         ! 
!    PURPOSE: Auxiliary routines used for memory allocs   !
!     and memory freeing                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module AuxiliaryRoutines
      contains 

       subroutine AllocateReal1DArray(var,st1,en1)
       use decomp_2d
       implicit none
       integer, intent(in) :: st1,en1
       real,allocatable,dimension(:),intent(inout) :: var
       integer :: alloc_stat, errorcode
  
       if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)
      
       if (alloc_stat /= 0) then
          call decomp_2d_abort(errorcode, &
               'Memory allocation failed when creating new arrays')
       end if
  
       var =0.0d0
  
       return
       end subroutine AllocateReal1DArray
      
!===========================================================================

       subroutine AllocateInt1DArray(var,st1,en1)
       use decomp_2d
       implicit none
       integer, intent(in) :: st1,en1
       integer,allocatable,dimension(:),intent(inout) :: var
       integer :: alloc_stat, errorcode
  
       if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)
      
       if (alloc_stat /= 0) then
          errorcode = 8
          call decomp_2d_abort(errorcode, &
               'Memory allocation failed when creating new arrays')
       end if
  
       var =0
  
       return
       end subroutine AllocateInt1DArray
  
!===========================================================================

       subroutine AllocateReal2DArray(var,st1,en1,st2,en2)
       use decomp_2d
       implicit none
       integer, intent(in) :: st1,en1,st2,en2
       real,allocatable,dimension(:,:),intent(inout) :: var
       integer :: alloc_stat, errorcode
  
       if (.not. allocated(var)) then 
         allocate(var(st1:en1,st2:en2), stat=alloc_stat)
       end if
      
       if (alloc_stat /= 0) then
          errorcode = 8
          call decomp_2d_abort(errorcode, &
               'Memory allocation failed when creating new arrays')
       end if
  
       var =0.0
  
       return
       end subroutine AllocateReal2DArray

!===========================================================================

       subroutine AllocateInt2DArray(var,st1,en1,st2,en2)
       use decomp_2d
       implicit none
       integer, intent(in) :: st1,en1,st2,en2
       integer,allocatable,dimension(:,:),intent(inout) :: var
       integer :: alloc_stat, errorcode
  
       if (.not. allocated(var)) then 
         allocate(var(st1:en1,st2:en2), stat=alloc_stat)
       end if
      
       if (alloc_stat /= 0) then
          errorcode = 8
          call decomp_2d_abort(errorcode, &
               'Memory allocation failed when creating new arrays')
       end if
  
       var =0
  
       return
       end subroutine AllocateInt2DArray

!===========================================================================

        subroutine AllocateReal3DArray(var,st1,en1,st2,en2,st3,en3)
        use decomp_2d
        implicit none
        integer, intent(in) :: st1,en1,st2,en2,st3,en3
        real,allocatable,dimension(:,:,:),intent(inout) :: var
        integer :: alloc_stat, errorcode
  
        if (.not. allocated(var)) then 
          allocate(var(st1:en1,st2:en2,st3:en3), stat=alloc_stat)
        end if
      
        if (alloc_stat /= 0) then
           errorcode = 8
           call decomp_2d_abort(errorcode, &
                'Memory allocation failed when creating new arrays')
        end if
  
        var =0.0
  
        return
        end subroutine AllocateReal3DArray

!===========================================================================

        subroutine DestroyReal1DArray(var)
        use decomp_2d
        implicit none
        real,allocatable,dimension(:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal1DArray

!===========================================================================

        subroutine DestroyInt1DArray(var)
        use decomp_2d
        implicit none
        integer,allocatable,dimension(:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyInt1DArray

!===========================================================================

        subroutine DestroyInt2DArray(var)
        use decomp_2d
        implicit none
        integer,allocatable,dimension(:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyInt2DArray

!===========================================================================

        subroutine DestroyReal2DArray(var)
        use decomp_2d
        implicit none
        real,allocatable,dimension(:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal2DArray

!===========================================================================

        subroutine DestroyReal3DArray(var)
        use decomp_2d
        implicit none
        real,allocatable,dimension(:,:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal3DArray
      
      end module AuxiliaryRoutines
