!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SlipBCs.F90                                    !
!    CONTAINS: subroutine SlipBCs                         !
!                                                         ! 
!    PURPOSE: Calculate slip/no-slip patterns on          !
!     cylinders from inputs in bou.in                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine SlipBCs
      use param
      implicit none
      integer :: jc,ic
      integer :: ntott, ntotz, ntot, acute, nzpoint
      integer :: i, j, k, ntpoint
      integer :: limit1, limit2, points_moved,truej
      integer, allocatable, dimension(:) :: thdispl
      real  :: th, d, m, m1, y1, x1
      real :: anglemoved, perdist, wint
      
!     Determines whether cylinder is free-slip or no-slip at a certain (jc,ic)
!     0 = free-slip, 1 = no-slip
!
!     Cannot set entire cylinder to be free-slip or there will be no forcing
!
      allocate(thdispl(1:nthm))

      acute =1
      icylsl=1 
      ocylsl=1 

      ntott = minstrt*multmin
      ntotz = minstrz*multmin
      ntot= ntott+ntotz

!!    Angular velocity: wint = virotint / rint

      perdist=2*pi/lamb

      wint = virotint / rint
      anglemoved = wint * time

      anglemoved=mod(anglemoved,perdist)
      points_moved=floor(anglemoved*dxt)

      do i=1,nthm
       thdispl(i)=mod(i+points_moved,nthm)
       if(thdispl(i).eq.0) thdispl(i)=nthm
      end do

!      write(*,*) points_moved,thdispl(1)

      if(ntot.ne.0) then
      icylsl=0

      ntpoint= nthm/ntott
      nzpoint= nzm/ntotz

      m1= float(minstrt)/float(minstrz)*float(nzm)/float(nthm)
      y1= float(nzm)/float(ntotz)
      x1= float(nthm)/float(ntott)

      if ((acute.eq.1).and.(ntotz.ne.0)) then
              limit1= -ntott
              limit2= ntotz
      else
              limit1= 0
              limit2= ntot
      endif

      do i= limit1, limit2
       do j=1,nthm
        truej=thdispl(j)
         if((ntotz.eq.0).and.(j-i*x1>0.00001).and. &
             (j-i*x1-x1/2.<0.0001)) then
           icylsl(:, truej)=1 
         else if (ntotz.ne.0) then
        do k= 1,nzm
         if((k-acute*m1*j-i*y1>0.00001).and. &
             (k-acute*m1*j-i*y1-y1/2.<0.00001)) then
          icylsl(k, truej)=1
         endif 
        enddo
        endif
       enddo
      enddo
      
       end if

!       write(*,*) 'matrix'
!       do j=1,nzm
!       write(*,*) icylsl(j,:)
!       end do

      return
      end
