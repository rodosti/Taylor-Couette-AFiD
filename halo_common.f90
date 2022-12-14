!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

    s1 = size(in,1)
    s2 = size(in,2)
    s3 = size(in,3)

    ! Calculate the starting index and ending index of output
    xs = 1
    xe = s1
    ys = 1 + level
    ye = s2 - level
    zs = 1 + level
    ze = s3  - level

    s2 = s2 - 2*level
    s3 = s3 - 2*level

    ! If needed, define MPI derived data type to pack halo data,
    ! then call MPI send/receive to exchange halo data

       ! *** north/south *** 
 
     tag_s = coord(1)
      if (coord(1)==dims(1)-1 .AND. periodic_y) then
         tag_n = 0
      else
         tag_n = coord(1) + 1
      end if
       icount = s3
       ilength = level * s1
       ijump = s1*(s2+2*level)

     if(tag_n.ne.tag_s) then
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            data_type, halo12, ierror)
       call MPI_TYPE_COMMIT(halo12, ierror)
       ! receive from south
       call MPI_IRECV(in(xs,ys-level,zs), 1, halo12, &
            neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
            requests(1), ierror)
       ! receive from north
       call MPI_IRECV(in(xs,ye+1,zs), 1, halo12, &
            neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
            requests(2), ierror)
       ! send to south
       call MPI_ISSEND(in(xs,ys,zs), 1, halo12, &
            neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
            requests(3), ierror)
       ! send to north
       call MPI_ISSEND(in(xs,ye-(level-1),zs), 1, halo12, &
            neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo12, ierror)
      else
!EP   MEMCPY when already in local memory to avoid indirect indexing
!EP   Performance difference unchecked
      do ilvl=1,level
       in(:,ye+ilvl,:) = in(:,ys+(ilvl-1),:)
       in(:,ys-ilvl,:) = in(:,ye-(ilvl-1),:)
      end do
     endif

       ! *** top/bottom ***
      tag_b = coord(2)
       if (coord(2)==dims(2)-1 .AND. periodic_z) then
          tag_t = 0
       else
         tag_t = coord(2) + 1
       end if

      ys = 1 
      s2 = s2 + 2*level

      icount = (s1 * s2) * level

     if(tag_t.ne.tag_b) then
       ! receive from bottom
       call MPI_IRECV(in(xs,ys,zs-level), icount, data_type, &
            neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
            requests(1), ierror)
       ! receive from top
       call MPI_IRECV(in(xs,ys,ze+1), icount, data_type, &
            neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
            requests(2), ierror)
       ! send to bottom
       call MPI_ISSEND(in(xs,ys,zs), icount, data_type, &
            neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
            requests(3), ierror)
       ! send to top
       call MPI_ISSEND(in(xs,ys,ze-(level-1)), icount, data_type, &
            neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
     else
      do ilvl=1,level
       in(:,:,ze+ilvl) = in(:,:,zs+(ilvl-1))
       in(:,:,zs-ilvl) = in(:,:,ze-(ilvl-1))
      end do
     endif


