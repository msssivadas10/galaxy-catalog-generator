module pair_counting_mod
    
    use omp_lib
    use iso_c_binding
    use spatial_hash_mod
    use pair_utils_mod
    implicit none

    private
    public :: count_pairs
    
contains

    subroutine count_pairs(pid, auto, periodic,                 &
                           box, nr, rbins, cnts, ncells,        &
                           npts1, grid_info1, grid_data1, pos1, &
                           npts2, grid_info2, grid_data2, pos2, &
                           nthreads, error_code                 &
        ) bind(c)
        !! Count the number of pairs between two sets of points, given their grid
        !! hash representation.

        integer(c_int64_t), intent(in), value :: pid
        !! Process ID
        
        integer(c_int), intent(in), value :: auto
        !! Flag for taking auto correlation (if both sets are same)
     
        integer(c_int), intent(in), value :: periodic
        !! Flag for using periodic boundary conditions
        
        integer(c_int64_t), intent(in), value :: nr        
        real(c_double)    , intent(in)        :: rbins(nr) 
        !! Distance bin edges
        
        integer(c_int64_t), intent(out) :: cnts(:)
        !! Pair counts array: must have size `nr-1`
        
        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid (same for both sets).
        
        type(boxinfo_t), intent(in) :: box
        !! Details about the space (both sets must be in the same box)
        
        ! -- Points and grid for set-1:
        integer(c_int64_t), intent(in), value :: npts1              !! Number of points (1)    
        real(c_double)    , intent(in)        :: pos1(3,npts1)      !! Position buffer (1)    
        integer(c_int64_t), intent(in)        :: grid_data1(npts1)  !! Grid representation (1)
        type(cinfo_t)     , intent(in)        :: grid_info1(ncells) !! Grid details (1)        
        
        ! -- Points and grid for set-2:
        integer(c_int64_t), intent(in), value :: npts2              !! Number of points (2)      
        real(c_double)    , intent(in)        :: pos2(3,npts2)      !! Position buffer (2)      
        integer(c_int64_t), intent(in)        :: grid_data2(npts2)  !! Grid representation (2)  
        type(cinfo_t)     , intent(in)        :: grid_info2(ncells) !! Grid details (2)          

        integer(c_int), intent(in), value :: nthreads
        !! Number of threads to use
        
        integer(c_int), intent(out) :: error_code
        !! Error code (0=success, 1=error)

        integer(c_int64_t), parameter :: chunk_size = 100

        character(256)     :: ifn
        real(c_double)     :: r2bins(nr)
        integer(c_int)     :: fi, iostat
        integer(c_int64_t) :: file_size_bytes, n_pairs, n_pairs_total, n_pairs_processed, &
                              j1, j2, start1, count1, start2, count2   
        integer(c_int64_t), allocatable :: pairs(:,:)

        error_code = 1

        call omp_set_num_threads(nthreads) ! set number of threads

        ! Precalculating cell pairs
        call enumerate_cell_pairs(pid, npts1, npts2, ncells, rbins(nr), box,  &
                                  periodic, grid_info1, grid_info2, nthreads, &
                                  error_code                                  &
        )
        if ( error_code /= 0 ) return 

        r2bins = rbins**2
        cnts   = 0_c_int64_t

        ! Calculating counts for uniformly populated cells: parallelised 
        ! over cell pairs:
        fi = 9 ! file unit for input
        write(ifn, '(i0,".cpstack.bin")') pid ! filename for input
        open(newunit=fi, file=ifn, access='stream', form='unformatted', &
             convert='little_endian', status='old', action='read'       &
        ) ! input file
        inquire(fi, size=file_size_bytes)

        allocate( pairs(2, chunk_size) )
        
        n_pairs_total = file_size_bytes / (2*c_sizeof(1_c_int64_t))
        if ( n_pairs_total > 0 ) then
            ! Parallel pair counting over cell pairs:
            n_pairs_processed = 0
            do while ( n_pairs_processed < n_pairs_total )

                ! Loading pairs as chunks
                n_pairs = min( chunk_size, n_pairs_total - n_pairs_processed )
                read(fi, iostat=iostat) pairs(1:2, 1:n_pairs)
                if ( iostat /= 0 ) exit ! EOF or error    

                ! Counting pairs
                call count_pairs_cells(pid, auto, periodic, n_pairs, pairs, &
                                       box, nr, r2bins, cnts, ncells,       &
                                       npts1, grid_info1, grid_data1, pos1, &
                                       npts2, grid_info2, grid_data2, pos2  &
                )
                n_pairs_processed = n_pairs_processed + n_pairs
                
            end do
        end if
        close(fi, status='delete') ! file is deleted on close
        deallocate( pairs )

        ! Calculating pair counts for pairs of dense cells. This is parallelised 
        ! over points. 
        fi = 9 ! file unit for input
        write(ifn, '(i0,".cpstack.spl.bin")') pid ! filename for input
        open(newunit=fi, file=ifn, access='stream', form='unformatted', &
             convert='little_endian', status='old', action='read'       &
        ) ! input file
        inquire(fi, size=file_size_bytes)
        n_pairs_total = file_size_bytes / (2*c_sizeof(1_c_int64_t))
        if ( n_pairs_total > 0 ) then
            ! Parallel pair counting over cell pairs:
            do

                read(fi, iostat=iostat) j1, j2
                if ( iostat /= 0 ) exit ! EOF or error    

                start1 = grid_info1(j1)%start
                count1 = grid_info1(j1)%count
                
                start2 = grid_info2(j2)%start
                count2 = grid_info2(j2)%count
                
                ! Counting pairs
                call count_pairs_parallel(pid, nthreads, periodic,                 &
                                          box%boxsize, nr, r2bins, cnts,           &
                                          npts1, start1, count1, grid_data1, pos1, & 
                                          npts2, start2, count2, grid_data2, pos2  &
                )
                
            end do
        end if
        close(fi, status='delete') ! file is deleted on close

        error_code = 0
        
    end subroutine count_pairs

    subroutine count_pairs_cells(pid, auto, periodic, npairs, pairs,  &
                                 box, nr, r2bins, cnts, ncells,       &
                                 npts1, grid_info1, grid_data1, pos1, &
                                 npts2, grid_info2, grid_data2, pos2  &
        ) bind(c)
        !! Count pairs with parallelization over cells.

        integer(c_int)    , intent(in)    :: auto, periodic
        integer(c_int64_t), intent(in)    :: pid, nr, ncells, npts1, npts2              
        integer(c_int64_t), intent(in)    :: npairs, pairs(2,npairs)
        integer(c_int64_t), intent(in)    :: grid_data1(npts1), grid_data2(npts2)  
        real(c_double)    , intent(in)    :: r2bins(nr), pos1(3,npts1), pos2(3,npts2)      
        type(cinfo_t)     , intent(in)    :: grid_info1(ncells), grid_info2(ncells) 
        type(boxinfo_t)   , intent(in)    :: box
        integer(c_int64_t), intent(inout) :: cnts(:)
        
    end subroutine count_pairs_cells

    subroutine count_pairs_sequential(periodic, boxsize, nr, r2bins, cnts,     &
                                      npts1, start1, count1, grid_data1, pos1, &
                                      npts2, start2, count2, grid_data2, pos2  &
        )
        !! Counting the pairs between sections of position buffers. This 
        !! is a sequential version of the broute force pair counting between 
        !! two cells.
        
        integer(c_int)    , intent(in)    :: periodic
        integer(c_int64_t), intent(in)    :: npts1, npts2, nr
        integer(c_int64_t), intent(in)    :: start1, start2, count1, count2
        integer(c_int64_t), intent(in)    :: grid_data1(npts1), grid_data2(npts2)
        real(c_double)    , intent(in)    :: pos1(3,npts1), pos2(3,npts2), boxsize(3)
        real(c_double)    , intent(in)    :: r2bins(nr)
        integer(c_int64_t), intent(inout) :: cnts(:) ! Must be of size nr-1

        integer(c_int64_t) :: j1, j2, i1, i2, bin_loc
        real(c_double)     :: dx(3), r2

        do j1 = start1, start1+count1-1

            ! Index of the point in position buffer
            i1 = grid_data1(j1) 
            
            do j2 = start2, start2+count2-1
                
                ! Index of the point in position buffer
                i2 = grid_data2(j2) 
                
                ! Check pair distance
                dx = pos2(1:3, i2) - pos1(1:3, i1) 
                if ( periodic /= 0 ) then
                    ! Correcting the coordinate distances for periodic boundary
                    dx = dx - boxsize*nint(dx / boxsize) 
                end if 
                r2      = sum( dx**2 ) ! squared distance between the points
                bin_loc = locate_bin( r2, r2bins, nr ) ! bin index
                if ( bin_loc > 0 ) then               
                    cnts(bin_loc) = cnts(bin_loc) + 1 ! add to bin
                end if
            
            end do
        end do

    end subroutine count_pairs_sequential

    subroutine count_pairs_parallel(pid, nthreads,                            & 
                                    periodic, boxsize, nr, r2bins, cnts,      &
                                    npts1, start1, count1, grid_data1, pos1,  &
                                    npts2, start2, count2, grid_data2, pos2   &
        )
        !! Counting the pairs between sections of position buffers. This 
        !! is a parallel version of the broute force pair counting between 
        !! two cells.
        
        integer(c_int)    , intent(in)    :: periodic, nthreads
        integer(c_int64_t), intent(in)    :: pid, npts1, npts2, nr
        integer(c_int64_t), intent(in)    :: start1, start2, count1, count2
        integer(c_int64_t), intent(in)    :: grid_data1(npts1), grid_data2(npts2)
        real(c_double)    , intent(in)    :: pos1(3,npts1), pos2(3,npts2), boxsize(3)
        real(c_double)    , intent(in)    :: r2bins(nr)
        integer(c_int64_t), intent(inout) :: cnts(:) ! Must be of size nr-1

        character(256)     :: fn
        integer(c_int)     :: tid, fu
        real(c_double)     :: dx(3), r2
        integer(c_int64_t) :: npairs, k, i1, i2, bin_loc
        integer(c_int64_t), allocatable :: local_cnts(:)

        npairs = count1*count2 ! Total number of possible pairs

        !$OMP  PARALLEL DEFAULT(SHARED) &
        !$OMP& PRIVATE(tid, k, i1, i2, bin_loc, r2, dx, local_cnts, fu, fn)

        allocate( local_cnts(nr) )

        !$OMP DO SCHEDULE(STATIC)
        do k = 0, npairs-1
            
            ! Index of the 1-st point in position buffer
            i1 = modulo(k, count1) + start1 ! cell index
            i1 = grid_data1(i1) 
            
            ! Index of the 2-nd point in position buffer
            i2 = modulo(k / count1, count2) + start2 ! cell index
            i2 = grid_data2(i2) 
            
            ! Check pair distance
            dx = pos2(1:3, i2) - pos1(1:3, i1) 
            if ( periodic /= 0 ) then
                ! Correcting the coordinate distances for periodic boundary
                dx = dx - boxsize*nint(dx / boxsize) 
            end if 
            r2      = sum( dx**2 ) ! squared distance between the points
            bin_loc = locate_bin( r2, r2bins, nr ) ! bin index
            if ( bin_loc > 0 ) then
                local_cnts(bin_loc) = local_cnts(bin_loc) + 1 ! add to bin
            end if

        end do
        !$OMP END DO
        
        ! Wrie data to temp file
        tid = omp_get_thread_num() + 1 ! thread ID 
        fu  = 10 + tid ! file unit for this thread
        write( fn, '(i0,".",i0,"cnts.bin")' ) pid, tid ! filename for this thread
        open(newunit=fu, file=fn, access='stream', form='unformatted',      &
             convert='little_endian', status='unknown', position='append',  &
             action='write'                                                 &
        )
        write(fu) local_cnts
        close(fu)

        deallocate( local_cnts )

        !$OMP END PARALLEL

        ! Accumulating counts
        allocate( local_cnts(nr) )
        do tid = 1, nthreads
            fu  = 10 + tid ! file unit for this thread
            write( fn, '(i0,".",i0,"cnts.bin")' ) pid, tid ! filename for this thread
            open(newunit=fu, file=fn, access='stream', form='unformatted',      &
                 convert='little_endian', status='unknown', position='append',  &
                 action='write'                                                 &
            )
            read(fu) local_cnts
            cnts(1:nr-1) = cnts(1:nr-1) + local_cnts(1:nr-1)
            close(fu, status='delete') ! file is deleted on close 
        end do
        deallocate( local_cnts )

    end subroutine count_pairs_parallel

! Helper functions:

    function locate_bin(val, edges, n) result(i)
        !! Locate the bin containing the given value using a binary search 
        !! on the bin edges array. 
        real(c_double)    , intent(in) :: val, edges(n)
        integer(c_int64_t), intent(in) :: n
        integer(c_int64_t) :: i, lo, hi, m
    
        ! Out-of-range check
        if (val < edges(1) .or. val >= edges(n)) then
            i = -1_c_int64_t
            return
        end if

        ! Binary search
        lo = 1_c_int64_t
        hi = n
        do while ( hi - lo > 1 )
            m = (lo + hi) / 2
            if ( val < edges(m) ) then
                hi = m
            else
                lo = m
            end if
        end do

        i = lo ! bin index such that edges(i) <= val < edges(i+1)
    
    end function locate_bin

end module pair_counting_mod