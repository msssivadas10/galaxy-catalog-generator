module pair_counting_mod
    !! Count pairs between two sets of points. 
    !! (TODO: testing)
    
    use omp_lib
    use iso_c_binding
    use spatial_hash_mod
    use pair_utils_mod
    implicit none

    private
    public :: count_pairs
    
contains

    subroutine count_pairs(pid, auto, periodic, box, nr, rbins, cnts, ncells, &
                           npts1, grid_info1, grid_data1, pos1,               &
                           npts2, grid_info2, grid_data2, pos2,               &
                           nthreads, error_code                               &
        ) bind(c)
        !! Count the number of pairs between two sets of points, given their grid
        !! hash representation.

        integer(c_int64_t), intent(in), value :: pid
        !! Process ID
        
        integer(c_int), intent(in), value :: auto
        !! Flag for taking auto correlation (if both sets are same)
     
        integer(c_int), intent(in), value :: periodic
        !! Flag for using periodic boundary conditions
        
        type(boxinfo_t), intent(in) :: box
        !! Details about the space (both sets must be in the same box)

        integer(c_int64_t), intent(in), value :: nr        
        real(c_double)    , intent(in)        :: rbins(nr) !! Distance bin edges
        integer(c_int64_t), intent(out)       :: cnts(nr)  !! Pair counts

        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid (same for both sets).

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

        character(256)     :: fn
        real(c_double)     :: r2bins(nr)
        integer(c_int)     :: fi, ios
        integer(c_int64_t) :: file_size, n_pairs, n_pairs_total, n_pairs_processed
        integer(c_int64_t) :: chunk_size  
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
        
        chunk_size = 1000
        allocate( pairs(2, chunk_size) )

        ! First, processing the 'normal' cell pairs stored in `pid.cplist.bin`.
        ! Pairs in these set are processed in parallel. 
        fi = 9 ! file unit for input
        write(fn, '(i0,".cplist.bin")') pid ! filename for input
        open(newunit=fi, file=fn, access='stream', form='unformatted', &
             convert='little_endian', status='old', action='read'      &
        ) ! input file
        inquire(fi, size=file_size)
        if ( file_size > 0 ) then
            
            ! Calculate the number if pairs: the file is a stream of int64 cell
            ! index pairs.
            n_pairs_total = file_size / (2*c_sizeof(1_c_int64_t))
            write(*,'("info: found ",i0," cell pairs in ",a," (filesize: ",i0," bytes)")') &
                n_pairs_total, trim(fn), file_size

            ! Loading pairs as chunks, and processing items in parallel...
            n_pairs_processed = 0
            do while ( n_pairs_processed < n_pairs_total )
                n_pairs = min( chunk_size, n_pairs_total - n_pairs_processed )
                read(fi, iostat=ios) pairs(:, 1:n_pairs)
                if ( ios /= 0 ) exit ! EOF or error 
                
                ! Counting pairs
                call distributed_pair_count_cells(pid, auto, periodic, n_pairs, pairs, &
                                                  box, nr, r2bins, cnts, ncells,       &
                                                  npts1, grid_info1, grid_data1, pos1, &
                                                  npts2, grid_info2, grid_data2, pos2, &
                                                  nthreads                             &
                )
                
                n_pairs_processed = n_pairs_processed + n_pairs
                write(*,'("info: processed ",i0," of ",i0," cell pairs (",f6.2,"%)")') & 
                    n_pairs_processed, n_pairs_total,                                  &
                    100*(dble(n_pairs_processed)/dble(n_pairs_total)) ! percentage completed
            end do
            
        end if
        close(fi, status='delete') ! file is deleted on close

        ! Second, special cell pairs (with much larger counts compared to average)
        ! are processed. For these, parellization is over the point pairs, as this 
        ! is much efficient than over cell pairs.
        fi = 9 ! file unit for input
        write(fn, '(i0,".cplist.spl.bin")') pid ! filename for input
        open(newunit=fi, file=fn, access='stream', form='unformatted', &
             convert='little_endian', status='old', action='read'      &
        ) ! input file
        inquire(fi, size=file_size)
        if ( file_size > 0 ) then
            
            ! Calculate the number if pairs: the file is a stream of int64 cell
            ! index pairs.
            n_pairs_total = file_size / (2*c_sizeof(1_c_int64_t))
            write(*,'("info: found ",i0," cell pairs in ",a," (filesize: ",i0," bytes)")') &
                n_pairs_total, trim(fn), file_size

            ! Loading pairs as chunks, and processing items in parallel...
            n_pairs_processed = 0
            do while ( n_pairs_processed < n_pairs_total )
                n_pairs = min( chunk_size, n_pairs_total - n_pairs_processed )
                read(fi, iostat=ios) pairs(:, 1:n_pairs)
                if ( ios /= 0 ) exit ! EOF or error 
                
                ! TODO: Counting pairs

                n_pairs_processed = n_pairs_processed + n_pairs
                write(*,'("info: processed ",i0," of ",i0," cell pairs (",f6.2,"%)")') & 
                    n_pairs_processed, n_pairs_total,                                  &
                    100*(dble(n_pairs_processed)/dble(n_pairs_total)) ! percentage completed
            end do
            
        end if
        close(fi, status='delete') ! file is deleted on close

        deallocate( pairs )
        error_code = 0
        
    end subroutine count_pairs

! Parallelized pair counting over cell pairs:

    subroutine distributed_pair_count_cells(pid, auto, periodic, n_pairs, pairs,  &
                                            box, nr, r2bins, cnts, ncells,        &
                                            npts1, grid_info1, grid_data1, pos1,  &
                                            npts2, grid_info2, grid_data2, pos2,  &
                                            nthreads                              &
        )
        !! Paralelly count pairs over a set of cell pairs.

        integer(c_int64_t), intent(in)    :: pid
        integer(c_int)    , intent(in)    :: auto, periodic, nthreads
        integer(c_int64_t), intent(in)    :: n_pairs, pairs(2,n_pairs)
        integer(c_int64_t), intent(in)    :: nr, ncells, npts1, npts2              
        real(c_double)    , intent(in)    :: r2bins(nr), pos1(3,npts1), pos2(3,npts2)
        integer(c_int64_t), intent(in)    :: grid_data1(npts1), grid_data2(npts2)  
        type(cinfo_t)     , intent(in)    :: grid_info1(ncells), grid_info2(ncells) 
        type(boxinfo_t)   , intent(in)    :: box
        integer(c_int64_t), intent(inout) :: cnts(nr)  

        integer(c_int)     :: tid
        integer(c_int64_t) :: p, j1, j2
        integer(c_int64_t), allocatable :: local_cnts(:)

        !$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(tid, p, j1, j2, local_cnts)

        allocate( local_cnts(nr) )
        local_cnts = 0_c_int64_t
        
        !$OMP DO SCHEDULE(static)
        do p = 1, n_pairs
            j1 = pairs(1, p)
            j2 = pairs(2, p)
            call count_pairs_serial(auto, periodic, j1, j2, box,         &
                                    nr, r2bins, local_cnts, ncells,      &
                                    npts1, grid_info1, grid_data1, pos1, &
                                    npts2, grid_info2, grid_data2, pos2  &
            )
        end do
        !$OMP END DO
        
        ! Temporarily save the counts
        tid = omp_get_thread_num() + 1 ! thread ID 
        call save_temporary_counts(pid, tid, nr, local_cnts) 
        deallocate( local_cnts )

        !$OMP END PARALLEL

        ! Load count data from each temporary file and update the global count 
        call merge_temporary_counts(pid, nthreads, nr, cnts)
        
    end subroutine distributed_pair_count_cells

    subroutine count_pairs_serial(auto, periodic, j1, j2,              &
                                  box, nr, r2bins, cnts, ncells,       &
                                  npts1, grid_info1, grid_data1, pos1, &
                                  npts2, grid_info2, grid_data2, pos2  &
        )
        !! Sequencially count pairs between points from two cells j1 and j2.

        integer(c_int)    , intent(in)    :: auto, periodic
        integer(c_int64_t), intent(in)    :: j1, j2, nr, ncells, npts1, npts2              
        real(c_double)    , intent(in)    :: r2bins(nr), pos1(3,npts1), pos2(3,npts2)
        integer(c_int64_t), intent(in)    :: grid_data1(npts1), grid_data2(npts2)  
        type(cinfo_t)     , intent(in)    :: grid_info1(ncells), grid_info2(ncells) 
        type(boxinfo_t)   , intent(in)    :: box
        integer(c_int64_t), intent(inout) :: cnts(nr) 

        integer(c_int64_t) :: i1, i2, i1_start, i2_start, i1_stop, i2_stop
        integer(c_int64_t) :: k1, k2, bin_loc
        real(c_double)     :: dx(3), r2

        ! Range of the block for cell j1
        i1_start = grid_info1(j1)%start
        i1_stop  = i1_start + grid_info1(j1)%count - 1
        
        ! Range of the block for cell j2
        i2_start = grid_info2(j2)%start
        i2_stop  = i2_start + grid_info2(j2)%count - 1

        do i1 = i1_start, i1_stop
            
            ! Actual index of the point in position buffer 
            k1 = grid_data1(i1) 
            
            do i2 = i2_start, i2_stop
                
                ! Actual index of the point in position buffer 
                k2 = grid_data2(i2) 

                ! If counting between same set of points, only count pairs with 
                ! k2 < k1, to avoid double counting.
                if ( auto /= 0 .and. k2 < k1 ) cycle 
                
                ! Check pair distance
                dx = pos2(1:3, k2) - pos1(1:3, k1) 
                if ( periodic /= 0 ) then
                    ! Correcting the coordinate distances for periodic boundary
                    dx = dx - box%boxsize*nint(dx / box%boxsize) 
                end if 
                r2      = sum( dx**2 ) ! squared distance between the points
                bin_loc = locate_bin( r2, r2bins, nr ) ! bin index
                if ( bin_loc > 0 ) then               
                    cnts(bin_loc) = cnts(bin_loc) + 1 ! add to bin
                end if

            end do
        
        end do
        
    end subroutine count_pairs_serial

! Parallelized pair counting over point pairs:

    subroutine sequential_pair_count_cells(pid, auto, periodic, n_pairs, pairs,  &
                                           box, nr, r2bins, cnts, ncells,        &
                                           npts1, grid_info1, grid_data1, pos1,  &
                                           npts2, grid_info2, grid_data2, pos2,  &
                                           nthreads                              &
        )
        !! Sequentially count pairs over a set of cell pairs. But, between each 
        !! cell pairs, pair cunting over point pairs will be in parallel.

        integer(c_int64_t), intent(in)    :: pid
        integer(c_int)    , intent(in)    :: auto, periodic, nthreads
        integer(c_int64_t), intent(in)    :: n_pairs, pairs(2,n_pairs)
        integer(c_int64_t), intent(in)    :: nr, ncells, npts1, npts2              
        real(c_double)    , intent(in)    :: r2bins(nr), pos1(3,npts1), pos2(3,npts2)
        integer(c_int64_t), intent(in)    :: grid_data1(npts1), grid_data2(npts2)  
        type(cinfo_t)     , intent(in)    :: grid_info1(ncells), grid_info2(ncells) 
        type(boxinfo_t)   , intent(in)    :: box
        integer(c_int64_t), intent(inout) :: cnts(nr)

        integer(c_int64_t) :: p, j1, j2

        do p = 1, n_pairs
            j1 = pairs(1, p)
            j2 = pairs(2, p)
            call count_pairs_parallel(pid, auto, periodic, j1, j2, box,    &
                                      nr, r2bins, cnts, ncells,            &
                                      npts1, grid_info1, grid_data1, pos1, &
                                      npts2, grid_info2, grid_data2, pos2, &
                                      nthreads                             &
            )
        end do
        
    end subroutine sequential_pair_count_cells

    subroutine count_pairs_parallel(pid, auto, periodic, j1, j2,         &
                                    box, nr, r2bins, cnts, ncells,       &
                                    npts1, grid_info1, grid_data1, pos1, &
                                    npts2, grid_info2, grid_data2, pos2, &
                                    nthreads                             &
        )
        !! Parallely count pairs between points from two cells j1 and j2.

        integer(c_int)    , intent(in)    :: auto, periodic, nthreads
        integer(c_int64_t), intent(in)    :: pid, j1, j2, nr, ncells, npts1, npts2              
        real(c_double)    , intent(in)    :: r2bins(nr), pos1(3,npts1), pos2(3,npts2)
        integer(c_int64_t), intent(in)    :: grid_data1(npts1), grid_data2(npts2)  
        type(cinfo_t)     , intent(in)    :: grid_info1(ncells), grid_info2(ncells) 
        type(boxinfo_t)   , intent(in)    :: box
        integer(c_int64_t), intent(inout) :: cnts(nr)
        
        real(c_double)     :: dx(3), r2
        integer(c_int)     :: tid
        integer(c_int64_t) :: i1, i2, i1_start, i2_start, i1_count, i2_count
        integer(c_int64_t) :: k_stop, k, bin_loc
        integer(c_int64_t), allocatable :: local_cnts(:)

        ! Range of the block for cell j1
        i1_start = grid_info1(j1)%start
        i1_count = grid_info1(j1)%count
        
        ! Range of the block for cell j2
        i2_start = grid_info2(j2)%start
        i2_count = grid_info2(j2)%count

        ! Total number of point pairs
        k_stop = i1_count * i2_count

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid, k, i1, i2, dx, r2, bin_loc)

        allocate( local_cnts(nr) )
        local_cnts = 0_c_int64_t
        
        !$OMP DO SCHEDULE(static)
        do k = 1, k_stop
            
            ! Index of the 1-st point in position buffer
            i1 = modulo(k, i1_count) + i1_start ! cell index
            i1 = grid_data1(i1) 
            
            ! Index of the 2-nd point in position buffer
            i2 = modulo(k / i1_count, i2_count) + i2_start ! cell index
            i2 = grid_data2(i2)

            ! If counting between same set of points, only count pairs with 
            ! k2 < k1, to avoid double counting.
            if ( auto /= 0 .and. i2 < i1 ) cycle 

            ! Check pair distance
            dx = pos2(1:3, i2) - pos1(1:3, i1) 
            if ( periodic /= 0 ) then
                ! Correcting the coordinate distances for periodic boundary
                dx = dx - box%boxsize*nint(dx / box%boxsize) 
            end if 
            r2      = sum( dx**2 ) ! squared distance between the points
            bin_loc = locate_bin( r2, r2bins, nr ) ! bin index
            if ( bin_loc > 0 ) then
                local_cnts(bin_loc) = local_cnts(bin_loc) + 1 ! add to bin
            end if
            
        end do
        !$OMP END DO

        ! Temporarily save the counts:
        tid = omp_get_thread_num() + 1 ! thread ID 
        call save_temporary_counts(pid, tid, nr, local_cnts) 
        deallocate( local_cnts )

        !$OMP END PARALLEL
        
        ! Load count data from each temporary file and update the global count 
        call merge_temporary_counts(pid, nthreads, nr, cnts)
        
    end subroutine count_pairs_parallel

! Helper functions:

    subroutine save_temporary_counts(pid, tid, nr, cnts)
        !! Write count data to a temporary file.
        integer(c_int64_t), intent(in) :: pid, nr
        integer(c_int)    , intent(in) :: tid
        integer(c_int64_t), intent(in) :: cnts(nr)  

        character(256) :: fn
        integer(c_int) :: fu

        fu  = 10 + tid ! file unit for this thread
        write(fn, '(i0,".",i0,".cnts.tmp")' ) pid, tid ! filename for this thread
        open(newunit=fu, file=fn, access='stream', form='unformatted',     &
            convert='little_endian', status='unknown', position='append',  &
            action='write'                                                 &
        )
        write(fu) cnts
        close(fu)
        
    end subroutine save_temporary_counts

    subroutine merge_temporary_counts(pid, nthreads, nr, cnts)
        !! Accumulate counts from temporary data.
        integer(c_int64_t), intent(in)    :: pid, nr
        integer(c_int)    , intent(in)    :: nthreads
        integer(c_int64_t), intent(inout) :: cnts(nr)  

        character(256) :: fn
        integer(c_int) :: fu, tid
        integer(c_int64_t), allocatable :: local_cnts(:)

        allocate( local_cnts(nr)  )

        fu = 10
        do tid = 1, nthreads
            write(fn, '(i0,".",i0,".tmp")' ) pid, tid ! filename for this thread
            open(newunit=fu, file=fn, access='stream', form='unformatted', &
                 convert='little_endian', status='old', action='read'      &
            )
            read(fu) local_cnts
            cnts = cnts + local_cnts
            close(fu, status='delete') ! file is deleted on close 
        end do
        
        deallocate( local_cnts  )
        
    end subroutine merge_temporary_counts

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