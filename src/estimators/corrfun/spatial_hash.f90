module spatial_hash_mod
    !! Spatial hashing for fast pair counting. Also, contains routines
    !! for spatial statistics calculation.
    
    use omp_lib
    use iso_c_binding
    implicit none

    private
    public :: build_grid_hash, get_grid_stats

    type, public, bind(c) :: cinfo_t
        !! A struct containing information about a grid cell. This, used 
        !! along with a sequence of point index blocks, can be used as an
        !! efficient spatial hash.    
        integer(c_int64_t) :: start !! Start index of a block
        integer(c_int64_t) :: count !! Size of the block  
    end type cinfo_t

    type, public, bind(c) :: boxinfo_t
        !! A struct for storing the details about the box 
        real(c_double)     ::  boxsize(3) !! Size of the box
        real(c_double)     ::   origin(3) !! Coordinates of the origin - lower left corner
        integer(c_int64_t) :: gridsize(3) !! Number of cells on each direction
        real(c_double)     :: cellsize(3) !! Size of the cell along each direction
    end type boxinfo_t

    type, public, bind(c) :: gstats_t
        !! A struct storing grid statistics
        integer(c_int64_t) :: count !! Number of cells
        integer(c_int64_t) :: min   !! Minimum value
        integer(c_int64_t) :: max   !! Maximum value
        integer(c_int64_t) :: empty !! Number of cells that are empty
        real(c_double)     :: avg   !! Average 
        real(c_double)     :: var   !! Variance
    end type gstats_t
    
contains

    subroutine build_grid_hash(pid, npts, pos, box, ncells, grid_info, &
                               index_list, nthreads, error_code        &
        ) bind(c)
        !! Calculate a grid spatial hash for the positions for fast and efficient 
        !! pair counting.

        integer(c_int64_t), intent(in), value :: pid
        !! A unique positive integer value (Used for safe file IO). 

        integer(c_int64_t), intent(in), value :: npts
        !! Number of points in the catalog

        real(c_double), intent(in) ::  pos(3, npts)
        !! Position

        type(boxinfo_t), intent(inout) :: box
        !! Details about the space

        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid - must be equal to `product(gridsize)`

        integer(c_int64_t), intent(out) :: index_list(npts)
        !! Indices of points, grouped by sorted cell index 

        type(cinfo_t), intent(out) :: grid_info(ncells)
        !! Grid data - start index and size of each cell group. 

        integer(c_int), intent(in), value :: nthreads
        !! Number of threads to use

        integer(c_int), intent(out) :: error_code
        !! Error code (0=success, 1=error)

        integer(c_int) :: fu, tid
        integer(c_int64_t) :: cell(3), i, j, p, j_start, j_end, chunk_size, rem
        integer(c_int64_t), allocatable :: cell_index(:), local_count(:)
        integer(c_int64_t), allocatable :: partial(:), offset(:)
        character(256)    , allocatable :: fn(:)

        error_code = 1
        if ( product(box%gridsize) /= ncells ) stop "mismatch in gridsize and ncells"
        
        call omp_set_num_threads(nthreads) ! set number of threads
        
        ! Setting up...
        box%cellsize = box%boxsize / box%gridsize 
        do j = 1, ncells
            grid_info(j)%count = 0_c_int64_t ! initialise all count to 0
            grid_info(j)%start = 0_c_int64_t ! initialise all start to 0
        end do
        ! Names for thread specific temp files
        allocate( fn(nthreads) )
        do tid = 1, nthreads
            write( fn(tid), '(i0,".",i0,".tmp")' ) pid, tid ! filename for this thread
        end do
        allocate( cell_index(npts) )

        ! -- Step 1 -- 
        ! Calculate the flattened index of the grid cell that contain the point. Also, 
        ! calculate the cell histogram - number of points in each cell. This part is 
        ! parallelised over multiple threads.
        ! 
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid, cell, local_count, i, j, fu)
        
        ! Allocate local histogram for this thread
        allocate( local_count(ncells) )
        local_count(:) = 0_c_int64_t
        
        !$OMP DO SCHEDULE(static)
        do i = 1, npts
            cell = int(( pos(:, i) - box%origin ) / box%cellsize) ! 3D cell index

            ! Apply bounds: all the cell indices must be within [0, gridsize-1].
            ! Any value outside that will be clipped. 
            cell = min(max([0_c_int64_t, 0_c_int64_t, 0_c_int64_t], cell), box%gridsize-1)

            ! Flattening the index (0 based)
            cell_index(i) = cell(1) + box%gridsize(1)*( cell(2) + box%gridsize(2)*cell(3) )
            
            ! Increment count for this cell
            j = cell_index(i) + 1
            local_count(j) = local_count(j) + 1

        end do
        !$OMP END DO

        ! For a thread-safe accumulatiion of the counts from different threads, each 
        ! thread will store the count to a temporary file, then update using that.
        tid = omp_get_thread_num() + 1 ! thread ID 
        fu  = 10 + tid ! file unit for this thread
        open(newunit=fu, file=fn(tid), access='stream', form='unformatted', &
             convert='little_endian', status='unknown', position='append',  &
             action='write'                                                 &
        )
        write(fu) local_count
        close(fu)
        
        deallocate( local_count  )

        !$OMP END PARALLEL

        ! Load count data from each temporary file and update the global count 
        allocate( local_count(ncells)  )
        fu = 10
        do tid = 1, nthreads
            open(newunit=fu, file=fn(tid), access='stream', form='unformatted', &
                 convert='little_endian', status='old', action='read'           &
            )
            read(fu) local_count
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
            do j = 1, ncells
                ! Accumulating counts: No need to sync this part, as each thread 
                ! works on a seperate block.
                grid_info(j)%count = grid_info(j)%count + local_count(j) 
            end do
            !$OMP END PARALLEL DO 
            close(fu, status='delete') ! file is deleted on close 
        end do
        deallocate( local_count  )

        ! -- Setp 2 -- 
        ! Prefix sum: calculating cell start index using counts. This is calculated in 
        ! parallel.
        
        allocate( partial(nthreads+1) )
        partial(:) = 0_c_int64_t

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid, j, j_start, j_end, chunk_size, rem)
        
        tid        = omp_get_thread_num()
        chunk_size = ncells / nthreads
        rem        = modulo(ncells, nthreads)
        j_start    = tid*chunk_size + min(tid, rem) + 1
        j_end      = j_start + chunk_size - 1
        if (tid < rem) j_end = j_end + 1

        grid_info(j_start)%start = 0_c_int64_t
        do j = j_start+1, j_end
            grid_info(j)%start = grid_info(j-1)%start + grid_info(j-1)%count
        end do
        partial(tid+2) = grid_info(j_end)%start + grid_info(j_end)%count
        !$OMP BARRIER

        ! Prefix sum over chunk totals:
        !$OMP SINGLE
        partial(1) = 1_c_int64_t
        do j = 2, nthreads+1
            partial(j) = partial(j) + partial(j-1)
        end do
        !$OMP END SINGLE

        ! Final start index:
        do j = j_start, j_end
            grid_info(j)%start = grid_info(j)%start + partial(tid+1)
        end do

        !$OMP END PARALLEL

        deallocate( partial )

        ! -- Step 3 --
        ! Scattering points into a sorted index array. This will create an array of 
        ! indices, where points in a specific cell are grouped together. Using this
        ! along with the grid_info, one can map points to its corresponding cells. 
        allocate( offset(ncells) )
        offset(:) = 0_c_int64_t

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, p)
        !$OMP DO
        do i = 1, npts
            j = cell_index(i) + 1
            !$OMP ATOMIC CAPTURE
            p         = offset(j)
            offset(j) = offset(j) + 1
            !$OMP END ATOMIC
            index_list( grid_info(j)%start + p ) = i
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        deallocate( offset )
        deallocate( cell_index )
        
        error_code = 0

    end subroutine build_grid_hash

    subroutine get_grid_stats(ncells, grid_info, nthreads, stats, error_code) bind(c)
        !! Calculate cell statistics  

        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid - must be equal to `product(gridsize)`

        type(cinfo_t), intent(in) :: grid_info(ncells)
        !! Grid data - start index and size of each cell group. 

        integer(c_int), intent(in), value :: nthreads
        !! Number of threads to use

        type(gstats_t), intent(inout) :: stats
        !! Grid statistics

        integer(c_int), intent(out) :: error_code
        !! Error code (0=success, 1=error)

        integer(c_int64_t) :: tid, j, n_total
        real(c_double)     :: delta1, delta2
        type(gstats_t), allocatable :: st(:)

        error_code = 1
        call omp_set_num_threads(nthreads)

        ! Get stats on blocks: average and variance are calculated using 
        ! Welford's online algorithm, per thread. These values are combined to
        ! get the actual values at the end. 
        ! Ref: <https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance>.
        allocate( st(nthreads) )
        do tid = 1, nthreads
            st(tid)%count =  0_c_int64_t
            st(tid)%min   = huge( st(tid)%min )
            st(tid)%max   = -1_c_int64_t
            st(tid)%empty =  0_c_int64_t
            st(tid)%avg   =  0._c_double
            st(tid)%var   =  0._c_double
        end do

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid, j, delta1, delta2)
        tid = omp_get_thread_num() + 1

        !$OMP DO SCHEDULE(static)
        do j = 1, ncells
            st(tid)%count = st(tid)%count + 1
            delta1        = grid_info(j)%count - st(tid)%avg
            st(tid)%avg   = st(tid)%avg + delta1 / dble( st(tid)%count )
            delta2        = grid_info(j)%count - st(tid)%avg
            st(tid)%var   = st(tid)%var + delta1*delta2
            st(tid)%min   = min( st(tid)%min, grid_info(j)%count )
            st(tid)%max   = max( st(tid)%max, grid_info(j)%count )
            if ( grid_info(j)%count < 1 ) st(tid)%empty = st(tid)%empty + 1
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        
        ! Merging the stats on blocks
        stats = st(1)
        do tid = 2, nthreads
            if ( st(tid)%count > 0 ) then
                n_total     = st(tid)%count + stats%count 
                delta1      = st(tid)%avg   - stats%avg
                stats%avg   = stats%avg + delta1 * dble(st(tid)%count) / dble(n_total)
                stats%var   = stats%var + st(tid)%var &
                                + delta1**2 * dble(stats%count*st(tid)%count) / dble(n_total)
                stats%count = n_total
                stats%empty = st(tid)%empty + stats%empty
            end if
            stats%min = min( stats%min, st(tid)%min )
            stats%max = max( stats%max, st(tid)%max )
        end do
        stats%var = stats%var / dble(stats%count) ! biased sample variance
        
        deallocate( st )
        error_code = 1

    end subroutine get_grid_stats
    
end module spatial_hash_mod