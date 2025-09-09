module corrfun_mod
    !! Pair counting and 2-point correlation function calculations.

    use omp_lib
    use iso_c_binding
    implicit none

    private
    public :: build_grid_hash

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

    subroutine count_pairs(auto, periodic, nr, rbins, counts, ncells, box, &
                           npts1, pos1, index_list1, grid_info1,           & ! set-1 
                           npts2, pos2, index_list2, grid_info2,           & ! set-2
                           error_code                                      & 
        ) bind(c)
        !! Count number of pairs between two sets of points.

        integer(c_int), intent(in), value :: auto
        !! Flag indicating pair counting between same sets of points
        
        integer(c_int), intent(in), value :: periodic
        !! Flag for using periodic boundary conditions

        integer(c_int64_t), intent(in), value :: nr   
        !! Size of distance bins array and count array

        real(c_double), intent(in) :: rbins(nr)
        !! Distance bin edges

        integer(c_int64_t), intent(out) :: counts(nr)
        !! Pair counts (last item will be 0 always).

        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid (same for both sets).

        type(boxinfo_t), intent(in) :: box
        !! Details about the space (both sets must be in the same box)

        integer(c_int64_t), intent(in), value :: npts1
        !! Number of points in set-1

        real(c_double), intent(in) :: pos1(3, npts1)
        !! Set-1 positions

        integer(c_int64_t), intent(in) :: index_list1(npts1)
        !! Sorted cell index list for set-1

        type(cinfo_t), intent(in) :: grid_info1(ncells)
        !! Grid specification for set-1

        integer(c_int64_t), intent(in) :: npts2
        !! Number of points in set-2

        real(c_double), intent(in) :: pos2(3, npts2)
        !! Set-2 positions

        integer(c_int64_t), intent(in) :: index_list2(npts2)
        !! Sorted cell index list for set-2

        type(cinfo_t), intent(in) :: grid_info2(ncells)
        !! Grid specification for set-2
        
        integer(c_int), intent(out) :: error_code
        !! Error code (0=success, 1=error)

        integer(c_int64_t) :: j1, cell1(3), delta(3), crange(3, 2)

        error_code = 1

        ! Number of cells to check on each direction 
        delta = ceiling( rbins(nr) / box%cellsize ) 

        do j1 = 1, ncells

            if ( grid_info1(j1)%count < 1 ) cycle ! Cell is empty
            
            ! Finding 3D cell index from flat index
            cell1(1) = modulo( j1, box%gridsize(1) ); cell1(3) = j1 / box%gridsize(1) 
            cell1(2) = modulo( cell1(3), box%gridsize(2) )
            cell1(3) = cell1(3) / box%gridsize(2)

            crange(:,1) = cell1 - delta
            crange(:,2) = cell1 + delta
            if ( periodic /= 0 ) then
                ! Using periodic boundary conditions: wrap around if the cell index 
                ! is outside range. 
                crange(:,1) = modulo( crange(:,1)-1, box%gridsize ) + 1
                crange(:,2) = modulo( crange(:,2)-1, box%gridsize ) + 1
            else
                ! No periodic wrapping: clip the value between 0 and gridsize
                crange(:,1) = max( min( 0_c_int64_t, crange(:,1) ), box%gridsize )
                crange(:,2) = max( min( 0_c_int64_t, crange(:,2) ), box%gridsize )
            end if

        end do

        error_code = 0
            
    end subroutine count_pairs

end module corrfun_mod
