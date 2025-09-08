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

    subroutine build_grid_hash(npts, pos, box, ncells, grid_info, &
                               index_list, nthreads, error_code   &
        ) bind(c)
        !! Calculate a grid spatial hash for the positions for fast and efficient 
        !! pair counting.

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

        integer(c_int64_t) :: cell(3), i, j, p
        integer(c_int64_t), allocatable :: cell_index(:), local_count(:), offset(:)

        error_code = 1
        if ( product(box%gridsize) /= ncells ) stop "mismatch in gridsize and ncells"
        
        ! Setting up...
        box%cellsize = box%boxsize / box%gridsize 
        allocate( cell_index(npts) )
        do i = 1, ncells
            grid_info(i)%count = 0_c_int64_t ! initialise all count to 0
        end do
        
        call omp_set_num_threads(nthreads) ! set number of threads

        ! -- Step 1 -- 
        ! Calculate the flattened index of the grid cell that contain the point. Also, 
        ! calculate the cell histogram - number of points in each cell. This part is 
        ! parallelised over multiple threads.
        
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(cell, local_count, i, j)
        
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

        !$OMP CRITICAL
        do j = 1, ncells 
            ! Accumulating counts from different threads... 
            grid_info(j)%count = grid_info(j)%count + local_count(j) 
        end do
        !$OMP END CRITICAL

        deallocate( local_count  )

        !$OMP END PARALLEL

        ! -- Setp 2 -- 
        ! Prefix sum: calculating cell start index using counts. This part is not 
        ! parallelised, as the number of cells will be usually much less than the 
        ! number of points. 
        grid_info(1)%start = 1_c_int64_t
        do j = 2, ncells
            grid_info(j)%start = grid_info(j-1)%start + grid_info(j-1)%count
        end do

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

end module corrfun_mod
