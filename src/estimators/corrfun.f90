module corrfun_mod
    !! Pair counting and 2-point correlation function calculations.
    use iso_c_binding
    implicit none

    ! private

    type, public, bind(c) :: cinfo_t
        integer(c_int64_t) :: start 
        integer(c_int64_t) :: count    
    end type cinfo_t
    
contains

    subroutine build_grid_hash(pos, npts, boxsize, origin, gridsize, ncells, &
                               grid_group, cellsize, grid_info, error_code   & 
        ) bind(c)
        !! Calculate a grid spatial hash for the positions for fast and 
        !! efficient pair counting.

        integer(c_int64_t), intent(in), value :: npts
        !! Number of points in the catalog

        real(c_double), intent(in) ::  pos(3, npts)
        !! Position

        real(c_double), intent(in) :: boxsize(3)
        !! Size of the bounding box for all positions
        
        real(c_double), intent(in) :: origin(3)
        !! Coordinates of the origin of the bounding box. All points must
        !! be within `[origin, origin + boxsize]` range.

        integer(c_int64_t), intent(in) :: gridsize(3)
        !! Size of the grid on each direction

        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid - must be equal to `product(gridsize)`

        integer(c_int64_t), intent(out) :: grid_group(npts)
        !! Indices of points, grouped by sorted cell index 

        type(cinfo_t), intent(out) :: grid_info(ncells)
        !! Grid data - start index and size of each cell group. 

        real(c_double), intent(out) :: cellsize(3)
        !! Size of a cell

        integer(c_int), intent(out) :: error_code
        !! Error code (0=success, 1=error)

        integer(c_int64_t) :: cell_index(npts), cell(3), i, j

        error_code = 1
        if ( product(gridsize) /= ncells ) stop "mismatch in gridsize and ncells"

        ! Calculate the flattened index of the grid cell that contain the 
        ! point.
        cellsize = boxsize / gridsize 
        do i = 1, npts
            cell = int(( pos(:, i) - origin ) / cellsize) ! 3D cell index

            ! Apply bounds: all the cell indices must be within [0, gridsize-1].
            ! Any value outside that will be clipped. 
            cell = min(max([0_c_int64_t, 0_c_int64_t, 0_c_int64_t], cell), gridsize-1)

            ! Flattening the index
            cell_index(i) = cell(1) + gridsize(1)*( cell(2) + gridsize(2)*cell(3) )
            grid_group(i) = i
        end do

        ! Sorting the point indices based on cell index. This will create 
        ! groups of points corresponding to same cells.
        call argsort(cell_index, grid_group, npts)
        
        ! Making the grid info: this is a summary of the grid groups, telling
        ! the index for the start of a points group and size of that group. 
        ! Each group correspond to a set of points in a certain cell.
        do i = 1, ncells
            grid_info(i)%start = 0; grid_info(i)%count = 0
        end do
        if ( npts > 0 ) then
            j = cell_index( grid_group(1) ) + 1
            grid_info(j)%start = 1 
            grid_info(j)%count = 1
            do i = 2, npts
                j = cell_index( grid_group(i) ) + 1
                if ( grid_info(j)%count == 0 ) grid_info(j)%start = i
                grid_info(j)%count = grid_info(j)%count + 1 
            end do
        end if

        error_code = 0
        
    end subroutine build_grid_hash

    subroutine argsort(a, idx, n)
        !! Sorting based on heapsort.
        integer(c_int64_t), intent(in)    :: n
        integer(c_int64_t), intent(in)    :: a(n)
        integer(c_int64_t), intent(inout) :: idx(n)
        integer(c_int64_t) :: i, t
      
        idx = [( i, i = 1, n )] ! Initialize indices
        do i = n/2, 1, -1
           call sift_down(i, n, a, idx)
        end do
        do i = n, 2, -1
           t = idx(1); idx(1) = idx(i); idx(i) = t  ! Swap root with last
           call sift_down(1_c_int64_t, i-1, a, idx) ! Restore heap
        end do

    end subroutine argsort

    subroutine sift_down(start, finish, a, idx)
        integer(c_int64_t), intent(in)    :: start, finish
        integer(c_int64_t), intent(in)    :: a(:)
        integer(c_int64_t), intent(inout) :: idx(:)
        integer(c_int64_t) :: root, child, swap, t
    
        root = start
        do while (2*root <= finish)
           child = 2*root
           swap  = root
           if ( a(idx(swap)) < a(idx(child)) ) swap = child
           if (child+1 <= finish) then
              if ( a(idx(swap)) < a(idx(child+1)) ) swap = child+1
           end if
           if (swap == root) then
              return
           else
              t    = idx(root); idx(root) = idx(swap); idx(swap) = t 
              root = swap
           end if
        end do
        
    end subroutine sift_down
    
end module corrfun_mod
