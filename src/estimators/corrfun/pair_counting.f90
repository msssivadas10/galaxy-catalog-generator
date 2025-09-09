module pair_counting_mod
    !! Pair counting and 2-point correlation function calculations.

    use omp_lib
    use iso_c_binding
    use spatial_hash_mod
    implicit none

    private
    public :: count_pairs
    
    integer(c_int), parameter :: POINT_PARALLEL = 2
    integer(c_int), parameter :: CELLS_PARALLEL = 1

contains

    subroutine get_parallel_strategy(auto, ncells, grid_info1, grid_info2, nthreads, &
                                     strategy, used_set                              &
        ) 
        !! Select an efficient strategy for parallelization, based on grid 
        !! statistics.

        integer(c_int), intent(in), value :: auto
        !! Auto correlation flag: 1=both sets are same.

        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid - must be equal to `product(gridsize)`

        ! integer(c_int64_t), intent(in), value :: npts1
        ! !! Number of points in set-1
        
        ! integer(c_int64_t), intent(in), value :: npts2
        ! !! Number of points in set-2

        type(cinfo_t), intent(in) :: grid_info1(ncells)
        !! Grid data (set-1) - start index and size of each cell group. 

        type(cinfo_t), intent(in) :: grid_info2(ncells)
        !! Grid data (set-2) - start index and size of each cell group. 

        integer(c_int), intent(in), value :: nthreads
        !! Number of threads to use

        integer(c_int), intent(out) :: strategy
        !! Strategy (1=over cells, 2=over points, -99=error)

        integer(c_int), intent(out) :: used_set
        !! Select the set to parallelize

        type(gstats_t) :: stats(2)
        integer(c_int) :: error_code, set
        real(c_double) :: f_empty(2), cv(2), r(2), score(2)

        strategy = -99

        ! Get stats for set-1
        call get_grid_stats(ncells, grid_info1, nthreads, stats(1), error_code)
        if ( error_code /= 0 ) return 
        
        ! Get stats for set-2
        if ( auto /= 1 ) then
            ! Both sets are same: reuse the stats
            stats(2) = stats(1)
        else
            call get_grid_stats(ncells, grid_info2, nthreads, stats(2), error_code)
            if ( error_code /= 0 ) return 
        end if

        ! Calculate some indicators
        do set = 1, 2
            ! Fraction of empty cells:
            f_empty(set) = dble(stats(set)%empty) / dble(stats(set)%count)
    
            ! Coefficient of variation:
            cv(set) = sqrt(stats(set)%var) / stats(set)%avg

            ! Max to mean ratio:
            r(set) = dble(stats(set)%max) / stats(set)%avg

        end do
        
        ! Weighted score based on these stats: set with lowest score is 
        ! selected:
        score = f_empty + cv + log(1._c_double + r)
        if ( score(1) <= score(2) ) then
            used_set = 1
        else
            used_set = 2
        end if

        ! Deciding strategy, based on the set used
        if ( f_empty(used_set) > 0.7_c_double ) then
            ! Set is more sparse: parallelise over points
            strategy = POINT_PARALLEL
        else if ( cv(used_set) < 0.5_c_double ) then
            ! Counts are more uniform among cells: parallelize over cells
            ! is the best strategy
            strategy = CELLS_PARALLEL
        else
            ! Default strategy is to parallelize over cells, unless there are 
            ! some cells with unusually large counts...
            if ( R(used_set) > 1000._c_double ) then
                strategy = POINT_PARALLEL
            else
                strategy = CELLS_PARALLEL
            end if 
        end if

        ! Look for cell wth count, much much larger than the average value. 
        ! For such cell pairs, a parallel loop over points is used, so that
        ! a thread will not get a huge number of points to look for...
        if ( strategy == POINT_PARALLEL ) then
            ! Parallelization already over points - no need to look for
            ! special cells, to parallelize over points.
            return
        end if
        
    end subroutine get_parallel_strategy

    subroutine convert_to_3d_index(gridsize, j, cell)
        !! Calculate the 3D cell indices from a flat index.
        integer(c_int64_t), intent(in)  :: gridsize(3), j
        integer(c_int64_t), intent(out) :: cell(3)
    
        ! Using `flat index = c1 + gs1(c2 + gs2*c3)`...
        cell(1) = modulo( j, gridsize(1) ); cell(3) = j / gridsize(1) 
        cell(2) = modulo( cell(3), gridsize(2) )
        cell(3) = cell(3) / gridsize(2)
        
    end subroutine convert_to_3d_index
    
    function convert_to_flat_index(gridsize, cell) result(j)
        !! Calculate the flat index from 3D cell indices.
        integer(c_int64_t), intent(in) :: gridsize(3), cell(3)
        integer(c_int64_t) :: j
    
        ! Using `flat index = c1 + gs1(c2 + gs2*c3)`...
        j = cell(1) + gridsize(1)*( cell(2) + gridsize(2)*cell(3) )

    end function convert_to_flat_index

    function apply_bc(cell, gridsize, periodic) result(cell_bc)
        !! Apply a periodic boundary condition to the cell index.
        integer(c_int)    , intent(in) :: periodic
        integer(c_int64_t), intent(in) :: gridsize(3), cell(3)
        integer(c_int64_t) :: cell_bc(3)
    
        if ( periodic /= 0 ) then
            ! Using periodic boundary conditions: wrap around if the cell  
            ! index is outside range. 
            cell_bc = modulo( cell-1, gridsize ) + 1
        else
            ! No periodic wrapping: clip the value between 0 and gridsize
            cell_bc = max( min( 0_c_int64_t, cell ), gridsize )
        end if
        
    end function apply_bc

    subroutine count_pairs_sequential(n1, pos1, n2, pos2, nr, rbins, counts)
        !! Count number of pairs from two sets of points in sequential 
        !! order.
        integer(c_int64_t), intent(in)  :: n1, n2, nr
        real(c_double)    , intent(in)  :: pos1(3,n1), pos2(3,n2), rbins(nr)
        integer(c_int64_t), intent(out) :: counts(nr)
    
        
    end subroutine count_pairs_sequential

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

        integer(c_int64_t) :: j1, j2, cell1(3), delta(3), cstart(3), cstop(3), &
                              c1, c2, c3

        error_code = 1

        ! Number of cells to check on each direction 
        delta = ceiling( rbins(nr) / box%cellsize ) 

        do j1 = 1, ncells

            if ( grid_info1(j1)%count < 1 ) cycle ! Cell is empty
            
            ! Finding 3D cell index from flat index
            call convert_to_3d_index(box%gridsize, j1, cell1)

            ! Find the neighbouring cells and count pairs
            cstart = apply_bc(cell1 - delta, box%gridsize, periodic)
            cstop  = apply_bc(cell1 + delta, box%gridsize, periodic)
            do c1 = cstart(1), cstop(1)
                do c2 = cstart(2), cstop(2)
                    do c3 = cstart(3), cstop(3)

                        ! Flat index of the neighbour cell
                        j2 = convert_to_flat_index(box%gridsize, [c1, c2, c3])
                        
                        ! Count number of paris
                        if ( grid_info1(j1)%count < 1 ) cycle ! Cell is empty

                    end do
                end do
            end do


        end do

        error_code = 0
            
    end subroutine count_pairs

end module pair_counting_mod
