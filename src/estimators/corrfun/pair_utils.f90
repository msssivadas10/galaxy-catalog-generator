module pair_utils_mod
    !! Helper functions for pair counting.

    use omp_lib
    use iso_c_binding
    use spatial_hash_mod
    implicit none
    
contains

    subroutine enumerate_cell_pairs(pid, npts1, npts2, ncells, rmax, box, periodic, &
                                    grid_info1, grid_info2, nthreads, error_code    &
        )
        !! List all possible cell pairs. Also, marks the pairs that needs special
        !! handling. Data are written to disk, to a file with name based on the 
        !! `pid` value.
        
        integer(c_int64_t), intent(in), value :: pid
        !! Process ID
        
        integer(c_int), intent(in), value :: periodic
        !! Flag for using periodic boundary conditions

        real(c_double), intent(in), value :: rmax
        !! Maximum distance bin value

        integer(c_int64_t), intent(in), value :: ncells
        !! Number of cells in the grid (same for both sets).

        type(boxinfo_t), intent(in) :: box
        !! Details about the space (both sets must be in the same box)

        integer(c_int64_t), intent(in), value :: npts1
        !! Number of points in set-1

        type(cinfo_t), intent(in) :: grid_info1(ncells)
        !! Grid specification for set-1

        integer(c_int64_t), intent(in) :: npts2
        !! Number of points in set-2

        type(cinfo_t), intent(in) :: grid_info2(ncells)
        !! Grid specification for set-2

        integer(c_int), intent(in), value :: nthreads
        !! Number of threads to use
        
        integer(c_int), intent(out) :: error_code
        !! Error code (0=success, 1=error)

        real(c_double), parameter :: F = 1000._c_double

        real(c_double)     :: navg1, navg2
        integer(c_int64_t) :: j1, j2, kcells, k
        integer(c_int64_t) :: delta(3), cell1(3), cell2(3), cstart(3), cstop(3)
        integer(c_int)     :: tid, fu, fp1, fp2, iostat, dense
        character(256)     :: fn

        error_code = 1
        if ( product(box%gridsize) /= ncells ) stop "mismatch in gridsize and ncells"
        
        ! Average cell counts
        navg1 = dble(npts1) / ncells ! for set-1
        navg2 = dble(npts2) / ncells ! for set-2
        
        ! Specify the number of neighbours to check on each direction. If the
        ! 3D cell index is`c[i]` for axis i, then the neighbour cell indices 
        ! runs from `c[i] - delta[i]` to `c[i] + delta[i]` (BC applied). 
        delta = ceiling( rmax / box%cellsize ) 

       !$OMP  PARALLEL DEFAULT(SHARED) &
       !$OMP& PRIVATE(tid, j1, j2, k, cell1, cell2, cstart, cstop, dense, fu, fn)

        tid = omp_get_thread_num() + 1 ! thread ID 
        fu  = 10 + tid ! file unit for this thread
        write( fn, '(i0,".",i0,".cpstack.tmp")' ) pid, tid ! filename for this thread
        open(newunit=fu, file=fn, access='stream', form='unformatted',      &
             convert='little_endian', status='unknown', position='append',  &
             action='write'                                                 &
        )

       !$OMP DO SCHEDULE(static)
        do j1 = 1, ncells

            if ( grid_info1(j1)%count < 1 ) cycle ! cell is empty
            
            ! Converting the flat cell index to 3D index
            cell1 = to_3d_index( j1, box%gridsize )

            ! Find the range of indices for the neighbouring cells. Periodic 
            ! wrapping is used if using a periodic BC. Otherwise, the values 
            ! are clipped to the range [0, gridsize-1].  
            cstart = cell1 - delta; call apply_bc(cstart, box%gridsize, periodic)
            cstop  = cell1 + delta; call apply_bc(cstop , box%gridsize, periodic)

            if ( grid_info1(j1)%count > F*navg1 ) then
                ! This cell contains much larger items, than an average cell have.
                ! To avoid any thread get stuck here, all cell pairs including this 
                ! one are kept for special handling.
                dense = 1_c_int
            else 
                dense = 0_c_int
            end if

            ! Walking through the neighbouring cells to make pairs:
            kcells = product(cstop - cstart + 1) ! number of neighbouring cells
            cell2  = cstart
            do k = 1, kcells

                ! Converting the 3D cell index to a flat index
                j2 = to_flat_index( cell2, box%gridsize )

                if ( grid_info2(j2)%count > 0 ) then 
                    if (( dense == 1 ) .or. ( grid_info2(j2)%count > F*navg2 )) then
                        ! Either this cell or the main cell is over-populated. 
                        write(fu) j1, j2, 1_c_int  ! cell pair with spacial handleing
                    else
                        write(fu) j1, j2, 0_c_int ! normal cell pair
                    end if
                end if
                
                ! Incrementing the cell indices
                cell2(1) = cell2(1) + 1
                if ( cell2(1) > cstop(1) ) then
                    cell2(1) = cstart(1)
                    cell2(2) =  cell2(2) + 1
                    if ( cell2(2) > cstop(2) ) then
                        cell2(2) = cstart(2)
                        cell2(3) =  cell2(3) + 1
                    end if
                end if

            end do

        end do
       !$OMP END DO

        close(fu)

       !$OMP END PARALLEL

        ! File for normal cell pair stack 
        fp1 = 8
        write( fn, '(i0,".cpstack.bin")' ) pid 
        open(newunit=fp1, file=fn, access='stream', form='unformatted',     &
             convert='little_endian', status='replace', position='append',  &
             action='write'                                                 &
        )
        
        ! File for special cell pair stack
        fp2 = 9
        write( fn, '(i0,".cpstack.spl.bin")' ) pid 
        open(newunit=fp2, file=fn, access='stream', form='unformatted',     &
             convert='little_endian', status='replace', position='append',  &
             action='write'                                                 &
        )
        
        ! Combining the individual stacks from threads
        do tid = 1, nthreads
            fu  = 10 + tid ! file unit for this thread
            write( fn, '(i0,".",i0,".cpstack.tmp")' ) pid, tid ! filename for this thread
            open(newunit=fu, file=fn, access='stream', form='unformatted', &
                 convert='little_endian', status='old', action='read'      &
            )
            do 
                read(fu, iostat=iostat) j1, j2, dense
                if ( iostat /= 0 ) exit ! EOF or error
                if ( dense == 1 ) then
                    write(fp2) j1, j2 ! to special (dense) pair stack
                else
                    write(fp1) j1, j2 ! to normal pair stack
                end if                 
            end do
            close(fu, status='delete') ! file is deleted on close 
        end do

        close(fp1)
        close(fp2)
        
        error_code = 0
        
    end subroutine enumerate_cell_pairs

! Helper functions

    function to_3d_index(j, gridsize) result(cell)
        !! Calculate the 3D cell indices from a flat index. This 3D index 
        !! follows 0-based indexing, but the flattened index is 1-based.
        integer(c_int64_t), intent(in)  :: gridsize(3), j
        integer(c_int64_t) :: cell(3)
    
        cell(1) = modulo( j-1, gridsize(1) ); cell(3) = (j-1) / gridsize(1) 
        cell(2) = modulo( cell(3), gridsize(2) )
        cell(3) = cell(3) / gridsize(2)
        
    end function to_3d_index
    
    function to_flat_index(cell, gridsize) result(j)
        !! Calculate the flat index from 3D cell indices. This 3D index 
        !! follows 0-based indexing, but the flattened index is 1-based.
        integer(c_int64_t), intent(in) :: gridsize(3), cell(3)
        integer(c_int64_t) :: j
    
        j = cell(1) + gridsize(1)*( cell(2) + gridsize(2)*cell(3) ) + 1

    end function to_flat_index

    subroutine apply_bc(cell, gridsize, periodic)
        !! Apply a boundary condition to the 3D cell index. This 3D index 
        !! follows 0-based indexing.
        integer(c_int)    , intent(in)    :: periodic
        integer(c_int64_t), intent(in)    :: gridsize(3)
        integer(c_int64_t), intent(inout) :: cell(3)
    
        if ( periodic /= 0 ) then
            ! Using periodic boundary conditions: wrap around if the cell  
            ! index is outside range. 
            cell = modulo( cell, gridsize )
        else
            ! No periodic wrapping: clip the value between 0 and gridsize
            cell = min( max( 0_c_int64_t, cell ), gridsize-1 )
        end if
        
    end subroutine apply_bc
    
end module pair_utils_mod