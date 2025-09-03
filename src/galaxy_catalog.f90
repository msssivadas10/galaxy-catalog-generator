module galaxy_catalog_mod
    !! An extension module to `halo_model_mod` module, containing methods
    !! for generating galaxy catalogs from a given halo catalog.
    
    use omp_lib
    use iso_c_binding
    use halo_model_mod
    use interpolate_mod
    use random_mod
    implicit none

    private
    public :: generate_galaxy_catalog
    
    type, public, bind(c) :: halodata_t
        !! A struct containing data about a halo
        integer(c_int64_t) :: id     !! Unique halo ID
        real(c_double)     :: pos(3) !! Halo position coordinates in Mpc
        real(c_double)     :: mass   !! Halo mass in Msun
    end type halodata_t
        
contains

    subroutine generate_galaxy_catalog(hcat, nhalos, sigma, ns, hmargs, &
                                       bbox, seed, fid, nthreads        &
        ) bind(c)
        !! Generate a galaxy catalog using the given halo catalog. 
        !! 
        !! The catalog will be saved to the file `{fid}.dat` in the working 
        !! directory. This binary file contains a string of records for each 
        !! galaxy, where each record consists of a parent halo ID (`int64`), 
        !! position coordinates (`float64[3]`), mass (`float64`) and typecode 
        !! (`char`, `'c'` for central and `'s'` for satellite). Little endian 
        !! byte order is used for portability.

        type(halodata_t), intent(in) :: hcat(nhalos)
        !! Halo catalog

        integer(c_int64_t), intent(in), value :: nhalos
        !! Number of halos 

        real(c_double), intent(in) :: sigma(3,ns)
        !! Natural spline data for calculating matter variance. This should 
        !! give the value of `log(sigma)=log(variance)/2` as a function of
        !! the natural log of halo mass in Msun.

        integer(c_int64_t), intent(in), value :: ns
        !! Size of the variance spline

        type(hmargs_t), intent(in) :: hmargs
        !! Halo model parameters

        real(c_double), intent(in) :: bbox(3,2)
        !! Bounding box for the halo positions. First row specifies the 
        !! minimum and second row maximum bound. 

        integer(c_int64_t), intent(in), value :: seed
        !! Seed value for random number generators

        integer(c_int64_t), intent(in), value :: fid
        !! A unique ID for generating output the filename storing galaxy 
        !! catalog data. Given a positive integer, the filename would be
        !! `{fid}.dat`. 
        ! 
        ! NOTE: This is a simple work-around to bypass the difficulties in 
        ! passing a string filename to this subroutine, when calling from 
        ! external C or python codes. In that case, a wrapper should be 
        ! written to handle the C string and convert it to fortran string.
        ! Using an ID and a file naming logic, the need for a seperate  
        ! wrapper can be dropped, and the user can do anything with the 
        ! generated file later. :)

        integer(c_int), intent(in), value :: nthreads
        !! Number of threads to use

        character(len=256) :: fn, tfn
        integer(c_int)     :: tid, fu, fo, ierr
        integer(c_int64_t) :: i, gbuf_size, ng, halo_id
        type(cgargs_t)     :: args 
        real(c_double)     :: gdata(4)
        real(c_double), allocatable :: gbuf(:,:) ! Galxy position and mass 

        ! Set number of threads
        call omp_set_num_threads(nthreads) 

        !$OMP PARALLEL PRIVATE(tid, tfn, fu, args, gbuf_size, gbuf, ng)
        
        tid = omp_get_thread_num() + 1 ! thread ID
        
        ! Opening a private temporary file for writing data from this thread. 
        ! These files will have a specific filename and unit ID based on the 
        ! thread ID. Data is saved in little endian binary format to avoid
        ! loss of precision.
        fu = 10 + tid ! file unit for this thread
        write(tfn, '(i0,".",i0,".tmp")') fid, tid ! filename for this thread
        open(newunit=fu, file=tfn, access='stream', form='unformatted',    &
             convert='little_endian', status='unknown', position='append', &
             action='write'                                                &
        )

        ! Initialising the random number generator. This RNG works on a state 
        ! private to the thread, with a seed offset by the main seed. So, the 
        ! generated RVs should be a different sequence on each thread.
        call pcg32_init(args%rstate, seed + 1000*tid) 

        ! Allocate galaxy data table: At first, an array that can hold a maximum 
        ! of 1024 galaxies are allocated on each thread. This is reallocated if 
        ! needed.
        gbuf_size = 1024
        allocate( gbuf(4,gbuf_size) )

        ! These are same for all halos
        args%boxsize(1:3) = bbox(1:3,2) - bbox(1:3,1) ! Boxsize
        args%offset(1:3)  = bbox(1:3,1) ! Offset or bottom-lower-left corner coordinates

        !$OMP DO SCHEDULE(static)
        do i = 1, nhalos

            ! Copy halo data to local args
            args%pos(1:3) = hcat(i)%pos(1:3)    ! position
            args%lnm      = log( hcat(i)%mass ) ! mass
            args%s        = exp( interpolate(args%lnm, ns, sigma) ) ! matter variance

            ! Setting up 
            call setup_catalog_generation(hmargs, args)

            ng = args%n_cen + args%n_sat ! total number of galaxies in this halo 
            if ( ng < 1 ) cycle

            ! Ensure there is enough space for storing all the expected galaxies. 
            ! Galaxy buffer is resized if needed.
            call ensure_capacity_thread(gbuf, gbuf_size, ng, 4_c_int64_t)
            
            ! Generating the galaxies
            call generate_galaxies( hmargs, args, ng, gbuf(:,1:ng) )

            ! Saving the galaxy data to the thread specific output file.
            write(fu) hcat(i)%id   ! halo unique ID
            write(fu) ng           ! number of galaxies in this block
            write(fu) gbuf(:,1:ng) ! galaxy data

        end do
        !$OMP END DO

        deallocate( gbuf )
        close(fu) ! closing the thread specific temp file

        !$OMP END PARALLEL 

        ! Merge data from all the temporary files to the specified output file
        ! in the correct format. That is, a binary file with each galaxy entry 
        ! corresponds to a parent halo ID (int64), galaxy type (character `C`
        ! for central and `S` for satellite), position coordinates and mass 
        ! (all float64).
        fu = 10 ! file unit for temporary outputs 
        fo = 11 ! file unit for main output
        write(fn, '(i0,".dat")') fid ! filename for main output
        open(newunit=fo, file=fn, access='stream', form='unformatted', &
             convert='little_endian', status='replace', action='write' &
        ) ! main output file 

        do tid = 1, nthreads
            write(tfn, '(i0,".",i0,".tmp")') fid, tid ! filename for this thread
            open(newunit=fu, file=tfn, access='stream', form='unformatted', &
                 convert='little_endian', status='old', action='read'       &
            ) ! temporary file

            do 
                ! Read halo ID and number of galaxies associated with this halo 
                read(fu, iostat=ierr) halo_id, ng 
                if ( ierr /= 0 ) exit ! end of file
                
                ! Load central galaxy data: always the first item in a block 
                ! corresponding to a halo. This is marked by the value 'c' in
                ! the output file. 
                read(fu, iostat=ierr) gdata(1:4)
                if ( ierr /= 0 ) exit
                write(fo) halo_id, gdata(1:4), 'c' 
                    
                ! Load satellite galaxy data. This is marked by the value 's'
                ! in the output file.
                do i = 2, ng
                    read(fu, iostat=ierr) gdata(1:4)
                    if ( ierr /= 0 ) exit
                    write(fo) halo_id, gdata(1:4), 's'
                end do
            end do

            ! Close temporary file: this will also delete the file
            close(fu, status='delete') 

        end do
        close(fo) ! close main output file
        
    end subroutine generate_galaxy_catalog

    subroutine ensure_capacity_thread(arr, current_size, needed_size, ncols)
        !! Safely grow a 2D buffer [ncols, capacity] for one thread

        real(c_double), allocatable, intent(inout) :: arr(:,:) ! shape (ncols, capacity)
        integer(c_int64_t), intent(inout) :: current_size
        integer(c_int64_t), intent(in)    :: needed_size, ncols
        real(c_double), allocatable :: tmp(:,:)
    
        if (needed_size > current_size) then
            ! Grow size exponentially or to needed size
            current_size = max(2*current_size, needed_size)
    
            allocate(tmp(ncols, current_size))
            if (allocated(arr)) tmp(:,1:size(arr,2)) = arr
            call move_alloc(tmp, arr)
        end if

    end subroutine ensure_capacity_thread
    
end module galaxy_catalog_mod