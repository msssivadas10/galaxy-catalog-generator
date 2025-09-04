module galaxy_catalog_mod
    !! An extension module to `halo_model_mod` module, containing methods
    !! for generating galaxy catalogs from a given halo catalog.
    
    use omp_lib
    use iso_c_binding
    use constants_mod
    use halo_model_mod
    use interpolate_mod
    use random_mod
    implicit none

    private
    public :: setup_catalog_generation, generate_galaxies, generate_galaxy_catalog
    
    type, public, bind(c) :: halodata_t
        !! A struct containing data about a halo
        integer(c_int64_t) :: id     !! Unique halo ID
        real(c_double)     :: pos(3) !! Halo position coordinates in Mpc
        real(c_double)     :: mass   !! Halo mass in Msun
    end type halodata_t

    type, public, bind(c) :: cgargs_t
        !! Struct containing various parameters for galaxy catalog
        !! generation (execpt halo model parameters). 

        real(c_double)     :: lnm    !! Natural log of halo mass in Msun
        real(c_double)     :: pos(3) !! Coordinates of the halo in Mpc units
        real(c_double)     :: lnr    !! Natural log of halo radius in Mpc
        real(c_double)     :: s      !! Matter variance corresponding to this halo mass
        real(c_double)     :: c      !! Concentration parameter for the halo
        integer(c_int64_t) :: n_cen  !! Number of central galaxies (0 or 1)
        integer(c_int64_t) :: n_sat  !! Number of satellite galaxies
        integer(c_int64_t) :: rstate !! Random number generator state

        real(c_double) :: boxsize(3)
        !! Size of the bounding box containing all the halos in the 
        !! simulation. Used for periodic wrapping galaxy position.

        real(c_double) :: offset(3)
        !! Coordinates of the bottom-lower-left corner of the bounding 
        !! box. Used for periodic wrapping galaxy position.

    end type
    
        
contains

    subroutine setup_catalog_generation(params, args) bind(c)
        !! Calculate various parameters for a galaxy catalog generation

        type(hmargs_t), intent(in) :: params 
        !! Halo model parameters

        type(cgargs_t), intent(inout) :: args
        !! Arguments for catalog generation

        real(c_double) :: rho_m, rho_h, p_cen, lam_sat

        rho_m = params%Om0 * ( critical_density_const * params%H0**2 ) ! Matter density at z=0 in Msun/Mpc^3 
        rho_h = rho_m ! Halo density (TODO: chek if the halo density is rho_m * self.Delta)

        ! Lagrangian radius (r) corresponding to halo mass
        args%lnr = ( args%lnm + log(3._c_double / (4*pi) / rho_h ) ) / 3._c_double ! r in Mpc

        ! Central galaxy count: this is drawn from a binomial distribution of n = 1, 
        ! so that the value is either 1 or 0. If the model uses a step function use 
        ! the average count as the actual count. 
        p_cen = central_count(params, args%lnm)
        if ( abs(params%sigma_m) < 1e-08_c_double ) then
            args%n_cen = int(p_cen, kind=c_int64_t)
        else
            args%n_cen = binomial_rv(args%rstate, 1_c_int64_t, p_cen)
        end if
        if ( args%n_cen < 1 ) return ! No central galaxies -> no satellites also 

        ! Satellite galaxy count: this drawn from a poisson distribution with the 
        ! calculated average.
        lam_sat    = satellite_count(params, args%lnm)
        args%n_sat = poisson_rv(args%rstate, lam_sat)
        
        if ( params%scale_shmf < exp(params%lnm_min - args%lnm) ) then 
            ! Halo mass correspond to an invalid satellite galaxy mass range: no satellites
            args%n_sat = 0
        end if

        ! Calculating halo concentration parameter
        args%c = halo_concentration(params, args%s)

    end subroutine setup_catalog_generation

    subroutine generate_galaxies(params, args, n, gdata) bind(c)
        !! Generate galaxy positions and mass.

        type(hmargs_t), intent(in) :: params 
        !! Halo model parameters

        type(cgargs_t), intent(inout) :: args
        !! Arguments for catalog generation

        integer(c_int64_t), intent(in), value :: n
        !! Size of the position and mass array: must be same total number 
        !! of galaxies i.e., `(args%n_cen+args%n_sat)`.

        real(c_double), intent(out) :: gdata(4,n)
        !! Galaxy positions (columns 1-3) and masses (column 4). 

        integer(c_int64_t) :: i
        real(c_double) :: m_halo, r_halo, c_halo, f, r, theta, phi, Ac, k1, k2, p

        if ( args%n_cen < 1 ) return ! No galaxies in this halo

        m_halo = exp(args%lnm) ! Halo mass in Msun
        r_halo = exp(args%lnr) ! Halo radius in Mpc
        c_halo = args%c        ! Halo concentration parameter

        ! Halo has a central galaxy: the position and mass of this galaxy is same
        ! as that of the parent halo.
        gdata(1:3,1) = args%pos(1:3) ! in Mpc
        gdata(4  ,1) = m_halo        ! in Msun

        if ( args%n_sat < 1 ) return ! No satellite galaxies in this halo
        do i = 1, args%n_sat
            
            ! Assigning random mass values to the satellite galaxies: These masses 
            ! are drown from a bounded pareto distribution, with bounds `m_min` and  
            ! `scale_shmf*m_halo`, and slope given by slope_shmf.
            !  
            ! NOTE: RVs are generated using inverse transform sampling 
            ! (<https://en.wikipedia.org/wiki/Pareto_distribution>)
            p  = -1._c_double / params%slope_shmf
            k1 = exp( (params%lnm_min - args%lnm)*params%slope_shmf )
            k2 = params%scale_shmf**params%slope_shmf
            f  = ( ( k2 - (k2 - k1) * uniform_rv(args%rstate) ) / (k1*k2) )**p ! m_sat / m_halo
            
            ! Generating random values corresponding to the distance of the galaxy 
            ! from the halo center. These RVs should follow a distribution matching 
            ! the NFW density profile of the halo. Sampling is done using the 
            ! inverse transformation method. 
            Ac    = uniform_rv(args%rstate)*( log(1 + c_halo) - c_halo / (1 + c_halo) )
            r     = (r_halo / c_halo) * nfw_c(Ac)
            theta = acos( 2*uniform_rv(args%rstate) - 1 ) ! -pi to pi
            phi   = 2*pi*uniform_rv(args%rstate)          !   0 to 2pi
            
            ! Satellite galaxy coordinates x, y, and z in Mpc
            gdata(1,i+1) = gdata(1,1) + r*sin(theta)*cos(phi)
            gdata(2,i+1) = gdata(2,1) + r*sin(theta)*sin(phi)
            gdata(3,i+1) = gdata(3,1) + r*cos(theta)
            ! Periodic wrapping of coordinates to stay within the bounding box
            call periodic_wrap( gdata(1:3,i+1), args%offset(1:3), args%boxsize(1:3) )
            
            ! Satellite galaxy mass in Msun
            gdata(4,i+1) = gdata(4,1) * f

        end do
        
    end subroutine generate_galaxies

    elemental subroutine periodic_wrap(x, offset, width)
        !! Periodicall wrap the value to the interval [offset, offset+width] 
        real(c_double), intent(inout) :: x
        real(c_double), intent(in)    :: offset, width

        ! If width <= 0, then no wrapping is done: this can be used for forcing 
        ! no wrapping, if needed...
        if ( width > 0. ) then
            ! Using `modulo` (for mathematical modulo) instead of `mod` function
            ! (remainder of a division)
            x = offset + modulo(x - offset, width)
        end if

    end subroutine periodic_wrap

    function nfw_c(a) result(res)
        !! Return the value of inverse to the NFW mass function. 

        real(c_double), intent(in) :: a
        real(c_double) :: res, x, p1, p2

        integer, parameter :: f8 = c_double

        ! Using an approximate piecewise rational function for inverting the 
        ! c-A(c) relation. This fit works good upto A(c) ~ 10, with less that  
        ! 10% error. For getting the actual value, one should invert the NFW 
        ! mass function, A(c) = log(1+c) - c/(1+c).  
        x = log10(a)
        if ( a < 1e-03_c_double ) then
            res = 0.50018962_f8*x + 0.15241388_f8
        else if ( a < 2._c_double ) then
            p1  = 2699.40545133_f8 + 2921.07235917_f8*x - 1162.90566455_f8*x**2 
            p2  = 3705.23701996_f8 - 3065.75405505_f8*x -   61.92662277_f8*x**2
            res = p1 / p2
        else
            p1  = 47.25501938_f8 + 32.98237791_f8*x - 38.26172387_f8*x**2  
            p2  = 66.29326139_f8 - 87.37718415_f8*x + 29.86558497_f8*x**2
            res = p1 / p2
        end if
        res = 10._c_double**res
        
    end function nfw_c

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