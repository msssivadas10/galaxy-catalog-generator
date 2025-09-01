module galaxy_catalog_mod
    
    use omp_lib
    use iso_c_binding
    ! use halo_model_mod
    ! use interpolate_mod
    implicit none
    
    type, public, bind(c) :: halodata_t
        !! Struct containing data about a halo
        integer(c_int64_t) :: id     !! Unique halo ID
        real(c_double)     :: pos(3) !! Halo position coordinates in Mpc
        real(c_double)     :: mass   !! Halo mass in Msun
    end type halodata_t

    type, public, bind(c) :: galdata_t
        !! Struct containing data about a galaxy
        integer(c_int64_t) :: parent_id !! Unique ID of the parent halo
        real(c_double)     :: pos(3)    !! Galaxy position coordinates in Mpc
        real(c_double)     :: mass      !! Galaxy mass in Msun 
        integer(c_int8_t)  :: gtype     !! Galaxy type: 1=central, 2=satellite 
    end type galdata_t
        
contains
    
end module galaxy_catalog_mod