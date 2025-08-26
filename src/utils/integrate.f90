module integrate_mod
    use iso_c_binding
    implicit none

    private
    public :: integrate, leggauss

    integer, parameter :: dp = c_double

    real(c_double), parameter :: K15(16) = [        &
    !   nodes               , weights               &
        0.000000000000000_dp, 0.209482141084728_dp, & !*
        0.207784955007898_dp, 0.204432940075298_dp, &
        0.405845151377397_dp, 0.190350578064785_dp, & !*
        0.586087235467691_dp, 0.169004726639267_dp, &
        0.741531185599394_dp, 0.140653259715525_dp, & !*
        0.864864423359769_dp, 0.104790010322250_dp, &
        0.949107912342759_dp, 0.063092092629979_dp, & !*
        0.991455371120813_dp, 0.022935322010529_dp  &
    ]
    !! Nodes and weigths for Kronrod-15 rule. Exact for polynomials upto 
    !! degree 29. Points are symmetric about x=0.
    
    real(c_double), public, parameter :: G7(8) = [  &
    !   nodes               , weights               &
        0.000000000000000_dp, 0.417959183673469_dp, & !*
        0.405845151377397_dp, 0.381830050505119_dp, & !*
        0.741531185599394_dp, 0.279705391489277_dp, & !*
        0.949107912342759_dp, 0.129484966168870_dp  & !*
    ]
    !! Nodes and weights of for Gauss-7 rule. Exact for polynomials upto 
    !! degree 13. Points are symmetric about x=0 and subset of Kronrod-15
    !! rule.

    type, bind(c) :: heapinfo_t
        integer(c_int64_t) :: size     !! Current size of the heap array
        integer(c_int64_t) :: itemsize !! Size of an item in the heap
        integer(c_int64_t) :: capacity !! Maximum size of the heap
    end type
    
contains

    !============================================================================!
    !         Integration of scalar function y = f(x) over [a, b]                !
    !============================================================================!

    subroutine integrate(f, a, b, args, abstol, reltol, maxiter, res, err, stat) bind(c)
        !! Calculate the integral of a scalar function f(x) over the interval [a, b]. 

        interface
            function f(x, args_) result(y) bind(c)
                !! Function to integrate
                import :: c_double, c_ptr
                real(c_double), value :: x
                type(c_ptr), value :: args_ 
                real(c_double) :: y
            end function
        end interface

        real(c_double), intent(in), value :: a
        !! Lower limit of integration

        real(c_double), intent(in), value :: b
        !! Upper limit of integration

        type(c_ptr), value :: args
        ! Other arguments to pass to the function

        real(c_double), intent(in), value :: abstol
        !! Absolute tolerance

        real(c_double), intent(in), value :: reltol
        !! Relative tolerance

        integer(c_int64_t), intent(in), value :: maxiter
        !! Maximum number of iterations for calculating integral

        real(c_double), intent(out) :: res
        !! Value of the integral of f over [a, b]

        real(c_double), intent(out) :: err
        !! Estimate of the error in integration

        integer(c_int), intent(out) :: stat
        !! Error code: 0=ok, 1=integral not converged

        integer(c_int64_t) :: iter
        real(c_double)     :: xa, xb, xm, I0, I1, I2, err0, err1, err2
        real(c_double)     :: idata(4) ! Interval data: [ -error, left, right, integral ]
        type(heapinfo_t)   :: hinfo
        real(c_double), allocatable :: heap(:, :)

        hinfo%size     = 0
        hinfo%itemsize = 4
        hinfo%capacity = 10*maxiter ! Heap capacity
        allocate( heap(hinfo%itemsize, hinfo%capacity) )

        ! Initial evaluation
        call G7K15(f, a, b, args, I0, err0)

        idata(:) = [ -err0, a, b, I0 ]
        call heap_push(heap, hinfo, idata)

        res  = I0
        err  = err0
        stat = 1
        do iter = 1, maxiter

            ! Stop if tolerance is met
            if ( err <= max(abstol, reltol*abs(res)) ) then
                stat = 0 
                exit
            endif

            ! Pop worst interval
            call heap_pop(heap, hinfo, idata)
            err0 = -idata(1)
            xa   =  idata(2)
            xb   =  idata(3)
            I0   =  idata(4)

            ! Refine it
            xm = 0.5_dp * (xa + xb)
            call G7K15(f, xa, xm, args, I1, err1)
            call G7K15(f, xm, xb, args, I2, err2)

            ! Update global sums
            res = res + (I1   + I2   - I0  ) ! replace old interval
            err = err + (err1 + err2 - err0)

            ! Push new intervals back
            idata(:) = [ -err1, xa, xm, I1 ]
            call heap_push(heap, hinfo, idata)
            idata(:) = [ -err2, xm, xb, I2 ]
            call heap_push(heap, hinfo, idata)
            
        end do

        deallocate(heap)
        
    end subroutine integrate

    subroutine G7K15(f, a, b, args, res, err) bind(c)
        !! Calculate the integral of f(x) over the interval [a, b] using G7-K15 rule.
        !! Value using K15 is returned and error is calculated with G7 rule.

        interface
            function f(x, args_) result(y) bind(c)
                !! Function to integrate
                import :: c_double, c_ptr
                real(c_double), value :: x
                type(c_ptr), value :: args_ 
                real(c_double) :: y
            end function
        end interface

        real(c_double), intent(in) :: a
        !! Lower limit of integration

        real(c_double), intent(in) :: b
        !! Upper limit of integration

        type(c_ptr), value :: args
        ! Other arguments to pass to the function

        real(c_double), intent(out) :: res
        !! Value of the integral of f over [a, b]

        real(c_double), intent(out) :: err
        !! Estimate of the error in integration

        real(c_double) :: xi, fx, scale
        integer(c_int) :: i

        scale = 0.5_dp * ( b - a )
        xi    = scale*( 1. + K15(1) ) + a
        fx    = f(xi, args)
        res   = K15(2) * fx ! store the value of integral using K15 rule
        err   =  G7(2) * fx ! store the value of integral using G7  rule
        do i = 2, 8
            xi  = scale*( 1. - K15(2*i-1) ) + a ! Left point
            fx  = f(xi, args)
            xi  = scale*( 1. + K15(2*i-1) ) + a ! Right point
            fx  = f(xi, args) + fx
            res = res + K15(2*i) * fx
            if ( mod(i, 2) == 1 ) err = err + G7(i+1) * fx ! Point in G7 rule
        end do
        res = scale * res            ! value of integral using K15 rule
        err = abs(res - scale * err) ! difference between G7 and K15 results
        
    end subroutine G7K15

    subroutine leggauss(n, x, w) bind(c)
        !! Generate Gauss-Legendre quadrature rule of order N for [-1, 1].

        integer(c_int), intent(in), value  :: n
        !! Order of the rule: number of nodes.

        real(c_double), intent(out) :: x(n)
        !! Nodes 

        real(c_double), intent(out) :: w(n)
        !! Weights

        real(c_double), parameter :: pi  = 3.141592653589793_c_double

        real(c_double) :: xj, xjo, pm, pn, ptmp
        integer(c_int) :: j, k

        ! -- For odd order, x = 0 is the {floor(n/2)+1}-th node
        if ( modulo(n, 2) == 1 ) then
            
            ! Calculating legendre polynomial P_n(0) using its reccurence relation
            xj = 0._c_double
            pm = 0._c_double
            pn = 1._c_double
            do k = 0, n-1
                ptmp = -k*pm / (k + 1._c_double)
                pm   = pn
                pn   = ptmp
            end do
        
            x(n/2 + 1) = 0._c_double
            w(n/2 + 1) = 2._c_double / (n*pm)**2 ! weight 
        end if

        ! -- Other nodes
        do j = 1, n/2

            ! Initial guess for j-th node (j-th root of n-th legendre polynomial)
            xj  = cos( (2*j - 0.5_c_double) * PI / (2*n + 1._c_double) ) 
            xjo = 100._c_double 
            do while ( abs(xj - xjo) > 1e-08_c_double )

                ! Calculating lenegedre polynomial P_n(xj) using its reccurence relation
                pm = 0._c_double
                pn = 1._c_double
                do k = 0, n-1
                    ptmp = ( (2*k + 1)*xj*pn - k*pm ) / (k + 1._c_double)
                    pm   = pn
                    pn   = ptmp
                end do
                
                ! Next estimate of the root
                xjo = xj
                xj  = xj - pn * (xj**2 - 1) / (n*xj*pn - n*pm)

            end do
            x(j)     = -xj
            w(j)     =  2*(1 - xj**2) / (n*xj*pn - n*pm)**2 !! weight for j-th node
            x(n-j+1) =  xj
            w(n-j+1) =  w(j)
            
        end do
    
    end subroutine leggauss

    !============================================================================!
    !                         Helper functions (private)                         !
    !============================================================================!

    subroutine heap_push(heap, info, item)
        !! Push an item into the heap.

        real(c_double), intent(inout) :: heap(:,:)
        !! Heap array

        type(heapinfo_t), intent(inout) :: info
        !! Details of the heap

        real(c_double), intent(in) :: item(:)
        !! Values to insert to the heap

        integer(c_int64_t) :: i, p
        real(c_double)     :: temp(info%itemsize)

        if ( info%size >= info%capacity ) stop 'error: heap full'
        info%size = info%size + 1
        i         = info%size
        heap(:,i) = item(:)

        ! Bubble up
        do while (i > 1)
            p = i / 2
            if (heap(1,i) < heap(1,p)) then
                temp(:)   = heap(:,i)
                heap(:,i) = heap(:,p)
                heap(:,p) = temp(:)
                i = p
            else
                exit
            end if
        end do
        
    end subroutine heap_push

    subroutine heap_pop(heap, info, item)
        !! Pop an item with lowet key from the heap. First item in each heap 
        !! item is used as the key.

        real(c_double), intent(inout) :: heap(:,:)
        !! Heap array

        type(heapinfo_t), intent(inout) :: info
        !! Details of the heap

        real(c_double), intent(out) :: item(:)
        !! Values to insert to the heap

        integer(c_int64_t) :: i, c
        real(c_double)     :: temp(info%itemsize)

        if (info%size <= 0) stop "error: heap empty"
        item(:)   = heap(:,1)         ! return min element
        heap(:,1) = heap(:,info%size) ! move last to root
        info%size = info%size - 1

        ! Bubble down
        i = 1
        do while (2*i <= info%size)
            c = 2*i
            if (c < info%size .and. heap(1,c+1) < heap(1,c)) c = c + 1
            if (heap(1,c) < heap(1,i)) then
                temp(:)   = heap(:,i)
                heap(:,i) = heap(:,c)
                heap(:,c) = temp(:)
                i = c
            else
                exit
            end if
        end do
        
    end subroutine heap_pop

end module integrate_mod

! program main
!     use iso_c_binding
!     use integrate_mod
!     implicit none
! 
!     real(c_double)    , parameter :: pi  = 3.141592653589793d0
!     real(c_double)    , parameter :: eps = 1.0d-08
!     integer(c_int64_t), parameter :: n   = 50
! 
!     real(c_double) :: res, err
!     integer(c_int) :: ier
!     type(c_ptr)    :: ctx
!   
!     ctx = c_null_ptr
!     call integrate(myfun, -10.0d0, 10.0d0, ctx, eps, eps, n, res, err, ier)
!     print *, 'Result =', res, ' Error =', err, ' Status:', ier
!     print *, sqrt(pi)
! 
! contains
!     function myfun(x, ctx) result(val) bind(c)
!         use iso_c_binding
!         real(c_double), value :: x
!         type(c_ptr), value :: ctx
!         real(c_double) :: val
!         val = exp(-x**2)
!     end function myfun
! 
! end program main


