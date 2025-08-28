module integrate_mod
    !! Routines for numerical integration.

    use iso_c_binding
    implicit none

    private
    public :: integrate, leggauss

    integer, parameter :: dp = c_double

    integer, parameter :: GKSIZE = 8
    !! Size of the Gauss-Kronrod rule array 

    real(c_double), parameter :: K15(2, 8) = reshape([ &
    !   nodes               , weights               
        0.000000000000000_dp, 0.209482141084728_dp, & !*
        0.207784955007898_dp, 0.204432940075298_dp, &
        0.405845151377397_dp, 0.190350578064785_dp, & !*
        0.586087235467691_dp, 0.169004726639267_dp, &
        0.741531185599394_dp, 0.140653259715525_dp, & !*
        0.864864423359769_dp, 0.104790010322250_dp, &
        0.949107912342759_dp, 0.063092092629979_dp, & !*
        0.991455371120813_dp, 0.022935322010529_dp  &
    ], shape(K15))
    !! Nodes and weigths for Kronrod-15 rule. Exact for polynomials upto 
    !! degree 29. Points are symmetric about x=0.
    
    real(c_double), parameter :: G7(2, 4) = reshape([ &
    !   nodes               , weights               
        0.000000000000000_dp, 0.417959183673469_dp, & !*
        0.405845151377397_dp, 0.381830050505119_dp, & !*
        0.741531185599394_dp, 0.279705391489277_dp, & !*
        0.949107912342759_dp, 0.129484966168870_dp  & !*
    ], shape(G7))
    !! Nodes and weights of for Gauss-7 rule. Exact for polynomials upto 
    !! degree 13. Points are symmetric about x=0 and subset of Kronrod-15
    !! rule.

contains

    subroutine integrate(f, a, b, args, abstol, reltol, maxiter, res, err, stat) bind(c)
        !! Calculate the integral of a scalar function f(x) over the interval [a, b]. 

        interface
            function f(x, args_) result(y) bind(c)
                !! Function to integrate
                import :: c_double, c_ptr
                real(c_double), value :: x
                type(c_ptr)   , value :: args_ 
                real(c_double) :: y
            end function
        end interface

        real(c_double), intent(in), value :: a
        !! Lower limit of integration

        real(c_double), intent(in), value :: b
        !! Upper limit of integration

        type(c_ptr), value :: args
        !! Other arguments to pass to the function

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
        
        integer(c_int64_t) :: iter, j
        real(c_double)     :: intg, intk, fval, scale
        real(c_double)     :: xa, xb, xm, I0, I1, I2, err0, err1, err2
        
        integer(c_int64_t) :: heap_size, heap_capacity
        real(c_double), allocatable :: heap(:, :)

        heap_size     = 0
        heap_capacity = 10*maxiter ! Heap capacity
        allocate( heap(4, heap_capacity) )

        ! Initial evaluation
        scale = 0.5_dp * (b - a)
        fval  = f(a + scale, args)
        intk  = fval * K15(2,1) 
        intg  = fval *  G7(2,1) 
        do j = 2, GKSIZE
            fval = f(a + scale * (1. - K15(1,j)) , args) + f(a + scale * (1. + K15(1,j)) , args)
            intk = intk + fval * K15(2,j)
            if ( mod(j, 2) == 1 ) intg = intg + fval * G7(2,(j+1)/2) ! Point also in G7 rule
        end do
        intk = scale * intk
        intg = scale * intg
        I0   = intk
        err0 = abs(intk - intg)
        call int_heap_push(heap, heap_size, heap_capacity, a, b, I0, err)

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
            call int_heap_pop(heap, heap_size, heap_capacity, xa, xb, I0, err0)
            
            xm = 0.5_dp * (xa + xb) 
            
            ! Refine on left interval
            scale = 0.5_dp * (xm - xa)
            fval  = f(xa + scale, args)
            intk  = fval * K15(2,1) 
            intg  = fval *  G7(2,1) 
            do j = 2, GKSIZE
                fval = f(xa + scale * (1. - K15(1,j)) , args) + f(xa + scale * (1. + K15(1,j)) , args)
                intk = intk + fval * K15(2,j)
                if ( mod(j, 2) == 1 ) intg = intg + fval * G7(2,(j+1)/2) ! Point also in G7 rule
            end do
            intk = scale * intk
            intg = scale * intg
            I1   = intk
            err1 = abs(intk - intg)
            call int_heap_push(heap, heap_size, heap_capacity, xa, xm, I1, err1) ! Push new interval back
            
            ! Refine on left interval
            scale = 0.5_dp * (xb - xm)
            fval  = f(xm + scale, args)
            intk  = fval * K15(2,1) 
            intg  = fval *  G7(2,1) 
            do j = 2, GKSIZE
                fval = f(xm + scale * (1. - K15(1,j)) , args) + f(xm + scale * (1. + K15(1,j)) , args)
                intk = intk + fval * K15(2,j)
                if ( mod(j, 2) == 1 ) intg = intg + fval * G7(2,(j+1)/2) ! Point also in G7 rule
            end do
            intk = scale * intk
            intg = scale * intg
            I2   = intk
            err2 = abs(intk - intg)
            call int_heap_push(heap, heap_size, heap_capacity, xm, xb, I2, err2) ! Push new interval back
            
            ! Update global sums
            res = res + (I1   + I2   - I0  ) ! replace old interval
            err = err + (err1 + err2 - err0)
            
        end do

        deallocate(heap)
        
    end subroutine integrate

    subroutine leggauss(n, x, w) bind(c)
        !! Generate Gauss-Legendre quadrature rule of order N for [-1, 1].

        integer(c_int), intent(in), value :: n    
        !! Order (number of nodes)
        
        real(c_double), intent(out) :: x(n) 
        !! Nodes 
        
        real(c_double), intent(out) :: w(n) 
        !! Weights

        real(c_double), parameter :: PI  = 3.141592653589793_dp
        !! Pi

        real(c_double) :: xj, xjo, pm, pn, ptmp
        integer(c_int) :: j, k

        ! If order is odd number, x = 0 is a node
        if ( modulo(n, 2) == 1 ) then
            ! Calculating legendre polynomial P_n(0) using its reccurence relation
            xj = 0.0_dp
            pm = 0.0_dp
            pn = 1.0_dp
            do k = 0, n-1
                ptmp = -k*pm / (k + 1.0_dp)
                pm   = pn
                pn   = ptmp
            end do
            x(n/2 + 1) = 0.0_dp
            w(n/2 + 1) = 2.0_dp / (n*pm)**2 ! weight 
        end if

        ! Other nodes (roots of the n-th legendre polynomial)
        do j = 1, n/2
            xj  = cos( (2*j - 0.5_dp)*PI / (2*n + 1.0_dp) ) ! initial guess for the root
            xjo = 100.0_dp ! or any large number > 1
            do while ( abs(xj - xjo) > 1e-08_dp )
                ! Calculating lenegedre polynomial P_n(xj) using its reccurence relation
                pm = 0.0_dp
                pn = 1.0_dp
                do k = 0, n-1
                    ptmp = ( (2*k + 1)*xj*pn - k*pm ) / (k + 1.0_dp)
                    pm   = pn
                    pn   = ptmp
                end do
                xjo = xj
                xj  = xj - pn * (xj**2 - 1) / (n*xj*pn - n*pm)
            end do
            x(j)     = -xj
            w(j)     =  2*(1 - xj**2) / (n*xj*pn - n*pm)**2 ! weight
            x(n-j+1) =  xj
            w(n-j+1) =  w(j)            
        end do
    
    end subroutine leggauss

    ! Helper functions (private):

    subroutine int_heap_push(heap, size, capacity, a, b, val, err)
        !! Push the integration result on a interval [a, b] to the interval heap.

        real(c_double), intent(inout) :: heap(4, capacity) 
        !! Heap array
        
        integer(c_int64_t), intent(inout) :: size
        !! Current size of the heap

        integer(c_int64_t), intent(in) :: capacity
        !! Maximum capacity of the heap

        real(c_double), intent(in) :: a
        !! Left end value of the interval

        real(c_double), intent(in) :: b
        !! Right end value of the interval

        real(c_double), intent(in) :: val
        !! Value of the integral in the interval

        real(c_double), intent(in) :: err
        !! Error estimate

        integer(c_int64_t) :: i, p
        real(c_double)     :: temp(4)

        if ( size >= capacity ) stop 'error: heap full'
        size       = size + 1
        i          = size
        heap(:, i) = [ -err, a, b, val ]

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
        
    end subroutine int_heap_push

    subroutine int_heap_pop(heap, size, capacity, a, b, val, err)
        !! Pop the interval with largest error from the heap.

        real(c_double), intent(inout) :: heap(4, capacity) 
        !! Heap array
        
        integer(c_int64_t), intent(inout) :: size
        !! Current size of the heap

        integer(c_int64_t), intent(in) :: capacity
        !! Maximum capacity of the heap

        real(c_double), intent(out) :: a
        !! Left end value of the interval

        real(c_double), intent(out) :: b
        !! Right end value of the interval

        real(c_double), intent(out) :: val
        !! Value of the integral in the interval

        real(c_double), intent(out) :: err
        !! Error estimate

        integer(c_int64_t) :: i, c
        real(c_double)     :: temp(4)

        if (size <= 0) stop "error: heap empty"
        
        ! Return min element
        err = -heap(1,1) 
        a   =  heap(2,1)
        b   =  heap(3,1)
        val =  heap(4,1)
        
        ! Move last to root
        heap(:,1) = heap(:,size) 
        size      = size - 1

        ! Bubble down
        i = 1
        do while (2*i <= size)
            c = 2*i
            if (c < size .and. heap(1,c+1) < heap(1,c)) c = c + 1
            if (heap(1,c) < heap(1,i)) then
                temp(:)   = heap(:,i)
                heap(:,i) = heap(:,c)
                heap(:,c) = temp(:)
                i = c
            else
                exit
            end if
        end do
        
    end subroutine int_heap_pop

end module integrate_mod

! program main
!     use iso_c_binding
!     use integrate_mod
!     implicit none
!     !
!     real(c_double)    , parameter :: PI  = 3.141592653589793_c_double
!     real(c_double)    , parameter :: eps = 1.0d-08
!     integer(c_int64_t), parameter :: n   = 50
!     real(c_double) :: res, err
!     integer(c_int) :: ier
!     type(c_ptr)    :: ctx
!     !   
!     ctx = c_null_ptr
!     call integrate(myfun, -10.0d0, 10.0d0, ctx, eps, eps, n, res, err, ier)
!     print *, 'Result =', res, ' Error =', err, ' Status:', ier
!     print *, sqrt(PI)
! contains
!     function myfun(x, ctx) result(val) bind(c)
!         use iso_c_binding
!         real(c_double), value :: x
!         type(c_ptr), value :: ctx
!         real(c_double) :: val
!         val = exp(-x**2)
!     end function myfun
! end program main
