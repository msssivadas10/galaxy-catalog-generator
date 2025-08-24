module quadutils
    use iso_c_binding
    implicit none

    private
    public :: integrate, G7K15, make_functable, integrate_functable

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

    subroutine make_functable(f, a, b, args, npts, rule, xy, stat) bind(c)
        !! Make a table of the values of function f(x) in the interval [a, b].
        !! 
        !! This table is used for calculating integral of the function f, which 
        !! can be expensive to compute. This function is evaluated at quadrature 
        !! nodes specified by the rule flag. 
        !!
        !! - rule=0: use Gauss-7 rule on an adaptive subdivision on interval [a, b]
        !! - rule=1: non-uniform adaptive grid, used for non-uniform Simpson's rule
        !! - rule=2: uniform grid for usual Simpson's rule integration    
        
        interface
            function f(x, args_) result(y) bind(c)
                !! Function to sample
                import :: c_double, c_ptr
                real(c_double), value :: x
                type(c_ptr)   , value :: args_ 
                real(c_double) :: y
            end function
        end interface

        real(c_double), intent(in), value :: a
        !! Lower limit 
        
        real(c_double), intent(in), value :: b
        !! Upper limit 

        type(c_ptr), value :: args
        !! Other arguments to pass to the function

        integer(c_int), intent(in), value :: rule
        !! Specify sampling rule (2=uniform, 1=non-uniform, 0=stacked G7).

        integer(c_int64_t), intent(in), value :: npts
        !! Number of points. If `rule=0`, this must be a multiple of 7.

        real(c_double), intent(out) :: xy(2, npts)
        !! Integration nodes and function values

        integer(c_int), intent(out) :: stat
        !! Error code: 0=ok, 1=error

        integer(c_int), parameter :: datsize = 7
        
        integer(c_int64_t) :: i, j, n, splits, idx(npts)
        real(c_double)     :: xa, xb, xm, xma, xmb, ya, yb, ym, yma, ymb, err, scale
        real(c_double)     :: xs(npts), ys(npts)
        real(c_double)     :: idata(7) ! Interval data: [ -error, left, value, right, value, mid, value ]    
        type(heapinfo_t)   :: hinfo
        real(c_double), allocatable :: heap(:, :) ! Priority queue for sub-intervals 
        ! The heap key is the error in the linear approximation of the function in 
        ! that interval. Intervals are divided if error is large.

        stat = 1
        if ( rule == 0 ) then
            if ( modulo(npts, 7) /= 0 ) stop "points must be multiple of 7 for rule=0"
            splits = npts / 7 ! Number of sub-intervals
        else
            splits = npts-1
        end if

        xa    = a
        ya    = f(xa, args)
        xs(1) = xa 
        ys(1) = ya 
        xb    = b
        yb    = f(xb, args)
        xs(2) = xb
        ys(2) = yb
        if ( splits > 1 ) then
            if ( rule == 2 ) then
                ! Sampling splits+1 points uniformly from the interval [a, b]
                
                scale = (b - a) / splits ! Stepsize
                
                xs(splits+1) = xs(2)
                ys(splits+1) = ys(2)
                do i = 2, splits
                    xm    = xs(i-1) + scale
                    xs(i) = xm
                    ys(i) = f(xm, args) 
                end do
            else
                ! Adaptive sampling of the interval into splits+1 points, based on 
                ! how the function changes
                
                hinfo%size     = 0
                hinfo%itemsize = 7
                hinfo%capacity = 10*(splits + 1) ! Heap capacity
                allocate( heap(hinfo%itemsize, hinfo%capacity) )

                xm       = 0.5_dp * ( xa + xb )
                ym       = f(xm, args)
                err      = abs( ym - 0.5_dp * ( ya + yb ) )
                idata(:) = [ -err, xa, ya, xb, yb, xm, ym ]
                call heap_push(heap, hinfo, idata)
                
                i = 2
                n = splits + 1 ! Number of end points
                do while ( i < n .and. hinfo%size > 0 )

                    ! Pop worst interval
                    call heap_pop(heap, hinfo, idata)
                    xa = idata(2)
                    ya = idata(3)
                    xb = idata(4)
                    yb = idata(5)
                    xm = idata(6)
                    ym = idata(7)

                    ! Add midpoint
                    i     = i + 1
                    xs(i) = xm
                    ys(i) = ym
                    
                    ! Subdivide into two new intervals
                    ! -- Left:
                    xma      = 0.5_dp * ( xa + xm )
                    yma      = f(xma, args)
                    err      = abs( yma - 0.5_dp * ( ya + ym ) )
                    idata(:) = [ -err, xa, ya, xm, ym, xma, yma ]
                    call heap_push(heap, hinfo, idata)

                    ! -- Right:
                    xmb      = 0.5_dp * ( xm + xb )
                    ymb      = f(xmb, args)
                    err      = abs( ymb - 0.5_dp * ( ym + yb ) )
                    idata(:) = [ -err, xm, ym, xb, yb, xmb, ymb ]
                    call heap_push(heap, hinfo, idata)

                end do

                deallocate( heap )   
                
                ! Sorting the samples based on x values
                call argsort(xs, idx, n)
                xy(1, 1:n) = xs(1:n) ! Using xy array as temp array for rearranging...
                xy(2, 1:n) = ys(1:n)
                do i = 1, n
                    xs(i) = xy(1, idx(i))
                    ys(i) = xy(2, idx(i))
                end do
                
            end if
        end if

        if ( rule /= 0 ) then
            ! Only the splits+1 points are required
            xy(1,1:npts) = xs(1:npts)
            xy(2,1:npts) = ys(1:npts)
            
            stat = 0
            return
        end if

        ! Calculate the G7 nodes and function values in each sub-intervals
        do i = 1, splits
            xa    = xs(i)
            scale = 0.5_dp * ( xs(i+1) - xa )
            do j = 1, 7
                if ( j > 3 ) then
                    xm = ( 1. + G7(2*j-7) ) * scale + xa
                else
                    xm = ( 1. - G7(9-2*j) ) * scale + xa
                end if
                xy(1, 7*i+j-7) = xm
                xy(2, 7*i+j-7) = f(xm, args)
            end do
        end do
        stat = 0
        
    end subroutine make_functable

    subroutine integrate_functable(x, y, rule, npts, res, stat) bind(c)
        !! Calculate the integral of a function using a table of the values of 
        !! function f(x) in the interval [a, b].

        integer(c_int), intent(in), value :: rule
        !! Integration rule. 0=stacked G7, 1=Simpson's rule with non-uniform spacing,
        !! 2=Simpson'r rule with uniform spacing.

        integer(c_int64_t), intent(in), value :: npts
        !! Number of points. If `rule=0`, this must be a multiple of 7.

        real(c_double), intent(in) :: x(npts)
        !! Integration nodes

        real(c_double), intent(in) :: y(npts)
        !! Function values at integration nodes

        real(c_double), intent(out) :: res
        !! Calculated value of the integral

        integer(c_int), intent(out) :: stat
        !! Error code: 0=ok, 1=error (if rule=0, but used rule=1 because npts % 7 != 0)

        integer(c_int64_t) :: i, nint
        real(c_double)     :: hl, hr, wl, wm, wr, scale

        stat = 1
        res  = 0.0_c_double
        if (( rule == 0 ) .and. ( modulo(npts, 7) == 0 )) then
            ! Integration using composite Gauss-7 rule
            do i = 1, npts, 7
                
                ! Calculate the scale factor for converting [-1, 1] to the subinterval
                ! [a, b], using the transformation x = scale * (1 + xg) + left
                scale = 0.5_c_double * ( x(i+6) - x(i) ) / G7(7) 

                res = res + scale * ( sum( y(i:i+3) * G7(8:2:-2) ) + sum( y(i+4:i+6) * G7(4:8:2) ) )
            end do
            stat = 0
        else
            ! Integration using Simpson's rule for uniform or non-uniform spacing
            
            nint = npts-1 ! Number of intervals
            
            do i = 1, nint-1, 2
                hl = x(i+1) - x(i)
                hr = x(i+2) - x(i+1)
                if ( abs(hr - hl) < 1e-08_c_double ) then
                    ! Use weights for uniform sample version
                    wl = 1._c_double
                    wm = 4._c_double
                    wr = 1._c_double
                else
                    ! Use weights for non-uniform sample version
                    wl = 2._c_double - hr / hl
                    wm = ( hl + hr )**2 / ( hl * hr )
                    wr = 2._c_double - hl / hr
                end if
                res = res + (hl + hr) * ( wl * y(i) + wm * y(i+1) + wr * y(i+2) )
            end do
            if ( modulo(nint, 2) == 1 ) then
                ! For odd number of intervals, the above calculation is done upto the second
                ! last interval and the last interval is handled seperately. 
                ! (see <https://en.wikipedia.org/wiki/Simpson%27s_rule>)

                hl = x(npts-1) - x(npts-2)
                hr = x(npts  ) - x(npts-1)
                if ( abs(hr - hl) < 1e-08_c_double ) then
                    ! Use weights for uniform sample version
                    wl = 1._c_double
                    wm = 4._c_double
                    wr = 1._c_double
                else
                    ! Use weights for non-uniform sample version
                    wl = hr*( 2*hr + 3*hl ) / (hl + hr)
                    wm = hr*(   hr + 3*hl ) / hl
                    wr = hr**3 / ( hl*(hr + hl) )
                end if
                res = res + (hl + hr) * ( wl * y(npts-2) + wm * y(npts-1) + wr * y(npts) )
            end if
            res = res / 6._c_double

            stat = 0
            if ( rule /= 0 ) stat = 0
        end if
        
    end subroutine integrate_functable

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

    subroutine argsort(arr, idx, n)
        !! Sort a float64 array using heapsort algorithm. Instead of sorting 
        !! the array, return the array of indices so that the array is sorted.
        
        integer(c_int64_t), intent(in) :: n
        !! Size of the array

        real(c_double), intent(in) :: arr(n)
        !! Array to sort

        integer(c_int64_t), intent(out) :: idx(n)
        !! Sorted index array

        integer(c_int64_t), parameter :: startpos = 1
        integer(c_int64_t) :: i, endpos, temp

        idx(1:n) = [( i, i = 1, n )] ! Initialize the index array

        ! Build heap (max-heap)
        do i = n/2, 1, -1
            call sift_down(arr, idx, i, n)
        end do

        ! Sort phase
        do endpos = n, 2, -1

            ! Swap root and last element
            temp        = idx(1)
            idx(1)      = idx(endpos)
            idx(endpos) = temp

            ! Restore heap property for reduced array
            call sift_down(arr, idx, startpos, endpos-1)

        end do
    end subroutine argsort

    subroutine sift_down(arr, idx, root, endpos)
        real(c_double), intent(in)  :: arr(:)
        integer(c_int64_t), intent(inout) :: idx(:)
        integer(c_int64_t), intent(in) :: root, endpos
        
        integer(c_int64_t) :: child, swappos
        integer(c_int64_t) :: temp

        swappos = root
        do
            child = 2*swappos
            if (child > endpos) exit  ! no children

            ! Pick larger child
            if (child < endpos) then
                if (arr( idx(child) ) < arr( idx(child+1) )) child = child + 1
            end if

            ! Compare parent with largest child
            if (arr( idx(swappos) ) < arr( idx(child) )) then
                temp         = idx(swappos)
                idx(swappos) = idx(child)
                idx(child)   = temp
                swappos      = child
            else
                exit
            end if
        end do

    end subroutine sift_down

end module quadutils

! program main
!     use iso_c_binding
!     use quadutils
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


