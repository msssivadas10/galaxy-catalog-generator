module functable
    !! WARNING: will be removed later
    use iso_c_binding
    implicit none

    private
    public :: make_functable

    integer, parameter :: dp = c_double

    type, bind(c) :: heapinfo_t
        integer(c_int64_t) :: size     !! Current size of the heap array
        integer(c_int64_t) :: itemsize !! Size of an item in the heap
        integer(c_int64_t) :: capacity !! Maximum size of the heap
    end type
    
contains

    subroutine make_functable(f, a, b, args, npts, unif, xy, stat) bind(c)
        !! Make a table of the values of function f(x) in the interval [a, b].

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

        integer(c_int), intent(in), value :: unif
        !! Specify uniform (1) or non-uniform (0) sampling.

        integer(c_int64_t), intent(in), value :: npts
        !! Number of points. If `rule=0`, this must be a multiple of 7.

        real(c_double), intent(out) :: xy(2, npts)
        !! Integration nodes and function values

        integer(c_int), intent(out) :: stat
        !! Error code: 0=ok, 1=error

        integer(c_int), parameter :: datsize = 7
        
        integer(c_int64_t) :: i, idx(npts)
        real(c_double)     :: xa, xb, xm, xma, xmb, ya, yb, ym, yma, ymb, err, scale
        real(c_double)     :: xs(npts), ys(npts)
        real(c_double)     :: idata(7) ! Interval data: [ -error, left, value, right, value, mid, value ]    
        type(heapinfo_t)   :: hinfo
        real(c_double), allocatable :: heap(:, :) ! Priority queue for sub-intervals 
        ! The heap key is the error in the linear approximation of the function in 
        ! that interval. Intervals are divided if error is large.

        stat   = 1

        xa    = a
        ya    = f(xa, args)
        xs(1) = xa 
        ys(1) = ya 
        xb    = b
        yb    = f(xb, args)
        xs(2) = xb
        ys(2) = yb
        if ( unif /= 1 ) then
            ! Sampling npts points uniformly from the interval [a, b]
            
            scale = (b - a) / (npts - 1) ! Stepsize
            
            xs(npts) = xs(2)
            ys(npts) = ys(2)
            do i = 2, npts-1
                xm    = xs(i-1) + scale
                xs(i) = xm
                ys(i) = f(xm, args) 
            end do
        else
            ! Adaptive sampling of the interval into npts points, based on 
            ! how the function changes
            
            hinfo%size     = 0
            hinfo%itemsize = 7
            hinfo%capacity = 10*npts ! Heap capacity
            allocate( heap(hinfo%itemsize, hinfo%capacity) )

            xm       = 0.5_dp * ( xa + xb )
            ym       = f(xm, args)
            err      = abs( ym - 0.5_dp * ( ya + yb ) )
            idata(:) = [ -err, xa, ya, xb, yb, xm, ym ]
            call heap_push(heap, hinfo, idata)
            
            i = 2
            do while ( i < npts .and. hinfo%size > 0 )

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
            call argsort(xs, idx, npts)
            xy(1, 1:npts) = xs(1:npts) ! Using xy array as temp array for rearranging...
            xy(2, 1:npts) = ys(1:npts)
            do i = 1, npts
                xs(i) = xy(1, idx(i))
                ys(i) = xy(2, idx(i))
            end do
            
        end if

        xy(1,1:npts) = xs(1:npts)
        xy(2,1:npts) = ys(1:npts)
        
    end subroutine make_functable

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
    
end module functable