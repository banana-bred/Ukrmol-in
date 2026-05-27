! Copyright 2019
!
! For a comprehensive list of the developers that contributed to these codes
! see the UK-AMOR website.
!
! This file is part of UKRmol-in (UKRmol+ suite).
!
!     UKRmol-in is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     UKRmol-in is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  UKRmol-in (in source/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.

!> \brief   Matrix cache module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module MatrixCache_module

    use precisn,              only: wp
    use const_gbl,            only: stdout
    use Timing_Module,        only: master_timer
    use MemoryManager_module, only: master_memory
    use consts_mpi_ci,        only: DEFAULT_EXPAND_SIZE

    implicit none

    public MatrixCache

    private

    type :: MatrixArray
        integer,  pointer :: ij(:,:)
        real(wp), pointer :: coefficient(:)
        integer           :: num_elems
        logical           :: constructed = .false.
    contains
       !procedure, public :: is_valid
        procedure, public :: destroy   => destroy_matrix_array
        procedure, public :: get       => get_from_array
        procedure, public :: insert    => insert_into_array
        procedure, public :: construct => construct_array
        procedure, public :: clear     => clear_array
    end type MatrixArray

    !> \brief   This handles the matrix elements and also expands the vector size if we have reached max capacity
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    type :: MatrixCache
        type(MatrixArray), allocatable :: matrix_arrays(:)
        integer :: n                                    !< Number of elements in the array of both the integral and coefficients
        logical :: constructed = .false.
        integer :: matrix_index
        integer :: num_matrix_chunks = 1
        integer :: max_capacity = 0                     !< The number of free slots in the array
        integer :: expand_size = DEFAULT_EXPAND_SIZE    !< How much we have to expand each
    contains
        procedure, public  :: construct
        procedure, public  :: insert_into_cache
        procedure, public  :: get_from_cache
        procedure, public  :: is_empty
        procedure, public  :: clear
        procedure, public  :: get_size
        procedure, public  :: expand_capacity
       !procedure, public  :: prune_threshold
        procedure, public  :: sort  => qsort_matelem
        procedure, public  :: print => print_matelems
        procedure, public  :: get_local_mat
        procedure, public  :: get_chunk_idx
        procedure, private :: quick_sort_function
        procedure, private :: partition_function
        procedure, private :: expand_array
        procedure, public  :: shrink_capacity
        procedure, public  :: clear_and_shrink
        procedure, private :: check_bounds
        procedure, public  :: destroy => destroy_cache
        procedure, private :: check_constructed
       !procedure, public  :: count_occurance
       !procedure, public  :: construct
    end type MatrixCache

contains

    subroutine check_constructed (this)
        class(MatrixCache) :: this

        if (.not. this % constructed) then
            write (stdout, "('Vector::constructed - Vector is not constructed')")
            stop "Vector::constructed - Vector is not constructed"
        end if

    end subroutine check_constructed


    integer function get_chunk_idx (this, idx)
        class(MatrixCache)  :: this
        integer, intent(in) :: idx
        integer             :: mat_idx, local_idx

        call this % get_local_mat(idx, mat_idx, local_idx)

        get_chunk_idx = mat_idx

    end function get_chunk_idx


    subroutine construct (this, expand_size_)
        class(MatrixCache)            :: this
        integer                       :: err
        integer, optional, intent(in) :: expand_size_

        if (present(expand_size_)) then
            this % expand_size = expand_size_
        else
            this % expand_size = DEFAULT_EXPAND_SIZE
        end if

        this % max_capacity = this % expand_size
        this % n = 0
        this % matrix_index = 1

        !Allocate the vectors
        allocate(this%matrix_arrays(this%matrix_index),stat=err)
        call master_memory%track_memory(storage_size(this%matrix_arrays)/8, size(this%matrix_arrays), err, &
                                        'MATRIXCACHE::MATRIXARRAY')

        call this % matrix_arrays(this % matrix_index) % construct(this % expand_size)
        !this%matrix_arrays(this%matrix_index)%num_elems = 0

        this % num_matrix_chunks = 1

        call this % clear

        if (err /= 0) then
            write (stdout, "('Matrix Cache::construct- arrays not allocated')")
            stop "Matrix Cache:: arrays not allocated"
        end if

        call this % clear
        this % constructed = .true.

    end subroutine construct


    subroutine insert_into_cache (this, i, j, coefficient)
        class(MatrixCache)   :: this
        integer,  intent(in) :: i, j
        real(wp), intent(in) :: coefficient
        integer :: mat_id, local_index

       !this % matrix_size = max(this % matrix_size, i, j)
        this % n = this % n + 1

        if (this % n > this % max_capacity) then
            call this % expand_array()
        end if

        call this % get_local_mat(this % n, mat_id, local_index)

        this % num_matrix_chunks = mat_id
        this % matrix_arrays(mat_id) % ij(1, local_index) = i
        this % matrix_arrays(mat_id) % ij(2, local_index) = j
        this % matrix_arrays(mat_id) % coefficient(local_index) = coefficient

        this % matrix_arrays(mat_id) % num_elems = this % matrix_arrays(mat_id) % num_elems + 1

    end subroutine insert_into_cache


    subroutine get_local_mat (this, idx, mat_id, local_index)
        class(MatrixCache)   :: this
        integer, intent(in)  :: idx
        integer, intent(out) :: mat_id
        integer, intent(out) :: local_index

        if (this % check_bounds(idx)) then
            mat_id = (idx - 1) / this % expand_size + 1
            !if(mat_id > this%num_matrix_chunks) then
                !write(stdout,"('Error, referencing a chunk larger than me! [idx,max_Capacity,expand_size,num_chunks,referenced_chunk] = ',5i8)") idx,this%max_capacity,this%expand_size,this%num_matrix_chunks,mat_id
                !stop "matid error"
            !endif
            local_index = idx - (mat_id - 1) * this % expand_size
        end if

    end subroutine get_local_mat


    subroutine get_from_cache (this, idx, i, j, coeff)
        class(MatrixCache)    :: this
        integer,  intent(in)  :: idx
        integer,  intent(out) :: i, j
        real(wp), intent(out) :: coeff

        integer :: mat_id, local_index

        if (this % check_bounds(idx)) then
            call this % get_local_mat(idx, mat_id, local_index)

            i = this % matrix_arrays(mat_id) % ij(1, local_index)
            j = this % matrix_arrays(mat_id) % ij(2, local_index)

            coeff = this % matrix_arrays(mat_id) % coefficient(local_index)
        end if

    end subroutine get_from_cache


    subroutine expand_array (this)
        class(MatrixCache) :: this
        type(MatrixArray)  :: temp_matrix(this % num_matrix_chunks + 1)
        integer            :: new_num_mats, err

        call this % check_constructed

       !call master_timer%start_timer("Expand array")
        temp_matrix(1:this % num_matrix_chunks) = this % matrix_arrays(1:this % num_matrix_chunks)
        call master_memory % free_memory(storage_size(this % matrix_arrays)/8, size(this % matrix_arrays))
        deallocate(this % matrix_arrays)

        this % max_capacity = this % max_capacity + this % expand_size

        new_num_mats = this % max_capacity / this % expand_size
        this % matrix_index = new_num_mats
        this % num_matrix_chunks = this % num_matrix_chunks + 1
        allocate(this % matrix_arrays(this % num_matrix_chunks), stat = err)
        call master_memory % track_memory(storage_size(this % matrix_arrays)/8, size(this % matrix_arrays), err, &
                                          'MATRIXCACHE::MATRIXARRAY')

        this % matrix_arrays(:) = temp_matrix(:)
        this % matrix_arrays(this % num_matrix_chunks) % num_elems = 0
        call this % matrix_arrays(this % num_matrix_chunks) % construct(this % expand_size)
       !call master_timer%stop_timer("Expand array")

    end subroutine expand_array


    subroutine destroy_matrix_array (this)
        class(MatrixArray) :: this

        if (.not. this % constructed) return
        if (associated(this % ij)) then
            call master_memory % free_memory(storage_size(this % ij)/8, size(this % ij))
            deallocate(this % ij)
            this % ij => null()
        end if
        if (associated(this % coefficient)) then
            call master_memory % free_memory(storage_size(this % coefficient)/8, size(this % coefficient))
            deallocate(this % coefficient)
            this % coefficient => null()
        end if
        this % constructed = .false.
        this % num_elems = 0

    end subroutine destroy_matrix_array


    subroutine construct_array (this, capacity)
        class(MatrixArray)  :: this
        integer, intent(in) :: capacity
        integer             :: err

        if (.not. this % constructed) then
            call master_memory % track_memory(storage_size(this % ij)/8, capacity * 2, 0, 'MATRIXCACHE::MATRIXARRAY::IJ')
            call master_memory % track_memory(storage_size(this % coefficient)/8, capacity, 0, 'MATRIXCACHE::MATRIXARRAY::COEFF')
            allocate(this % ij(2, capacity), this % coefficient(capacity), stat = err)
            this % constructed = .true.
        end if

    end subroutine construct_array


    subroutine clear_array (this)
        class(MatrixArray) :: this

        this % num_elems = 0

    end subroutine clear_array


    subroutine insert_into_array (this, i, j, coeff)
        class(MatrixArray)   :: this
        integer,  intent(in) :: i, j
        real(wp), intent(in) :: coeff

        this % num_elems = this % num_elems + 1
        this % ij(this % num_elems, 1) = i
        this % ij(this % num_elems, 2) = j
        this % coefficient(this % num_elems) = coeff

    end subroutine insert_into_array


    subroutine get_from_array(this, idx, i, j, coeff)
        class(MatrixArray), intent(in)  :: this
        integer,            intent(in)  :: idx
        integer,            intent(out) :: i, j
        real(wp),           intent(out) :: coeff

        if (idx > this % num_elems) stop "Matrix Array access segfault!"

        i = this % ij(1, idx)
        j = this % ij(2, idx)
        coeff = this % coefficient(idx)

    end subroutine get_from_array


    subroutine expand_capacity (this, capacity)
        class(MatrixCache)  :: this
        integer, intent(in) :: capacity

        do while (this % max_capacity < capacity)
            call this % expand_array
        end do

    end subroutine expand_capacity


    logical function check_bounds (this, i)
        class(MatrixCache), intent(in) :: this
        integer,            intent(in) :: i

        if (i <= 0 .or. i > this % n) then
            write (stdout, "('MatrixCache::check_bounds - Out of Bounds access', 2i4)") i, this % n
            stop "MatrixCache::check_bounds - Out of Bounds access"
            check_bounds = .false.
        else
            check_bounds = .true.

        end if

    end function check_bounds


    logical function is_empty (this)
        class(MatrixCache), intent(in) :: this

        is_empty = (this % n == 0)

    end function is_empty


    subroutine clear (this)
        class(MatrixCache) :: this
        integer            :: ido

        if (allocated(this % matrix_arrays)) this % matrix_arrays(:) % num_elems = 0
        this % n = 0

    end subroutine clear


    integer function get_size (this)
        class(MatrixCache), intent(in) :: this

        get_size = this % n

    end function get_size


    subroutine shrink_capacity (this)
        class(MatrixCache) :: this
        integer            :: ido, total_size, shrink_size

        type(MatrixArray), allocatable :: temp_matrix(:)

        if (allocated(this % matrix_arrays)) then
            total_size = size(this % matrix_arrays)
            shrink_size = min(this % n / this % expand_size + 1, this % num_matrix_chunks)

            if (shrink_size == this % num_matrix_chunks) return

            do ido = shrink_size + 1, total_size
                call this % matrix_arrays(ido) % destroy
            end do

            allocate(temp_matrix(shrink_size))
            call master_memory % track_memory(storage_size(temp_matrix)/8, size(temp_matrix), 0, 'MATRIXCACHE::SHRINK::TEMP')
            temp_matrix(1:shrink_size) = this % matrix_arrays(1:shrink_size)
            call master_memory % free_memory(storage_size(this % matrix_arrays)/8, size(this % matrix_arrays))
            deallocate(this % matrix_arrays)

            allocate(this % matrix_arrays(shrink_size))
            call master_memory % track_memory(storage_size(this % matrix_arrays)/8, size(this % matrix_arrays), 0, &
                                              'MATRIXCACHE::SHRINK::MATRIXARRAYS')
            this % matrix_arrays(1:shrink_size) = temp_matrix(1:shrink_size)
            this % max_capacity = shrink_size * this % expand_size
            this % num_matrix_chunks = shrink_size
            call master_memory % free_memory(storage_size(temp_matrix)/8, size(temp_matrix))
            deallocate(temp_matrix)
        end if

    end subroutine shrink_capacity


    subroutine clear_and_shrink (this)
        class(MatrixCache) :: this

        call this % clear
        call this % shrink_capacity

    end subroutine clear_and_shrink


    subroutine destroy_cache (this)
        class(MatrixCache) :: this
        integer            :: num_arrays, ido

        if (allocated(this % matrix_arrays)) then
            num_arrays = size(this % matrix_arrays)

            do ido = 1, num_arrays
                call this % matrix_arrays(ido) % destroy
            end do

            call master_memory % free_memory(storage_size(this % matrix_arrays)/8, size(this % matrix_arrays))
            deallocate(this % matrix_arrays)
        end if

        this % constructed = .false.

    end subroutine destroy_cache


    subroutine print_matelems (this)
        class(MatrixCache) :: this
        integer            :: labels(8), ido, arrs, jdo, elm
        real(wp)           :: coeff

        write (stdout, "('Outputting Matrix elements....')")

        this % matrix_index = this % n / this % expand_size + 1
        elm = 0

        do ido = 1, this % num_matrix_chunks
            write (stdout, *) this % matrix_arrays(ido) % num_elems
            do jdo = 1, this % matrix_arrays(ido) % num_elems
                write (stdout, "(2i8,' -- ',D16.8)") this % matrix_arrays(ido) % ij(1:2,jdo), &
                                                     this % matrix_arrays(ido) % coefficient(jdo)
            end do
        end do

    end subroutine print_matelems


    subroutine qsort_matelem (this)
        class(MatrixCache) :: this

        call this % check_constructed
        if (this % n <= 1) return
        call this % quick_sort_function(1, this % n)
       !call QsortMatrixElement_ji(this % matrix_elements(1:this % n))

    end subroutine qsort_matelem


    recursive subroutine quick_sort_function (this, start_idx, end_idx)
        class(MatrixCache)  :: this
        integer, intent(in) :: start_idx, end_idx
        integer :: mat_size, iq

        mat_size = end_idx - start_idx + 1

        if (mat_size > 1) then
            call this % partition_function(iq, start_idx, end_idx)
            call this % quick_sort_function(start_idx, iq - 1)
            call this % quick_sort_function(iq, end_idx)
        end if

    end subroutine quick_sort_function


    subroutine partition_function (this, marker, start_idx, end_idx)
        class(MatrixCache)   :: this
        integer, intent(in)  :: start_idx, end_idx
        integer, intent(out) :: marker
        integer  :: i, j, temp_ij(2), x_ij(2), A_ij(2), mat_id_i, mat_id_j, local_index_i, local_index_j
        real(wp) :: temp_coeff, x_coeff, A_coeff

        call this % get_from_cache(start_idx, x_ij(1), x_ij(2), x_coeff)

        i = start_idx - 1
        j = end_idx + 1

        do
            j = j - 1
            do
                call this % get_from_cache(j, A_ij(1), A_ij(2), A_coeff)
                if (A_ij(1) <= x_ij(1)) exit
                j = j - 1
            end do
            i = i + 1
            do
                call this % get_from_cache(i, A_ij(1), A_ij(2), A_coeff)
                if (A_ij(1) >= x_ij(1)) exit
                i = i + 1
            end do
            if (i < j) then
                ! exchange A(i) and A(j)
                call this % get_local_mat(i, mat_id_i, local_index_i)
                call this % get_local_mat(j, mat_id_j, local_index_j)

                temp_ij(:) = this % matrix_arrays(mat_id_i) % ij(:,local_index_i)
                temp_coeff = this % matrix_arrays(mat_id_i) % coefficient(local_index_i)

                this % matrix_arrays(mat_id_i) % ij(:,local_index_i) &
                    = this % matrix_arrays(mat_id_j) % ij(:,local_index_j)
                this % matrix_arrays(mat_id_i) % coefficient(local_index_i) &
                    = this % matrix_arrays(mat_id_j) % coefficient(local_index_j)

                this % matrix_arrays(mat_id_j) % ij(:,local_index_j) = temp_ij(:)
                this % matrix_arrays(mat_id_j) % coefficient(local_index_j) = temp_coeff
            else if (i == j) then
                marker = i + 1
                return
            else
                marker = i
                return
            end if
        end do

    end subroutine partition_function

    ! Recursive Fortran 95 quicksort routine adapted to
    ! sort matrix elements into ascending either ij or ji
    ! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
    ! Based on algorithm from Cormen et al., Introduction to Algorithms,
    ! 1997 printing

    ! Made F conformant by Walt Brainerd

    !subroutine construct(this,i,j,k,l,inttype,positron,coefficent)
    !    class(MatrixElement),intent(inout)    ::    this
    !    integer,intent(in)            :: i,j,k,l,positron,inttype
    !    real(wp),intent(in)            :: coefficent
    !
    !    this%coefficient = coefficient
    !    call pack4ints64(i,j,k,l,this%integral_label)
    !    this%positron = positron
    !
    !    if(inttype /= ONE_ELECTRON_INTEGRAL .or. inttype /= TWO_ELECTRON_INTEGRAL) then
    !        write(stdout,"('Error integral type of ',i4,' is not 1 or 2 electron')") inttype
    !        stop "Integral label is neither 1 or 2 electron"
    !    else
    !        this%integral_type = inttype
    !    endif
    !
    !
    !
    !
    !end subroutine

end module MatrixCache_module
