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

!> \brief   Symbolic module
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module handles the storage of symbolic elements produced from slater rules
!> when integrating configuration state functions
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial documentation Version.
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!> \note 22/02/2019 - Jakub Benda: Removed dependency on C++ std::map.
!>
module Symbolic_module

    use const_gbl,              only: stdout
    use consts_mpi_ci,          only: DEFAULT_INTEGRAL_THRESHOLD, NIDX
    use containers,             only: bstree
    use iso_c_binding,          only: c_loc, c_ptr, c_f_pointer
    use mpi_gbl,                only: mpi_reduceall_max, mpi_mod_rotate_arrays_around_ring, mpi_xermsg
    use precisn,                only: longint, wp
    use MemoryManager_module,   only: master_memory
    use Parallelization_module, only: grid => process_grid
    use Utility_module,         only: lexicographical_compare, pack_ints, unpack_ints

    implicit none

    public SymbolicElementVector

    private

    !> \brief This class handles the storage symbolic elements
    !>
    !> This class handles the storage symbolic elements and also expands the vector size if we have reached max capacity
    !> Additionaly, uses a binary search tree to perform a binary search on integrals labels during insertion to ensure no
    !> repeating elements this is automatic and removes the compression stage in the original SCATCI code.
    !>
    type, extends(bstree) :: SymbolicElementVector
        private
        integer(longint), allocatable :: electron_integral(:,:)  !< The packed integral storage
        real(wp),         allocatable :: coeffs(:)               !< The coefficients of the integral
        integer                       :: n = 0                   !< Number of elements in the array of both the integral and coefficients
        logical                       :: constructed = .false.   !< Whether this class has been constructed
        real(wp)                      :: threshold  = 0.0        !< The integral threshold
        integer                       :: max_capacity = 0        !< The number of free slots in the array
        integer                       :: expand_size = 10        !< How much we have to expand each
    contains
        ! type-bound procedures used by the bstree interface
        procedure, public  :: compare => bstree_compare

        ! own type-bound procedures
        procedure, public  :: construct
        procedure, public  :: insert_ijklm_symbol
        procedure, public  :: insert_symbol
        procedure, public  :: remove_symbol_at
        procedure, public  :: is_empty
        procedure, public  :: clear
        procedure, public  :: check_same_integral
        procedure, public  :: get_integral_label
        procedure, public  :: synchronize_symbols
        procedure, public  :: modify_coeff
        procedure, public  :: get_coefficient
        procedure, public  :: get_coeff_and_integral
        procedure, public  :: get_size
        procedure, public  :: add_symbols
        procedure, public  :: print => print_symbols
        procedure, private :: expand_array
        procedure, private :: check_bounds
        procedure, public  :: destroy
        procedure, private :: check_constructed
    end type  SymbolicElementVector

contains

    !> \brief   Compare two integral index sets
    !> \authors J Benda
    !> \date    2019
    !>
    !> This is the predicate (order-defining function) used by the binary search tree. Given two pointers into the
    !> electron integral array, it returns 0 when the corresponding integral index sets are equal, -1 when the
    !> first one is (lexicographically) less, and +1 when the first one is (lexicographically) greater.
    !>
    !> When either of the indices is non-positive, the dummy value pointed to by the optional argument is used instead.
    !>
    !> \param[in] this  Symbolic vector containing the electron_integral storage.
    !> \param[in] i     Position of the first integral index set.
    !> \param[in] j     Position of the first integral index set.
    !> \param[in] data  Integral index set used for comparing with sets being processed.
    !>
    !> \returns -1/0/+1 as a lexicographical spaceship operator applied on the index sets.
    !>
    integer function bstree_compare (this, i, j, data) result (verdict)

        class(SymbolicElementVector), intent(in) :: this
        integer,                      intent(in) :: i, j
        type(c_ptr),        optional, intent(in) :: data

        integer(longint)          :: ii(NIDX), jj(NIDX)
        integer(longint), pointer :: kk(:)

        ! get the default index set to compare with (if provided)
        if (present(data)) call c_f_pointer(data, kk, (/ NIDX /))

        ! obtain the index sets to compare
        if (i <= 0) then; ii = kk; else; ii = this % electron_integral(:,i); end if
        if (j <= 0) then; jj = kk; else; jj = this % electron_integral(:,j); end if

        ! compare the index sets using lexicographical comparison
        verdict = lexicographical_compare(ii, jj)

    end function bstree_compare


    !> \brief   A simple class to check if this has been properly constructed
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine check_constructed (this)
        class(SymbolicElementVector) :: this
        if (.not. this % constructed) then
            write (stdout, "('Vector::constructed - Vector is not constructed')")
            stop "Vector::constructed - Vector is not constructed"
        end if
    end subroutine check_constructed


    !> \brief   Constructs the class
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Constructs the class by allocating space for the electron integrals.
    !>
    !> \param[inout] this       Vector object to initialize.
    !> \param[in] threshold     The mininum integral coefficient threshold we should store
    !> \param[in] initial_size  Deprecated, doesnt do anything
    !>
    subroutine construct (this, threshold, initial_size)
        class(SymbolicElementVector)   :: this
        real(wp), optional, intent(in) :: threshold
        integer,  optional, intent(in) :: initial_size
        integer :: err

        if (present(threshold)) then
            this % threshold = threshold
        else
            this % threshold = DEFAULT_INTEGRAL_THRESHOLD
        end if

        this % expand_size = 100
        this % max_capacity = this % expand_size

        !Allocate the vectors
        allocate(this % electron_integral(NIDX, this % max_capacity), this % coeffs(this % max_capacity), stat = err)
        call master_memory % track_memory(storage_size(this % electron_integral)/8, size(this % electron_integral), &
                                          err, 'SYMBOLIC::ELECTRONINT')
        call master_memory % track_memory(storage_size(this % coeffs)/8, size(this % coeffs), err, 'SYMBOLIC::ELECTRONCOEFF')
        if (err /= 0) then
            write (stdout, "('SymbolicVector::construct- arrays not allocated')")
            stop "SymbolicVector arrays not allocated"
        end if

        !Do an intial clear
        call this % clear

        !This has been succesfully constructed
        this % constructed = .true.
    end subroutine construct


    !> \brief   A function to check whether the same integral label exists
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this    Vector object to query.
    !> \param[in]  integral  Packed integral label
    !> \param[out] idx       The index to the electron integral array, -1 if not found, positive otherwise
    !>
    !> \result check_same_integral True: we found one the same one. False: there isn't one
    !>
    logical function check_same_integral (this, integral, idx) result (found)

        class(SymbolicElementVector), intent(in)  :: this
        integer(longint), target,     intent(in)  :: integral(NIDX)
        integer,                      intent(out) :: idx

        call this % check_constructed

        idx   = this % locate(-1, c_loc(integral))
        found = idx > 0

    end function check_same_integral


    !> \brief   Insert unpacked integral labels
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This subroutine is used to insert unpacked integral labels. Wraps insert_symbol by automatically packing.
    !>
    !> \param[inout] this      Vector object to update.
    !> \param[in] i,j,k,l      The integral labels
    !> \param[in] m            Positron label
    !> \param[in] coeff        The integral coefficient
    !> \param[in] check_same_  Deprecated, orginally used to debug but now does nothing
    !>
    subroutine insert_ijklm_symbol (this, i, j, k, l, m, coeff, check_same_)
        class(SymbolicElementVector)  :: this
        integer,           intent(in) :: i, j, k, l, m
        real(wp),          intent(in) :: coeff
        logical, optional, intent(in) :: check_same_
        integer(longint) :: integral_label(NIDX)
        logical          :: check_same

        if (present(check_same_)) then
            check_same = check_same_
        else
            check_same = .true.
        end if

        call pack_ints(i, j, k, l, m, 0, 0, 0, integral_label)
        call this % insert_symbol(integral_label, coeff, check_same)
    end subroutine insert_ijklm_symbol


    !> \brief   Insert a packed integral symbol into the class
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This subroutine is used to insert the integral into the class, it performs a check to see if it exists, if not then
    !> it will insert (and expand if needed) into the integral array. Otherwise it will simply add the coefficients to found
    !> integral index.
    !>
    !> \param[inout] this        Vector object to update.
    !> \param[in] integral_label Packed integral label
    !> \param[in] coeff          The integral coefficient
    !> \param[in] check_same_    Deprecated, orginally used to debug but now does nothing
    !>
    subroutine insert_symbol (this, integral_label, coeff, check_same_)
        class(SymbolicElementVector)  :: this
        integer(longint),  intent(in) :: integral_label(NIDX)
        real(wp),          intent(in) :: coeff
        logical, optional, intent(in) :: check_same_
        integer :: idx = 0, check_idx
        logical :: check_same

        if (present(check_same_)) then
            check_same = check_same_
        else
            check_Same = .true.
        end if

        ! Filter small coefficients
        if (abs(coeff) < this % threshold) return

        if (.not.check_same) then

            this % n = this % n + 1                                     ! update number of elements
            idx = this % n                                              ! the idx is the last element

            if (this % n > this % max_capacity) then                    ! if we've exceeded capacity then expand
                call this % expand_array()
            end if

            this % electron_integral(:,idx) = integral_label(:)         ! place into the array
            this % coeffs(idx) = 0.0_wp                                 ! zero the coefficient since we can guareentee the compiler will do so
            call this % insert(idx)                                     ! insert idx into the binary tree

        !Check if we have the same integral, if we do get the index
        else if (.not. this % check_same_integral(integral_label, idx)) then

            this % n = this % n + 1                                     ! update number of elements
            idx = this % n                                              ! the idx is the last element

            if (this % n > this % max_capacity) then                    ! if we've exceeded capacity then expand
                call this % expand_array()
            end if

            this % electron_integral(:,idx) = integral_label(:)         ! place into the array
            this % coeffs(idx) = 0.0_wp                                 ! zero the coefficient since we can guareentee the compiler will do so
            call this % insert(idx)                                     ! insert index into the binary tree
        end if

        !Now we have an index, lets add in the coeffcient
        this % coeffs(idx) = this % coeffs(idx) + coeff
        !write(111,"('COEFF :',2es16.8)") coeff,this%coeffs(idx)

    end subroutine insert_symbol


    !> \brief   Array expansion subroutine
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> It simply copies the coeffcieint and integral array into a temporary one, reallocates a new array increased by expand_size
    !> and then recopies the elements and updates max capacity.
    !>
    subroutine expand_array(this)
        class(SymbolicElementVector)  :: this
        integer(longint), allocatable :: temp_integral(:,:)
        real(wp),         allocatable :: temp_coeffs(:)
        integer :: temp_capacity,ido,err

        !Of course lets check if we are constructed
        call this % check_constructed

        !Now the tempcapacity is the current max capacity
        temp_capacity = this % max_capacity

        allocate(temp_integral(NIDX, this % max_capacity), temp_coeffs(this % max_capacity))

        !Copy over the elements into the temporary array
        do ido = 1, temp_capacity
            temp_integral(:,ido) = this % electron_integral(:,ido)
        end do

        temp_coeffs(1:temp_capacity) = this % coeffs(1:temp_capacity)

        call master_memory % free_memory(storage_size(this % electron_integral)/8, size(this % electron_integral))
        call master_memory % free_memory(storage_size(this % coeffs)/8, size(this % coeffs))
        !deallocate our old values
        deallocate(this % electron_integral)
        deallocate(this % coeffs)

        !Increase the max capacity
        this % max_capacity = this % max_capacity + this % expand_size
        !Allocate our new array with the new max capacity
        allocate(this % electron_integral(NIDX, this % max_capacity), this % coeffs(this % max_capacity), stat = err)
        call master_memory % track_memory(storage_size(this % electron_integral)/8, size(this % electron_integral), &
                                          err, 'SYMBOLIC::EXP::ELECTRONINT')
        call master_memory % track_memory(storage_size(this % coeffs)/8, size(this % coeffs), err, 'SYMBOLIC::EXP::ELECTRONCOEFF')

        !Zero the elements
        this % electron_integral(:,:) = -1
        this % coeffs(:) = 0.0_wp

        do ido = 1, temp_capacity
            !Copy our old values back
            this % electron_integral(:,ido) = temp_integral(:,ido)
        end do
        this % coeffs(1:temp_capacity) = temp_coeffs(1:temp_capacity)
        deallocate(temp_integral, temp_coeffs)

    end subroutine expand_array


    !> \brief   This inserts one symbolic vector into another scaled by a coefficient
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this    The symbolic vector to modify.
    !> \param[in]    rhs     The symbolic vector to add into this class.
    !> \param[in]    alpha_  The scaling value to be applied to each coefficient.
    !>
    subroutine add_symbols (this, rhs, alpha_)
        class(SymbolicElementVector)             :: this
        class(SymbolicElementVector), intent(in) :: rhs
        real(wp), optional,           intent(in) :: alpha_
        real(wp) :: alpha
        integer  :: ido

        if (present(alpha_)) then
            alpha = alpha_
        else
            alpha = 1.0_wp
        end if

        !Loop through each rhs symbolic and add it to me scaled by alpha
        do ido = 1, rhs % n
            call this % insert_symbol(rhs % electron_integral(:,ido), rhs % coeffs(ido) * alpha)
        end do
    end subroutine add_symbols


    subroutine synchronize_symbols (this)
        class(SymbolicElementVector) :: this
        real(wp),         allocatable         :: coeffs(:)
        integer(longint), allocatable, target :: labels(:,:)
        integer(longint), pointer             :: label_ptr(:)
        integer(longint)                      :: my_num_symbols, largest_num_symbols, procs_num_of_symbols
        integer                               :: ido, proc_id, ierr

        if (grid % gprocs <= 1) then
            return
        end if

        my_num_symbols = this % get_size()

        call mpi_reduceall_max(my_num_symbols, largest_num_symbols, grid % gcomm)

        if (largest_num_symbols == 0) return

        call master_memory % track_memory(storage_size(labels)/8, int(largest_num_symbols * NIDX), 0, 'SYMBOLIC::MPISYNCH::LABELS')
        call master_memory % track_memory(storage_size(coeffs)/8, int(largest_num_symbols), 0, 'SYMBOLIC::MPISYNCH::COEFFS')
        allocate(labels(NIDX,largest_num_symbols), coeffs(largest_num_symbols), stat = ierr)

        label_ptr(1:largest_num_symbols*NIDX) => labels(:,:)
        labels = 0
        coeffs = 0

        if (my_num_symbols > 0) then
            labels(:,1:my_num_symbols) = this % electron_integral(:,1:my_num_symbols)
            coeffs(1:my_num_symbols) = this % coeffs(1:my_num_symbols)
        end if

        procs_num_of_symbols = my_num_symbols

        do proc_id = 1, grid % gprocs - 1
            call mpi_mod_rotate_arrays_around_ring(procs_num_of_symbols, label_ptr, coeffs, largest_num_symbols, grid % gcomm)
            do ido = 1, procs_num_of_symbols
                call this % insert_symbol(labels(:,ido), coeffs(ido))
            end do
        end do

        call master_memory % free_memory(storage_size(labels)/8, size(labels))
        call master_memory % free_memory(storage_size(coeffs)/8, size(coeffs))
        deallocate(labels, coeffs)

    end subroutine synchronize_symbols


    !> \brief   Range check
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Simply checks the index and wheter it exists within the array.
    !>
    !> \param[inout] this  Vector object to query.
    !> \param[in]    i     The index we wish to access
    !>
    logical function check_bounds (this, i)
        class(SymbolicElementVector) :: this
        integer, intent(in) :: i

        if (i <= 0 .or. i > this % n) then
            write (stdout, "('SymbolicVector::check_bounds - Out of Bounds access')")
            stop "SymbolicVector::check_bounds - Out of Bounds access"
            check_bounds = .false.
        else
            check_bounds = .true.
        end if

    end function check_bounds


    !> \brief  Removes an integral and coefficient
    !> \author A Al-Refaie
    !> \date   2017
    !>
    !> Never been used and pretty much useless.
    !>
    subroutine remove_symbol_at (this, idx)
        class(SymbolicElementVector) :: this
        integer, intent(in)          :: idx
        integer                      :: ido

        if (this % check_bounds(idx)) then
            do ido = idx + 1, this % n
                this % electron_integral(:,ido-1) = this % electron_integral(:,ido)
                this % coeffs(ido-1)              = this % coeffs(ido)
            end do
            this % electron_integral(:, this % n) = 0
            this % coeffs(this % n) = 0
            this % n = this % n - 1
        end if

    end subroutine remove_symbol_at


    !> \brief   Emptiness check
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    !>
    !> Simply returns whether we are storing integrals and coeffs or not.
    !>
    !> \result True if empty or False if we have any symbols.
    !>
    logical function is_empty (this)
        class(SymbolicElementVector) :: this

        is_empty = (this % n == 0)

    end function is_empty


    !> \brief   Clear array
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    !>
    !> Clears our array (not really but it essential resets the symbol counter to zero, which is way quicker).
    !>
    subroutine clear (this)
        class(SymbolicElementVector) :: this

        this % n = 0
        call this % bstree % destroy

    end subroutine clear


    !> \brief   Update coeff with contribution
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Add contribution to coeff at given index.
    !>
    subroutine modify_coeff (this, idx, coeff)
        class(SymbolicElementVector) :: this
        integer, intent(in) :: idx
        real(wp) :: coeff

        this % coeffs(idx) = this % coeffs(idx) + coeff

    end subroutine modify_coeff


    !> \brief   Get integral label at specific index
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    !>
    function get_integral_label (this, idx)
        class(SymbolicElementVector)   :: this
        integer, intent(in)            :: idx
        integer(longint)               :: get_integral_label(NIDX)

        if (this % check_bounds(idx)) get_integral_label(:) = this % electron_integral(:,idx)

    end function get_integral_label


    !> \brief   Get coefficient at specific index
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    !>
    real(wp) function get_coefficient (this, idx)
        class(SymbolicElementVector) :: this
        integer, intent(in)          :: idx

        if (this % check_bounds(idx)) get_coefficient = this % coeffs(idx)

    end function get_coefficient


    !> \brief   Get both label and coeffcient at specific index
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    !>
    subroutine get_coeff_and_integral (this, idx, coeff, label)
        class(SymbolicElementVector), intent(in)  :: this
        integer,                      intent(in)  :: idx
        real(wp),                     intent(out) :: coeff
        integer(longint),             intent(out) :: label(NIDX)

        if (this % check_bounds(idx)) then
            label = this % electron_integral(:,idx)
            coeff = this % coeffs(idx)
        end if

    end subroutine get_coeff_and_integral


    !> \brief   Returns the number of symbolic elements stored
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    integer function get_size(this)
        class(SymbolicElementVector) :: this
        get_size = this % n
    end function get_size


    !> \brief   Cleans up the class by deallocating arrays
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    !>
    subroutine destroy (this)
        class(SymbolicElementVector) :: this

        if (allocated(this % electron_integral)) then
            call master_memory % free_memory(storage_size(this % electron_integral)/8, size(this % electron_integral))
            deallocate(this % electron_integral)
        end if
        if (allocated(this % coeffs)) then
            call master_memory % free_memory(storage_size(this % coeffs)/8, size(this % coeffs))
            deallocate(this % coeffs)
        end if

        call this % bstree % destroy

        this % constructed = .false.

    end subroutine destroy


    !> \brief   Print currently stored symbols
    !> \authors A Al-Refaie
    !> \date    2017
    !> \public
    subroutine print_symbols (this)
        class(SymbolicElementVector) :: this
        integer  :: labels(8), ido
        real(wp) :: coeff

        if (.not. this % is_empty()) write (stdout, "('Outputting symbolic elements....')")
        do ido = 1, this % n
            call unpack_ints(this % electron_integral(:,ido), labels)
            write (stdout, "(5i4,' -- ',es14.3)") labels(1:5),this % coeffs(ido)
        end do

    end subroutine print_symbols

end module Symbolic_module
