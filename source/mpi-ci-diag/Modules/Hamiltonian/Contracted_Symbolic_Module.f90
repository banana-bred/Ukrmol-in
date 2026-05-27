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
!> This module handles the storage of symbolic elements produced from Slater rules
!> when integrating configuration state functions.
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial documentation version
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!> \note 22/02/2019 - Jakub Benda: Removed dependency on C++ std::map.
!>
module Contracted_Symbolic_module

    use const_gbl,              only: stdout
    use consts_mpi_ci,          only: DEFAULT_INTEGRAL_THRESHOLD, NIDX
    use containers,             only: bstree
    use iso_c_binding,          only: c_loc, c_ptr, c_f_pointer
    use mpi_gbl,                only: myrank, mpiint, mpi_mod_allgather, mpi_mod_rotate_cfp_arrays_around_ring, &
                                      mpi_mod_rotate_int_arrays_around_ring, mpi_xermsg, mpi_reduceall_max
    use precisn,                only: longint, wp
    use MemoryManager_module,   only: master_memory
    use Parallelization_module, only: grid => process_grid
    use Symbolic_module,        only: SymbolicElementVector
    use Timing_Module,          only: master_timer
    use Utility_module,         only: lexicographical_compare, pack_ints, unpack_ints

    implicit none

    public ContractedSymbolicElementVector

    private

    !> \brief   This class handles the storage symbolic elements
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This class handles the storage symbolic elements and also expands the vector size if we have reached max capacity
    !> Additionaly, uses a binary search tree to perform a binary search on integrals labels during insertion to ensure no
    !> repeating elements this is automatic and removes the compression stage in the original SCATCI code.
    !>
    type, extends(bstree) :: ContractedSymbolicElementVector
        private

        integer(longint), allocatable :: electron_integral(:,:)   !< The packed integral storage
        real(wp),         allocatable :: coeffs(:,:,:)            !< The coefficients of the integral

        integer  :: n = 0                   !< Number of elements in the array of both the integral and coefficients
        logical  :: constructed = .false.   !< Whether this class has been constructed
        real(wp) :: threshold  = 0.0        !< The integral threshold
        integer  :: num_states_1 = 0
        integer  :: num_states_2 = 0
        integer  :: max_capacity = 0        !< The number of free slots in the array
        integer  :: expand_size = 100       !< How much we have to expand each
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
        procedure, public  :: estimate_synchronize_cost
        procedure, public  :: modify_coeff
        procedure, public  :: get_coefficient
        procedure, public  :: get_coeff_and_integral
        procedure, public  :: get_size
        procedure, public  :: get_num_targets_sym1
        procedure, public  :: get_num_targets_sym2
        procedure, public  :: add_symbols
        procedure, public  :: reduce_symbols
        procedure, public  :: print => print_symbols
       !procedure, public  :: prune_threshold
       !procedure, public  :: sort_symbols
        procedure, private :: expand_array
        procedure, private :: check_bounds
        procedure, public  :: destroy
        procedure, private :: check_constructed
        procedure, private :: synchronize_symbols_II
       !procedure, public  :: count_occurance
       !procedure, public  :: construct
    end type  ContractedSymbolicElementVector

    class(ContractedSymbolicElementVector), pointer :: to_be_synched

contains

    !> \brief   Compare two integral index sets
    !> \authors J Benda
    !> \date    2019
    !>
    !> This is the predicate (order-defining function) used by the binary search tree. Given two pointers into the
    !> electron integral array, it returns 0 when the corresponding integral index sets are equal, -1 when the
    !> first one is (lexicographically) less, and +1 when the first one is (lexicographically) greater.
    !>
    !> When either of the indices is non-positive, the dummy value stored in the dummy argument `data` is used instead of
    !> the corresponding index set.
    !>
    !> \param[in] this  Symbolic vector containing the reference electron_integral storage.
    !> \param[in] i     Position of the first integral index set.
    !> \param[in] j     Position of the first integral index set.
    !> \param[in] data  Integral index set to which compare the set being processed.
    !>
    !> \returns -1/0/+1 as a lexicographical spaceship operator applied on the index sets.
    !>
    integer function bstree_compare (this, i, j, data) result (verdict)

        class(ContractedSymbolicElementVector), intent(in) :: this
        integer,                      intent(in) :: i, j
        type(c_ptr),        optional, intent(in) :: data

        integer(longint)          :: ii(NIDX), jj(NIDX)
        integer(longint), pointer :: kk(:)

        ! get the default index set to compare with (if provided)
        if (present(data)) call c_f_pointer(data, kk, (/ NIDX /))

        ! obtain pointers to the index sets
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
        class(ContractedSymbolicElementVector) :: this
        if (.not. this % constructed) then
            write (stdout, "('Vector::constructed - Vector is not constructed')")
            stop "Vector::constructed - Vector is not constructed"
        end if

    end subroutine check_constructed


    !> \brief   Constructs the class by allocating space for the electron integrals.
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout]  this      Vector object to update.
    !> \param[in] n1            Number of states 1.
    !> \param[in] n2            Number of states 2.
    !> \param[in] threshold     The mininum integral coefficient threshold we should store
    !> \param[in] initial_size  Deprecated, doesnt do anything
    !>
    subroutine construct (this, n1, n2, threshold, initial_size)
        class(ContractedSymbolicElementVector) :: this
        real(wp), optional, intent(in)         :: threshold
        integer,  optional, intent(in)         :: initial_size

        integer :: err, n1, n2

        if (present(threshold)) then
            this % threshold = threshold
        else
            this % threshold = DEFAULT_INTEGRAL_THRESHOLD
        end if

        this % expand_size = 100
        this % max_capacity = this % expand_size
        this % num_states_1 = n1
        this % num_states_2 = n2

        !Allocate the vectors
        allocate(this % electron_integral(NIDX, this % max_capacity), this % coeffs(n1, n2, this % max_capacity), stat = err)
        call master_memory % track_memory (storage_size(this % electron_integral)/8, &
                                           size(this % electron_integral), err, 'CONSYMBOLIC::ELECTRONINT')
        call master_memory % track_memory(storage_size(this % coeffs)/8, &
                                          size(this % coeffs), err, 'CONSYMBOLIC::ELECTRONCOEFF')
        if (err /= 0) then
            call mpi_xermsg('Contracted_Symbolic_module', 'construct', 'SymbolicVector arrays not allocated', 1, 1)
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

        class(ContractedSymbolicElementVector), intent(in)  :: this
        integer(longint),  target,              intent(in)  :: integral(NIDX)
        integer,                                intent(out) :: idx

        call this % check_constructed

        idx   = this % locate(-1, c_loc(integral))
        found = idx > 0

    end function check_same_integral


    !> \brief Insert unpacked integral labels.
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Wraps insert_symbol by automatically packing.
    !>
    !> \param[inout]  this     Vector object to update.
    !> \param[in] i,j,k,l      The integral labels
    !> \param[in] m            Positron label
    !> \param[in] coeffs       The integral coefficient
    !> \param[in] check_same_  Deprecated, orginally used to debug but now does nothing
    !>
    subroutine insert_ijklm_symbol (this, i, j, k, l, m, coeffs, check_same_)
        class(ContractedSymbolicElementVector) :: this
        integer,  intent(in)            :: i, j, k, l, m
        real(wp), intent(in)            :: coeffs(this % num_states_1, this % num_states_2)
        logical,  intent(in), optional  :: check_same_

        integer(longint) :: integral_label(NIDX)
        logical          :: check_same

        if (present(check_same_)) then
            check_same = check_same_
        else
            check_same = .true.
        end if

        call pack_ints(i, j, k, l, m, 0, 0, 0, integral_label)
        call this%insert_symbol(integral_label, coeffs, check_same)

    end subroutine insert_ijklm_symbol


    !> \brief   Insert a packed integral symbol into the class
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This subroutine is used to insert the integral into the class, it performs a check to see if it exists, if not then
    !> it will insert (and expand if needed) into the integral array. Otherwise it will simply add the coefficients
    !> to found integral index.
    !>
    !> \param[inout]  this        Vector object to update.
    !> \param[in] integral_label  Packed integral label
    !> \param[in] coeffs          The integral coefficient
    !> \param[in] check_same_     Deprecated, orginally used to debug but now does nothing
    !>
    subroutine insert_symbol (this, integral_label, coeffs, check_same_)
        class(ContractedSymbolicElementVector) :: this
        integer(longint), intent(in)           :: integral_label(NIDX)
        real(wp),         intent(in)           :: coeffs(this % num_states_1, this % num_states_2)
        logical,          intent(in), optional :: check_same_

        integer :: idx = 0, check_idx
        logical :: check_same

        if (present(check_same_)) then
            check_same = check_same_
        else
            check_Same = .true.
        end if

        ! Filter small coefficients
        !if(abs(coeff) < this%threshold) return

        if (.not. this % check_same_integral(integral_label, idx)) then
                !If not then update number of elements
                this % n = this % n + 1
                !The idx is the last element
                idx = this % n
                !If we've exceeded capacity then expand
                if (this % n > this % max_capacity) then
                    call this % expand_array()
                end if
                !Place into the array
                this % electron_integral(:,idx) = integral_label(:)
                !Zero the coefficient since we can guareentee the compiler will do so
                this % coeffs(:,:,idx) = 0.0_wp
                !Insert it into the binary tree
                call this % insert(idx)
        end if

        !Now we have an index, lets add in the coeffcient
        this % coeffs(:,:,idx) = this % coeffs(:,:,idx) + coeffs(:,:)

    end subroutine insert_symbol


    !> \brief   This is the array expansion subroutine.
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> It simply copies the coeffcieint and integral array into a temporary one, reallocates a new array increased
    !> by expand_size and then recopies the elements and updates max capacity.
    !>
    subroutine expand_array (this)
        class(ContractedSymbolicElementVector) :: this
        integer(longint), allocatable :: temp_integral(:,:)
        real(wp),         allocatable :: temp_coeffs(:,:,:)
        integer :: temp_capacity, ido, err

        !write(stdout,"('Size of arrays = ',5i12)") this%max_capacity,this%num_states_1,this%num_states_2,size(this%electron_integral),size(this%coeffs)

        !Of course lets check if we are constructed
        call this % check_constructed

        !Now the tempcapacity is the current max capacity
        temp_capacity = this % max_capacity

        allocate(temp_integral(NIDX, this % max_capacity), &
                 temp_coeffs(this % num_states_1, this % num_states_2, this % max_capacity), stat = err)
        call master_memory % track_memory(storage_size(temp_integral)/8, size(temp_integral), err, 'CONSYMBOLIC::EXP::TEMPINT')
        call master_memory % track_memory(storage_size(temp_coeffs)/8, size(temp_coeffs), err, 'CONSYMBOLIC::EXP::TEMPCOEFF')

        !Copy over the elements into the temporary array
        do ido = 1, temp_capacity
            temp_integral(:,ido) = this % electron_integral(:,ido)
            temp_coeffs(:,:,ido) = this % coeffs(:,:,ido)
        end do

        call master_memory % free_memory(storage_size(this % electron_integral)/8, size(this % electron_integral))
        call master_memory % free_memory(storage_size(this % coeffs)/8, size(this % coeffs))

        !deallocate our old values
        deallocate(this % electron_integral)
        deallocate(this % coeffs)

        !Increase the max capacity
        this % max_capacity = this % max_capacity + this % expand_size

        !Allocate our new array with the new max capacity
        allocate(this % electron_integral(NIDX, this % max_capacity), &
                 this % coeffs(this % num_states_1, this % num_states_2, this % max_capacity), stat = err)
        call master_memory % track_memory (storage_size(this % electron_integral)/8, &
                                           size(this % electron_integral), err,'CONSYMBOLIC::EXP::ELECTRONINT_EXP')
        call master_memory % track_memory(storage_size(this % coeffs)/8, &
                                          size(this % coeffs), err, 'CONSYMBOLIC::EXP::ELECTRONCOEFF_EXP')
        !Zero the elements
        this % electron_integral = -1
        this % coeffs = 0.0_wp

        do ido = 1, temp_capacity
            !Copy our old values back
            this % electron_integral(:,ido) = temp_integral(:,ido)
            this % coeffs(:,:,ido)          = temp_coeffs(:,:,ido)
        end do

        call master_memory % free_memory(storage_size(temp_integral)/8, size(temp_integral))
        call master_memory % free_memory(storage_size(temp_coeffs)/8, size(temp_coeffs))

        deallocate(temp_integral, temp_coeffs)

    end subroutine expand_array


    !> \brief   Inserts one symbolic vector into another scaled by
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this   Vector object to update.
    !> \param[in]    rhs    The symbolic vector to add into this class
    !> \param[in]    alpha  The scaling value to be applied to each coefficient
    !>
    subroutine    add_symbols (this, rhs, alpha)
        class(ContractedSymbolicElementVector)   :: this
        class(SymbolicElementVector), intent(in) :: rhs
        real(wp)         :: alpha(this % num_states_1, this % num_states_2), int_coeff
        real(wp)         :: coeffs(this % num_states_1, this % num_states_2)
        integer          :: ido
        integer(longint) :: label(NIDX)

        coeffs = 0._wp

        !Loop through each rhs symbolic and add it to me scaled by alpha
        do ido = 1, rhs % get_size()
            call rhs % get_coeff_and_integral(ido, int_coeff, label)
            coeffs = int_coeff * alpha
            call this % insert_symbol(label, coeffs)
        end do

    end subroutine add_symbols


    subroutine reduce_symbols (this, rhs)
        class(ContractedSymbolicElementVector)             :: this
        class(ContractedSymbolicElementVector), intent(in) :: rhs
        integer          :: ido

        !Loop through each rhs symbolic and add it to me scaled by alpha
        do ido = 1, rhs % get_size()
            call this % insert_symbol(rhs % electron_integral(:,ido), rhs % coeffs(:,:,ido))
        end do

    end subroutine reduce_symbols


    subroutine synchronize_symbols (this)
        class(ContractedSymbolicElementVector) :: this

        integer :: my_num_symbols, largest_num_symbols, ierr
        integer :: ido, proc_id, jdo

        integer(longint), allocatable, target :: labels(:,:)
        integer(longint), pointer             :: label_ptr(:)
        integer(longint)                      :: procs_num_of_symbols_int, procs_num_of_symbols_coeff, size_labels, size_coeffs

        real(wp), allocatable, target :: coeffs(:,:,:)
        real(wp), pointer             :: coeffs_ptr(:)

        if (grid % gprocs <= 1) then
            return
        end if

        call master_timer % start_timer("Symbol Synchronize")

        my_num_symbols = this % get_size()
        call mpi_reduceall_max(my_num_symbols, largest_num_symbols, grid % lcomm)

        if (largest_num_symbols == 0) then
            call master_timer % stop_timer("Symbol Synchronize")
            return
        end if

        call master_memory % track_memory (storage_size(labels)/8, largest_num_symbols * NIDX, 0, 'CONSYMBOLIC::MPISYNCH::LABELS')
        call master_memory % track_memory (storage_size(coeffs)/8, largest_num_symbols * this % num_states_1 * this % num_states_2,&
                                            0, 'CONSYMBOLIC::MPISYNCH::COEFFS')

        allocate(labels(NIDX, largest_num_symbols), &
                 coeffs(this % num_states_1, this % num_states_2, largest_num_symbols), &
                 stat = ierr)

        size_labels = size(labels, kind = longint)
        size_coeffs = size(coeffs, kind = longint)

        label_ptr(1:size_labels) => labels(:,:)
        coeffs_ptr(1:size_coeffs) => coeffs(:,:,:)

        labels = 0
        coeffs = 0

        if (my_num_symbols > 0) then
            labels(:,1:my_num_symbols) = this % electron_integral(:,1:my_num_symbols)
            do ido = 1, my_num_symbols
                coeffs(1:this%num_states_1,1:this%num_states_2,ido) = this % coeffs(1:this%num_states_1,1:this%num_states_2,ido)
            end do
        end if

        procs_num_of_symbols_int = my_num_symbols * NIDX
        procs_num_of_symbols_coeff = my_num_symbols * this % num_states_1 * this % num_states_2
        do proc_id = 1, grid % lprocs - 1
            call mpi_mod_rotate_int_arrays_around_ring(procs_num_of_symbols_int, label_ptr, size_labels, grid % lcomm)
            call mpi_mod_rotate_cfp_arrays_around_ring(procs_num_of_symbols_coeff, coeffs_ptr, size_coeffs, grid % lcomm)

            do ido = 1, procs_num_of_symbols_int / NIDX
                call this % insert_symbol(labels(:,ido), coeffs(1:this%num_states_1,1:this%num_states_2,ido))
            end do
        end do

        call master_timer % stop_timer("Symbol Synchronize")
        call master_memory % free_memory(storage_size(labels)/8, size(labels))
        call master_memory % free_memory(storage_size(coeffs)/8, size(coeffs))

        deallocate(labels, coeffs)

    end subroutine synchronize_symbols


    subroutine synchronize_symbols_II (this)
        class(ContractedSymbolicElementVector), target :: this

        integer(longint), allocatable, target :: labels(:,:)
        integer(longint), pointer             :: label_ptr(:)

        real(wp), allocatable, target :: coeffs(:,:,:)
        real(wp), pointer             :: coeffs_ptr(:)

        integer :: my_num_symbols, largest_num_symbols, procs_num_of_symbols_int, procs_num_of_symbols_coeff, ido, proc_id, jdo
        integer :: ierr, local_communicator, global_communicator, odd_even

        integer(kind=mpiint) :: master_comm, global_rank, global_nprocs, error, proc, tag = 1

        if (grid % gprocs <= 1) then
            return
        end if

    end subroutine synchronize_symbols_II


    !> \brief   Simply checks the index and wheter it exists within the array
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this  Vector object to update.
    !> \param[in]    i     The index we wish to access.
    !>
    logical function check_bounds (this, i)
        class(ContractedSymbolicElementVector) :: this
        integer, intent(in)                    :: i

        if (i <= 0 .or. i > this % n) then
            write (stdout, "('SymbolicVector::check_bounds - Out of Bounds access')")
            stop "SymbolicVector::check_bounds - Out of Bounds access"
            check_bounds = .false.
        else
            check_bounds = .true.
        end if

    end function check_bounds


    !> \brief   Removes an integral and coefficient, never been used and pretty much useless.
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine remove_symbol_at (this, idx)
        class(ContractedSymbolicElementVector) :: this
        integer, intent(in) :: idx
        integer             :: ido

        if (this % check_bounds(idx)) then
            do ido = idx + 1, this % n
                this % electron_integral(:, ido - 1) = this % electron_integral(:, ido)
                !this%coeffs(ido-1)        = this%coeffs(ido)
            end do
            this % electron_integral(:, this % n) = 0
           !this % coeffs(this % n) = 0
            this % n = this % n - 1
        end if

    end subroutine remove_symbol_at


    !> @brief   Simply returns whether we are storing integrals and coeffs or not
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \result True if empty or False if we have any symbols
    !>
    logical function is_empty (this)

        class(ContractedSymbolicElementVector) :: this

        is_empty = (this % n == 0)

    end function is_empty


    !> \brief   Clear our array (not really but it essentialy resets the symbol counter to zero which is way quicker).
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine clear (this)

        class(ContractedSymbolicElementVector) :: this

        this % n = 0
        call this % bstree % destroy

    end subroutine clear


    subroutine modify_coeff (this, idx, coeff)
        class(ContractedSymbolicElementVector) :: this
        integer, intent(in) :: idx
        real(wp) :: coeff
        !this%coeffs(idx) = this%coeffs(idx)+ coeff
    end subroutine modify_coeff


    !> @brief   Get integral label at specific index
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    function get_integral_label (this, idx)
        class(ContractedSymbolicElementVector) :: this
        integer, intent(in)                    :: idx
        integer(longint)                       :: get_integral_label(NIDX)

        if (this % check_bounds(idx)) get_integral_label(:) = this % electron_integral(:,idx)

    end function get_integral_label


    !> \brief   Get coefficient at specific index
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    real(wp) function get_coefficient(this, idx, n1, n2)
        class(ContractedSymbolicElementVector) :: this
        integer, intent(in) :: idx, n1, n2

        !if(this%check_bounds(idx)==.true.) get_coefficient=this%coeffs(idx)
        get_coefficient = this % coeffs(n1, n2, idx)

    end function get_coefficient


    !> \brief   Get both label and coeffcient at specific index
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine get_coeff_and_integral (this, idx, coeff, label)
        class(ContractedSymbolicElementVector), intent(in) :: this
        integer,          intent(in)  :: idx
        real(wp),         intent(out) :: coeff(this % num_states_1, this % num_states_2)
        integer(longint), intent(out) :: label(:)

        if (this % check_bounds(idx)) then
            label = this % electron_integral(:,idx)
            coeff = this % coeffs(:,:,idx)
        end if

    end subroutine get_coeff_and_integral


    !> \brief   Returns the number of symbolic elements stored
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer function get_size (this)
        class(ContractedSymbolicElementVector) :: this
        get_size = this % n
    end function get_size


    integer function get_num_targets_sym1 (this)
        class(ContractedSymbolicElementVector) :: this
        get_num_targets_sym1 = this % num_states_1
    end function get_num_targets_sym1


    integer function get_num_targets_sym2 (this)
        class(ContractedSymbolicElementVector) :: this
        get_num_targets_sym2 = this % num_states_2
    end function get_num_targets_sym2


    !> \brief   Cleans up the class by deallocating arrays
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine destroy (this)
        class(ContractedSymbolicElementVector) :: this

        if (allocated(this % electron_integral)) then
            call master_memory % free_memory(storage_size(this % electron_integral)/8, size(this % electron_integral))
            deallocate(this % electron_integral)
        end if
        if(allocated(this % coeffs)) then
            call master_memory % free_memory(storage_size(this % coeffs)/8, size(this % coeffs))
            deallocate(this % coeffs)
        end if

        call this % bstree % destroy

        this % constructed = .false.

    end subroutine destroy


    !> \brief   Print currently stored symbols
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine print_symbols (this)
        class(ContractedSymbolicElementVector) :: this
        integer  :: labels(8), ido
        real(wp) :: coeff

        if (.not. this % is_empty()) write (stdout, "('Outputting symbolic elements....')")

        do ido = 1, this % n
            call unpack_ints(this % electron_integral(:,ido), labels)
            write (1923 + myrank, "(5i4,' -- ',5es14.3)") labels(1:5), this % coeffs(:,:,ido)
        end do

    end subroutine print_symbols


    logical function estimate_synchronize_cost (this)
        class(ContractedSymbolicElementVector) :: this
        integer(longint) :: memory_cost
        integer          :: my_num_symbols, largest_num_symbols
        integer          :: gathered_bool, global_bool(grid % gprocs)

        if (grid % gprocs <= 1) then
            estimate_synchronize_cost = .true.
            return
        end if

        my_num_symbols = this % get_size()
        gathered_bool = 0
        estimate_synchronize_cost = .true.

        call mpi_reduceall_max(my_num_symbols, largest_num_symbols, grid % gcomm)

        memory_cost = 2 * largest_num_symbols * (2 + this % num_states_1 * this % num_states_2) * 8

        if (memory_cost >= master_memory % get_scaled_available_memory(0.75_wp)) gathered_bool = 1

        global_bool = 0

        call mpi_mod_allgather(gathered_bool, global_bool, grid % gcomm)

        if (sum(global_bool) /= 0) then
            estimate_synchronize_cost = .false.
        end if

    end function estimate_synchronize_cost

end module Contracted_Symbolic_module
