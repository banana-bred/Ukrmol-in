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

!> \brief   SLEPc matrix module
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module is built in only when SLEPc library is available.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module SLEPCMatrix_module

    use const_gbl,                only: stdout
    use consts_mpi_ci,            only: MAT_DENSE, MAT_SPARSE
    use precisn,                  only: wp
    use mpi_gbl,                  only: mpi_mod_allgather, mpi_mod_scan_sum, mpi_xermsg
    use DistributedMatrix_module, only: DistributedMatrix
    use MatrixCache_module,       only: MatrixCache
    use MemoryManager_module,     only: master_memory
    use Parallelization_module,   only: grid => process_grid
    use Timing_Module,            only: master_timer
    use slepcsys,                 only: SlepcInitialize
    use petsc

#include <finclude/petsc.h>

    implicit none

    public SLEPCMatrix, initialize_slepc

    private

    type :: CSRFormat
        PetscInt                 :: row
        PetscInt,    allocatable :: col(:)
        PetscScalar, allocatable :: coeff(:)
        PetscInt                 :: num_cols
    contains
        procedure, public :: construct => construct_csr
        procedure, public :: sort      => sort_csr
        final :: destroy_csr
    end type CSRFormat

    type, extends(DistributedMatrix) :: SLEPCMatrix

        type(MatrixCache)  :: slepc_cache
        integer            :: start_row
        integer            :: end_row
        integer            :: local_size
        integer            :: mem_track
        Mat, pointer       :: A

        ! These are required to ensure that both PETSC and the local MPI configuration behave
        PetscInt           :: start_row_PETSC
        PetscInt           :: end_row_PETSC
        PetscInt           :: local_size_PETSC
        PetscInt           :: matrix_dimension_PETSC

        !> These are required for preallocation
        PetscInt, allocatable :: diagonal_nonzero(:)
        PetscInt, allocatable :: offdiagonal_nonzero(:)

        !This stores which row belongs to which process
        integer :: diagonal_start,diagonal_end
        logical :: compressed = .false.
        logical :: mat_created = .false.

    contains

        procedure, public   :: print                => print_SLEPC
        procedure, public   :: setup_diag_matrix    => initialize_struct_SLEPC
        procedure, public   :: get_matelem_self     => get_matelem_SLEPC
        procedure, public   :: clear_matrix         => clear_SLEPC
        procedure, public   :: destroy_matrix       => destroy_SLEPC
        procedure, public   :: finalize_matrix_self => finalize_matrix_SLEPC

        procedure, public   :: insert_into_diag_matrix => insert_into_hard_cache
        procedure, private  :: compress_cache_to_csr_format
        procedure, private  :: insert_csr_into_hard_cache
        procedure, public   :: print_nonzeros
        procedure, public   :: get_PETSC_matrix
        procedure, private  :: create_PETSC_mat
        procedure, private  :: destroy_PETSC_mat

    end type SLEPCMatrix

contains

    !> \brief   Initialize SLEPc
    !> \authors J Benda
    !> \date    2019
    !>
    !> SLEPc needs to be initialized by all processes at once. Originally, it was being initialized in the SLEPc matrix init
    !> routine, but since there are now multiple concurrent diagonalizations and not all of them need to use SLEPc, this needs
    !> to be separated out.
    !>
    subroutine initialize_slepc

        PetscErrorCode :: ierr

        call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
        write (stdout, "('SLEPc initialization status = ',I0)") ierr
        if (ierr /= 0) then
            call mpi_xermsg('SLEPCMatrix_module', 'initialize_slepc', 'Failed to initialize SLEPc', int(ierr), 1)
        end if

    end subroutine initialize_slepc


    subroutine construct_csr (this, num_cols)
        class(CSRFormat)    :: this
        integer, intent(in) :: num_cols

        if (allocated(this % coeff)) deallocate(this % coeff)
        if (allocated(this % col)) deallocate(this % col)

        allocate(this % coeff(num_cols))
        allocate(this % col(num_cols))

    end subroutine construct_csr


    subroutine destroy_csr (this)
        type(CSRFormat) :: this

        if (allocated(this % col)) deallocate(this % col)
        if (allocated(this % coeff)) deallocate(this % coeff)

    end subroutine destroy_csr


    subroutine sort_csr (this)
        class(CSRFormat) :: this

        call QsortCSR (this % col, this % coeff)

    end subroutine sort_csr


    function get_PETSC_Matrix (this) result(m)
        class(SLEPCMatrix) :: this
        Mat, pointer       :: m

        m => this % A

    end function get_PETSC_Matrix


    subroutine create_PETSC_mat (this, matrix_size, matrix_type)
        class(SLEPCMatrix)   :: this
        PetscInt, intent(in) :: matrix_size
        integer,  intent(in) :: matrix_type
        integer              :: ierr, non_zero_guess

        if (this % mat_created) then
            stop "PETSC matrix already created!!"
        endif

        allocate(this % A)

        if (matrix_type == MAT_DENSE) then

        else if (matrix_type == MAT_SPARSE) then

        end if

        this % mat_created = .true.

    end subroutine create_PETSC_mat


    subroutine destroy_PETSC_mat (this)
        class(SLEPCMatrix) :: this
        PetscErrorCode     :: ierr

        if (this % mat_created) then
            call MatDestroy(this % A, ierr)
            deallocate(this % A)
            this % A => null()
            call master_memory % free_memory(this % mem_track, 1)
            this % mem_track = 0
            this % mat_created = .false.
        end if

    end subroutine destroy_PETSC_mat


    subroutine Shuffle (a)
        PetscInt, intent(inout) :: a(:)

        integer :: i, randpos, temp
        real    :: r

        do i = size(a), 2, -1
            call random_number(r)
            randpos = int(r * i) + 1
            temp = a(randpos)
            a(randpos) = a(i)
            a(i) = temp
        end do

    end subroutine Shuffle


    subroutine initialize_struct_SLEPC (this, matrix_size, matrix_type, block_size)

        class(SLEPCMatrix)  :: this
        integer, intent(in) :: matrix_size, matrix_type, block_size

        PetscInt            :: non_zero_guess, petsc_row, to_insert, one = 1, dummy_int, matrix_size_PETSC
        PetscErrorCode      :: ierr
        real(wp)            :: mem_usage_begin, mem_usage_end, ratio, dummy_real
        integer             :: ido, total_vals, per_elm, l_update, c_update, indicies, total_elems

        integer,     allocatable :: mem_use_total(:)
        PetscInt,    allocatable :: test_insertion_ij(:)
        PetscScalar, allocatable :: test_insertion_val(:)

        !If we are dealing with a dense matrix then we don't bother since SLEPC automatically handles this
        ierr = 0

        call this % slepc_cache % construct
        call this % destroy_PETSC_mat

        this % matrix_dimension_PETSC = this % matrix_dimension
        write (stdout, "('Mat_size = ',i8)") this % matrix_dimension_PETSC

        call this % create_PETSC_mat(this % matrix_dimension_PETSC, matrix_type)

        matrix_size_PETSC = matrix_size
        this % local_size_PETSC = PETSC_DECIDE

        !Lets find our local dimension size
        call PetscSplitOwnership(int(grid % gcomm, kind(PETSC_COMM_WORLD)), &
                                 this % local_size_PETSC, this % matrix_dimension_PETSC, ierr)

        !Convert to our local integer
        this % local_size = this % local_size_PETSC

        write (stdout, "('Local_size = ',i8)") this % local_size
        write (stdout, "('THRESHOLD = ',es16.8)") this % threshold

        !Get the end index
        call mpi_mod_scan_sum(this % local_size, this % end_row, grid % gcomm)
        !Get the start index
        this % start_row = this % end_row - this % local_size

        !store the rest just in case petsc needs it
        this % start_row_PETSC = this % start_row
        this % end_row_PETSC = this % end_row
        write (stdout, "('start,end = ',2i8)") this % start_row, this % end_row

        !Dense doesn';t require much
        if (this % matrix_type == MAT_DENSE) then
            call PetscMemoryGetCurrentUsage(mem_usage_begin, ierr)
            call MatCreateDense(grid % gcomm, PETSC_DECIDE, PETSC_DECIDE, matrix_size_PETSC, matrix_size_PETSC, &
                                PETSC_NULL_SCALAR_ARRAY, this % A, ierr)
            call PetscMemoryGetCurrentUsage(mem_usage_end, ierr)

            write (stdout, "('Matrix uses ',f8.3,' GB of memory')") (mem_usage_end - mem_usage_begin) / 1024**3
            this % mem_track = mem_usage_end - mem_usage_begin

            call master_memory % track_memory(this % mem_track, 1, 0, 'PETSC MATRIX')
            return
        end if

        !Now we need to allocate what is needed
        if (allocated(this % diagonal_nonzero)) deallocate(this % diagonal_nonzero)
        if (allocated(this % offdiagonal_nonzero)) deallocate(this % offdiagonal_nonzero)
        !The final procedure should take care of all the arrays
        !This will store our non-zeros for the local matrix
        allocate(this % diagonal_nonzero(this % local_size), stat = ierr)
        allocate(this % offdiagonal_nonzero(this % local_size), stat = ierr)

        this % diagonal_nonzero = 0
        this % offdiagonal_nonzero = 0

        total_elems = 0

        this % diagonal_start = this % start_row_PETSC
        this % diagonal_end = min(this % start_row_PETSC + this % local_size, this % matrix_dimension)

        do ido = 1, this % local_size
            total_vals = this % matrix_dimension - this % start_row - ido + 1
            this % diagonal_nonzero(ido) = this % diagonal_end - this % diagonal_start + 1 - ido
            this % offdiagonal_nonzero(ido) = total_vals - this % diagonal_nonzero(ido)

            !If we're above the block size then we use the 90% sparse guess
            if ((this % start_row + ido - 1) >= block_size) then
                !Now we are dealing with sparsity
                !Lets compute the ratio between the diagonal and off diagonal
                ratio = real(this % diagonal_nonzero(ido)) / real(total_vals)
                !Reduce the number of values to 30% of max
                total_vals = real(total_vals) * 0.1
                !Now distribute
                !There should be at least one in the diagonal
                !The off diagonal will be larger so we will give the diagonal a bit more leeway
                this % diagonal_nonzero(ido) = max(1, int(real(total_vals) * ratio))
                this % offdiagonal_nonzero(ido) = int(real(total_vals) * (1.0 - ratio))
            end if
            total_elems = total_elems + this % diagonal_nonzero(ido) + this % offdiagonal_nonzero(ido)
            !write(stdout,"(5i8)") this%start_row,this%end_row,total_vals,this%diagonal_nonzero(ido),this%offdiagonal_nonzero(ido)
        end do

        call PetscMemoryGetCurrentUsage(mem_usage_begin, ierr)

        this % continuum_counter = 0
        this % L2_counter = 0

        call MatCreateSBAIJ(grid % gcomm, one, PETSC_DECIDE, PETSC_DECIDE, matrix_size_PETSC, matrix_size_PETSC, &
                            PETSC_DEFAULT_INTEGER, this % diagonal_nonzero, PETSC_DEFAULT_INTEGER, this % offdiagonal_nonzero, &
                            this % A, ierr)

        !call MatCreate(PETSC_COMM_WORLD,this%A,ierr)
        !CHKERRQ(ierr)
        !call MatSetType(this%A,MATSBAIJ,ierr)
        !CHKERRQ(ierr)
        !call MatSetSizes(this%A,PETSC_DECIDE,PETSC_DECIDE,matrix_size,matrix_size,ierr)
        !CHKERRQ(ierr)
                !Now we set the preallocation
        !if(nprocs > 1) then
        !    call MatMPISBAIJSetPreallocation(this%A,1,PETSC_DECIDE,this%diagonal_nonzero,PETSC_DECIDE,this%offdiagonal_nonzero,ierr)
        !else
        !    call MatSeqSBAIJSetPreallocation(this%A,1,PETSC_DECIDE,this%diagonal_nonzero,ierr)
        !endif
        !call MatSetUp(this%A,ierr)
        !call MatXAIJSetPreallocation(this%A, 1,PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,this%diagonal_nonzero, this%offdiagonal_nonzero,ierr)
        if (ierr /= 0) then
            call mpi_xermsg('SLEPCMatrix_module', 'initialize_struct_SLEPC', 'Failed to create matrix', int(ierr), 1)
        end if
        !call MatSetUp(this%A,ierr)

        call MatSetOption(this % A, MAT_IGNORE_OFF_PROC_ENTRIES,    PETSC_TRUE,  ierr)
        call MatSetOption(this % A, MAT_NO_OFF_PROC_ENTRIES,        PETSC_TRUE,  ierr)
        call MatSetOption(this % A, MAT_IGNORE_LOWER_TRIANGULAR,    PETSC_TRUE,  ierr)
        call MatSetOption(this % A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)

        deallocate(this % diagonal_nonzero, this % offdiagonal_nonzero)

        allocate(test_insertion_ij(this % matrix_dimension), test_insertion_val(this % matrix_dimension))

        do ido = 1, this % matrix_dimension
            test_insertion_ij(ido)  = ido - 1
            test_insertion_val(ido) = 1.0_wp
        end do

        !--------We use this to embed the non-zero structce into the block to improve performance when inserting values
        call master_timer % start_timer("Insertion time")
        do ido = 1, this % local_size
            total_vals = this % matrix_dimension - this % start_row - ido + 1
            petsc_row = ido - 1
            dummy_int = this % start_row + ido
            if ((dummy_int - 1) < block_size) then
                to_insert = size(test_insertion_ij(dummy_int:this % matrix_dimension))
                !write(stdout,*) ido,to_insert,dummy_int,total_vals
                !call MatSetValuesBlocked(this%A,1,ido,csr(ido)%num_cols,csr(ido)%col,csr(ido)%coeff,INSERT_VALUES,ierr)
                call MatSetValuesBlocked(this % A, one, [dummy_int - one], to_insert, &
                                         test_insertion_ij(dummy_int : this % matrix_dimension), &
                                         test_insertion_val(dummy_int : this % matrix_dimension), &
                                         INSERT_VALUES, ierr)
            else
                exit
            end if
        end do
        call MatAssemblyBegin(this % A, MAT_FLUSH_ASSEMBLY, ierr)
        call MatAssemblyEnd(this % A, MAT_FLUSH_ASSEMBLY, ierr)

        !Now lets clear just in case we have any issues
        call MatZeroEntries(this % A, ierr)


        deallocate(test_insertion_ij, test_insertion_val)

        write (stdout, "('PETSC matrix creation completed!')")

        call master_timer % stop_timer("Insertion time")

        call PetscMemoryGetCurrentUsage(mem_usage_end, ierr)

        per_elm = (storage_size(test_insertion_ij) + storage_size(test_insertion_val)) / 8

        write (stdout, "('Matrix uses ',f0.0,' B of memory')") mem_usage_end - mem_usage_begin
        write (stdout, "('Estimated uses ',i0,' B of memory')") total_elems * per_elm

        this % mem_track = mem_usage_end - mem_usage_begin

        allocate(mem_use_total(grid % gprocs))

        call mpi_mod_allgather(total_elems, mem_use_total, grid % gcomm)

        total_elems = sum(mem_use_total)
        total_elems = total_elems / grid % gprocs
        this % mem_track = total_elems

        deallocate(mem_use_total)

        write (stdout, "('Average uses ',i0,' B of memory')") total_elems*per_elm

        call master_memory % track_memory(per_elm, total_elems, 0, 'PETSC MATRIX')

    end subroutine initialize_struct_SLEPC


    !TODO: Write a 'compressed' version of this
    subroutine get_matelem_SLEPC (this, idx, i, j, coeff)
        class(SLEPCMatrix)    :: this
        integer,  intent(in)  :: idx
        integer,  intent(out) :: i, j
        real(wp), intent(out) :: coeff

        i = 0
        j = 0
        coeff = 0.0
        !call this%matrix_cache%get_from_cache(idx,i,j,coeff)

    end subroutine get_matelem_SLEPC


    subroutine print_nonzeros (this)
        class(SLEPCMatrix) :: this
        integer            ::    ido

        do ido = 1, this % local_size
            write (stdout, "('dnnz = ',i8,'onnz = ',i8)") this % diagonal_nonzero(ido), this % offdiagonal_nonzero(ido)
        end do

    end subroutine print_nonzeros


    logical function compress_cache_to_csr_format (this, matrix_cache, csr_matrix)
        class(SLEPCMatrix)                        :: this
        type(CSRFormat), allocatable, intent(out) :: csr_matrix(:)
        class(MatrixCache)                        :: matrix_cache

        integer  :: ido, jdo, start_idx, end_idx, num_elems, row_count, row_idx, first_row
        integer  :: current_chunk, current_row, next_row, last_chunk, i, j
        real(wp) :: coeff

        compress_cache_to_csr_format = .false.

        !Lets not bother if its a DENSE matrix
        if (this % matrix_type == MAT_DENSE) return

        !sort the matrix cache in ascending row order
        call matrix_cache % sort

        current_row = -1
        row_count = 0
        do ido = 1, matrix_cache % get_size()
            call matrix_cache % get_from_cache(ido, i, j, coeff)
            if (ido == 1) first_row = i
            if (i /= current_row) then
                row_count = row_count + 1
                current_row = i
            end if
        end do

        call master_timer % start_timer("CSR conversion")

        num_elems = matrix_cache % get_size()

        if (row_count == 0) then
            call master_timer % stop_timer("CSR conversion")
            return
        end if

        if (allocated(csr_matrix)) deallocate(csr_matrix)

        !allocate the number of local rows
        allocate(csr_matrix(row_count))

        csr_matrix(:) % num_cols = 0

        next_row = -1
        last_chunk = 1
        start_idx = 1
        end_idx = 99
        ido = 1
        row_idx = 1
        current_row = first_row

        !Loop through number of elements
        do while (ido <= num_elems)
            row_count = 0

            !counting phase
            do jdo = ido, num_elems
                call matrix_cache % get_from_cache(jdo, i, j, coeff)

                !Check row and threshold count
                if (i == current_row) then
                    if (abs(coeff) > this % threshold) row_count = row_count + 1
                    end_idx = jdo
                else
                    next_row = i
                    exit
                end if
            end do
            csr_matrix(row_idx) % row = current_row
            csr_matrix(row_idx) % num_cols = row_count

            !Do we have anything to store?
            if (row_count > 0) then
                !allocation phase
                call csr_matrix(row_idx) % construct(row_count)

                !now we loop through and store the elements
                do jdo = 1, row_count
                    !Get the element
                    call matrix_cache % get_from_cache(jdo + start_idx - 1, i, j, coeff)
                    !add it to the list
                    csr_matrix(row_idx) % col(jdo) = j
                    csr_matrix(row_idx) % coeff(jdo) = coeff
                    current_chunk = matrix_cache % get_chunk_idx(jdo + start_idx - 1)
                    !also if we've moved to the next chunk then delete the lastone
                    if (last_chunk /= current_chunk) then
                        call matrix_cache % matrix_arrays(last_chunk) % destroy
                        last_chunk = current_chunk
                    end if
                end do
            end if

            !If we've hit the last element previously then we can leave
            ido = end_idx + 1
            start_idx = end_idx + 1
            row_idx = row_idx + 1
            current_row  = next_row
        end do

        compress_cache_to_csr_format = .true.

        !After this we can delete anything left over
        call matrix_cache % destroy
        call matrix_cache % construct
        call master_timer % stop_timer("CSR conversion")

    end function compress_cache_to_csr_format


    !> \brief Inserts an element into the hard storage which is considered the final location before diagonalization
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> It also checks wherther the element exists within the aloowed range and tells us if it was successfully inserted.
    !>
    logical function insert_into_hard_cache (this, row, column, coefficient)
        class(SLEPCMatrix)   :: this
        integer,  intent(in) :: row, column
        real(wp), intent(in) :: coefficient
        PetscInt             :: i, j
        PetscErrorCode       :: ierr

        !if(row==column) call this%store_diagonal(row,coefficient)

        i = row - 1
        j = column - 1

        if (abs(coefficient) < this % threshold) return

        !Dense has slightly different rules for this
        if (this % matrix_type == MAT_DENSE) then
            if (i >= this % start_row .or. i < this % end_row ) then
                call MatSetValue(this % A, i, j, coefficient, INSERT_VALUES, ierr)
            end if
            if (j >= this % start_row .or. j < this % end_row) then
                call MatSetValue(this % A, j, i, coefficient, INSERT_VALUES, ierr)
                insert_into_hard_cache = .true.
            end if
            return
        end if

        if (j < this % start_row .or. j >= this % end_row) then
            insert_into_hard_cache = .false.
            return
        end if

        insert_into_hard_cache = .true.
        !write(stdout,"(2i8,1)") row,column
        !call this%slepc_cache%insert_into_cache(column-1,row-1,coefficient)

        call MatSetValue(this % A, j, i, coefficient, INSERT_VALUES, ierr)

        this % n = this % n + 1

    end function insert_into_hard_cache


    subroutine insert_csr_into_hard_cache (this, csr)
        class(SLEPCMatrix)            :: this
        class(CSRFormat),  intent(in) :: csr(:)
        integer        :: size_csr
        PetscInt       :: row, column, ido, one = 1
        PetscErrorCode :: ierr

        if (this % matrix_type == MAT_DENSE) return

        size_csr = size(csr)

        do ido = 1, size_csr
            !call csr(ido)%sort
            row = csr(ido) % row
            if (row < this % start_row .or. row >= this % end_row) cycle
            if (csr(ido) % num_cols == 0) cycle
            !write(stdout,"('Inserting row ',2i8)") row,csr(ido)%num_cols

            call MatSetValuesBlocked(this % A, one, [csr(ido) % row], csr(ido) % num_cols, &
                                     csr(ido) % col, csr(ido) % coeff, INSERT_VALUES, ierr)

            this % n = this % n + csr(ido) % num_cols
        end do

    end subroutine insert_csr_into_hard_cache


    subroutine finalize_matrix_SLEPC (this)
        class(SLEPCMatrix) :: this
        PetscErrorCode     :: ierr

        if (.not. this % mat_created) then
            stop "Finalizing matrix that doesn't exist!!!"
        else
            write (stdout, "('Finalizing SLEPC matrix')")
            call master_timer % start_timer("Matrix assembly")
            call MatAssemblyBegin(this % A, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(this % A, MAT_FINAL_ASSEMBLY, ierr)
            !call MatSetOption(this%A,MAT_SYMMETRIC ,PETSC_TRUE,ierr)
            !call MatSetOption(this%A,MAT_HERMITIAN ,PETSC_TRUE,ierr)
            !call MatSetOption(this%A,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr)
            call master_timer % stop_timer("Matrix assembly")
            write (stdout, "('Finished assembly')")
        end if

    end subroutine finalize_matrix_SLEPC


    subroutine print_SLEPC (this)
        class(SLEPCMatrix) :: this

        write (stdout, "('-------TEMP CACHE---------')")
        call this % temp_cache % print

    end subroutine print_SLEPC


    subroutine clear_SLEPC (this)
        class(SLEPCMatrix) :: this

        this % n = 0

    end subroutine clear_SLEPC


    subroutine destroy_SLEPC (this)
        class(SLEPCMatrix) :: this

        call this % clear
        call this % destroy_PETSC_mat

    end subroutine destroy_SLEPC

    recursive subroutine QsortCSR (A, coeff)
        PetscInt,    intent(inout), dimension(:) :: A
        PetscScalar, intent(inout), dimension(:) :: coeff
        integer :: iq

        if (size(A) > 1) then
            call Partition(A, coeff, iq)
            call QsortCSR(A(:iq-1), coeff(:iq-1))
            call QsortCSR(A(iq:), coeff(iq:))
        end if

    end subroutine QsortCSR


    subroutine Partition (A, coeff, marker)
        PetscInt,    intent(inout), dimension(:) :: A
        PetscScalar, intent(inout), dimension(:) :: coeff
        integer, intent(out) :: marker
        integer     :: i, j
        PetscInt    :: temp_A
        PetscScalar :: temp_coeff
        PetscInt    :: x_A      ! pivot point

        x_A = A(1)
        i = 0
        j = size(A) + 1

        do
            j = j - 1
            do
                if (A(j) <= x_A) exit
                j = j - 1
            end do
            i = i + 1
            do
                if (A(j)>= x_A) exit
                i = i + 1
            end do
            if (i < j) then
                ! exchange A(i) and A(j)
                temp_A     = A(i);     A(i)     = A(j);     A(j)     = temp_A
                temp_coeff = coeff(i); coeff(i) = coeff(j); coeff(j) = temp_coeff
            else if (i == j) then
                marker = i + 1
                return
            else
                marker = i
                return
            end if
        end do

    end subroutine Partition

end module SLEPCMatrix_module
