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

!> \brief   Distributed matrix module
!> \authors A Al-Refaie
!> \date    2017
!>
!> Provides DistributedMatrix used by parallel diagonalizers. Other specialized matrix types
!> are based on this. See \ref SCALAPACKMatrix_module::SCALAPACKMatrix and
!> \ref SLEPCMatrix_module::SLEPCMatrix.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module DistributedMatrix_module

    use const_gbl,              only: stdout
    use consts_mpi_ci,          only: MAT_DENSE
    use precisn,                only: longint, wp
    use mpi_gbl,                only: nprocs, mpi_xermsg, mpi_reduceall_min, mpi_reduceall_max, &
                                      mpi_mod_rotate_arrays_around_ring, mpi_reduceall_inplace_sum_cfp, mpi_reduceall_sum_cfp
    use BaseMatrix_module,      only: BaseMatrix
    use MatrixCache_module,     only: MatrixCache
    use MemoryManager_module,   only: master_memory
    use Parallelization_module, only: grid => process_grid
    use Timing_Module,          only: master_timer

    implicit none

    public DistributedMatrix

    private

    !> \brief   Distributed matrix type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Distributed matrix is used for construction of the Hamiltonian in distributed (MPI) mode.
    !>
    type, abstract, extends(BaseMatrix) :: DistributedMatrix
        !>This cache holds the temprorary matrix values before they are sent to the relevant process
        type(MatrixCache) :: temp_cache

        real(wp) :: memory_scale = 0.75_wp

        integer  :: continuum_counter
        integer  :: L2_counter
        integer  :: start_continuum_update
        integer  :: start_L2_update

    contains

        procedure, public   :: print => print_Distributed
        procedure, public   :: update_continuum => update_continuum_Distributed
        procedure, public   :: update_pure_L2   => update_L2_Distributed

        procedure, public   :: initialize_struct_self  => initialize_struct_Distributed

        procedure, public   :: setup_diag_matrix

        procedure, public   :: construct_self      => construct_mat_Distributed
        procedure, public   :: insert_matelem_self => insert_matelem_Distributed
        procedure, public   :: get_matelem_self    => get_matelem_Distributed
        procedure, public   :: clear_self          => clear_Distributed
        procedure, public   :: destroy_self        => destroy_Distributed
        procedure, public   :: finalize_matrix     => finalize_Distributed
        procedure, public   :: finalize_matrix_self
        procedure, public   :: destroy_matrix
        procedure, public   :: clear_matrix

        procedure, public   :: insert_into_diag_matrix
        procedure, private  :: convert_temp_cache_to_array
        procedure, private  :: update_counter
    end type

contains

    subroutine construct_mat_Distributed (this)
        class(DistributedMatrix) :: this

        write (stdout, "('Constructing Distributed matrix')")

        !Initialize the temporary cache
        call this % temp_cache % construct

        !Clear them both
        call this % temp_cache % clear

    end subroutine construct_mat_Distributed


    subroutine initialize_struct_Distributed (this, matrix_size, matrix_type, block_size)
        class(DistributedMatrix) :: this
        integer, intent(in)      :: matrix_size, matrix_type, block_size
        integer                  :: ifail

        call this % temp_cache % clear

        this % memory_scale = 0.75_wp
        call this % setup_diag_matrix(matrix_size, matrix_type, block_size)

        this % continuum_counter = 0
        this % L2_counter = 0

        !this is an estimate based on the number of L2*contiuum functions
        !We will use a a rough estimate which wuill be the dimension of the matrix *100
        !Lets figure out when we need to perform a continuum update

        call this % update_counter

    end subroutine initialize_struct_Distributed


    subroutine insert_matelem_Distributed (this, i, j, coefficient, class, thresh)
        class(DistributedMatrix) :: this
        integer,  intent(in)     :: i, j, class
        real(wp), intent(in)     :: coefficient, thresh
        logical :: dummy
        integer :: row, column

        if (nprocs <= 1) then
            dummy = this % insert_into_diag_matrix(i, j, coefficient)
            return
        end if

        if (class /= 8 .and. class /= 2  .and. this % matrix_type /= MAT_DENSE) then
            call this % temp_cache % insert_into_cache(i, j, coefficient)
        else
            !Is it within the threshold
            if (abs(coefficient) < thresh) return

            !Does it belong to me?
            if (this % insert_into_diag_matrix(i, j, coefficient)) return

            !Otherwise we insert it into the temporary cache
            call this % temp_cache % insert_into_cache(i, j, coefficient)
        end if

    end subroutine insert_matelem_Distributed


    subroutine get_matelem_Distributed (this, idx, i, j, coeff)
        class(DistributedMatrix) :: this
        integer,  intent(in)     :: idx
        integer,  intent(out)    :: i, j
        real(wp), intent(out)    :: coeff

        i = -1
        j = -1
        coeff = 0

    end subroutine get_matelem_Distributed


    subroutine setup_diag_matrix (this, matrix_size, matrix_type, block_size)
        class(DistributedMatrix) :: this
        integer, intent(in)      :: matrix_size, matrix_type, block_size
    end subroutine setup_diag_matrix


    !> This inserts an element into the hard storage which is considered the final location before diagonalization
    !> It also checks wherther the element exists within the aloowed range and tells us if it was successfully inserted
    logical function insert_into_diag_matrix (this, row, column, coefficient)
        class(DistributedMatrix) :: this
        integer,  intent(in)     :: row, column
        real(wp), intent(in)     :: coefficient

        insert_into_diag_matrix = .false.
        call mpi_xermsg('DistributedMatrix_module', 'insert_into_diag_matrix', 'Not implemented', 1, 1)

    end function insert_into_diag_matrix


    subroutine update_continuum_Distributed (this, force_update)
        class(DistributedMatrix) :: this
        logical, intent(in)      :: force_update

        real(wp), allocatable    :: matrix_coeffs(:)
        real(wp)                 :: coeff

        integer :: number_of_chunks, num_elements, nelms_chunk, num_elems, i, j, ido, jdo, ierr, error, length, temp
        logical :: dummy

        if (this % temp_cache % is_empty()) return

        if (nprocs <= 1) return

        this % continuum_counter = this % continuum_counter + 1

        if (this % continuum_counter < this % start_continuum_update .and. .not. force_update) return

        this % continuum_counter = 0

        call master_timer % start_timer("Update Continuum")

        !Otherwise we will need to to an MPI reduce on all elements
        !Here we make the assumption that all the processes have the same elements (which they should in this case)

        number_of_chunks = this % temp_cache % num_matrix_chunks

        write (stdout, "('Number of elements to reduce = ',i8)") this % temp_cache % get_size()

        !Loop through each chunk and reduce

        do ido = 1, number_of_chunks
            !Get the number f elements in the matrix chunk
            nelms_chunk = this % temp_cache % matrix_arrays(ido) %  num_elems
            allocate(matrix_coeffs(nelms_chunk))

            !Reduce them all
            !call mpi_reduce_inplace_sum_cfp(matrix_coeffs,nelms_chunk)
            !For some strange reason MPI_IN_PLACE does not work here -_-

            call mpi_reduceall_sum_cfp(this % temp_cache % matrix_arrays(ido) % coefficient(1:nelms_chunk), &
                                       matrix_coeffs, nelms_chunk, grid % gcomm)

            this % temp_cache % matrix_arrays(ido) % coefficient(1:nelms_chunk) = matrix_coeffs(1:nelms_chunk)

            deallocate(matrix_coeffs)
        end do

        num_elems = this % temp_cache % get_size()

        !Insert what is in the temp cache
        do ido = 1, num_elems
            call this % temp_cache % get_from_cache(ido, i, j, coeff)
            dummy = this % insert_into_diag_matrix(i, j, coeff)
        end do

        !Clear it
        call this % temp_cache % clear_and_shrink
        call master_timer % stop_timer("Update Continuum")
        call this % update_counter

    end subroutine update_continuum_Distributed


    subroutine convert_temp_cache_to_array (this, matrix_ij, matrix_coeffs)

        class(DistributedMatrix)          :: this
        integer(longint),  intent(inout)  :: matrix_ij(:,:)
        real(wp),          intent(inout)  :: matrix_coeffs(:)

        integer :: ido, i, j

        write (stdout, "('Converting cache to array')")
        do ido = 1, this % temp_cache % n
            call this % temp_cache % get_from_cache(ido, i, j, matrix_coeffs(ido))
            matrix_ij(ido, 1) = i
            matrix_ij(ido, 2) = j
        end do
        write (stdout, "('done')")
        call this % temp_cache % clear_and_shrink

    end subroutine convert_temp_cache_to_array


    subroutine update_L2_Distributed (this, force_update, count_)
        class(DistributedMatrix)     :: this
        logical, intent(in)          :: force_update
        real(wp), allocatable        :: matrix_coeffs(:)
        integer, optional            :: count_

        integer(longint), allocatable, target :: matrix_ij(:,:)
        integer(longint), pointer             :: mat_ptr(:)
        integer(longint)                      :: my_num_of_elements, procs_num_of_elements, largest_num_of_elems

        integer  :: count_amount, ido, proc_id, i, j, ierr
        logical  :: dummy
        real(wp) :: coeff

        count_amount = 1

        if (present(count_)) count_amount = count_

        this % L2_counter = this % L2_counter + count_amount

        if (this % L2_counter < this % start_L2_update .and. .not. force_update) return

        my_num_of_elements = this % temp_cache % get_size()

        if (nprocs <= 1) then
            do ido = 1, my_num_of_elements
                call this % temp_cache % get_from_cache(ido, i, j, coeff)
                dummy = this % insert_into_diag_matrix(i, j, coeff)
            end do

            call this % temp_cache % clear
            return
        end if

        call mpi_reduceall_max(my_num_of_elements, largest_num_of_elems, grid % gcomm)
        !No point if we aint full
        if (largest_num_of_elems < this % start_L2_update .and. .not. force_update) then
            this % L2_counter = largest_num_of_elems
            return
        end if

        this % L2_counter = 0

        call this % temp_cache % shrink_capacity

        !Lets start off with a much simpler method and move on to something more complex
        !First off we find the largest_number of elements between procs

        write (stdout, "('The largest number of elements is ',i10,' mine is ',i10)") largest_num_of_elems, my_num_of_elements
        flush(stdout)
        !Now we allocate the needed space

        if (largest_num_of_elems == 0) return
        write (stdout, "('Starting L2 update')")
        flush(stdout)

        call master_timer % start_timer("Update L2")
        allocate(matrix_ij(largest_num_of_elems, 2), matrix_coeffs(largest_num_of_elems), stat = ierr)

        call master_memory % track_memory(storage_size(matrix_ij)/8,     size(matrix_ij),     ierr, "DIST::L2UPDATE::MAT_IJ")
        call master_memory % track_memory(storage_size(matrix_coeffs)/8, size(matrix_coeffs), ierr, "DIST::L2UPDATE::MAT_COEFFS")

        if (ierr /=0 ) then
            write (stdout, "('Memory allocation error during update!')")
            stop "Memory allocation error"
        end if

        matrix_ij = 0
        matrix_coeffs = 0.0

        mat_ptr(1:largest_num_of_elems*2) => matrix_ij(:,:)

        call this % convert_temp_cache_to_array(matrix_ij(:,:), matrix_coeffs(:))
        call this % temp_cache%clear_and_shrink

        procs_num_of_elements = my_num_of_elements

        !whats nice is that we can guarentee that the processes caches will not contain coefficients that belong to itself
        do proc_id = 1, grid % gprows * grid % gpcols - 1
            call mpi_mod_rotate_arrays_around_ring(procs_num_of_elements, mat_ptr, &
                                                   matrix_coeffs, largest_num_of_elems, grid % gcomm)
            do ido = 1, procs_num_of_elements
                dummy = this % insert_into_diag_matrix(int(matrix_ij(ido,1)), int(matrix_ij(ido,2)), matrix_coeffs(ido))
            end do
        end do

        call master_memory % free_memory(storage_size(matrix_ij)/8, size(matrix_ij))
        call master_memory % free_memory(storage_size(matrix_coeffs)/8, size(matrix_coeffs))

        !Free the space
        deallocate(matrix_ij,matrix_coeffs)

        !Reduce the temporary cache to free more space
        call this % temp_cache % clear_and_shrink
        call master_timer % stop_timer("Update L2")

        write (stdout, "('Finished L2 update')")
        flush(stdout)

        call this % update_counter

    end subroutine update_L2_Distributed


    subroutine clear_matrix (this)
        class(DistributedMatrix) :: this
    end subroutine clear_matrix


    subroutine finalize_matrix_self (this)
        class(DistributedMatrix) :: this
    end subroutine finalize_matrix_self


    subroutine finalize_Distributed (this)
        class(DistributedMatrix) :: this

        call mpi_reduceall_inplace_sum_cfp(this % diagonal, this % matrix_dimension, grid % gcomm)
        call this % finalize_matrix_self

    end subroutine finalize_Distributed


    subroutine destroy_matrix (this)
        class(DistributedMatrix) :: this
    end subroutine destroy_matrix


    subroutine print_Distributed (this)
        class(DistributedMatrix) :: this

        write (stdout, "('-------TEMP CACHE---------')")

        call this % temp_cache % print

    end subroutine print_Distributed


    subroutine clear_Distributed (this)
        class(DistributedMatrix) :: this

        call this % temp_cache % clear_and_shrink
        call this % clear_matrix

    end subroutine clear_Distributed


    subroutine destroy_Distributed (this)
        class(DistributedMatrix) :: this

        call this % clear
        call this % temp_cache % destroy
        call this % destroy_matrix

    end subroutine destroy_Distributed


    subroutine update_counter (this)
        class(DistributedMatrix) :: this
        integer  :: ifail, per_elm, dummy_int, c_update, l_update
        real(wp) :: dummy_real

        this % continuum_counter = 0
        this % L2_counter = 0

        !this is an estimate based on the number of L2*contiuum functions
        !We will use a a rough estimate which wuill be the dimension of the matrix *100
        !Lets figure out when we need to perform a continuum update

        per_elm = kind(dummy_int) * 2 + kind(dummy_real) + 4

        c_update = master_memory % get_scaled_available_memory(this % memory_scale) / (per_elm * 2)
        l_update = master_memory % get_scaled_available_memory(this % memory_scale) / (per_elm * 2)

        call mpi_reduceall_min(c_update, this % start_continuum_update, grid % gcomm)
        call mpi_reduceall_min(l_update, this % start_L2_update, grid % gcomm)

        call master_memory % print_memory_report

        write (stdout, "(2i16,' updates will occur at continuum = ',i12,' and L2 = ',i12)") &
            c_update, l_update, this % start_continuum_update, this % start_L2_update

    end subroutine update_counter

end module DistributedMatrix_module
