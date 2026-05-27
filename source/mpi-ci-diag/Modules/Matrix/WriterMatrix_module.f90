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

!> \brief   Write matrix module
!> \authors A Al-Refaie
!> \date    2017
!>
!> Defines the WriterMatrix type which is responsible for writing the Hamiltonian matrix
!> to a disk file. This matrix type is used by all serial diagonalizers.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module WriterMatrix_module

    use const_gbl,                only: stdout
    use mpi_gbl,                  only: master
    use precisn,                  only: wp
    use scatci_routines,          only: WRTEM
    use DistributedMatrix_module, only: DistributedMatrix
    use MatrixCache_module,       only: MatrixCache
    use Options_module,           only: Options
    use Parallelization_module,   only: grid => process_grid

    implicit none

    public WriterMatrix

    private

    !> \brief   Matrix associated with a disk drive
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This type is used for construction of the final Hamiltonian matrix. It allows continuous
    !> insertion of matrix elements, which are periodically flushed (as individual "records")
    !> to the matrix disk file unit in order to maintain reasonable consumption of memory.
    !>
    type, extends(DistributedMatrix) :: WriterMatrix
        !> This cache holds the temprorary matrix values before they are sent to the relevant process.
        !> This cache holds what each process requires for SCALAPACK diagonalization.
        type(MatrixCache) :: write_cache

        integer :: num_elements_per_record
        integer :: matrix_unit
    contains
        procedure, public :: print                   => print_Writer
        procedure, public :: finalize_matrix_self    => finalize_matrix_Writer
        procedure, public :: set_options             =>  set_options_writer
        procedure, public :: get_matelem_self        => get_matelem_Writer
        procedure, public :: clear_matrix            => clear_Writer
        procedure, public :: destroy_matrix          => destroy_Writer
        procedure, public :: write_cache_to_unit
        procedure, public :: insert_into_diag_matrix => insert_into_write_cache
        procedure, public :: get_matrix_unit
    end type WriterMatrix

contains

    !> \brief   Return Hamiltonian disk file unit
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer function get_matrix_unit (this)
        class(WriterMatrix) :: this

        get_matrix_unit = this % matrix_unit

    end function get_matrix_unit


    !> \brief   Initialize the matrix cache using the Options object
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine set_options_writer (this, option_val)
        class(WriterMatrix)        :: this
        class(Options), intent(in) :: option_val

        this % matrix_unit = option_val % hamiltonian_unit
        this % num_elements_per_record = option_val % num_matrix_elements_per_rec

        !we set the expansion size to the number of records + 1
        call this % write_cache % construct(this % num_elements_per_record + 1)

    end subroutine set_options_writer


    !> \brief   Not implemented
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine get_matelem_Writer (this, idx, i, j, coeff)
        class(WriterMatrix)   :: this
        integer,  intent(in)  :: idx
        integer,  intent(out) :: i, j
        real(wp), intent(out) :: coeff

        i = -1
        j = -1
        coeff = -1.0

        stop "get element not yet implemented"

    end subroutine get_matelem_Writer


    !> \brief   Write matrix chunk to file
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Unless the element cache is empty, this subroutine flushes the elements in cache to the Hamiltonian
    !> disk file unit as a new record. The matrix element cache is empty on return from this subroutine.
    !>
    subroutine write_cache_to_unit (this)
        class(WriterMatrix) :: this
        integer             :: record_count

        if (this % write_cache % is_empty()) return

        record_count = this % write_cache % get_size()

        !Write it all out
        call WRTEM(this % matrix_unit, record_count, this % num_elements_per_record, &
                   this % write_cache % matrix_arrays(1) % ij(1:2,1:record_count), &
                   this % write_cache % matrix_arrays(1) % coefficient(1:record_count))
        call this % write_cache % clear

    end subroutine write_cache_to_unit


    !> \brief   Insert matrix element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This inserts an element into the hard storage which is considered the final location before diagonalization.
    !> It also checks wherther the element exists within the allowed range and tells us if it was successfully inserted.
    !> Elements smaller than a specific threshold will be ignored. When the limit size of the cache is reached,
    !> it is emptied (appended) to the Hamiltonian disk file.
    !>
    logical function insert_into_write_cache (this, row, column, coefficient)
        class(WriterMatrix)  :: this
        integer,  intent(in) :: row, column
        real(wp), intent(in) :: coefficient
        integer              :: record_count

        if (grid % grank /= master) then
            insert_into_write_cache = .false.
            return
        end if

        if (row == column) then
            if (row > size(this % diagonal)) then
                write (stdout, "('Row greater than size of diagonal!!!!',2i12)") row, size(this % diagonal)
                stop "Error"
            end if
            this % diagonal(row) = coefficient
        end if

        insert_into_write_cache = .true.
        if (abs(coefficient) < this % threshold) return

        call this % write_cache % insert_into_cache(row, column, coefficient)
        record_count = this % write_cache % get_size()

        !Update number of elements
        this % n = this % n + 1
        insert_into_write_cache = .true.

        if (record_count == this % num_elements_per_record) then
            !Write it all out
            call this % write_cache_to_unit
        end if

    end function insert_into_write_cache


    !> \brief   Write Hamiltonian to disk
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Uses SCATCI subroutine WRTEM to write Hamiltonian matrix to its disk file unit.
    !> Also, appends a dummy finalizing record at the end of the file.
    !>
    subroutine finalize_matrix_Writer (this)
        class(WriterMatrix) :: this
        integer  :: dummy_ij(2, this % num_elements_per_record), record_count = 0
        real(wp) :: dummy_coeff(this % num_elements_per_record)

        if (grid % grank /= master) return

        write (stdout, "('finishing unit')")
        call this % write_cache_to_unit
        record_count = 0
        call WRTEM(this % matrix_unit, record_count, this % num_elements_per_record, dummy_ij, dummy_coeff)
        write (stdout, "('done')")

    end subroutine finalize_matrix_Writer


    !> \brief   Print matrix to standard output
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine print_Writer (this)
        class(WriterMatrix) :: this

        write (stdout, "('-------TEMP CACHE---------')")
        call this % temp_cache % print
        !write(stdout,"('-------HARD CACHE---------')")
        !call this%matrix_cache%print

    end subroutine print_Writer


    !> \brief   Clear matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Clears matrix element cache.
    !>
    subroutine clear_Writer (this)
        class(WriterMatrix) :: this

        call this % write_cache % clear

    end subroutine clear_Writer


    !> \brief   Destroy matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine destroy_Writer (this)
        class(WriterMatrix) :: this

        call this % clear
        call this % write_cache % destroy

    end subroutine destroy_Writer

end module WriterMatrix_module
