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

!> \brief   ScaLAPACK matrix module
!> \authors A Al-Refaie, J Benda
!> \date    2017 - 2019
!>
!> Only compiled in when ScaLAPACK library is available during configuration step. Uses a compile-time macro for the BLAS/LAPACK
!> integer, which is long integer by default, or controllable by the compiler flags -Dblas64int/-Dblas32int.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!> \note 24/02/2019 - Jakub Benda: Enable multiple concurrent diagonalizations.
!>
module SCALAPACKMatrix_module

    use blas_lapack_gbl,          only: blasint
    use cdenprop_defs,            only: get_scalapack_max_block_size!, descinit
    use const_gbl,                only: stdout
    use mpi_gbl,                  only: mpi_xermsg
    use precisn,                  only: longint, wp
    use DistributedMatrix_module, only: DistributedMatrix
    use MemoryManager_module,     only: master_memory
    use Parallelization_module,   only: grid => process_grid

    implicit none

    public SCALAPACKMatrix

    private

    type, extends(DistributedMatrix) :: SCALAPACKMatrix
        !>This cache holds the temprorary matrix values before they are sent to the relevant process
        !>This cache holds what each process requires for SCALAPACK diagobnalization
        real(wp), allocatable :: a_local_matrix(:,:)
        integer(blasint) :: mat_dimen
        integer(blasint) :: scal_block_size
        integer(blasint) :: local_row_dimen, local_col_dimen
        integer(blasint) :: descr_a_mat(50), descr_z_mat(50)
        integer(blasint) :: lda
        logical  :: store_mat = .false.
    contains
        procedure, public  :: print                   => print_SCALAPACK
        procedure, public  :: am_i_involved
        procedure, public  :: setup_diag_matrix       => initialize_struct_SCALAPACK
        procedure, public  :: get_matelem_self        => get_matelem_SCALAPACK
        procedure, public  :: clear_matrix            => clear_SCALAPACK
        procedure, public  :: destroy_matrix          => destroy_SCALAPACK
        procedure, public  :: insert_into_diag_matrix => insert_into_local_matrix
        procedure, private :: is_this_me
    end type SCALAPACKMatrix

contains

    !> \brief   Is this process part of the BLACS grid
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This function returns "true" when part of the matrix is allocated to this process.
    !> It returns false if this process is not part of the BLACS grid at all.
    !>
    logical function am_i_involved (this)

        class(SCALAPACKMatrix) :: this

        am_i_involved = this % store_mat

    end function am_i_involved


    !> \brief   Initialize the type
    !> \authors A Al-Refaie, J Benda
    !> \date    2017 - 2022
    !>
    !> This subroutine also prepares the context for BLACS (ie. the MPI communicator layout).
    !> The number of processes is split into a grid of `nprow` times `npcol` processes, where
    !> both `nprow` and `npcol` are as close as possible to the square root of the number of
    !> processes.
    !>
    subroutine initialize_struct_SCALAPACK (this, matrix_size, matrix_type, block_size)

        class(SCALAPACKMatrix) :: this
        integer, intent(in)    :: matrix_size, matrix_type, block_size

        integer  :: ifail, per_elm, dummy_int, c_update, l_update
        real(wp) :: dummy_real

        !Blas variables
        integer(blasint)           :: ierr, ido, nb, info, nproc, nprow, npcol, myrow, mycol
        integer(blasint)           :: scalapack_min_block_size, scalapack_max_block_size
        integer(blasint), external :: numroc

        this % mat_dimen = matrix_size

        nprow = grid % gprows
        npcol = grid % gpcols
        myrow = grid % mygrow
        mycol = grid % mygcol

        nproc = nprow * npcol

        write (stdout, "('nproc = ',i4,'matdimen = ',i4,'nrow = ',i4,' ncol = ',i4,'myrow = ',i4,' mycol = ',i4)") &
            nproc, this % mat_dimen, nprow, npcol, myrow, mycol

        write (stdout, "('context = ',i4,'nprocs = ',i4,'matdimen = ',i8,'nrow = ',i8,' ncol = ',i8,'myrow = ',i8,' mycol = ',i8)")&
            grid % gcntxt, nproc, this % mat_dimen, nprow, npcol, myrow, mycol

        this % store_mat = .true.

        if (.not. this % am_i_involved()) return

        scalapack_min_block_size = 1
        scalapack_max_block_size = get_scalapack_max_block_size()

        this % scal_block_size = min(int(this % mat_dimen) / nprow, this % matrix_dimension / npcol)
        this % scal_block_size = min(this % scal_block_size, scalapack_max_block_size)
        this % scal_block_size = max(this % scal_block_size, scalapack_min_block_size)

        this % local_row_dimen = numroc(this % mat_dimen, this % scal_block_size, myrow, 0, nprow)
        this % local_col_dimen = numroc(this % mat_dimen, this % scal_block_size, mycol, 0, npcol)

        this % lda = max(1_blasint, this % local_row_dimen)

        write (stdout, "('block_size = ',i4,' local_row_size = ',i8,' local_col_size = ',i8,' lda = ',i8)") &
            this % scal_block_size, this % local_row_dimen, this % local_col_dimen, this % lda

        call descinit(this % descr_a_mat, this % mat_dimen, this % mat_dimen, this % scal_block_size, this % scal_block_size, &
                      0_blasint, 0_blasint, grid % gcntxt, this % lda, info)
        if (info /= 0) then
            call mpi_xermsg('SCALAPACKMatrix_module', 'initialize_struct_SCALAPACK', &
                            'Error in getting ScaLAPACK description for A', int(info), 1)
        end if

        call descinit(this % descr_z_mat, this % mat_dimen, this % mat_dimen, this % scal_block_size, this % scal_block_size, &
                      0_blasint, 0_blasint, grid % gcntxt, this % lda, info)
        if (info /= 0) then
            call mpi_xermsg('SCALAPACKMatrix_module', 'initialize_struct_SCALAPACK', &
                            'Error in getting ScaLAPACK description for Z', int(info), 1)
        end if

        if (allocated(this % a_local_matrix)) then
            call master_memory % free_memory(storage_size(this % a_local_matrix)/8, size(this % a_local_matrix))
            deallocate(this % a_local_matrix)
        end if

        allocate(this % a_local_matrix(this % local_row_dimen, this % local_col_dimen), stat = ifail)
        call master_memory % track_memory(storage_size(this % a_local_matrix)/8, size(this % a_local_matrix), ifail, &
                                          'SCALAPACKMATRIX::A_LOCAL_MATRIX')

        this % a_local_matrix = 0
        this % n = real(this % local_row_dimen) * real(this % local_col_dimen)

    end subroutine initialize_struct_SCALAPACK


    !> \brief   Check workitem association with current process
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    logical function is_this_me (this, proc_row, proc_col)

        class(SCALAPACKMatrix)       :: this
        integer(blasint), intent(in) :: proc_row, proc_col

        is_this_me = (proc_row == grid % mygrow) .and. &
                     (proc_col == grid % mygcol)

    end function is_this_me


    !> \brief   Retrieve matrix element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine get_matelem_SCALAPACK (this, idx, i, j, coeff)

        class(SCALAPACKMatrix) :: this
        integer,  intent(in)   :: idx
        integer,  intent(out)  :: i, j
        real(wp), intent(out)  :: coeff

        integer(blasint)           :: i_loc, j_loc, proc_row, proc_col
        integer(blasint), external :: indxl2g

        i_loc = idx / this % local_row_dimen  + 1
        j_loc = mod(idx, int(this % local_col_dimen)) + 1

        i = indxl2g(i_loc, this % scal_block_size, grid % mygrow, 0, grid % gprows)
        j = indxl2g(j_loc, this % scal_block_size, grid % mygcol, 0, grid % gpcols)

        coeff = this % a_local_matrix(i_loc, j_loc)

    end subroutine get_matelem_SCALAPACK


    !> \brief   Insert new element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This inserts an element into the hard storage which is considered the final location before diagonalization.
    !> It also checks wherther the element exists within the aloowed range and tells us if it was successfully inserted.
    !>
    logical function insert_into_local_matrix (this, row, column, coefficient)

        class(SCALAPACKMatrix) :: this
        integer,  intent(in)   :: row, column
        real(wp), intent(in)   :: coefficient

        integer(blasint) :: proc_row, proc_col, i_loc, j_loc, blas_row, blas_col

        blas_row = row
        blas_col = column

        if (row == column) call this % store_diagonal(row, coefficient)
        if (.not. this % am_i_involved()) return

        !Figure out which proc it belongs to and the local matrix index
        call infog2l(blas_row, blas_col, this % descr_a_mat,     &
                     grid % gprows, grid % gpcols, &
                     grid % mygrow, grid % mygcol, &
                     i_loc, j_loc, proc_row, proc_col)

        if (this % is_this_me(proc_row, proc_col)) then
            this % a_local_matrix(i_loc, j_loc) = coefficient
            insert_into_local_matrix = .true.
        else
            insert_into_local_matrix = .false.
        end if

    end function insert_into_local_matrix


    !> \brief   Print matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine print_SCALAPACK(this)
        class(SCALAPACKMatrix) :: this

        write (stdout, "('-------TEMP CACHE---------')")
        call this % temp_cache % print
        !write(stdout,"('-------HARD CACHE---------')")
        !call this%matrix_cache%print

    end subroutine print_SCALAPACK


    !> \brief   Clear matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine clear_SCALAPACK (this)
        class(SCALAPACKMatrix) :: this

        if (allocated(this % a_local_matrix)) this % a_local_matrix = 0

    end subroutine clear_SCALAPACK


    !> \brief   Destroy matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine destroy_SCALAPACK(this)
        class(SCALAPACKMatrix) :: this

        if (allocated(this % a_local_matrix)) then
            call master_memory % free_memory(storage_size(this % a_local_matrix)/8, size(this % a_local_matrix))
            if (allocated(this % a_local_matrix)) deallocate(this % a_local_matrix)
        end if

    end subroutine destroy_SCALAPACK

end module SCALAPACKMatrix_module
