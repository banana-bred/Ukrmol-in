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

!> \brief   Diagonalizer type using ELPA backend
!> \authors J Benda
!> \date    2019
!>
!> This module will be built only if ELPA (and ScaLAPACK) is available.
!>
module ELPADiagonalizer_module

    use Diagonalizer_module,          only: BaseDiagonalizer
    use BaseIntegral_module,          only: BaseIntegral
    use blas_lapack_gbl,              only: blasint
    use const_gbl,                    only: stdout
    use mpi_gbl,                      only: mpi_xermsg, nprocs, myrank
    use Options_module,               only: Options
    use Parallelization_module,       only: grid => process_grid
    use precisn_gbl,                  only: wp
    use ELPAMatrix_module,            only: ELPAMatrix
    use SCALAPACKDiagonalizer_module, only: SCALAPACKDiagonalizer
    use SCALAPACKMatrix_module,       only: SCALAPACKMatrix

    use iso_c_binding, only: c_int
    use elpa,          only: elpa_t, elpa_allocate, elpa_deallocate, elpa_init, elpa_uninit, ELPA_OK, ELPA_SOLVER_2STAGE

    implicit none

    integer, parameter :: elpaint = c_int

    !> \brief   Diagonalizer class
    !> \authors J Benda
    !> \date    2019
    !>
    !> ELPA diagonalizer class is based on SCALAPACKDiagonalizer and shares with its
    !> parent class everything (including processing of solutions) except the very backend
    !> used for the symmetric matrix diagonalization (subroutine diagonalize_backend).
    !>
    type, extends(SCALAPACKDiagonalizer) :: ELPADiagonalizer
    contains
        procedure, public :: diagonalize_backend => diagonalize_backend_ELPA
    end type ELPADiagonalizer

contains

    !> \brief   Initialize the ELPA library
    !> \authors J Benda
    !> \date    2019
    !>
    !> Set upt the ELPA library for use. This needs to be called once at the beginning of the
    !> program. Initialization of the library selects a specific version of the API to use
    !> when communicating with the library. Currently, the API of version 20181112 is used.
    !>
    subroutine initialize_elpa

        integer(elpaint), parameter :: ELPA_API = 20181112
        integer(elpaint) :: success

        success = elpa_init(ELPA_API)
        call elpa_assert('initialize_elpa', 'ELPA API version not supported', success)

    end subroutine initialize_elpa


    !> \brief   Diagonalization kernel
    !> \authors J Benda
    !> \date    2019
    !>
    !> Perform the diagonalization of a symmetric cyclically block-distributed ScaLAPACK matrix.
    !> Use the ELPA library. It is assumed that the library has been initialized prior to this
    !> call. The eigenvectors will be stored in `z_mat`, eigenvalues in `w`, both of which will
    !> be allocated by this subroutine.
    !>
    subroutine diagonalize_backend_ELPA (this, matrix_elements, num_eigenpair, z_mat, w, all_procs, option, integrals)

        class(ELPADiagonalizer)               :: this
        class(SCALAPACKMatrix), intent(in)    :: matrix_elements
        class(BaseIntegral),    intent(in)    :: integrals
        type(Options),          intent(in)    :: option
        integer,                intent(in)    :: num_eigenpair
        logical,                intent(in)    :: all_procs
        real(wp), allocatable,  intent(inout) :: z_mat(:,:), w(:)

        integer                :: ierr
        integer(elpaint)       :: success, nb, loc_r, loc_c, mat_dimen, myrow, mycol, nev, comm
        class(elpa_t), pointer :: elpa

        write (stdout, "('Diagonalization will be done with ELPA')")

        comm  = grid % gcomm
        myrow = grid % mygrow
        mycol = grid % mygcol
        nev   = num_eigenpair

        mat_dimen = matrix_elements % get_matrix_size()
        loc_r     = matrix_elements % local_row_dimen
        loc_c     = matrix_elements % local_col_dimen
        nb        = matrix_elements % scal_block_size

        ! allocate the eigenvectors
        allocate(z_mat(loc_r, loc_c), w(mat_dimen), stat = ierr)
        if (ierr /= 0) then
            call mpi_xermsg('ELPADiagonalizer_module', 'diagonalize_backend_ELPA', 'Memory allocation error.', ierr, 1)
        end if

        ! create the ELPA object
        elpa => elpa_allocate(success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to create ELPA object.', success)

        ! set matrix-related ELPA parameters
        call elpa % set('na', mat_dimen, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "na" for ELPA.', success)
        call elpa % set('nev', mat_dimen, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "nev" for ELPA.', success)
        call elpa % set('local_nrows', loc_r, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "local_nrows" for ELPA.', success)
        call elpa % set('local_ncols', loc_c, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "local_ncols" for ELPA.', success)
        call elpa % set('nblk', nb, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "nblk" for ELPA.', success)
        call elpa % set('mpi_comm_parent', comm, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "mpi_comm_parent" for ELPA.', success)
        call elpa % set('process_row', myrow, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "process_row" for ELPA.', success)
        call elpa % set('process_col', mycol, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set "process_col" for ELPA.', success)

        ! finalize the setup
        success = elpa % setup()
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to set up the ELPA object.', success)

        ! use 2-stage solver by default (i.e. when not given in environment variable)
        call get_environment_variable('ELPA_DEFAULT_solver', status = ierr)
        if (ierr > 0) then
            ! environment variable ELPA_DEFAULT_solver does not exist (= 1) or environment variables are not supported at all (= 2)
            call elpa % set("solver", ELPA_SOLVER_2STAGE, success)
            call elpa_assert('diagonalize_backend_ELPA', 'Failed to switch ELPA to the 2-stage solver.', success)
        end if

        ! run the solver
        call elpa % eigenvectors(matrix_elements % a_local_matrix, w, z_mat, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Diagonalization using ELPA failed.', success)

        ! free the memory, possible getting ready for another diagonalization
        call elpa_deallocate(elpa, success)
        call elpa_assert('diagonalize_backend_ELPA', 'Failed to deallocate the ELPA object.', success)

    end subroutine diagonalize_backend_ELPA


    !> \brief   Check ELPA error code
    !> \authors J Benda
    !> \date    2019
    !>
    !> Verify that the given ELPA error code is equal to ELPA_OK. If not, print the given error message
    !> and terminate the program.
    !>
    subroutine elpa_assert (sub, msg, err)

        character(len=*), intent(in) :: sub, msg
        integer(elpaint), intent(in) :: err

        if (err /= ELPA_OK) then
            call mpi_xermsg('ELPADiagonalizer_module', sub, msg, int(err), 1)
        end if

    end subroutine elpa_assert

end module ELPADiagonalizer_module
