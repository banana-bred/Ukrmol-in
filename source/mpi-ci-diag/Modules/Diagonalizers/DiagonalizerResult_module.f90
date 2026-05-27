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

!> \brief   Diagonalizer result data object
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module DiagonalizerResult_module

    use cdenprop_defs,       only: CIVect
    use consts_mpi_ci,       only: SAVE_CONTINUUM_COEFFS
    use BaseIntegral_module, only: BaseIntegral
    use Options_module,      only: Options

    implicit none

    public DiagonalizerResult

    private

    !> \brief   Output from diagonalization
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Base class for \ref SolutionHandler_module::SolutionHandler and \ref Target_RMat_CI_module::Target_RMat_CI.
    !> Can be used as a storage for the eigensolution when that is to be passed to CDENPROP.
    !>
    type, abstract :: DiagonalizerResult
        type(CIVect) :: ci_vec  !< Data structure for keeping the solution of the eigensystem.

        integer :: vector_storage = SAVE_CONTINUUM_COEFFS
        integer :: continuum_dimen
    contains
        procedure(handle_eigenvalues), deferred, public :: handle_eigenvalues
        procedure(handle_eigenvector), deferred, public :: handle_eigenvector
        procedure(export_eigenvalues), deferred, public :: export_eigenvalues
        procedure :: export_header
        procedure :: write_header
        procedure :: finalize_solutions
    end type DiagonalizerResult

    abstract interface
        !> \brief   Main build routine of the hamiltonian
        !> \authors A Al-Refaie
        !> \date    2017
        !>
        !> All build must be done within this routine in order to be used by MPI-SCATCI.
        !>
        !> \param[out] matrix_elements  Resulting matrix elements from the build.
        !>
        subroutine handle_eigenvalues (this, eigenvalues, diagonals, num_eigenpairs, vec_dimen)
            use precisn, only : wp
            import :: DiagonalizerResult
            class(DiagonalizerResult) :: this
            integer,  intent(in)      :: num_eigenpairs, vec_dimen
            real(wp), intent(in)      :: eigenvalues(num_eigenpairs), diagonals(vec_dimen)
        end subroutine handle_eigenvalues

        !> \brief   Main build routine of the hamiltonian
        !> \authors A Al-Refaie
        !> \date    2017
        !>
        !> All build must be done within this routine in order to be used by MPI-SCATCI.
        !>
        !> \param[out] matrix_elements  Resulting matrix elements from the build.
        !>
        subroutine handle_eigenvector (this, eigenvector, vec_dimen)
            use precisn, only : wp
            import :: DiagonalizerResult
            class(DiagonalizerResult) :: this
            integer,  intent(in)      :: vec_dimen
            real(wp), intent(in)      :: eigenvector(vec_dimen)
        end subroutine handle_eigenvector

        !> \brief   build routine of the hamiltonian
        !> \authors A Al-Refaie
        !> \date    2017
        !>
        !> All build must be done within this routine in order to be used by MPI-SCATCI.
        !>
        !> \param[out] matrix_elements  Resulting matrix elements from the build.
        !>
        subroutine export_eigenvalues (this, eigenvalues, diagonals, num_eigenpairs, vec_dimen, ei, iphz)
            use precisn, only : wp
            import :: DiagonalizerResult
            class(DiagonalizerResult) :: this
            integer,  intent(in)      :: num_eigenpairs, vec_dimen
            real(wp), intent(in)      :: eigenvalues(num_eigenpairs), diagonals(vec_dimen)
            real(wp), allocatable     :: ei(:)
            integer,  allocatable     :: iphz(:)
        end subroutine export_eigenvalues
    end interface

contains

    !> \brief   Save eigenvector information internally
    !> \authors J Benda
    !> \date    2019
    !>
    !> This dummy function is expected to be overridden by derived types that support passing
    !> the eigensolution to CDENPROP. This is used in \ref SolutionHandler_module::SolutionHandler
    !> after the main diagonalization.
    !>
    subroutine export_header (this, option, integrals)
        class(DiagonalizerResult)       :: this
        class(Options),      intent(in) :: option
        class(BaseIntegral), intent(in) :: integrals
    end subroutine export_header


    !> \brief   Save eigenvector information to disk
    !> \authors J Benda
    !> \date    2019
    !>
    !> This dummy function is expected to be overridden by derived types that support writing
    !> the eigensolution to disk file. This is used in \ref SolutionHandler_module::SolutionHandler
    !> after the main diagonalization.
    !>
    subroutine write_header (this, option, integrals)
        class(DiagonalizerResult)       :: this
        class(Options),      intent(in) :: option
        class(BaseIntegral), intent(in) :: integrals
    end subroutine write_header


    !> \brief   Finalize saving eigenvector information to disk
    !> \authors J Benda
    !> \date    2019
    !>
    !> This dummy function is expected to be overridden by derived types that support writing
    !> the eigensolution to disk file. This is used in \ref SolutionHandler_module::SolutionHandler
    !> after the main diagonalization.
    !>
    subroutine finalize_solutions (this, option)
        class(DiagonalizerResult)       :: this
        class(Options),      intent(in) :: option
    end subroutine finalize_solutions

end module DiagonalizerResult_module
