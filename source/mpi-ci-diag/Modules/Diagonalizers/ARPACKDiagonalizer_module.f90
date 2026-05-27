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

!> \brief   Diagonalizer type using Arpack backend
!> \authors A Al-Refaie
!> \date    2017
!>
!> Requires Arpack routine \c mkarp. This module is only included in the build if ARPACK_LIBRARIES are
!> given on the CMake command line.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module ARPACKDiagonalizer_module

    use precisn,                   only: longint, wp
    use const_gbl,                 only: stdout
    use BaseMatrix_module,         only: BaseMatrix
    use BaseIntegral_module,       only: BaseIntegral
    use WriterMatrix_module,       only: WriterMatrix
    use Diagonalizer_module,       only: BaseDiagonalizer
    use DiagonalizerResult_module, only: DiagonalizerResult
    use Options_module,            only: Options
    use mpi_gbl,                   only: myrank, master, mpi_xermsg
    use m_4_mkarp,                 only: mkarp, default_arpack_max_iter
    use consts_mpi_ci,             only: PASS_TO_CDENPROP

    implicit none

    type, extends(BaseDiagonalizer) :: ARPACKDiagonalizer
    contains
        procedure, public  :: diagonalize => diagonalize_ARPACK
        procedure, private :: diagonalize_writermatrix
    end type ARPACKDiagonalizer

contains

    subroutine diagonalize_ARPACK (this, matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)
        class(ARPACKDiagonalizer)       :: this
        class(DiagonalizerResult)       :: dresult
        class(BaseMatrix),   intent(in) :: matrix_elements
        class(BaseIntegral), intent(in) :: integrals
        type(Options),       intent(in) :: option
        integer,             intent(in) :: num_eigenpair
        logical,             intent(in) :: all_procs

        integer  :: max_iterations
        real(wp) :: max_tolerance

        !Only the first rank can do it
        if (myrank /= master) return

        max_iterations = option % max_iterations
        max_tolerance  = option % max_tolerance

        if (max_iterations < 0) then
            max_iterations = MAX(num_eigenpair * 50, default_arpack_max_iter)
        end if

        select type(matrix_elements)
            type is (WriterMatrix)
                call this % diagonalize_writermatrix(matrix_elements, num_eigenpair, dresult, max_iterations, max_tolerance, &
                                                     option, integrals)
            class is (BaseMatrix)
                call mpi_xermsg('ARPACKDiagonalizer_module', 'diagonalize_ARPACK', &
                                'Only WriterMatrix format can use Arpack', 1, 1)
        end select

    end subroutine diagonalize_ARPACK


    subroutine diagonalize_writermatrix (this, matrix_elements, num_eigenpair, dresult, max_iterations, max_tolerance, &
                                         option, integrals)
        class(ARPACKDiagonalizer)       :: this
        class(DiagonalizerResult)       :: dresult
        type(WriterMatrix),  intent(in) :: matrix_elements
        class(BaseIntegral), intent(in) :: integrals
        type(Options),       intent(in) :: option
        integer,             intent(in) :: num_eigenpair
        integer,             intent(in) :: max_iterations
        real(wp),            intent(in) :: max_tolerance

        !for now these are default but will be included i the input

        integer :: matrix_unit
        integer :: matrix_size
        integer :: i1, i2

        integer  :: arpack_numeig
        integer  :: arpack_maxit
        real(wp) :: arpack_maxtol
        integer  :: num_matrix_elements_per_record, num_elems

        real(wp), allocatable :: eigenvalues(:), eigenvector(:)

        write (stdout, "('Diagonalization done with ARPACK')")
        write (stdout, "('Parameters:')")
        write (stdout, "('N: ',i8)") matrix_elements % get_matrix_size()
        write (stdout, "('Requested # of eigenpairs',i8)") num_eigenpair

        arpack_numeig = num_eigenpair
        arpack_maxit = max_iterations
        arpack_maxtol = max_tolerance
        matrix_size = matrix_elements % get_matrix_size()
        num_elems = matrix_elements % get_size()
        matrix_unit = matrix_elements % get_matrix_unit()

        allocate(eigenvalues(num_eigenpair), eigenvector(matrix_size))

        call mkarp(matrix_size, arpack_numeig, arpack_maxit, arpack_maxtol, matrix_unit)

        ! 3. Opening the produced file "eigensystem_arpack_file.bin" file.
        !    This file contains the produced eigenpairs and it is
        !    written to the hard disk.

        open (unit = 11, file = "___tmp_eigensys", status = "old", action = "read", form = "unformatted")

        ! 4. Reading the eigenvalues.

        do i1 = 1, num_eigenpair
            read (unit = 11) eigenvalues(i1)
        end do

        if (iand(option % vector_storage_method, PASS_TO_CDENPROP) /= 0) then
            call dresult % export_header(option, integrals)
            if (allocated(dresult % ci_vec % CV)) deallocate(dresult % ci_vec % CV)
            allocate (dresult % ci_vec % CV(matrix_size, num_eigenpair))
            dresult % ci_vec % CV_is_scalapack = .false.
        end if
        call dresult % write_header(option, integrals)
        call dresult % handle_eigenvalues(eigenvalues, matrix_elements % diagonal, num_eigenpair, matrix_size)

        ! 5. Reading the eigenvectors.

        do i1 = 1, num_eigenpair
            do i2 = 1, matrix_size
                read (unit = 11) eigenvector(i2)
            end do
            call dresult % handle_eigenvector(eigenvector, matrix_size)
            if (iand(option % vector_storage_method, PASS_TO_CDENPROP) /= 0) then
                dresult % ci_vec % CV(1:matrix_size, i1) = eigenvector(1:matrix_size)
            end if
            eigenvector = 0
        end do

        ! 6. Closing the unit and deleting the file
        !    "___tmp_eigensys".

        close (unit = 11, status = "delete")

        if (iand(option % vector_storage_method, PASS_TO_CDENPROP) /= 0) then
            call dresult % export_eigenvalues(eigenvalues(1:num_eigenpair), matrix_elements % diagonal, &
                                              num_eigenpair, matrix_size, dresult % ci_vec % ei, dresult % ci_vec % iphz)
            write (stdout, '(/," Eigensolutions exported into the CIVect structure.",/)')
        end if

        deallocate(eigenvalues, eigenvector)

    end subroutine diagonalize_writermatrix

end module ARPACKDiagonalizer_module
