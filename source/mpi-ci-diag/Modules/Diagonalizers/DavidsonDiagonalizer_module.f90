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

!> \brief   Diagonalizer type using Davidson backend
!> \authors A Al-Refaie
!> \date    2017
!>
!> This type is always available. Uses the SCATCI routine \c mkdvm.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module DavidsonDiagonalizer_module

    use precisn,                   only: wp
    use const_gbl,                 only: stdout
    use BaseIntegral_module,       only: BaseIntegral
    use BaseMatrix_module,         only: BaseMatrix
    use WriterMatrix_module,       only: WriterMatrix
    use Diagonalizer_module,       only: BaseDiagonalizer
    use DiagonalizerResult_module, only: DiagonalizerResult
    use Options_module,            only: Options
    use mpi_gbl,                   only: master, myrank, mpi_mod_bcast, mpi_xermsg
    use consts_mpi_ci,             only: PASS_TO_CDENPROP

    implicit none

    type, extends(BaseDiagonalizer) :: DavidsonDiagonalizer
    contains
        procedure, public  :: diagonalize => diagonalize_davidson
        procedure, private :: diagonalize_writermatrix
    end type DavidsonDiagonalizer

contains

    subroutine diagonalize_davidson(this, matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)
        class(DavidsonDiagonalizer)     :: this
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
            max_iterations = max(num_eigenpair * 40, 500)
        end if

        print *, 'max_iterations', max_iterations
        print *, 'max_tolerance', max_tolerance

        select type(matrix_elements)
            type is (WriterMatrix)
                call this % diagonalize_writermatrix(matrix_elements, num_eigenpair, dresult, all_procs, &
                                                     max_iterations, max_tolerance, option, integrals)
            class is (BaseMatrix)
                call mpi_xermsg('DavidsonDiagonalizer_module', 'diagonalize_davidson', &
                                'Only WriterMatrix format can use Davidson', 1, 1)
        end select

    end subroutine diagonalize_davidson


    subroutine diagonalize_writermatrix (this, matrix_elements, num_eigenpair, dresult, all_procs, max_iterations, max_tolerance, &
                                         option, integrals)
        class(DavidsonDiagonalizer)     :: this
        class(DiagonalizerResult)       :: dresult
        type(WriterMatrix),  intent(in) :: matrix_elements
        class(BaseIntegral), intent(in) :: integrals
        type(Options),       intent(in) :: option
        integer,             intent(in) :: num_eigenpair
        logical,             intent(in) :: all_procs
        integer,             intent(in) :: max_iterations
        real(wp),            intent(in) :: max_tolerance

        !for now these are default but will be included i the input

        real(wp) :: CRITC = 1E-10_wp
        real(wp) :: CRITR = 1E-8_wp
        real(wp) :: ORTHO = 1E-7_wp

        integer :: matrix_unit
        integer :: matrix_size
        integer :: num_matrix_elements_per_record, num_elems

        real(wp), allocatable :: hamiltonian(:), diagonal(:), eig(:)
        real(wp), pointer     :: vec_ptr(:)

        integer :: mpi_num_eig, mpi_num_vecs, ido, ierr

        real(wp), allocatable :: eigen(:)
        real(wp), allocatable :: vecs(:,:)

        character(len=4), dimension(30) :: NAMP
        integer,          dimension(10) :: NHD
        integer,          dimension(20) :: KEYCSF, NHE
        real(kind=wp),    dimension(41) :: DTNUC

        write (stdout, "('Diagonalization done with Davidson')")
        write (stdout, "('Parameters:')")
        write (stdout, "('N: ',i8)") matrix_elements % get_matrix_size()
        write (stdout, "('Requested # of eigenpairs',i8)") num_eigenpair

        allocate(eigen(num_eigenpair), vecs(matrix_elements % get_matrix_size(), num_eigenpair))

        eigen = 0
        vecs = 0
        matrix_size = 0

        if (myrank == master) then
            !Get our matrix unit
            matrix_unit = matrix_elements % get_matrix_unit()

            !Lets follow scatcis example even though we can get these information from the matrix class itself
            rewind matrix_unit
            read (matrix_unit) matrix_size, num_matrix_elements_per_record, NHD, NAMP, NHE, DTNUC

            num_elems = matrix_elements % get_size()
            print *, 'matrix_size', matrix_size
            print *, 'num_elems', num_elems

            ! allocate(hamiltonian(matrix_size*num_eigenpair),diagonal(matrix_size),eig(matrix_size))
            allocate(diagonal(matrix_size))
            call mkdvm(vecs, eigen, diagonal, matrix_size, num_eigenpair, stdout, matrix_unit, &
                       num_matrix_elements_per_record, num_elems, max_tolerance, CRITC, CRITR, ORTHO, &
                       max_iterations, ierr)

            if (ierr /= 0) then
                call mpi_xermsg('DavidsonDiagonalizer_module', 'DavidsonDiagonalizer', &
                                'Davidson diagonalization failed', ierr, 1)
            end if
        end if

        if (all_procs) then
            call mpi_mod_bcast(matrix_size, master)
            mpi_num_vecs = matrix_size * num_eigenpair
            mpi_num_eig = num_eigenpair
            call mpi_mod_bcast(eigen,master)
            do ido = 1, num_eigenpair
                call mpi_mod_bcast(vecs(:,ido), master)
            end do
        end if

        if (iand(option % vector_storage_method, PASS_TO_CDENPROP) /= 0) then
            call dresult % export_header(option, integrals)
        end if
        call dresult % write_header(option, integrals)
        call dresult % handle_eigenvalues(eigen, matrix_elements % diagonal, num_eigenpair, matrix_size)

        do ido = 1, num_eigenpair
            call dresult % handle_eigenvector(vecs(:,ido), matrix_size)
        end do

        ! store eigensolutions for the post-processing stage
        if (iand(dresult % vector_storage, PASS_TO_CDENPROP) /= 0) then
            !Eigenvectors:
            if (allocated(dresult % ci_vec % CV)) deallocate(dresult % ci_vec % CV)
            call move_alloc(vecs, dresult % ci_vec % CV)
            dresult % ci_vec % CV_is_scalapack = .false.

            !Eigenvalues and phases:
            call dresult % export_eigenvalues(eigen(1:num_eigenpair), matrix_elements % diagonal, &
                                              num_eigenpair, matrix_size, dresult % ci_vec % ei, dresult % ci_vec % iphz)

            write (stdout, '(/," Eigensolutions exported into the CIVect structure.",/)')
        else
            deallocate (vecs)
        end if

        deallocate(eigen)

    end subroutine diagonalize_writermatrix

end module DavidsonDiagonalizer_module
