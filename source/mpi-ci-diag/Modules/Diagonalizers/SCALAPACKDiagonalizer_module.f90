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

!> \brief   Diagonalizer type using ScaLAPACK backend
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module will be built only if ScaLAPACK is available. Uses the ScaLAPACK routine \c pdsyevd and several supporting
!> ScaLAPACK and BLACS routines. Uses a compile-time macro for the BLAS/LAPACK integer, which is long integer by default,
!> or controllable by the compiler flags -Dblas64int/-Dblas32int.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module SCALAPACKDiagonalizer_module

    use blas_lapack_gbl,           only: blasint
    use const_gbl,                 only: stdout
    use consts_mpi_ci,             only: PASS_TO_CDENPROP, SAVE_CONTINUUM_COEFFS, SAVE_L2_COEFFS
    use mpi_gbl,                   only: master, myrank, mpi_xermsg
    use precisn,                   only: wp
    use BaseIntegral_module,       only: BaseIntegral
    use BaseMatrix_module,         only: BaseMatrix
    use Diagonalizer_module,       only: BaseDiagonalizer
    use DiagonalizerResult_module, only: DiagonalizerResult
    use Options_module,            only: Options
    use Parallelization_module,    only: grid => process_grid
    use SCALAPACKMatrix_module,    only: SCALAPACKMatrix
    use Timing_module,             only: master_timer

    implicit none

    public SCALAPACKDiagonalizer

    private

    type, extends(BaseDiagonalizer) :: SCALAPACKDiagonalizer
    contains
        procedure, public :: diagonalize => diagonalize_SCALAPACK
        procedure, public :: diagonalize_backend => diagonalize_backend_SCALAPACK
        procedure, public :: process_solution
    end type SCALAPACKDiagonalizer

contains

    subroutine diagonalize_SCALAPACK (this, matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)

        class(SCALAPACKDiagonalizer)    :: this
        class(DiagonalizerResult)       :: dresult
        class(BaseMatrix),   intent(in) :: matrix_elements
        class(BaseIntegral), intent(in) :: integrals
        type(Options),       intent(in) :: option
        integer,             intent(in) :: num_eigenpair
        logical,             intent(in) :: all_procs

        real(wp), allocatable :: z_mat(:,:), w(:)
        integer :: ierr

        write (stdout, "('Utilizing Optimized SCALAPACK matrix')")
        write (stdout, "('N: ',I0)") matrix_elements % get_matrix_size()
        write (stdout, "('Requested # of eigenpairs ',I0)") num_eigenpair

        ! call the diagonalization backend
        select type(matrix_elements)
            class is (SCALAPACKMatrix)
                call this % diagonalize_backend(matrix_elements, num_eigenpair, z_mat, w, all_procs, option, integrals)
                call this % process_solution(matrix_elements, num_eigenpair, z_mat, w, dresult, all_procs, option, integrals)
        end select

    end subroutine diagonalize_SCALAPACK


    subroutine diagonalize_backend_SCALAPACK (this, matrix_elements, num_eigenpair, z_mat, w, all_procs, option, integrals)

        class(SCALAPACKDiagonalizer)       :: this
        class(SCALAPACKMatrix), intent(in) :: matrix_elements
        class(BaseIntegral),    intent(in) :: integrals
        type(Options),          intent(in) :: option
        integer,                intent(in) :: num_eigenpair
        logical,                intent(in) :: all_procs
        real(wp), allocatable,  intent(inout) :: z_mat(:,:), w(:)

        real(wp),         allocatable :: work(:), pdsyevd_lwork(:), pdormtr_lwork(:)
        integer(blasint), allocatable :: iwork(:)
        integer(blasint)              :: lwork, liwork, npcol, mat_dimen, info, loc_r, loc_c, one = 1
        integer                       :: ierr

        write (stdout, "('Diagonalization will be done with SCALAPACK')")

        npcol     = grid % gpcols
        loc_r     = matrix_elements % local_row_dimen
        loc_c     = matrix_elements % local_col_dimen
        mat_dimen = matrix_elements % get_matrix_size()

        ! allocate the eigenvectors
        allocate(z_mat(loc_r, loc_c), w(mat_dimen), stat = ierr)
        if (ierr /= 0) then
            call mpi_xermsg('SCALAPACKDiagonalizer_module', 'diagonalize_backend_SCALAPACK', 'Memory allocation error.', ierr, 1)
        end if

        ! Calculating the workspace size is quite a bit of alchemy (no pun intended). The ScaLAPACK subroutine `pdsyevd` has some
        ! well-defined requirements which it correctly returns on a workspace query. However, when eigenvectors are desired,
        ! `pdsyevd` eventually calls `pdormtr` and requirements of that subroutine are not fully included in the "optimal" work
        ! size returned by `pdsyevd`. From the inspection of the ScaLAPACK code, `pdsyevd` will always require 2*mat_dimen elements
        ! on top of any `pdormtr` workspace. So, we query both subroutines for the real workspace size and use the larger of
        !
        !       lwork[pdsyevd]
        !
        ! and
        !
        !       lwork[pdormtr] + 2*mat_dimen

        liwork = 7*mat_dimen + 8*npcol + 2

        allocate(pdsyevd_lwork(1), pdormtr_lwork(1), iwork(liwork))

        call pdsyevd('V', 'L', mat_dimen, matrix_elements % a_local_matrix, one, one, matrix_elements % descr_a_mat, w, z_mat, &
                     one, one, matrix_elements % descr_z_mat, pdsyevd_lwork, -one, iwork, liwork, info)
        call pdormtr('L', 'L', 'N', mat_dimen, mat_dimen, matrix_elements % a_local_matrix, one, one, &
                     matrix_elements % descr_a_mat, w, z_mat, one, one, matrix_elements % descr_z_mat, pdormtr_lwork, -one, info)

        lwork = ceiling(max(pdsyevd_lwork(1), pdormtr_lwork(1) + 2 * mat_dimen))

        allocate(work(lwork), stat = info)

        if (info == 0) then
            write (stdout, "('ScaLAPACK workspace size: lwork = ',I0,', liwork = ',I0)") lwork, liwork
        else
            call mpi_xermsg('SCALAPACKDiagonalizer_module', 'DoSCALAPACKMat', 'Memory allocation error', int(info), 1)
        end if

        ! And now the distributed diagonalization itself.

        call pdsyevd('V', 'L', mat_dimen, matrix_elements % a_local_matrix, one, one, matrix_elements % descr_a_mat, w, z_mat, &
                     one, one, matrix_elements % descr_z_mat, work, lwork, iwork, liwork, info)

        if (info == 0) then
            write (stdout, "(/'Diagonalization finished successfully!')")
        else
            call mpi_xermsg('SCALAPACKDiagonalizer_module', 'DoSCALAPACKMat', 'PDSYEVD returned non-zero code', int(info), 1)
        end if

    end subroutine diagonalize_backend_SCALAPACK


    subroutine process_solution (this, matrix_elements, num_eigenpair, z_mat, w, dresult, all_procs, option, integrals)

        class(SCALAPACKDiagonalizer)       :: this
        class(DiagonalizerResult)          :: dresult
        class(SCALAPACKMatrix), intent(in) :: matrix_elements
        class(BaseIntegral),    intent(in) :: integrals
        type(Options),          intent(in) :: option
        integer,                intent(in) :: num_eigenpair
        logical,                intent(in) :: all_procs
        real(wp), allocatable,  intent(inout) :: z_mat(:,:), w(:)

        real(wp), allocatable :: global_vector(:)

        integer(blasint) :: myrow, mycol, nprow, npcol, blacs_context, i, j, iprow, ipcol
        integer(blasint) :: mat_dimen, nb, loc_r, loc_c, jdimen, idimen

        blacs_context = grid % gcntxt
        myrow         = grid % mygrow
        mycol         = grid % mygcol
        nprow         = grid % gprows
        npcol         = grid % gpcols

        mat_dimen     = matrix_elements % get_matrix_size()
        loc_r         = matrix_elements % local_row_dimen
        loc_c         = matrix_elements % local_col_dimen
        nb            = matrix_elements % scal_block_size

        if (iand(dresult % vector_storage, SAVE_CONTINUUM_COEFFS) /= 0 .and. &
            iand(dresult % vector_storage, SAVE_L2_COEFFS) == 0) then
            write (stdout, '(/," Only the continuum part of the eigenvectors will be saved to disk.")')
            write (stdout, '(" Index of the last continuum configuration: ",i10)') dresult % continuum_dimen
        end if

        ! save the results to disk
        allocate(global_vector(mat_dimen))
        call blacs_barrier(blacs_context, 'a')

        ! store Hamiltonian information into the CDENPROP vector
        if (iand(option % vector_storage_method, PASS_TO_CDENPROP) /= 0) then
            call dresult % export_header(option, integrals)
        end if

        ! write Hamiltonian information into the output disk file (the current context master will do this)
        if (grid % grank == master) then
            call dresult % write_header(option, integrals)
        end if

        ! save eigenvalues by master of the context (to disk), or by everyone (to internal arrays)
        if (all_procs .or. grid % grank == master) then
            call dresult % handle_eigenvalues(w(1:num_eigenpair), matrix_elements % diagonal, num_eigenpair, int(mat_dimen))
        end if

        ! reconstruct and save the eigenvectors to the disk file
        call master_timer % start_timer("Collect eigenvectors")
        do jdimen = 1, num_eigenpair

            ! skip collection of eigenvectors when not required
            if (iand(dresult % vector_storage, SAVE_L2_COEFFS) == 0 .and. &
                iand(dresult % vector_storage, SAVE_CONTINUUM_COEFFS) == 0 .and. &
                .not. all_procs) then
                if (grid % grank == master) then
                    call dresult % handle_eigenvector(global_vector, int(mat_dimen))  ! a dummy call to write empty record
                end if
                cycle
            end if

            ! clear the local copy of the eigenvector
            global_vector = 0

            ! set non-zero elements of global_vector by processes that own those elements
            do idimen = 1, mat_dimen
                call infog2l(idimen, jdimen, matrix_elements % descr_z_mat, nprow, npcol, myrow, mycol, i, j, iprow, ipcol)
                if (myrow == iprow .and. mycol == ipcol) then
                    global_vector(idimen) = z_mat(i,j)
                end if
            end do

            ! collect the eigenvector...
            if (all_procs) then
                ! ... to all members of the context and let all of them 'handle' the eigenvector
                call dgsum2d(blacs_context, 'all', ' ', mat_dimen, 1_blasint, global_vector, mat_dimen, -1_blasint, -1_blasint)
                call dresult % handle_eigenvector(global_vector, int(mat_dimen))
            else
                ! ... to master of the context and let it 'handle' the eigenvector
                call dgsum2d(blacs_context, 'all', ' ', mat_dimen, 1_blasint, global_vector, mat_dimen, 0_blasint, 0_blasint)
                if (grid % grank == master) then
                    call dresult % handle_eigenvector(global_vector, int(mat_dimen))
                end if
            end if

        end do
        call master_timer % stop_timer("Collect eigenvectors")

        ! fill cdenprop data structure for further processing of eigenvectors, if required
        if (iand(dresult % vector_storage, PASS_TO_CDENPROP) /= 0) then
            !Distributed eigenvectors and the descriptor:
            if (allocated(dresult % ci_vec % CV)) deallocate(dresult % ci_vec % CV)
            call move_alloc(z_mat, dresult % ci_vec % CV)
            dresult % ci_vec % blacs_context   = blacs_context
            dresult % ci_vec % local_row_dimen = loc_r
            dresult % ci_vec % local_col_dimen = loc_c
            dresult % ci_vec % scal_block_size = nb
            dresult % ci_vec % myrow           = myrow
            dresult % ci_vec % mycol           = mycol
            dresult % ci_vec % nprow           = nprow
            dresult % ci_vec % npcol           = npcol
            dresult % ci_vec % mat_dimen       = matrix_elements % get_matrix_size()
            dresult % ci_vec % descr_CV_mat    = matrix_elements % descr_z_mat
            dresult % ci_vec % lda             = matrix_elements % lda
            dresult % ci_vec % CV_is_scalapack = .true.

            dresult % ci_vec % mat_dimen_r = dresult % ci_vec % mat_dimen
            dresult % ci_vec % mat_dimen_c = dresult % ci_vec % mat_dimen

            !Eigenvalues and phases:
            call dresult % export_eigenvalues(w(1:num_eigenpair), matrix_elements % diagonal, num_eigenpair, &
                                                int(mat_dimen), dresult % ci_vec % ei, dresult % ci_vec % iphz)

            write (stdout, '(/," Eigensolutions exported into the CIVect structure.",/)')
        end if

        ! finish writing the solution file and let other processes access it
        if (grid % grank == master) then
            call dresult % finalize_solutions(option)
        end if

        if (allocated(z_mat)) deallocate(z_mat)

        call blacs_barrier(blacs_context, 'a')

    end subroutine process_solution

end module SCALAPACKDiagonalizer_module
