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

!> \brief   Diagonalizer type using LAPACK backend
!> \authors A Al-Refaie
!> \date    2017
!>
!> This type is always available. Uses the LAPACK routine \c dsyevd. Uses a compile-time macro for the BLAS/LAPACK integer,
!> which is long integer by default, or controllable by the compiler flags -Dblas64int/-Dblas32int.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module LapackDiagonalizer_module

    use blas_lapack_gbl,           only: blasint
    use precisn,                   only: longint, wp
    use const_gbl,                 only: stdout
    use BaseIntegral_module,       only: BaseIntegral
    use BaseMatrix_module,         only: BaseMatrix
    use WriterMatrix_module,       only: WriterMatrix
    use Diagonalizer_module,       only: BaseDiagonalizer
    use DiagonalizerResult_module, only: DiagonalizerResult
    use Options_module,            only: Options
    use mpi_gbl,                   only: master, myrank, nprocs, mpi_mod_bcast, mpi_mod_barrier
    use consts_mpi_ci,             only: PASS_TO_CDENPROP, SAVE_CONTINUUM_COEFFS, SAVE_L2_COEFFS

    implicit none

    type, extends(BaseDiagonalizer) :: LapackDiagonalizer
    contains
        procedure         :: diagonalize => choose_matelem
        procedure         :: diagonalize_generic
        procedure, public :: diagonalize_writer
       !procedure, public :: choose_matelem
    end type

contains

    subroutine choose_matelem (this, matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)
        class(LapackDiagonalizer)       :: this
        class(DiagonalizerResult)       :: dresult
        class(BaseMatrix),   intent(in) :: matrix_elements
        class(BaseIntegral), intent(in) :: integrals
        type(Options),       intent(in) :: option
        integer,             intent(in) :: num_eigenpair
        logical,             intent(in) :: all_procs

        select type(matrix_elements)
            type is (WriterMatrix)
                call this % diagonalize_writer(matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)
            !class is (BaseMatrix)
                !call this%diagonalize_generic(matrix_elements,num_eigenpair,result,all_procs)
            !class is (BaseMatrix)
                !call this%DoNonSLEPCMat(matrix_elements,num_eigenpair,eigen,vecs,maxit,maxtol)
        end select
    end subroutine choose_matelem


    subroutine diagonalize_writer (this, matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)
        class(LapackDiagonalizer)       :: this
        class(DiagonalizerResult)       :: dresult
        class(WriterMatrix), intent(in) :: matrix_elements
        class(BaseIntegral), intent(in) :: integrals
        type(Options),       intent(in) :: option
        integer,             intent(in) :: num_eigenpair
        logical,             intent(in) :: all_procs
        real(wp)                        :: EONE
        real(wp), dimension(1)          :: eshift
        real(wp), allocatable, target   :: hmt(:,:), eig(:)
        integer                         :: matrix_unit
        integer                         :: matrix_size, error, ido, jdo
        integer                         :: num_matrix_elements_per_record, num_elems

        character(len=4), dimension(30) :: NAMP
        integer,          dimension(10) :: NHD
        integer,          dimension(20) :: KEYCSF, NHE
        real(kind=wp),    dimension(41) :: DTNUC

        eshift = 0

        write (stdout, "('Diagonalization done with LAPACK')")
        write (stdout, "('Parameters:')")
        write (stdout, "('N: ',i8)") matrix_elements % get_matrix_size()
        write (stdout, "('Requested # of eigenpairs',i8)") num_eigenpair

        allocate(eig(matrix_elements % get_matrix_size()))
        matrix_size = 0

        if (myrank == master) then
            !Get our marix unit
            matrix_unit = matrix_elements % get_matrix_unit()

            !Lets follow scatcis example even though we can get
            !these information from the matrix class itself
            rewind matrix_unit
            read (matrix_unit) matrix_size, num_matrix_elements_per_record, NHD, NAMP, NHE, DTNUC

            num_elems = matrix_elements % get_size()

            allocate(hmt(matrix_size,matrix_size))

            call HMAT(hmt, 1, 1, eshift, 0, matrix_size, stdout, matrix_unit, num_matrix_elements_per_record, eone)
            call QLDIAG(matrix_size, hmt, eig)

            if (iand(option % vector_storage_method, PASS_TO_CDENPROP) /= 0) then
                call dresult % export_header(option, integrals)
            end if
            call dresult % write_header(option, integrals)
            call dresult % handle_eigenvalues(eig(1:num_eigenpair), matrix_elements % diagonal, num_eigenpair, matrix_size)

            do ido = 1, num_eigenpair
                call dresult % handle_eigenvector(hmt(:,ido), matrix_size)
            end do
        end if

        call mpi_mod_barrier(error)

        if (iand(dresult % vector_storage, SAVE_CONTINUUM_COEFFS) /= 0 .and. &
            iand(dresult % vector_storage, SAVE_L2_COEFFS) == 0) then
            write (stdout, '(/," Only the continuum part of the eigenvectors will be saved to disk.")')
            write (stdout, '(" Index of the last continuum configuration: ",i10)') dresult % continuum_dimen
        end if

        if (all_procs) then
            if (myrank /= master) then
                matrix_size = matrix_elements%get_matrix_size()
                allocate(hmt(matrix_size,matrix_size))
            end if

            call mpi_mod_barrier(error)
            call mpi_mod_bcast(eig, master)

            if (myrank /= master) then
                call dresult % handle_eigenvalues(EIG(1:num_eigenpair), matrix_elements % diagonal, num_eigenpair, matrix_size)
            end if

            do ido = 1, num_eigenpair
                call mpi_mod_bcast(hmt(1:matrix_size,ido), master)
                if (myrank /= master) call dresult % handle_eigenvector(hmt(1:matrix_size,ido), matrix_size)
            end do
        end if

        if (iand(dresult % vector_storage, PASS_TO_CDENPROP) /= 0) then
            !Eigenvectors:
            if (allocated(dresult % ci_vec % CV)) deallocate(dresult % ci_vec % CV)
            call move_alloc(hmt, dresult % ci_vec % CV)
            dresult % ci_vec % CV_is_scalapack = .false.

            !Eigenvalues and phases:
            call dresult % export_eigenvalues(eig(1:num_eigenpair), matrix_elements % diagonal, &
                                              num_eigenpair, matrix_size, dresult % ci_vec % ei, dresult % ci_vec % iphz)

            write (stdout, '(/," Eigensolutions exported into the CIVect structure.",/)')
        end if

        deallocate(eig)

    end subroutine diagonalize_writer


    subroutine diagonalize_generic (this, matrix_elements, num_eigenpair, eigen, vecs, all_procs)
        class(LapackDiagonalizer)        :: this
        class(BaseMatrix), intent(in)    :: matrix_elements
        integer,           intent(in)    :: num_eigenpair
        real(wp),          intent(inout) :: eigen(:)
        real(wp),          intent(inout) :: vecs(:,:)
        logical,           intent(in)    :: all_procs

        real(wp)                      :: coeff
        real(wp), pointer             :: pointer_mat(:,:)
        real(wp), allocatable, target :: matrix(:)
        real(wp), allocatable         :: t_eigen(:), work(:)
        integer(blasint), allocatable :: iwork(:)
        integer(blasint)              :: matrix_size, nelems, lwork, liwork
        integer                       :: num_elements, ido, jdo, i, j, ierror, info

        !Since we're dealing with lapack we dont need to worry about the optional arguments

        matrix_size = matrix_elements % get_matrix_size()

        !Write out the parameters
        write (stdout, "('Diagonalization done with LAPACK')")
        write (stdout, "('Parameters:')")
        write (stdout, "('N: ',i8)") matrix_size
        write (stdout, "('Requested # of eigenpairs',i8)") num_eigenpair
        write (stdout, "('Number of matrix elements',i8)") matrix_elements % get_size()

        nelems = int(matrix_size, longint) * int(matrix_size, longint)

        if (nprocs > 1) then
            stop "Serial diagonalizer used in MPI, use SCALAPACK instead"
        end if

        !Lets allocate our matrix, this potentially uses a lot of memory
        allocate(matrix(nelems), stat = ierror)
        pointer_mat(1:matrix_size,1:matrix_size) => matrix

        if (ierror /= 0) then
            stop "LapackDiagonalizer --matrix-- Out of memory"
        end if

        allocate(t_eigen(matrix_size), stat = ierror)

        if (ierror /= 0) then
            stop "LapackDiagonalizer --eigen_-- Out of memory"
        end if

        num_elements = matrix_elements % get_size()
        t_eigen(:) = 0.0
        matrix = 0.0
        eigen(:) = 0.0

        !lets fill the matrix with the coefficients
        do ido = 1, num_elements
            call matrix_elements % get_matrix_element(ido, i, j, coeff)
            pointer_mat(i,j) = pointer_mat(i,j) + coeff
            !if (matrix_size == 475) write(91,"(2i8,D16.8)") i, j, val(ido)
        end do

        lwork = -1
        liwork = -1

        !Our starting work array
        allocate(work(1), iwork(1))

        write (stdout, "('Diagonalizing......')", advance = 'no')

        call dsyevd('V', 'L', matrix_size, matrix,  matrix_size, t_eigen, WORK, lwork,iWORK, LIWORK, INFO)

        if (info /= 0) then
            write (stdout, "(' dsyev returned ',i8)") info
            stop 'lapack_dsyev - dsyev failed'
        end if

        LWORK = int(WORK(1))
        LIWORK = int(IWORK(1))

        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))

        call dsyevd('V', 'L', matrix_size, matrix,  matrix_size, t_eigen, WORK, lwork, iWORK, LIWORK, INFO)

        if (info /= 0) then
            write (stdout, "(' dsyev returned ',i8)") info
            stop 'lapack_dsyev - dsyev failed'
        end if

        write (stdout, "('done!')")

        eigen(1:num_eigenpair) = t_eigen(1:num_eigenpair)

        do ido = 1, num_eigenpair
            do jdo = 1, matrix_size
                vecs(jdo,ido)=pointer_mat(jdo,ido)
            end do
        end do

        !Clean up
        deallocate(matrix)
        deallocate(t_eigen)
        deallocate(work)

    end subroutine diagonalize_generic

end module LapackDiagonalizer_module
