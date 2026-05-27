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

!> \brief   Diagonalizer type using SLEPc backend
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module will be built only if SLEPc is available.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module SLEPCDiagonalizer_module

    use blas_lapack_gbl,           only: blasint
    use const_gbl,                 only: stdout
    use consts_mpi_ci,             only: PASS_TO_CDENPROP
    use cdenprop_defs,             only: CIvect
    use mpi_gbl,                   only: master, mpi_xermsg
    use precisn,                   only: longint, wp

    use BaseMatrix_module,         only: BaseMatrix
    use BaseIntegral_module,       only: BaseIntegral
    use Diagonalizer_module,       only: BaseDiagonalizer
    use DiagonalizerResult_module, only: DiagonalizerResult
    use Options_module,            only: Options
    use Parallelization_module,    only: grid => process_grid
    use Timing_module,             only: master_timer
    use SLEPCMatrix_module,        only: SLEPCMatrix

    use petsc
    use slepceps

#include <finclude/petsc.h>
#include <finclude/slepceps.h>

    implicit none

    public SLEPCDiagonalizer

    private

    type, extends(BaseDiagonalizer) :: SLEPCDiagonalizer
    contains
        procedure, public :: diagonalize => diagonalize_slepc
       !procedure, public :: DoNonSLEPCMat
        procedure, public :: UseSLEPCMat
        procedure, public :: SelectEPS
    end type

contains

    subroutine UseSLEPCMat (this, matrix_elements, mat_out)
        class(SLEPCDiagonalizer)          :: this
        class(SLEPCMatrix), intent(in)    :: matrix_elements
        Mat, pointer,       intent(inout) :: mat_out

        write(stdout,"('--------Optimized SLEPC Matrix Format chosen------')")

        mat_out => matrix_elements % get_PETSC_matrix()

        if (.not. associated(mat_out)) then
            stop "Matrix doesn;t exist for SLEPC, are you sure you created it?"
        end if

    end subroutine UseSLEPCMat


    subroutine SelectEPS (this, eps)
        class(SLEPCDiagonalizer) :: this
        EPS, intent(inout)       :: eps
        PetscErrorCode           :: ierr

        write(stdout,"('--------KRYLOVSCHUR used as Diagonalizer------')")

        !     ** Create eigensolver context
        call EPSSetType(eps,EPSKRYLOVSCHUR,ierr)
        !call EPSSetTrueResidual(eps,PETSC_TRUE,ierr);
        !call EPSLanczosSetReorthog(eps,EPS_LANCZOS_REORTHOG_SELECTIVE,ierr)
        if (ierr /= 0) then
            call mpi_xermsg('SLEPCDiagonalizer_module', 'SelectEPS', 'Failed to select eigensolver', int(ierr), 1)
        end if

    end subroutine SelectEPS


    subroutine diagonalize_slepc (this, matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)
        class(SLEPCDiagonalizer)        :: this
        class(DiagonalizerResult)       :: dresult
        class(BaseMatrix),   intent(in) :: matrix_elements
        class(BaseIntegral), intent(in) :: integrals
        type(Options),       intent(in) :: option
        integer,             intent(in) :: num_eigenpair
        logical,             intent(in) :: all_procs

        integer               :: mat_dimen, number_solved
        integer(blasint)      :: grid_master_context
        real(wp), allocatable :: diagonal_temp(:), eigen(:)

        type(CIvect) :: eigenvector

        !-----SLEPC values
        Mat, pointer         :: A
        EPS                  :: eps
        EPSType              :: tname
        ST                   :: st_eps
        EPSConvergedReason   :: conv_reason
        PetscReal            :: maxtol,tol, error
        PetscScalar          :: kr, ki, current_shift
        PetscScalar, pointer :: xx_v(:)
        Vec                  :: xr, xi, f_xr
        PetscInt             :: nev, maxit, its, nconv, matrix_size, ido, n
        PetscErrorCode       :: ierr
        VecScatter           :: vsctx
        KSP                  :: ksp
        PC                   :: pc

        maxit = option % max_iterations
        maxtol = option % max_tolerance

        if (maxit < 0) maxit = PETSC_DEFAULT_INTEGER

        matrix_size = matrix_elements % get_matrix_size()
        mat_dimen = matrix_size

        write (stdout, "()")
        write (stdout, "('-----------WARNING----------------------- ')")
        write (stdout, "('When dealing with symmetries that are degenerate SLEPC may give incorrect values')")
        write (stdout, "('YOU HAVE BEEN WARNED')")
        write (stdout, "('---------------------------------------- ')")
        write (stdout, "('N: ',i8)") matrix_size
        write (stdout, "('Requested # of eigenpairs',i8)") num_eigenpair
        write (stdout, "('Maximum iterations',i8)") maxit
        write (stdout, "('Maximum tolerance',es9.2)") maxtol
        write (stdout, "()")

        select type(matrix_elements)
            type is (SLEPCMatrix)
                call this % UseSLEPCMat(matrix_elements, A)
            class is (BaseMatrix)
                stop "Matrix format not yet implemented for SLEPC"
            !class is (BaseMatrix)
                !call this%DoNonSLEPCMat(matrix_elements,num_eigenpair,eigen,vecs,maxit,maxtol)
        end select

        call MatCreateVecs(A, xr, xi, ierr)

        if (all_procs) then
            call VecScatterCreateToAll(xr, vsctx, f_xr, ierr)
        else
            call VecScatterCreateToZero(xr, vsctx, f_xr, ierr)
        end if

        call EPSCreate(int(grid % gcomm, kind(PETSC_COMM_WORLD)), eps, ierr)
        call this % SelectEPS(eps)
        call EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr)
        call EPSSetProblemType(eps, EPS_HEP, ierr)

!#if defined(PETSC_HAVE_MUMPS)
!            if(num_eigenpair > 20) then
!                call EPSGetST(eps,st_eps,ierr);
!                call STSetType(st_eps,STSINVERT,ierr)
!                call EPSSetST(eps,st_eps,ierr)
!                call EPSSetWhichEigenpairs(eps,EPS_ALL,ierr);
!                   call EPSSetInterval(eps,-90,-80,ierr);
!                call STGetKSP(st_eps,ksp,ierr);
!                call KSPSetType(ksp,KSPPREONLY,ierr);
!                call KSPGetPC(ksp,pc,ierr);
!
!                call PCSetType(pc,PCCHOLESKY,ierr);
!
!                call EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE,ierr)
!                call PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS,ierr)
!
!                 call PetscOptionsInsertString("-mat_mumps_icntl_13 1 -mat_mumps_icntl_24 1 -mat_mumps_cntl_3 1e-12",ierr);
!
!            else
!                call EPSSetWhichEigenpairs(eps,    EPS_SMALLEST_REAL ,ierr);
!            endif
!#else
!                call EPSSetWhichEigenpairs(eps,    EPS_SMALLEST_REAL ,ierr);
!#endif

        !call MatView(A,    PETSC_VIEWER_DRAW_WORLD,ierr)

        !** Set operators. In this case, it is a standard eigenvalue problem

        !Create the ST
        n = num_eigenpair
        call EPSSetDimensions(eps, n, PETSC_DECIDE, PETSC_DECIDE, ierr)
        !  call EPSSetType(eps,EPSJD,ierr)
        call EPSSetTolerances(eps, maxtol, maxit, ierr)
        call master_timer % start_timer("SLEPC solve")
        !call EPSView(eps,    PETSC_VIEWER_STDOUT_WORLD,ierr)
        call EPSSolve(eps, ierr)
        if (ierr /= 0) then
            call mpi_xermsg('SLEPCDiagonalizer_module', 'SelectEPS', 'Failed to select eigensolver', int(ierr), 1)
        end if
        call master_timer % stop_timer("SLEPC solve")
        !call EPSGetConvergedReason(eps,conv_reason,ierr)

        ! write(stdout,"('Converged reason ',i8)") conv_reason
        call EPSGetIterationNumber(eps, its, ierr)
        call EPSGetTolerances(eps, tol, maxit, ierr)

        if (grid % grank == master) then
            write (stdout, *)
            write (stdout, "(('/ Stopping condition: tol=',1P,E11.4,', maxit=',I0,', stopped at it=',I0))") tol, maxit, its
        end if

        call EPSGetConverged(eps, nconv, ierr)

        if (nconv < num_eigenpair) then
            write (stdout, *) 'Not all requested eigenpairs have converged!!!!'
            write (stdout, *) 'Only ', nconv, ' have converged against ', num_eigenpair, ' requested'
            stop "Diagoanlization not completed"
        end if

        allocate(eigen(num_eigenpair))

        do ido = 0, min(nconv, n) - 1
            ! ** Get converged eigenpairs: i-th eigenvalue is stored in kr
            ! ** (real part) and ki (imaginary part)
            call EPSGetEigenpair(eps, ido, kr, ki, xr, xi, ierr)
            eigen(ido + 1) = PetscRealPart(kr)
        end do

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
            call dresult % handle_eigenvalues(eigen, matrix_elements % diagonal, num_eigenpair, mat_dimen)
        end if

        ! when the post-processing stage is required, set up the needed ScaLAPACK matrices
        if (iand(dresult % vector_storage, PASS_TO_CDENPROP) /= 0) then

            ! split off a 1-by-1 (i.e. master only) sub-grid
            grid_master_context = grid % gcntxt
#if defined(usempi) && defined(scalapack)
            call blacs_gridinit(grid_master_context, 'R', 1_blasint, 1_blasint)
#endif

            ! set up ScaLAPACK storage for a single eigenvector (master-only ScaLAPACK matrix)
            eigenvector % blacs_context = grid_master_context
            eigenvector % CV_is_scalapack = .true.
            call eigenvector % init_CV(mat_dimen, 1)

            ! set up ScaLAPACK storage for all eigenvectors, distributed among processes in the current BLACS grid
            dresult % ci_vec % blacs_context   = grid % gcntxt
            dresult % ci_vec % CV_is_scalapack = .true.
            call dresult % ci_vec % init_CV(mat_dimen, num_eigenpair)

        end if

        do ido = 0, min(nconv, n) - 1
            ! ** Get converged eigenpairs: i-th eigenvalue is stored in kr
            ! ** (real part) and ki (imaginary part)
            call EPSGetEigenpair(eps, ido, kr, ki, xr, xi, ierr)

            ! ** Compute the relative error associated to each eigenpair
            !call EPSComputeError(eps,ido,EPS_ERROR_RELATIVE,error,ierr)
            !write(stdout,*) PetscRealPart(kr), error
            !160      format (1P,'   ',E12.4,'       ',E12.4)
            !eigen(ido+1) = PetscRealPart(kr)
            ! write(stdout,*)eigen(ido+1)

            !Grab the vectors
            call VecScatterBegin(vsctx, xr, f_xr, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call VecScatterEnd(vsctx, xr, f_xr, INSERT_VALUES, SCATTER_FORWARD, ierr)
            if (all_procs .or. grid % grank == master) then
                call VecGetArrayRead(f_xr, xx_v, ierr)
                call dresult % handle_eigenvector(xx_v(1:mat_dimen), mat_dimen)
            end if
            if (iand(dresult % vector_storage, PASS_TO_CDENPROP) /= 0) then
                if (grid % grank == master) then
                    eigenvector % CV(1:mat_dimen, 1) = xx_v(1:mat_dimen)
                end if
                call dresult % ci_vec % redistribute(eigenvector, grid % gcntxt, mat_dimen, 1, 1, 1, 1, ido + 1)
            end if
            if (all_procs .or. grid % grank == master) then
                call VecRestoreArrayRead(f_xr, xx_v, ierr)
            end if
        end do

        if (iand(dresult % vector_storage, PASS_TO_CDENPROP) /= 0) then
            call dresult % export_eigenvalues(eigen(1:num_eigenpair), matrix_elements % diagonal, num_eigenpair, &
                                              mat_dimen, dresult % ci_vec % ei, dresult % ci_vec % iphz)
            write (stdout, '(/," Eigensolutions exported into the CIVect structure.",/)')
        end if

        ! finish writing the solution file and let other processes access it
        if (grid % grank == master) then
            call dresult % finalize_solutions(option)
        end if

        call EPSDestroy(eps, ierr)
        call VecDestroy(xr, ierr)
        call VecDestroy(xi, ierr)

        deallocate(eigen)

    end subroutine diagonalize_slepc

end module SLEPCDiagonalizer_module
