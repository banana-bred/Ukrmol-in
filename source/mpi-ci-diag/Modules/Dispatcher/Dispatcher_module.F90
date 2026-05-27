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

!> \brief   Dispatcher module
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module provides the automatic selection of the correct classes needed by MPI SCATCI.
!> These classes are the Integrals, Matricies and Diagonalizers
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial revision.
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Dispatcher_module

    use precisn,                          only: longint, wp
    use const_gbl,                        only: stdout
    use mpi_gbl,                          only: nprocs, mpi_xermsg
    use consts_mpi_ci
    use Options_module,                   only: Options
    use BaseMatrix_module,                only: BaseMatrix
    use WriterMatrix_module,              only: WriterMatrix
    use Diagonalizer_module,              only: BaseDiagonalizer
    use LAPACKDiagonalizer_module,        only: LapackDiagonalizer
    use DavidsonDiagonalizer_module,      only: DavidsonDiagonalizer
#ifdef arpack
    use ARPACKDiagonalizer_module,        only: ARPACKDiagonalizer
#endif
    use MatrixElement_module,             only: MatrixElementVector
    use Timing_Module,                    only: master_timer
    use BaseIntegral_module,              only: BaseIntegral
    use SWEDEN_module,                    only: SWEDENIntegral
    use ALCHEMY_module,                   only: ALCHEMYIntegral
    use UKRMOL_module,                    only: UKRMOLIntegral
   !use BaseCorePotential_module,         only: BaseCorePotential
   !use MOLPROCorePotential_module,       only: MOLPROCorePotential
#ifdef slepc
    use SLEPCMatrix_module,               only: SLEPCMatrix, initialize_slepc
    use SLEPCDiagonalizer_module,         only: SLEPCDiagonalizer
   !use SLEPCDavidsonDiagonalizer_module, only: SLEPCDavidsonDiagonalizer
#endif
#ifdef scalapack
    use SCALAPACKMatrix_module,           only: SCALAPACKMatrix
    use SCALAPACKDiagonalizer_module,     only: SCALAPACKDiagonalizer
#endif
#ifdef useelpa
    use ELPAMatrix_module,                only: ELPAMatrix
    use ELPADiagonalizer_module,          only: ELPADiagonalizer, initialize_elpa
#endif

    implicit none

    public DispatchMatrixAndDiagonalizer, DispatchIntegral, initialize_libraries

    private

contains

    !> \brief   Initialize libraries used by MPI-SCATCI
    !> \authors J Benda
    !> \date    2019
    !>
    !> Currently only SLEPc needs global initialization.
    !>
    subroutine initialize_libraries
#if defined(usempi) && defined(slepc)
        call initialize_slepc
#endif
#if defined(usempi) && defined(scalapack) && defined(useelpa)
        call initialize_elpa
#endif
    end subroutine initialize_libraries


    !> \brief   This subroutine determines which matrix format/diagonalizer pair to be used by SCATCI in build diagonalization
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> It determines what type of diagonalizer based on the matrix size and number of eigenpairs
    !> and whether to use the Legacy SCATCI diagonalizers for nprocs = 1 / -Dnompi option or to use an MPI one-
    !>
    !> \param[in]  diag_choice          From Options, used to force a particular diagonalizer.
    !> \param[in]  force_serial         Whether to use serial method even in distributed run.
    !> \param[in]  matrix_size          Size of Matrix. Used with number_of_eigenpairs to select the correct matrix format.
    !> \param[in]  number_of_eigenpairs No. of eigensolutions requested.
    !> \param[out] matrix               The selected matrix format.
    !> \param[out] diagonalizer         The selected diagonalizer.
    !> \param[in]  integral             A constructed and loaded integral class. This is required to write out the matrix header.
    !> \param[in]  matrix_io            The IO unit for the matrix header.
    !>
    subroutine DispatchMatrixAndDiagonalizer (diag_choice, force_serial, matrix_size, number_of_eigenpairs, matrix, diagonalizer, &
                                              integral, matrix_io)

        use Parallelization_module, only: grid => process_grid

        class(BaseMatrix),       pointer, intent(out) :: matrix
        class(BaseDiagonalizer), pointer, intent(out) :: diagonalizer
        class(BaseIntegral),              intent(in)  :: integral

        integer, intent(in) :: diag_choice, force_serial
        integer, intent(in) :: matrix_size, number_of_eigenpairs
        integer, intent(in) :: matrix_io
        integer             :: final_choice

        !Now if we forced a diagonalizer then use that, other wise use the rules
        if (diag_choice /= SCATCI_DECIDE) then
            final_choice = diag_choice
        else
            if (number_of_eigenpairs > real(matrix_size) * 0.2) then ! More the 20% of eigenvectors requires a dense diagonalizer
                final_choice = DENSE_DIAG
            else if (number_of_eigenpairs <= 3) then ! Less than three, we use davidson
                final_choice = DAVIDSON_DIAG
            else !Inbetween we use an iterative library like ARPACK or SLEPC
                final_choice = ITERATIVE_DIAG
            end if
        end if

#ifdef usempi
        !If we are dealing with an MPI diagonalization then we select the appropriate version of the diagonalizer
        if (nprocs > 1 .and. force_serial == 0) then
            select case (final_choice)
                case (DENSE_DIAG)                ! The MPI DENSE diagonalization is SCALAPACK
#ifdef scalapack
#ifdef useelpa
                    write (stdout, "('ELPA chosen as diagonalizer of choice. Hohoho!')")
                    if (grid % gprows > matrix_size .or. grid % gpcols > matrix_size) then
                        write (stdout, "('Too many processes for ELPA (BLACS grid size > matrix size). Defaulting to SCALAPACK.')")
                        write (stdout, "(3(a,i0))") 'nprocs = ', grid%gprows, '×', grid%gpcols, ' and matrix_size = ', matrix_size
                        allocate(SCALAPACKMatrix::matrix)
                        allocate(SCALAPACKDiagonalizer::diagonalizer)
                    else
                        allocate(ELPAMatrix::matrix)
                        allocate(ELPADiagonalizer::diagonalizer)
                    end if
#else
                    write (stdout, "('ScaLAPACK chosen as diagonalizer of choice. Nice!')")
                    allocate(SCALAPACKMatrix::matrix)
                    allocate(SCALAPACKDiagonalizer::diagonalizer)
#endif
#else
                    call mpi_xermsg('Dispatcher_module', 'DispatchMatrixAndDiagonalizer', &
                                    'ScaLAPACK needed for dense diagonalization!', 1, 1)
#endif
                case (DAVIDSON_DIAG) !Use the SCATCI Davidson diagonalizer
                    write (stdout, "('Davidson chosen as diagonalizer of choice. Sweet!')")
#ifdef slepc
                    write (stdout, "('SLEPc chosen as diagonalizer of choice. Amazing!')")
                    allocate(SLEPCMatrix::matrix)
                    allocate(SLEPCDiagonalizer::diagonalizer)
#elif defined(scalapack)
#ifdef useelpa
                    write (stdout, "('Compiled without SLEPc, defaulting to ELPA')")
                    if (grid % gprows > matrix_size .or. grid % gpcols > matrix_size) then
                        write (stdout, "('Too many processes for ELPA (BLACS grid size > matrix size). Defaulting to SCALAPACK.')")
                        write (stdout, "(3(a,i0))") 'nprocs = ', grid%gprows, '×', grid%gpcols, ' and matrix_size = ', matrix_size
                        allocate(SCALAPACKMatrix::matrix)
                        allocate(SCALAPACKDiagonalizer::diagonalizer)
                    else
                        allocate(ELPAMatrix::matrix)
                        allocate(ELPADiagonalizer::diagonalizer)
                    end if
#else
                    write (stdout, "('Compiled without SLEPc, defaulting to ScaLAPACK')")
                    allocate(SCALAPACKMatrix::matrix)
                    allocate(SCALAPACKDiagonalizer::diagonalizer)
#endif
#else
                    call mpi_xermsg('Dispatcher_module', 'DispatchMatrixAndDiagonalizer', &
                                    'SLEPc or ScaLAPACK needed for distributed diagonalization!', 1, 1)
#endif
                case (ITERATIVE_DIAG)                ! ITERATIVE uses SLEPC
#ifdef slepc
                    write (stdout, "('SLEPc chosen as diagonalizer of choice. Amazing!')")
                    allocate(SLEPCMatrix::matrix)
                    allocate(SLEPCDiagonalizer::diagonalizer)
#elif defined(scalapack)
#ifdef useelpa
                    write (stdout, "('Compiled without SLEPc, defaulting to ELPA')")
                    if (grid % gprows > matrix_size .or. grid % gpcols > matrix_size) then
                        write (stdout, "('Too many processes for ELPA (BLACS grid size > matrix size). Defaulting to SCALAPACK.')")
                        write (stdout, "(3(a,i0))") 'nprocs = ', grid%gprows, '×', grid%gpcols, ' and matrix_size = ', matrix_size
                        allocate(SCALAPACKMatrix::matrix)
                        allocate(SCALAPACKDiagonalizer::diagonalizer)
                    else
                        allocate(ELPAMatrix::matrix)
                        allocate(ELPADiagonalizer::diagonalizer)
                    end if
#else
                    write (stdout, "('Compiled without SLEPc, defaulting to ScaLAPACK')")
                    allocate(SCALAPACKMatrix::matrix)
                    allocate(SCALAPACKDiagonalizer::diagonalizer)
#endif
#else
                    call mpi_xermsg('Dispatcher_module', 'DispatchMatrixAndDiagonalizer', &
                                    'SLEPc or ScaLAPACK needed for distributed diagonalization!', 1, 1)
#endif
                case default
                    write (stdout, "('That is not a diagonalizer I would have chosen tbh, lets not do it')")
                    call mpi_xermsg('Dispatcher_module', 'DispatchMatrixAndDiagonalizer', &
                                    'Invalid diagonalizer flag chosen', 1, 1)
            end select
            return
        end if
#endif

        select case (final_choice)
            case (DENSE_DIAG)  !Use our rolled LAPACK diagonalizer
                write (stdout, "('LAPACK chosen as diagonalizer of choice. Nice!')")
                !allocate(MatrixElementVector::matrix)
                call integral % write_matrix_header(matrix_io, matrix_size)
                allocate(WriterMatrix::matrix)
                allocate(LAPACKDiagonalizer::diagonalizer)
            case (DAVIDSON_DIAG) !Use the SCATCI Davidson diagonalizer
                write (stdout, "('Davidson chosen as diagonalizer of choice. Sweet!')")
                call integral % write_matrix_header(matrix_io, matrix_size)
                allocate(WriterMatrix::matrix)
                allocate(DavidsonDiagonalizer::diagonalizer)
            case (ITERATIVE_DIAG) !Use the SCATCI ARPACK diagonalizer
#ifdef arpack
                write (stdout, "('ARPACK chosen as diagonalizer of choice. Amazing!')")
                call integral % write_matrix_header(matrix_io, matrix_size)
                allocate(WriterMatrix::matrix)
                allocate(ARPACKDiagonalizer::diagonalizer)
#else
                write (stdout, "('ARPACK not available - using LAPACK. Sorry!')")
                call integral % write_matrix_header(matrix_io, matrix_size)
                allocate(WriterMatrix::matrix)
                allocate(LAPACKDiagonalizer::diagonalizer)
#endif
            case default
                write (stdout, "('That is not a diagonalizer I would have chosen tbh, lets not do it')")
                stop "Invalid diagonalizer flag chosen"
        end select

    end subroutine DispatchMatrixAndDiagonalizer


    !> \brief   This subroutine dispatches the correct integral class based on simple parameters
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[in]  sym_group_flag        Which symmetry group, this determines either ALCHEMY or SWEDEN/UKRMOL+.
    !> \param[in]  ukrmolplus_integrals  For D_2h_, determines either SWEDEN or UKRMOL+.
    !> \param[out] integral              The resultant integral object.
    !>
    subroutine DispatchIntegral (sym_group_flag, ukrmolplus_integrals, integral)
        class(BaseIntegral), pointer, intent(out) :: integral
        integer, intent(in) :: sym_group_flag
        logical, intent(in) :: ukrmolplus_integrals

        if (sym_group_flag /= SYMTYPE_D2H) then
            allocate(ALCHEMYIntegral::integral)
        else
            if (ukrmolplus_integrals) then
                allocate(UKRMOLIntegral::integral)
            else
                allocate(SWEDENIntegral::integral)
            end if
        end if

    end subroutine DispatchIntegral


    !subroutine DispatchECP(ecp_type,ecp_filename,vals_defined,num_targ_symm,num_targ_states_symm,spatial_symm,spin_symm,ecp)
    !    integer,intent(in)                ::    ecp_type
    !    character(len=line_len),intent(in)        ::    ecp_filename
    !    logical,intent(in)                ::    vals_defined
    !    integer,intent(in)                ::    num_targ_symm,num_targ_states_symm(num_targ_symm),spatial_symm(num_targ_symm),spin_symm(num_targ_symm)
    !    class(BaseCorePotential),pointer,intent(out)    ::    ecp
    !
    !    if(ecp_type == ECP_TYPE_NULL) then
    !        ecp => null()
    !        return
    !    else if(vals_defined == .true.) then
    !        select case(ecp_type)
    !            case(ECP_TYPE_MOLPRO)
    !                allocate(MOLPROCorePotential::ecp)
    !            case default
    !                stop "Unrecognized ECP type"
    !        end select
    !
    !        call ecp%construct(num_targ_symm,num_targ_states_symm,spatial_symm,spin_symm)
    !
    !        call ecp%parse_ecp(ecp_filename)
    !    else
    !        write(stdout,"('TARGET SPATIAL SYMMETRY AND MULTIPLICITY NOT DEFINED IN INPUT')")
    !        stop "TARGET SPATIAL SYMMETRY AND MULTIPLICITY NOT DEFINED IN INPUT"
    !    endif
    !
    !end subroutine DispatchECP

end module Dispatcher_module
