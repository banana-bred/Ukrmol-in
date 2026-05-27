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

!> \brief   Target Rmat CI module
!> \authors A Al-Refaie
!> \date    2017
!>
!> Handles the storage of the target coefficients for the contracted Hamiltonian.
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial Revision
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Target_RMat_CI_module

    use precisn,                   only: wp
    use const_gbl,                 only: stdout
    use consts_mpi_ci,             only: SYMTYPE_D2H
    use scatci_routines,           only: CIRMAT
    use Options_module,            only: Options, Phase
    use MemoryManager_module,      only: master_memory
    use DiagonalizerResult_module, only: DiagonalizerResult
   !use BaseCorePotential_module

    implicit none

    public Target_RMat_CI, read_ci_mat

    private

    !> \brief   This class handles the storage of target CI coefficients.
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    type, extends(DiagonalizerResult) :: Target_RMat_CI
        integer  :: set         !< Which ci set this belongs to
        integer  :: nstat       !< Number of eigenvalues
        integer  :: num_csfs    !< Size of the eigenvectors
        real(wp) :: core_energy !< The associated core energy (for printing)
        integer  :: current_eigenvector = 0

        real(wp), allocatable :: eigenvectors(:,:)  !< The ci coefficients
        real(wp), allocatable :: eigenvalues(:)     !< The associated eigenvalues
        type(Phase), pointer  :: phase              !< Not used, a vestigal type

       !class(BaseCorePotential), pointer :: effective_core_potential => null()
    contains
        procedure, public  :: initialize      => initialize_cirmat
        procedure, public  :: print           => print_cirmat
        procedure, public  :: get_coefficient => get_coefficient_ci
        procedure, public  :: print_vecs
        procedure, public  :: export_eigenvalues
        procedure, public  :: handle_eigenvalues => store_eigenvalues
        procedure, public  :: handle_eigenvector => store_eigenvector
        procedure, public  :: modify_diagonal
        procedure, public  :: modify_L2_diagonal
    end type Target_RMat_CI

contains

    !> \brief   Sets up and allocates eigenvalues and eigenvectors
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this      Target CI object to initialize.
    !> \param[in] set          Which CI set this belongs to
    !> \param[in] nstat        Number of eigenpairs/target states
    !> \param[in] num_csfs     Size of the eigenvectors/number of ci coefficients per target state
    !> \param[in] phase_       Not used
    !> \param[in] core_energy  Core energy used for eigenvalue printing
    !>
    subroutine initialize_cirmat (this, set, nstat, num_csfs, phase_, core_energy)
        class(Target_RMat_CI)    :: this
        integer,     intent(in)  :: set, nstat, num_csfs
        real(wp),    intent(in)  :: core_energy
        type(Phase), intent(in),  target   :: phase_
       !class(BaseCorePotential), pointer  :: ecp

        integer :: ifail

        !Set relevant values
        this % set = set
        this % nstat = nstat
        this % num_csfs = num_csfs
        this % core_energy = core_energy
        this % phase => phase_

       !if (associated(ecp)) then
       !    this % effective_core_potential => ecp
       !end if

        !allocate the necessary arrays and track the memory usage
        allocate(this % eigenvectors(this % num_csfs, this % nstat), this % eigenvalues(this % nstat), stat = ifail)
        call master_memory % track_memory(storage_size(this % eigenvectors)/8, &
                                          size(this % eigenvectors), ifail, 'TARGET_RMAT_CI::EIGENVECTORS')
        call master_memory % track_memory(storage_size(this % eigenvalues)/8,  &
                                          size(this % eigenvalues),  ifail, 'TARGET_RMAT_CI::EIGENVALUES')

        !Zero the arrays
        this % eigenvectors(:,:) = 0.0_wp
        this % eigenvalues(:) = 0.0_wp

    end subroutine initialize_cirmat


    !> \brief   Prints the eigen-energies of the target states stored
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine print_cirmat (this)
        class(Target_RMat_CI) :: this
        integer               :: I, J

        write (stdout, '(/," SET",I4,4x,A)') &
            this % set
        write (stdout, '(/," NOCSF=",I5,4x,"NSTAT=",I5,4x,"EN   =",F20.10)') &
            this % num_csfs, this % nstat, this % core_energy
        write (stdout, '(/," EIGEN-ENERGIES",/,(16X,5F20.10))') &
            (this % eigenvalues(I) + this % core_energy, I = 1, this % nstat)
        write (stdout, '(  " STATE",I3,":",I3," transformation vectors selected")') &
            this % set, this % nstat

    end subroutine print_cirmat


    !> \brief   Print vectors to standard output
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine print_vecs (this)
        class(Target_RMat_CI) :: this
        integer               :: ido, jdo

        do ido = 1, this % nstat
            write (stdout, *) ' Target state ', this % set, ido
            do jdo = 1, this % num_csfs
                write (stdout, *) this % eigenvectors(jdo, ido)
            end do
        end do

    end subroutine print_vecs


    !> \brief   Get the specific coefficient for a target state and configuration
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[in] this              Target CI object to query.
    !> \param[in] target_state      Which target state to get from.
    !> \param[in] configuration_idx Which particular coefficient to get.
    !>
    !> \result The resultant coefficient.
    !>
    real(wp) function get_coefficient_CI (this, target_state, configuration_idx)
        class(Target_RMat_CI) :: this
        integer, intent(in)   :: target_state, configuration_idx

        get_coefficient_CI = this % eigenvectors(configuration_idx, target_state)

    end function  get_coefficient_CI


    !> \brief   Reads the target coefficients from file
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[in] option  An initialized Options object.
    !> \param[in] cirmats An allocated and initialized array of target CI objects.
    !>
    subroutine read_ci_mat (option, cirmats)
        class(Options),        intent(in)    :: option
        class(Target_RMat_CI), intent(inout) :: cirmats(:)

        real(wp), allocatable :: eigenvectors(:)
        real(wp), allocatable :: eigenvalues(:)

        integer :: ido, jdo, kdo, ii, jj, ci_size, nstat, ifail

        !Allocate the temporary arrays for eigenvalues and eienvectors
        allocate(eigenvalues(option % total_num_tgt_states), eigenvectors(option % total_num_ci_target), stat = ifail)
        call master_memory % track_memory(storage_size(eigenvectors)/8, size(eigenvectors), ifail, &
                                          'TARGET_RMAT_CI::READ::EIGENVECTORS')
        call master_memory % track_memory(storage_size(eigenvalues)/8, size(eigenvalues), ifail, &
                                          'TARGET_RMAT_CI::READ::EIGENVALUES')

        !Use SCATCI to read them
        call CIRMAT(option % num_target_sym,        &
                    option % ci_target_switch,      &
                    option % ci_set_number,         &
                    eigenvectors,                   &
                    eigenvalues,                    &
                    option % ci_sequence_number,    &
                    option % num_ci_target_sym,     &
                    stdout,                         &
                    option % num_target_state_sym,  &
                    option % phase_index,           &
                    option % num_continuum,         &
                    option % sym_group_flag)

        ii = 0
        jj = 0

        !Seperate out all the eigenvectors and values into distinct target symmetries and states.
        do ido = 1, option % num_target_sym
            ci_size = cirmats(ido) % num_csfs
            nstat   = cirmats(ido) % nstat
            do jdo = 1, nstat
                jj = jj + 1
                cirmats(ido) % eigenvalues(jdo) = eigenvalues(jj)
                cirmats(ido) % eigenvectors(1:ci_size,jdo) = eigenvectors(ii+1:ii+ci_size)
                ii = ii + ci_size
            end do
        end do

        !Stop tracking this memory
        call master_memory % free_memory(storage_size(eigenvectors)/8, size(eigenvectors))
        call master_memory % free_memory(storage_size(eigenvalues)/8,  size(eigenvalues))

        !Cleanup
        deallocate(eigenvalues)
        deallocate(eigenvectors)

    end subroutine read_ci_mat


    !> \brief   Save eigenvalues from diagonalizer
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Copy eigenvalues obtained by the diagonalizer to internal storage of Target_RMat_CI.
    !>
    subroutine store_eigenvalues (this, eigenvalues, diagonals, num_eigenpairs, vec_dimen)
        class(Target_RMat_CI) :: this
        integer,  intent(in)  :: num_eigenpairs, vec_dimen
        real(wp), intent(in)  :: eigenvalues(num_eigenpairs), diagonals(vec_dimen)

        this % eigenvalues = eigenvalues

       !if (associated(this % effective_core_potential)) then
       !       call this % effective_core_potential % modify_target_energies(this % set, this % eigenvalues)
       !end if

    end subroutine store_eigenvalues


    !> \brief   Save eigenvector from diagonalizer
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Copy eigenvector obtained by the diagonalizer to internal storage of Target_RMat_CI.
    !>
    subroutine store_eigenvector (this, eigenvector, vec_dimen)
        class(Target_RMat_CI) :: this
        integer,  intent(in)  :: vec_dimen
        real(wp), intent(in)  :: eigenvector(vec_dimen)

        this % current_eigenvector = this % current_eigenvector + 1
        this % eigenvectors(:, this % current_eigenvector) = eigenvector(:)

    end subroutine store_eigenvector


    !> \brief   To be overriden by derived types
    !> \authors J Benda
    !> \date    2018
    !>
    !> This subroutine body was added so that Target_RMat_CI does not need to be `abstract`.
    !> Otherwise there were problems with compilation using non-Intel compilers.
    !>
    subroutine export_eigenvalues (this, eigenvalues, diagonals, num_eigenpairs, vec_dimen, ei, iphz)
        use precisn, only: wp
        class(Target_RMat_CI)  :: this
        integer,  intent(in)   :: num_eigenpairs, vec_dimen
        real(wp), intent(in)   :: eigenvalues(num_eigenpairs), diagonals(vec_dimen)
        real(wp), allocatable  :: ei(:)
        integer,  allocatable  :: iphz(:)
    end subroutine export_eigenvalues


    !> \brief   Dummy subroutine
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Has no effect.
    !>
    subroutine modify_diagonal (this, value)
        class(Target_RMat_CI), intent(in)    :: this
        real(wp),              intent(inout) :: value

       !if (associated(this % effective_core_potential)) then
       !    call this % effective_core_potential % modify_ecp_diagonal(value)
       !end if

    end subroutine modify_diagonal


    !> \brief   Dummy subroutine
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Has no effect.
    !>
    subroutine modify_L2_diagonal (this, value)
        class(Target_RMat_CI), intent(in)    :: this
        real(wp),              intent(inout) :: value

       !if (associated(this % effective_core_potential)) then
       !       call this % effective_core_potential % modify_ecp_diagonal_L2(value)
       !end if

    end subroutine modify_L2_diagonal

end module Target_RMat_CI_module
