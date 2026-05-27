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

!> \brief   Target CI Hamiltonian module
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module handles the computation of multiple target symmetry hamiltonians in a single run.
!>
!> \note 01/02/2017 - Ahmed Al-Refaie: Initial revision.
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module CI_Hamiltonian_module

    use const_gbl,              only: stdout
    use consts_mpi_ci,          only: NORMAL_CSF, MATRIX_EVALUATE_DIFF, MATRIX_EVALUATE_FULL, NO_PURE_TARGET_INTEGRAL, MAT_DENSE
    use precisn,                only: wp
    use BaseMatrix_module,      only: BaseMatrix
    use Hamiltonian_module,     only: BaseHamiltonian
    use Options_module,         only: Options
    use Parallelization_module, only: grid => process_grid
    use Symbolic_module,        only: SymbolicElementVector
    use Utility_module,         only: compute_total_box, box_index_to_ij

    implicit none

    public Target_CI_Hamiltonian

    private

    !> \brief   Hamiltonian type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This class differs from the Uncontracted as it allows for the computation of
    !> multiple targets of differing symmetries in a single run.
    !>
    type, extends(BaseHamiltonian) :: Target_CI_Hamiltonian
        integer :: target_symmetry  !< The target symmetry hamiltonian to compute
        integer :: csf_skip         !< The CSF index stride
    contains
        procedure, public :: initialize => initialize_target_CI_hamiltonian
        procedure, public :: build_hamiltonian => build_target_CI_hamiltonian_fast
    end type Target_CI_Hamiltonian

contains

    !> \brief   Sets up the hamiltonian to build for a specific target symmetry
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this         hamiltonian object to update.
    !> \param[in] target_symmetry The target symmetry hamiltonian to build
    !>
    subroutine initialize_target_CI_hamiltonian (this, target_symmetry)

        class(Target_CI_Hamiltonian) :: this
        integer, intent(in) :: target_symmetry

        this % target_symmetry = target_symmetry
        if (this % options % csf_type == NORMAL_CSF) then
            this % csf_skip = this % options % num_continuum_orbitals_target(target_symmetry)
        else
            this % csf_skip = 2
        end if

        this % initialized = .true.
    end subroutine initialize_target_CI_hamiltonian


    !> \brief   Builds the Target Hamiltonian of a specific target symmetry
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this           Hamiltonian object to update.
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine build_target_CI_hamiltonian (this, matrix_elements)

        class(Target_CI_Hamiltonian)     :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        type(SymbolicElementVector)      :: symbolic_elements, ref_symbolic_elements
        integer  :: starting_index, num_csfs, ido, jdo, csf_a_idx, csf_b_idx, pzero_, num_elements
        real(wp) :: mat_coeff, element_one

        !Just a quick check in the old diagonal flag
        if (this % diagonal_flag == NO_PURE_TARGET_INTEGRAL) then
            this % diagonal_flag = MATRIX_EVALUATE_FULL
        else
            this % diagonal_flag = MATRIX_EVALUATE_DIFF
        end if

        this % phase_flag = 0

        !Get the total number of CSFs within the symmetry
        num_csfs = this % options % num_ci_target_sym(this % target_symmetry)

        !Initialize the matrix structure as DENSE with no block structure
        call matrix_elements % initialize_matrix_structure(num_csfs, MAT_DENSE, num_csfs)

        !Construct the symbolic elements and clear
        call symbolic_elements % construct
        call symbolic_elements % clear
        call ref_symbolic_elements % construct

        starting_index = 1
        !Get the start of the CSFs within the target symmetry
        do ido = 2, this % target_symmetry
            starting_index = starting_index + this % options % num_ci_target_sym(ido - 1) * this % csf_skip
        end do

        this % orbitals % MFLG = 0

        !Do the (1,1) diagonal --------------------------------------------------------!
        csf_a_idx = starting_index

        !Everyone must have this
        call this % slater_rules(this % csfs(csf_a_idx), this % csfs(csf_a_idx), ref_symbolic_elements, 0, .false.)

        !Evaluate the symbols for the integral for the matrix element (1,1)
        if (.not. ref_symbolic_elements % is_empty()) then
            mat_coeff = this % evaluate_integrals(ref_symbolic_elements, 0)
            this % element_one = mat_coeff
            !Lets have the master keep it, it will eventually reach everyone else if they need it
            if (grid % mygrow == 0 .and. grid % mygcol == 0) call matrix_elements % insert_matrix_element(1, 1, mat_coeff)
            call matrix_elements % update_pure_L2(.false.)
        end if
        !write(stdout,*)this%element_one

        !Clear the symbols
        call symbolic_elements % clear

        !------------Do the rest of diagonals---------------------------!
        do ido = 2, num_csfs
            call matrix_elements % update_pure_L2(.false.)
            csf_a_idx = csf_a_idx + this % csf_skip
            call this % slater_rules(this % csfs(csf_a_idx), this % csfs(csf_a_idx), symbolic_elements, 0)
            if (symbolic_elements % is_empty()) cycle

            if (this % diagonal_flag == 0) call symbolic_elements % add_symbols(ref_symbolic_elements, -1.0_wp)
            mat_coeff = this%evaluate_integrals(symbolic_elements,0)
            if (this % diagonal_flag == 0) mat_coeff = mat_coeff + this % element_one
            num_elements = num_elements + 1

            call matrix_elements % insert_matrix_element(ido, ido, mat_coeff)
            call symbolic_elements % clear
        end do

        this % orbitals % MFLG = 1

        csf_a_idx = starting_index - this % csf_skip
        !Do the off-idagonal
        do ido = 1, num_csfs - 1
            !Get the first CSf index
            csf_a_idx = csf_a_idx + this % csf_skip
            !Start the second CSF index
            csf_b_idx = csf_a_idx
            !Loop thorugh the lower triangular
            do jdo = ido + 1, num_csfs
                !Compute the second CSf index
                csf_b_idx = csf_b_idx + this % csf_skip
                !Perform an update if necessary
                call matrix_elements % update_pure_L2(.false.)
                !Slater the CSF
                call this % slater_rules(this % csfs(csf_a_idx), this % csfs(csf_b_idx), symbolic_elements, 1)
                !If empty, skip
                if (symbolic_elements % is_empty()) cycle

                num_elements = num_elements + 1
                !Evaulate the matrix coefficient
                mat_coeff = this % evaluate_integrals(symbolic_elements, 0)

                !Insert into the matrix
                call matrix_elements % insert_matrix_element(jdo, ido, mat_coeff)
                !Clear
                call symbolic_elements % clear
            end do
        end do
        !Force an update to clear the matrix cache if needed
        call matrix_elements % update_pure_L2(.true.)

        !Cleanup
        call symbolic_elements % destroy
        call ref_symbolic_elements % destroy

        !Finalize if needed
        call matrix_elements % finalize_matrix

    end subroutine build_target_CI_hamiltonian


    !> \brief   Builds the Target Hamiltonian of a specific target symmetry fast, currently not working
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \todo Fix this so the results are correct.
    !>
    !> \param[inout] this           Hamiltonian object to update.
    !> \param[out] matrix_elements  A BaseMatrix object, stores the calculated elements.
    !>
    subroutine build_target_CI_hamiltonian_fast (this, matrix_elements)

        class(Target_CI_Hamiltonian)     :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        type(SymbolicElementVector)      :: symbolic_elements, ref_symbolic_elements
        integer  :: starting_index, num_csfs, ido, jdo, csf_a_idx, csf_b_idx, pzero_, loop_idx, total_vals
        integer  :: num_elements,begin_csf,loop_skip,my_idx,loop_ido,m_flag
        real(wp) :: mat_coeff, element_one

        !Just a quick check in the old diagonal flag
        if (this % diagonal_flag == NO_PURE_TARGET_INTEGRAL) then
            this % diagonal_flag = MATRIX_EVALUATE_FULL
        else
            this % diagonal_flag = MATRIX_EVALUATE_DIFF
        end if

        this % phase_flag = 0
        num_elements = 0

        !Get the total number of CSFs within the symmetry
        num_csfs = this % options % num_ci_target_sym(this % target_symmetry)

        !Initialize the matrix structure as DENSE with no block structure
        call matrix_elements % initialize_matrix_structure(num_csfs, MAT_DENSE, num_csfs)

        !Construct the symbolic elements and clear
        call symbolic_elements % construct
        call symbolic_elements % clear
        call ref_symbolic_elements % construct

        starting_index = 1
        !Get the start of the CSFs within the target symmetry
        do ido = 2, this % target_symmetry
            starting_index = starting_index + this % options % num_ci_target_sym(ido - 1) * this % csf_skip
        end do

        m_flag = 0

        !Do the (1,1) diagonal --------------------------------------------------------!
        csf_a_idx = starting_index

        !Everyone must have this
        call this % slater_rules(this % csfs(csf_a_idx), this % csfs(csf_a_idx), ref_symbolic_elements, m_flag, .false.)

        !Evaluate the symbols for the integral for the matrix element (1,1)
        if (.not. ref_symbolic_elements % is_empty()) then
            mat_coeff = this % evaluate_integrals(ref_symbolic_elements, 0)
            this % element_one = mat_coeff
            !Lets have the master keep it, it will eventually reach everyone else if they need it
            if (grid % mygrow == 0 .and. grid % mygcol == 0) call matrix_elements % insert_matrix_element(1, 1, mat_coeff)
            call matrix_elements % update_pure_L2(.false.)
        end if
        !write(stdout,*)this%element_one

        !Clear the symbols
        call symbolic_elements % clear

        !------------Do the rest of diagonals---------------------------!
        do ido = 2,num_csfs
            call matrix_elements % update_pure_L2(.false.)
            csf_a_idx = csf_a_idx + this % csf_skip

            call this % slater_rules(this % csfs(csf_a_idx), this % csfs(csf_a_idx), symbolic_elements, 0)
            if (symbolic_elements % is_empty()) cycle

            if (this % diagonal_flag == 0) call symbolic_elements % add_symbols(ref_symbolic_elements, -1.0_wp)
            mat_coeff = this % evaluate_integrals(symbolic_elements, 0)
            if (this % diagonal_flag == 0) mat_coeff = mat_coeff + this % element_one
            num_elements = num_elements + 1

            call matrix_elements % insert_matrix_element(ido, ido, mat_coeff)
            call symbolic_elements % clear
        end do

        begin_csf = starting_index - this % csf_skip
        total_vals = compute_total_box(num_csfs, num_csfs)
        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        do loop_ido = 1, total_vals, loop_skip
            loop_idx = loop_ido + my_idx
            call matrix_elements % update_pure_L2(.false.)
            if (loop_idx > total_vals) cycle

            call box_index_to_ij(loop_idx, num_csfs, ido, jdo)
            if (ido <= jdo) cycle

            csf_a_idx = begin_csf + ido * this % csf_skip
            csf_b_idx = begin_csf + jdo * this % csf_skip

            call this % slater_rules(this % csfs(csf_a_idx), this % csfs(csf_b_idx), symbolic_elements, 1, .false.)

            if (symbolic_elements % is_empty()) cycle

            num_elements = num_elements + 1
            mat_coeff = this % evaluate_integrals(symbolic_elements, this % options % phase_correction_flag)

            call matrix_elements % insert_matrix_element(ido, jdo, mat_coeff)
            call symbolic_elements % clear
        end do

        call matrix_elements % update_pure_L2(.true.)
        call symbolic_elements % destroy
        call  matrix_elements % finalize_matrix

        write (stdout, "('Num of elements = ',i0)") num_elements

    end subroutine build_target_CI_hamiltonian_fast

end module CI_Hamiltonian_module
