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

!> \brief   Uncontracted Hamiltonian module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Uncontracted_Hamiltonian_module

    use const_gbl,              only: stdout
    use consts_mpi_ci,          only: MAT_DENSE
    use precisn,                only: wp
    use BaseMatrix_module,      only: BaseMatrix
    use Hamiltonian_module,     only: BaseHamiltonian
    use Parallelization_module, only: grid => process_grid
    use Symbolic_module,        only: SymbolicElementVector
    use Utility_module,         only: compute_total_box, box_index_to_ij

    implicit none

    public Uncontracted_Hamiltonian

    private

    type, extends(BaseHamiltonian) :: Uncontracted_Hamiltonian
    contains
        procedure, public :: build_hamiltonian => build_Uncontracted_hamiltonian
    end type Uncontracted_Hamiltonian

contains

    !> \brief   Initialize the type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine build_Uncontracted_hamiltonian (this, matrix_elements)
        class(Uncontracted_Hamiltonian)  :: this
        class(BaseMatrix), intent(inout) :: matrix_elements
        type(SymbolicElementVector)      :: symbolic_elements, ref_symbolic_elements
        integer  :: starting_index, num_csfs, ido, jdo, num_continuum, csf_a_idx, csf_b_idx
        integer  :: num_elements, loop_skip, my_idx, loop_ido, total_vals, m_flag, diag
        real(wp) :: mat_coeff, element_one

        diag = 0
        num_continuum = this % options % last_continuum
        num_csfs = this % options % num_csfs
        num_elements = 1
        this % orbitals % MFLG = 0

        call matrix_elements % initialize_matrix_structure(num_csfs, MAT_DENSE, num_csfs)
        !get the first element and slater it regardless of idiag flag
        call ref_symbolic_elements % construct
        call this % slater_rules(this % csfs(1), this % csfs(1), ref_symbolic_elements, diag, .false.)

        element_one = this % evaluate_integrals(ref_symbolic_elements, 0)

        call symbolic_elements % construct
        call symbolic_elements % clear

        loop_skip = max(1, grid % gprocs)
        my_idx = max(grid % grank, 0)

        total_vals = compute_total_box(num_csfs, num_csfs)
        do loop_ido = 1, total_vals, loop_skip

            call matrix_elements % update_pure_L2(.false.)
            ido = loop_ido + my_idx
            if (ido > total_vals) cycle

            call box_index_to_ij(ido, num_csfs, csf_a_idx, csf_b_idx)

            if (csf_a_idx == csf_b_idx) then
                m_flag = 0
            else
                m_flag = 1
            end if
            if (csf_a_idx < csf_b_idx) cycle

            call this % slater_rules(this % csfs(csf_a_idx), this % csfs(csf_b_idx), symbolic_elements, m_flag, .false.)

            if (m_flag == 0 .and. this % diagonal_flag == 0) then
                call symbolic_elements % add_symbols(ref_symbolic_elements, -1.0_wp)
            end if

            mat_coeff = this % evaluate_integrals(symbolic_elements, this % options % phase_correction_flag)
            if (this % diagonal_flag == 0 .and. this % orbitals % MFLG == 0) mat_coeff = mat_coeff + element_one
            call matrix_elements % insert_matrix_element(csf_a_idx, csf_b_idx, mat_coeff)
            num_elements = num_elements + 1

            call symbolic_elements % clear

        end do

        !Flush
        call matrix_elements % update_pure_L2(.true.)

        !call matrix_elements%print
        call symbolic_elements % destroy
        call matrix_elements % finalize_matrix

        write (stdout, "('Num of elements = ',i0)") num_elements

    end subroutine build_Uncontracted_hamiltonian

end module Uncontracted_Hamiltonian_module
