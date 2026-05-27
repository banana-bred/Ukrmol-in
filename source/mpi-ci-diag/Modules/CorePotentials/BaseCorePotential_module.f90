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

!> \brief   Base type for core potentials
!> \authors A Al-Refaie
!> \date    2017
!>
!> Used by specialized core potential types, eg. MOLPROCorePotential.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module BaseCorePotential_module

    use precisn, only: wp

    implicit none

    type, abstract :: BaseCorePotential
        real(wp), allocatable :: target_energies(:,:)
        real(wp)              :: L2_ECP_deficit
        integer               :: num_target_sym
        integer, allocatable  :: num_target_per_symmetries(:)
        integer, allocatable  :: spatial_symmetry(:)
        integer, allocatable  :: spin_symmetry(:)
        real(wp)              :: ground_energy
    contains
        procedure, public                   :: construct => construct_core
        procedure(generic_ecp),    deferred :: parse_ecp
        procedure(generic_input),  deferred :: parse_input
        procedure(generic_energy), deferred :: modify_target_energies
        procedure, public                   :: modify_ecp_diagonal
        procedure, public                   :: modify_ecp_diagonal_L2
    end type BaseCorePotential

    abstract interface
        !> \brief   Main build routine of the hamiltonian
        !> \authors A Al-Refaie
        !> \date    2017
        !> \public
        !>
        !> All build must be done within this routine in order to be used by MPI-SCATCI.
        !>
        !> \param[out] matrix_elements Resulting matrix elements from the build
        !>
        subroutine generic_ecp (this, filename)
            import :: BaseCorePotential
            class(BaseCorePotential) :: this
            character(len=*)         :: filename
        end subroutine generic_ecp
    end interface

    abstract interface
        !> \brief   Main build routine of the hamiltonian
        !> \authors A Al-Refaie
        !> \date    2017
        !> \public
        !>
        !> All build must be done within this routine in order to be used by MPI-SCATCI.
        !>
        !> \param[out] matrix_elements Resulting matrix elements from the build
        !>
        subroutine generic_input(this)
            import :: BaseCorePotential
            class(BaseCorePotential) :: this
        end subroutine generic_input
    end interface

    abstract interface
        subroutine generic_energy (this, target_sym, target_energy)
            import :: BaseCorePotential
            class(BaseCorePotential) :: this
            integer, intent(in)      :: target_sym
            real(wp)                 :: target_energy(:)
        end subroutine generic_energy
    end interface

contains

    subroutine construct_core (this, num_target_sym, num_target_per_symmetries, spatial_symmetry, spin_symmetry)
        class(BaseCorePotential) :: this
        integer, intent(in) :: num_target_sym
        integer, intent(in) :: num_target_per_symmetries(num_target_sym), &
                               spatial_symmetry(num_target_sym), &
                               spin_symmetry(num_target_sym)

        this % num_target_sym = num_target_sym
        allocate(this % num_target_per_symmetries(num_target_sym))

        this % num_target_per_symmetries = num_target_per_symmetries
        allocate(this % target_energies(maxval(this % num_target_per_symmetries), this % num_target_sym))

        allocate(this % spatial_symmetry(num_target_sym))
        allocate(this % spin_symmetry(num_target_sym))

        this % spatial_symmetry = spatial_symmetry
        this % spin_symmetry = spin_symmetry

        this % target_energies = 0.0
        this % L2_ECP_deficit =  0.0
    end subroutine construct_core


    subroutine modify_ecp_diagonal (this, value)
        class(BaseCorePotential), intent(in) :: this
        real(wp), intent(inout) :: value

        value = value - this % ground_energy
        value = value + this % target_energies(1,1)
    end subroutine modify_ecp_diagonal


    subroutine modify_ecp_diagonal_L2 (this, value)
        class(BaseCorePotential),intent(in) :: this
        real(wp), intent(inout) :: value

        value = value + this % L2_ECP_deficit
    end subroutine modify_ecp_diagonal_L2

end module BaseCorePotential_module
