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

!> \brief   Base integral module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module BaseIntegral_module

    use consts_mpi_ci,          only: NAME_LEN_MAX, NIDX
    use mpi_gbl,                only: master, myrank
    use precisn,                only: longint, wp
    use Options_module,         only: Options
    use Parallelization_module, only: grid => process_grid
    use Utility_module,         only: unpack_ints

    implicit none

    public BaseIntegral

    type, abstract :: BaseIntegral
        integer, allocatable :: orbital_mapping(:)
        integer              :: num_orbitals
        real(wp)             :: core_energy
        integer              :: positron_flag
        logical              :: quantamoln_flag

        !>Matrix header information
        integer              :: matrix_size
        integer              :: nhe(20), nhd(10)
        integer              :: num_symmetries
        integer              :: number_of_matrix_records
        integer              :: nnuc
        real(wp)             :: DTNUC(41)

        character(NAME_LEN_MAX) :: name
    contains
        procedure(generic_finalize),   deferred :: finalize_self
        procedure(generic_initialize), deferred :: initialize_self
        procedure(generic_load),       deferred :: load_integrals
        procedure(generic_get),        deferred :: get_integral_ijklm
        procedure(generic_geometries), deferred :: write_geometries
        procedure(generic_destroy),    deferred :: destroy_integrals

        procedure, public :: initialize          => initialize_base
        procedure, public :: finalize            => finalize_base
        procedure, public :: write_matrix_header => base_write_header
        procedure, public :: get_core_energy
        procedure, public :: get_num_nuclei
        procedure, public :: get_integralf
    end type BaseIntegral

    abstract interface
        subroutine generic_initialize (this, option)
            use Options_module, only: Options
            import :: BaseIntegral
            class(BaseIntegral)        :: this
            class(Options), intent(in) :: option
        end subroutine generic_initialize
    end interface

    abstract interface
        subroutine generic_finalize (this)
            use Options_module, only: Options
            import :: BaseIntegral
            class(BaseIntegral)        :: this
        end subroutine generic_finalize
    end interface

    abstract interface
        subroutine generic_load (this, iounit)
            import :: BaseIntegral
            class(BaseIntegral) :: this
            integer, intent(in) :: iounit
        end subroutine generic_load
    end interface

    abstract interface
        function generic_get (this, i, j, k, l, m) result(coeff)
            use precisn, only : wp
            import :: BaseIntegral
            class(BaseIntegral) :: this
            integer, intent(in) :: i, j, k, l, m
            real(wp)            :: coeff
        end function generic_get
    end interface

    abstract interface
        subroutine generic_geometries (this, iounit)
            use precisn, only : wp
            import :: BaseIntegral
            class(BaseIntegral), intent(in) :: this
            integer,             intent(in) :: iounit
        end subroutine generic_geometries
    end interface

    abstract interface
        subroutine generic_destroy (this)
            import :: BaseIntegral
            class(BaseIntegral) :: this
        end subroutine generic_destroy
    end interface

contains

    !> \brief   Initialize the class
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine initialize_base (this, option, num_orbitals, mapping)

        class(BaseIntegral)          :: this
        class(Options),   intent(in) :: option
        integer,          intent(in) :: num_orbitals
        integer,          intent(in) :: mapping(num_orbitals)

        this % num_orbitals = num_orbitals

        allocate(this % orbital_mapping(num_orbitals))

        this % orbital_mapping(:) = mapping(:)
        this % positron_flag = option % positron_flag
        this % quantamoln_flag = option % QuantaMolN

        !>Store needed matrix header information
        this % num_symmetries = option % num_syms

        this % number_of_matrix_records = option % num_matrix_elements_per_rec
        this % name = option % name

        call this % initialize_self(option)

    end subroutine initialize_base


    !> \brief   Finalize the class
    !> \authors J Benda
    !> \date    2019
    !>
    !> Deallocate memory used by the integrals storage.
    !>
    subroutine finalize_base (this)

        class(BaseIntegral) :: this

        call this % finalize_self

        if (allocated(this % orbital_mapping)) then
            deallocate (this % orbital_mapping)
        end if

    end subroutine finalize_base


    !> \brief   Write data file header
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine base_write_header (this, matrix_io, matrix_size)

        class(BaseIntegral), intent(in) :: this
        integer,             intent(in) :: matrix_io, matrix_size

        if (grid % grank == master) then
            open (unit = matrix_io, form = 'unformatted')
            write (matrix_io) matrix_size, this % number_of_matrix_records , 0, matrix_size, &
                              0, this % num_symmetries, 0, 0, 0, 0, this % nnuc, 0, &
                              this % NAME, this % NHE, this % DTNUC
        end if

    end subroutine base_write_header


    !> \brief   Retrieve integral value
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    function get_integralf (this, integral) result(coeff)

        class(BaseIntegral)          :: this
        integer(longint), intent(in) :: integral(NIDX)
        real(wp) :: coeff
        integer  :: label(8)

        call unpack_ints(integral, label)
        coeff = this % get_integral_ijklm(label(1), label(2), label(3), label(4), label(5))

    end function get_integralf


    !> \brief   Retrieve integral value
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    real(wp) function get_core_energy (this)

        class(BaseIntegral) :: this

        get_core_energy = this % core_energy

    end function get_core_energy


    !> \brief   Retrieve integral value
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer function get_num_nuclei (this)

        class(BaseIntegral) :: this

        get_num_nuclei = this % nnuc

    end function get_num_nuclei

end module BaseIntegral_module
