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

!> \brief   UKRmol+ integral module
!> \authors A Al-Refaie
!> \date    2017
!>
!> The module makes use of the UKRmol+ integral library (GBTOlib).
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module UKRMOL_module

    use precisn,             only: longint, wp
    use const_gbl,           only: stdout
    use BaseIntegral_module, only: BaseIntegral
    use Options_module,      only: Options

    implicit none

    public UKRMOLIntegral

    private

    type, extends(BaseIntegral) :: UKRMOLIntegral
        integer  :: num_csf
        integer  :: num_one_electron_integrals
        integer  :: num_two_electron_integrals
        integer  :: symmetry_type
        integer  :: matrix_io_unit
        integer  :: num_symmetry
        real(wp) :: exotic_mass
        integer  :: hamiltonian_unit
        integer,          allocatable :: num_orbitals_sym(:)
        real(wp),         allocatable :: xnuc(:), ynuc(:), znuc(:), charge(:)
        character(len=8), allocatable, dimension(:) :: CNAME
    contains
        procedure, public :: initialize_self    => initialize_UKRMOL
        procedure, public :: finalize_self      => finalize_UKRMOL
        procedure, public :: load_integrals     => load_UKRMOL
        procedure, public :: get_integral_ijklm => get_integral_UKRMOL
        procedure, public :: write_geometries   => write_geometries_UKRMOL
        procedure, public :: destroy_integrals  => destroy_integral_UKRMOL
    end type

contains

    !> \brief   Initialize the type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine initialize_UKRMOL (this, option)
        class(UKRMOLIntegral)      :: this
        class(Options), intent(in) :: option

        write (stdout, "('Using UKRMOL+ integrals')")

        this % num_csf = option % num_csfs
        this % num_symmetry = option % num_syms

        allocate(this % num_orbitals_sym(this % num_symmetry))

        this % num_orbitals_sym(:) = option % num_electronic_orbitals_sym(:)
        this % hamiltonian_unit = option % hamiltonian_unit
        this % exotic_mass = option % exotic_mass
        this % matrix_io_unit = option % hamiltonian_unit
        this % name = option % name

    end subroutine initialize_UKRMOL


    !> \brief   Finalize the type
    !> \authors J Benda
    !> \date    2019
    !>
    !> Deallocate memory used by the integral storage.
    !>
    subroutine finalize_UKRMOL (this)

        use ukrmol_interface_gbl, only: destroy_ukrmolp

        class(UKRMOLIntegral) :: this

        if (allocated(this % num_orbitals_sym)) deallocate (this % num_orbitals_sym)
        if (allocated(this % xnuc))   deallocate (this % xnuc)
        if (allocated(this % ynuc))   deallocate (this % ynuc)
        if (allocated(this % znuc))   deallocate (this % znuc)
        if (allocated(this % charge)) deallocate (this % charge)
        if (allocated(this % cname))  deallocate (this % cname)

        call destroy_ukrmolp

    end subroutine finalize_UKRMOL


    !> \brief   Read UKRmol+ integrals
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Uses calls to GBTOlib to perform reading itself.
    !>
    subroutine load_UKRMOL (this, iounit)
        use ukrmol_interface_gbl, only: read_ukrmolp_ints, get_geom

        class(UKRMOLIntegral) :: this
        integer, intent(in)   :: iounit
        integer :: nalm = 0, lembf, matrix_size = 0,dummy_int

        this % num_one_electron_integrals = -1
        this % num_two_electron_integrals = -1

        call read_ukrmolp_ints(this % matrix_io_unit, iounit, this % number_of_matrix_records, &
                               this % num_one_electron_integrals, this % num_two_electron_integrals, &
                               this % num_csf, stdout, this % symmetry_type, this % num_symmetry, this % num_orbitals_sym, &
                               this % positron_flag, 1.0_wp, this % name, &
                               nalm, .false.)

        !Get geometries
        call get_geom(this % nnuc, this % cname, this % xnuc, this % ynuc, this % znuc, this % charge)

        !Read header files locally to store to memory
        rewind (this % matrix_io_unit)
        read (this % matrix_io_unit) matrix_size, this % number_of_matrix_records , dummy_int, matrix_size, &
                                     dummy_int, this % num_symmetries, dummy_int, dummy_int, dummy_int, dummy_int, this % nnuc, &
                                     dummy_int, this % NAME, this % NHE, this % DTNUC
        this % core_energy = this % DTNUC(1)
        close (this % matrix_io_unit)

    end subroutine load_UKRMOL


    !> \brief   Get integral value
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Obtain the integral value from GBTOlib.
    !>
    function get_integral_UKRMOL (this, i, j, k, l, m) result(coeff)
        use ukrmol_interface_gbl, only: get_integral

        class(UKRMOLIntegral) :: this
        integer, intent(in)   :: i, j, k, l, m
        integer  :: a, b, c, d
        real(wp) :: coeff

        if (i == 0) then
            a = this % orbital_mapping(j)
            b = this % orbital_mapping(l)
            coeff = get_integral(a, b, 0, 0, m)
        else
            a = this % orbital_mapping(i)
            b = this % orbital_mapping(j)
            c = this % orbital_mapping(k)
            d = this % orbital_mapping(l)
            coeff = get_integral(a, b, c, d, m)
        end if

    end function get_integral_UKRMOL


    !> \brief   Write geometries to file
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine write_geometries_UKRMOL (this, iounit)

        class(UKRMOLIntegral), intent(in) :: this
        integer,               intent(in) :: iounit
        integer :: ido, ifail, i

        do ido = 1, this % nnuc
            write (iounit) this % cname(ido), this % xnuc(ido), this % ynuc(ido), this % znuc(ido), this % charge(ido)
        end do

        write (stdout, 120) (this % cname(i), this % xnuc(i), this % ynuc(i), this % znuc(i), this % charge(i), i = 1, this % nnuc)
        120 FORMAT(/,' Nuclear data     X         Y         Z       Charge',/,(3x,a8,2x,4F10.6))

    end subroutine write_geometries_UKRMOL


    !> \brief   Destroy this type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine destroy_integral_UKRMOL (this)
        use ukrmol_interface_gbl, only: destroy_ukrmolp

        class(UKRMOLIntegral) :: this

        !Where is the destory integral or cleanup function?
        !MADE ONE!
        call destroy_ukrmolp()

    end subroutine destroy_integral_UKRMOL

end module UKRMOL_module
