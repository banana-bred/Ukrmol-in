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
!> This module handles construction of the spin orbital table and of generating symbolic elements from determinants
!>
!> \note 01/02/2017 - Ahmed Al-Refaie: Initial revision.
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Orbital_module

    use const_gbl,       only: stdout
    use consts_mpi_ci,   only: SYMTYPE_CINFV, SYMTYPE_DINFH, SYMTYPE_D2H
    use mpi_gbl,         only: mpi_xermsg
    use precisn,         only: wp
    use scatci_routines, only: MKORBS, PMKORBS
    use Symbolic_Module, only: SymbolicElementVector

    implicit none

    private

    public SpinOrbital, OrbitalTable

    !> \brief   This type holds a single spin orbital
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Basic building block of \ref OrbitalTable.
    !>
    type :: SpinOrbital
        integer :: orbital_number   !< MN
        integer :: gerude           !< MG
        integer :: m_quanta         !< MM
        integer :: spin             !< MS
        integer :: positron_flag    !< MPOS
        integer :: electron_idx     !< Electron number for the references
    end type SpinOrbital


    !> \brief   This class generates the molecular and spin orbitals, stores them and generates symblic elements from determinants
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Consists of an array of objects \ref SpinOrbital. Provides access functions for retrieval of properties
    !> of individual orbitals.
    !>
    type :: OrbitalTable
        logical :: initialized = .false.    !< Whether the class has been initialized
        integer :: total_num_spin_orbitals  !< The total number of spin orbitals
        integer :: total_num_orbitals       !< The total number of molecular orbitals
        integer :: symmetry_type            !< Which symmetry group we are dealing with
        integer :: num_symmetries           !< How many symmetries we have
        integer :: positron_flag            !< Whether we have exotic matter
        integer :: MFLG = 0
        integer,           allocatable :: mcorb(:), mcon(:), nsrb(:)  !< MCORB - Orbital mapping?
        type(SpinOrbital), allocatable :: spin_orbitals(:)            !< Our table of spin orbitals
        integer,           allocatable :: orbital_map(:)              !< Mapping of determinants to spin orbitals
    contains
        procedure, public    ::    initialize => initialize_table
        procedure, public    ::    construct => construct_table
        procedure, public    ::    get_orbital_number
        procedure, public    ::    get_spin
        procedure, public    ::    get_gerude
        procedure, public    ::    get_electron_number
        procedure, public    ::    compute_electron_index
        procedure, public    ::    check_max_mcon_in_determinants
        procedure, public    ::    evaluate_IJKL_and_coeffs
        procedure, public    ::    get_mcon
        procedure, public    ::    get_minimum_mcon
        procedure, public    ::    get_two_minimum_mcon
        procedure, public    ::    add_positron
        procedure, private   ::    evaluate_case_one
        procedure, private   ::    evaluate_case_two
        procedure, private   ::    evaluate_case_three
        procedure, private   ::    evaluate_case_other
    end type OrbitalTable

contains

    !> \brief   Basic initialization of the data structure
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Stores given arguments and allocates memory for all spin-orbitals.
    !>
    !> \param[inout] this   Orbital object to update.
    !> \param[in] nsrb      Number of spin-orbitals.
    !> \param[in] norb      Number of orbitals.
    !> \param[in] symtype   Symmetry type (group class). Expects named constant from \ref consts_mpi_ci.
    !> \param[in] nsym      Number of symmetries (irreducible representations).
    !> \param[in] positron_flag  Non-zero value indicates positron.
    !>
    subroutine initialize_table (this, nsrb, norb, symtype, nsym, positron_flag)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: nsrb, nsym, norb, symtype, positron_flag

        this % total_num_spin_orbitals = nsrb
        this % total_num_orbitals      = norb
        this % symmetry_type           = symtype
        this % num_symmetries          = nsym
        this % positron_flag           = positron_flag

        if (allocated(this % spin_orbitals))    deallocate (this % spin_orbitals)
        if (allocated(this % mcon))             deallocate (this % mcon)
        if (allocated(this % nsrb))             deallocate (this % nsrb)
        if (allocated(this % mcorb))            deallocate (this % mcorb)
        if (allocated(this % orbital_map))      deallocate (this % orbital_map)

        allocate(this % spin_orbitals(this % total_num_spin_orbitals))
        allocate(this % mcon(this % total_num_spin_orbitals))
        allocate(this % nsrb(this % total_num_spin_orbitals))
        allocate(this % mcorb(this % total_num_orbitals))
        allocate(this % orbital_map(this % total_num_orbitals))

        this % mcorb(:) = 0

    end subroutine initialize_table


    !> \brief   Define all spin-orbitals
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Assumes that \ref initialize_table has been called before.
    !>
    !> Sets all attributes of the individual spin-orbital data structures. This is actually a simple
    !> wrapper around the original SCATCI subroutines MKORBS and PMKORBS,
    !> where the construction of the spin-orbitals is carried out.
    !>
    !> On exit from this subroutine, all spin-orbitals contained in this type have correctly defined attributes
    !> (quantum numbers).
    !>
    subroutine construct_table (this, num_orbital_target_sym, num_orbital_target_sym_dinf, &
                                num_orbitals, num_elec_orbitals, num_orbitals_congen)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: num_orbital_target_sym(this % num_symmetries),        &
                                              num_orbital_target_sym_dinf(this % num_symmetries),   &
                                              num_orbitals(this % num_symmetries),                  &
                                              num_elec_orbitals(this % num_symmetries),             &
                                              num_orbitals_congen(this % num_symmetries)

        integer :: MN(this % total_num_spin_orbitals), MG(this % total_num_spin_orbitals), &
                   MM(this % total_num_spin_orbitals), MS(this % total_num_spin_orbitals), MPOS(this % total_num_spin_orbitals)
        integer :: dummy_iposit = 0, dummy_lusme = 0, ido

        if (this % positron_flag /= 0 .and. this % symmetry_type == SYMTYPE_D2H) then
            call PMKORBS(num_orbitals, num_elec_orbitals, num_orbitals_congen, this % num_symmetries,    &
                         MN, MG, MM, MS, this % mcon, this % mcorb, this % total_num_orbitals,           &
                         this % total_num_spin_orbitals, this % orbital_map, MPOS, this % positron_flag, &
                         this % symmetry_type, dummy_lusme, stdout)
        else
            call MKORBS(this % num_symmetries, MN, MG, MM, MS, this % total_num_orbitals,      &
                        this % total_num_spin_orbitals, this % orbital_map, MPOS, this % mcon, &
                        this % mcorb, this % positron_flag, num_orbital_target_sym,            &
                        num_orbital_target_sym_dinf, this % symmetry_type, dummy_lusme, stdout)
        end if

        do ido = 1, this % total_num_spin_orbitals
            this % spin_orbitals(ido) % orbital_number = MN(ido)
            this % spin_orbitals(ido) % gerude         = MG(ido)
            this % spin_orbitals(ido) % m_quanta       = MM(ido)
            this % spin_orbitals(ido) % spin           = MS(ido)
            this % spin_orbitals(ido) % positron_flag  = MPOS(ido)
            this % spin_orbitals(ido) % electron_idx   = 0
        end do

    end subroutine construct_table


    !> \brief   Assign electrons to reference spin-orbitals
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Assign consecutive electron index to all stored spin-orbitals of the reference determinant.
    !>
    subroutine compute_electron_index (this, num_electrons, reference_determinants)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: num_electrons
        integer,             intent(in)    :: reference_determinants(num_electrons)

        integer :: ido

        do ido = 1, num_electrons
            this % spin_orbitals(reference_determinants(ido)) % electron_idx = ido
        end do

    end subroutine compute_electron_index


    integer function check_max_mcon_in_determinants (this, n, determinants)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: determinants(:), n

        integer :: i

        check_max_mcon_in_determinants = 0
        do i = 1, n
            check_max_mcon_in_determinants = max(check_max_mcon_in_determinants, this % mcon(determinants(i)))
        end do

    end function check_max_mcon_in_determinants


    !> \brief   This compares the determinants and generates the proper coefficents and ijklm values needed for the symbols
    !>
    subroutine evaluate_IJKL_and_coeffs (this, dtrs, coeff, symmetry_type, symbol, flag)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: dtrs(4), symmetry_type, flag
        real(wp),            intent(in)    :: coeff
        class(SymbolicElementVector)       :: symbol

        type(SpinOrbital) :: P,R            !Determinants A
        type(SpinOrbital) :: Q,S            !Determinants B
        integer  :: KA, KB, lpositron = 0
        integer  :: NSFA, NSFB, MLA, MLB    !I'm not 100% what these do alpha/beta?
        real(wp) :: SIGN

        if (this % positron_flag /= 0) lpositron = this % add_positron(dtrs, 1, 4)

        !Get our spin orbitals
        P = this % spin_orbitals(dtrs(1))
        Q = this % spin_orbitals(dtrs(2))
        R = this % spin_orbitals(dtrs(3))
        S = this % spin_orbitals(dtrs(4))

        !If there are opposing spins
        if(P % spin + Q % spin == 1 .or. R % spin + S % spin == 1 .or. P % positron_flag /= Q % positron_flag) then
            NSFA = 0
        else
            NSFA = 1
            MLA = 0
            if (symmetry_type < SYMTYPE_D2H) then
                KA = MAX(dtrs(1), dtrs(2))
                KB = MAX(dtrs(3), dtrs(4))
                if (KA < KB) then
                    if (R % m_quanta * S % m_quanta < 0) mla = 1
                else
                    if(P % m_quanta * Q % m_quanta < 0) mla = 1
                end if
            end if
        end if

        !If there are opposing spins or if there is a single positron
        if (P % spin + S % spin == 1 .or. Q % spin + R % spin == 1 .or. P % positron_flag /= S % positron_flag) then
            NSFB = 0
        else
            NSFB = 1
            MLB = 0
            if (symmetry_type < SYMTYPE_D2H) then
                KA = MAX(dtrs(1), dtrs(4))
                KB = MAX(dtrs(3), dtrs(2))
                if (KA < KB) then
                    if (Q % m_quanta * R % m_quanta < 0) mlb = 1
                else
                    if (P % m_quanta * S % m_quanta < 0) mlb = 1
                end if
            end if
        end if

        !No exchange so ignore
        if (NSFA == 0 .and. NSFB == 0) return

        SIGN = 1_wp
        if (lpositron == 1) SIGN = -1_wp

        if (flag == 1) then
            if (nsfa /= 0) then
                call this % evaluate_case_other(P % orbital_number, Q % orbital_number, R % orbital_number, &
                                                S % orbital_number, mla, 0, SIGN * coeff, symbol)
            end if
            if (nsfb /= 0) then
                call this % evaluate_case_other(P % orbital_number, S % orbital_number, R % orbital_number, &
                                                Q % orbital_number, mlb, 1, SIGN * coeff, symbol)
            end if
        else if (P % orbital_number == Q % orbital_number .and. R % orbital_number == S % orbital_number) then

            call this % evaluate_case_one(P % orbital_number, R % orbital_number, nsfa, nsfb, mla, mlb, SIGN * coeff, symbol)

        else if (P % orbital_number == Q % orbital_number .or. R % orbital_number == S % orbital_number) then
            if (nsfa /= 0) then
                call this % evaluate_case_other(P % orbital_number, Q % orbital_number, R % orbital_number, &
                                                S % orbital_number, mla, 0, SIGN * coeff, symbol)
            end if
            if (nsfb /= 0) then
                call this % evaluate_case_other(P % orbital_number, S % orbital_number, R % orbital_number, &
                                                Q % orbital_number, mlb, 1, SIGN * coeff, symbol)
            end if
        else if (P % orbital_number == R % orbital_number .and. Q % orbital_number == S % orbital_number) then

            call this % evaluate_case_three(P % orbital_number, R % orbital_number, nsfa, nsfb, mla, mlb, SIGN * coeff, symbol)

        else if (P % orbital_number == R % orbital_number) then
            if (nsfa /= 0) then
                call this % evaluate_case_other(P % orbital_number, Q % orbital_number, R % orbital_number, &
                                                S % orbital_number, mla, 0, SIGN * coeff, symbol)
            end if
            if (nsfb /= 0) then
                call this % evaluate_case_other(P % orbital_number, S % orbital_number, R % orbital_number, &
                                                Q % orbital_number, mlb, 1, SIGN * coeff, symbol)
            end if
        else if(P % orbital_number == S % orbital_number .and. Q % orbital_number == R % orbital_number) then

            call this % evaluate_case_two(P % orbital_number, R % orbital_number, nsfa, nsfb, mla, mlb, SIGN * coeff, symbol)

        else
            if (nsfa /= 0) then
                call this % evaluate_case_other(P % orbital_number, Q % orbital_number, R % orbital_number, &
                                                S % orbital_number, mla, 0, SIGN * coeff, symbol)
            end if
            if (nsfb /= 0) then
                call this % evaluate_case_other(P % orbital_number, S % orbital_number, R % orbital_number, &
                                                Q % orbital_number, mlb, 1, SIGN * coeff, symbol)
            end if
        end if

    end subroutine evaluate_IJKL_and_coeffs


    subroutine evaluate_case_one (this, p, r, nsfa, nsfb, mla, mlb, coeff, symbol)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: p, r, NSFA, NSFB, MLA, MLB
        real(wp),            intent(in)    :: coeff
        class(SymbolicElementVector)       :: symbol

        integer :: I, J

        if (P < R) then
            I = R
            J = P
        else
            I = P
            J = R
        end if

        if (nsfa /= 0) then
            if (mla == 0) then
                call symbol % insert_ijklm_symbol(I, I, J, J, 1, coeff)
            else
                call symbol % insert_ijklm_symbol(I, I, J, J, 0, coeff)
            end if
        end if

        if (nsfb /= 0) then
            if (mlb == 0) then
                call symbol % insert_ijklm_symbol(I, J, I, J, 1, -coeff)
            else
                call symbol % insert_ijklm_symbol(I, J, I, J, 0, -coeff)
            end if
        end if

    end subroutine evaluate_case_one


    subroutine evaluate_case_two (this, p, q, nsfa, nsfb, mla, mlb, coeff, symbol)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: p, q, NSFA, NSFB, MLA, MLB
        real(wp),            intent(in)    :: coeff
        class(SymbolicElementVector)       :: symbol

        integer :: I, J

        if (P < Q) then
            I = Q
            J = P
        else
            I = P
            J = Q
        end if

        if (nsfa /= 0) then
            if (mla == 0) then
                call symbol % insert_ijklm_symbol(I, J, I, J, 1, coeff)
            else
                call symbol % insert_ijklm_symbol(I, J, I, J, 0, coeff)
            end if
        end if

        if (nsfb /= 0) then
            if (mlb == 0) then
                call symbol % insert_ijklm_symbol(I, I, J, J, 1, -coeff)
            else
                call symbol % insert_ijklm_symbol(I, I, J, J, 0, -coeff)
            end if
        end if

    end subroutine evaluate_case_two


    subroutine evaluate_case_three (this, p, q, nsfa, nsfb, mla, mlb, coeff, symbol)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: p, q, NSFA, NSFB, MLA, MLB
        real(wp),            intent(in)    :: coeff
        class(SymbolicElementVector)       :: symbol

        integer :: I, J

        if (P < Q) then
            I = Q
            J = P
        else
            I = P
            J = Q
        end if

        if (nsfa /= 0) then
            if (mla == 0) then
                call symbol % insert_ijklm_symbol(I, J, I, J, 1, coeff)
            else
                call symbol % insert_ijklm_symbol(I, J, I, J, 0, coeff)
            end if
        end if

        if (nsfb /= 0) then
            if (mlb == 0) then
                call symbol % insert_ijklm_symbol(I, J, I, J, 1, -coeff)
            else
                call symbol % insert_ijklm_symbol(I, J, I, J, 0, -coeff)
            end if
        end if

    end subroutine evaluate_case_three


    subroutine evaluate_case_other (this, p, q, r, s, mla, ms, coeff, symbol)
        class(OrbitalTable), intent(inout) :: this
        integer,             intent(in)    :: p, q, r, s, mla, ms
        real(wp),            intent(in)    :: coeff
        class(SymbolicElementVector)       :: symbol

        real(wp) :: cfd
        integer  :: I, J, K, L

        if (P < Q) then
            I = Q
            J = P
        else
            I = P
            J = Q
        end if
        if (R < S) then
            K = S
            L = R
        else
            K = R
            L = S
        end if

        cfd = coeff
        if (ms /= 0) cfd = -cfd

        if (I > K .or. (I == K .and. J > L)) then
            call symbol % insert_ijklm_symbol(I, J, K, L, mla, cfd)
        else
            call symbol % insert_ijklm_symbol(K, L, I, J, mla, cfd)
        end if

    end subroutine evaluate_case_other


    !> Simple function to get orbital number for a specific spin orbital
    integer function get_orbital_number (this, spin_orbital)

        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: spin_orbital

        if (spin_orbital > this % total_num_spin_orbitals .or. spin_orbital < 1) then
            write (stdout, "('Orbital selected: ',2i12)") spin_orbital, this % total_num_spin_orbitals
            flush(stdout)
            call mpi_xermsg('Orbital_module', 'get_orbital_number', 'Selected orbital not within range', 1, 1)
        end if

        get_orbital_number = this % spin_orbitals(spin_orbital) % orbital_number

    end function get_orbital_number


    !> Simple function to get orbital number for a specific spin orbital
    integer function get_spin (this, spin_orbital)
        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: spin_orbital

        if (spin_orbital > this % total_num_spin_orbitals .or. spin_orbital < 1) then
            write (stdout, "('Orbital selected: ',2i12)") spin_orbital, this % total_num_spin_orbitals
            flush(stdout)
            call mpi_xermsg('Orbital_module', 'get_spin', 'Selected orbital not within range', 1, 1)
        end if

        get_spin = this % spin_orbitals(spin_orbital) % spin

    end function get_spin


    integer function get_gerude (this, spin_orbital)
        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: spin_orbital

        if (spin_orbital > this % total_num_spin_orbitals .or. spin_orbital < 1) then
            write (stdout, "('Orbital selected: ',2i12)") spin_orbital, this % total_num_spin_orbitals
            flush(stdout)
            call mpi_xermsg('Orbital_module', 'get_gerude', 'Selected orbital not within range', 1, 1)
        end if

        get_gerude = this % spin_orbitals(spin_orbital) % gerude

    end function get_gerude


    integer function get_electron_number (this, spin_orbital)
        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: spin_orbital

        if (spin_orbital > this % total_num_spin_orbitals .or. spin_orbital < 1) then
            write (stdout, "('Orbital selected: ',2i12)") spin_orbital, this % total_num_spin_orbitals
            flush(stdout)
            call mpi_xermsg('Orbital_module', 'get_electron_number', 'Selected orbital not within range', 1, 1)
        end if

        get_electron_number = this % spin_orbitals(spin_orbital) % electron_idx

    end function get_electron_number


    integer function get_minimum_mcon (this, determinants)
        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: determinants(4)

        !get_minimum_mcon = min(this%get_two_minimum_mcon(determinants(1),determinants(2)),  &
        !           &    this%get_two_minimum_mcon(determinants(3),determinants(4)))

        get_minimum_mcon = min(this % mcon(determinants(1)), &
                               this % mcon(determinants(2)), &
                               this % mcon(determinants(3)), &
                               this % mcon(determinants(4)))

    end function get_minimum_mcon


    integer function get_mcon (this, spin_orbital)
        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: spin_orbital

        if (spin_orbital > this % total_num_spin_orbitals .or. spin_orbital < 1) then
            write (stdout, "('Orbital selected: ',2i12)") spin_orbital, this % total_num_spin_orbitals
            flush(stdout)
            call mpi_xermsg('Orbital_module', 'get_mcon', 'Selected orbital not within range', 1, 1)
        end if

        get_mcon = this % mcon(spin_orbital)

    end function get_mcon


    integer function get_two_minimum_mcon (this, determinant_one, determinant_two)
        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: determinant_one, determinant_two

        get_two_minimum_mcon = min(this % mcon(determinant_one), this % mcon(determinant_two))

    end function get_two_minimum_mcon


    integer function add_positron (this, determinants, ia, ib)
        class(OrbitalTable), intent(in) :: this
        integer,             intent(in) :: determinants(4), ia, ib
        integer :: icount, ido, ip

        add_positron = 0
        icount = 0

        do ido = ia, ib
            ip = determinants(ido)
            if (this % spin_orbitals(ip) % positron_flag /= 0) icount = icount + 1
        end do

        if (icount == 0) return
        if (icount == 2) then
            add_positron = 1
            return
        end if

        write (stdout, "(' ICOUNT ',i4,' NDTC = ',4i4)") icount, determinants(:)
        flush(stdout)
        call mpi_xermsg('Orbital_module', 'add_positron', 'Invalid count', 1, 1)

    end function add_positron

end module Orbital_module
