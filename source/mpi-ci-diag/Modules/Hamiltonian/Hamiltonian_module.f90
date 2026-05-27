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

!> \brief   Base Hamiltonian module
!> \authors A Al-Refaie
!> \date    2017
!>
!> Provides an abstract hamiltonian class that has built in functionalities for all inherited hamiltonians.
!>
!> \todo Rewrite the idiag = 0 matrix evaluation into a more easily handled form and chenge variable naming to be less confusing.
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial documentation version
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Hamiltonian_module

    use precisn,                only: longint, wp
    use BaseIntegral_module,    only: BaseIntegral
    use consts_mpi_ci,          only: NIDX
    use CSF_module,             only: CSFObject, CSFOrbital
    use Options_module,         only: Options
    use Orbital_module,         only: OrbitalTable, SpinOrbital
    use Parallelization_module, only: grid => process_grid
    use Symbolic_module,        only: SymbolicElementVector
    use Utility_module,         only: unpack_ints

    implicit none

    public BaseHamiltonian

    !> \brief   This is an abstract class that contains the majority of functionality required to construct hamiltonians
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This class provides an abstraction of a lot of the components required to build the hamiltonian.
    !> For example, the user does not need to worry about the specific implementation of Slaters rules or evaluating integrals
    !> this class wraps these features for you and should allow one to implement hamiltonians closer to those described in papers.
    !> This is not a Matrix class. BaseHamiltonain deals with building the matrix whilst BaseMatrix deals with storing the matrix.
    !>
    type, abstract :: BaseHamiltonian
        class(OrbitalTable), pointer :: orbitals  !< Our orbitals required to generate symblic elements
        class(Options),      pointer :: options   !< Scatci program settings
        class(BaseIntegral), pointer :: integral  !< The integrals we are using
        class(CSFObject),    pointer :: csfs(:)   !< Our configuration state functions

        integer  :: NFLG = 0
        integer  :: diagonal_flag
        integer  :: positron_flag           !< Positron aware flag
        integer  :: phase_flag              !< whether to evaluate integrals whilst dealing with phase
        logical  :: constructed = .false.   !< Has the hamiltonain been constructed
        logical  :: initialized = .false.   !< Has the hamiltonian been initialized
        integer  :: job_id = 0              !< Whose job it is to (soon to be deprecated)
        integer  :: number_of_integrals = 0 !< How many integrals have been evaluated?
        real(wp) :: element_one = 0.0       !< First element for idiag = 0
        real(wp) :: enh_factor              !< Enhancement factor set in options
        integer  :: enh_factor_type         !< Enhancement factor type set in options

        type(SymbolicElementVector) :: reference_symbol   !< Symbols for idiag = 0
    contains
        !Constructor
        procedure, public  :: construct => construct_base_hamiltonian
        procedure(generic_build), deferred :: build_hamiltonian

        !-----------------------Private procedures-----------------!
        procedure, public  :: slater_rules
       !procedure, public  :: obey_slater_rules
        procedure, public  :: evaluate_integrals
        procedure, public  :: evaluate_integrals_singular
       !procedure, public  :: store_first_element

        !-------Slater rule functions------------------!
        procedure, private :: handle_two_pair
        procedure, private :: handle_one_pair
        procedure, private :: handle_same
        procedure, public  :: my_job
    end type BaseHamiltonian

    abstract interface
        !> \brief   Main build routine of the hamiltonian
        !> \authors A Al-Refaie
        !> \date    2017
        !>
        !> All build must be done within this routine in order to be used by MPI-SCATCI.
        !> \param[out] matrix_elements  Resulting matrix elements from the build.
        !>
        subroutine generic_build(this,matrix_elements)
            use BaseMatrix_module, only: BaseMatrix
            import :: BaseHamiltonian
            class(BaseHamiltonian)           :: this
            class(BaseMatrix), intent(inout) :: matrix_elements
        end subroutine generic_build
    end interface

contains

    !> \brief   Constructs the hamiltonain object and provides all the various componenets required for building
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this   Hamiltonian object to update.
    !> \param[in] option    An Options class containing run constants
    !> \param[in] csfs      An array of configuration state functions
    !> \param[in] orbitals  An initialized Orbitals class
    !> \param[in] integral  A BaseIntegral class for evaluating symbols
    !>
    subroutine construct_base_hamiltonian (this, option, csfs, orbitals, integral)
        class(BaseHamiltonian),      intent(inout) :: this
        class(Options),      target, intent(in)    :: option
        class(CSFObject),    target, intent(in)    :: csfs(:)
        class(OrbitalTable), target, intent(in)    :: orbitals
        class(BaseIntegral), target, intent(in)    :: integral

        !Assign our member pointers for each class
        this % options  => option
        this % csfs     => csfs
        this % orbitals => orbitals
        this % integral => integral
        this % phase_flag = this % options % phase_correction_flag
        this % positron_flag = this % options % positron_flag
        this % enh_factor = this % options % enh_factor
        this % enh_factor_type = this % options % enh_factor_type

        this % job_id = 0
        this % diagonal_flag = min(this % options % matrix_eval, 2)
        this % constructed = .true.

    end subroutine construct_base_hamiltonian


    !> \brief   Performs Slater rules and returns symbolic elements
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Computes the Slater rules and adds them to the symbolic elements class that is provided.
    !> Additionaly has an automatic MPI mode by usage of the job_ = .true. keyword.
    !> Whilst not faster than manually unrolling a loop, it provides a quick and easy way to fairly
    !> effectively parallelize the loop. However this is not OpenMP compatible
    !>
    !> \param[inout] this       Hamiltonian object to query.
    !> \param[in]  csf_a        First configuration state function to slater
    !> \param[in]  csf_b        Second configuration state function to slater
    !> \param[out] sym_elements Resulting smbolic elements
    !> \param[in]  flag         Whether we are doing the diagonal or not
    !> \param[in]  job_         Whether to use the 'easy' MPI split method
    !>
    subroutine slater_rules (this, csf_a, csf_b, sym_elements, flag, job_)
        class(BaseHamiltonian),       intent(in)    :: this
        class(CSFObject),             intent(in)    :: csf_a, csf_b
        class(SymbolicElementVector), intent(inout) :: sym_elements
        integer,           intent(in)  :: flag
        logical, optional, intent(in)  :: job_
        logical  :: job
        integer  :: diag, num_differences, dtr_diff(4), dtr_diff_fast(4), dtr_idx_a, dtr_idx_b
        real(wp) :: coeff

        !Check for job parameter
        if (present(job_)) then
            job = job_
        else
            job = .true.
        end if

        !Clear the symbolic elements ahead of time
        call sym_elements % clear()

        !Check if we should be doing this
        if (job) then
            !If not then return
            if (.not. this % my_job()) return
        end if

        !Check if we are not the same
        if (csf_a % id /= csf_b % id) then
            !If we are not then quickly check if we need to do a slater calculation
            if (csf_a % check_slater(csf_b) > 4) return
            !if(this%obey_slater_rules(csf_a,csf_b) > 4) return
        end if

        !Loop through the first set of determinants
        do dtr_idx_a = 1, csf_a % num_orbitals
            !Loop through the second set of determinants
            do dtr_idx_b = 1, csf_b % num_orbitals

                !Use the fast slater rules method to check the number of differences and get the differing determinants
                call csf_a % orbital(dtr_idx_a) % compare_excitations_fast(csf_b % orbital(dtr_idx_b), &
                                                                           this % options % num_electrons, &
                                                                           coeff, &
                                                                           num_differences, &
                                                                           dtr_diff)

                !Depending on the number of differences, we handle it appropriately
                select case (num_differences)
                    case(3:)
                        cycle
                    case(2)
                        call this % handle_two_pair(dtr_diff, coeff, sym_elements, flag)
                    case(1)
                        call this % handle_one_pair(dtr_diff, coeff, sym_elements, csf_a % orbital(dtr_idx_a), flag)
                    case(0)
                        call this % handle_same(dtr_diff, coeff, sym_elements, csf_a % orbital(dtr_idx_a), flag)
                end select

            end do
        end do

    end subroutine slater_rules


    !> \brief   Handles if there are two pairs of spin orbitals that differ
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this       Hamiltonian object to query.
    !> \param[in]  dtrs         A 4 element array of slater determinants that differ.
    !> \param[in]  coeff        The computed coefficients from the Slater rules.
    !> \param[out] sym_element  Resulting symbolic elements.
    !> \param[in]  flag         Whether we are doing the diagonal or not.
    !>
    subroutine handle_two_pair (this, dtrs, coeff, sym_element, flag)
        class(BaseHamiltonian)       ::    this
        class(SymbolicElementVector) ::    sym_element
        integer,  intent(in)         ::    dtrs(4)
        real(wp), intent(in)         ::    coeff
        integer,  intent(in)         :: flag

        if (this % diagonal_flag> 1 .and. this % orbitals % get_minimum_mcon(dtrs) > 0) return

        !Evaluate for symbolic elements
        call this % orbitals % evaluate_IJKL_and_coeffs(dtrs, coeff, 0, sym_element, flag)

    end subroutine handle_two_pair


    !> \brief   Handles if there are one pair of spin orbitals that differ
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this       Hamiltonian object to query.
    !> \param[in]  dtrs         A 4 element array of slater determinants
    !> \param[in]  coeff        The computed coefficients from the Slater rules
    !> \param[in]  csf          The first CSFs determinant
    !> \param[out] sym_element  Resulting symbolic elements
    !> \param[in]  flag         Whether we are doing the diagonal or not
    !>
    subroutine handle_one_pair (this, dtrs, coeff, sym_element, csf, flag)
        class(BaseHamiltonian)           :: this
        class(SymbolicElementVector)     :: sym_element
        class(CSFOrbital), intent(in)    :: csf
        integer,           intent(inout) :: dtrs(4)
        real(wp),          intent(in)    :: coeff
        integer,           intent(in)    :: flag

        type(SpinOrbital) :: so1, so2

        integer :: ref_dtrs(this % options % num_electrons)
        integer :: num_electrons, N, M, ido, ia1, ia2, lnadd
        integer :: nchk !Whether its continuum maybe?

        !Reset the positron flag
        lnadd = 0

        !Get number of electrons
        num_electrons = this % options % num_electrons

        !Copt the reference determinant to a seperate array
        !ref_dtrs(:) = this%options%reference_dtrs(:)

        !Loop through the reference and replace the spin orbitals with ones from the CSF
        !do ido=1,csf%num_replaced_dtrs
        !    N = csf%replacing(ido)
        !    M = this%orbitals%get_electron_number(N)
        !    ref_dtrs(M) = csf%replaced(ido)
        !enddo

        call csf % get_determinants(ref_dtrs, num_electrons)

        !Not sure what this does exactly
        if (this % diagonal_flag > 1 .and. this % orbitals % get_two_minimum_mcon(dtrs(3), dtrs(4)) > 0) then
            nchk = 1
        else
            nchk = 0
        end if

        !Loop through all the spin orbitals
        do ido = 1, num_electrons
            !If we have the same spin orbitals then skip
            if (ref_dtrs(ido) == dtrs(3)) cycle
            !Again not quite sure what this does
            if (nchk == 1 .and. this % orbitals % get_mcon(ref_dtrs(ido)) > 0) cycle
            !Place the reference into our determinant array
            dtrs(1) = ref_dtrs(ido)
            dtrs(2) = ref_dtrs(ido)
            !Now evaluate for symbols
            call this % orbitals % evaluate_IJKL_and_coeffs(dtrs, coeff, 0, sym_element, flag)
        end do

        !If this run contains positrons, then add them in if necessary
        if (this % positron_flag /= 0) lnadd = this % orbitals % add_positron(dtrs, 3, 4)
        if (nchk == 1) return

        !Get our original slater replacements
        so1 = this % orbitals % spin_orbitals(dtrs(3))
        so2 = this % orbitals % spin_orbitals(dtrs(4))

        !Check if they posses the same spin and m quanta
        if (so1 % m_quanta /= so2 % m_quanta .or. so1 % spin /= so2 % spin) return
        if (so1 % orbital_number < so2 % orbital_number) then
            ia1 = so2 % orbital_number
            ia2 = so1 % orbital_number
        else
            ia1 = so1 % orbital_number
            ia2 = so2 % orbital_number
        end if

        !if they do then add a 0I0J + P symbol
        call sym_element % insert_ijklm_symbol(0, ia1, 0, ia2, lnadd, coeff)

    end subroutine handle_one_pair


    !> \brief   Handles if there are no differences between spin orbitals
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this       Hamiltonian object to query.
    !> \param[in]  dtrs         A 4 element array of slater determinants
    !> \param[in]  coeff        The computed coefficients from the Slater rules
    !> \param[in]  csf          The first CSFs determinant
    !> \param[out] sym_element  Resulting symbolic elements
    !> \param[in]  flag         Whether we are doing the diagonal or not
    !>
    subroutine handle_same (this, dtrs, coeff, sym_element, csf, flag)
        class(BaseHamiltonian)           :: this
        class(SymbolicElementVector)     :: sym_element
        class(CSFOrbital), intent(in)    :: csf
        integer,           intent(inout) :: dtrs(4)
        integer                          :: pdtrs(4)
        real(wp),          intent(in)    :: coeff
        integer,           intent(in)    :: flag

        type(SpinOrbital) :: so1, so2

        integer :: ref_dtrs(this % options % num_electrons)
        integer :: num_electrons, N, M, ido, jdo, ia1, ia2, lnadd
        integer :: nchk   !Whether its continuum maybe?

        !Reset the positron flag
        lnadd = 0

        !Get number of electrons
        num_electrons = this % options % num_electrons

        !Copt the reference determinant to a seperate array
        !ref_dtrs(:) = this%options%reference_dtrs(:)

        !Loop through the reference and replace the spin orbitals with ones from the CSF
        !do ido=1,csf%num_replaced_dtrs
        !    N = csf%replacing(ido)
        !    M = this%orbitals%get_electron_number(N)
        !    ref_dtrs(M) = csf%replaced(ido)
        !enddo

        call csf % get_determinants(ref_dtrs, num_electrons)

        do ido = 2, num_electrons
            dtrs(1) = ref_dtrs(ido)
            dtrs(2) = dtrs(1)
            do jdo = 1, ido - 1
                dtrs(3) = ref_dtrs(jdo)
                dtrs(4) = dtrs(3)
                if (this % diagonal_flag > 1 .and. this % orbitals % get_two_minimum_mcon(dtrs(1), dtrs(3)) > 0) cycle
                !Now evaluate for symbols
                call this % orbitals % evaluate_IJKL_and_coeffs(dtrs, coeff, 0, sym_element, flag)
            end do
        end do

        if (this % NFLG /= 0) this % NFLG = 2
        do ido = 1, num_electrons
            N = ref_dtrs(ido)
            if (this % diagonal_flag > 1 .and. this % orbitals%get_mcon(N) /= 0) cycle

            M = this % orbitals % get_orbital_number(N)
            if (this % positron_flag /= 0) then
              pdtrs(1) = N
              pdtrs(2) = N
              lnadd = this % orbitals % add_positron(pdtrs, 1, 2)
            endif
            call sym_element % insert_ijklm_symbol(0, M, 0, M, lnadd, coeff)
        end do

    end subroutine handle_same


    !> \brief   Evaluates a single integral from labels and also checks for dummy orbitals
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this   Hamiltonian object to update.
    !> \param[in] label     A packed integral label
    !> \param[in] coeff     The computed coefficients from the Slater rules
    !> \param[in] dummy_orb Check wheter to ignore the integral if it contains a dummy orbital of this value
    !>
    !> \return Evaluated integral
    !>
    real(wp) function evaluate_integrals_singular (this, label, coeff, dummy_orb)
        USE ukrmol_interface_gbl, ONLY : GET_ORB_ENERGY

        class(BaseHamiltonian), intent(inout) :: this
        integer,                intent(in)    :: dummy_orb
        real(wp),               intent(in)    :: coeff
        integer(longint),       intent(in)    :: label(NIDX)

        integer :: ido, jdo, num_elements, lwd(8), i, orb_n, orb_num, &
                   belec_orb, velec_orb, posi_orb, enh_factor_type

        real(wp):: integralf, orb_diff, oe, mp2_corr, x1

        real(kind=wp), dimension(4) :: occu_orbe, virt_orbe
        logical :: with_pos


        evaluate_integrals_singular = 0.0

        if (coeff == 0.0_wp) return

        !Search for the dummy orbital
        if (this % phase_flag > 0) then
            call unpack_ints(label, lwd)
            do jdo = 1, 4
                if (LWD(jdo) == this % phase_flag) return
            end do
        else if (this % phase_flag == 0) then
            call unpack_ints(label, lwd)
            do jdo = 1, 4
                if (lwd(jdo) <= 0) cycle
                if (this % orbitals % mcorb(LWD(jdo)) == 0) return
            end do
        end if

        integralf = this % integral % get_integralf(label)
        call unpack_ints(label,lwd)

        !~~~~~~~~~~~ENHANCEMENT FACTOR BIT ~~~~~~~~~~~~~~~~
        with_pos=.False.
        if (lwd(1) > 0) then
          do i=1, 4
            if (this % orbitals % mcorb(LWD(i)) .ne. lwd(i)) then
              with_pos=.True.
            endif
          enddo
        end if

        if (this%enh_factor_type .gt. 1) then
          if (with_pos) then
            belec_orb = 0
            velec_orb = 0
            posi_orb = 0

            do orb_n=1, 4
              virt_orbe(orb_n) = 0.0
              occu_orbe(orb_n) = 0.0

              orb_num = this % orbitals % mcorb(LWD(orb_n)) !orb index

              oe = GET_ORB_ENERGY(orb_num)
              if (oe .gt. 0.0) then !virtual or continuum orbital
                virt_orbe(orb_n) = oe
                if (this % orbitals % mcorb(LWD(orb_n)) .ne. lwd(orb_n)) then
                  velec_orb = velec_orb + 1
                else
                  posi_orb = posi_orb + 1
                endif
              else !bound orbital
                occu_orbe(orb_n) = oe
                if (this % orbitals % mcorb(LWD(orb_n)) .ne. lwd(orb_n)) then
                  belec_orb = belec_orb + 1
                else
                  posi_orb = posi_orb + 1
                endif
              endif
            enddo
            orb_diff = 0.0
            if (belec_orb .eq. 1 .and. velec_orb .eq. 1) then
              !Energy diff between bound and virtual/continuum
              !orbitals
              orb_diff = sum(occu_orbe(1:4)) - sum(virt_orbe(1:4))!Will always be negative
              x1 = integralf
            elseif (this % orbitals % mcorb(LWD(1)) .eq. this % orbitals % mcorb(LWD(2)) .and. &
                    this % orbitals % mcorb(LWD(3)) .eq. this % orbitals % mcorb(LWD(4)) .and. &
                    this % orbitals % mcorb(LWD(1)) .eq. this % orbitals % mcorb(LWD(3)) .and. &
                    velec_orb .ne. 0) then
              !Electron and positron occupying the same orbital -
              !virtual Ps
              orb_diff= 1.0 !Orbital energy that all particles are in.
              x1 = 1
            endif

            if (orb_diff .ne. 0.0) then !include MP2 correction
              if (integralf .ge. 0) then
                mp2_corr = ( (integralf*x1) / orb_diff) * this%enh_factor
              else
                mp2_corr = (-(integralf*x1) / orb_diff) * this%enh_factor
              endif
              integralf = integralf + mp2_corr
            endif
          endif
        else
          if (with_pos) then
            integralf = integralf * this%enh_factor !use provided average enhancement factor. Apply to all (positron) 2-particle integrals
          endif
        endif
        !~~~~~~~~~~~~~END ENHANCEMENT FACTOR BIT~~~~~~~~~~~~~~~~~~~~

        evaluate_integrals_singular = coeff * integralf

        !$OMP ATOMIC
        this % number_of_integrals = this % number_of_integrals + 1

    end function evaluate_integrals_singular


    !> \brief   Evaluates all integrals from a list of symbolic elements
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this           Hamiltonian object to update.
    !> \param[in] symbolic_elements Symbolic elements to be evaluated.
    !> \param[in] dummy_orb         Check wheter to ignore the integral if it contains a dummy orbital of this value.
    !>
    !> \return Total evaluated integrals from symbols.
    !>
    real(wp) function evaluate_integrals (this, symbolic_elements, dummy_orb)
        class(BaseHamiltonian),       intent(inout) :: this
        class(SymbolicElementVector), intent(in)    :: symbolic_elements
        integer,                      intent(in)    :: dummy_orb

        integer          :: ido, num_elements
        integer(longint) :: label(NIDX)
        real(wp)         :: coeff

        num_elements = symbolic_elements % get_size()

        evaluate_integrals = 0.0_wp
        if (num_elements == 0) return
        do ido = 1, num_elements
            call symbolic_elements % get_coeff_and_integral(ido, coeff, label)
            evaluate_integrals = evaluate_integrals + this % evaluate_integrals_singular(label, coeff, dummy_orb)
        end do

    end function evaluate_integrals


    !> \brief   Used by the easy MPI parallelization of slater loops
    logical function my_job (this)
        class(BaseHamiltonian) :: this

        my_job = .true.

        if (grid % gprocs == 1 .or. grid % grank == -1) return

        if (this % job_id /= grid % grank) then
            my_job = .false.
        end if

        this % job_id = this % job_id + 1
        if (this % job_id == grid % gprocs) then
            this % job_id = 0
        end if

    end function my_job


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \deprecated No longer used.
    !>
    integer function fast_slat (rep1, rep2, fin1, fin2, num_rep1, num_rep2, spin_orbs, tot_num_spin_orbs)
        integer, intent(in)    :: num_rep1, rep1(num_rep1), fin1(num_rep1)
        integer, intent(in)    :: num_rep2, rep2(num_rep2), fin2(num_rep2), tot_num_spin_orbs
        integer, intent(inout) :: spin_orbs(tot_num_spin_orbs)

        integer :: i_orbital
        integer :: i_dtrs
        integer :: i_spin_orbital
        integer :: orb_add
        integer :: orb_remove

        do i_dtrs = 1, num_rep1
            orb_remove = (rep1(i_dtrs) + 1) / 2
            orb_add = (fin1(i_dtrs) + 1) / 2
            spin_orbs(orb_remove) = spin_orbs(orb_remove) - 1
            spin_orbs(orb_add) = spin_orbs(orb_add) + 1
        end do

        do i_dtrs = 1, num_rep2
            orb_remove = (rep2(i_dtrs) + 1) / 2
            orb_add = (fin2(i_dtrs) + 1) / 2
            spin_orbs(orb_remove) = spin_orbs(orb_remove) + 1
            spin_orbs(orb_add) = spin_orbs(orb_add) - 1
        end do

        fast_slat = 0

        do i_dtrs = 1, num_rep1
            orb_remove = (rep1(i_dtrs) + 1) / 2
            orb_add = (fin1(i_dtrs) + 1) / 2
            fast_slat = fast_slat + spin_orbs(orb_remove)
            spin_orbs(orb_remove) = 0
            fast_slat = fast_slat + spin_orbs(orb_add)
            spin_orbs(orb_add) = 0
        end do

        do i_dtrs = 1, num_rep2
            orb_remove = (rep2(i_dtrs) + 1) / 2
            orb_add = (fin2(i_dtrs) + 1) / 2
            fast_slat = fast_slat + spin_orbs(orb_remove)
            spin_orbs(orb_remove) = 0
            fast_slat = fast_slat + spin_orbs(orb_add)
            spin_orbs(orb_add) = 0
        end do

    end function fast_slat

end module Hamiltonian_module
