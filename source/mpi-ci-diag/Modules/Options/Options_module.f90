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

!> \brief   Options module
!> \authors A Al-Refaie, J Benda
!> \date    2017 - 2019
!>
!> This module reads input namelists and can be used to pass program run information to other modules.
!> MPI-SCATCI is mostly compatible with SCATCI, except for omission of the following keywords:
!>  - &INPUT: nctgt, gucont, ncont, nobc, thrhm, lusme
!>  - &CINORN: nfta, npflg, itgt, nopvec, ntgt, notgt, ncipfg
!>
!> These keywords are ignored in the input. Additionally, MPI-SCATCI understands the following extra
!> keywords:
!>  - &CINORN: memp, forse, exrc, vecstore, ecp, targmul, targspace, luprop
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial version.
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!> \note 24/01/2019 - Jakub Benda: Read multiple stacked input configurations.
!>
!> \todo Sort the read namelist inputs by `nciset` to avoid error in positioning on output!
!>
module Options_module

    use const_gbl,              only: line_len, stdin, stdout
    use consts_mpi_ci,          only: mxint, NO_CI_TARGET, NORMAL_CSF, GENERATE_TARGET_WF, DIAGONALIZE, DEFAULT_ARCHER_MEMORY, &
                                      SCATCI_DECIDE, ECP_TYPE_NULL, SAVE_ALL_COEFFS, NO_DIAGONALIZATION, DIAGONALIZE_NO_RESTART, &
                                      NORMAL_PHASE, NAME_LEN_MAX, DO_CI_TARGET_CONTRACT, PROTOTYPE_CSF, SYMTYPE_DINFH
    use mpi_gbl,                only: master, mpiint, mpi_xermsg
    use precisn,                only: longint, wp
    use scatci_routines,        only: RDNFTO, movep, mkpt, ham_pars_io
    use MemoryManager_module,   only: master_memory
    use Parallelization_module, only: grid => process_grid

    implicit none

    public OptionsSet, Options, Phase, ham_pars_io_wrapper

    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    type Phase
        integer              :: size_phaze
        integer, allocatable :: phase_index(:)
    end type Phase


    !> \brief   Calculation setup for a single Hamiltonian diagonalization
    !> \authors A Al-Refaie, J Benda
    !> \date    2017 - 2019
    !>
    !> Corresponds to data read from a single pair of &INPUT and &CINORN namelists.
    !> For default values see the subroutines \ref read_all_input and \ref read_all_cinorn.
    !>
    type :: Options
        !--------------SCATCI-INPUTS----------------------------!
        integer          :: megul                       !< MEGUL
        integer          :: ci_target_flag              !< ICITG
        integer          :: csf_type                    !< IEXPC
        integer          :: ci_target_switch            !< NFTG
        real(wp)         :: exotic_mass                 !< SCALEM
        real(wp)         :: enh_factor                  !< ENH_FACTOR
        integer          :: enh_factor_type             !< ENH_FACTOR_TYPE
        integer          :: diagonalization_flag        !< ICIDG
        integer          :: integral_unit               !< NFTI
        integer          :: hamiltonian_unit            !< NFTE
        integer          :: CI_output_unit              !< NFTW
        integer          :: num_matrix_elements_per_rec !< LEMBF
        integer          :: phase_correction_flag       !< NCORB
        integer          :: num_eigenpairs              !< NSTAT
        logical          :: use_UKRMOL_integrals        !< UKRMOLP_INTS
        logical          :: QuantaMolN                  !< QMOLN
        integer          :: integral_ordering           !< IORD
        integer          :: print_dataset_heading       !< NPCVC
        integer          :: output_ci_set_number        !< NCISET
        integer(longint) :: total_memory                !< MEMP
        integer          :: diagonalizer_choice         !< IGH
        integer          :: use_SCF                     !< SCFUSE
        integer          :: force_serial                !< FORSE
        integer          :: print_flags(2)              !< NPFLG

        !> Old_name:Name given to CSF
        character(len=NAME_LEN_MAX) :: name
        character(len=NAME_LEN_MAX) :: diag_name

        !------------Header information on CSF------------------------------!
        integer  :: lambda = 0                      !< MGVN       - Lambda or IRR for CSF
        real(wp) :: spin = 0.0_wp                   !< ISS        - Spin value
        real(wp) :: spin_z = 0.0_wp                 !< ISZ        - Z-Spin
        real(wp) :: refl_quanta = 0.0_wp            !< R          - Reflection quanta Sigma+(1.0) Sigma-(-1.0)
        real(wp) :: pin = 0.0_wp                    !< PIN        - ???
        integer  :: tot_num_orbitals = 0            !< NORB       - Total # of Orbitals
        integer  :: tot_num_spin_orbitals = 0       !< NSRB       - Total # spin orbitals
        integer  :: num_csfs = 0                    !< NOCSF      - # of CSFs
        integer  :: num_electrons = 0               !< NELT       - # of electrons
        integer  :: length_coeff = 0                !< LCDOF      - Length of corresponding coefficients array
        integer  :: matrix_eval = -1                !< IDIAG      - Diagonalization flag
        integer  :: num_syms = 0                    !< NSYM       - # of symmetries
        integer  :: sym_group_flag = 0              !< SYMTYP     - Point Group Flag
        integer  :: length_dtrs = 0                 !< LNDOF      - Length of packed slater determinants
        integer  :: positron_flag = 0               !< IPOSIT     - Positron flag
        integer  :: output_flag(6)                  !< NPFLAG     - Output flag for SCATCI ( or SPEEDY)
        real(wp) :: threshold = 0.0_wp              !< thresh     -
        integer  :: num_continuum = 0               !< nctarg     - Nunber of continuum functions
        integer  :: maximum_num_target_states = 0   !< NMMX       -  Max number of states
        integer  :: total_num_ci_target = 0         !< NCITOT     - Total number of CI components
        integer  :: total_num_target_states = 0     !< NTGT       - total number target states
        integer  :: num_target_sym  = 0             !< ntgtsym    - Number of target symmetries
        integer  :: last_continuum = 0              !< ncont      - Number of last CSF containing a continuum orbital
        integer  :: num_expanded_continuum_CSF = 0  !< Ncont2     - number of expanded continuum functions
        integer  :: seq_num_last_continuum_csf = 0  !< lcont      - Sequence number of last_CSF
        integer  :: continuum_block_size = 0
        integer  :: num_L2_CSF = 0                  !<            - Number of L2 functions
        integer  :: actual_num_csfs = 0             !< MXCSF      - Number of actual CSFS
        integer  :: contracted_mat_size = 0         !< MOCSF      - Contracted matrix size
        integer  :: total_num_tgt_states = 0        !< NTGT       - Total number of starget states
        integer  :: verbosity = 1
        integer  :: exchange_flag = 0               !<NOEX        - Switch off exchange between projectile and target electrons

        !---------------allocatable header information----------------!
        integer,     allocatable :: num_orbitals_sym(:)                     !< NOB    - Number of orbitals per symmetry
        integer,     allocatable :: num_electronic_orbitals_sym(:)          !< NOBE   - Number of electronic orbitals per symmetry
        integer,     allocatable :: num_positron_orbitals_sym(:)            !< NOBP   - Number of positron orbitals per symmetry
        integer,     allocatable :: num_unused_orbitals_sym(:)              !< NOBV   - Number of positron orbitals per symmetry
        integer,     allocatable :: reference_dtrs(:)                       !< NDTRF  - Reference determinants
        integer,     allocatable :: num_dtr_csf_out(:)                      !< NODO   - Number of determinants in CSF in output
        type(Phase), allocatable :: ci_phase(:)                             !< IPHZ   - Index of Phase factor
        integer,     allocatable :: phase_index(:)
        integer,     allocatable :: num_orbitals_sym_dinf(:)                !< NOBL   - The D-inf-h version of NOB
        integer,     allocatable :: num_target_orbitals_sym_congen(:)       !< NOB0   - Number of target orbitals for each symmetry from congen
        integer,     allocatable :: num_target_orbitals_sym_dinf_congen(:)  !< NOB0L  - Number of target orbitals for each symmetry from congen
        integer,     allocatable :: num_target_orbitals_sym(:)              !< NOBC   - Number of target orbitals for each symmetry (and in a sense their starting point)
        integer,     allocatable :: num_ci_target_sym(:)                    !< NCTGT  - Number of CI components of each target symmetry.
        integer,     allocatable :: num_continuum_orbitals_target(:)        !< NOTGT  - Number of continuum orbitals associated with each target state.
        integer,     allocatable :: num_target_state_sym(:)                 !< NUMTGT - Number of target states for each symmetries
        integer,     allocatable :: lambda_continuum_orbitals_target(:)     !< MCONT  - Lambda value ('symmetry') of the continuum orbitals associated with each target state.
        integer,     allocatable :: gerade_sym_continuum_orbital_target(:)  !< GUCONT - Gerade/ungerade symmetry of the continuum orbitals associated with each target state.
        integer,     allocatable :: orbital_sequence_number(:)              !< kpt    - Pointer to orbital sequence
        integer,     allocatable :: ci_set_number(:)                        !< NTGTF  - Set number on CI unit
        integer,     allocatable :: ci_sequence_number(:)                   !< NTGTS  - Set number on CI unit
        real(wp),    allocatable :: energy_shifts(:)                        !< ESHIFT - Energy shifts

        !----------------------miscellaneous variables
        integer,     allocatable :: temp_orbs(:)

        !ECP variables
        integer,     allocatable  :: target_multiplicity(:)
        integer,     allocatable  :: target_spatial(:)
        logical                   :: multiplicity_defined = .false.
        logical                   :: spatial_defined = .false.
        logical                   :: all_ecp_defined = .false.
        integer                   :: ecp_type = ECP_TYPE_NULL
        character(len=line_len*5) :: ecp_filename
        integer                   :: exclude_rowcolumn = -1

        ! Eigenvector output write order
        integer(mpiint) :: preceding_writer = -1    !< Rank of process that will write the previous dataset into the same file (or -1 if none).
        integer(mpiint) :: following_writer = -1    !< Rank of process that will write the following dataset into the same file (or -1 if none).

        real(wp) :: eigenvec_conv = 1d-10
        real(wp) :: residual_conv = 1d-8
        real(wp) :: orthog_trigger = 1d-7
        real(wp) :: max_tolerance = 1d-12
        integer  :: max_iterations = -1

        integer  :: vector_storage_method = SAVE_ALL_COEFFS
    contains
        procedure, public  :: read => read_input
        procedure, private :: read_congen
        procedure, private :: compute_variables
        procedure, private :: compute_expansions
        procedure, private :: check_sanity_congen_header
        procedure, private :: check_sanity_complete
        procedure, private :: arrange_phase
        procedure, public  :: do_ci_contraction
    end type Options


    !> \brief   MPI-SCATCI outer interface
    !> \authors J Benda
    !> \date    2024
    !>
    !> Extract boundary amplitudes, and/or transition multipoles, and/or write RMT data file.
    !>
    !> MPI-SCATCI outer interface uses target properties for construction of
    !>  - inner region transition dipoles (properties read from unit `lutdip`)
    !>  - outer region multipole potentials (properties read from unit `lutarg`)
    !>
    type :: InterfaceOptions
        integer,       allocatable :: idtarg(:)  ! denprop-to-congen state order permutation
        integer,       allocatable :: itspin(:)  ! spin multiplicity of each target spin-symmetry used in CI expansion
        logical  :: write_amp = .false.  ! write boundary amplitudes (swinterf "fort.21")
        logical  :: write_dip = .false.  ! write transition dipole moments (cdenprop_target "fort.624")
        logical  :: write_rmt = .false.  ! write RMT input file (rmt_interface "molecular_data")
        logical  :: all_props = .true.   ! whether to calculate IRR-with-itself properties and energy-sort states (.true. for RMT)
        logical  :: auto_nvo  = .false.  ! Automatically determine contracted virtual orbitals
        integer  :: luprop    = -1       ! orbital properties file unit to read if different than NFTI
        integer  :: lupropw   = 624      ! output unit for CDENPROP properties
        integer  :: lutdip    = 24       ! target properties (dipoles) for use in CDENPROP (typically with ISW = 1)
        integer  :: lutarg    = 24       ! target properties for use in SWINTERF (typically with ISW = 0)
        integer  :: isw       = 0        ! calculate electronic and nuclear (= 0) or electronic only (= 1) propeties in CDENPROP
        integer  :: nfdm                 ! number of extra R-matrix sub-spheres for RMT (excluding the main R-matrix sphere)
        integer  :: luchan    = 10       ! output unit with outer channel descriptions
        integer  :: lurmt     = 21       ! output unit with boundary amplitudes
        real(wp) :: rmatr                ! radius of the R-matrix sphere
        real(wp) :: delta_r              ! spacing between the R-matrix sub-spheres
        character(len=1) :: cform = 'F'  ! formatted or unformatted channel description file
        character(len=1) :: rform = 'U'  ! formatted or unformatted R-matrix file
    end type InterfaceOptions


    !> \brief   MPI-SCATCI input
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> Input configuration of the calculation as a set of \ref Options types. Typically,
    !> there is one \ref Options instance for each Hamiltonian (irreducible representation)
    !> to diagonalize. Also contains data from `&outer_interface` namelists.
    !>
    type :: OptionsSet
        type(Options), allocatable :: opts(:)
        type(InterfaceOptions), allocatable :: interface_opts(:)
    contains
        procedure, public  :: read => read_all_namelists
        procedure, private :: read_all_input
        procedure, private :: read_all_cinorn
        procedure, private :: read_outer_interface
        procedure, public  :: setup_write_order
    end type OptionsSet

contains

    !> \brief   Read all relevant namelists from the input
    !> \authors A Al-Refaie, J Benda
    !> \date    2017 - 2019
    !>
    !> Finds out whether the program received the input file path via the command line; if not, it will be read
    !> from the standard input. Then reads all &INPUT namelists, followed by all &CINORN namelists. The data read
    !> from the namelists are used to initialize objects of type \ref Options.
    !>
    subroutine read_all_namelists (this)

        class(OptionsSet), intent(inout) :: this

        integer             :: num_arguments, io, err, num_input, num_cinorn, i
        character(len=800)  :: input_filename

        num_arguments = command_argument_count()

        if (num_arguments > 1) then
            write (stdout, "('Only one argument allowed')")
            write (stdout, "('Command line must be:')")
            write (stdout, "('./mpi-scatci <input filename>')")
            call mpi_xermsg('Options_module', 'read_all_namelists', 'Too many command-line arguments', 1, 1)
        end if

        if (num_arguments == 1) then
            call get_command_argument(1, input_filename)
            input_filename = trim(input_filename)

            open (newunit = io, action = 'read', status = 'old', file = input_filename, iostat = err)

            if (err /= 0) then
                call mpi_xermsg('Options_module', 'read_all_namelists', 'Could not open input file ' // trim(input_filename), 2, 1)
            end if
        else
            io = stdin
        end if

        if (allocated(this % opts)) deallocate(this % opts)

        num_input  = this % read_all_input(io)
        num_cinorn = this % read_all_cinorn(io)

        call this % read_outer_interface(io)

        if (io /= stdin) close (io)

        ! we require at least one valid input
        if (num_input == 0) then
            call mpi_xermsg('Options_module', 'read_all_namelists', 'Missing &INPUT namelist in the input file', 3, 1)
        end if

        ! we require equal number of &INPUT and &CINORN namelists
        if (num_input /= num_cinorn) then
            call mpi_xermsg('Options_module', 'read_all_namelists', 'Unequal number of &INPUT and &CINORN namelists', 4, 1)
        end if

        ! let each Options instance finalize itself
        do i = 1, num_input
            call this % opts(i) % read
        end do

    end subroutine read_all_namelists


    !> \brief   Read all &INPUT namelists from file unit
    !> \authors J Benda
    !> \date    2019 - 2022
    !>
    !> Reads all &INPUT namelists in the order present in the source file and stores the data
    !> in the internal array (reallocated whenever necessary).
    !>
    !> \param[inout] this  OptionSet object to update.
    !> \param[in]    io    File unit to process.
    !>
    !> \return  Number of &INPUT namelists read from the file.
    !>
    integer function read_all_input (this, io) result (i)

        use iso_fortran_env, only: iostat_end

        class(OptionsSet), intent(inout) :: this
        integer,           intent(in)    :: io

        type(Options), allocatable :: tmp_opts(:)
        character(NAME_LEN_MAX)    :: name

        logical :: qmoln, ukrmolp_ints
        integer :: icidg, icitg, iexpc, idiag, nfti, nfte, lembf, megul, nftg, ntgsym, ncont, ncorb, iord, scfuse, lusme
        integer :: io_stat, ierr, nset, nrec, nnuc, nocsf, nstat, mgvn, nelt, k, j, iposit, enh_factor_type, NOEX

        character(1000)      :: io_msg
        real(wp)             :: scalem, thrhm, spin, spinz, e0, enh_factor
        integer, allocatable :: numtgt(:), notgt(:), nctgt(:), ntgtf(:), ntgts(:), mcont(:), gucont(:), nobc(:)

        namelist /INPUT/ icidg, icitg, iexpc, idiag, nfti, nfte, lembf, megul, nftg, ntgsym, numtgt, notgt, nctgt, ntgtf, ntgts, &
                         mcont, gucont, ncont, ncorb, nobc, iord, name, scalem, scfuse, thrhm, qmoln, ukrmolp_ints, lusme, iposit, &
                         enh_factor, enh_factor_type, NOEX

        i = 0

        allocate (numtgt(mxint), notgt(mxint), nctgt(mxint), ntgtf(mxint), ntgts(mxint), mcont(mxint), gucont(mxint), &
                  nobc(mxint), stat = ierr)

        if (ierr /= 0) then
            call mpi_xermsg('Options_module', 'read_all_input', 'Input arrays allocation failure', 1, 1)
        end if

        rewind (io, iostat = io_stat)

        namelist_loop: do

            ! set namelist entries to their default values
            icidg     = DIAGONALIZE
            icitg     = NO_CI_TARGET
            iexpc     = NORMAL_CSF
            idiag     = -1
            nfti      = 16
            nfte      = 26
            lembf     = 5000
            megul     = 13
            nftg      = GENERATE_TARGET_WF
            ntgsym    = 0
            numtgt(:) = 0
            notgt(:)  = 0
            nctgt(:)  = 1
            ntgtf(:)  = 0
            ntgts(:)  = 0
            mcont(:)  = 0
            ncorb     = NORMAL_PHASE
            iord      = 1
            name      = ''
            scalem    = 0
            scfuse    = 0
            iposit    = 0
            qmoln     = .false.
            enh_factor = 1.0_wp
            enh_factor_type = 0
            ukrmolp_ints = .true.
            NOEX         = 0

            ! read the next &INPUT namelist
            read (io, nml = INPUT, iostat = io_stat, iomsg = io_msg)
            if (io_stat == iostat_end) then
                exit namelist_loop  ! no more namelists in the file...
            end if
            if (io_stat /= 0) then
                call mpi_xermsg('Options_module', 'read_all_input', &
                                'Error when reading &INPUT namelist: ' // trim(io_msg), i + 1, 1)
            end if

            ! seems to be a valid namelist
            i = i + 1

            ! re-allocate array of Options structures
            if (allocated(this % opts)) then
                call move_alloc(this % opts, tmp_opts)
            end if
            allocate (this % opts(i), stat = ierr)
            if (ierr /= 0) then
                call mpi_xermsg('Options_module', 'read_all_input', 'Input structures allocation failure', 2, 1)
            end if
            if (allocated(tmp_opts)) then
                this % opts(1:size(tmp_opts)) = tmp_opts(:)
                deallocate (tmp_opts)
            end if

            ! prevent overflows in below assigments (but should be assured by length check of the namelist data already)
            if (ntgsym > mxint) then
                call mpi_xermsg('Options_module', 'read_all_input', 'NTGSYM exceeds built-in limit MXINT', 3, 1)
            end if
            this % opts(i) % total_num_tgt_states = sum(numtgt(1:ntgsym))
            if (this % opts(i) % total_num_tgt_states > mxint) then
                call mpi_xermsg('Options_module', 'read_all_input', 'Number of target states exceeds built-in limit MXINT', 4, 1)
            end if

            ! copy data from the namelist to the structure
            this % opts(i) % diagonalization_flag        = icidg
            this % opts(i) % ci_target_flag              = icitg
            this % opts(i) % csf_type                    = iexpc
            this % opts(i) % matrix_eval                 = idiag
            this % opts(i) % integral_unit               = nfti
            this % opts(i) % hamiltonian_unit            = nfte
            this % opts(i) % num_matrix_elements_per_rec = lembf
            this % opts(i) % megul                       = megul
            this % opts(i) % ci_target_switch            = nftg
            this % opts(i) % num_target_sym              = ntgsym

            if (ntgsym > 0) then
                call allocate_integer_array(this % opts(i) % num_target_state_sym, ntgsym)
                call allocate_integer_array(this % opts(i) % num_continuum_orbitals_target, ntgsym)
                call allocate_integer_array(this % opts(i) % lambda_continuum_orbitals_target, ntgsym)
                call allocate_integer_array(this % opts(i) % ci_set_number, this % opts(i) % total_num_tgt_states)
                call allocate_integer_array(this % opts(i) % ci_sequence_number, this % opts(i) % total_num_tgt_states)
                this % opts(i) % num_target_state_sym             = numtgt(1:ntgsym)
                this % opts(i) % num_continuum_orbitals_target    = notgt(1:ntgsym)
                this % opts(i) % ci_set_number                    = ntgtf(1:this % opts(i) % total_num_tgt_states)
                this % opts(i) % ci_sequence_number               = ntgts(1:this % opts(i) % total_num_tgt_states)
                this % opts(i) % lambda_continuum_orbitals_target = mcont(1:ntgsym)
            end if

            this % opts(i) % phase_correction_flag       = ncorb
            this % opts(i) % integral_ordering           = iord
            this % opts(i) % name                        = name
            this % opts(i) % exotic_mass                 = scalem
            this % opts(i) % use_SCF                     = scfuse
            this % opts(i) % positron_flag               = iposit
            this % opts(i) % enh_factor                  = enh_factor
            this % opts(i) % enh_factor_type             = enh_factor_type
            this % opts(i) % QuantaMolN                  = qmoln
            this % opts(i) % use_UKRMOL_integrals        = ukrmolp_ints
            this % opts(i) % exchange_flag               = NOEX

            !if ( this % opts(i) % positron_flag /= 0 ) then
            !    call mpi_xermsg('Options_module', 'read_all_input', &
            !                    'Positrons are not yet implemented for mpi-scatci.', i, 1)
            !end if

            ! also read spin-symmetry information about the target states from the target CI expansion coefficients file (if given)
            if (nftg > 0) then
                allocate (this % opts(i) % target_spatial(this % opts(i) % total_num_tgt_states), &
                          this % opts(i) % target_multiplicity(this % opts(i) % total_num_tgt_states))
                do j = 1, this % opts(i) % total_num_tgt_states
                    if (this % opts(i) % ci_set_number(j) <= 0) then
                        call mpi_xermsg('Options_module', 'read_all_input', &
                                        'Non-zero NFTG requires setting also non-zero NTGTF/NTGTS.', i, 1)
                    end if
                    call movep(nftg, this % opts(i) % ci_set_number(j), ierr, 0, 0)
                    if (ierr /= 0) then
                        call mpi_xermsg('Options_module', 'read_all_input', &
                                        'CI data not found in unit', nftg, 1)
                    end if
                    read (nftg)  ! dummy read to skip header
                    read (nftg) nset, nrec, name, nnuc, nocsf, nstat, mgvn, spin, spinz, nelt, e0
                    this % opts(i) % target_spatial(j) = mgvn
                    this % opts(i) % target_multiplicity(j) = nint(2 * spin + 1)
                end do
            end if

        end do namelist_loop

    end function read_all_input


    !> \brief   Read all &CINORN namelists from file unit
    !> \authors J Benda
    !> \date    2019
    !>
    !> Reads all &CINORN namelists in the order present in the source file and stores the data
    !> in the internal array (reallocated whenever necessary).
    !>
    !> \param[inout] this  OptionSet object to update.
    !> \param[in]    io    File unit to process.
    !>
    !> \return  Number of &CINORN namelists read from the file.
    !>
    integer function read_all_cinorn (this, io) result (i)

        use iso_fortran_env, only: iostat_end
        use const_gbl,       only: set_verbosity_level

        class(OptionsSet), intent(inout) :: this
        integer,           intent(in)    :: io

        type(Options), allocatable :: tmp_opts(:)

        integer  :: nfta, nftw, nciset, npflg(2), npcvc, itgt, nopvec, ntgt, notgt, ncipfg, nkey, large, thrprt
        integer  :: keycsf, nstat, igh, maxiter, forse, exrc, vecstore, ecp, targmul(mxint), targspace(mxint)
        integer  :: ierr, io_stat, nshifts, ishift, verbosity
        real(wp) :: critc, critr, ortho, crite, memp, eshift(mxint)

        character(256)          :: io_msg
        character(NAME_LEN_MAX) :: name

        namelist /CINORN/ nfta, nftw, nciset, npflg, npcvc, name, itgt, nopvec, ntgt, notgt, eshift, ncipfg, nkey, large, thrprt, &
                          keycsf, nstat, igh, critc, critr, ortho, crite, maxiter, memp, forse, exrc, vecstore, ecp, &
                          targmul, targspace, verbosity

        i = 0

        rewind (io, iostat = io_stat)

        namelist_loop: do

            ! set namelist entries to their default values
            nftw      = 25
            nciset    = 1
            npcvc     = 1
            name      = ''
            nstat     = 0
            igh       = SCATCI_DECIDE
            critc     = 1e-10_wp
            critr     = 1e-08_wp
            ortho     = 1e-07_wp
            crite     = 1e-12_wp
            maxiter   = -1
            memp      = DEFAULT_ARCHER_MEMORY
            forse     = 0
            exrc      = -1
            vecstore  = SAVE_ALL_COEFFS
            ecp       = ECP_TYPE_NULL
            targmul   = 0
            targspace = 0
            eshift    = 0.0_wp
            npflg     = 0
            verbosity = 1

            ! read the next &CINORN namelist
            read (io, nml = CINORN, iostat = io_stat, iomsg = io_msg)
            if (io_stat == iostat_end) then
                exit namelist_loop  ! no more namelists in the file...
            end if
            if (io_stat /= 0) then
                call mpi_xermsg('Options_module', 'read_all_input', &
                                'Error when reading &CINORN namelist: ' // trim(io_msg), i + 1, 1)
            end if

            ! seems to be a valid namelist
            i = i + 1

            ! re-allocate array of Options structures as needed
            if (.not. allocated(this % opts) .or. size(this % opts) < i) then
                if (allocated(this % opts)) call move_alloc(this % opts, tmp_opts)
                allocate (this % opts(i), stat = ierr)
                if (ierr /= 0) then
                    call mpi_xermsg('Options_module', 'read_all_cinorn', 'Input structures allocation failure', 2, 1)
                end if
                if (allocated(tmp_opts)) this % opts(1:size(tmp_opts)) = tmp_opts(:)
            end if

            ! energy shifts: count them (up to the last non-zero specified), allocate a permanent array, and copy data
            nshifts = 0
            do ishift = 1, size(eshift)
                if (eshift(ishift) /= 0) nshifts = ishift
            end do
            if (nstat > 0 .and. nshifts > nstat) then
                call mpi_xermsg('Options_module', 'read_all_cinorn', &
                                'Too many energy shifts given for the given number of eigenpairs', nstat, 1)
            end if
            call allocate_real_array(this % opts(i) % energy_shifts, nshifts)
            if (nshifts > 0) then
                this % opts(i) % energy_shifts = eshift(1:nshifts)
            end if

            ! copy data from the namelist to the structure
            this % opts(i) % CI_output_unit        = nftw
            this % opts(i) % output_ci_set_number  = nciset
            this % opts(i) % print_dataset_heading = npcvc
            this % opts(i) % diag_name             = name
            this % opts(i) % num_eigenpairs        = nstat
            this % opts(i) % diagonalizer_choice   = igh
            this % opts(i) % eigenvec_conv         = critc
            this % opts(i) % residual_conv         = critr
            this % opts(i) % orthog_trigger        = ortho
            this % opts(i) % max_tolerance         = crite
            this % opts(i) % max_iterations        = maxiter
            this % opts(i) % print_flags           = npflg

            ! MPI-SCATCI specific parameters
            this % opts(i) % total_memory          = int(memp * 1024_longint**3, longint)
            this % opts(i) % force_serial          = forse
            this % opts(i) % exclude_rowcolumn     = exrc
            this % opts(i) % vector_storage_method = vecstore
            this % opts(i) % ecp_type              = ecp
            this % opts(i) % verbosity             = verbosity

            call set_verbosity_level(verbosity)

        end do namelist_loop

    end function read_all_cinorn


    !> \brief   Read contents of the OUTER interface namelist
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> Reads setup(s) of the OUTER interfaces(s).
    !>
    !> \param[inout] this  OptionSet object to update.
    !> \param[in]    io    File unit to process.
    !>
    subroutine read_outer_interface (this, io)

        use iso_fortran_env, only: iostat_end

        class(OptionsSet), intent(inout) :: this
        integer,           intent(in)    :: io

        type(InterfaceOptions), allocatable :: options(:)
        type(InterfaceOptions)              :: option

        logical  :: write_amp, write_dip, write_rmt, all_props, auto_nvo, is_scattering
        integer  :: luprop, lupropw, lutdip, lutarg, luchan, lurmt, nfdm, io_stat, ntarg, idtarg(mxint), itarget_spin(mxint), isw, n
        real(wp) :: delta_r, rmatr

        character(1)   :: cform, rform
        character(256) :: io_msg
        character(900) :: err_msg

        namelist /OUTER_INTERFACE/ write_amp, write_dip, write_rmt, all_props, luprop, nfdm, delta_r, rmatr, ntarg, idtarg, &
                                   lutdip, lutarg, lupropw, luchan, lurmt, cform, rform, isw, auto_nvo, itarget_spin

        rewind (io, iostat = io_stat)

        namelist_loop: do

            write_amp = .true.
            write_dip = .false.
            write_rmt = .false.
            all_props = .true.
            auto_nvo  = .false.
            lutdip    = 24
            lutarg    = 24
            luprop    = -1
            isw       = 0
            nfdm      = 18
            delta_r   = 0.08_wp
            rmatr     = -1.0_wp
            ntarg     = 0
            idtarg    = 0
            itarget_spin = 2
            lupropw   = 624
            luchan    = 10
            lurmt     = 21
            cform     = 'F'
            rform     = 'U'

            read (io, nml = OUTER_INTERFACE, iostat = io_stat, iomsg = io_msg)
            if (io_stat == iostat_end) then
                exit namelist_loop  ! no more namelists in the file...
            end if
            if (io_stat /= 0) then
                call mpi_xermsg('Options_module', 'read_outer_interface', &
                                'Error when reading &OUTER_INTERFACE namelist: ' // trim(io_msg), 1, 1)
            end if
            if (rmatr < 0 .and. (write_amp .or. write_rmt)) then
                call mpi_xermsg('Options_module', 'read_outer_interface', 'Missing "rmatr" in &outer_interface namelist.', 1, 1)
            end if

            option % write_amp = write_amp
            option % write_dip = write_dip
            option % write_rmt = write_rmt
            option % all_props = all_props
            option % auto_nvo  = auto_nvo
            option % lutdip    = lutdip
            option % lutarg    = lutarg
            option % luprop    = luprop
            option % itspin    = itarget_spin
            option % isw       = isw
            option % nfdm      = nfdm
            option % delta_r   = delta_r
            option % rmatr     = rmatr
            option % lupropw   = lupropw
            option % luchan    = luchan
            option % lurmt     = lurmt
            option % cform     = cform
            option % rform     = rform

            if (ntarg > 0) then
                option % idtarg = idtarg(1:ntarg)
                if (any(option % idtarg <= 0)) then
                    call mpi_xermsg('Options_module', 'read_outer_interface', &
                                    'IDTARG must be provided for all NTARG target states.', 1, 1)
                end if
            else if (option % write_amp .or. option % write_rmt) then
                if (any(this % opts(:) % ci_target_switch == 0)) then
                    call mpi_xermsg('Options_module', 'read_outer_interface', &
                                    'When NFTG = 0, outer interface requires setting NTARG and IDTARG.', 1, 1)
                end if
            end if

            ! determine if this is a CI-contracted (scattering) run
            is_scattering = any([(this % opts(n) % do_ci_contraction(), n = 1, size(this % opts))])

            ! test availability of target (full molecular) properties
            if (write_amp .or. write_rmt .or. (write_dip .and. is_scattering)) then
                open (lutarg, status = 'old', action = 'read', form = 'formatted', iostat = io_stat, iomsg = io_msg)
                if (io_stat /= 0) then
                    write (err_msg, '(a,i0,a)') 'Missing target properties (lutarg = ', lutarg, '):'//new_line(' ')//' *  '//io_msg
                    call mpi_xermsg('Options_module', 'read_outer_interface', err_msg, 1, 1)
                end if
                close (lutarg)
            end if

            ! test availability target transition dipoles
            if (write_rmt .or. (write_dip .and. is_scattering)) then
                open (lutdip, status = 'old', action = 'read', form = 'formatted', iostat = io_stat, iomsg = io_msg)
                if (io_stat /= 0) then
                    write (err_msg, '(a,i0,a)') 'Missing target dipoles (lutdip = ', lutdip, '):'//new_line(' ')//' *  '//io_msg
                    call mpi_xermsg('Options_module', 'read_outer_interface', err_msg, 2, 1)
                end if
                close (lutdip)
            end if

            call move_alloc(this % interface_opts, options)
            n = 0
            if (allocated(options)) then
                n = size(options)
            end if
            allocate (this % interface_opts(n + 1))

            if (n > 0) then
                this % interface_opts(1 : n) = options
            end if
            this % interface_opts(n + 1) = option

        end do namelist_loop

    end subroutine read_outer_interface


    !> \brief   Read and check input files
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine read_input (this)
        class(Options), intent(inout) :: this

        call this % read_congen

        call this % compute_variables
        call this % compute_expansions
        call this % check_sanity_complete

        call master_memory % construct(this % total_memory)
        call master_memory % print_memory_report

    end subroutine read_input


    subroutine check_sanity_congen_header (this)
        class(Options), intent(inout) :: this

        !   write(stdout,"(' ',A6,': ',A)") 'SCATCI',this%name
        !       write(stdout,16) this%lambda,this%spin,this%num_electrons,this%spin_z,this%matrix_eval,this%refl_quanta,this%threshold,&
        !         & this%num_syms,this%MEGUL,this%num_csfs,this%output_flag(1:6)
        !       16 FORMAT(' MGVN =',I10,5X,'S    =',F5.1,/,' NELT =',I10,5X,'SZ   =' &
        !              &       ,F5.1,/,' IDIAGT=',I10,5X,'R    =',F5.1,/, &
        !             &       ' THRES=',D10.1,5X,'NSYM =',I5,/,' MEGUL=',I3,', NOCSF=' &
        !             &       ,I10,/,' NPFLG',I10,5I4)

        if (this % num_csfs <= 0) THEN
            write (stdout, "('Options::check_sanity_congen_header - No CSFS detected...stopping')")
            stop 'No of CSFS = 0'
        end if

    end subroutine check_sanity_congen_header


    !> \brief   Print a subset of options read from the input
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Contrary to its name does not actually check anything.
    !>
    subroutine check_sanity_complete (this)
        class(Options), intent(inout) :: this

        write (stdout, *)
        write (stdout, "('---------------SCATCI run options------------------')")
        write (stdout, "('Name of run is ', a)") trim(this % name)
        write (stdout, "('CONGEN unit ', i4)") this % megul
        write (stdout, "('Hamiltonian unit ', i4)") this % hamiltonian_unit
        write (stdout, "('Number of hamiltonian records to write ', i4)") this % num_matrix_elements_per_rec
        write (stdout, "('Integral unit ', i4)") this % integral_unit

        write (stdout, "('Integral Type - ')", advance = 'no')
        if (this % sym_group_flag <= 1) then
            write (stdout, "('ALCHEMY')")
        else if (this % sym_group_flag == 2 .and. this % use_UKRMOL_integrals) then
            write (stdout, "('UKRMOL+')")
        else
            write (stdout, "('SWEDEN')")
        end if

        if (this%enh_factor_type .eq. 1) then
          write(stdout,'(/,"Using an enhancement factor of ",f10.7)') this%enh_factor
        elseif (this%enh_factor_type .eq. 2) then
          this%enh_factor = 1.0
          write(stdout,'(/,"Computing MP2 enhancement factor")')
        elseif (this%enh_factor_type .eq. 3) then
          write(stdout,'(/,"Computing enhanced MP2 enhancement factor using ",f10.7)')this%enh_factor
        else
          this%enh_factor = 1.0
        endif

        write (stdout, "('CI target Contraction? ')", advance = 'no')
        if (this % ci_target_flag == NO_CI_TARGET) then
            write (stdout, "('No')")
        else
            write (stdout, "('Yes')")
            write (stdout, "('Generate CI vectors? ')", advance = 'no')
            if (this % ci_target_switch == GENERATE_TARGET_WF) then
                write (stdout, "('Yes')")
            else
                write (stdout, "('No, We are reading from unit ',i4)") this % ci_target_switch
            end if
        end if

        write (stdout, "('Do we diagonalize? ')", advance = 'no')
        if (this % diagonalization_flag == NO_DIAGONALIZATION) then
            write (stdout, "('No')")
        else if (this % diagonalization_flag == DIAGONALIZE) then
            write (stdout, "('Yes with restart for ARPACK')")
        else if (this % diagonalization_flag == DIAGONALIZE_NO_RESTART) then
            write (stdout, "('Yes with NO restart for ARPACK')")
        endif

        write (stdout, "('Num eigenpairs = ',i10)") this % num_eigenpairs
        write (stdout, *)
        write (stdout, "('   -----------Molecule options------------    ')")
        write (stdout, "('Number of electrons = ', i10)") this % num_electrons
        write (stdout, "('Number of molecular orbitals = ', i10)") this % tot_num_orbitals
        write (stdout, "('Number of spin orbitals = ', i10)") this % tot_num_spin_orbitals
        write (stdout, *)
        write (stdout, "('   -----------CSF options------------    ')")
        write (stdout, "('Number of configuration state functions = ', i10)") this % num_csfs
        write (stdout, "('Number of continuum functions = ', i10)") this % last_continuum
        this % num_L2_CSF = this % num_csfs - this % last_continuum
        write (stdout, "('Number of L2 functions = ', i10)") this % num_L2_CSF
        write (stdout, *)
        write (stdout, "('   -----------Orbital/Symmetry options------------    ')")
        write (stdout, "('Number of symmetries = ', i4)") this % num_syms
        write (stdout, "('Number of orbitals per symmetry = ', 20i4)") this % num_orbitals_sym(:)
        write (stdout, "('Number of electronic orbitals per symmetry = ', 20i4)") this % num_electronic_orbitals_sym(:)
        write (stdout, "('Number of exotic orbitals per symmetry = ', 20i4)") this % num_positron_orbitals_sym(:)

        write (stdout, "('Use ECP? ')", advance = 'no')
        if (this % ecp_type == ECP_TYPE_NULL) then
            write (stdout, "('No')")
        else
            write (stdout, "('Yes')")
        end if

    end subroutine check_sanity_complete


    !> \brief   Reads all input information
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Reads saved setup from the congen output binary file.
    !>
    !> \todo Switch from old name list way into newer more robust input reading (use input.f90 like in TROVE?)
    !>
    subroutine read_congen (this)

        class(Options), intent(inout) :: this

        !Dummy variables for RDNFTO
        integer :: ido, nobep, ntgcon, idiagt
        character(len=NAME_LEN_MAX) :: congen_name

        !Begin reading first round of header information
        rewind (this % megul)
        read (this % megul) congen_name, this % lambda, this % spin, this % spin_z, this % refl_quanta, this % pin,       &
                            this % tot_num_orbitals, this % tot_num_spin_orbitals, this % num_csfs, this % num_electrons, &
                            this % length_coeff,idiagt, this % num_syms, this % sym_group_flag, this % length_dtrs,       &
                            this % output_flag, this % threshold, this % num_continuum, ntgcon

        if (trim(this % name) == '') this % name = congen_name
        if (this % matrix_eval < 0) this % matrix_eval = idiagt
        if (ntgcon > 0) this % num_target_sym = ntgcon

        call this % check_sanity_congen_header()

        !allocate header arrays
        call allocate_integer_array(this % num_orbitals_sym,                    this % num_syms)
        call allocate_integer_array(this % num_electronic_orbitals_sym,         this % num_syms)
        call allocate_integer_array(this % num_positron_orbitals_sym,           this % num_syms)
        call allocate_integer_array(this % num_unused_orbitals_sym,             this % num_syms)
        call allocate_integer_array(this % reference_dtrs,                      this % num_electrons)
        call allocate_integer_array(this % num_dtr_csf_out,                     this % num_csfs)
        call allocate_integer_array(this % phase_index,                         this % num_continuum)
        call allocate_integer_array(this % num_orbitals_sym_dinf,               this % num_syms * 2)
        call allocate_integer_array(this % num_target_orbitals_sym_congen,      this % num_syms)
        call allocate_integer_array(this % num_target_orbitals_sym,             this % num_syms)
        call allocate_integer_array(this % num_target_orbitals_sym_dinf_congen, this % num_syms * 2)
        call allocate_integer_array(this % num_ci_target_sym,                   this % num_target_sym)
        call allocate_integer_array(this % temp_orbs,                           this % num_target_sym)
        call allocate_integer_array(this % num_continuum_orbitals_target,       this % num_target_sym)
        call allocate_integer_array(this % lambda_continuum_orbitals_target,    this % num_target_sym)
        call allocate_integer_array(this % gerade_sym_continuum_orbital_target, this % num_target_sym)
        call allocate_integer_array(this % num_target_state_sym,                this % num_target_sym)
        call allocate_integer_array(this % orbital_sequence_number,             this % num_csfs)

        call RDNFTO(this % megul, this % num_orbitals_sym, this % num_target_orbitals_sym_congen, this % num_orbitals_sym_dinf,  &
                    this % num_target_orbitals_sym_dinf_congen, this % num_syms, this % reference_dtrs, this % num_electrons,    &
                    this % num_dtr_csf_out, this % num_csfs, this % phase_index, this % num_continuum, this % num_ci_target_sym, &
                    this % temp_orbs, this % lambda_continuum_orbitals_target, this % gerade_sym_continuum_orbital_target,       &
                    this % num_target_sym, this % num_target_sym, this % num_electronic_orbitals_sym,                            &
                    this % num_positron_orbitals_sym, this % num_unused_orbitals_sym, this % positron_flag , this % exchange_flag)

        do ido = 1, this % num_syms
            if (this % num_positron_orbitals_sym(ido) == 0) then
                this % num_electronic_orbitals_sym(ido) = this % num_orbitals_sym(ido)
                write (stdout, *) 'NOB(i)=', this % num_orbitals_sym(ido)
                write (stdout, *) 'setting NOBE(i)=', this % num_electronic_orbitals_sym(ido)
            end if
            nobep = this % num_electronic_orbitals_sym(ido) + this % num_positron_orbitals_sym(ido)
            if (nobep /= this % num_orbitals_sym(ido)) then
                write (stdout, *) 'ERROR on input:'
                write (stdout, *) 'not: NOB(i)=NOBE(i)+NOBP(i)'
                write (stdout, *) 'i=', ido
                write (stdout, *) 'NOB(i)=', this % num_orbitals_sym(ido)
                write (stdout, *) 'NOBE(i)=', this % num_electronic_orbitals_sym(ido)
                write (stdout, *) 'NOBP(i)=', this % num_positron_orbitals_sym(ido)
            end if
            if (this % num_unused_orbitals_sym(ido) == 0) then
                this % num_unused_orbitals_sym(ido) = this % num_target_orbitals_sym_congen(ido)
            end if
            if (this % num_target_orbitals_sym_congen(ido) > this % num_orbitals_sym(ido)) then
                write (stdout, *) 'ERROR on input:'
                write (stdout, *) 'NOB0(i) > NOB(i)'
                write (stdout, *) 'i=', ido
                write (stdout, *) 'NOB(i)=', this % num_orbitals_sym(ido)
                write (stdout, *) 'NOB0(i)=', this % num_target_orbitals_sym_congen(ido)
            end if
        end do

        ! Allocate the phases for each target symmetry (must be allocated as it is used in main program)
        allocate(this % ci_phase(this % num_target_sym))

        ! Set up the phases
        !call this % arrange_phase(this % phase_index)

    end subroutine read_congen


    subroutine arrange_phase (this, temp_phase)
        class(Options), intent(inout) :: this
        integer,        intent(in)    :: temp_phase(:)
        integer                       :: ido, ci_num

        ci_num = 0

        if (sum(this % num_ci_target_sym) == this % num_continuum) then
            do ido = 1, this % num_target_sym
                this % ci_phase(ido) % size_phaze = this % num_ci_target_sym(ido)
                allocate(this % ci_phase(ido) % phase_index(this % ci_phase(ido) % size_phaze))
                this % ci_phase(ido) % phase_index(1:this % ci_phase(ido) % size_phaze) = &
                    temp_phase(ci_num + 1 : ci_num + this % ci_phase(ido) % size_phaze)
                ci_num = ci_num + this % ci_phase(ido) % size_phaze
            end do
        else if (this % num_continuum == 0) then
            return
        else
            stop "[Options] Problem with phase not matching sum of NCTGT"
        end if

    end subroutine arrange_phase


    !> \brief   Set up write order flags
    !> \authors J Benda
    !> \date    2019
    !>
    !> To prevent mutual overwriting of the common output file, the MPI group masters, who write the eigenvectors to disk,
    !> need to communicate, so that the next set written by the next process commences only once the previous write is
    !> finished.
    !>
    subroutine setup_write_order (this)

        class(OptionsSet), intent(inout) :: this

        integer :: i, j, k, n

        ! set up preceding/following relations for writing datasets in the same files
        do i = 1, size(this % opts)

            n = 0   ! number of datasets written to the same file before this process
            k = 0   ! which setup index corresponds to the last write to the same file

            do j = 1, i - 1
                if (this % opts(j) % CI_output_unit == this % opts(i) % CI_output_unit) then
                    n = n + 1
                    k = j
                end if
            end do

            ! if automatic numbering of datasets is used, get the true unit number for this setup
            if (this % opts(i) % output_ci_set_number == 0) then
                this % opts(i) % output_ci_set_number = n + 1
            end if

            ! if there is a previous dataset in the same file, link them to avoid lock/overwriting
            if (k > 0) then
                if (this % opts(i) % output_ci_set_number <= this % opts(k) % output_ci_set_number) then
                    call mpi_xermsg('Options_module', 'read_all_namelists', &
                                    'The entries "nciset" must form ascending sequence.', i, 1)
                end if
                this % opts(k) % following_writer = grid % group_master_world_rank(grid % which_group_is_work(i))
                this % opts(i) % preceding_writer = grid % group_master_world_rank(grid % which_group_is_work(k))
            end if

        end do

    end subroutine setup_write_order


    !> \brief   Computes additional variables required by scatci
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine compute_variables (this)
        class(Options), intent(inout) :: this
        integer :: itgtsym, isym, I, num_orb_temp

        if (this % csf_type == PROTOTYPE_CSF) then
            this % num_target_orbitals_sym = 0
            do itgtsym = 1, this % num_target_sym

                !If we're working with D_inf_h then we need to 'double' the symmtery
                if (this % sym_group_flag /= SYMTYPE_DINFH) then
                    isym = this % lambda_continuum_orbitals_target(itgtsym) + 1
                else
                    isym = (this % lambda_continuum_orbitals_target(itgtsym) * 2) + 1
                    if (mod(this % lambda_continuum_orbitals_target(itgtsym), 2) == 0 .and. &
                        this % gerade_sym_continuum_orbital_target(itgtsym) == -1) isym = isym + 1
                    if (mod(this % lambda_continuum_orbitals_target(itgtsym), 2) == 1 .and. &
                        this % gerade_sym_continuum_orbital_target(itgtsym) ==  1) isym = isym + 1
                end if

                !If there is already a value here then don't bother (impossible since the array must be allocated)
                if (this % num_target_orbitals_sym(isym) > 0) cycle

                num_orb_temp = this % num_orbitals_sym(isym) - this % num_continuum_orbitals_target(itgtsym)

                if (num_orb_temp < this % num_target_orbitals_sym_dinf_congen(isym)) then
                    write (stdout, '(/,A)') 'Error: Incompatible number of target orbitals between CONGEN and SCATCI'
                    write (stdout, '(7X,*(A,I0))') 'Target type ', itgtsym, ' with MGVN=', isym - 1
                    write (stdout, '(7X,*(A,I0))') 'Total number of orbitals for MGVN=', isym - 1, ': ', &
                        this % num_orbitals_sym(isym)
                    write (stdout, '(7X,*(A,I0))') 'Number of continuum orbitals associated with target type: ', &
                        this % num_continuum_orbitals_target(itgtsym)
                    write (stdout, '(7X,*(A,I0))') 'Number of target orbitals for MGVN=', isym - 1, ' from CONGEN: ', &
                        this % num_target_orbitals_sym_dinf_congen(isym)

                    call mpi_xermsg('Options_module', 'compute_variables', &
                                    'Incompatible number of target orbitals. Wrong CONGEN file or NOTGT?', itgtsym, 1)
                end if

                this % num_target_orbitals_sym(isym) = num_orb_temp

            end do

            write (stdout, 1010) this % num_target_sym, (this % num_continuum_orbitals_target(i), i = 1, this % num_target_sym)
            1010 format(/,' Number of target symmetries in expansion,   NTGSYM =',I5,/, &
                          ' Number of continuum orbitals for each state, NOTGT =',20I5,/,(' ',20I5))
            write (stdout, 1020) (this % num_ci_target_sym(i), i = 1, this % num_target_sym)
            1020 format(  ' Number of CI components for each state,      NCTGT =',20I10,/,('  ',20I5))
            write (stdout, 1025) (this % num_target_state_sym(i), i = 1, this % num_target_sym)
            1025 format(  ' Number of target states of each symmetry,   NUMTGT =',20I5,/,('  ',20I5))
            write (stdout, 1030) (this % lambda_continuum_orbitals_target(i), i = 1, this % num_target_sym)
            1030 format(  ' Continuum M projection  for each state,      MCONT =',20I5,/,('  ',20I5))
            if (this % gerade_sym_continuum_orbital_target(1) /= 0) then
                write (stdout, 1040) (this % gerade_sym_continuum_orbital_target(i), i = 1, this % num_target_sym)
            end if
            1040 format(  ' Continuum G/U symmetry  for each state,     GUCONT =',20I5,/,('  ',20I5))
            write (stdout, *)

        else if (this % num_continuum > 0) then
            this % num_continuum_orbitals_target(1:this%num_target_sym) = this % temp_orbs(1:this%num_target_sym)
        end if

        this % all_ecp_defined = this % multiplicity_defined .and. this % spatial_defined

        write (stdout, "(' NOBC =',I10,20I4)") (this % num_target_orbitals_sym(I), I = 1, this % num_syms)
        write (stdout, "(' NOB  =',I10,20I4)") (this % num_orbitals_sym_dinf(I), I = 1, this % num_syms)
        write (stdout, "(' NOB0 =',I10,20I4)") (this % num_target_orbitals_sym_dinf_congen(I), I = 1 , this % num_syms)

    end subroutine compute_variables


    subroutine compute_expansions (this)
        class(Options), intent(inout) :: this
        integer                       :: lusme = 0

        if (this % csf_type /= NORMAL_CSF .or. this % ci_target_flag /= NO_CI_TARGET) then

            call mkpt(this % orbital_sequence_number,               & ! kpt
                      this % num_orbitals_sym_dinf,                 & ! nob
                      this % num_target_orbitals_sym,               & ! nobc
                      this % sym_group_flag,                        & ! symtype
                      this % num_continuum_orbitals_target,         & ! notgt
                      this % num_ci_target_sym,                     & ! nctgt
                      this % lambda_continuum_orbitals_target,      & ! mcont
                      this % gerade_sym_continuum_orbital_target,   & ! gucont
                      this % contracted_mat_size,                   & ! mocsf
                      this % num_csfs,                              & ! nocsf0
                      this % actual_num_csfs,                       & ! mxcsf
                      this % last_continuum,                        & ! ncont
                      this % num_expanded_continuum_CSF,            & ! ncont2
                      this % num_target_sym,                        & ! ntgsym
                      this % total_num_target_states,               & ! ntgt
                      this % num_target_state_sym,                  & ! numtgt
                      this % total_num_ci_target,                   & ! ncitot
                      this % maximum_num_target_states,             & ! nummx
                      this % csf_type,                              & ! iexpc
                      this % ci_target_flag,                        & ! icitg
                      lusme)

            1050 format(' SCATCI will expand NOCSF =',i10,' prototype CSFs',/15x, &
                        'into MXCSF =',i10,' actual configurations',/15x, &
                        'and  MOCSF =',i10,' dimension final Hamiltonian')
            write (stdout, 1050) this % num_csfs, this % actual_num_csfs, this % contracted_mat_size

        else

            this % contracted_mat_size = this % num_csfs
            this % actual_num_csfs = this % num_csfs
            this % total_num_target_states = this % num_target_sym

        end if

        write (stdout, "(/' Number of last continuum CSF, NCONT =',i7)") this % last_continuum

        this % seq_num_last_continuum_csf = this % contracted_mat_size - (this % num_csfs - this % last_continuum)
        this % continuum_block_size =  this % seq_num_last_continuum_csf
        if (this % num_eigenpairs == 0) this % num_eigenpairs = this % contracted_mat_size

    end subroutine compute_expansions


    logical function do_ci_contraction (this)
        class(Options), intent(inout) :: this

        do_ci_contraction = (this % ci_target_flag == DO_CI_TARGET_CONTRACT)

    end function do_ci_contraction


    subroutine allocate_integer_array (array, num_elms)
        integer, allocatable :: array(:)
        integer              :: num_elms, err

        if (.not. allocated(array)) then
            allocate(array(num_elms), stat = err)
            call master_memory % track_memory(storage_size(array)/8, size(array), err, 'OPTION::INT_ARRAY')
            if (err /= 0) then
                call mpi_xermsg('Options_module', 'allocate_integer_array', 'Error in allocating one of the arrays', 1, 1)
            end if
            array = 0
        end if

    end subroutine allocate_integer_array


    subroutine allocate_real_array (array, num_elms)
        real(wp), allocatable :: array(:)
        integer              :: num_elms, err

        if (.not. allocated(array)) then
            allocate(array(num_elms), stat = err)
            call master_memory % track_memory(storage_size(array)/8, size(array), err, 'OPTION::REAL_ARRAY')
            if (err /= 0) then
                call mpi_xermsg('Options_module', 'allocate_real_array', 'Error in allocating one of the arrays', 1, 1)
            end if
            array = 0
        end if

    end subroutine allocate_real_array


    !> \brief   Read/write information about Hamiltonian
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Wrapper around SCATCI's ham_pars_io subroutine.
    !>
    !> Reads or writes the file 'ham_data' containing the parameters needed on the call to dgmain.
    !> If `write_ham_pars` is true, then the file is written interpretting the arguments as input.
    !> If `write_ham_pars` is false, then the file is read in interpretting the argumens as output.
    !>
    !> MPI-SCATCI uses only the writing mode, and only when the diagonalization is not undertaken.
    !>
    subroutine ham_pars_io_wrapper (option, write_ham_pars, nelms)

        class(Options), intent(inout) :: option
        logical,        intent(in)    :: write_ham_pars
        integer,        intent(inout) :: nelms
        integer                       :: size_phase

        size_phase = size(option % phase_index)

        if (grid % grank /= master) return

        call ham_pars_io(write_ham_pars,                            &
                         option % seq_num_last_continuum_csf,       &
                         option % num_target_sym,                   &
                         option % num_target_state_sym,             &
                         option % num_continuum_orbitals_target,    &
                         option % lambda_continuum_orbitals_target, &
                         option % hamiltonian_unit,                 &
                         nelms,                                     &
                         option % sym_group_flag,                   &
                         option % spin,                             &
                         option % spin_z,                           &
                         option % num_electrons,                    &
                         option % lambda,                           &
                         option % integral_unit,                    &
                         option % phase_index,                      &
                         size_phase,                                &
                         option % use_UKRMOL_integrals)

     end subroutine ham_pars_io_wrapper

end module Options_module
