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
! ------------------------------------------------------------------------------
!   Modified by Josh Forer <j.forer@posteo.net> 2026:
!     Work around gcc-fortran 16 and -fdefault-integer-8 (DTIO default 8-byte
!     integers are not compatible)  by routing gcc-fortran builds through the
!     existing intel_bug_workaround path
! ------------------------------------------------------------------------------

!> \brief   Further processing of the diagonalization results
!> \authors J Benda
!> \date    2019 - 2023
!>
!> Contains subroutines that can be used to extract further interesting data from the eigenvectors.
!> In particular, this means extraction of boundary amplitudes for outer region codes (see \ref outer_interface)
!> and calculation of multipole moments, mainly for use in the time-dependent R-matrix code (see \ref cdenprop_properties).
!>
!> This module currently supports only the UKRmol+ integral format.
!>
module Postprocessing_module

    use blas_lapack_gbl, only: blasint
    use cdenprop_defs,   only: CIvect
    use precisn,         only: shortint, wp

    implicit none

    private

    public postprocess

    ! data types used by RMT and expected in the 'molecular_data' file
    integer, parameter :: rmt_int  = shortint
    integer, parameter :: rmt_real = wp

    !> \brief   SWINTERF-inspired post-processing object
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> This type extracts and stores outer region channels, as well as boundary amplitudes
    !> evaluated from the solutions.
    !>
    type :: OuterInterface
        integer :: maxchan = 10000  !< Maximal number of channels per irreducible representation. TODO replace by actual number
        integer :: maxnmo           !< Maximal number of continuum molecular orbitals.
        integer :: ntarg            !< Number of outer region states.
        integer :: nchan            !< Number of outer region channels.
        integer :: ismax            !< Largest angular momentum transfer considered (a.k.a. lamax).
        logical :: auto_nvo         !< Automatically determine contracted virtual orbitals.

        integer  :: nfdm            !< Number of radii at which to evaluate the boundary amplitudes.
        real(wp) :: rmatr           !< R-matrix radius.

        real(wp),     allocatable :: r_points(:)  !< Boundary amplitude evaluation radii.
        type(CIvect), allocatable :: wamp(:)      !< Boundary amplitudes (per channel, per eigenstate) for each evaluation radius.

        integer,  allocatable :: idtarg(:)        !< Denprop-to-congen state order permutation.
        integer,  allocatable :: irrchl(:)        !< Irreducible representation of the channel state.
        integer,  allocatable :: ichl(:)          !< Index of target per channel.
        integer,  allocatable :: lchl(:)          !< Projectile angular momentum per channel.
        integer,  allocatable :: mchl(:)          !< Projection of angular momentum per channel.
        integer,  allocatable :: qchl(:)          !< Projectile charge per channel.
        real(wp), allocatable :: echl(:)          !< Projectile energy per channel.
        real(wp), allocatable :: a(:,:)           !< Long-range coupling coefficients.
        integer,  allocatable :: core_charges(:)  !< Number of electrons, on each atom, represented by ECPs.
    contains
        procedure, public  :: init   => init_outer_interface
        procedure, public  :: deinit => deinit_outer_interface
        procedure, public  :: setup_amplitudes

        procedure, public  :: extract_data
        procedure, private :: get_channel_info
        procedure, private :: get_channel_couplings
        procedure, private :: get_boundary_data

        procedure, public  :: write_data
        procedure, private :: write_channel_info
        procedure, private :: write_boundary_data
    end type OuterInterface

contains

    !> \brief   Full post-processing pass
    !> \authors J Benda
    !> \date    2019
    !>
    !> Call the post-processing sequence, starting with evaluation of inner-region dipole moments, followed by
    !> extraction of information useful for the outer region.
    !>
    !> \param[in]  SCATCI_input    Input SCATCI namelist data for all symmetries.
    !> \param[in]  solutions       Eigenvectors and related stuff calculated by the diagonalizer.
    !>
    subroutine postprocess (SCATCI_input, solutions)

        use class_molecular_properties_data, only: molecular_properties_data
        use Options_module,                  only: OptionsSet
        use SolutionHandler_module,          only: SolutionHandler

        type(OptionsSet),                   intent(in)    :: SCATCI_input
        type(SolutionHandler), allocatable, intent(inout) :: solutions(:)

        integer :: iitf

        if (allocated(SCATCI_input % interface_opts)) then
            call redistribute_solutions(SCATCI_input, solutions)
            do iitf = 1, size(SCATCI_input % interface_opts)
                block
                    type(molecular_properties_data) :: properties_all
                    call cdenprop_properties(SCATCI_input, iitf, solutions, properties_all)
                    call outer_interface(SCATCI_input, iitf, solutions, properties_all)
                end block
            end do
        end if

    end subroutine postprocess


    !> \brief   Redistribute solutions from groups to everyone
    !> \authors J Benda
    !> \date    2019
    !>
    !> This subroutine needs to be called before using solutions in other subroutines of this module. The diagonalizations
    !> may have run in MPI groups and so the eigenvectors are not evenly distributed among all processes, which would cause
    !> comlications in communication. This subroutine will redistribute the eigenvectors from groups to the MPI world group.
    !>
    !> \param[in]  SCATCI_input    Input SCATCI namelist data for all symmetries.
    !> \param[in]  solutions       Eigenvectors and related stuff calculated by the diagonalizer.
    !>
    subroutine redistribute_solutions (SCATCI_input, solutions)

        use consts_mpi_ci,          only: PASS_TO_CDENPROP, NO_DIAGONALIZATION
        use mpi_gbl,                only: mpiint, mpi_mod_bcast
        use Options_module,         only: OptionsSet
        use SolutionHandler_module, only: SolutionHandler
        use Parallelization_module, only: grid => process_grid

        type(OptionsSet),                   intent(in)    :: SCATCI_input
        type(SolutionHandler), allocatable, intent(inout) :: solutions(:)

        type(SolutionHandler) :: oldsolution
        integer(mpiint)       :: srcrank
        integer               :: i, nnuc, nstat

        ! nothing to do when the solutions are already uniformly distributed everywhere
        if (grid % sequential) return

        ! set up the uniformly distributed solutions
        do i = 1, size(SCATCI_input % opts)
            if (iand(SCATCI_input % opts(i) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                SCATCI_input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then

                ! 1. Share metadata for the solutions from the master of group that found that solution to everyone.

                srcrank = grid % group_master_world_rank(grid % which_group_is_work(i))

                call mpi_mod_bcast(solutions(i) % core_energy,  srcrank)
                call mpi_mod_bcast(solutions(i) % vec_dimen,    srcrank)

                call mpi_mod_bcast(solutions(i) % ci_vec % nstat,           srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % mgvn,            srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % s,               srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % sz,              srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % nnuc,            srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % e0,              srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % CV_is_scalapack, srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % mat_dimen_r,     srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % mat_dimen_c,     srcrank)

                nnuc  = solutions(i) % ci_vec % nnuc
                nstat = solutions(i) % ci_vec % nstat

                if (.not. allocated(solutions(i) % ci_vec % ei)) then
                    allocate (solutions(i) % ci_vec % ei(nstat))
                end if

                call mpi_mod_bcast(solutions(i) % ci_vec % cname(1:nnuc),   srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % charge(1:nnuc),  srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % xnuc(1:nnuc),    srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % ynuc(1:nnuc),    srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % znuc(1:nnuc),    srcrank)
                call mpi_mod_bcast(solutions(i) % ci_vec % ei(1:nstat),     srcrank)

                ! 2. Backup the original solutions, which are distributed over the MPI group.

                call oldsolution % destroy()
                oldsolution = solutions(i)

                if (.not. grid % is_my_group_work(i)) then
                    oldsolution % ci_vec % blacs_context = -1
                    oldsolution % ci_vec % descr_CV_mat = -1

                    ! Initialize dummy submatrices for ranks that are not members of the original group.
                    ! This is needed for formal correctness; gfortran's -fcheck=all otherwise chokes
                    ! when passing these submatrices to pdgemr2d in the redistribution below.
                    allocate (oldsolution % ci_vec % cv(0, 0))
                end if

                ! 3. Redistribute the eigenvectors over whole MPI world.

                solutions(i) % ci_vec % blacs_context   = grid % wcntxt
                solutions(i) % ci_vec % nprow           = grid % wprows
                solutions(i) % ci_vec % npcol           = grid % wpcols

                call solutions(i) % ci_vec % init_CV(int(oldsolution % ci_vec % mat_dimen_r), &
                                                     int(oldsolution % ci_vec % mat_dimen_c))
                call solutions(i) % ci_vec % redistribute(oldsolution % ci_vec)

            end if
        end do

    end subroutine redistribute_solutions


    !> \brief   Evaluate multipoles in CDENPROP
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> Calculate properties (multipole moments) of the solutions using CDENPROP library. Only those symmetries
    !> that have `igh = 1, vecstore = 3` set will be used. This subroutine also assumes that setups for all
    !> symmetries are compatible, in particular:
    !>   - use the same molecular integral format
    !>   - use the same CI targets (in case of scattering diagonalization)
    !>
    !> \param[in]  SCATCI_input    Input SCATCI namelist data for all symmetries.
    !> \param[in]  iitf            Which outer_interface namelist to consider.
    !> \param[in]  solutions       Eigenvectors and related stuff calculated by the diagonalizer.
    !> \param[out] properties_all  Table of properties calculated by CDENPROP.
    !>
    subroutine cdenprop_properties (SCATCI_input, iitf, solutions, properties_all)

        use class_molecular_properties_data, only: molecular_properties_data
        use class_namelists,                 only: cdenprop_namelists
        use consts_mpi_ci,                   only: PASS_TO_CDENPROP, NO_DIAGONALIZATION, NO_CI_TARGET
        use cdenprop_procs,                  only: cdenprop_drv
        use cdenprop_defs,                   only: ir_max
        use mpi_gbl,                         only: master
        use Options_module,                  only: OptionsSet
        use SolutionHandler_module,          only: SolutionHandler
        use Timing_module,                   only: master_timer

        type(OptionsSet),                   intent(in)    :: SCATCI_input
        integer,                            intent(in)    :: iitf
        type(SolutionHandler), allocatable, intent(inout) :: solutions(:)
        type(molecular_properties_data),    intent(inout) :: properties_all

        type(cdenprop_namelists) :: namelist_ij
        type(CIvect)             :: oldblock

        integer :: i, j, ntgsym, ninputs, maxpole, nblock, mxstat

        ninputs = size(SCATCI_input % opts)
        maxpole = 2

        ! for now, skip evaluation of properties if RMT data not required
        if (.not. SCATCI_input % interface_opts(iitf) % write_dip .and. &
            .not. SCATCI_input % interface_opts(iitf) % write_rmt) return

        ! prepare free slots for the distributed dense property matrices
        call properties_all % clean
        call properties_all % preallocate_property_blocks(ninputs**2 * (maxpole + 1)**2)
        properties_all % dense_blocks(:) % mat_dimen = 0

        ! calculate properties for all pairs of symmetries
        sym_i_loop: do i = 1, ninputs
            if (iand(SCATCI_input % opts(i) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                SCATCI_input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then

                sym_j_loop: do j = i, ninputs
                    if (iand(SCATCI_input % opts(j) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                        SCATCI_input % opts(j) % diagonalization_flag /= NO_DIAGONALIZATION) then

                        ! only calculate self-overlaps when required
                        if (i == j .and. .not. SCATCI_input % interface_opts(iitf) % all_props) cycle

                        ! virtual CDENPROP namelist
                        call namelist_ij % init(2)

                        ! assign values obtained from SCATCI namelists
                        namelist_ij % lucsf(1:2)   = (/ SCATCI_input % opts(i) % megul, &
                                                        SCATCI_input % opts(j) % megul /)
                        namelist_ij % nciset(1:2)  = (/ SCATCI_input % opts(i) % output_ci_set_number, &
                                                        SCATCI_input % opts(j) % output_ci_set_number /)
                        namelist_ij % nstat        = 0      ! use all states
                        namelist_ij % lupropw      = 0      ! no output needed now, will be done later
                        namelist_ij % lucivec      = 0      ! no CI units needed, will be given as arguments
                        namelist_ij % numtgt       = 0      ! clear target states; subset of the array populated below
                        namelist_ij % lutdip       = SCATCI_input % interface_opts(iitf) % lutdip
                        namelist_ij % ukrmolp_ints = SCATCI_input % opts(i) % use_UKRMOL_integrals
                        namelist_ij % luprop       = merge(SCATCI_input % opts(i) % integral_unit, &
                                                           SCATCI_input % interface_opts(iitf) % luprop, &
                                                           SCATCI_input % interface_opts(iitf) % luprop < 0)
                        namelist_ij % isw          = SCATCI_input % interface_opts(iitf) % isw

                        ntgsym = SCATCI_input % opts(i) % num_target_sym

                        ! following depend on whether this is target or scattering calculation
                        if (SCATCI_input % opts(i) % ci_target_flag == NO_CI_TARGET) then
                            namelist_ij % max_multipole    = 2  ! dipole + quadrupole
                        else
                            namelist_ij % max_multipole    = 1  ! dipole only
                            namelist_ij % numtgt(1:ntgsym) = SCATCI_input % opts(i) % num_target_state_sym(1:ntgsym)
                            namelist_ij % itarget_spin(1:ntgsym) = SCATCI_input % interface_opts(iitf) % itspin(1:ntgsym)
                        end if

                        ! use BLACS context of the first vector for communication (they are all of the same shape anyway)
                        solutions(j) % ci_vec % blacs_context   = solutions(i) % ci_vec % blacs_context
                        solutions(j) % ci_vec % descr_CV_mat(2) = solutions(i) % ci_vec % descr_CV_mat(2)

                        ! get the properties in CDENPROP
                        call master_timer % start_timer("Cdenprop")
                        call cdenprop_drv(solutions(i) % ci_vec, solutions(j) % ci_vec, namelist_ij, properties_all)
                        call master_timer % stop_timer("Cdenprop")
                        call master_timer % report_timers

                        ! cleanup
                        call namelist_ij % dealloc

                    end if
                end do sym_j_loop

            end if
        end do sym_i_loop

        ! prefer SWINTERF-compatible format of the property file (so that it can be, e.g., read back to CDENPROP_ALL)
        properties_all % swintf_format = .true.

        ! write calculated properties
        if (SCATCI_input % interface_opts(iitf) % write_dip .and. properties_all % no_symmetries > 0) then
            if (SCATCI_input % interface_opts(iitf) % all_props) call properties_all % sort_by_energy
            call properties_all % write_properties(SCATCI_input % interface_opts(iitf) % lupropw, &
                                                   SCATCI_input % opts(1) % use_UKRMOL_integrals)
        end if

        ! inflate all property matrices to uniform size for easy output for RMT, which wants it that way
        if (SCATCI_input % interface_opts(iitf) % write_rmt) then
            nblock = properties_all % no_blocks
            mxstat = max(maxval(properties_all % dense_blocks(1:nblock) % mat_dimen_r), &
                         maxval(properties_all % dense_blocks(1:nblock) % mat_dimen_c))
            do i = 1, nblock
                if (properties_all % dense_blocks(i) % mat_dimen_r /= mxstat .or. &
                    properties_all % dense_blocks(i) % mat_dimen_c /= mxstat) then
                        oldblock = properties_all % dense_blocks(i)
                        call properties_all % dense_blocks(i) % init_CV(mxstat, mxstat)
                        properties_all % dense_blocks(i) % CV = 0
                        call properties_all % dense_blocks(i) % redistribute(oldblock)
                        call oldblock % final_CV
                end if
                properties_all % dense_blocks(i) % mat_dimen = max(properties_all % dense_blocks(i) % mat_dimen_r, &
                                                                   properties_all % dense_blocks(i) % mat_dimen_c)
            end do
        end if

    end subroutine cdenprop_properties


    !> \brief   Extract data needed by outer-region codes
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> This subroutine does the following:
    !>  - write channel information to unit LUCHAN
    !>  - write boundary amplitudes to unit LURMT
    !>  - write RMT data to file `molecular_data`
    !>
    !> Essentially, it should cover the work of SWINTERF and RMT_INTERFACE. At the moment, only the symmetries
    !> marked with `vecstore` with the second least-significant bit set (i.e., 2, 5 or 7) are processed in the outer interface.
    !>
    !> The outer interface has an automatic support for the contracted virtual orbitals. When the number of continuum orbitals
    !> specified in GBTOlib (scatci_integrals) is different than the number of CI target continuum orbitals specified in the input
    !> for MPI-SCATCI, the difference (i.e. target orbitals included as fake continuum) are automatically recognized as contracted
    !> virtual orbitals. This is equivalent to the NVO infrastructure in SWINTERF. If more fine-grained control over calculation
    !> of NVO is desired, use the legacy SCATCI with SWINTERF instead.
    !>
    !> \param[in]  SCATCI_input      Input SCATCI namelist data for all symmetries.
    !> \param[in]  iitf              Which outer_interface namelist to use.
    !> \param[in]  solutions         Eigenvectors and related stuff calculated by the diagonalizer.
    !> \param[out] inner_properties  Table of (N+1) properties calculated by CDENPROP.
    !>
    subroutine outer_interface (SCATCI_input, iitf, solutions, inner_properties)

        use class_molecular_properties_data, only: molecular_properties_data
        use consts_mpi_ci,                   only: PASS_TO_CDENPROP, NO_DIAGONALIZATION
        use Options_module,                  only: OptionsSet
        use SolutionHandler_module,          only: SolutionHandler
        use Timing_module,                   only: master_timer

        type(OptionsSet),                   intent(in) :: SCATCI_input
        integer,                            intent(in) :: iitf
        type(SolutionHandler), allocatable, intent(in) :: solutions(:)
        type(molecular_properties_data),    intent(in) :: inner_properties

        type(molecular_properties_data)   :: target_properties
        type(OuterInterface), allocatable :: intf(:)

        integer :: i, nsym, nset

        nset = 0                          ! number of sets in the output files fort.10 and fort.21
        nsym = size(SCATCI_input % opts)  ! number of inputs (symmetries) that we dealt with

        if (.not. SCATCI_input % interface_opts(iitf) % write_amp .and. &
            .not. SCATCI_input % interface_opts(iitf) % write_rmt) return

        allocate (intf(nsym))

        call target_properties % read_properties(SCATCI_input % interface_opts(iitf) % lutarg, 1)

        call master_timer % start_timer("Outer interface")

        ! open the outer interface files for writing
        if (SCATCI_input % interface_opts(iitf) % cform == 'F') then
            open (unit = SCATCI_input % interface_opts(iitf) % luchan, form = 'formatted', action = 'write')
        else
            open (unit = SCATCI_input % interface_opts(iitf) % luchan, form = 'unformatted', action = 'write')
        end if
        if (SCATCI_input % interface_opts(iitf) % rform == 'F') then
            open (unit = SCATCI_input % interface_opts(iitf) % lurmt, form = 'formatted', action = 'write')
        else
            open (unit = SCATCI_input % interface_opts(iitf) % lurmt, form = 'unformatted', action = 'write')
        end if

        ! run SWINTERF for all symmetries
        do i = 1, nsym
            call intf(i) % init(SCATCI_input % interface_opts(iitf))
            if (iand(SCATCI_input % opts(i) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                SCATCI_input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then
                nset = nset + 1
                call intf(i) % setup_amplitudes(SCATCI_input % opts(i))
                call intf(i) % extract_data(SCATCI_input % opts(i), target_properties, solutions(i))
                if (SCATCI_input % interface_opts(iitf) % write_amp) then
                    call intf(i) % write_data(i, SCATCI_input, iitf, target_properties, solutions(i))
                end if
            end if
        end do

        close (unit = SCATCI_input % interface_opts(iitf) % luchan)
        close (unit = SCATCI_input % interface_opts(iitf) % lurmt)

        if (SCATCI_input % interface_opts(iitf) % write_rmt) then
            call write_rmt_data(SCATCI_input, iitf, inner_properties, solutions, intf)
        end if

        call master_timer % stop_timer("Outer interface")

    end subroutine outer_interface


    !> \brief   Allocate memory for channel information
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> \param[in] this  Interface object to initialize.
    !> \param[in] input Input SCATCI namelist data for all symmetries.
    !>
    subroutine init_outer_interface (this, input)

        use mpi_gbl,        only: mpi_xermsg
        use Options_module, only: InterfaceOptions

        class(OuterInterface),  intent(inout) :: this
        type(InterfaceOptions), intent(in)    :: input

        integer :: i, ierr

        this % ntarg  = 0
        this % nchan  = 0
        this % ismax  = 0
        this % maxnmo = 0

        this % nfdm     = merge(input % nfdm, 0, input % write_rmt)
        this % rmatr    = input % rmatr
        this % auto_nvo = input % auto_nvo

        allocate (this % r_points(this % nfdm + 1))

        do i = 1, this % nfdm + 1
            this % r_points(i) = this % rmatr - (this % nfdm + 1 - i) * input % delta_r
        end do

        allocate (this % irrchl(this % maxchan), &
                  this % ichl(this % maxchan),   &
                  this % lchl(this % maxchan),   &
                  this % mchl(this % maxchan),   &
                  this % qchl(this % maxchan),   &
                  this % echl(this % maxchan), stat = ierr)

        if (ierr /= 0) then
            call mpi_xermsg('OuterInterface', 'init_outer_interface', 'Failed to allocate memory for channel data.', 1, 1)
        end if

        if (allocated(input % idtarg)) then
            this % ntarg  = size(input % idtarg)
            this % idtarg = input % idtarg
        end if

    end subroutine init_outer_interface


    !> \brief   Release memory held by this object
    !> \authors J Benda
    !> \date    2019
    !>
    !> \param[in] this  Interface object to finalize.
    !>
    subroutine deinit_outer_interface (this)
#if defined(usempi) && defined(scalapack)
        use cdenprop_defs, only: blacs_gridexit
#endif
        class(OuterInterface), intent(inout) :: this

        integer :: i

        ! release the no-longer-needed BLACS context
        if (allocated(this % wamp)) then
            do i = 1, size(this % wamp)
                call this % wamp(i) % final_CV
#if defined(usempi) && defined(scalapack)
                if (i == 1 .and. this % wamp(this % nfdm + 1) % blacs_context >= 0) then
                    call blacs_gridexit(this % wamp(this % nfdm + 1) % blacs_context)
                end if
#endif
            end do
            deallocate (this % wamp)
        end if

        if (allocated(this % irrchl))   deallocate (this % irrchl)
        if (allocated(this % ichl))     deallocate (this % ichl)
        if (allocated(this % lchl))     deallocate (this % lchl)
        if (allocated(this % mchl))     deallocate (this % mchl)
        if (allocated(this % qchl))     deallocate (this % qchl)
        if (allocated(this % echl))     deallocate (this % echl)
        if (allocated(this % a))        deallocate (this % a)
        if (allocated(this % r_points)) deallocate (this % r_points)

        this % nchan  = 0
        this % ismax  = 0
        this % maxnmo = 0
        this % nfdm   = 0

    end subroutine deinit_outer_interface


    !> \brief   Initialize raw boundary amplitudes (in integral library)
    !> \authors J Benda
    !> \date    2019
    !>
    !> \param[in] this  Interface object to read.
    !> \param[in] opts  SCATCI input namelist data.
    !>
    subroutine setup_amplitudes (this, opts)

        use Options_module,       only: Options
        use ukrmol_interface_gbl, only: read_ukrmolp_basis, eval_amplitudes, get_ecp_core_charges

        class(OuterInterface), intent(inout) :: this
        type(Options),         intent(in)    :: opts

        if (opts % use_UKRMOL_integrals) then
            call read_ukrmolp_basis(opts % integral_unit)
            call eval_amplitudes(this % rmatr, .true.)
            call get_ecp_core_charges(this % core_charges)
        end if

    end subroutine setup_amplitudes


    !> \brief   Get interface data
    !> \authors J Benda
    !> \date    2019 - 2023
    !>
    !> This mirros the functionality of SWITERF (UKRmol-out), particularly:
    !>  - reads the target property file (*fort.24*)
    !>  - assembles the list of outer region channels
    !>  - evaluates amplitudes of the eigenstates at R-matrix boundary and, if requested, in its vicinity for use in RMT
    !>
    !> All data are then stored internally for writing.
    !>
    !> \param[in] this      Interface object to read.
    !> \param[in] opts      SCATCI input namelist data.
    !> \param[in] prop      Target properties.
    !> \param[in] solution  Eigenvectors and related stuff calculated by the diagonalizer.
    !>
    subroutine extract_data (this, opts, prop, solution)

        use class_molecular_properties_data, only: molecular_properties_data
        use mpi_gbl,                         only: master, myrank
        use Options_module,                  only: Options
        use SolutionHandler_module,          only: SolutionHandler

        class(OuterInterface),           intent(inout) :: this
        type(molecular_properties_data), intent(in)    :: prop
        type(Options),                   intent(in)    :: opts
        type(SolutionHandler),           intent(in)    :: solution

        integer,  parameter :: ismax  = 2
        real(wp), parameter :: alpha0 = 0
        real(wp), parameter :: alpha2 = 0

        logical :: use_pol
        integer :: lamax

        lamax = ismax
        use_pol = .false.  ! if set to .true., then channel_couplings will take into account the polarizability

        if (alpha0 /= 0 .or. alpha2 /= 0) then
           use_pol = .true.
           lamax = max(ismax, 3)  ! polarizability (1/r^4) corresponds to lambda = 3
        end if

        this % ismax = lamax

        call this % get_channel_info(opts, prop)
        call this % get_boundary_data(opts, prop, solution)

        ! only the master process writes the outer region potential coefficients; other do not need to have them
        if (myrank == master) then
            call this % get_channel_couplings(lamax, prop % no_states, prop, alpha0, alpha2, use_pol)
        end if

        ! redefine magnetic quantum number by multiplying in the projectile charge
        this % mchl(1:this % nchan) = this % mchl(1:this % nchan) * this % qchl(1:this % nchan)

    end subroutine extract_data


    !> \brief   Write interface data
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> Updates the channel and R-matrix dump files (fort.10 and fort.21, respectively) by adding/overwriting a data set.
    !>
    !> \param[in] this      Interface object to read.
    !> \param[in] isym      Output set number for this irreducible representation.
    !> \param[in] opts      SCATCI input namelist data.
    !> \param[in] iitf      Which outer_interface namelist to use.
    !> \param[in] prop      Target properties.
    !> \param[in] solution  Eigenvectors and related stuff calculated by the diagonalizer.
    !>
    subroutine write_data (this, isym, opts, iitf, prop, solution)

        use class_molecular_properties_data, only: molecular_properties_data
        use Options_module,                  only: OptionsSet
        use SolutionHandler_module,          only: SolutionHandler

        class(OuterInterface),           intent(in) :: this
        type(molecular_properties_data), intent(in) :: prop
        type(OptionsSet),                intent(in) :: opts
        type(SolutionHandler),           intent(in) :: solution
        integer,                         intent(in) :: isym, iitf

        call this % write_channel_info(isym, &
                                       opts % opts(isym), &
                                       opts % interface_opts(iitf) % luchan, &
                                       opts % interface_opts(iitf) % cform, &
                                       prop)

        call this % write_boundary_data(isym, &
                                        opts % opts(isym), &
                                        opts % interface_opts(iitf) % lurmt, &
                                        opts % interface_opts(iitf) % rform, &
                                        prop, &
                                        solution)

    end subroutine write_data


    !> \brief   Assemble the list of outer channels
    !> \authors J Benda
    !> \date    2019
    !>
    !> Adapted from SWCHANL (UKRmol-out). Reads the target information from `fort.24` and, with the knowledge
    !> of the (N+1)-electron setup, determines the individual channels.
    !>
    !> \param[inout] this  Interface object to update.
    !> \param[in]    opts  SCATCI namelist information for this irreducible representation.
    !> \param[in]    prop  Target properties.
    !>
    subroutine get_channel_info (this, opts, prop)

        use algorithms,                      only: findloc
        use class_molecular_properties_data, only: molecular_properties_data
        use const_gbl,                       only: stdout
        use global_utils,                    only: mprod
        use mpi_gbl,                         only: mpi_xermsg
        use Options_module,                  only: Options
        use precisn,                         only: wp
        use ukrmol_interface_gbl,            only: ukp_preamp

        class(OuterInterface),           intent(inout) :: this
        type(Options),                   intent(in)    :: opts
        type(molecular_properties_data), intent(in)    :: prop

        integer  :: i, j, k, mgvn, spin, tgt_symm, tgt_mgvn, tgt_spin, nch, prj_mgvn, congen_tgt_id, denprop_tgt_id, isym
        real(wp) :: tgt_enrg, ebase

        ebase  = prop % energies(1)
        mgvn   = opts % lambda
        spin   = nint(2 * opts % spin + 1)

        this % maxnmo = sum(opts % num_orbitals_sym)
        this % nchan  = 0

        ! construct default IDTARG, i.e. denprop (energy-order) index for each congen (symmetry-order) target state
        if (.not. allocated(this % idtarg)) then
            this % ntarg = opts % total_num_target_states
            allocate (this % idtarg(this % ntarg))
            this % idtarg = -1
            congen_tgt_id = 0
            congen_symm_loop: do j = 1, size(opts % num_target_state_sym)
                congen_targ_loop: do k = 1, opts % num_target_state_sym(j)
                    congen_tgt_id = congen_tgt_id + 1
                    denprop_tgt_id = 0
                    tgt_mgvn = opts % target_spatial(congen_tgt_id)
                    tgt_spin = opts % target_multiplicity(congen_tgt_id)
                    denprop_targ_loop: do i = 1, prop % no_states
                        isym = prop % energies_index(2, i)
                        if (prop % symmetries_index(1, isym) == tgt_mgvn .and. &
                            prop % symmetries_index(2, isym) == tgt_spin) then
                            denprop_tgt_id = denprop_tgt_id + 1
                        end if
                        if (denprop_tgt_id == k) then
                            this % idtarg(congen_tgt_id) = i
                            exit denprop_targ_loop
                        end if
                    end do denprop_targ_loop
                end do congen_targ_loop
            end do congen_symm_loop
            if (any(this % idtarg == -1)) then
                call mpi_xermsg('Postprocessing_module', 'get_channel_info', 'Property file not compatible with target symmetries&
                                &given in the SCATCI input file.', 1, 1)
            end if
        end if

        write (stdout, '(/,1x,"Target state indices sorted by energy (IDTARG): ")')
        write (stdout, '(  1x,20I4)') this % idtarg

        do i = 1, prop % no_states

            ! skip states not referenced in IDTARG
            if (findloc(this % idtarg, i, 1) == 0) cycle

            ! retrieve target data from the properties structure
            tgt_symm = prop % energies_index(2, i)
            tgt_mgvn = prop % symmetries_index(1, tgt_symm)
            tgt_spin = prop % symmetries_index(2, tgt_symm)
            tgt_enrg = prop % energies(i)

            ! angular momentum composition rule -- must allow for projectile spin
            if (abs(spin - tgt_spin) /= 1) cycle

            ! get projectile total symmetry IRR from the multiplication table
            prj_mgvn = mprod(tgt_mgvn + 1, mgvn + 1, 0, stdout)

            ! retrieve channel descriptions from the integral library
            if (opts % use_UKRMOL_integrals) then
                call ukp_preamp(prj_mgvn, this % nchan + 1, this % lchl, this % mchl, this % qchl, nch, this % maxnmo)
            else
                call mpi_xermsg('Postprocessing_module', 'get_channel_info', &
                                'Interface for SWINTERF integrals not implemented.', 1, 1)
            end if

            ! fill in the remaining channel identifiers
            do j = 1, nch
                this % ichl(this % nchan + j) = i
                this % echl(this % nchan + j) = 2 * (tgt_enrg - ebase)
                this % irrchl(this % nchan + j) = prj_mgvn
            end do

            ! update total number of channels
            this % nchan = this % nchan + nch

        end do

        write (stdout, '(/," Channel Target  l  m  q       Irr  Energy   Energy(eV)")')

        do i = 1, this % nchan
            write (stdout, '(I5,I8,I5,2I3,2x,A4,I3,2F10.6)') i, this % ichl(i), this % lchl(i), this % mchl(i), this % qchl(i), &
                                                             'N/A', this % irrchl(i), this % echl(i), this % echl(i) / 0.073500D0
        end do

    end subroutine get_channel_info


    !> \brief   Evaluate boundary amplitudes for propagation and RMT
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> Calculates the contribution to each channel amplitude per eigenstate. The master process retrieves the raw boundary
    !> amplitudes from the integral library, then distributes them in the BLACS context of the solution vector, where the
    !> multiplication takes place.
    !>
    !> The boundary amplitudes evaluated at the R_matrix sphere (but not those evaluated inside) are then collected to master
    !> process, who will need to write them to the record-based swinterf output file.
    !>
    !> \todo This subroutine contains a potential hidden memory scaling obstacle, which is the evaluation of the boundary
    !>       amplitudes in GBTOlib. Even though GBTOlib internally stores just a compact array of orbital boundary amplitudes
    !>       (dimension Npw × Norb), its interface only allows retrieval of the channel boundary amplitudes (Nchan × Norb).
    !>       For calculations with large numbers of channels (= target states), this is a significantly larger matrix, which
    !>       should be distributed. However, this is not possible with the current GBTOlib interface.
    !>
    !> \param[in]  this       Interface object to update.
    !> \param[in]  opts       SCATCI namelist information for this irreducible representation.
    !> \param[in]  prop       Results from denprop.
    !> \param[in]  solution   Eigenvectors and related stuff calculated by the diagonalizer.
    !>
    subroutine get_boundary_data (this, opts, prop, solution)

        use algorithms,                      only: findloc
        use class_molecular_properties_data, only: molecular_properties_data
#if defined(usempi) && defined(scalapack)
        use cdenprop_defs,                   only: blacs_get, blacs_gridinit
#endif
        use const_gbl,                       only: stdout
        use global_utils,                    only: mprod
        use mpi_gbl,                         only: master, myrank, mpi_xermsg
        use Options_module,                  only: Options
        use SolutionHandler_module,          only: SolutionHandler
        use ukrmol_interface_gbl,            only: eval_amplitudes, ukp_readamp

        class(OuterInterface),           intent(inout)  :: this
        type(Options),                   intent(in)     :: opts
        type(molecular_properties_data), intent(in)     :: prop
        type(SolutionHandler),           intent(in)     :: solution

        type(CIvect) :: orb_amps, dist_orb_amps

        integer,  allocatable :: idchl(:), nvo(:), tgcontirr(:)
        real(wp), allocatable :: ampls(:)

        integer :: i, j, k, ir, ncontmo(8), iprnt = 0, nchan, nstat, norbs, morbs, ncnt, offs, totirr, tgtirr

        nchan = this % nchan
        nstat = solution % ci_vec % nstat
        morbs = maxval(opts % num_orbitals_sym)  ! safe upper bound on maximal number of continuum orbitals per channel...
        norbs = dot_product(opts % num_target_state_sym, opts % num_continuum_orbitals_target)  ! ... and total

        allocate(ampls(morbs), idchl(nchan), tgcontirr(this % ntarg), nvo(opts % num_syms))

        ! set up energy order channel map
        k = 0
        do j = 1, this % ntarg
            do i = 1, nchan
                if (this % ichl(i) == this % idtarg(j)) then
                    k = k + 1
                    idchl(k) = i
                end if
            end do
            totirr = opts % lambda + 1
            tgtirr = prop % symmetries_index(1, prop % energies_index(2, j)) + 1
            tgcontirr(j) = mprod(totirr, tgtirr, 0, stdout)
        end do

        write (stdout, '(/,1x,"Channel target state energy map (IDCHL): ")')
        write (stdout, '(  1x,20I5)') idchl

        ! set up storage for the raw boundary amplitudes (master-only ScaLAPACK matrix)
        orb_amps % nprow = 1
        orb_amps % npcol = 1
        orb_amps % blacs_context   = -1
        orb_amps % CV_is_scalapack = solution % ci_vec % CV_is_scalapack
#if defined(usempi) && defined(scalapack)
        if (orb_amps % CV_is_scalapack) then
            call blacs_get(-1_blasint, 0_blasint, orb_amps % blacs_context)
            call blacs_gridinit(orb_amps % blacs_context, 'R', 1_blasint, 1_blasint)
        end if
#endif
        call orb_amps % init_CV(nchan, norbs)   ! WARNING: this can be a fairly big matrix, but cannot be distributed

        ! set up storage for the evaluated boundary amplitudes
        if (.not. allocated(this % wamp)) then
            allocate (this % wamp(this % nfdm + 1))
        end if

        ! set up cyclic distribution descriptors
        this % wamp(:) % nprow = solution % ci_vec % nprow
        this % wamp(:) % npcol = solution % ci_vec % npcol
        this % wamp(:) % blacs_context = solution % ci_vec % blacs_context
        this % wamp(:) % CV_is_scalapack = solution % ci_vec % CV_is_scalapack

        ! allocate the matrices
        do ir = 1, this % nfdm + 1
            call this % wamp(ir) % init_CV(nchan, nstat)
        end do

        ! for all evaluation radii
        do ir = 1, this % nfdm + 1

            ! get raw boundary amplitudes of the continuum orbitals
            if (myrank == master) then
                ! retrieve orbital amplitudes from the integral library
                call eval_amplitudes(this % r_points(ir), .false.)
                call ukp_readamp(orb_amps % CV, this % nchan, this % irrchl, this % lchl, this % mchl, this % qchl, ncontmo, &
                                 opts % lambda_continuum_orbitals_target, iprnt)

                ! calculate number of contracted virtual orbitals (i.e. orbitals that are not continuum orbitals according
                ! to GBTOlib, but they are specified as continuum in SCATCI anyway)
                if (this % auto_nvo) then
                    nvo(1 : opts % num_syms) = opts % num_orbitals_sym(1 : opts % num_syms) &
                                             - opts % num_target_orbitals_sym(1 : opts % num_syms) &
                                             - ncontmo(1 : opts % num_syms)
                else
                    nvo = 0
                end if

                if (ir == 1) then
                    write (stdout, '(/,1x,"Virtual orbitals per CI target continuum spin-symmetry: ")')
                    write (stdout, '(*(I4))') pack(nvo, [ (findloc(tgcontirr, i, 1) /= 0, i = 1, opts % num_syms) ])
                    write (stdout, '(/,1x,"Virtual orbitals per CI target state continuum in energy order (NVO): ")')
                    write (stdout, '(20I4)') nvo(tgcontirr)
                end if

                if (any(nvo(tgcontirr) < 0)) then
                    call mpi_xermsg('Postprocessing_module', 'get_boundary_data', &
                                    'Some NVOs are negative. This looks like incorrect setup of continuum orbitals.', 1, 1)
                end if

                ! re-align the orbital amplitude data to their absolute positions in the list of all virt+cont orbitals:
                !
                !              1st target used in channels...........2nd target used in channels...........etc....
                !    orb_amps: [virtual orbitals][continuum orbitals][virtual orbitals][continuum orbitals]etc...
                do i = 1, nchan
                    j = idchl(i)            ! channel index in CONGEN ordering of target states
                    k = this % irrchl(j)    ! partial wave irreducible representation (1, 2, ...)
                    if (i == 1) then
                        ! initial offset for the first channel: place evaluated continuum orbital boundary amplitudes
                        ! at positions after the last virtual orbital
                        offs = nvo(k)
                    else if (this % ichl(j) /= this % ichl(idchl(i - 1))) then
                        ! if the target changed between the preceding and the current channel, increment the offset by the number
                        ! of continuum orbitals of the previous state and the number of virtual orbitals of the current state
                        offs = offs + ncnt + nvo(k)
                    end if
                    ncnt = ncontmo(k)
                    ampls(1:ncnt) = orb_amps % CV(j, 1:ncnt)
                    orb_amps % CV(j, 1:norbs) = 0
                    orb_amps % CV(j, offs + 1:offs + ncnt) = ampls(1:ncnt)
                end do
            end if

            ! multiply the raw boundary amplitudes with the eigenvectors
            if (.not. solution % ci_vec % CV_is_scalapack) then
                call this % wamp(ir) % A_B_matmul(orb_amps, solution % ci_vec, 'N', 'N')
            else
                ! distribute raw boundary amplitudes over solution context
                dist_orb_amps % nprow = solution % ci_vec % nprow
                dist_orb_amps % npcol = solution % ci_vec % npcol
                dist_orb_amps % blacs_context   = solution % ci_vec % blacs_context
                dist_orb_amps % CV_is_scalapack = solution % ci_vec % CV_is_scalapack
                call dist_orb_amps % init_CV(nchan, norbs)
                call dist_orb_amps % redistribute(orb_amps)

                ! multiply distributed matrices
                call this % wamp(ir) % A_B_matmul(dist_orb_amps, solution % ci_vec, 'N', 'N')
            end if

            ! change normalization convention to the RMT one (but not for the sufrace boundaries)
            if (ir < this % nfdm + 1) then
                this % wamp(ir) % CV = sqrt(2.0_wp) * this % wamp(ir) % CV
            end if

        end do

    end subroutine get_boundary_data


    !> \brief   Write the channel list to disk
    !> \authors J Benda
    !> \date    2019 - 2023
    !>
    !> Adapted from WRITCH (UKRmol-out). Writes the information into the selected set
    !> of the channel output file *fort.10*.
    !>
    !> \param[in]  this       Interface object to update.
    !> \param[in]  nchset     Output data set in the file unit.
    !> \param[in]  opts       SCATCI namelist information for this irreducible representation.
    !> \param[in]  luchan     Unit number for channel info output.
    !> \param[in]  cform      'F'/'U' for formatted or unformatted output.
    !> \param[in]  prop       Results from denprop.
    !>
    subroutine write_channel_info (this, nchset, opts, luchan, cform, prop)

        use algorithms,                      only: findloc
        use class_molecular_properties_data, only: molecular_properties_data
        use consts,                          only: amu
        use consts_mpi_ci,                   only: SYMTYPE_D2H
        use mpi_gbl,                         only: myrank, master
        use Options_module,                  only: Options

        class(OuterInterface), intent(in) :: this
        integer,               intent(in) :: nchset, luchan
        type(Options),         intent(in) :: opts
        character(len=1),      intent(in) :: cform

        type(molecular_properties_data), intent(in) :: prop

        integer, parameter :: keych = 10

        integer :: i, k, nvd, nrec, ninfo, ndata, nvib, ndis, stot, gutot, ntarg, ivtarg(1), iv(1), ion, mtarg, starg, gutarg, nnuc
        integer :: isymm
        real(wp) :: r, rmass, etarg

        if (myrank /= master) return    ! the file is written by rank 0 only

        nnuc  = prop % no_nuclei
        nvib  = 0
        ndis  = 0
        gutot = SYMTYPE_D2H  ! to indicate polyatomic code

        if (opts % use_UKRMOL_integrals) then
           ion   = nint(sum(prop % charge(1:nnuc) - this % core_charges(1:nnuc))) - (opts % num_electrons - 1)
        else
           ion   = nint(sum(prop % charge(1:nnuc))) - (opts % num_electrons - 1)
        endif

        r     = 0
        ntarg = this % ntarg
        stot  = nint(2 * opts % spin + 1)
        rmass = sum(1 / prop % mass(1:nnuc), prop % mass(1:nnuc) > 0)

        if (rmass > 0) then
            rmass = amu / rmass
        end if

        nvd   = nvib + ndis
        nrec  = 3 + ntarg + this % nchan + nvd
        ninfo = 1
        ndata = nrec - ninfo

        if (cform == 'F') then
            write (luchan, '(16I5)') keych, nchset, nrec, ninfo, ndata
            write (luchan, '(A80)') opts % diag_name
            write (luchan, '(16I5)') ntarg, nvib, ndis, this % nchan
            write (luchan, '(4I5,2D20.12)') opts % lambda, stot, gutot, ion, r, rmass

            k = 0
            do i = 1, prop % no_states
                if (findloc(this % idtarg, i, 1) /= 0) then
                    k = k + 1
                    isymm  = prop % energies_index(2, i)
                    mtarg  = prop % symmetries_index(1, isymm)
                    starg  = prop % symmetries_index(2, isymm)
                    gutarg = 0
                    etarg  = prop % energies(i)
                    write (luchan, '(4I5,2D20.12)') k, mtarg, starg, gutarg, etarg
                end if
            end do
            do i = 1, this % nchan
                write (luchan, '(4I5,2D20.12)') i, this % ichl(i), this % lchl(i), this % mchl(i), this % echl(i)
            end do
            do i = 1, nvd
                write (luchan, '(16I5)') i, ivtarg(i), iv(i)
            end do
        else
            write (luchan) keych, nchset, nrec, ninfo, ndata
            write (luchan) opts % diag_name(1:80)
            write (luchan) ntarg, nvib, ndis, this % nchan
            write (luchan) opts % lambda, stot, gutot, ion, r, rmass

            k = 0
            do i = 1, prop % no_states
                if (findloc(this % idtarg, i, 1) > 0) then
                    k = k + 1
                    isymm  = prop % energies_index(2, i)
                    mtarg  = prop % symmetries_index(1, isymm)
                    starg  = prop % symmetries_index(2, isymm)
                    gutarg = 0
                    etarg  = prop % energies(i)
                    write (luchan) k, mtarg, starg, gutarg, etarg
                end if
            end do
            do i = 1, this % nchan
                write (luchan) i, this % ichl(i), this % lchl(i), this % mchl(i), this % echl(i)
            end do
            do i = 1, nvd
                write (luchan) i, ivtarg(i), iv(i)
            end do
        end if

    end subroutine write_channel_info


    !> \brief   Write the R-matrix amplitudes to disk
    !> \authors J Benda
    !> \date    2019 - 2025
    !>
    !> Writes the information into the selected set of the R-matrix output file (*fort.21*).
    !> Adapted from WRITRM (UKRmol-out). The boundary amplitudes are written using derived-type I/O, so that if they are
    !> distributed, the data are gathered during the writing one column by another. This is to avoid the need to fetch
    !> the complete boundary amplitude matrix, which would negatively affect memory scaling.
    !>
    !> There is a bug (CMPLRLIBS-35389) in the Intel oneAPI compiler*, which prevents use of the derived-type I/O with large
    !> datasets. To overcome this bug, a temporary workaround is introduced. In case of Intel compiler the amplitudes are actually
    !> all collected by the master and written at once using a standard defined I/O. This construction can be removed
    !> once the compiler is fixed.
    !>
    !> *) https://community.intel.com/t5/Intel-Fortran-Compiler/Unformatted-DTIO-failure-with-gt-2GiB-records/m-p/1705503
    !>
    !> \param[in]  this       Interface object to update.
    !> \param[in]  nrmset     Output data set in the file unit.
    !> \param[in]  opts       SCATCI namelist information for this irreducible representation.
    !> \param[in]  lurmt      Unit number for boundary amplitudes output.
    !> \param[in]  rform      'F'/'U' for formatted or unformatted output.
    !> \param[in]  prop       Results from denprop.
    !> \param[in]  solution   Eigenvectors and related stuff calculated by the diagonalizer.
    !>
    subroutine write_boundary_data (this, nrmset, opts, lurmt, rform, prop, solution)

        use class_molecular_properties_data, only: molecular_properties_data
        use consts,                          only: amu
        use consts_mpi_ci,                   only: SYMTYPE_D2H
        use mpi_gbl,                         only: myrank, master
        use precisn_gbl,                     only: wp_bytes
        use Options_module,                  only: Options
        use SolutionHandler_module,          only: SolutionHandler

        class(OuterInterface),           intent(in) :: this
        integer,                         intent(in) :: nrmset, lurmt
        type(Options),                   intent(in) :: opts
        character(len=1),                intent(in) :: rform
        type(molecular_properties_data), intent(in) :: prop
        type(SolutionHandler),           intent(in) :: solution

        integer, parameter :: keyrm = 11
        character(len=10)  :: io_msg
        type(CIvect)       :: wamp

        integer  :: i, j, k, nchsq, nrec, ninfo, ndata, ntarg, nvib, ndis, nchan, stot, gutot, ion, ibut, iex, nocsf, npole, nnuc
        integer  :: ismax, nstat, io_stat, v_list(0)
        logical  :: intel_bug_workaround
        real(wp) :: r, rmass, rmatr

        nnuc  = prop % no_nuclei
        iex   = 0
        ibut  = 0
        nocsf = 0
        npole = 0
        ismax = 2
        rmatr = this % rmatr
        ntarg = this % ntarg
        nstat = solution % ci_vec % nstat
        nvib  = 0
        ndis  = 0
        nchan = this % nchan
        stot  = nint(2 * opts % spin + 1)
        gutot = SYMTYPE_D2H  ! to indicate polyatomic code

        if (opts % use_UKRMOL_integrals) then
           ion   = nint(sum(prop % charge(1:nnuc) - this % core_charges(1:nnuc))) - (opts % num_electrons - 1)
        else
           ion   = nint(sum(prop % charge(1:nnuc))) - (opts % num_electrons - 1)
        endif

        r     = 0
        nchsq = nchan * (nchan + 1) / 2
        nrec  = 4
        rmass = sum(1 / prop % mass(1:nnuc), prop % mass(1:nnuc) > 0)

        if (rmass > 0) then
            rmass = amu / rmass
        end if

        if (rform == 'F') then
            if (ismax > 0) nrec = nrec + (nchsq * ismax + 3) / 4
            nrec = nrec + (nstat + 3) / 4
            nrec = nrec + (nchan * nstat + 3) / 4
            if (npole > 0) nrec = nrec + (nocsf * npole + 3) / 4
            if (ibut > 0) nrec = nrec + (3 * nchan + 3) / 4
            if (abs(ibut) > 1) nrec = nrec + 1 + (nchan + 3) / 4 + (iex + 3) / 4 + (nchan * iex + 3) / 4
        else
            nrec = nrec + 2
            if (ismax > 0) nrec = nrec + 1
            if (npole > 0) nrec = nrec + 1
            if (ibut > 0) nrec = nrec + 1
            if (abs(ibut) > 1) nrec = nrec + 4
        end if

        ninfo = 1
        ndata = nrec - ninfo
        intel_bug_workaround = .false.

#ifdef __INTEL_COMPILER
        ! auxiliary master-only copy of the boundary amplitudes to work around Intel oneAPI bug CMPLRLIBS-35389;
        ! used for matrix size > 1 GiB = 2^30 B = 2^27 FP64 values
        if (nchan*nstat > 2**30/wp_bytes) then
            intel_bug_workaround = .true.
            wamp % nprow = 1
            wamp % npcol = 1
            wamp % blacs_context = -1
            wamp % CV_is_scalapack = this % wamp(this % nfdm + 1) % CV_is_scalapack
#if defined(usempi) && defined(scalapack)
            if (wamp % CV_is_scalapack) then
                call blacs_get(-1_blasint, 0_blasint, wamp % blacs_context)
                call blacs_gridinit(wamp % blacs_context, 'R', 1_blasint, 1_blasint)
            end if
#endif
            call wamp % init_CV(nchan, nstat)
        end if
#endif

#ifdef __GFORTRAN__
        ! Work around GCC PR fortran/122957: as of GCC 16, gfortran rejects the
        ! combination of -fdefault-integer-8 and user-defined derived-type I/O
        ! (the corresponding miscompile, GCC #108680, was the underlying defect).
        ! Force the same master-side auxiliary-copy path used by the Intel branch
        ! above, which avoids the DTIO write entirely.
        intel_bug_workaround = .true.
        wamp % nprow = 1
        wamp % npcol = 1
        wamp % blacs_context = -1
        wamp % CV_is_scalapack = this % wamp(this % nfdm + 1) % CV_is_scalapack
#if defined(usempi) && defined(scalapack)
        if (wamp % CV_is_scalapack) then
            call blacs_get(-1_blasint, 0_blasint, wamp % blacs_context)
            call blacs_gridinit(wamp % blacs_context, 'R', 1_blasint, 1_blasint)
        end if
#endif
        call wamp % init_CV(nchan, nstat)
#endif

        if (rform == 'F') then
            if (myrank == master) then
                write (lurmt, '(10I7)') keyrm, nrmset, nrec, ninfo, ndata
                write (lurmt, '(A80)') opts % diag_name
                write (lurmt, '(16I5)') ntarg, nvib, ndis, nchan
                write (lurmt, '(4I5,2D20.13)') opts % lambda, stot, gutot, ion, r, rmass
                write (lurmt, '(4I5,2D20.13)') ismax, nstat, npole, ibut, rmatr
                !if (abs(ibut) > 1) write(lurmt, '(i10,d20.13,i5)') nocsf, ezero, iex
                if (ismax > 0) write(lurmt, '(4D20.13)') ((this % a(i,k),i=1,nchsq),k=1,ismax)
                write (lurmt, '(4D20.13)') (solution % ci_vec % ei(i) + solution % core_energy,i=1,nstat)
                if (intel_bug_workaround) then
                    call wamp % redistribute(this%wamp(this%nfdm + 1), this%wamp(this%nfdm + 1)%blacs_context)
                    write (lurmt, '(4D20.13)') wamp % CV(1:nchan, 1:nstat)
#ifndef __GFORTRAN__
                else
                    write (lurmt, '(dt(4,20,13))') this % wamp(this % nfdm + 1)  ! derived-type I/O (formatted)
#endif
                end if
                !if (npole > 0) write (lurmt, '(4D20.13)') ((solution % ci_vec % cv(i,j),i=1,nocsf),j=1,npole)
                !if (ibut > 0) write (lurmt, '(4D20.13)') ((bcoef(i,j),i=1,3),j=1,nchan)
                !if (abs(ibut) > 1) write (lurmt, '(4D20.13)') (sfac(j),j=1,nchan)
                !if (abs(ibut) > 1) write (lurmt, '(4D20.13)') (ecex(j),j=1,iex)
                !if (abs(ibut) > 1) write (lurmt, '(4D20.13)') ((rcex(i,j),i=1,nchan),j=1,iex)
            else
                if (intel_bug_workaround) then
                    ! a dummy call from non-master tasks to collaborate with redistribution to master above
                    call wamp % redistribute(this%wamp(this%nfdm + 1), this%wamp(this%nfdm + 1)%blacs_context)
#ifndef __GFORTRAN__
                else
                    ! a dummy call from non-master tasks to collaborate with master's derived-type I/O above
                    call this % wamp(this % nfdm + 1) % formatted_write(-1, '', v_list, io_stat, io_msg)
#endif
                end if
            end if
        else
            if (myrank == master) then
                write(lurmt) keyrm, nrmset, nrec, ninfo, ndata
                write(lurmt) opts % diag_name(1:80)
                write(lurmt) ntarg, nvib, ndis, nchan
                write(lurmt) opts % lambda, stot, gutot, ion, r, rmass
                write(lurmt) ismax, nstat, npole, ibut, rmatr
                !if (abs(ibut) > 1) write(lurmt) nocsf, ezero, iex
                if (ismax > 0) write(lurmt) ((this % a(i,k),i=1,nchsq),k=1,ismax)
                write (lurmt) (solution % ci_vec % ei(i) + solution % core_energy,i=1,nstat)
                if (intel_bug_workaround) then
                    call wamp % redistribute(this%wamp(this%nfdm + 1), this%wamp(this%nfdm + 1)%blacs_context)
                    write (lurmt) wamp % CV(1:nchan, 1:nstat)
#ifndef __GFORTRAN__
                else
                    write (lurmt) this % wamp(this % nfdm + 1)  ! derived-type I/O (unformatted)
#endif
                end if
                !if (npole > 0) write (lurmt) ((solution % ci_vec % cv(i,j),i=1,nocsf),j=1,npole)
                !if (ibut > 0) write (lurmt) ((bcoef(i,j),i=1,3),j=1,nchan)
                !if (abs(ibut) > 1) write (lurmt) (sfac(j),j=1,nchan)
                !if (abs(ibut) > 1) write (lurmt) (ecex(j),j=1,iex)
                !if (abs(ibut) > 1) write (lurmt) ((rcex(i,j),i=1,nchan),j=1,iex)
            else
                if (intel_bug_workaround) then
                    ! a dummy call from non-master tasks to collaborate with redistribution to master above
                    call wamp % redistribute(this%wamp(this%nfdm + 1), this%wamp(this%nfdm + 1)%blacs_context)
#ifndef __GFORTRAN__
                else
                    ! a dummy call from non-master tasks to collaborate with master's derived-type I/O above
                    call this % wamp(this % nfdm + 1) % unformatted_write(-1, io_stat, io_msg)
#endif
                end if
            end if
        end if

    end subroutine write_boundary_data


    !> \brief   Evaluate long-range channel coupling coefficients
    !> \authors Z Masin
    !> \date    2014
    !>
    !> This is a completely rewritten version of the original SWINTERF routine SWASYMC. The coupling potentials for the outer region
    !> calculation are defined as:
    !>
    !> \f[
    !>     V_{ij}(r)=\sum_{\lambda=0}^{\infty}\frac{1}{r^{\lambda+1}}\times \nonumber \\
    !>     \times\underbrace{\sum_{m=-\lambda}^{\lambda}\langle{\cal Y}_{l_{i},m_{i}}\vert {\cal {Y}}_{\lambda ,m}\vert{\cal Y}_{l_{j},m_{j}}\rangle
    !>      \sqrt{\frac{4\pi}{2\lambda+1}}\underbrace{\sqrt{\frac{4\pi}{2\lambda+1}}\left( T_{ij}^{\lambda m} -
    !>      \langle \Phi_{i}\vert\Phi_{j}\rangle
    !>      \sum_{k=1}^{Nuclei}Z_{k}{\cal Y}_{\lambda,m}(\hat{\mathbf{R}}_{k})R_{k}^{\lambda}\right)}_{Q_{ij}^{\lambda m}}}_{a_{ij\lambda}} .
    !> \f]
    !>
    !> This routine calculates the coefficients \f$ a_{ij\lambda} \f$ using the information for all channels (i,j) and the target
    !> permanent and transition multipole moments \f$ Q_{ij}^{\lambda m} \f$. The main difference from SWASYMC is that here
    !> the coupling coefficients for the real spherical harmonics are calculated independently of the molecular orientation,
    !> symmetry and the lambda value. Tests were performed for pyrazine (\f$ D_{2h} \f$), uracil (\f$ C_s \f$) and water (\f$ C_{2v} \f$).
    !> The inclusion of the additional polarization potential is also possible, but note that only the spherical part is taken into
    !> account at the moment (see the comment for the variable `alpha2`).
    !>
    !> \note 06/02/2019 - Jakub Benda: Copied over from UKRmol-out.
    !>
    !> \param[inout] this      Interface object to update.
    !> \param[in] ismax        Maximum lambda value for the multipole potential contribution.
    !> \param[in] ntarg        Total number of target states.
    !> \param[in] target_properties   Denprop data.
    !> \param[in] alpha0       Spherical part of the ground state target polarizability.
    !> \param[in] alpha2       Non-spherical part of the ground state target polarizability. Note that this has not been implemented
    !>                         in this subroutine yet since it is not clear to me what is the convention defining this value in the old code.
    !> \param[in] use_pol      If .true. then the coefficients for the polarization potential will be constructed.
    !>
    subroutine get_channel_couplings (this, ismax, ntarg, target_properties, alpha0, alpha2, use_pol)

        use coupling_obj_gbl,                only: couplings_type
        use class_molecular_properties_data, only: molecular_properties_data
        use mpi_gbl,                         only: mpi_xermsg
        use precisn,                         only: wp

        class(OuterInterface),           intent(inout) :: this
        type(molecular_properties_data), intent(in)    :: target_properties

        integer,  intent(in)    :: ismax, ntarg
        real(wp), intent(in)    :: alpha0, alpha2
        logical,  intent(in)    :: use_pol

        real(wp), parameter :: fourpi  = 12.5663706143591729538505735331180115_wp
        real(wp), parameter :: roneh   =  0.70710678118654752440084436210484903_wp !sqrt(0.5)
        real(wp), parameter :: rothree =  0.57735026918962576450914878050195745_wp !1/sqrt(3)

        ! The Gaunt coefficients for the real spherical harmonics defined below are needed to express the xx, yy, zz angular
        ! behaviour in terms of the real spherical harmonics.

        ! x^2 = x*x ~ X_{11}*X_{11} = x2_X00*X_{00} + x2_X20*X_{20} + x2_X22*X_{22}
        real(wp), parameter :: x2_X00 =  0.282094791773878E+00_wp
        real(wp), parameter :: x2_X20 = -0.126156626101008E+00_wp
        real(wp), parameter :: x2_X22 =  0.218509686118416E+00_wp
        ! y^2 = y*y ~ X_{1-1}*X_{1-1} = y2_X00*X_{00} + y2_X20*X_{20} + y2_X22*X_{22}
        real(wp), parameter :: y2_X00 =  0.282094791773878E+00_wp
        real(wp), parameter :: y2_X20 = -0.126156626101008E+00_wp
        real(wp), parameter :: y2_X22 = -0.218509686118416E+00_wp
        ! z^2 = z*z ~ X_{10}*X_{10} = z2_X00*X_{00} + z2_X20*X_{20}
        real(wp), parameter :: z2_X00 =  0.282094791773878E+00_wp
        real(wp), parameter :: z2_X20 =  0.252313252202016E+00_wp

        type(couplings_type) :: couplings

        real(wp), allocatable :: prop(:,:,:)

        real(wp) :: cpl, sph_cpl, fac
        integer  :: l1, m1, q1, l2, m2, q2, no_cpl, lqt, isq, it1, it2
        integer  :: ch_1, ch_2, lambda, mlambda, lmin, lmax
        logical  :: use_alpha2

        if (use_pol .and. ntarg > 1) then
            write (*, '("WARNING: adding polarization potential while more than one target state is present in the outer region.")')
        end if

        ! get a convenient dense matrix set with the properties
        call target_properties % decompress(1, ismax, prop)

        ! find size parameters based on polarizability setup
        use_alpha2 = (use_pol .and. alpha2 /= 0.0_wp)
        lmax = merge(max(ismax, 3), ismax, use_pol)     ! maximal angular momentum transfer (polarizability corresponds to lambda=3)
        no_cpl = this % nchan * (this % nchan + 1) / 2  ! total number of unique combinations of the scattering channels

        ! allocate memory for the coupling coefficients
        if (allocated(this % a)) deallocate (this % a)
        allocate (this % a(no_cpl, lmax))               ! indices (:, 0) are never needed, though
        this % a(:,:) = 0.0_wp

        do ch_1 = 1, this % nchan

            l1 = this % lchl(ch_1)  ! the l, |m| values correspond to the l,|m| values of the real spherical harmonics
            m1 = this % mchl(ch_1)
            q1 = this % qchl(ch_1)  ! q is sign(m) or 0 if m=0.
            m1 = m1 * q1            ! determine m
            it1 = this % ichl(ch_1) ! sequence number of the target state corresponding to this channel

            do ch_2 = 1, ch_1

                l2 = this % lchl(ch_2)  ! the l, |m| values correspond to the l,|m| values of the real spherical harmonics
                m2 = this % mchl(ch_2)
                q2 = this % qchl(ch_2)  ! q is sign(m) or 0 if m=0.
                m2 = m2 * q2            ! determine m
                it2 = this % ichl(ch_2) ! sequence number of the target state corresponding to this channel

                ! use the selection rules for the real spherical harmonics to determine the range of lambda values
                ! which may give non-zero real Gaunt coefficients
                call couplings % bounds_rg(l1, l2, m1, m2, lmin, lmax)

                !don't include potentials with lambda > ismax
                if (lmax > ismax) lmax = ismax

                ! lambda = 0 would correspond to the monopole contribution, i.e. taking into account the total (perhaps nonzero)
                ! molecular charge. This is taken into account separately in RSOLVE.
                if (lmin == 0) lmin = 2

                ! linear index corresponding to the current combination of the channels (ch_1,ch_2).
                lqt = ch_1 * (ch_1 - 1) / 2 + ch_2

                ! Loop over the multipole moments which may be non-zero: we loop in steps of two since the selection rules
                ! for the r.s.h. imply that only if the sum of all L values is even the coupling may be non-zero.
                do lambda = lmin, lmax, 2

                    ! see the formula for V_{ij}: effectively this converts the inner region property from solid harmonic
                    ! normalization to spherical harmonic normalization as needed by the Leg. expansion.
                    fac = sqrt(fourpi / (2 * lambda + 1.0_wp))

                    do mlambda = -lambda, lambda

                        ! rgaunt: the coupling coefficient for the real spherical harmonics (l1,m1,l2,m2,lambda,mlambda);
                        ! The factor 2.0_wp converts the units of the inner region
                        ! properties from Hartree to Rydberg since that's the energy
                        ! unit used in the outer region.
                        cpl = 2.0_wp * fac * couplings % rgaunt(l1, lambda, l2, m1, mlambda, m2)

                        ! linear index corresponding to the current (lambda,mlambda) values. lambda*lambda is the number
                        ! of all (lambda1,mlambda1) combinations for lambda1 = lambda-1.
                        isq = lambda * lambda + lambda + mlambda

                        ! increment the value of the coefficient for the multipole coupling potential of order lambda between
                        ! the target states it1 and it2
                        this % a(lqt, lambda) = this % a(lqt, lambda) + cpl * prop(it1, it2, isq)

                    end do !mlambda

                end do !lambda

                ! add polarizabilities for the ground state channels, but only if use_pol == .true.
                if (use_pol .and. it1 == 1 .and. it2 == 1) then

                    ! the radial dependence of the polarization potential is r^{-4} which corresponds to lambda = 3
                    lambda = 3

                    ! obviously, the spherical polarizability (l=0) couples only the channels with identical angular behaviour
                    if (l1 == l2 .and. m1 == m2) then
                        sph_cpl = 1.0_wp
                        write(*, '("Spherical part of the polarization potential will be added to the &
                                   &channel (target state,l,m), coefficient value: (",3i5,") ",e25.15)') &
                                   it1, l1, m1, -sph_cpl * alpha0
                        this % a(lqt, lambda) = this % a(lqt, lambda) - sph_cpl * alpha0  ! add the spherical part of the polarizability
                    else
                        sph_cpl = 0.0_wp
                    end if

                    ! Add the non-spherical part of the polarizability; the polarization potential here corresponds to:
                    ! alpha_{2}*C_{i}*r_{-4}, where the angular behaviour of C_{i} = xy, xz, yz, xx, yy, zz
                    ! Note that we assume below that the polarizability tensor has the same values (alpha2) for all components.
                    ! Also note that as long as l1,m1 and l2,m2 belong to the same IR, only the totally symmetric components
                    ! of the pol. tensor (xx,yy,zz) may contribute.
                    if (use_alpha2) then

                        write(*, '("Non-spherical part of the polarization potential &
                                   &(target state,l1,m2,l2,m2): ",5i5)') it1, l1, m1, l2, m2

                        ! x^2 component
                        cpl = x2_X00 * couplings % rgaunt(l1, 0, l2, m1, 0, m2) + &
                              x2_X20 * couplings % rgaunt(l1, 2, l2, m1, 0, m2) + &
                              x2_X22 * couplings % rgaunt(l1, 2, l2, m1, 2, m2)
                        write(*, '("x^2 coefficient: ",e25.15)') -cpl * alpha2
                        this % a(lqt, lambda) = this % a(lqt, lambda) - cpl * alpha2

                        ! y^2 component
                        cpl = y2_X00 * couplings % rgaunt(l1, 0, l2, m1, 0, m2) + &
                              y2_X20 * couplings % rgaunt(l1, 2, l2, m1, 0, m2) + &
                              y2_X22 * couplings % rgaunt(l1, 2, l2, m1, 2, m2)
                        write(*, '("y^2 coefficient: ",e25.15)') -cpl * alpha2
                        this % a(lqt, lambda) = this % a(lqt, lambda) - cpl * alpha2

                        ! z^2 component
                        cpl = z2_X00 * couplings % rgaunt(l1, 0, l2, m1, 0, m2) + &
                              z2_X20 * couplings % rgaunt(l1, 2, l2, m1, 0, m2)
                        write(*, '("z^2 coefficient: ",e25.15)') -cpl * alpha2
                        this % a(lqt, lambda) = this % a(lqt, lambda) - cpl * alpha2

                        ! xz component
                        cpl = couplings % rgaunt(l1, 2, l2, m1, 1, m2) * rothree
                        write(*, '("xz coefficient: ",e25.15)') -cpl * alpha2
                        this % a(lqt, lambda) = this % a(lqt, lambda) - cpl * alpha2

                        ! yz component
                        cpl = couplings % rgaunt(l1, 2, l2, m1, -1, m2) * rothree
                        write(*, '("yz coefficient: ",e25.15)') -cpl * alpha2
                        this % a(lqt, lambda) = this % a(lqt, lambda) - cpl * alpha2

                        ! xy component
                        cpl = couplings % rgaunt(l1, 2, l2, m1, -2, m2) * rothree
                        write(*, '("xy coefficient: ",e25.15)') -cpl * alpha2
                        this % a(lqt, lambda) = this % a(lqt, lambda) - cpl * alpha2

                  end if  ! use_alpha2

               end if  ! use_pol

            end do  ! ch_2

         end do  ! ch_1

    end subroutine get_channel_couplings


    !> \brief   Compose the RMT molecular input data file
    !> \authors J Benda
    !> \date    2019 - 2024
    !>
    !> Based on `generate_data_for_rmt` from the original program `rmt_interface` (UKRmol-out).
    !>
    !> This subroutine writes the input data file used by molecular version of RMT. The file is always called 'molecular_data',
    !> it is a binary stream file and it contains the following data:
    !>  - (For m = -1, 0, 1) s, s1, s2, iidip(1:s), ifdip(1:s), dipsto(1:s1,1:s2, 1:s): Number of distinct symmetry pairs, maximal
    !>    number of states in bra symmetry maximal number of states in ket symmetry, list of bra symmetries, list of ket symmetries,
    !>    (N + 1)-electron state transition dipole blocks (padded by zeros to s1 x s2).
    !>  - ntarg: Number of targets (N-electron ions).
    !>  - (For m = -1, 0, 1) crlv: Dipole transition matrix for N-electron states
    !>  - n_rg, rg, lm_rg: Number of real Gaunt coeffcients, array of them, their angular momentum quantum numbers.
    !>  - nelc, nz, lrang2, lamax, ntarg, inast, nchmx, nstmx, lmaxp1: Number of electrons, nuclear charge, maximal angular momentum,
    !>    maximal angular momentum transfer, number of target states, number of irreducible representations, maximal number of channels
    !>    in all irreducible representations, maximal number of eigenstates in all irreducible representations, maximal angular
    !>    momentum.
    !>  - rmatr, bbloch: R-matrix radius, Bloch b-coeffcient.
    !>  - etarg, ltarg, starg: Ionic state eigen-energies, their irreducible representations (one-based), and spins.
    !>  - nfdm, delta_r: Number of finite difference points inside the R-matrix sphere, their spacing.
    !>  - r_points: List of the finite difference points inside the R-matrix sphere.
    !>  - (For all target symmetries) lrgl, nspn, npty, nchan, mnp1, nconat, l2p, m2p, eig, wamp, cf, ichl, s1, s2, wamp2:
    !>    Irreducible representation (one-based), spin multiplicity, parity (not used, always zero), number of continuum channels,
    !>    number of (N + 1)-electron states, number of continuum channels per target (ionic) state, angular momentum per channel,
    !>    angular momentum projection per channel, eigenenergies of the (N + 1)-electron states (+ core energy), contributions of the
    !>    (N + 1)-electron eigenstates to every channel (boundary amplitudes), long-range multipole coupling coefficients,
    !>    target state index per continuum channel, row size of wamp2, col size of wamp2, wamp evaluated at inner-region
    !>    finite difference points.
    !>
    !> \todo At the moment, the subroutine assumes that all inner states have the same spin.
    !>
    !> \warning The subroutine makes use of MPI I/O for writing the distributed arrays. In the present, straightforward implementation
    !>          this limits the local portion of the distributed matrices to 2**31 elements (i.e. 16 GiB for 8-byte real arrays).
    !>          If this became a problem, one would need to rewrite the subroutine to write everything in one go as a custom MPI
    !>          datatype, whose building routines accept long integers. However, 16 GiB per core per distributed matrix seems quite
    !>          acceptable for now.
    !>
    !> \param[in] input              Input SCATCI namelist data for all symmetries.
    !> \param[in] iitf               Which outer_interface namelist to use.
    !> \param[in] inner_properties   Results from cdenprop.
    !> \param[in] solutions          Results from diagonalization.
    !> \param[in] intf               Extracted channel and boundary data.
    !>
    subroutine write_rmt_data (input, iitf, inner_properties, solutions, intf)

        use algorithms,                      only: findloc, insertion_sort
        use class_molecular_properties_data, only: molecular_properties_data
        use const_gbl,                       only: stdout
        use consts_mpi_ci,                   only: PASS_TO_CDENPROP, NO_DIAGONALIZATION
        use mpi_gbl,                         only: myrank, master, mpiint, mpi_xermsg, mpi_mod_file_open_write, &
                                                   mpi_mod_file_close, mpi_mod_file_set_size, mpi_mod_file_write, mpi_mod_barrier
        use Options_module,                  only: OptionsSet
        use Timing_module,                   only: master_timer
        use SolutionHandler_module,          only: SolutionHandler

        type(OptionsSet),                   intent(in) :: input
        integer,                            intent(in) :: iitf
        type(molecular_properties_data),    intent(in) :: inner_properties
        type(SolutionHandler), allocatable, intent(in) :: solutions(:)
        type(OuterInterface),  allocatable, intent(in) :: intf(:)

        real(rmt_real),   allocatable :: cf(:,:), crlv(:,:), etarg(:), rg(:), eig(:), r_points(:)
        integer(rmt_int), allocatable :: iidip(:), ifdip(:), lm_rg(:,:), ltarg(:), starg(:), nconat(:), l2p(:), m2p(:), ichl(:), &
                                         propblocks(:)
        integer,          allocatable :: ions(:)

        type(molecular_properties_data) :: target_properties

        type(CIvect)     :: wmat
        real(rmt_real)   :: bbloch = 0, prop, rmatr, delta_r
        integer          :: i, j, k, l, m, irr_map(8,8), u, v, nsymm, ind, ip, irr1, irr2, it1, it2, nnuc, &
                            ierr, nblock, sym1, sym2, ref
        integer(mpiint)  :: fh
        integer(rmt_int) :: s1, s2, ntarg, n_rg, nelc, nz, lamax, inast, nchmx, nstmx, lmaxp1, nchan, lrgl, nspn, npty, &
                            mnp1, mxstat, lrang2, lmax, nfdm

        call target_properties % read_properties(input % interface_opts(iitf) % lutdip, 1)

        ! number of N + 1 irrs for which we have inner solution and dipoles
        inast   = count(iand(input % opts(:) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                        input % opts(:) % diagonalization_flag /= NO_DIAGONALIZATION)

        ref     = maxloc(intf(:) % ntarg, 1)                ! irr with the first non-zero number of outer region states
        nnuc    = target_properties % no_nuclei             ! number of nuclei in the molecule
        nfdm    = input % interface_opts(iitf) % nfdm       ! number of boundary amplitudes evaluation radii *inside* R-matrix sphere
        delta_r = input % interface_opts(iitf) % delta_r    ! spacing between the evaluation distances
        rmatr   = input % interface_opts(iitf) % rmatr      ! radius of the inner region
        nsymm   = size(input % opts(:))                     ! number of N + 1 irrs for which we have inputs (= size(intf(:)))
        ntarg   = intf(ref) % ntarg                         ! number of target ("ionic") states
        nelc    = input % opts(1) % num_electrons - 1       ! number of electrons of the target (ion)
        mxstat  = maxval(inner_properties % dense_blocks(1:inner_properties % no_blocks) % mat_dimen)  ! largest dipole block size

        ! avoid the routine altogether when there are no data for any reason
        if (nsymm == 0 .or. inast == 0 .or. ntarg == 0 .or. mxstat <= 0) then
            write (stdout, *) 'Nothing to do in RMT interface'
            return
        end if

        nz = sum(target_properties % charge(1:nnuc) - intf(ref) % core_charges(1:nnuc))   ! number of protons in the molecule

        lrang2  = 0         ! angular momentum limit
        lamax   = 0         ! angular momentum limit
        nchmx   = 0         ! maximal number of channels across irreducible representations
        nstmx   = 0         ! maximal number of inner region eigenstates across representations

        call master_timer % start_timer("RMT interface")

        ! get energy-sorted list of ionic states for use in outer region (only a subset of all if IDTARG is incomplete)
        ions = intf(ref) % idtarg
        call insertion_sort(ions)

        ! open the output file and reset its contents
        call mpi_mod_file_open_write('molecular_data', fh, ierr)
        call mpi_mod_file_set_size(fh, 0)

        allocate (propblocks(inner_properties % no_blocks))

        write (stdout, '(/,"Writing RMT data file")')
        write (stdout, '(  "=====================")')
        write (stdout, '(/,"  Inner region dipole block max size: ",I0)') mxstat
        write (stdout, '(/,"  Ionic states used in outer region: ")')
        write (stdout, '(/,"    ",20I5)') ions

        ! store all (N + 1)-electron dipole blocks
        do m = -1, 1

            ! find out which IRR pairs yield non-zero matrix elements of Y[1,m]
            irr_map = 0
            nblock = 0
            do ip = 1, inner_properties % no_blocks
                if (inner_properties % dense_index(3, ip) == 1 .and. &
                    inner_properties % dense_index(4, ip) == m) then
                    sym1 = inner_properties % dense_index(1, ip)
                    sym2 = inner_properties % dense_index(2, ip)
                    irr1 = inner_properties % symmetries_index(1, sym1) + 1
                    irr2 = inner_properties % symmetries_index(1, sym2) + 1
                    nblock = nblock + 1
                    propblocks(nblock) = ip
                    irr_map(irr1, irr2) = nblock
                end if
            end do

            ! get some free memory
            if (allocated(iidip)) deallocate (iidip)
            if (allocated(ifdip)) deallocate (ifdip)
            allocate (iidip(nblock))
            allocate (ifdip(nblock))
            iidip(:) = -1
            ifdip(:) = -1

            ! compose sequence of bra and ket symmetry indices
            do ip = 1, nblock
                sym1 = inner_properties % dense_index(1, propblocks(ip))
                sym2 = inner_properties % dense_index(2, propblocks(ip))
                irr1 = inner_properties % symmetries_index(1, sym1) + 1
                irr2 = inner_properties % symmetries_index(1, sym2) + 1
                iidip(ip) = irr1
                ifdip(ip) = irr2
            end do

            ! write dipole metadata
            call mpi_mod_file_write(fh, int(nblock, rmt_int))
            call mpi_mod_file_write(fh, mxstat)
            call mpi_mod_file_write(fh, mxstat)
            call mpi_mod_file_write(fh, iidip, nblock)
            call mpi_mod_file_write(fh, ifdip, nblock)

            ! collective write of the ScaLAPACK matrices to the stream file
            do ip = 1, nblock
                call inner_properties % dense_blocks(propblocks(ip)) % stream_write(fh)
            end do
        end do

        ! set up angular momentum limits
        do i = 1, nsymm
            if (iand(input % opts(i) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then
                if (intf(i) % nchan > 0) then
                    lrang2 = max(int(lrang2), maxval(intf(i) % lchl(1:intf(i) % nchan)))
                    lamax  = max(int(lamax), intf(i) % ismax)
                    nchmx  = max(int(nchmx), intf(i) % nchan)
                end if
                nstmx = max(int(nstmx), solutions(i) % vec_dimen)
            end if
        end do
        lrang2 = lrang2 + 1
        lmaxp1 = lrang2

        ! number of target (ion) states
        call mpi_mod_file_write(fh, ntarg)

        ! allocate memory for auxiliary arrays (dummy for non-master processes)
        if (myrank == master) then
            allocate (etarg(ntarg), ltarg(ntarg), starg(ntarg), nconat(ntarg), l2p(nchmx), m2p(nchmx), eig(nstmx), &
                      cf(nchmx, nchmx), ichl(nchmx), crlv(ntarg, ntarg), r_points(nfdm + 1), stat = ierr)
        else
            allocate (etarg(0), ltarg(0), starg(0), nconat(0), l2p(0), m2p(0), eig(0), rg(0), lm_rg(0, 0), &
                      cf(0, 0), ichl(0), crlv(0, 0), r_points(0), stat = ierr)
        end if
        if (ierr /= 0) then
            call mpi_xermsg('Postprocessing_module', 'write_rmt_data', 'Memory allocation failure.', 1, 1)
        end if

        ! set up a distributed matrix for reshaping the boundary amplitudes
        wmat % nprow = 1
        wmat % npcol = 1
        wmat % blacs_context = intf(1) % wamp(nfdm + 1) % blacs_context
        wmat % CV_is_scalapack = intf(1) % wamp(nfdm + 1) % CV_is_scalapack
        call wmat % init_CV(int(nchmx), int(nstmx))

        ! target data
        if (myrank == master) then
            do i = 1, ntarg
                etarg(i) = target_properties % energies(ions(i))
                ltarg(i) = target_properties % symmetries_index(1, target_properties % energies_index(2, ions(i))) + 1
                starg(i) = target_properties % symmetries_index(2, target_properties % energies_index(2, ions(i)))
            end do
        end if

        ! finite difference points inside (and at) the R-matrix sphere
        if (myrank == master) then
            do i = 1, nfdm + 1
                r_points(i) = rmatr - (nfdm + 1 - i) * delta_r
            end do
        end if

        ! store all N-electron dipole blocks
        do m = -1, 1
            if (myrank == master) then
                crlv(:,:) = 0
                do ip = 1, target_properties % non_zero_properties
                    if (target_properties % properties_index(3, ip) == 1 .and. &
                        target_properties % properties_index(4, ip) == m) then
                        ! get state information for this property
                        it1  = target_properties % properties_index(1, ip)
                        it2  = target_properties % properties_index(2, ip)
                        prop = target_properties % properties(ip)
                        ! check that these two states are requested in the outer region
                        it1 = findloc(ions, it1, 1)
                        it2 = findloc(ions, it2, 1)
                        if (it1 > 0 .and. it2 > 0) then
                            crlv(it1, it2) = prop
                            crlv(it2, it1) = prop
                        end if
                    end if
                end do
            end if

            call mpi_mod_file_write(fh, crlv, int(ntarg), int(ntarg))
        end do

        ! generate the coupling coefficients needed in RMT to construct the laser-target and laser-electron asymptotic potentials
        lmax = 0
        do i = 1, nsymm
            if (iand(input % opts(i) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then
                if (intf(i) % nchan > 0) then
                    lmax = max(int(lmax), maxval(intf(i) % lchl(1:intf(i) % nchan)))
                end if
            end if
        end do
        if (myrank == master) then
            call generate_couplings(lmax, n_rg, rg, lm_rg)
        end if

        call mpi_mod_file_write(fh, n_rg)
        call mpi_mod_file_write(fh, rg, int(n_rg))
        call mpi_mod_file_write(fh, lm_rg, 6, int(n_rg))
        call mpi_mod_file_write(fh, nelc)
        call mpi_mod_file_write(fh, nz)
        call mpi_mod_file_write(fh, lrang2)
        call mpi_mod_file_write(fh, lamax)
        call mpi_mod_file_write(fh, ntarg)
        call mpi_mod_file_write(fh, inast)
        call mpi_mod_file_write(fh, nchmx)
        call mpi_mod_file_write(fh, nstmx)
        call mpi_mod_file_write(fh, lmaxp1)
        call mpi_mod_file_write(fh, rmatr)
        call mpi_mod_file_write(fh, bbloch)
        call mpi_mod_file_write(fh, etarg, int(ntarg))
        call mpi_mod_file_write(fh, ltarg, int(ntarg))
        call mpi_mod_file_write(fh, starg, int(ntarg))
        call mpi_mod_file_write(fh, nfdm)
        call mpi_mod_file_write(fh, delta_r)
        call mpi_mod_file_write(fh, r_points, int(nfdm + 1))

        ! write channel information (per symmetry)
        do i = 1, nsymm
            if (iand(input % opts(i) % vector_storage_method, PASS_TO_CDENPROP) /= 0 .and. &
                input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then

                if (myrank == master) then
                    lrgl  = 1 + input % opts(i) % lambda
                    nspn  = int(2 * input % opts(i) % spin + 1, rmt_int)
                    npty  = 0
                    nchan = intf(i) % nchan
                    mnp1  = solutions(i) % vec_dimen

                    do j = 1, ntarg
                        nconat(j) = count(intf(i) % ichl(1:nchan) == j)
                    end do

                    l2p(:)  = 0
                    m2p(:)  = 0
                    ichl(:) = 0

                    do j = 1, nchan
                        l2p(j) = intf(i) % lchl(j)
                        m2p(j) = intf(i) % mchl(j)
                        ichl(j) = intf(i) % ichl(j)
                    end do
                end if

                call mpi_mod_file_write(fh, lrgl)
                call mpi_mod_file_write(fh, nspn)
                call mpi_mod_file_write(fh, npty)
                call mpi_mod_file_write(fh, nchan)
                call mpi_mod_file_write(fh, mnp1)
                call mpi_mod_file_write(fh, nconat, int(ntarg))
                call mpi_mod_file_write(fh, l2p, int(nchmx))
                call mpi_mod_file_write(fh, m2p, int(nchmx))

                if (myrank == master) then
                    eig(:) = 0
                    eig(1:mnp1) = solutions(i) % ci_vec % ei(1:mnp1) + solutions(i) % core_energy
                end if

                call mpi_mod_file_write(fh, eig, int(nstmx))

                ! resize the boundary amplitudes matrix to nchmx by nstmx
                wmat % CV = 0
                call wmat % redistribute(intf(i) % wamp(nfdm + 1), intf(i) % wamp(nfdm + 1) % blacs_context)
                wmat % CV = sqrt(2.0_wp) * wmat % CV    ! switch to RMT convention
                call wmat % stream_write(fh)

                do l = 1, lamax
                    if (myrank == master) then
                        cf(:,:) = 0
                        if (l <= intf(i) % ismax) then
                            do k = 1, nchan
                                do j = 1, nchan
                                    u = max(j, k)
                                    v = min(j, k)
                                    ind = u * (u - 1) / 2 + v
                                    cf(j, k) = intf(i) % a(ind, l)
                                end do
                            end do
                        end if
                    end if
                    call mpi_mod_file_write(fh, cf, int(nchmx), int(nchmx))
                end do

                call mpi_mod_file_write(fh, ichl, int(nchmx))

                s1 = intf(i) % wamp(1) % mat_dimen_r
                s2 = intf(i) % wamp(1) % mat_dimen_c

                call mpi_mod_file_write(fh, s1)
                call mpi_mod_file_write(fh, s2)

                ! collective write of the ScaLAPACK matrix to the stream file
                do j = 1, nfdm
                    call intf(i) % wamp(j) % stream_write(fh)
                end do

            end if
        end do

        call mpi_mod_file_close(fh, ierr)

        call master_timer % stop_timer("RMT interface")

    end subroutine write_rmt_data


    !> \brief   Evaluate angular couplings for outer region of RMT
    !> \authors Z Masin
    !> \date    2017
    !>
    !> Precompute the necessary Gaunt coefficients.
    !>
    !> \param[in]  maxl   Limit on angular quantum number.
    !> \param[out] n_rg   Number of calculated Gaunt coefficients.
    !> \param[out] rg     Values of the non-zero Gaunt coefficients.
    !> \param[out] lm_rg  Integer sextets l1,m1,l2,m2,l3,m3 for each of the non-zero Gaunt coefficients.
    !>
    subroutine generate_couplings (maxl, n_rg, rg, lm_rg)

        use coupling_obj_gbl, only: couplings_type

        integer(rmt_int),              intent(in)    :: maxl
        integer(rmt_int),              intent(out)   :: n_rg
        real(rmt_real),   allocatable, intent(inout) :: rg(:)
        integer(rmt_int), allocatable, intent(inout) :: lm_rg(:,:)

        type(couplings_type) :: cpl

        integer  :: l1, l2, l3, m1, m2, m3, err
        real(wp) :: tmp

        call cpl % prec_cgaunt(int(maxl))

        l3 = 1

        ! How many non-zero dipole couplings there are
        n_rg = 0
        do l1 = 0, maxl
            do m1 = -l1, l1
                do l2 = 0, maxl
                    do m2 = -l2, l2
                        do m3 = -l3, l3
                            tmp = cpl % rgaunt(l1, l2, l3, m1, m2, m3)
                            if (tmp /= 0) n_rg = n_rg + 1
                        end do !m3
                    end do !m2
                end do !l2
            end do !m1
        end do !l1

        if (allocated(rg)) deallocate (rg)
        if (allocated(lm_rg)) deallocate (lm_rg)
        n_rg = max(n_rg, 1_rmt_int)  ! make sure that even if all couplings are zero this routine allocates rg, lm_rg
        allocate (rg(n_rg), lm_rg(6, n_rg), stat = err)
        rg = 0
        lm_rg = -2

        ! Save them
        n_rg = 0
        do l1 = 0, maxl
            do m1 = -l1, l1
                do l2 = 0, maxl
                    do m2 = -l2, l2
                        do m3 = -l3, l3
                            tmp = cpl % rgaunt(l1, l2, l3, m1, m2, m3)
                            if (tmp /= 0) then
                                n_rg = n_rg + 1
                                rg(n_rg) = tmp
                                lm_rg(1:6, n_rg) = (/ l1, m1, l2, m2, l3, m3 /)
                            end if
                        end do !m3
                    end do !m2
                end do !l2
            end do !m1
        end do !l1

   end subroutine generate_couplings

end module Postprocessing_module
