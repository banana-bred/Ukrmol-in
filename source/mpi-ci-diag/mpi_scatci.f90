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

!> \brief   MPI-SCATCI main program
!> \authors A Al-Refaie, J Benda
!> \date    2017 - 2022
!>
!> MPI-SCATCI diagonalizes Hamiltonians for one or more irreducible representations. The input for each
!> is read from the single input file as the pair of standard SCATCI namelists &input and &cinorn.
!> If there are fewer MPI processes than Hamiltonians to diagonalize, each process is given several
!> Hamiltonians to diagonalize sequentially. If there are more processes than Hamiltonians, some or all
!> Hamiltonians will be diagonalized in a distributed regime. All processes that participate on diagonalization
!> of the same matrix are part of a "diagonalization group" which has separate BLACS context and their own
!> MPI communicator.
!>
!> The eigenvectors are kept in the memory and, if desired, are then used to calculate transition dipole
!> moments and further derived data.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!> \note 25/01/2019 - Jakub Benda: Extension to support multiple subsequent diagonalizations.
!> \note 30/01/2019 - Jakub Benda: Evaluation of properties using libcdenprop, including cross-symmetry ones.
!> \note 22/05/2022 - Jakub Benda: Reimplemented reading of eigensolutions from disk.
!>
program mpi_scatci

    use consts_mpi_ci,                   only: NO_DIAGONALIZATION, MAIN_HAMILTONIAN, TARGET_HAMILTONIAN, NO_CI_TARGET, &
                                               READ_FROM_FILE
    use const_gbl,                       only: stdout
    use global_utils,                    only: print_ukrmol_header
    use precisn,                         only: wp, longint
    use mpi_gbl,                         only: mpi_mod_start, mpi_mod_finalize, master, myrank, mpi_mod_print_info
    use Options_module,                  only: OptionsSet, ham_pars_io_wrapper
    use Orbital_module,                  only: OrbitalTable
    use CSF_module,                      only: CSFObject, CSFOrbital, CSFManager
    use CI_Hamiltonian_module,           only: Target_CI_Hamiltonian
    use BaseMatrix_module,               only: BaseMatrix
    use BaseIntegral_module,             only: BaseIntegral
    use Diagonalizer_module,             only: BaseDiagonalizer
    use MemoryManager_module,            only: master_memory
    use Parallelization_module,          only: process_grid
    use Postprocessing_module,           only: postprocess
    use Uncontracted_Hamiltonian_module, only: Uncontracted_Hamiltonian
    use Contracted_Hamiltonian_module,   only: Contracted_Hamiltonian
    use Timing_module,                   only: master_timer
    use Target_RMat_CI_module,           only: Target_RMat_CI, read_ci_mat
    use Dispatcher_module,               only: DispatchIntegral, DispatchMatrixAndDiagonalizer, initialize_libraries
    use WriterMatrix_module,             only: WriterMatrix
    use SolutionHandler_module,          only: SolutionHandler

    implicit none

    type(OptionsSet)                        ::      SCATCI_input
    type(OrbitalTable)                      ::      Orbitals
    type(Uncontracted_Hamiltonian)          ::      enrgms
    type(Contracted_Hamiltonian)            ::      enrgmx
    type(CSFManager)                        ::      Configuration_manager
    class(BaseIntegral),          pointer   ::      Integrals           => null()
    class(Target_CI_Hamiltonian), pointer   ::      tgt_ci_hamiltonian  => null()
    class(BaseDiagonalizer),      pointer   ::      diagonalizer        => null()
    class(BaseMatrix),            pointer   ::      matrix_elements     => null()
    type(SolutionHandler),    allocatable   ::      solutions(:)
    type(CSFObject),          allocatable   ::      CSFS(:)
    type(Target_RMat_CI),     allocatable   ::      ci_rmat(:)
    real(wp),                 allocatable   ::      test_eig(:), test_vecs(:,:)
    integer                                 ::      num_mat_elms, i, j
    real(wp)                                ::      usr_ef
    integer                                 ::      usr_ef_type

    logical, parameter :: sequential_diagonalizations = .false.
    logical, parameter :: master_writes_to_stdout = .true.
    logical, parameter :: allow_shared_memory = .true.

    call mpi_mod_start(master_writes_to_stdout, allow_shared_memory)
    call master_timer % initialize()
    if (myrank == master) call print_ukrmol_header(stdout)
    call mpi_mod_print_info(stdout)

    ! Initialize all libraries used by MPI-SCATCI that need it (e.g. SLEPc)
    call initialize_libraries

    ! Read options (for all symmetries)
    call SCATCI_input % read
    allocate (solutions(size(SCATCI_input % opts)))

    ! Set up the MPI/BLACS groups
    call process_grid % setup(size(SCATCI_input % opts), sequential_diagonalizations)
    call process_grid % summarize
    if (.not. process_grid % sequential) call SCATCI_input % setup_write_order

    call master_timer % start_timer("Wall Time")

    ! Process all symmetries
    symmetry_loop: do i = 1, size(SCATCI_input % opts)

        if (.not. process_grid % is_my_group_work(i)) cycle

        call Orbitals % initialize(SCATCI_input % opts(i) % tot_num_spin_orbitals,    &
                                   SCATCI_input % opts(i) % tot_num_orbitals,         &
                                   SCATCI_input % opts(i) % sym_group_flag,           &
                                   SCATCI_input % opts(i) % num_syms,                 &
                                   SCATCI_input % opts(i) % positron_flag)
        call Orbitals % construct(SCATCI_input % opts(i) % num_orbitals_sym,                    &
                                  SCATCI_input % opts(i) % num_target_orbitals_sym_dinf_congen, &
                                  SCATCI_input % opts(i) % num_orbitals_sym,                    &
                                  SCATCI_input % opts(i) % num_electronic_orbitals_sym,         &
                                  SCATCI_input % opts(i) % num_target_orbitals_sym_congen)

        ! Compute the electron number for each reference orbital
        call Orbitals % compute_electron_index(SCATCI_input % opts(i) % num_electrons, &
                                               SCATCI_input % opts(i) % reference_dtrs)

        ! Construct our csf manager
        call Configuration_manager % initialize(SCATCI_input % opts(i), Orbitals)

        ! Construct our CSFS
        call master_timer % start_timer("Construct CSFs")
        call Configuration_manager % create_csfs(CSFS,                      &
                SCATCI_input % opts(i) % orbital_sequence_number,           &
                SCATCI_input % opts(i) % num_ci_target_sym,                 &
                SCATCI_input % opts(i) % lambda_continuum_orbitals_target)
        call master_timer % stop_timer("Construct CSFs")
        call master_memory % print_memory_report

        ! Read molecular integrals (assuming all symmetries use the same!)
        if (.not. associated(Integrals)) then
            call DispatchIntegral(SCATCI_input % opts(i) % sym_group_flag, &
                                  SCATCI_input % opts(i) % use_UKRMOL_integrals, &
                                  Integrals)
            call Integrals % initialize(SCATCI_input % opts(i), &
                                        Orbitals % total_num_orbitals, &
                                        Orbitals % orbital_map)
            call Integrals % load_integrals(SCATCI_input % opts(i) % integral_unit)
            call master_memory % print_memory_report
        end if

        ! if only reading already calculated eigensolutions (for processing in CDENPROP or outer interface)
        if (SCATCI_input % opts(i) % diagonalization_flag == READ_FROM_FILE) then
            call solutions(i) % construct(SCATCI_input % opts(i))
            call solutions(i) % export_header(SCATCI_input % opts(i), Integrals)
            call solutions(i) % read_from_file(SCATCI_input % opts(i))
            deallocate(CSFS)
            call Configuration_manager % finalize
            cycle
        end if

        write (stdout, *) '  Load CI target and construct matrix elements'

        ! Before we do anything else let us get our ECP if we need them
       !call DispatchECP(SCATCI_input % opts(i) % ecp_type,             &
       !                 SCATCI_input % opts(i) % ecp_filename,         &
       !                 SCATCI_input % opts(i) % all_ecp_defined,      &
       !                 SCATCI_input % opts(i) % num_target_sym,       &
       !                 SCATCI_input % opts(i) % num_target_state_sym, &
       !                 SCATCI_input % opts(i) % target_spatial,       &
       !                 SCATCI_input % opts(i) % target_multiplicity, ecp)

        if (SCATCI_input % opts(i) % do_ci_contraction()) then
            usr_ef = SCATCI_input % opts(i) % enh_factor
            usr_ef_type = SCATCI_input % opts(i) % enh_factor_type
            SCATCI_input % opts(i) % enh_factor = 1.0
            SCATCI_input % opts(i) % enh_factor_type = 0
            ! allocate the CIRmats
            allocate(ci_rmat(SCATCI_input % opts(i) % num_target_sym))
            do j = 1, SCATCI_input % opts(i) % num_target_sym
                call ci_rmat(j) % initialize(j,                                                &
                                             SCATCI_input % opts(i) % num_target_state_sym(j), &
                                             SCATCI_input % opts(i) % num_ci_target_sym(j),    &
                                             SCATCI_input % opts(i) % ci_phase(j),             &
                                             Integrals % get_core_energy())
            end do

            if (SCATCI_input % opts(i) % ci_target_switch > 0) then
                call read_ci_mat(SCATCI_input % opts(i), ci_rmat)
            else
                do j = 1, SCATCI_input % opts(i) % num_target_sym
                    allocate(Target_CI_Hamiltonian::tgt_ci_hamiltonian)
                    call DispatchMatrixAndDiagonalizer(SCATCI_input % opts(i) % diagonalizer_choice,    &
                                                       SCATCI_input % opts(i) % force_serial,           &
                                                       SCATCI_input % opts(i) % num_ci_target_sym(j),   &
                                                       ci_rmat(j) % nstat,                              &
                                                       matrix_elements,                                 &
                                                       diagonalizer,                                    &
                                                       Integrals,                                       &
                                                       SCATCI_input % opts(i) % hamiltonian_unit)

                    call matrix_elements % construct(TARGET_HAMILTONIAN + j)
                    call matrix_elements % set_options(SCATCI_input % opts(i))

                    call tgt_ci_hamiltonian % construct(SCATCI_input % opts(i), CSFS, Orbitals, Integrals)
                    call tgt_ci_hamiltonian % initialize(j)

                    call master_timer % start_timer("Build CI Hamiltonian")
                    call tgt_ci_hamiltonian % build_hamiltonian(matrix_elements)
                    call master_timer % stop_timer("Build CI Hamiltonian")
                    call master_timer % report_timers

                    call master_timer % start_timer("Diagonalization")
                    call diagonalizer % diagonalize(matrix_elements, ci_rmat(j) % nstat, ci_rmat(j), .true., &
                                                    SCATCI_input % opts(i), Integrals)
                    call master_timer % stop_timer("Diagonalization")

                    call ci_rmat(j) % print
                    call master_timer % report_timers

                    call matrix_elements % destroy
                    if (process_grid % grank == master) close (SCATCI_input % opts(i) % hamiltonian_unit)
                    deallocate(matrix_elements, diagonalizer)
                    deallocate(tgt_ci_hamiltonian)
                end do
            end if
            SCATCI_input % opts(i) % enh_factor = usr_ef
            SCATCI_input % opts(i) % enh_factor_type = usr_ef_type
        end if

        if (SCATCI_input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then
            call DispatchMatrixAndDiagonalizer(SCATCI_input % opts(i) % diagonalizer_choice,  &
                                               SCATCI_input % opts(i) % force_serial,         &
                                               SCATCI_input % opts(i) % contracted_mat_size,  &
                                               SCATCI_input % opts(i) % num_eigenpairs,       &
                                               matrix_elements,                               &
                                               diagonalizer,                                  &
                                               Integrals,                                     &
                                               SCATCI_input % opts(i) % hamiltonian_unit)
            call matrix_elements % construct(MAIN_HAMILTONIAN)
            call matrix_elements % set_options(SCATCI_input % opts(i))
        else
            call Integrals % write_matrix_header(SCATCI_input % opts(i) % hamiltonian_unit, &
                                                 SCATCI_input % opts(i) % contracted_mat_size)
            allocate(WriterMatrix::matrix_elements)
            call matrix_elements % construct(MAIN_HAMILTONIAN)
            call matrix_elements % set_options(SCATCI_input % opts(i))
        end if

        call matrix_elements % exclude_row_column(SCATCI_input % opts(i) % exclude_rowcolumn)

        ! We are doing a scattering calculation
        if (SCATCI_input % opts(i) % do_ci_contraction()) then
            call enrgmx % construct(SCATCI_input % opts(i), CSFS, Orbitals, Integrals)
            call enrgmx % initialize(ci_rmat)
            call master_timer % start_timer("C-Hamiltonian Build")
            call enrgmx % build_hamiltonian(matrix_elements)
            call master_timer % stop_timer("C-Hamiltonian Build")
            deallocate (ci_rmat)
        else
            call enrgms % construct(SCATCI_input % opts(i), CSFS, Orbitals, Integrals)
            call master_timer % start_timer("Target-Hamiltonian Build")
            call enrgms % build_hamiltonian(matrix_elements)
            call master_timer % stop_timer("Target-Hamiltonian Build")
        end if

        ! Lets deallocate all the csfs as well
        deallocate(CSFS)

        ! Clear the integrals, we dont need them anymore
        call Configuration_manager % finalize

        call master_memory % print_memory_report
        write (stdout, "('Matrix N is ',i8)") matrix_elements % get_matrix_size()
        write (stdout, "('Num elements is ',i14)") matrix_elements % get_size()
        call master_timer % report_timers

        if (SCATCI_input % opts(i) % diagonalization_flag /= NO_DIAGONALIZATION) then
            SCATCI_input % opts(i) % num_eigenpairs = min(SCATCI_input % opts(i) % num_eigenpairs, &
                                                          matrix_elements % get_matrix_size())
            call solutions(i) % construct(SCATCI_input % opts(i))
            call master_timer % start_timer("Diagonalization")
            call diagonalizer % diagonalize(matrix_elements,                          &
                                            SCATCI_input % opts(i) % num_eigenpairs,  &
                                            solutions(i),                             &
                                            .false.,                                  &
                                            SCATCI_input % opts(i),                   &
                                            Integrals)
            call master_timer % stop_timer("Diagonalization")
            call master_timer % report_timers
        else
            write (stdout, '(/,"Diagonalization not selected")')
            num_mat_elms = matrix_elements % get_size()

            call master_timer % start_timer("Matrix Save")
           !call matrix_elements % save(SCATCI_input % hamiltonian_unit, SCATCI_input % num_matrix_elements_per_rec)
            call master_timer % stop_timer("Matrix Save")

            if (process_grid % grank == master) call ham_pars_io_wrapper(SCATCI_input % opts(i), .true., num_mat_elms)

            write (stdout, '(/,"Parameters written to: ham_data")')
        end if

        if (process_grid % grank == master) close (SCATCI_input % opts(i) % hamiltonian_unit)

        call matrix_elements % destroy

        deallocate (diagonalizer)
        deallocate (matrix_elements)

    end do symmetry_loop

    ! use solutions in outer interface etc.
    call postprocess(SCATCI_input, solutions)
    call master_timer % stop_timer("Wall Time")
    call master_timer % report_timers

    ! release (possibly shared) memory occupied by the integral storage
    if (associated(Integrals)) then
        call Integrals % finalize
        nullify (Integrals)
    end if

    ! cleanly exit MPI
    call mpi_mod_finalize

end program mpi_scatci
