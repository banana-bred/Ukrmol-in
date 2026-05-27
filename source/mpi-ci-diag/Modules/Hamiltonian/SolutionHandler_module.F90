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

!> \brief   Solution handler module
!> \authors A Al-Refaie
!> \date    2017
!>
!> Manages the eigevalue and eigenvector output. These are passed to individual subroutines
!> by the caller; solution handler itself contains just a little basic information.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module SolutionHandler_module

    use blas_lapack_gbl,           only: blasint
    use cdenprop_defs,             only: CIvect, maxnuc
    use const_gbl,                 only: stdout
    use consts_mpi_ci,             only: SYMTYPE_D2H, SAVE_CONTINUUM_COEFFS, SAVE_L2_COEFFS, MAX_NEIG, NAME_LEN_MAX
    use precisn,                   only: wp
    use mpi_gbl,                   only: master, mpi_xermsg, mpi_mod_isend, mpi_mod_recv
    use params,                    only: ccidata, c8stars, cblank, cpoly
    use scatci_routines,           only: MOVEP, CIVIO, MOVEW
    use BaseIntegral_module,       only: BaseIntegral
    use DiagonalizerResult_module, only: DiagonalizerResult
    use Options_module,            only: Options
    use Parallelization_module,    only: grid => process_grid
    use SWEDEN_module,             only: SWEDENIntegral
    use UKRMOL_module,             only: UKRMOLIntegral

    implicit none

    !> \brief   Solution writer
    !> \authors A Al-Refaie, J Benda
    !> \date    2017 - 2019
    !>
    !> Provides a comfortable interface to eigenvector disk output.
    !> By default writes the *fort.25* file.
    !>
    type, extends(DiagonalizerResult) :: SolutionHandler
        integer               :: io_unit = 25
        integer               :: symmetry_type
        integer               :: num_eigenpairs
        real(wp), allocatable :: energy_shifts(:)
        real(wp)              :: core_energy
        real(wp)              :: energy_shift
        integer, allocatable  :: phase(:)
        integer               :: size_phase
        integer               :: vec_dimen
        integer               :: nciset
        integer               :: current_eigenvector = 0
        integer               :: print_all_eigs = 0
    contains
        procedure, public :: construct
        procedure, public :: write_header       => write_header_sol
        procedure, public :: export_header      => export_header_sol
        procedure, public :: export_eigenvalues
        procedure, public :: shift_eigenvalues
        procedure, public :: handle_eigenvalues => write_eigenvalues
        procedure, public :: handle_eigenvector => write_eigenvector
        procedure, public :: finalize_solutions => finalize_solutions_sol
        procedure, public :: read_from_file
        procedure, public :: destroy
    end type SolutionHandler

contains

    !> \brief   Set up solution handler
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Initialize internal parameters, most notably the number of the output unit.
    !>
    subroutine construct (this, opts)

        class(SolutionHandler)      :: this
        class(Options),  intent(in) :: opts

        integer :: n

        this % symmetry_type   = opts % sym_group_flag
        this % io_unit         = opts % ci_output_unit
        this % continuum_dimen = opts % seq_num_last_continuum_csf
        this % vec_dimen       = opts % contracted_mat_size
        this % num_eigenpairs  = opts % num_eigenpairs
        this % vector_storage  = opts % vector_storage_method
        this % print_all_eigs  = opts % print_flags(2)

        if (.not.(allocated(this % energy_shifts))) then
            allocate(this % energy_shifts(this % num_eigenpairs))
            n = min(this % num_eigenpairs, size(opts % energy_shifts))
            this % energy_shifts = 0.0_wp
            this % energy_shifts(1:n) = opts % energy_shifts(1:n)
        end if

    end subroutine construct


    !> \brief   Release resources
    !> \authors J Benda
    !> \date    2019
    !>
    !> Release all allocated memory (phases and the CDENPROP data vector).
    !>
    subroutine destroy (this)

        class(SolutionHandler), intent(inout) :: this

        if (allocated(this % energy_shifts)) then
            deallocate (this % energy_shifts)
        end if

        if (allocated(this % phase)) then
            deallocate (this % phase)
        end if

        call this % ci_vec % final_CV

    end subroutine destroy


    !> \brief   Start writing eigenvector record
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Opens or rewinds the output file, positions it to the beginning of the requested record
    !> and writes basic information about the solution, excluding phase information, eigenvalues
    !> and eigenvectors.
    !>
    subroutine write_header_sol (this, option, integrals)

        class(SolutionHandler)          :: this
        class(Options),      intent(in) :: option
        class(BaseIntegral), intent(in) :: integrals
        integer                         :: Nset = 1
        integer                         :: print_header
        integer                         :: nnuc, nrec, i, buf(1), n = 1, tag = 0
        integer                         :: ierror, nth, nhd(10)

        !first lets quickly grab the core_energy while we're here
        this % core_energy = integrals % get_core_energy()

        if (grid % grank /= master) return
        this % nciset = option % output_ci_set_number
        if (option % output_ci_set_number < 0) return

        ! If output unit dataset for these solutions follows after another dataset in the same unit, we need to wait
        ! for a signal from master of the MPI group responsible for writing of that dataset first, so that we do not
        ! access the unit concurrently.
        if (option % preceding_writer >= 0) then
            call mpi_mod_recv(option % preceding_writer, tag, buf, n)
        end if

        this % size_phase = min(size(option % phase_index), this % vec_dimen)

        if (.not.(allocated(this % phase))) then
            allocate(this % phase(this % vec_dimen))
            this % phase = 0
            this % phase(1:this % size_phase) = option % phase_index(1:this % size_phase)
        end if

        Nset = option % output_ci_set_number
        nth = Nset

        print_header = option % print_dataset_heading

        if (this % symmetry_type >= SYMTYPE_D2H) then

            ! Move to the correct position in the file. The subroutine movep can, in theory, change the value of "nth". However,
            ! OptionsSet should have made sure that the values of output_ci_set_number are such that this will never happen.
            ! (I.e. that we are writing always the next unit in line.) Otherwise we would risk locking of the MPI groups and/or
            ! overwriting some MPI group results by other MPI groups.
            call movep(this % io_unit, nth, ierror, print_header, stdout)

            if (ierror /= 0) then
                call mpi_xermsg('SolutionHandler_module', 'write_header_sol', &
                                ' ERROR POSITIONING FOR OUTPUT OF CI COEFFICIENTS', 1, 1)
            end if

            Nset = nth
            nnuc = integrals % get_num_nuclei()
            nrec = nnuc + this % num_eigenpairs + 1

            write (this % io_unit) c8stars, cblank, cblank, ccidata
            write (this % io_unit) Nset, nrec, option % diag_name(1:120), nnuc,  option % contracted_mat_size, &
                                   this % num_eigenpairs, option % lambda, option % spin, option % spin_z, &
                                   option % num_electrons, this % core_energy, option % num_target_sym, &
                                   (option % lambda_continuum_orbitals_target(i), &
                                    option % num_continuum_orbitals_target(i), i = 1, option % num_target_sym)

            call integrals % write_geometries(this % io_unit)

        else

            ! Move to the correct position in the file
            call movew(this % io_unit, nth, ierror, print_header, stdout)

            if (ierror /= 0) then
                call mpi_xermsg('SolutionHandler_module', 'write_header_sol', &
                                ' ERROR POSITIONING FOR OUTPUT OF CI COEFFICIENTS', 2, 1)
            end if

            nhd( :) = 0
            nhd( 2) = option % contracted_mat_size
            nhd( 3) = this % num_eigenpairs
            nhd( 4) = option % num_syms
            nhd( 8) = option % contracted_mat_size
            nhd( 9) = integrals % get_num_nuclei()
            nhd(10) = merge(1, 0, this % num_eigenpairs < option % contracted_mat_size)

            ! We are using alchemy integrals
            write (this % io_unit) nth, nhd, option % diag_name(1:120), integrals % NHE, integrals % DTNUC

        end if

        write (stdout, "(' CI data will be stored as set number',I3)") Nset

    end subroutine write_header_sol


    !> \brief   Write basic information into a CI vector
    !> \authors A Al-Refaie, J Benda
    !> \date    2017 - 2019
    !>
    !> Store general data related to the eigenvectors (like global quantum numbers, nuclear information
    !> and number of eigenvectors) to the provided CDENPROP data vector. This is used when passing
    !> data from MPI-SCATCI directly to the CDENPROP library.
    !>
    subroutine export_header_sol (this, option, integrals)

        class(SolutionHandler)           :: this
        class(Options),      intent(in)  :: option
        class(BaseIntegral), intent(in)  :: integrals

        integer :: err, i

        !first lets quickly grab the core_energy while we're here
        this % core_energy = integrals % get_core_energy()
        this % size_phase = min(size(option % phase_index), this % vec_dimen)

        if (.not.(allocated(this % phase))) then
            allocate(this % phase(this % vec_dimen), stat = err)
            if (err /= 0) then
                call mpi_xermsg('SolutionHandler_module', 'export_header_sol', &
                                'Error allocating internal phase array', 1, 1)
            end if
            this % phase = 0
            this % phase(1:this % size_phase) = option % phase_index(1:this % size_phase)
        end if

        if (this % symmetry_type >= SYMTYPE_D2H) then
            this % ci_vec % Nset  = option % output_ci_set_number
            this % ci_vec % nrec  = integrals % get_num_nuclei() + this % num_eigenpairs + 1
            this % ci_vec % name  = option % diag_name
            this % ci_vec % nnuc  = integrals % get_num_nuclei()
            this % ci_vec % nocsf = option % contracted_mat_size
            this % ci_vec % nstat = this % num_eigenpairs
            this % ci_vec % mgvn  = option % lambda
            this % ci_vec % s     = option % spin
            this % ci_vec % sz    = option % spin_z
            this % ci_vec % nelt  = option % num_electrons
            this % ci_vec % e0    = this % core_energy

            if (this % ci_vec % nnuc > maxnuc) then
                call mpi_xermsg('SolutionHandler_module', 'export_header_sol', &
                                'Increase the value of maxnuc in cdenprop_defs and recompile', 2, 1)
            end if

            select type (integrals)
                type is (UKRMOLIntegral)
                    this % ci_vec % cname (1:this % ci_vec % nnuc) = integrals % cname (1:this % ci_vec % nnuc)
                    this % ci_vec % charge(1:this % ci_vec % nnuc) = integrals % charge(1:this % ci_vec % nnuc)
                    this % ci_vec % xnuc  (1:this % ci_vec % nnuc) = integrals % xnuc  (1:this % ci_vec % nnuc)
                    this % ci_vec % ynuc  (1:this % ci_vec % nnuc) = integrals % ynuc  (1:this % ci_vec % nnuc)
                    this % ci_vec % znuc  (1:this % ci_vec % nnuc) = integrals % znuc  (1:this % ci_vec % nnuc)
                type is (SWEDENIntegral)
                    this % ci_vec % cname (1:this % ci_vec % nnuc) = integrals % cname (1:this % ci_vec % nnuc)
                    this % ci_vec % charge(1:this % ci_vec % nnuc) = integrals % charge(1:this % ci_vec % nnuc)
                    this % ci_vec % xnuc  (1:this % ci_vec % nnuc) = integrals % xnuc  (1:this % ci_vec % nnuc)
                    this % ci_vec % ynuc  (1:this % ci_vec % nnuc) = integrals % ynuc  (1:this % ci_vec % nnuc)
                    this % ci_vec % znuc  (1:this % ci_vec % nnuc) = integrals % znuc  (1:this % ci_vec % nnuc)
                class default
                    call mpi_xermsg('SolutionHandler_module', 'export_header_sol', &
                                    'Geometry export implemented only for UKRMOLIntegral, SWEDENIntegral', 3, 1)
            end select

            write (stdout, '(/,"CI VECTOR HEADER DATA:")')
            write (stdout, '("----------------------")')
            write (stdout, *) this % ci_vec % nset,  &
                              this % ci_vec % nrec,  &
                              this % ci_vec % NAME,  &
                              this % ci_vec % nnuc,  &
                              this % ci_vec % nocsf, &
                              this % ci_vec % nstat, &
                              this % ci_vec % mgvn,  &
                              this % ci_vec % s,     &
                              this % ci_vec % sz,    &
                              this % ci_vec % nelt,  &
                              this % ci_vec % e0
            do i = 1, this % ci_vec % nnuc
                write (stdout, *) this % ci_vec % cname(i), &
                                  this % ci_vec % xnuc(i),  &
                                  this % ci_vec % ynuc(i),  &
                                  this % ci_vec % znuc(i),  &
                                  this % ci_vec % charge(i)
            end do
            write (stdout, '("----------------------")')
        else
            call mpi_xermsg('SolutionHandler_module', 'export_header_sol', &
                            'Output of ALCHEMY header not implemented yet', 4, 1)
        end if

        write (stdout, '(" CI data will be stored in memory for use in subsequent subroutines")')

    end subroutine export_header_sol


    !> \brief   Shift eigenvalues according to the NAMELIST input ESHIFT
    !> \authors T Meltzer
    !> \date    2020
    !>
    !> Shift eigenvalues according to the NAMELIST input ESHIFT. This subroutine is called only
    !> when some of the energy_shifts is non-zero.
    !>
    subroutine shift_eigenvalues (this, output_eigenvalues, num_eigenpairs)

        class(SolutionHandler)  :: this
        integer,  intent(in)    :: num_eigenpairs
        real(wp), intent(inout) :: output_eigenvalues(num_eigenpairs)
        integer :: ido, jdo, neig

        if (grid % grank /= master) return
        if (this % nciset < 0) return

        output_eigenvalues = output_eigenvalues + this % energy_shifts

        neig = this % num_eigenpairs
        if (this % print_all_eigs == 0) then
            neig = min(this % num_eigenpairs, MAX_NEIG)
        end if

        write (stdout, '("ENERGY SHIFTS TO BE APPLIED")')
        do ido = 1, neig, 5
            write (stdout, "(5F20.10)") (this % energy_shifts(jdo), jdo = ido, min(ido + 4, neig))
        end do

    end subroutine shift_eigenvalues


    !> \brief   Save Hamiltonian eigen-energies to disk file
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Continues the output operation commenced in \ref write_header_sol by writing the phase information,
    !> eigenvalues and diagonals(?).
    !>
    subroutine write_eigenvalues (this, eigenvalues, diagonals, num_eigenpairs, vec_dimen)

        class(SolutionHandler) :: this
        integer,  intent(in)   :: num_eigenpairs, vec_dimen
        real(wp), intent(in)   :: eigenvalues(num_eigenpairs), diagonals(vec_dimen)
        real(wp)               :: output_eigenvalues(num_eigenpairs)
        integer :: ido, jdo, neig

        if (grid % grank /= master) return
        if (this % nciset < 0) return

        output_eigenvalues = eigenvalues

        neig = this % num_eigenpairs
        if (this % print_all_eigs == 0) then
            neig = min(this % num_eigenpairs, MAX_NEIG)
        end if

        if (sum(abs(this % energy_shifts)) .gt. 0_wp) then
            write (stdout, '("ORIGINAL EIGEN-ENERGIES")')
            do ido = 1, neig, 5
                write (stdout, "(16X,5F20.10)") (output_eigenvalues(jdo) + this % core_energy, jdo = ido, min(ido + 4, neig))
            end do
            call this % shift_eigenvalues(output_eigenvalues, num_eigenpairs)
        end if

        write (stdout, '("EIGEN-ENERGIES")')
        do ido = 1, neig, 5
            write (stdout, "(16X,5F20.10)") (output_eigenvalues(jdo) + this % core_energy, jdo = ido, min(ido + 4, neig))
        end do
        write (this % io_unit) this % phase, output_eigenvalues, diagonals

    end subroutine write_eigenvalues


    !> \brief   Write eigevalues and phases to supplied arrays
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Copies (a subset of) eigenvalues and phases stored in the calling objecst to supplied
    !> arrays. This is used when preparing diagonalization data for processing in CDENPROP.
    !>
    subroutine export_eigenvalues (this, eigenvalues, diagonals, num_eigenpairs, vec_dimen, ei, iphz)

        class(SolutionHandler) :: this
        integer,  intent(in)   :: num_eigenpairs, vec_dimen
        real(wp), intent(in)   :: eigenvalues(num_eigenpairs), diagonals(vec_dimen)
        real(wp), allocatable  :: ei(:)
        real(wp)               :: output_eigenvalues(num_eigenpairs)
        integer,  allocatable  :: iphz(:)
        integer                :: ido, jdo, err

        output_eigenvalues = eigenvalues

        if (sum(abs(this % energy_shifts)) .gt. 0_wp) then
            call this % shift_eigenvalues(output_eigenvalues, num_eigenpairs)
        end if

        if (allocated(ei)) deallocate(ei)
        allocate(ei, source = output_eigenvalues, stat = err)
        if (err /= 0) then
            call mpi_xermsg('SolutionHandler_module', 'export_eigenvalues', 'Error allocating ei', size(output_eigenvalues), 1)
        end if

        if (allocated(iphz)) deallocate(iphz)
        allocate(iphz, source = this % phase, stat = err)
        if (err /= 0) then
            call mpi_xermsg('SolutionHandler_module', 'export_eigenvalues', 'Error allocating iphz', size(this % phase), 1)
        end if

    end subroutine export_eigenvalues


    !> \brief   Save Hamiltonian eigen-vectors to disk file
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Completes the output operation commenced in \ref write_header_sol by writing the eigenvectors
    !> (one per a call).
    !>
    subroutine write_eigenvector (this, eigenvector, vec_dimen)

        class(SolutionHandler) :: this
        integer,  intent(in)   :: vec_dimen
        real(wp), intent(in)   :: eigenvector(vec_dimen)
        integer                :: start_dimen, end_dimen, ido, j

        if (grid % grank /= master) return
        if (this % nciset < 0) return

        this % current_eigenvector = this % current_eigenvector + 1

        start_dimen = 1
        end_dimen = this % vec_dimen

        if (iand(this % vector_storage, SAVE_L2_COEFFS) == 0) then
            end_dimen = this % continuum_dimen
        end if
        if (iand(this % vector_storage, SAVE_CONTINUUM_COEFFS) == 0) then
            start_dimen = this % continuum_dimen + 1
        end if

        write (this % io_unit) this % current_eigenvector, (eigenvector(j), j = start_dimen, end_dimen)

    end subroutine write_eigenvector


    !> \brief   Import eigensolution from file
    !> \authors J Benda
    !> \date    2022 - 2025
    !>
    !> Instead of calculating the eigenpairs, simple read them from a file written earlier.
    !> Store them in a distributed CIvect structure for processing in CDENPROP.
    !>
    subroutine read_from_file (this, option)

        class(SolutionHandler)          :: this
        class(Options),      intent(in) :: option
        integer                         :: print_header, ierror, nset, nrec, nnuc, nth, nocsf, nstat, mgvn, nelt, i, j, m
        integer(blasint)                :: one = 1
        real(wp)                        :: e0, s, sz
        real(wp), allocatable           :: vec(:), ei(:), phase(:)
        character(len=NAME_LEN_MAX)     :: diag_name

        nset = option % output_ci_set_number
        nth = nset
        print_header = option % print_dataset_heading

        this % ci_vec % CV_is_scalapack = .true.
        this % ci_vec % blacs_context   = grid % gcntxt
        this % ci_vec % myrow           = grid % mygrow
        this % ci_vec % mycol           = grid % mygcol
        this % ci_vec % nprow           = grid % gprows
        this % ci_vec % npcol           = grid % gpcols

        call this % ci_vec % init_CV(this % vec_dimen, this % num_eigenpairs)

        allocate (this % ci_vec % ei(this % num_eigenpairs), ei(this % num_eigenpairs))
        allocate (this % ci_vec % iphz(this % num_eigenpairs), phase(this % num_eigenpairs))
        allocate (vec(this % vec_dimen))

        if (grid % grank == master) then
            call movep(this % io_unit, nth, ierror, print_header, stdout)

            if (ierror /= 0) then
                call mpi_xermsg('SolutionHandler_module', 'read_from_file', 'Failed to read eigensolution', 1, 1)
            end if

            ! read header to get the number of nuclei
            read (this % io_unit)
            read (this % io_unit) nset, nrec, diag_name, nnuc, nocsf, nstat, mgvn, s, sz, nelt, e0

            ! skip nuclear information
            do i = 1, nnuc
                read (this % io_unit)
            end do

            ! read the phases and eigenenergies
            read (this % io_unit) phase, ei
        else
            phase = 0
            ei = 0
        end if

#if defined(usempi) && defined(scalapack)
        ! broadcast the data to other tasks
        call dgsum2d(grid % gcntxt, 'all', ' ', this % vec_dimen, one, phase, this % vec_dimen, -one, -one)
        call dgsum2d(grid % gcntxt, 'all', ' ', this % vec_dimen, one, ei, this % vec_dimen, -one, -one)
#endif

        this % ci_vec % iphz = phase
        this % ci_vec % ei = ei

        ! read all requested eigenvectors
        do j = 1, this % num_eigenpairs
            vec = 0
            ! master writes the eigenvector
            if (grid % grank == master) then
                read (this % io_unit) m, vec
            end if
#if defined(usempi) && defined(scalapack)
            ! master broadcasts the eigenvector
            call dgsum2d(grid % gcntxt, 'all', ' ', this % vec_dimen, one, vec, this % vec_dimen, -one, -one)
#endif
            ! all processes independently include the elements into their storage
            do i = 1, this % vec_dimen
                call this % ci_vec % set_CV_element(vec(i), i, j)
            end do
        end do

    end subroutine read_from_file


    !> \brief   Finalize writing the eigenvector disk file
    !> \authors J Benda
    !> \date    2019
    !>
    !> To avoid concurrent, conflicting accesses to the disk file in case that multiple datasets are to be written to it,
    !> by different processes, MPI-SCATCI uses a semaphore system. Once a process finishes writing its results into a dataset,
    !> it notifies the next process in line that wants to write to the same file unit (but into the subsequent dataset).
    !>
    subroutine finalize_solutions_sol (this, option)

        class(SolutionHandler)          :: this
        class(Options),      intent(in) :: option

        integer :: buf(1), tag = 0, n = 1

        flush(this % io_unit)
        close(this % io_unit)

        if (option % following_writer >= 0) then
            call mpi_mod_isend(option % following_writer, buf, tag, n)
        end if

    end subroutine finalize_solutions_sol

end module SolutionHandler_module
