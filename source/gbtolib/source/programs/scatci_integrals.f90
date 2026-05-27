! Copyright 2020
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
program scatci_integrals

   use iso_fortran_env,      only: iostat_end
   use const_gbl,            only: line_len, sym_op_nam_len, stdout, set_verbosity_level
   use atomic_basis_gbl,     only: atomic_orbital_basis_obj, bto_shell_data_obj
   use free_scattering_gbl,  only: free_scattering
   use integral_storage_gbl, only: integral_options_obj, integral_storage_obj, p2d_array_obj
   use molecular_basis_gbl,  only: molecular_orbital_basis_obj, sym_ortho_io
   use mpi_gbl,              only: mpi_xermsg, mpi_mod_start, mpi_mod_print_info, mpi_mod_finalize
   use precisn_gbl,          only: cfp
   use symmetry_gbl,         only: symmetry_obj

   implicit none

   !Limits on the size of the namelist arrays below.
   integer, parameter :: max_exp = 20
   integer, parameter :: max_cont_l = 15

   !MPI initialization parameters
   logical, parameter :: allow_shared_memory = .false.  ! not yet properly implemented

   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
   type(molecular_orbital_basis_obj) :: molecular_orbital_basis
   type(BTO_shell_data_obj) :: BTO_shell_data
   type(integral_options_obj) :: integral_options
   type(integral_storage_obj), target :: atomic_1el_integral_storage, atomic_2el_integral_storage, molecular_integral_storage
   type(p2d_array_obj), target :: ao_1el_integrals, ao_2el_integrals
   type(sym_ortho_io) :: sym_ortho(8)

   integer :: err, j, k, number_of_continuum_shells, no_cont_exps(0:max_cont_l)
   integer :: n_cont_cgto, n_cont_bto, n_cont_cgto_orbs(8), n_cont_bto_orbs(8)
   integer :: n_PCO_cgto, number_of_PCO_shells, n_PCO_cgto_orbs(8), n_tot_orbs(8)=0 ! DDL Number of PCO CGTO functions, num of PCO atomic shells, and num of PCO orbs per irr, and num total orbitals per irr
   logical :: include_btos
   logical, allocatable :: function_sym(:), to_delete(:)
   logical :: include_PCOs ! DDL PCO inclusion flag
   real(kind=cfp), allocatable :: overlap_matrix(:,:), amplitudes(:,:)
   real(kind=cfp) :: PCO_exponents(1:max_exp,0:max_cont_l) = 0.0_cfp ! DDL value of PCO gto exponents

   !namelist variables
   character(len=line_len) :: molden_file = '', mo_integrals_file_name = './moints', ao_integrals_file_name = './aoints'
   character(len=line_len) :: scratch_directory = '', basis_input = '', can_mo_integrals_file_name = './moints_canonical'
   integer :: nob(1:8) = 0, select_orbitals_by = 1, no_sym_op = 0, min_l = -1, max_l = -1, nE = -1, verbosity = 1
   integer :: bspline_indices(2,0:max_cont_l) = -1, min_bspline_l = -1, max_bspline_l = -1, bspline_order = -1, no_bsplines = -1
   integer :: two_particle_integrals_method = 1  ! 1 = transformation of AO integrals, 2 = direct computation from Poission equation
   integer :: mixed_ints_method = -1, max_l_legendre_1el = -1, max_l_legendre_2el = -1
   integer :: molecular_2el_algorithm = 0  ! 0 = auto, 1 = sparse, other = dense
   integer :: alpha_or_beta = 1  ! 0 = use both spins, 1 = use only alpha orbitals, 2 = use only beta orbitals
   integer :: num_PCOs(0:max_cont_l) = 0, min_PCO_l=-1, max_PCO_l=-1 ! DDL number of PCOs for a given l, min and max l for PCOs
   integer :: HF_max_iters = 10
   integer :: max_property_l = 2 ! max l for Property integrals (for multipole expansion of density)
   integer :: aoints_index_scheme = 1  ! 1 = standard, 2 = compact
   character(len=sym_op_nam_len) :: sym_op(1:3) = (/'  ','  ','  '/)
   real(kind=cfp) :: del_thrs(1:8) = (/-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp/)
   real(kind=cfp) :: exponents(1:max_exp,0:max_cont_l) = 0.0_cfp
   real(kind=cfp) :: a = -1.0_cfp, min_energy = 0.0_cfp, max_energy = 1.0_cfp, bspline_grid_start = -1.0_cfp
   real(kind=cfp) :: max_ijrs_size = -1.0_cfp, delta_r1 = 0.25_cfp
   real(kind=cfp) :: PCO_alpha0(0:max_cont_l) = -1.0_cfp, PCO_beta(0:max_cont_l) = -1.0_cfp, PCO_gto_thrs(0:max_cont_l)=-1.0_cfp ! DDL PCO parameters alpha0 and beta, and PCO to continuum<F3> gto thresholds
   real(kind=cfp) :: PCO_del_thrs(1:8) = -1.0_cfp ! DDL PCO deletion thresholds per irr
   real(kind=cfp) :: dipole_damp_factor = 0.0_cfp ! when non-zero the dipole properties are calculated with the radial part of the dipole operator: r * exp(-dipole_damp_factor*r)
   real(kind=cfp) :: HF_convergence = 1E-10_cfp
   logical :: use_ylm_scratch = .true., use_aoints_scratch = .false., use_transformation_scratch = .false.
   logical :: run_free_scattering = .false., save_ao_integrals_to_disk = .false., do_two_particle_integrals = .true.
   logical :: two_p_continuum = .false., use_spherical_cgto_alg = .true., redirect_master = .true.
   logical :: check_target_target_orbital_overlaps = .true., check_target_continuum_orbital_overlaps = .true.
   logical :: check_continuum_continuum_orbital_overlaps = .true., calc_radial_densities = .false.
   logical :: preorthogonalize_continuum = .false., print_1el_ints = .false., print_2el_ints = .false.
   logical :: ortho_continuum_against_all_tgt_orbs = .false., qmoln = .false. ! added back qmoln logical for controlling QEC functionality (BC)
   logical :: reortho_target_orbitals = .true. ! controls execution of Gramm-Schmidt reorthogonalization of the MOLDEN target orbitals
   logical :: construct_canonical_continuum = .false.
   logical :: keep_ao_integrals = .false.
   logical :: only_construct_fock_blocks = .false.
   logical :: canonize_virtuals_instead_of_continuum = .false.
   logical :: calculate_bbb_ints = .false. !must be .true. in case of 2-el integrals by Poisson method
   logical :: calculate_overlap_ints = .true., calculate_kinetic_energy_ints = .true., calculate_property_ints = .true.
   !if 1-el Hamiltonian ints. are required, then nuclear attractrion integrals must be selected too. This is done automatically below.
   logical :: calculate_nuclear_attraction_ints = .true., calculate_one_el_hamiltonian_ints = .true.

   namelist /target_data/ molden_file, nob, no_sym_op, sym_op, a, select_orbitals_by, alpha_or_beta
   namelist /pco_data/ PCO_alpha0, PCO_beta, num_PCOs, min_PCO_l, max_PCO_l,PCO_gto_thrs,PCO_del_thrs  ! DDL PCO input parameters
   namelist /continuum_data/ del_thrs, exponents, min_l, max_l, run_free_scattering, min_energy, max_energy, nE, bspline_indices, &
                             min_bspline_l, max_bspline_l, bspline_grid_start, bspline_order, no_bsplines
   namelist /process_control/ max_ijrs_size, ao_integrals_file_name, mo_integrals_file_name, save_ao_integrals_to_disk, &
                              do_two_particle_integrals, use_spherical_cgto_alg, check_target_target_orbital_overlaps, &
                              check_target_continuum_orbital_overlaps, check_continuum_continuum_orbital_overlaps, &
                              two_p_continuum, mixed_ints_method, max_l_legendre_1el, max_l_legendre_2el, calc_radial_densities, &
                              preorthogonalize_continuum, print_1el_ints, print_2el_ints, scratch_directory, delta_r1, &
                              verbosity, redirect_master, ortho_continuum_against_all_tgt_orbs, molecular_2el_algorithm, qmoln, &
                              basis_input, dipole_damp_factor, construct_canonical_continuum, keep_ao_integrals, HF_convergence, &
                              HF_max_iters, only_construct_fock_blocks, canonize_virtuals_instead_of_continuum, &
                              aoints_index_scheme, use_transformation_scratch, &
                              calculate_bbb_ints, max_property_l, two_particle_integrals_method, use_ylm_scratch, &
                              use_aoints_scratch, calculate_overlap_ints, calculate_kinetic_energy_ints, calculate_property_ints, &
                              calculate_nuclear_attraction_ints, calculate_one_el_hamiltonian_ints


      !this routine is contained in this program at the bottom: read-in the namelist variables and check for sanity of the input
      !idat = 5
      call process_namelist

      !initialize MPI and output redirections
      call mpi_mod_start(.not.redirect_master, allow_shared_memory)
      call set_verbosity_level(verbosity)
      call mpi_mod_print_info(stdout)

      !generate atomic and molecular basis sets from the input, or read them from file; also evaluates atomic 1-el integrals
      call setup_basis_sets

      !save the atomic basis and the molecular orbital basis to disk
      call atomic_orbital_basis%write(mo_integrals_file_name)
      call molecular_orbital_basis%write(mo_integrals_file_name)

      !
      ! TRANSFORM THE 1-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE MOLECULAR ORBITALS:
      !

      !Describe where the transformed AO->MO integrals will be stored :
      err = molecular_integral_storage%init(disk=mo_integrals_file_name)
      if (err /= 0) then
         call mpi_xermsg('main', 'main', 'error initializing the target molecular_integral_storage', err, 1)
      end if

      molecular_orbital_basis%ao_integral_storage => atomic_1el_integral_storage !point to the storage for the atomic integrals
      call molecular_orbital_basis%one_electron_integrals(molecular_integral_storage,integral_options)
     

      !
      ! RUN FREE-POTENTIAL SCATTERING:
      !

      if (run_free_scattering .and. number_of_continuum_shells > 0) then
         call free_scattering(molecular_integral_storage, molecular_orbital_basis, a, min_energy, max_energy, nE)
      end if

      !
      ! CALCULATE THE 2-ELECTRON INTEGRALS:
      !

      if (construct_canonical_continuum) then
         ! always calculate 2-electron integrals when constructing canonical continuum
         call calculate_2el_integrals
         call molecular_orbital_basis%construct_canonical_continuum(atomic_1el_integral_storage, atomic_2el_integral_storage, &
                         molecular_integral_storage, integral_options, molecular_2el_algorithm, n_cont_bto, n_cont_cgto, &
                         mo_integrals_file_name, can_mo_integrals_file_name)
         call atomic_1el_integral_storage%final
      else
         ! conditionally calculate 2-electron integrals, but release atomic storage first for maximum memory
         call atomic_1el_integral_storage%final
         if (do_two_particle_integrals) then
            call calculate_2el_integrals
         end if
      endif

      !
      ! CALCULATE RADIAL CHARGE DENSITIES FOR ALL ORBITALS:
      !

      if (calc_radial_densities .and. integral_options % a > 0.0_cfp) then
         call molecular_orbital_basis % radial_charge_density(integral_options % a, 0.0_cfp, &
              integral_options % a, 0.1_cfp, .true., amplitudes)
      end if

      if (do_two_particle_integrals) then
         call atomic_2el_integral_storage%final
      end if


      call molecular_integral_storage%final

      call mpi_mod_finalize

contains

   subroutine process_namelist

      use pco_gbl, only: generate_PCO_exponents

      integer :: i, j, idat, eof
         
 
         call open_inp_file(idat)

         read (idat, nml = target_data)

         if (molden_file == '') then
            call mpi_xermsg('test', 'process_namelist', 'No Molden input file specified.', 1, 1)
         end if

         if (no_sym_op < 0) then
            call mpi_xermsg('test', 'process_namelist', 'The number of symmetry operations cannot be smaller than 0.', 2, 1)
         end if

         if (all(nob == 0)) then
            call mpi_xermsg('test', 'process_namelist', 'Target orbitals (NOB) must not be all equal to zero.', 3, 1)
         end if

         if (select_orbitals_by <= 0 .or. select_orbitals_by > 3) then
            call mpi_xermsg('test', 'process_namelist', 'Error in input: select_orbitals_by was out of the range [1;3].', 4, 1)
         end if

         include_PCOs=.false.
         number_of_PCO_shells=0
         ! Is there PCO data in namelist pco_data?
         read(idat,nml=pco_data,iostat=eof)
         rewind(idat)
         if(eof .eq. 0 ) then ! pco_data namelist read in, PCO Data check
            if (min_PCO_l < 0 .and. max_PCO_l >= 0 .or. max_PCO_l < 0 .and. min_PCO_l >= 0) then
               call mpi_xermsg('test', 'process_namelist', &
                               'When one of min_PCO_l/max_PCO_l is given the other one must be specified too.', 17, 1)
            else if (min_PCO_l >= 0 .and. max_PCO_l < min_PCO_l) then
               call mpi_xermsg('test', 'process_namelist', 'max_PCO_l must not be smaller than min_PCO_l.', 18, 1)
            else if (min_PCO_l >= 0 .and. max_PCO_l >= min_PCO_l) then
               do i=min_PCO_l, max_PCO_l
                  if(PCO_alpha0(i).le.0.0) call mpi_xermsg('test', 'process_namelist', &
                                                           'PCO_alpha0 components must all be greater than zero.', 19, 1)
                  if(PCO_beta(i).le.0.0) call mpi_xermsg('test', 'process_namelist', &
                                                         'PCO_beta components must all be greater than zero.', 20, 1)
                  if(num_PCOs(i).le.0.0) call mpi_xermsg('test', 'process_namelist', &
                                                         'num_PCOs components must all be greater than zero.', 21, 1)
                  number_of_PCO_shells=number_of_PCO_shells + num_PCOs(i)
               end do
               if( number_of_PCO_shells > 0 ) then
                   include_PCOs=.true.
                   call generate_PCO_exponents( min_PCO_l,max_PCO_l,maxval(num_PCOs),num_PCOs(min_PCO_l:max_PCO_l), &
                        PCO_alpha0(min_PCO_l:max_PCO_l),PCO_beta(min_PCO_l:max_PCO_l), &
                        PCO_exponents(1:maxval(num_PCOs),min_PCO_l:max_PCO_l), &
                        PCO_gto_thrs(min_PCO_l:max_PCO_l))
               else
                   call mpi_xermsg('test', 'process_namelist', &
                                   'No PCO shells contained in PCO data.', 22, 0)
               endif
            end if
         else if (eof .eq. iostat_end) then
            !write(stdout,'("No Pseudo-Continuum Orbital data name list, pco_data, not using PCOs. ",i0)') eof
         else ! If there is an issue not related to end of file then check with a non iostat call to read the nml
            read(idat,nml=pco_data)
         end if


         read(idat,nml=continuum_data)

         if (min_l < 0 .and. max_l >= 0 .or. max_l < 0 .and. min_l >= 0) then
            call mpi_xermsg('test', 'process_namelist', &
                            'When one of min_l/max_l is given the other one must be specified too.', 5, 1)
         end if

         if (min_l >= 0 .and. max_l < min_l) then
            call mpi_xermsg('test', 'process_namelist', 'max_l must not be smaller than min_l.', 4, 1)
         end if

         no_cont_exps(:) = 0
         number_of_continuum_shells = 0

         if (min_l .ge. 0 .and. max_l .ge. min_l) then
            do i=min_l,max_l
               j = 0
               do
                  j = j + 1
                  if (exponents(j,i) .le. 0.0_cfp) exit
               enddo
               j = j - 1

               if (j .le. 0) call mpi_xermsg('test','process_namelist','No exponents given in L shell < max_L.',6,1)
               no_cont_exps(i) = j
               number_of_continuum_shells = number_of_continuum_shells + j
            enddo !i
         endif

         include_btos = .false.

         if (min_bspline_l < 0 .and. max_bspline_l >= 0 .or. max_bspline_l < 0 .and. min_bspline_l >= 0) then
            call mpi_xermsg('test', 'process_namelist', &
                            'When one of min_bspline_l/max_bspline_l is given the other one must be specified too.', 7, 1)
         end if

         if (min_bspline_l >= 0 .and. max_bspline_l < min_bspline_l) then
            call mpi_xermsg('test', 'process_namelist', 'max_bspline_l must not be smaller than min_bspline_l.', 8, 1)
         end if

         if (min_bspline_l .ge. 0 .and. max_bspline_l .ge. min_bspline_l) then
            do i=min_bspline_l,max_bspline_l
               if (bspline_indices(2,i) < bspline_indices(1,i)) then
                  call mpi_xermsg('test', 'process_namelist', 'Incompatible values of bspline_indices.', 9, 1)
               end if
               if (bspline_indices(1,i) <= 0 .or. bspline_indices(2,i) <= 0) then
                  call mpi_xermsg('test', 'process_namelist', 'Values in bspline_indices are out of range.', 10, 1)
               end if
               number_of_continuum_shells = number_of_continuum_shells + bspline_indices(2,i)-bspline_indices(1,i)+1
               if (bspline_indices(2,i)-bspline_indices(1,i)+1 > 0) include_btos = .true.
            enddo !i
         endif

         if (include_btos) then
            if (bspline_grid_start < 0.0_cfp) then
               call mpi_xermsg('test', 'process_namelist', 'The starting point of the radial BTO grid is < 0.0_cfp.', 11, 1)
            end if
            if (bspline_order <= 0) then
               call mpi_xermsg('test', 'process_namelist', 'The order of the radial B-splines is <= 0.', 12, 1)
            end if
            if (no_bsplines <= 0) then
               call mpi_xermsg('test', 'process_namelist', 'The number of radial B-splines is <= 0.', 13, 1)
            end if
         endif

         read(idat,nml=process_control)

         if (ao_integrals_file_name == '') then
            call mpi_xermsg('test', 'process_namelist', &
                            'The input variable ao_integrals_file_name is empty but must not be empty.', 14, 1)
         end if

         if (mo_integrals_file_name == '') then
            call mpi_xermsg('test', 'process_namelist', &
                            'The input variable mo_integrals_file_name is empty but must not be empty.', 15, 1)
         end if

         if (ao_integrals_file_name == mo_integrals_file_name) then
            call mpi_xermsg('test', 'process_namelist', &
                            'The input variables ao/mo_integrals_file_name must be different.', 16, 1)
         end if

         close(idat)

   end subroutine process_namelist


   subroutine open_inp_file(idat_out)

      implicit none
      integer, intent(out) :: idat_out

      character(len=1024) :: filename_inp
      logical :: ex
      integer :: length, stat

         ! prefer input filename from command line, default to 'inp' in the current working directory
         if (command_argument_count() > 0) then
            call get_command_argument(1, filename_inp, length, stat)
            if (stat /= 0) then
               call mpi_xermsg('scatci_integrals', 'process_namelist', 'Failed to get filename from command line.', 1, 1)
            end if
         else
            filename_inp = 'inp'
         end if

         ! check that the input file exists
         inquire (file = trim(filename_inp), exist = ex)
         if (.not. ex) then
            call mpi_xermsg('scatci_integrals', 'open_inp_file', 'The file "' // trim(filename_inp) // '" does not exist!', 1, 1)
         end if

         ! open the input file
         open (newunit = idat_out, file = trim(filename_inp), form = 'formatted', status = 'old', action = 'read')

   end subroutine open_inp_file


   subroutine setup_basis_sets

      if (basis_input == '') then
         call construct_basis_sets
      else
         call read_basis_sets
      end if

   end subroutine setup_basis_sets


   subroutine read_basis_sets

      ! read atomic & molecular basis
      call atomic_orbital_basis%read(basis_input)
      molecular_orbital_basis % ao_basis => atomic_orbital_basis
      call molecular_orbital_basis%read(basis_input)

      ! also perform 1el integrals (otherwise done in construct_basis_sets)
      call setup_integral_options
      call calculate_1el_integrals

   end subroutine read_basis_sets


   subroutine construct_basis_sets

      use atomic_basis_gbl,       only: CGTO_shell_data_obj
      use basis_data_generic_gbl, only: orbital_data_obj
      use molden_gbl,             only: molden_input_obj
      use pco_gbl,                only: rm_cont_gt_pco
      use symmetry_gbl,           only: geometry_obj, determine_pg
      use const_gbl,              only: level1

      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)
      type(orbital_data_obj),    allocatable :: orbital_data(:)

      type(molden_input_obj) :: molden_input
      type(geometry_obj)     :: geometry
      type(symmetry_obj)     :: symmetry_data

      integer :: io, ntot(8), nob_ortho(8), number_of_target_shells, pg, n_sym, n, i
      logical :: ortho_against_all_tgt_orbs
      integer :: idat

      !check existence of output files
      call abort_if_file_exists(trim(ao_integrals_file_name))
      call abort_if_file_exists(trim(mo_integrals_file_name))

      !open and analyze the Molden input file
      io = 1
      call molden_input%init(molden_file, io, alpha_or_beta, determine_pg(no_sym_op, sym_op))
      call molden_input%get_mo_num(ntot)
      write (stdout, '(/,"The number of orbitals per symmetry on the Molden file is: ",8i5,/)') ntot(1:8)

      !read in: the target geometry, CGTO basis data and all molecular orbitals.
      call molden_input%read(geometry%nucleus, CGTO_shell_data, orbital_data)

      if (ortho_continuum_against_all_tgt_orbs .and. (number_of_continuum_shells > 0)) then
         !The continuum will be orthogonalized against the full set of target orbitals.
         ortho_against_all_tgt_orbs = .true.
         nob_ortho = ntot
      else
         !The continuum will be orthogonalized as usual against the number of nob(:) orbitals from each symmetry.
         ortho_against_all_tgt_orbs = .false.
         nob_ortho = nob
      end if

      !close the Molden file
      call molden_input%final()

      geometry%no_nuc = size(geometry%nucleus)
      number_of_target_shells = size(CGTO_shell_data)

      !set-up the symmetry input
      geometry%no_sym_op = no_sym_op
      geometry%sym_op(1:3) = sym_op(1:3)
      geometry%use_symmetry = .true.

      !if this is a scattering or pseudo-state calculation add the scattering centre to the list of nuclei.
      if (number_of_PCO_shells > 0 .or. number_of_continuum_shells > 0) call geometry%add_scattering_centre

      err = symmetry_data%init(geometry)
      if (err /= 0) call mpi_xermsg('scatci_integrals', 'main', 'Symmetry initialization failed.', err, 1)

      !Get the ID for this point group symmetry
      pg = symmetry_data%get_pg()

      !Get the number of IRRs for this point group symmetry
      n_sym = symmetry_data%get_no_irrep(pg)

      !Set the point group symmetry ID for orbital sets in all symmetries
      orbital_data(1:n_sym)%point_group = pg

      !Keep only the requested number of target molecular orbitals:
      do i = 1, n_sym
         if (nob(i) > 0 .and. i > size(orbital_data)) then
            print *, i, nob(i)
            call mpi_xermsg('scatci_integrals', 'main', &
                            'Target molecular orbitals for the requested symmetry are not on the Molden file.', 1, 1)
         end if
         if (.not.(ortho_against_all_tgt_orbs)) then
            ! Keep only nob(i) orbitals in this symmetry: the criterion for selection of the orbitals is either orbital index
            ! (select_orbitals_by==1) or orbital energy (select_orbitals_by==2) or occupation (select_orbitals_by==3)
            call orbital_data(i)%keep_first_n_orbitals(nob(i), select_orbitals_by)
         end if
      end do

      !
      ! ADD THE ATOMIC SHELLS TO THE BASIS:
      !

      !todo temp
      !number_of_target_shells = 0

      if (include_PCOs .and. maxval(no_cont_exps ) > 0) then
         call rm_cont_gt_pco(PCO_gto_thrs(min_PCO_l:max_PCO_l), &
                             exponents, min_l, max_l, no_cont_exps, maxval(no_cont_exps), number_of_continuum_shells, &
                             PCO_exponents, min_PCO_l, max_PCO_l, num_PCOs, maxval(num_PCOs))
      end if

      !READ THE ECPs from the namelists &ECP
      call open_inp_file(idat)

      do i = 1, geometry%no_nuc
         call geometry%nucleus(i)%ecp%read_ecp_nml(idat)
      enddo

      n = number_of_target_shells + number_of_continuum_shells + number_of_PCO_shells
      err = atomic_orbital_basis%init(n, geometry)

      !add the target CGTO shells
      do i = 1, number_of_target_shells
         call atomic_orbital_basis%add_shell(CGTO_shell_data(i))
      end do

      !if constructing canonical continuum, make the core basis better fock-orthogonal
      if (construct_canonical_continuum) then
         call restricted_hartree_fock(geometry,n_sym,CGTO_shell_data,orbital_data)
      endif

      ! Add the PCOs to the atomic_orbital_basis
      n_PCO_cgto = 0
      if (include_PCOs) then
          call atomic_orbital_basis % add_cms_gtos(min_PCO_l, max_PCO_l, maxval(num_PCOs), num_PCOs(min_PCO_l:max_PCO_l), &
                                                   PCO_exponents(1:maxval(num_PCOs),min_PCO_l:max_PCO_l), n_PCO_cgto, .false.)
      end if

      ! Add the continuum CGTOs to the atomic_orbital_basis
      n_cont_cgto = 0
      if (maxval(no_cont_exps) > 0) then
          call atomic_orbital_basis % add_cms_gtos(min_l, max_l, maxval(no_cont_exps), no_cont_exps(min_l:max_l), &
                                                   exponents(1:maxval(no_cont_exps),min_l:max_l), n_cont_cgto, .true.)
      end if

      !add the BTO continuum functions
      n_cont_bto = 0
      if (include_btos) then

         !Set-up the grid of radial B-splines
         call BTO_shell_data % bspline_grid % init_grid(bspline_grid_start, a, -1.0_cfp, bspline_order, no_bsplines)

         do j = min_bspline_l, max_bspline_l
            if (j < 0) exit

            BTO_shell_data%l = j
            BTO_shell_data%number_of_functions = 2*BTO_shell_data%l + 1
            BTO_shell_data%non_zero_at_boundary = .true.

            if (bspline_grid_start > 0.0_cfp) then
               if (bspline_indices(1,j) < BTO_shell_data%bspline_grid%ind_0_der) then
                  print *, j, bspline_indices(1,j), BTO_shell_data%bspline_grid%ind_0_der
                  call mpi_xermsg('scatci_integrals', 'main', &
                                  'B-spline with a non-zero derivative at r1 has been included in the basis.', 1, 1)
               end if
            end if

            do i=bspline_indices(1,j),bspline_indices(2,j)
               BTO_shell_data%bspline_index = i

               call atomic_orbital_basis%add_shell(BTO_shell_data)
               n_cont_bto = n_cont_bto + BTO_shell_data%number_of_functions
            end do !i
         end do !j
      end if

      !
      ! ADD THE MOLECULAR ORBITALS TO THE BASIS:
      !

      molecular_orbital_basis%ao_basis => atomic_orbital_basis
      err = molecular_orbital_basis%init(n_sym,geometry)
      if (err .ne. 0) stop "error initializing molecular orbital basis"

      i = number_of_target_shells
      j = number_of_target_shells + number_of_PCO_shells
      k = number_of_target_shells + number_of_PCO_shells + number_of_continuum_shells
      n_PCO_cgto_orbs = 0
      if ( include_PCOs ) then ! Need to calculate the number of PCO orbs per irr
          do n=1,n_sym
              call atomic_orbital_basis % get_symmetry_flags(n,function_sym,i+1,j) ! Mark those AOs with IRR=n which are between i+1 and j, which represent PCOs functions.
              call orbital_data(n) % add_cms_continuum(function_sym) ! Now add the PCO functions for this symmetry to the orbital data set
              n_PCO_cgto_orbs(n) = count( function_sym )           ! For inclusion of PCOs count the number of orbitals in PCO index range and with IRR=n.
          end do
      end if
      n_cont_cgto_orbs = 0
      n_cont_bto_orbs  = 0
      do n=1,n_sym
          call atomic_orbital_basis % get_symmetry_flags(n,function_sym,j+1,k) ! Mark those AOs with IRR=ni between j+1 and k, which represent continuum functions.
          call orbital_data(n) % add_cms_continuum(function_sym) ! Add the initial set of continuum functions for this symmetry to the orbital data set
          call molecular_orbital_basis % add_shell(orbital_data(n))
          n_cont_cgto_orbs(n) = count(function_sym(atomic_orbital_basis % n_target_fns &
                                      : atomic_orbital_basis % n_target_fns + n_cont_cgto))
          n_cont_bto_orbs(n)  = count(function_sym(atomic_orbital_basis % number_of_functions - n_cont_bto &
                                      : atomic_orbital_basis % number_of_functions))
      end do

      ! Always print the orbital table from SCATCI_INTEGRALS
      call molecular_orbital_basis%print_energy_sorted_orbital_table(level1)
      call setup_integral_options
      call calculate_1el_integrals

      !
      ! ORBITAL ORTHOGONALIZATION:
      !

      call atomic_orbital_basis%assemble_overlap_matrix(ao_1el_integrals, overlap_matrix)

      do i = 1, n_sym

         !GS orthogonalize the target orbitals among themselves using the AO overlap integrals calculated before.
         if (nob_ortho(i) > 0) then
            write(stdout,'("0.0 Starting the orthogonalisation for symmetry ", I5)') i !Added a list of check points
            call molecular_orbital_basis % orthogonalize(gramm_schmidt = .true., overlap_matrix = overlap_matrix, symmetry = i, &
                 active_start = 1, active_end = nob_ortho(i), check_overlaps = .true.)
         end if

         !Pseudo-orbital initial orthogonalisation
         k = nob_ortho(i)+n_PCO_cgto_orbs(i)
         if (n_PCO_cgto_orbs(i) .gt. 0 .and. nob_ortho(i) .gt. 0 ) then
            write(stdout,'("1.0 Inital one-by-one orthogonalisation of PCOs ", I5)') i !Checkpoint
            call molecular_orbital_basis%orthogonalize(gramm_schmidt_one_by_one=.true.,overlap_matrix=overlap_matrix,symmetry=i,&
                 active_start =nob_ortho(i)+1,active_end=k,&
                 passive_start=1,             passive_end=nob_ortho(i),check_overlaps=.true.)
         end if

         !Symmetrically orthogonalise the pseudo orbitals among themselves.
         if (n_PCO_cgto_orbs(i) > 0) then
            write(stdout,'("1.1 Symmetric orthogonalisation of PCOs ", I5)') i !Checkpoint
            if (allocated(sym_ortho(i) % to_delete)) sym_ortho(i) % to_delete = .false.
            sym_ortho(i)%del_thrs = PCO_del_thrs
            call molecular_orbital_basis % orthogonalize(&
                 symmetric = .true., sym_ortho_data = sym_ortho(i), overlap_matrix = overlap_matrix, symmetry = i, &
                 active_start = nob_ortho(i) + 1, active_end = k, check_overlaps = .true.)
            call molecular_orbital_basis%delete_orbitals(i,sym_ortho(i)%to_delete)
            n_PCO_cgto_orbs(i)=n_PCO_cgto_orbs(i)-count(sym_ortho(i)%to_delete)
            k = nob_ortho(i)+n_PCO_cgto_orbs(i)
         end if

         !GS orthogonalize the continuum orbitals one-by-one against the target orbitals + PCOs
         j=molecular_orbital_basis%get_number_of_orbitals(i)
         if (j > k .and. k > 0) then
            write(stdout,'("2.0 One-by-one orthogonalisation of continuum ", I5)') i !Checkpoint
            call molecular_orbital_basis%orthogonalize(gramm_schmidt_one_by_one=.true.,overlap_matrix=overlap_matrix,symmetry=i,&
                 active_start =k+1,active_end=j, &
                 passive_start=1,  passive_end=k,check_overlaps=.true.)
         end if

         !If required orthogonalize the CGTO and BTO continuum separately first.
         if (preorthogonalize_continuum) then
            !BTOs: symmetric orthogonalization
            if (n_cont_bto_orbs(i) > 0) then
               write(stdout,'("2.1 Symmetric orthogonalisation of BTO continuum ", I5)') i !Checkpoint
               if (allocated(sym_ortho(i) % to_delete)) then
                  sym_ortho(i) % to_delete = .false.
               end if
               sym_ortho(i)%del_thrs = del_thrs
               call molecular_orbital_basis % orthogonalize(&
                    symmetric = .true., sym_ortho_data = sym_ortho(i), overlap_matrix = overlap_matrix, symmetry = i, &
                    active_start = k + n_cont_cgto_orbs(i) + 1, active_end = j, check_overlaps = .true.)
               call molecular_orbital_basis%delete_orbitals(i,sym_ortho(i)%to_delete)
            end if
            !CGTOs: symmetric orthogonalization
            if (n_cont_cgto_orbs(i) > 0) then
               write(stdout,'("2.2 Symmetric orthogonalisation of CGTO continuum ", I5)') i !Checkpoint
               if (allocated(sym_ortho(i) % to_delete)) then
                  sym_ortho(i) % to_delete = .false.
               end if
               sym_ortho(i)%del_thrs = del_thrs
               call molecular_orbital_basis % orthogonalize(&
                    symmetric = .true., sym_ortho_data = sym_ortho(i), overlap_matrix = overlap_matrix, symmetry = i, &
                    active_start = k + 1, active_end = k + n_cont_cgto_orbs(i), check_overlaps = .true.)
               call molecular_orbital_basis%delete_orbitals(i,sym_ortho(i)%to_delete)
            end if
         end if

         !ALL: symmetrically orthogonalize the BTO and CGTO continuum orbitals among themselves
         j = molecular_orbital_basis%get_number_of_orbitals(i)
         if (j > k) then
            write (stdout, '("2.3 Symmetric orthogonalisation of total continuum ", I5)') i ! DDL checkpoint
            if (allocated(sym_ortho(i) % to_delete)) then
               sym_ortho(i) % to_delete = .false.
            end if
            sym_ortho(i)%del_thrs = del_thrs
            call molecular_orbital_basis % orthogonalize(&
                 symmetric = .true., sym_ortho_data = sym_ortho(i), overlap_matrix = overlap_matrix, symmetry = i, &
                 active_start = k + 1, active_end = j, check_overlaps = .true.)
            call molecular_orbital_basis%delete_orbitals(i,sym_ortho(i)%to_delete)
         end if

         if (ortho_against_all_tgt_orbs .and. (ntot(i) > nob(i))) then
            if (allocated(to_delete)) deallocate(to_delete)
            j = molecular_orbital_basis%get_number_of_orbitals(i)
            allocate(to_delete(j),stat=err)
            if (err /= 0) call mpi_xermsg('main','main','error allocating to_delete array', err, 1)

            !Delete all target orbitals that we don't want to include in the integral transformation
            to_delete = .false.
            to_delete(nob(i)+1:ntot(i)) = .true.
            call molecular_orbital_basis%delete_orbitals(i, to_delete)
         end if

         ! DDL Orbital number write out including PCOs, for scripts.
         n_tot_orbs(i) = molecular_orbital_basis%get_number_of_orbitals(i)
         write (stdout, '("3.0 Completed orthogonalisation of symmetry ", I5)') i
         write (stdout, '("3.1 Orbitals for symmetry", I3, ", for target, PCO, continuum", 3I5)') &
                       i, nob(i), n_PCO_cgto_orbs(i),n_tot_orbs(i)-nob(i)-n_PCO_cgto_orbs(i)

      end do !i

      ! DDL Orbital number write out including PCOs, for scripts.
      write (stdout, '("Final number of orbitals for target    ",/, 8I5)') nob(:)
      write (stdout, '("Final number of orbitals for PCOs      ",/, 8I5)') n_PCO_cgto_orbs(:)
      write (stdout, '("Final number of orbitals for TGT+PCOs  ",/, 8I5)') n_PCO_cgto_orbs(:) + nob(:)
      write (stdout, '("Final number of orbitals for continuum ",/, 8I5)') n_tot_orbs(:) - nob(:) - n_PCO_cgto_orbs(:)
      write (stdout, '("Final number of orbitals for total     ",/, 8I5)') n_tot_orbs(:)

      call molecular_orbital_basis%delete_small_coefficients
      call molecular_orbital_basis%print

      ! Small change added by BC to output orbital data needed by QEC in a useful file.
      if (qmoln) then
          call molecular_orbital_basis%write_qec_orbital_table
      end if

      call molecular_orbital_basis%print_orbitals

   end subroutine construct_basis_sets


   subroutine setup_integral_options

      !If 1-electron Hamiltionian integrals are required the nuclear attraction ints. must be selected too.
      if (calculate_one_el_hamiltonian_ints) calculate_nuclear_attraction_ints = .true.

      integral_options%a = a !-1.0_cfp
      integral_options%max_ijrs_size = max_ijrs_size
      integral_options%calculate_overlap_ints = calculate_overlap_ints
      integral_options%calculate_kinetic_energy_ints = calculate_kinetic_energy_ints
      integral_options%calculate_property_ints = calculate_property_ints
      integral_options%max_property_l = max(2,max_property_l)
      integral_options%calculate_bbb_ints = calculate_bbb_ints
      integral_options%calculate_nuclear_attraction_ints = calculate_nuclear_attraction_ints
      integral_options%calculate_one_el_hamiltonian_ints = calculate_one_el_hamiltonian_ints
      integral_options%use_spherical_cgto_alg = use_spherical_cgto_alg
      integral_options%mixed_ints_method = mixed_ints_method
      integral_options%max_l_legendre_1el = max_l_legendre_1el
      integral_options%max_l_legendre_2el = max_l_legendre_2el
      integral_options%scratch_directory = scratch_directory
      integral_options%delta_r1 = delta_r1
      !integral_options%tol = 0.0_cfp !1e-8
      integral_options%print_integrals = print_1el_ints
      integral_options%dipole_damp_factor = dipole_damp_factor
      integral_options%two_p_continuum = two_p_continuum
      integral_options%keep_ao_integrals = keep_ao_integrals
      integral_options%HF_convergence = HF_convergence
      integral_options%HF_max_iters = HF_max_iters
      integral_options%only_construct_fock_blocks = only_construct_fock_blocks
      integral_options%canonize_virtuals_instead_of_continuum = canonize_virtuals_instead_of_continuum
      integral_options%use_ylm_scratch = use_ylm_scratch
      integral_options%use_aoints_scratch = use_aoints_scratch
      integral_options%use_transformation_scratch = use_transformation_scratch
      integral_options%aoints_index_scheme = aoints_index_scheme

   end subroutine setup_integral_options


   subroutine calculate_1el_integrals

      !describe where the AO integrals will be stored
      err = atomic_1el_integral_storage % init(memory=ao_1el_integrals) !(disk='./oneel')
      if (err /= 0) then
         call mpi_xermsg('main', 'main', 'error initializing the target atomic_integral_storage', err, 1)
      end if

      ! perform the integration
      call atomic_orbital_basis % one_electron_integrals(atomic_1el_integral_storage, integral_options)

   end subroutine calculate_1el_integrals


   subroutine calculate_2el_integrals

      call molecular_integral_storage%final

      integral_options%calculate_two_el_ints = .true.
      integral_options%print_integrals = print_2el_ints

      !describe where the AO integrals will be stored 
      err = atomic_2el_integral_storage%init(memory=ao_2el_integrals) !(disk='./oneel')
      if (err /= 0) then
         call mpi_xermsg('main', 'main', 'error initializing the target atomic_integral_storage', err, 1)
      endif

      !
      ! CALCULATE THE ATOMIC 2-ELECTRON INTEGRALS (unless Poisson equation method was chosen):
      !

      if (two_particle_integrals_method == 1) then
         call atomic_orbital_basis%two_electron_integrals(atomic_2el_integral_storage,integral_options)
      end if

      !
      ! TRANSFORM THE 2-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE MOLECULAR ORBITALS:
      !

      !describe where the MO integrals will be stored
      err = molecular_integral_storage%init(disk=mo_integrals_file_name) !(memory=mo_integrals)
      if (err /= 0) then
         call mpi_xermsg('main','main','error initializing the target molecular_integral_storage', err, 1)
      endif

      molecular_orbital_basis%ao_integral_storage => atomic_2el_integral_storage !point to the storage for the atomic integrals
      if (two_particle_integrals_method /= 1) then
         call molecular_orbital_basis % two_electron_integrals_poisson(molecular_integral_storage, integral_options)
      else if (molecular_2el_algorithm == 1 .or. (molecular_2el_algorithm == 0 .and. n_cont_bto > 0 .and. n_cont_cgto == 0)) then
         call molecular_orbital_basis % two_electron_integrals_sparse(molecular_integral_storage, integral_options)
      else
         call molecular_orbital_basis % two_electron_integrals(molecular_integral_storage, integral_options)
      end if

   end subroutine calculate_2el_integrals


   !> \brief  Terminate program if file exists
   !> \author J Benda
   !> \date   2020
   !>
   !> Stops the program if the provided (relative or absolute) path refers to an already existing file.
   !>
   subroutine abort_if_file_exists (filename)

      character(len=*), intent(in) :: filename
      logical :: fileexists

      fileexists = .false.
      inquire (file = filename, exist = fileexists)
      if (fileexists) then
         call mpi_xermsg('scatci_integrals', 'abort_if_file_exists', 'File "' // filename // '" already exists.', 1, 1)
      end if

   end subroutine abort_if_file_exists

   subroutine restricted_hartree_fock(geometry,n_sym,CGTO_shell_data,orbital_data)
      use symmetry_gbl, only : geometry_obj
      use basis_data_generic_gbl , only : orbital_data_obj
      use atomic_basis_gbl, only : CGTO_shell_data_obj
      use mpi_gbl, only: nprocs

      implicit none
      type(geometry_obj), intent(in) :: geometry
      integer, intent(in) :: n_sym
      type(CGTO_shell_data_obj), intent(in) :: CGTO_shell_data(:)
      type(orbital_data_obj), allocatable :: orbital_data(:)

      integer :: i, number_of_target_shells
      type(atomic_orbital_basis_obj), target :: atomic_orbital_basis_loc
      type(molecular_orbital_basis_obj) :: molecular_orbital_basis_loc
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data_loc(:)

      if (nprocs > 1) then
         call mpi_xermsg('scatci_integrals', 'restricted_hartree_fock', 'Not implemented for MPI mode.', 1, 1)
      endif

      number_of_target_shells = size(CGTO_shell_data)
      allocate(CGTO_shell_data_loc(number_of_target_shells))
      CGTO_shell_data_loc = CGTO_shell_data

      err = atomic_orbital_basis_loc%init(number_of_target_shells,geometry)
      if (err .ne. 0) stop "error initializing atomic orbital basis"

      do i=1,number_of_target_shells
         call atomic_orbital_basis_loc % add_shell(CGTO_shell_data_loc(i))
      enddo

      molecular_orbital_basis_loc%ao_basis => atomic_orbital_basis_loc
      err = molecular_orbital_basis_loc%init(n_sym,geometry)
      if (err .ne. 0) stop "error initializing molecular orbital basis"

      do i=1,n_sym
         call molecular_orbital_basis_loc % add_shell(orbital_data(i))
      enddo

      call setup_integral_options
      call molecular_orbital_basis_loc % solve_roothan_equations(integral_options,orbital_data)

      err = molecular_orbital_basis_loc % final()
      if (err .ne. 0) stop "error finalizing molecular orbital basis"

      deallocate(CGTO_shell_data_loc)

   end subroutine restricted_hartree_fock

end program scatci_integrals
