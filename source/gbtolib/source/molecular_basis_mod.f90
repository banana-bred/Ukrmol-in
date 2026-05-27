! Copyright 2019
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
module molecular_basis_gbl
   use file_mapping_gbl, only: file_mapping, map
   use precisn_gbl
   use const_gbl, only: line_len, no_header, stdout, level1, level2, level3
   use symmetry_gbl
   use integral_storage_gbl
   use mpi_gbl
   use atomic_basis_gbl
   use basis_data_generic_gbl
   use mpi_memory_gbl
 
   implicit none
 
   private
 
   !Only the objects themselves are visible from outside of this module.
   public molecular_orbital_basis_obj, orbital_data_obj, sym_ortho_io

   type, extends(basis_data_generic_obj) :: molecular_orbital_basis_obj
       !> Pointer to the AO basis object in terms of which this MO set is specified. This pointer must be associated before calling init.
       class(atomic_orbital_basis_obj), pointer :: ao_basis => null()
       !> This pointer must be associated by the user before calling one_electron_integrals or two_electron_integrals to the set of atomic integrals to be transformed.
       class(integral_storage_obj), pointer :: ao_integral_storage => null()
       !> Orbital data sets (typically one for each IRR). This array will have size number_of_shells but only during initialization, i.e. while the orbital sets are being added.
       !> This array is discarded following initialization.
       type(orbital_data_obj), allocatable, private :: orbital_data(:)
       !> Array which maps the absolute index of the orbital in the basis to the index of the orbital within the orbital_data array (absolute_to_relative).
       !> relative_to_absolute is the inverse of absolute_to_relative.
       !> absolute_to_relative(2,absolute_index) = IRR of the orbital with absolute index absolute_index.
       !> absolute_to_relative(1,absolute_index) = index of the orbital within the set of orbitals with the same IRR.
       integer, allocatable, private :: absolute_to_relative(:,:), relative_to_absolute(:,:)
       !> Point group symmetry of the orbitals in this set.
       integer :: pg = -1
       !> Number of irreducible representations.
       integer :: no_irr = 0
       !> so2mo_range: indices of the first and the last MO to which a given AO contributes. mo2so_range is the opposite for the MOs.
       !> These arrays have a meaning identical to the one explained in:
       !> S.Yamamoto, U. Nagashima, CPC 166 (2005) 58-65.
       integer, allocatable :: so2mo_range(:,:), mo2so_range(:,:)
       !> Arrays used to compute the index for the 2-particle symmetric (AB|CD) integrals. For details see add_function.
       integer, allocatable :: block_offset(:), sym_offset(:)
       !> For each orbital in the basis this list contains the value .true./.false. depending on whether the orbital is continuum or not.
       logical, private, allocatable :: is_continuum(:)
       !> Indices of the 2-electron integrals for each type of integral. Determined by two_electron_integrals.
       !> The pointer attribute is necessary for shared memory usage.
       integer, pointer :: ijkl_indices(:,:) => null()
       !> Used in case of shared memory allocation of ijkl_indices. It stores the MPI window that ijkl_indices is stored in.
       !> If it is -1 then we are in the non-shared mode for ijkl_indices.
       integer, private :: shared_window_ijkl_indices  = -1
       !> Index of the last integral stored in ijkl_indices.
       integer :: ind_ijkl_integral = 0
       !> Used only during initialization while the orbital sets are being added.
       logical, private, allocatable :: sets_added(:)
       !> Set to .true. following input of all basis functions for which space has been allocated.
       logical, private :: initialized = .false.
       !> Set to .true. following a call to init.
       logical, private :: init_called = .false.
       !> Auxiliary helper structure for virtual memory mapping.
       type(file_mapping) :: mmap
   contains
       !> Allocates space for a given number and type of shells of basis functions.
       procedure :: init
       !> Initialize the object by reading everything from the given integrals file and read all MO integrals.
       !> The input is the path to the moints file and the atomic basis object to be pointed at.
       !> On output the molecular basis is initialized, the integral arrays and the integral options objects
       !> are initialized with the data read-in from the moints file. It is assumed that the moints file includes
       !> 1-electron and 2-electron integrals.
       procedure :: read_basis_and_integrals
       !> Finalizes the basis set.
       procedure :: final
       !> Adds data for one shell into the basis set.
       procedure :: add_shell
       !> Prints the orbital basis set data to the stdout unit.
       procedure :: print => print_molecular_orbital_basis_obj
       !> Prints the orbital basis set data to a file. Uses by QEC controlled
       !> with a logical in the scatci_integrals program
       procedure :: write_qec_orbital_table
       !> Prints the orbital basis set data to any unit number
       procedure :: print_molecular_orbital_info
       !> Prints the orbital coefficients to the stdout unit.
       procedure :: print_orbitals
       !> Prints an orbital table showing the orbitals sorted in energy.
       procedure :: print_energy_sorted_orbital_table
       !> Calculates the values in arrays so2mo, mo2so, absolute_to_relative, relative_to_absolute, block_offset, sym_offset, is_continuum.
       procedure, private :: determine_auxiliary_indices
       !> Performs Gramm-Schmidt or Symmetric orthogonalization. The type of orthogonalization is selected using the logical
       !> parameters 'gramm_schmidt' or 'symmetric'. By default for G-S orthogonalization all orbitals are orthogonalized
       !> starting from the orbital with index 1 in the basis set. By default for symmetric orthogonalization all orbitals
       !> are orthogonalized. The range of 'active' and 'passive' orbitals can be selected by specifying the optional integers
       !> active_start, active_end and passive_start, passive_end which specify the range of indices for the orbitals to orthogonalize
       !> and not to orthogonalize. The format for the indices is (/num,sym/) where num is the (external) number of the orbital
       !> and sym is its symmetry. For symmetric orthogonalization the set of 'passive' orbitals is used only to check at the end
       !> that all 'active' orbitals are orthogonal to the 'passive' orbitals. The orthogonalization requires the AO overlap integrals array. 
       !> If symmetric orthogonalization is required then the data structure of type sym_ortho_io must be also present on input.
       !> It specifies the deletion thresholds for each symmetry. On output it contains the list of orbitals to delete
       !> (i.e. those that didn't pass the deletion threshold criterion) in the form of a logical array which marks the orbitals
       !> for deletion.
       procedure :: orthogonalize
       !> Deletes specified orbitals of a given symmetry. The orbitals to delete are marked in an input logical array which must
       !> have size equal to the number of orbitals in the given symmetry.
       procedure :: delete_orbitals
       !> Calculates and stores 1-electron integrals for all pairs of shells in the basis. The atomic integrals to be transformed
       !> are input via the type-bound pointer ao_integral_storage.
       procedure :: one_electron_integrals
       !> Calculates and stores 2-electron integrals for all pairs of shells in the basis. The atomic integrals to be transformed
       !> are input via the type-bound pointer ao_integral_storage.
       procedure :: two_electron_integrals
       !> Variant of "two_electron_integrals" useful for B-spline-only continuum basis. The subroutine makes use of sparsity of
       !> the two-electron integral matrix.
       procedure :: two_electron_integrals_sparse
       !> Transforms two of four atomic orbital indices to molecular ones.
       procedure :: transform_two_indices
       !> Retrieves a block of atomic integrals from the integral array.
       procedure :: fetch_atomic_integrals_block
       !> Moves temporary integral arrays to this object's data arrays and stores the integrals to disk.
       procedure :: finalize_two_electron_integrals_sparse
       !> Calculates indices for 1- or 2-electron integrals given their type and the number of basis function pairs/quartets. The
       !> two_p_continuum input variable is not used for MO integrals, the indexing method is the same whether we calculate the two
       !> particle continuum integrals or not.
       procedure :: integral_index
       !> Returns the name of the basis set.
       procedure :: get_basis_name
       !> Returns the name of the i-th shell in the basis set.
       procedure :: get_shell_name
       !> Returns the shell data for the i-th shell in the basis set.
       procedure :: get_shell_data
       !> Returns an array containing data for all shells of CGTOs in the basis.
       procedure :: get_all_orbital_sets
       !> Returns the matrix containing coefficients for all orbitals in the basis.
       procedure :: get_orbital_coefficient_matrix
       !> Returns the value of initialized.
       procedure :: is_initialized
       !> Same as for atomic_orbital_basis_obj but here we use the symmetry information for the orbitals.
       procedure :: get_continuum_flags
       !> Returns the number of orbitals in a given symmetry.
       procedure :: get_number_of_orbitals
       !> Returns the index within its own symmetry of a given orbital.
       procedure :: get_index_within_symmetry
       !> Returns the symmetry of a given orbital.
       procedure :: get_orbital_symmetry
       !> Given the pair of numbers num,sym it returns the absolute index of the orbital within the whole orbital set.
       procedure :: get_absolute_index
       !> Calculates orbital amplitudes for all continuum channels.
       procedure :: calculate_amplitudes
       !> Calculates radial charge densities of all orbitals in the basis. If the input value of rmat_radius is > 0.0_cfp then
       !> the continuum functions will be normalized to the R-matrix radius rmat_radius. If rmat_radius .le. 0 then normalization
       !> of the continuum functions will not be done. This is useful in case the orbital set corresponds to the Dyson orbitals obtained from CDENPROP: in this
       !> case the orbital coefficients already include the continuum normalization factors and therefore rmat_radius must be set to < 0.0_cfp.
       procedure :: radial_charge_density => orbital_radial_charge_density
       !> Deletes orbital coefficients with magnitude smaller than thrs_orb_cf.
       procedure :: delete_small_coefficients
       !> Writes to disk the array ijkl_indices and the value ind_ijkl_integral. As usual, only master writes its own array.
       procedure :: write_ijkl_indices
       !> Reads from the disk the array ijkl_indices and the value ind_ijkl_integral. The reading is done by the master task which also perform redistribution to other processes.
       !> If shared-memory MPI is used then each NODE keeps only one copy of the array this%ijkl_indices, otherwise the array is kept by every MPI task.
       procedure :: read_ijkl_indices
       !> Evaluates a given orbital (specified by its absolute index) at a set of points in space. The sign of the orbital
       !> at the corresponding points is output too in a separate array. No normalization of the continuum functions is performed
       !> so make sure you're either using atomic basis set whose functions have been normalized to the required R-matrix radius or 
       !> use orbital coefficients which include the continuum normalization factors (as is the case for the Dyson orbitals produced by CDENPROP).
       procedure :: eval_orbital
       !> Constructs canonical continuum
       procedure :: construct_canonical_continuum
       !> Construct a block of Fock matrix specified by absolute indices of two molecular orbital basis functions
       procedure, private :: construct_fock_matrix_block
       !> Calculates overlap integrals for a given range of molecular orbitals. It reports the results to standard output.
       !> On input the routine requires the storage for the 1-electron atomic integrals.
       procedure, private :: check_overlaps
       !> Obtains atomic and molecular basis and solves iteratively the roothan equations
       procedure :: solve_roothan_equations
       !> Calculates and stores 2-electron integrals for all pairs of shells in the basis using the Poisson equation.
       procedure :: two_electron_integrals_poisson
       !> Last stage of the Poisson equation method - populates the integral array and writes it to disk.
       procedure :: finalize_two_electron_integrals_poisson
       !> Internal subrotine used in two_electron_integrals_poisson
       !> Basically the same as free_scattering in one symmetry
       !> Used to store R-matrix amplitudes and energies for free Hamiltonian (Poisson eq.) into type diag_ham_obj
       procedure, private :: kinetic_energy_eigenstates

end type molecular_orbital_basis_obj

   !> \class <diag_ham_obj>
   !> The results of kintetic_energy_eigenstates are stored in this type
   type diag_ham_obj
      !> Number of channels in the symmetry
      integer :: nchan
      !> Number of functions
      integer :: n
      !> Ang mom l, azim mom m of the symmetry
      integer, allocatable :: chan_lm(:,:)
      !> On output, diagonalized hamiltonian is stored here
      real(kind=cfp), allocatable :: ham(:,:)
      !> On output, the R-matrix amplitudes are stored in w and energies in e
      real(kind=cfp), allocatable :: w(:,:), e(:)
   end type diag_ham_obj

   !> \class <sym_ortho_io>
   !> This data structure is used for input/output by the method orthogonalize of the molecular_orbital_basis_obj.
   type sym_ortho_io
      !> On input: deletion threshold for linearly dependent orbitals.
      real(kind=cfp) :: del_thrs(1:8) = (/-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp/)
      !> On output: number of functions 'deleted' in each symmetry.
      logical, allocatable :: to_delete(:)
   contains
      !> Checks that the input data is OK. We check: del_thrs.
      procedure :: check => check_sym_ortho_io
   end type sym_ortho_io

 contains
 
   function init(this,n,geometry)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: n
      class(geometry_obj), intent(in) :: geometry
      integer :: init
 
         write(level3,'("--------->","molecular_orbital_basis_obj:init")')
 
         init = 0
 
         if (this%initialized) then
            init = this%final()
            if (init .ne. 0) then
               call xermsg ('molecular_orbital_basis_obj', 'init', &
                            'Finalization has failed. See molecular_orbital_basis_obj%final for details.', init, 0)
               return
            endif
         endif
 
         if (n < 0) then
            init = 1
            call xermsg('molecular_orbital_basis_obj','init','On input the value of n was out of range.',init,0)
            return
         endif
 
         init = this%symmetry_data%init(geometry)
         if (init .ne. 0) then
            call xermsg ('molecular_orbital_basis_obj', 'init', &
                         'Symmetry initialization failed. See symmetry_obj%init for details.', init, 0)
            return
         endif

         this%pg = this%symmetry_data%get_pg()
         this%no_irr = this%symmetry_data%get_no_irrep(this%pg)

         !Note that in this case we don't try to exit the routine with an error
         !code. On run-time unassociated pointer can cause all kinds of mess and
         !this may make it hard to track the source of the problem if the error
         !codes from this routine are not properly processed.
         if (.not. associated(this % ao_basis)) then
            call xermsg ('molecular_orbital_basis_obj', 'init', &
                         'The ao_basis pointer to the AO basis set has not been associated. Fatal error.', 1, 1)
         end if

         !todo test that the atomic and the molecular symmetry data are the same!!!

         if (n /= this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'init', &
                         'The number n on input must be equal to the number of symmetries.', 2, 1)
         end if

         allocate(this%orbital_data(this%no_irr),this%sets_added(this%no_irr),stat=init)
         if (init .ne. 0) then
            call xermsg('molecular_orbital_basis_obj','init','Memory allocation of this%orbital_data has failed.',init,0)
            return
         endif

         this%sets_added = .false.
         this%init_called = .true.
         this%number_of_functions = 0
         this%number_of_shells = 0

         write(level3,'("<---------","molecular_orbital_basis_obj:init")')
 
   end function init

   function read_basis_and_integrals(this, file_name, atomic_orbital_basis, &
                                     one_electron_integrals, two_electron_integrals, options)

      use const_gbl, only: ijkl_indices_header, one_p_sym_ints, two_p_sym_ints
      use integral_storage_gbl

      implicit none

      integer :: read_basis_and_integrals
      class(molecular_orbital_basis_obj) :: this
      character(len=line_len), intent(in) :: file_name
      type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
      type(p2d_array_obj), target :: one_electron_integrals, two_electron_integrals
      type(integral_options_obj) :: options

      integer :: err, lunit, pos
      logical, parameter :: tgt_is_local = .true.
      type(integral_storage_obj) :: storage
      type(integral_storage_obj), target :: tgt_storage
      type(data_header_obj) :: header

         write(level3,'("--------->","molecular_orbital_basis_obj:read_basis_and_integrals")')

         read_basis_and_integrals = 0

         if (this%is_initialized()) then
            err = this%final()
            if (err /= 0) call xermsg ('molecular_orbital_basis_obj', 'read_basis_and_integrals', &
                                       'Finalization failed.', err, 1)
         endif

         !READ-IN ALL BASIS SETS AND ORBITAL DATA

         call atomic_orbital_basis%read(file_name)
         this%ao_basis => atomic_orbital_basis
         call this%read(file_name)

         err = storage%init(disk=file_name)
         if (err /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'read_basis_and_integrals', &
                         'Error initializing the integrals storage object.', err, 1)
         end if

         !READ-IN THE AUXILIARY ARRAY FOR 2-EL INTEGRAL INDEXING

         !We need to get the full header corresponding to the ijkl_indices.
         err = storage % integral_file % get_header_containing (header, &
                                                                this% get_basis_name(), &
                                                                ijkl_indices_header)
         if (err == 0) then
            call mpi_mod_bcast(header%first_record, master)
            lunit = storage%integral_file%get_unit_no()
            call mpi_mod_bcast(lunit, master)
            call this%read_ijkl_indices(lunit, file_name, header%first_record, pos)
         end if

         !READ-IN THE 1-ELECTRON INTEGRALS INTO ONE_ELECTRON_INTEGRALS

         err = tgt_storage%init(memory=one_electron_integrals)
         tgt_storage%integral_file%identifier = storage%integral_file%identifier
         if (err /= 0) call xermsg('molecular_orbital_basis_obj', 'read_basis_and_integrals', &
                                   'Error initializing tgt_storage.',err,1)

         !we need to get the full header corresponding to the transformed 1-electron integrals.
         err = storage % integral_file % get_header_containing(header, &
                                                               this % get_basis_name(), &
                                                               one_p_sym_ints)
         if (err /= 0) call xermsg('molecular_orbital_basis_obj', 'read_basis_and_integrals', &
                                   'Error locating the transformed 1-electron integrals; &
                                   &see data_header_obj%get_header_containing for details.',err,1)

         call tgt_storage%read(storage,header%name,options,tgt_is_local)
         call tgt_storage%final

         write(level1,'(/,"One electron integrals read-in.")')

         !READ-IN THE 2P INTEGRALS
         err = tgt_storage%init(memory=two_electron_integrals)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj', 'read_basis_and_integrals', &
                                     'Error initializing tgt_storage.',err,1)

         !we need to get the full header corresponding to the transformed 2-electron integrals.
         err = storage % integral_file % get_header_containing(header, &
                                                               this % get_basis_name(), &
                                                               two_p_sym_ints)
         if (err /= 0) call xermsg('molecular_orbital_basis_obj', 'read_basis_and_integrals', &
                                   'Error locating the transformed 2-electron integrals; &
                                   &see data_header_obj%get_header_containing for details.',err,1)

         call tgt_storage%read(storage,header%name,options,tgt_is_local)
         call tgt_storage%final

         write(level2,'(/,"Two electron integrals read-in.")')

         call storage%final

         write(level3,'("<---------","molecular_orbital_basis_obj:read_basis_and_integrals")')

   end function read_basis_and_integrals
 
   function final(this)
      use iso_c_binding, only: c_associated
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer :: final

      integer :: err
 
         write(level3,'("--------->","molecular_orbital_basis_obj:final")')

         final = 0
 
         this%number_of_shells = 0
         this%number_of_functions = 0
         if (allocated(this%orbital_data)) deallocate(this%orbital_data)
         if (allocated(this%absolute_to_relative)) deallocate(this%absolute_to_relative)
         if (allocated(this%relative_to_absolute)) deallocate(this%relative_to_absolute)
         if (allocated(this%sets_added)) deallocate(this%sets_added)
         if (this%shared_window_ijkl_indices /= -1) then
            call mpi_memory_deallocate_integer_2dim(this%ijkl_indices,size(this%ijkl_indices),this%shared_window_ijkl_indices)
            this%shared_window_ijkl_indices = -1
         else if (map .and. c_associated(this % mmap % ptr)) then
            call this % mmap % finalize
         else
            if (associated(this%ijkl_indices)) deallocate(this%ijkl_indices)
         endif
         this%ijkl_indices => null()

         this%init_called = .false.

         this%pg = -1
         this%no_irr = 0

         deallocate(this%so2mo_range,this%mo2so_range,this%block_offset,this%sym_offset,this%is_continuum,stat=err)
         if (err .ne. 0) final = 1
        
         this%initialized = .false.
         !WE MUST NOT NULLIFY THE POINTER HERE: IF WE DO IT HERE THEN READING OF
         !THE BASIS WILL FAIL SINCE THE POINTER WILL NOT BE ASSOCIATED WHEN
         !ADDING SETS OF ORBITALS. (FINALIZATION OF THE BASIS IS ALWAYS DONE
         !BEFORE READING ANY BASIS).
         !nullify(this%ao_basis)
 
         write(level3,'("<---------","molecular_orbital_basis_obj:final")')
 
   end function final
 
   subroutine add_shell(this,shell_data)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(shell_data_obj), intent(inout) :: shell_data

      type(orbital_data_obj), allocatable :: temp_orbital_data(:)
      integer :: err, i

         write(level3,'("--------->","molecular_orbital_basis_obj:add_shell")')

         if (.not. this % init_called) then
            call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                         'Attempt to call add_shell before calling init.', 1, 1)
         end if
 
         if (this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                         'All orbital sets for which space has been allocated have already been supplied.', 2, 1)
         end if

         call shell_data%normalize
         call shell_data%print

         select type (orbital_set => shell_data)
            type is (orbital_data_obj)

               !No duplicities allowed and symmetries of all orbital sets must be consistent:
               if (this % pg /= orbital_set % point_group) then
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'All sets of orbitals in the basis must belong to the same point group symmetry.', 3, 1)
               end if

               if (orbital_set%irr <= 0 .or. orbital_set%irr > this % no_irr) then
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The IRR of the orbital set on input is out of range.', 4, 1)
               end if

               if (this % sets_added(orbital_set % irr)) then
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The basis set already contains orbitals for this symmetry.', 5, 1)
               end if

               if (orbital_set%number_of_coefficients .ne. this%ao_basis%number_of_functions) then
                  print *,orbital_set%number_of_coefficients,this%ao_basis%number_of_functions
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The number of orbital coefficients is not compatible with the AO basis set.', 6, 1)
               endif

               !Resize the basis of orbitals:
               call move_alloc(this%orbital_data,temp_orbital_data)
               this%number_of_shells = this%number_of_shells + 1
               allocate(this%orbital_data(this%number_of_shells),stat=err)
               if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'add_shell', 'Memory allocation 1 failed.',err, 1)

               !Copy the data for the previous sets of orbitals:
               do i=1,this%number_of_shells-1
                  this%orbital_data(i) = temp_orbital_data(i)
               enddo !i
               if (this%number_of_shells > 1) deallocate(temp_orbital_data)

               !Add the new orbital_set data to the list:
               this%orbital_data(this%number_of_shells) = orbital_set
               write(level3,'("Orbital set of type orbital_data_obj has been added to the molecular basis.")')
               this%sets_added(orbital_set%irr) = .true.
            class default
               call xermsg ('molecular_orbital_basis_obj', 'add_shell', 'The shell type must be orbital_data_obj.',1, 1)
         end select

         this%number_of_functions = this%number_of_functions + shell_data%number_of_functions
 
         if (this%number_of_shells .eq. this%no_irr) then

            write(level2,'(/,"Orbitals for all symmetries have been supplied. &
                             &Generating indices and analyzing orbital coefficients...")')

            call this%determine_auxiliary_indices
   
            this%initialized = .true.

         endif
 
         write(level3,'("<---------","molecular_orbital_basis_obj:add_shell")')
 
   end subroutine add_shell

   !> \brief   Subroutine to write the orbital information to file. This is
   !>          needed by QEC to easily get the number of continuum orbitals per
   !>          symmetry irrep. Writes the data to file named orbital.data.
   !> \authors Bridgette Cooper
   !> \date    October 2019
   subroutine write_qec_orbital_table(this)
      implicit none
      integer :: ierr, unit_number
      class(molecular_orbital_basis_obj) :: this

         unit_number = 42018
         open(unit=unit_number, file='orbital.data', status='replace', form='FORMATTED', &
              iostat=ierr)
         if (ierr .ne. 0) then
            call xermsg('molecular_basis_mod','write_qec_orbital_table', &
                        'Error opening the file to write orbital data to', &
                        ierr, 1)
         end if
         call this%print_molecular_orbital_info(unit_number)
         close(unit=unit_number) 

   end subroutine write_qec_orbital_table
   
   !> \brief   Subroutine to print the molecular orbital information to the
   !>          standard output.
   !> \authors Zdenek Masin. Modified by Bridgette Cooper
   !> \date    October 2019
   subroutine print_molecular_orbital_basis_obj(this)
      implicit none
      class(molecular_orbital_basis_obj) :: this

         call this%print_molecular_orbital_info(level3)

   end subroutine print_molecular_orbital_basis_obj

   !> \brief  General function for printing the molecular orbital basis
   !>         information to any unit number. Called by print_molecular_basis_obj
   !>         with stdout as the unit number.
   !> \author Zdenek Masin. Modified by Bridgette Cooper
   !> \date   October 2019
   subroutine print_molecular_orbital_info(this, unit_number)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer :: unit_number !> The unit number to write this information to

      integer :: i, n_tgt(this%no_irr), n_cnt(this%no_irr)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'print_molecular_orbital_info', &
                         'The object has not been initialized or not all orbitals have been read-in.', 1, 1)
         end if

         write(unit_number,'(/,"--------->","molecular_orbital_basis_obj:print",/)')

         write(unit_number,'("Point-group symmetry identifier: ",i0)') this%pg
         write(unit_number,'("Number of irreducible representations: ",i0)') this%no_irr
         write(unit_number,'("Number of molecular orbitals in each irreducible representation: ")') 
         write(unit_number,'(8(i0,1x))') this%orbital_data(1:this%no_irr)%number_of_functions

         n_tgt = 0
         n_cnt = 0
         do i=1,this%number_of_functions
            if (this%is_continuum(i)) then
               n_cnt(this%absolute_to_relative(2,i)) = n_cnt(this%absolute_to_relative(2,i)) + 1
            else
               n_tgt(this%absolute_to_relative(2,i)) = n_tgt(this%absolute_to_relative(2,i)) + 1
            endif
         enddo !i

         write(unit_number,'("Number of target orbitals: ")')
         write(unit_number,'(8(i0,1x))') n_tgt(1:this%no_irr)

         write(unit_number,'("Number of continuum orbitals: ")')
         write(unit_number,'(8(i0,1x))') n_cnt(1:this%no_irr)

         write(unit_number,'("Name of the associated AO basis: ",a)') trim(this%ao_basis%get_basis_name())
         write(unit_number,'("Number of AO basis functions for each irreducible representation: ")') 
         write(unit_number,'(8(i0,1x))') this%orbital_data(1:this%no_irr)%number_of_coefficients

         write(unit_number,'(/,"Symmetries and indices of the orbitals:")')

         write(unit_number,'("Index within symmetry, Orbital symmetry, Overall index, Is continuum")')
         do i=1,this%number_of_functions
            write(unit_number,'(3(i0,1x),1X,l1)') this%absolute_to_relative(1:2,i), i, this%is_continuum(i)
         enddo

         write(unit_number,'("<---------","done:molecular_orbital_basis_obj:print")')

   end subroutine print_molecular_orbital_info

   subroutine print_orbitals(this)
      use const_gbl, only: orbs_line
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: i, j, k, m, mx, symmetry
      real(kind=cfp) :: cf_tmp(orbs_line)
!      real(kind=cfp), allocatable :: cf(:,:)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'print_orbitals', &
                         'The object has not been initialized or not all orbitals have been read-in.', 1, 1)
         end if

         write(level3,'(/,"--------->","molecular_orbital_basis_obj:print_orbitals",/)')

         write(level3,'(/,"so2mo_range:")')
         mx = this%ao_basis%number_of_functions
         do i=1,mx
            write(level3,'("AO ",i0," MO start, end: ",i0,1x,i0)') i, this%so2mo_range(1:2,i)
         enddo

         write(level3,'(/,"mo2so_range:")')
         do i=1,this%number_of_functions
            write(level3,'("MO ",i0," AO start, end: ",i0,1x,i0)') i, this%mo2so_range(1:2,i)
         enddo

         write(level3,'(/,10X,"Orbital coefficients follow")')

!         call this%get_orbital_coefficient_matrix(cf)
!         write(stdout,'(/,10X,"Orbital coefficients follow")')
!         k = 0
!         do i=1,this%number_of_functions/orbs_line
!            write(stdout,'(/,10X,50(i,2X))') (j,j=k+1,k+orbs_line)
!            do j=1,mx
!               write(stdout,'(i,50e25.15)') j, cf(j,k+1:k+orbs_line)
!            enddo
!            k = k + orbs_line
!         enddo
!
!         m = mod(this%number_of_functions,orbs_line)
!         if (m > 0) then
!            write(stdout,'(/,10X,50(i,2X))') (j,j=k+1,k+m)
!            do j=1,mx
!               write(stdout,'(i,50e25.15)') j, cf(j,k+1:k+m)
!            enddo
!         endif

         do symmetry=1,this%no_irr
            write(level3,'(/,10X,"Symmetry: ",i4)') symmetry
            k = 0
            do i=1,this%orbital_data(symmetry)%number_of_functions/orbs_line
               write(level3,'(/,10X,50(i0,2X))') (this%relative_to_absolute(j,symmetry),j=k+1,k+orbs_line)
               do j=1,mx
                  cf_tmp(1:orbs_line) = this%orbital_data(symmetry)%coefficients(j,k+1:k+orbs_line)
                  write(level3,'(i0,50e25.15)') j, cf_tmp(1:orbs_line) !this%orbital_data(symmetry)%coefficients(j,k+1:k+orbs_line)
               enddo
               k = k + orbs_line
            enddo
   
            m = mod(this%orbital_data(symmetry)%number_of_functions,orbs_line)
            if (m > 0) then
               write(level3,'(/,10X,50(i0,2X))') (this%relative_to_absolute(j,symmetry),j=k+1,k+m)
               do j=1,mx
                  cf_tmp(1:m) = this%orbital_data(symmetry)%coefficients(j,k+1:k+m)
                  write(level3,'(i0,50e25.15)') j, cf_tmp(1:m) !this%orbital_data(symmetry)%coefficients(j,k+1:k+m)
               enddo
            endif
         enddo !symmetry

         write(level3,'("<---------","done:molecular_orbital_basis_obj:print_orbitals")')

   end subroutine print_orbitals

   subroutine print_energy_sorted_orbital_table(this, print_level_override)
      use common_obj_gbl, only: print_orbital_table
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in), optional :: print_level_override

      integer :: i, j, n, err, n_tgt, num, print_level
      integer, allocatable :: num_sym(:,:)
      real(kind=cfp), allocatable :: energies(:)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'print_energy_sorted_orbital_table', &
                         'The object has not been initialized or not all orbitals have been read-in.', 1, 1)
         end if

         write(level3,'(/,"--------->","molecular_orbital_basis_obj:print_energy_sorted_orbital_table",/)')

         if (present(print_level_override)) then
            print_level = print_level_override
         else
            print_level = level2
         endif

         write(print_level,'(/,10X,"Continuum orbitals will be ignored.")')

         n_tgt = this%number_of_functions - count(this%is_continuum)

         allocate(energies(n_tgt),num_sym(2,n_tgt),stat=err)
         if (err /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'print_energy_sorted_orbital_table', 'Memory allocation failed.', err, 1)
         end if

         n = 0
         do i=1,this%number_of_shells !over all symmetries
            do j=1,this%orbital_data(i)%number_of_functions
               num = this%relative_to_absolute(j,i)
               if (this%is_continuum(num)) cycle
               n = n + 1
               energies(n) = this%orbital_data(i)%energy(j)
               num_sym(1,n) = j
               num_sym(2,n) = i
            enddo !i
         enddo

         call print_orbital_table(energies,num_sym,n_tgt,this%number_of_shells,.true.,print_level)

         write(level3,'("<---------","done:molecular_orbital_basis_obj:print_energy_sorted_orbital_table")')

   end subroutine print_energy_sorted_orbital_table
 
   subroutine one_electron_integrals(this,integral_storage,integral_options)
      use const_gbl
      use omp_lib, only: omp_get_wtime
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(integral_options_obj), intent(in) :: integral_options
      class(integral_storage_obj), intent(inout) :: integral_storage
 
      !Input/output of the calculated integrals:
      type(integral_storage_obj) :: ao_integrals_disk
      type(integral_options_obj) :: ao_int_opt
      integer :: number_of_integrals, lunit, first_record, current_pos, last_record, d1, d2, no_blocks
      integer :: p, q, i, j, ij, err, no_ao, no_mo, int_type
      integer, allocatable :: int_index(:,:), ind(:)
      real(kind=cfp), allocatable :: cf(:,:), cf_t(:,:), iq(:,:), iq_t(:,:), ao_int(:)
      real(kind=cfp) :: mo_int
      real(kind=wp) :: t
      type(p2d_array_obj), target :: integral_src, integral_tgt !we really need two of these in case disk-to-disk AO-MO run is required
      type(p2d_array_obj), pointer :: ao_integrals, mo_integrals
      logical, parameter :: ao_is_local = .true.
      integer, parameter :: number_of_blocks = 0
      character(len=line_len), allocatable :: column_descriptor(:)
      character(len=line_len) :: ao_header, mo_header
 
         call mpi_mod_barrier(err)
         t = omp_get_wtime()
 
         write(level3,'("--------->","molecular_orbital_basis_obj:one_electron_integrals")')
 
         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', 'The basis set has not been initialized.', 1, 1)
         end if

         if (this % number_of_functions == 0) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                         'Number of molecular orbitals on input is zero.', 2, 1)
         end if

         if (.not.associated(this%ao_integral_storage) .or. .not.associated(this%ao_basis)) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                         'On input at least one of this%ao_integral_storage, this%ao_basis have not been associated.', 4, 1)
         endif

         !Header for the AO integrals that we're looking for
         ao_header = this%ao_integral_storage%contruct_header_string(this%ao_basis%get_basis_name(),one_electron_ints)

         !Header for the MO ntegrals that will be calculated
         mo_header = integral_storage%contruct_header_string(this%get_basis_name(),one_p_sym_ints)

         !In this section we associate ao_integrals which contains the input AO integrals with the appropriate source.
         !In case the AO integrals are stored in memory then we point directly to the array
         !holding them. If the AO integrals are on disk then we load them into the local array 'integral' and set the pointer ao_integrals to that.
         !At the moment 1p integral transform using SHARED input AO integrals is not supported.
         if (this%ao_integral_storage%in_memory()) then
            ao_integrals => this%ao_integral_storage%integral
            if (this%ao_integral_storage%data_header%name .ne. ao_header) then
               call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                            'The AO integrals on input are not compatible with the AO basis set for the MO basis set.', 5, 1)
            endif
            write(level2,'("AO integrals with header: ",a)') ao_header
            if (ao_integrals % have_offsets()) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'The AO integrals on input must be LOCAL, i.e. not shared accross processes.', 6, 1)
            end if
         endif

         !load the AO integrals into memory as LOCAL arrays (if they are stored on disk)
         if (this%ao_integral_storage%on_disk()) then
            write(level2,'("Loading AO integrals from the disk...")')

            err = ao_integrals_disk%init(memory=integral_src)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', 'Memory allocation 3 failed.', err, 1)
            end if

            call ao_integrals_disk%read(this%ao_integral_storage,ao_header,ao_int_opt,ao_is_local)
            ao_integrals => ao_integrals_disk%integral !this points to the local array 'integral_src'
            write(level2,'("AO integrals with header: ",a)') ao_header
            write(level2,'("...done")')
         endif

         !BEYOND THIS POINT ao_integrals POINTS TO AN ARRAY CONTAINING THE AO INTEGRALS TO BE TRANSFORMED

         call ao_integrals%get_array_dim(d1,d2,no_blocks) !This gives the dimensions of the AO integrals array
         call ao_integrals%get_column_descriptor(column_descriptor)

         write(level2,'("On input there is ",i0," types of AO integrals")') d2
         write(level2,'("Number of AO integrals of each type: ",i0)') d1

         !note that instead of loading the AO integrals we can calculate them now since we have the pointer to the AO basis set (this%ao_basis) and the AO integral routine...

         number_of_integrals = this%number_of_functions*(this%number_of_functions+1)/2 !total number of integrals to calculate; note that this does not include symmetry reduction that we can acheive if we know the symmetry of the 1p operator.

         !Allocate the output arrays if we request the output to be stored in memory.
         if (integral_storage%in_memory()) then
            integral_storage%data_header%name = mo_header

            mo_integrals => integral_storage%integral
            !We allocate space for a non-indexed (that is purely local) array with d2 columns and number_of_integrals rows.
            !The columns in the mo_integrals array correspond to the types of AO integrals we have on input.
            err = mo_integrals%init(number_of_integrals,d2,number_of_blocks,column_descriptor)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'Array initialization failed; see p2d_array_obj%init.', err, 1)
            end if
         endif

         !If we request the output to be stored on disk then start a new record on the data file that will contain the integrals.
         !We also allocate temporary storage for the transformed integrals.
         if (integral_storage%on_disk()) then

            !temporary storage for the integrals
            mo_integrals => integral_tgt
            !We allocate space for a non-indexed (that is purely local) array with d2 columns and number_of_integrals rows.
            !The columns in the mo_integrals array correspond to the types of AO integrals we have on input.
            err = mo_integrals%init(number_of_integrals,d2,number_of_blocks,column_descriptor)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'Array initialization 2 failed; see p2d_array_obj%init.', err, 1)
            end if

            lunit = integral_storage%integral_file%get_unit_no()             !unit that is associated to the file opened
            first_record = integral_storage%integral_file%start_record(mo_header) !position, within the data file, of the first record available for the integral data
         endif

         !BEYOND THIS POINT mo_integrals POINTS TO AN ARRAY CONTAINING THE TRANSFORMED MO INTEGRALS

         !1-PARTICLE INTEGRAL TRANSFORM STARTS HERE: the AO integrals are accessed through the pointer ao_integrals; the MO integrals are accessed through the pointer mo_integrals
!
!--------The orbitals may have changed so recalculate the so2mo,mo2so indices.
!
         call this%determine_auxiliary_indices

         no_ao = this%ao_basis%number_of_functions !total number of AOs
         no_mo = this%number_of_functions !total number of MOs

         allocate(cf_t(no_mo,no_ao),int_index(1:2,no_ao),ind(no_ao),iq(this%number_of_functions,no_ao), &
                  iq_t(no_ao,this%number_of_functions),ao_int(no_ao),stat=err)
         if (err /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', 'Memory allocation 4 failed.', err, 1)
         end if

         !Copy the orbital coefficients to one array: this relies on the fact that the molecular orbitals are indexed symmetry by symmetry.
         call this%get_orbital_coefficient_matrix(cf)

         !transpose the MO coefficient matrix: we use it in the first step where we iterate over the MOs.
         cf_t = transpose(cf)

         !iterate over all types of AO integrals (distribute types cyclically among parallel tasks)
         do int_type = myrank + 1, d2, nprocs

            write(level1,'("Transforming AO integral type: ",a," ...")') trim(adjustl(column_descriptor(int_type)))

            !iterate over all AO integrals (p|O|q) and accumulate their contributions to (i|O|q) where i is MO.
            iq(:,:) = 0.0_cfp
            do p=1,no_ao
               !construct the list of AO indices corresponding to all unique pairs (pq) of the AOs.
               do q=1,p
                  int_index(1,q) = p
                  int_index(2,q) = q
               enddo !q
   
               !calculate indices of the AO integrals corresponding to all unique pairs (pq) of the AOs.
               ind(1:p) = this%ao_basis%integral_index(column_descriptor(int_type),int_index(1:2,1:p),&
                                                       integral_options%two_p_continuum)

               !load the AO integrals (p|O|q) with GLOBAL indices 'ind' into ao_int.
               do j=1,p
                  ao_int(j) = ao_integrals%a(ind(j),int_type)
               enddo

               !transform the first index: p->i 
               do j=1,p
                  q = int_index(2,j)
                  if (ao_int(j) .ne. 0.0_cfp) then
                     if (p .eq. q) then
                        do i = this%so2mo_range(1,p), this%so2mo_range(2,p)
                           iq(i,q) = iq(i,q) + ao_int(j)*cf_t(i,p)
                        enddo !i
                     else !we assume the integral is symmetric: (p|O|q) = (q|O|p) so we calculate at once contributions of (p|O|q) to iq(i,q) and iq(i,p)
                        do i = this%so2mo_range(1,p), this%so2mo_range(2,p)
                           iq(i,q) = iq(i,q) + ao_int(j)*cf_t(i,p)
                        enddo !i
                        do i = this%so2mo_range(1,q), this%so2mo_range(2,q)
                           iq(i,p) = iq(i,p) + ao_int(j)*cf_t(i,q)
                        enddo !i
                     endif
                  endif
               enddo !j
            enddo !p

            iq_t = transpose(iq)

            !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,ij,mo_int,q) SHARED(this,iq_t,cf,integral_options,mo_integrals,int_type)
            !$OMP DO SCHEDULE(DYNAMIC)
            do i=1,this%number_of_functions
               do j=1,i
                  !the ij indexing is equivalent to the index function
                  ij = i*(i-1)/2+j

                  mo_int = sum(iq_t(this%mo2so_range(1,j):this%mo2so_range(2,j),i) &
                               * cf(this%mo2so_range(1,j):this%mo2so_range(2,j),j))

                  !The line above is the equivalent of:
                  !mo_int = 0.0_cfp
                  !do q=this%mo2so_range(1,j), this%mo2so_range(2,j)
                  !   mo_int = mo_int + iq_t(q,i)*cf(q,j)
                  !enddo !q

                  if (abs(mo_int) < integral_options%tol) then
                     mo_int = 0.0_cfp
                  else
                     mo_integrals%a(ij,int_type) = mo_int
                  endif

               enddo !j
            enddo !i
            !$OMP END DO
            !$OMP END PARALLEL

            write(level1,'("...done")')

         enddo !int_type

         deallocate(cf_t,int_index,iq,ind)

         nullify(ao_integrals)
         !FROM NOW ON ao_integrals is nullified

         ! synchronize integral types
         do int_type = 1, d2
            call mpi_reduce_inplace_sum_cfp(mo_integrals%a(:,int_type), number_of_integrals)
         end do

         !If requested print the non-zero integrals
         if (integral_options%print_integrals) then
            call mo_integrals%print(.true.)
         endif

         !dump all integrals to disk and close the record
         if (integral_storage%on_disk()) then

            write(level2,'("Saving integrals to disk...")')

            !The first record are the integral options.
            call integral_options%write(lunit,first_record,current_pos)

            ! The second record are the ordered integrals: we write them using a method whose choice
            ! depends on whether each process keeps the copy of the full integral array or not.
            i = master
            call mo_integrals%write(lunit,current_pos,last_record,i)
         
            !Every process closes the record so that they all keep identical header information.
            call integral_storage%integral_file%close_record(mo_header,first_record,last_record)
               
            err = mo_integrals%final()
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'one_electron_integrals', &
                             'Deallocation of the temporary integral array failed.', 5, 1)
            end if

            nullify(mo_integrals)
            !FROM NOW ON mo_integrals is nullified

            write(level2,'("...done")')

         else
 
            nullify(mo_integrals)
            !FROM NOW ON mo_integrals is nullified

         endif
 
         write(level1,'(a,f0.1)') 'One-electron integrals transformed in [s]: ', omp_get_wtime() - t
         write(level3,'("<---------","molecular_orbital_basis_obj:one_electron_integrals")')
 
         call mpi_mod_barrier(err)
 
   end subroutine one_electron_integrals


    !> \brief   Transform atomic 2-electron integrals to molecular ones
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Calculates the molecular-orbital 2-electron integrals (ij|kl) from the known atomic-orbital 2-electron integrals [pq|rs].
    !> The algorithm proceeds in several steps:
    !>  1. Expand all atomic integrals to a working array, including all those redundant in the sense of index symmetries.
    !>  2. Sort the integral array so that all integrals with common p,q are clustered together, resulting in a sequence of
    !>     sparse matrices whose elements are indexed by atomic indices r,s.
    !>  3. Transform both indices of these sparse matrices to molecular ones by application of the expansion coefficient matrices,
    !>     yielding object [pq|kl).
    !>  4. Sort the integral array again, but now group together elements with the same k,l, resulting in a sequence of sparse
    !>     matrices whose elements are indexed by atomic indices p,q. Symbolically: (kl|pq].
    !>  5. Transform both indices p,q by application of the coefficient matrices, yielding (kl|ij).
    !>  6. Finally, sort the integrals to satisfy expectations of the library.
    !>
    !> In the present implementation, steps 1, 2 and 3 are fused to avoid large memory use; only a few blocks are expanded
    !> at a time.
    !>
    !> Some work can be saved with the knowledge of index and orbital symmetries:
    !>  - In the step 1 discard CCCC and CCCT types if not needed.
    !>  - In the step 3 it is possible to calculate the numbers [pq|kl) just for k >= l, due to symmetry k <-> l. In addition
    !>    to this, the CCCC and CCCT combinations can be ignored if two electrons in continuum are not required.
    !>  - In the step 5 it is possible to used the same, ignoring anything else than i >= j, and one can skip also all i,j
    !>    pairs, whose combined orbital symmetry is different than the combined orbital symmetry of the pair k,l. This is
    !>    due to the fact that the product of the two pairs i,j and k,l must be totally symmetric. Furthermore, one can
    !>    skip all i,j pairs that violate [ij] >= [kl], due to the symmetry i,j <-> k,l. Here [ij] is the standard triangular
    !>    multi-index function. Again, as before, CCCC and CCCT combination can be skipped.
    !>
    !> Note that this subroutine uses rectangular indexing function and transforms the indices to the triangular ones only
    !> at the very end.
    !>
    subroutine two_electron_integrals_sparse(this, integral_storage, integral_options)

        use omp_lib,                  only: omp_get_wtime
        use sparse_block_storage_gbl, only: sparse_block_storage_obj

        class(molecular_orbital_basis_obj)         :: this
        class(integral_options_obj), intent(in)    :: integral_options
        class(integral_storage_obj), intent(inout) :: integral_storage

        type(p2d_array_obj),     pointer     :: ao_integrals
        type(sparse_block_storage_obj)       :: mo_integrals
        character(len=line_len), allocatable :: column_descriptor(:)
        character(len=:),        allocatable :: scratch_prefix

        real(kind=cfp), allocatable :: cf(:,:), cft(:,:), Cv(:)
        integer,        allocatable :: Cp(:), Cj(:)

        real(kind=cfp) :: tolerance, max_buffer_size
        real(kind=wp)  :: start_t, end_t
        integer        :: d1, d2, no_blocks_ao, no_ao, no_mo, Rn
        logical        :: skip2ec, ao_is_local

        write(level3,'("--------->","molecular_orbital_basis_obj:two_electron_integrals_sparse")')

        start_t = omp_get_wtime()

        ao_integrals => this % ao_integral_storage % integral

        ! use scratch if specified, but only when `use_transformation_scratch` was requested
        scratch_prefix = integral_options % scratch_directory
        if (.not. integral_options % use_transformation_scratch) scratch_prefix = ''

        ! figure out if the AO integrals stored in memory are distributed (multiple MPI tasks) or local (single MPI task)
        ao_is_local = .false.
        if (this % ao_integral_storage % in_memory()) then
           ao_is_local = .not.(ao_integrals % have_offsets())
        endif

        call ao_integrals % get_column_descriptor(column_descriptor)
        call ao_integrals % get_array_dim(d1, d2, no_blocks_ao)

        no_ao = this % ao_basis % number_of_functions  ! total number of AOs
        no_mo = this % number_of_functions             ! total number of MOs

        write(level1, '("On input there is ",I0," types of AO integrals")') d2
        write(level1, '("Total number of AO integrals of each type: ",I0)') d1
        !write(level2, '("Every process keeps the full set of AO integrals.")')
        write(level1, '("Number of atomic orbitals: ",I0)') no_ao
        write(level1, '("Number of molecular orbitals: ",I0)') no_mo

        ! construct coefficient matrix C(no_ao,no_mo)
        call this % get_orbital_coefficient_matrix(cf)
        call sparse_from_dense(Cp, Cj, Cv, cf)
        write(level1, '("Coefficient matrix density: ",I0," %")') (Cp(no_mo + 1) - 1) * 100 / (no_ao * no_mo)

        max_buffer_size = integral_options % max_ijrs_size

        ! prepare the block storage for semitransformed integrals
        call mo_integrals % reset(no_ao*no_ao, scratch_prefix)

        ! transform the two summation indices r,s to k,l
        tolerance = 0
        skip2ec = (this % ao_basis % n_cont_fns > 0) .and. (.not. integral_options % two_p_continuum)

        call this % transform_two_indices(no_ao, no_ao, no_mo, mo_integrals, Cp, Cj, Cv, &
                                          1, tolerance, skip2ec, max_buffer_size, &
                                          integral_options % aoints_index_scheme, &
                                          integral_options % keep_ao_integrals)

        ! invert indexing (ie. make p,q the summation indices)
        call invert_indexing(mo_integrals%Rn, mo_integrals%Ri, no_ao, no_mo)

        ! sort the integral array (cluster by k,l)
        call sort_intermediate_integrals(mo_integrals%Rn, mo_integrals%Ri, mo_integrals%Rv, no_ao * no_ao, &
                                         integral_options%max_ijrs_size, scratch_prefix)

        ! in the current (new) scheme all AO integrals will have been used by now
        ! so we can neglect those that are too small in the final step.
        tolerance = integral_options % tol

        ! prepare the block storage for fully transformed integrals
        call mo_integrals % reset(no_mo*no_mo, scratch_prefix)

        ! transform the two summation indices p,q to i,j
        call this % transform_two_indices(no_mo, no_ao, no_mo, mo_integrals, Cp, Cj, Cv, &
                                          2, tolerance, skip2ec, max_buffer_size, &
                                          integral_options % aoints_index_scheme, &
                                          integral_options % keep_ao_integrals)

        ! reindex the array to the expected triangular format
        call rect_index_to_tri_index(mo_integrals%Rn, mo_integrals%Ri, no_mo, 1)

        ! sort the integral array in accord with the triangular indexing
        call sort_final_integrals(mo_integrals%Rn, mo_integrals%Ri, mo_integrals%Rv)

        ! store integrals to disk and finalize the calculation process
        call this % finalize_two_electron_integrals_sparse(mo_integrals%Rn, mo_integrals%Ri, mo_integrals%Rv, &
                                                           integral_storage, integral_options, column_descriptor)

        end_t = omp_get_wtime()

        write(level1,'("two_electron_integrals_sparse took [s]: ",F25.15)') end_t - start_t
        write(level3,'("<---------","molecular_orbital_basis_obj:two_electron_integrals_sparse")')

    end subroutine two_electron_integrals_sparse


    !> \brief   Copy atomic 2-electron integrals to a sparse array
    !> \authors Jakub Benda, Zdenek Masin
    !> \date    2018, 2020 - 2025
    !>
    !> Retrieve non-zero [pq|rs] atomic integrals for given block index (ie. for given p,q), restricted to the upper
    !> triangle (r <= s). Write them into the supplied (pre-allocated) arrays Ri(:), Rv(:) as elements of sparse 4-index tensor;
    !> Ri(:) will contain rectangular zero-based multi-indices of the non-zero elements, Rv(:) will contain the values
    !> of the integrals.
    !>
    !> If the optional parameter `block_offset` is specified, implying prefetching in MPI mode, both triangles are fetched.
    !>
    !> \param this         Reference to the parent type.
    !> \param ao_integrals Pointer to the p2d_array_obj structure where the AO integrals are stored
    !> \param iblk         Zero-based rectangular multi-index of the block to retrive.
    !> \param rbeg         Starting index for arrays Ri, Rv. In case block_offset is present this value
    !>                     has the meaning of the offset for storage of the integrals in the Rv array
    !>                     which is gradually accumulated in the main loop over iblk.
    !> \param Rn           On return, number of elements written into the arrays Rv and Ri.
    !> \param Ri           Zero-based multi-index for each element in Rv.
    !> \param Rv           Array of non-zero atomic two-electron integrals.
    !> \param skip2ec      Whether to skip [CC|CT] and [CC|CC] integrals (mostly .TRUE.).
    !> \param index_scheme AO integral indexing scheme (1 = standard, 2 = compact).
    !> \param block_offset If present then the integrals from ao_integrals are copied to the block_offset,
    !>                     and Rv arrays in the blocked storage format used in case of nprocs > 1.
    !>                     In this case it is required that the ao_integrals are also stored in
    !>                     the block format. The array Ri is not touched at all (may be unallocated on input).
    !>
    subroutine fetch_atomic_integrals_block(this, ao_integrals, iblk, rbeg, Rn, Ri, Rv, skip2ec, index_scheme, block_offset)

        use atomic_basis_gbl, only: compact_integral_index
        use gto_routines_gbl, only: index_1p_continuum, find_mapping
        use const_gbl,        only: two_el_ints_prefetched

        class(molecular_orbital_basis_obj) :: this

        type(p2d_array_obj), pointer :: ao_integrals
        integer,                    intent(in)    :: iblk, index_scheme
        integer,                    intent(inout) :: rbeg
        integer,                    intent(out)   :: Rn
        real(kind=cfp),             intent(inout) :: Rv(:,:)
        integer,                    intent(inout) :: Ri(:,:)
        logical,                    intent(in)    :: skip2ec
        integer,       allocatable, optional      :: block_offset(:)


        integer :: no_ao, no_T_ao, last_CT_fn, n_prec_ints, n_TT_pairs, p, q, r, s, u, v, pq, pqm, rs, rsm, I
        integer :: ptype, qtype, rtype, stype, n_shell_rs
        logical :: two_p_continuum
        ! needed for indexing in case of distributed (MPI) calculation
        logical :: ao_is_local, source_ao_integrals_prefetched
        logical :: is_CCTT
        integer :: ind_orig(4), n(4), map(4), n_map(3), r_shell, s_shell, p_shell, q_shell, j, k, l, quartet_ind, ind_triang
        integer :: p_rel, q_rel, r_rel, s_rel, ind_start_r, ind_start_s, ind_start_p, ind_start_q, i1, i2, i3, i4
        integer :: n_shell_r, n_shell_s, n_shell_p, n_shell_q
        integer :: block_offset_source, output_storage_format
        character(len=line_len), allocatable :: column_descriptor(:)

        if (.not.associated(ao_integrals)) then
           call xermsg ('molecular_orbital_basis_obj', 'fetch_atomic_integrals_block', &
                        'On input the pointer to the AO integrals has not been associated.', 1, 1)
        endif

        no_ao       = this % ao_basis % number_of_functions
        n_prec_ints = this % ao_basis % n_prec_ints
        n_TT_pairs  = this % ao_basis % n_TT_pairs
        last_CT_fn  = this % ao_basis % last_CT_fn

        ! the last atomic T orbital
        no_T_ao = this % ao_basis % n_target_fns

        ! figure out if the AO integrals stored in memory are distributed (multiple MPI tasks) or local (single MPI task)
        ao_is_local = .false.
        if (this % ao_integral_storage % in_memory()) then
           ao_is_local = .not.(ao_integrals % have_offsets())
        else
           call xermsg ('molecular_orbital_basis_obj', 'fetch_atomic_integrals_block', &
                        'Reading of AO integrals from disk not implemented.', 2, 1)
        endif

        if (present(block_offset)) then
           ! The AO integrals will be copied into the Ri, Rv arrays in the block format
           ! with the global offset for storage in Rv given by rbeg on input.
           if (ao_is_local) then
              call xermsg ('molecular_orbital_basis_obj', 'fetch_atomic_integrals_block', &
                           'block form requires AO integrals on input stored in the block format.', 3, 1)
           endif
        endif

        ! Figure out the format of the source AO integrals
        call ao_integrals % get_column_descriptor (column_descriptor)

        if (column_descriptor(1) == two_el_ints_prefetched) then
           source_ao_integrals_prefetched = .true.
        else
           source_ao_integrals_prefetched = .false.
        endif

        if (present(block_offset)) then
           ! block format in blocks by the pq-index
           output_storage_format = 1
        else
           ! sparse format (index/value) pair
           output_storage_format = 0
        endif

        block_offset_source = -1
        Rn = 0
        p = 1 + iblk / no_ao
        q = 1 + mod(iblk, no_ao)

        ptype = merge(0, 1, p <= no_T_ao)  ! T (= 0) or C (= 1)
        qtype = merge(0, 1, q <= no_T_ao)  ! T (= 0) or C (= 1)

        !pq = max(p,q) * (max(p,q) - 1) / 2 + min(p,q)

        do r = 1, no_ao

            rtype = merge(0, 1, r <= no_T_ao)  ! T (= 0) or C (= 1)

            ! skip CCCT and CCCC combinations
            if (skip2ec .and. ptype + qtype + rtype == 3) then
                exit
            end if

            do s = merge(r, 1, output_storage_format == 0), no_ao  ! return only non-redundant half (unless prefetching for MPI)

                stype = merge(0, 1, s <= no_T_ao)  ! T (= 0) or C (= 1)
                !rs = max(r,s) * (max(r,s) - 1) / 2 + min(r,s)

                ! skip TCCC, CTCC and CCTC combinations
                if (skip2ec .and. ptype + qtype + rtype + stype == 3) then
                    exit
                end if

                ! get mapped indices
                rsm = this%ao_basis%function_pair_ordered_index(r,s)
                pqm = this%ao_basis%function_pair_ordered_index(p,q)
   
                ! compactify the overall non-redundant multi-index
                u = max(pqm, rsm) ; v = min(pqm, rsm) ; I = u * (u - 1) / 2 + v
                ind_triang = I

                ! AO integral indexing depends on whether all AO integrals are kept by each MPI task (are local) or not
                if (ao_is_local) then

                   if (index_scheme == 2) then
                       ! compact integral indexing
                       I = compact_integral_index(p, q, r, s)
                   else if (skip2ec) then
                       ! special indexing for CCTT class when 2p continuum is disabled
                       if (v <= n_TT_pairs .and. u > last_CT_fn) I = n_prec_ints + v + n_TT_pairs * (u - last_CT_fn - 1)
                       if (u <= n_TT_pairs .and. v > last_CT_fn) I = n_prec_ints + u + n_TT_pairs * (v - last_CT_fn - 1)
                   end if
   
                else ! AO integrals are distributed among the MPI tasks

                   ! special indexing of the shell quartets for CCTT class when 2p continuum is disabled
                   is_CCTT = .false.
                   if (skip2ec) then
                       if (v <= n_TT_pairs .and. u > last_CT_fn) is_CCTT = .true.
                       if (u <= n_TT_pairs .and. v > last_CT_fn) is_CCTT = .true.
                   end if

                   if (source_ao_integrals_prefetched) then

                      ! The integrals from ao_integrals % a array were stored at the prefetching step (i.e. when this routine was called with block_offset).
                      ! The integrals are stored in blocks by pq-index.
                      block_offset_source = ao_integrals % block_offset(iblk+1)
                      if (block_offset_source == -1) cycle

                      I = block_offset_source + rsm

                   else

                      ! index of the shell the p,q,r,s function is part of
                      p_shell = this % ao_basis % indices_to_shells(1,p)
                      q_shell = this % ao_basis % indices_to_shells(1,q)
                      r_shell = this % ao_basis % indices_to_shells(1,r)
                      s_shell = this % ao_basis % indices_to_shells(1,s)
   
                      i = max(p_shell,q_shell)
                      j = min(p_shell,q_shell)
                      k = max(r_shell,s_shell)
                      l = min(r_shell,s_shell)

                      !quartet_ind = index of the quartet of shells to which the [pq|rs] integral belongs
                      quartet_ind = index_1p_continuum(this % ao_basis % ordered_shell_pairs, &
                                             i, j, k, l, is_CCTT, &
                                             this % ao_basis % last_CT_sh, &
                                             this % ao_basis % n_prec_sh, &
                                             this % ao_basis % n_TT_sh_pairs)
   
                      block_offset_source = ao_integrals % block_offset(quartet_ind)
   
                      if (block_offset_source == -1) then
                         ! the requested integral is not kept by this task
                         cycle
                      endif

                      ! starting index of the functions in the shell the p,q,r,s function is part of
                      ind_start_p = this % ao_basis % shell_descriptor(4,p_shell)
                      ind_start_q = this % ao_basis % shell_descriptor(4,q_shell)
                      ind_start_r = this % ao_basis % shell_descriptor(4,r_shell)
                      ind_start_s = this % ao_basis % shell_descriptor(4,s_shell)
   
                      ! number of functions in the shell of which the p,q,r,s function is part of
                      n_shell_p = this % ao_basis % shell_descriptor(5,p_shell)
                      n_shell_q = this % ao_basis % shell_descriptor(5,q_shell)
                      n_shell_r = this % ao_basis % shell_descriptor(5,r_shell)
                      n_shell_s = this % ao_basis % shell_descriptor(5,s_shell)
   
                      ! index of the p,q,r,s function within the shell it belongs to
                      p_rel = this % ao_basis % indices_to_shells(2,p)
                      q_rel = this % ao_basis % indices_to_shells(2,q)
                      r_rel = this % ao_basis % indices_to_shells(2,r)
                      s_rel = this % ao_basis % indices_to_shells(2,s)
   
                      !Each function p,q,r,s corresponds to a given shell. In the section below we permute the p,q,r,s indices
                      !to the standard order (determined by the starting indices of the functions in each shell)
                      !which was used to save the integrals (see atomic_orbital_basis_obj % two_electron_integrals).
                      ind_orig(1:4) = (/ind_start_p,ind_start_q,ind_start_r,ind_start_s/)
                      n(1:4) = (/n_shell_p,n_shell_q,n_shell_r,n_shell_s/)
   
                      !Map the order of the corresponding p,q,r,s shells to the order in which the integrals within the shell are saved.
                      call find_mapping(ind_orig,n,n_map,map)
   
                      !Permute the relative indices into the order in which the integrals were saved
                      ind_orig(1:4) = (/p_rel,q_rel,r_rel,s_rel/)
                      i1 = ind_orig(map(1))
                      i2 = ind_orig(map(2))
                      i3 = ind_orig(map(3))
                      i4 = ind_orig(map(4))

                      !Compute the index of the integral [pq|rs] within its own quartet of shells and add to it
                      !the offset for the corresponding quartet of shells.
                      I = block_offset_source + i1 + n_map(1)*(i2-1) + n_map(2)*(i3-1) + n_map(3)*(i4-1)

                   endif

                endif

                ! sanity check - the index must not exceed the size of the atomic integral storage
                if (I > size(ao_integrals % a, 1)) then
                    write (stdout, '("Inconsistency: Index overflow!")')
                    write (stdout, '("- Number of AOs: ",I0)') no_ao
                    write (stdout, '("- n_prec_ints: ",I0)') n_prec_ints
                    write (stdout, '("- n_TT_pairs: ",I0)') n_TT_pairs
                    write (stdout, '("- last_CT_fn: ",I0)') last_CT_fn
                    write (stdout, '("- p, q, r, s: ",I0,1x,I0,1x,I0,1x,I0)') p, q, r, s
                    write (stdout, '("- u, v, I: ",I0,1x,I0,1x,I0)') u, v, I
                    flush (stdout)
                    stop 1
                end if

                select case (output_storage_format)

                case (1)

                   ! Store this integral into the dense array in blocks by the pq-index: the rs-part is the same for all rs-indices.
                   ! I.e. this is the "prefetched" storage format. It is assumed that the block_offset is preset externally!
                   j = block_offset(iblk+1) + rsm

                   if (j > size(Rv,1)) then
                      write(stdout,'("error -out of bounds",4i4,2i10)') p,q,r,s,j,I
                      flush(stdout)
                   endif
                   Rv(j, 1) = ao_integrals % a(I, 1)

                   ! On the prefetching stage Rn has the meaning of the largest rs-index encountered
                   Rn = max(Rn, rsm)

                case (0)

                   ! store this integral into the sparse array
                   if (I > 0) then
                      if (ao_integrals % a(I, 1) /= 0) then
                         Rn = Rn + 1
                         Ri(Rn + rbeg - 1, 1) = (((p-1) * no_ao + (q-1)) * no_ao + (r-1)) * no_ao + (s-1)  ! zero-based multi-index
                         Rv(Rn + rbeg - 1, 1) = ao_integrals % a(I, 1)
                      end if
                   end if

                end select

            end do ! s

        end do  ! r

    end subroutine fetch_atomic_integrals_block


    !> \brief   Extract non-zero elements from a dense matrix
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Given a dense matrix A, copy its non-zero elements into the array Sv, with corresponding CSC indices in Sp and Sj.
    !>
    !> \param Sp  Index of the first element in Sj/Sv in the given column. If two consecutive elements in Sp have the same
    !>            value, it means that the column has zero length.
    !> \param Sj  Row index of the corresponding element in Sv.
    !> \param Sv  Non-zero elements of A.
    !> \param A   Dense matrix (column major storage).
    !>
    subroutine sparse_from_dense(Sp, Sj, Sv, A)

        real(kind=cfp), allocatable, intent(in)  :: A(:,:)
        real(kind=cfp), allocatable, intent(out) :: Sv(:)
        integer,        allocatable, intent(out) :: Sp(:), Sj(:)

        integer :: row, col, ierr, m, n, pos, nnz

        m   = size(A, 1)     ! leading dimension of the matrix
        n   = size(A, 2)     ! other dimension
        nnz = count(A /= 0)  ! number of non-zero elements

        allocate(Sp(n + 1), Sj(nnz), Sv(nnz), stat = ierr)
        if (ierr /= 0) call xermsg ('molecular_orbital_basis_obj', 'sparse_from_dense', 'Memory allocation failure.', 1, 1)

        pos = 0
        Sp(1) = 1

        do col = 1, n
            do row = 1, m
                if (A(row,col) /= 0) then
                    pos = pos + 1
                    Sj(pos) = row
                    Sv(pos) = A(row,col)
                end if
            end do
            Sp(col + 1) = pos + 1
        end do

    end subroutine sparse_from_dense


    !> \brief   Convert sparse triangle to a sparse rectangular matrix
    !> \authors J Benda
    !> \date    2025
    !>
    !> Unpacks a sparse triangular matrix in COO format (Ri, Rv) to a rectangular symmetric sparse matrix in CSC format
    !> (Wp, Wj, Wv).
    !>
    !> \param[in]  nAO Number of atomic orbitals (dimension of the matrix).
    !> \param[in]  Ri  COO row-column multi-index indexing non-zero elements in the sparse triangle.
    !> \param[in]  Rv  Non-zero elements in sparse triangle matrix.
    !> \param[out] Wp  CSC column pointers.
    !> \param[out] Wj  CSC row indices of non-zero elements in the sparse matrix.
    !> \param[out] Wv  Non-zero elements in output sparse matrix.
    !>
    subroutine symmetrize_sparse_triangle(nAO, Ri, Rv, Wp, Wj, Wv, rbeg, rend0)

        use sort_gbl, only: heap_sort_int_float

        integer,   intent(in)    :: nAO, rbeg, rend0
        integer,   intent(in)    :: Ri(:)
        real(cfp), intent(in)    :: Rv(:)
        integer,   intent(inout) :: Wp(:), Wj(:)
        real(cfp), intent(inout) :: Wv(:)

        integer :: ab, c, cd, d, p, rcur, rend

        rend = rend0

        ! mirror off-diagonal entries, append at the end of the provided storage
        do rcur = rbeg, rend0
            ! decode the multi-index
            cd = mod(Ri(rcur), nAO * nAO)
            ab = Ri(rcur) - cd
            c = 1 + cd / nAO
            d = 1 + mod(cd, nAO)
            ! copy this non-zero element to output arrays
            Wj(rcur - rbeg + 1) = Ri(rcur)
            Wv(rcur - rbeg + 1) = Rv(rcur)
            ! append non-diagonal element at end
            if (c /= d) then
                rend = rend + 1
                Wj(rend - rbeg + 1) = ab + (d - 1)*nAO + c - 1
                Wv(rend - rbeg + 1) = Rv(rcur)
            end if
        end do

        ! sort by COO multi-index
        call heap_sort_int_float(rend - rbeg + 1, Wj, Wv)

        ! convert to CSC format
        p = 1
        do rcur = rbeg, rend
            ! decode the multi-index
            cd = mod(Wj(rcur - rbeg + 1), nAO * nAO)
            c = 1 + cd / nAO
            d = 1 + mod(cd, nAO)
            ! fill all column pointers up to the current column
            do while (p <= c)
                Wp(p) = rcur - rbeg + 1
                p = p + 1
            end do
            ! set the row index for this element
            Wj(rcur - rbeg + 1) = d
        end do

        ! finalize the column pointer array (point beyond the storage)
        Wp(p:nAO + 1) = (rend + 1) - rbeg + 1

    end subroutine symmetrize_sparse_triangle


    !> \brief   Transform both indices of each block of a large sparse matrix by the (sparse) coefficient matrix
    !> \authors Jakub Benda, Zdenek Masin
    !> \date    2018, 2020
    !>
    !> In the first step of the transformation (ie. [pq|rs] -> [pq|kl)), the input arrays Ri, Rv are ignored and the sparse
    !> blocks Rb of the 4-index tensor R are retrieved directly from the atomic integrals storage. Ri, Rv are then allocated
    !> and filled with the transformed data.
    !>
    !> In the second step, the sparse blocks Rb are truly read from the input arrays Ri, Rv. These are then destroyed and
    !> reallocated to the proper size before returning the transformed data.
    !>
    !> Transformation of every block amounts to two sparse matrix multiplications (with some restrictions on non-zero
    !> elements). The result of each such transformation is a subset of a dense block. Each transformed block is stored in
    !> a separate array in a linked list data structure to avoid frequent reallocations of a growing array of elements.
    !> These are merged to Ri, Rv at the very end of the transformation step.
    !>
    !> Note that this subroutine deliberately transforms only half of all blocks. Specifically, the input in the first
    !> step are [pq|rs] and only blocks with p >= q are considered by this function. Similarly, the subroutine
    !> `fetch_atomic_integrals_block` called from here also provides only triangles. The triangles are symmetrized
    !> prior to the transformation of each block. On exit from this subroutine we have semi-transformed
    !> integrals [pq|kl), once again restricted do k >= l due to a triangular reduction in `transform_one_index`.
    !> Due to the subsequent sorting of the semitransformed integrals between the first and the second step of the sparse
    !> transformation, the inputs to the second step will then be the integrals (kl|pq], where only k >= l and p >= q are
    !> non-zero. This is perfectly compatible with the requirements above, so no integrals are calculated needlessly.
    !>
    !> \todo The merging could be parallelized, too, if each block also stored its offset.
    !>
    !> \param Nbk Number of blocks along the diagonal (= nAO in the first step, nMO in the second step).
    !> \param nAO Number of atomic orbitals (dimension of the index to transform).
    !> \param nMO Number of molecular orbitals (dimension of the index after transformation).
    !> \param mo_integrals    Storage object for working with the (semi)transformed two-electron integrals.
    !> \param Cp  Column pointers in Ci (ie. where the coefficients corresponding to i-th MO start in Cv).
    !> \param Cj  Row indices corresponding to elements in Cv.
    !> \param Cv  Structurally non-zero elements of C (the MO <- AO coefficient matrix).
    !> \param step            First or second step of the transformation.
    !> \param tolerance       Minimal absolute value of integral to keep it.
    !> \param skip2ec         Whether to skip CCCT and CCCC combinations.
    !> \param max_buffer_size Maximum size (in MiB) of the auxiliary buffer array: used only in MPI mode.
    !> \param index_scheme    AO integrals indexing scheme (1 = standard, 2 = compact).
    !> \param keep_ao_integrals Logical switch to preserve or not the array of AO integrals.
    !>
    subroutine transform_two_indices(this, Nbk, nAO, nMO, mo_integrals, Cp, Cj, Cv, step, tolerance, skip2ec, max_buffer_size, &
                                     index_scheme, keep_ao_integrals)

        use const_gbl            , only: abel_prod_tab, two_el_ints_prefetched
        use omp_lib              , only: omp_get_num_threads, omp_get_thread_num, omp_get_wtime
        use sparse_block_storage_gbl, only: sparse_block_storage_obj
        use special_functions_gbl, only: ipair

        implicit none

        class(molecular_orbital_basis_obj) :: this
        type(sparse_block_storage_obj), target, intent(inout) :: mo_integrals

        real(kind=cfp), allocatable, intent(in)    :: Cv(:)
        integer,        allocatable, intent(in)    :: Cp(:), Cj(:)
        integer,                     intent(in)    :: Nbk, nAO, nMO, step, index_scheme
        logical,                     intent(in)    :: skip2ec
        real(kind=cfp),              intent(in)    :: tolerance, max_buffer_size
        logical,                     intent(in)    :: keep_ao_integrals

        real(wp) :: start_t, end_t, start_prefetch, end_prefetch
        real(kind=cfp) :: x
        type(p2d_array_obj), pointer         :: ao_integrals => null()
        type(p2d_array_obj), target          :: ao_integrals_prefetched
        character(len=line_len), allocatable :: column_descriptor(:)
        real(kind=cfp), allocatable  :: Vv(:,:), Wv(:,:)
        real(kind=cfp), pointer      :: Rv(:,:)
        integer,        allocatable  :: Vp(:), Vj(:,:), Wq(:), Wj(:,:), Rp(:), Rj(:,:), blocks(:)
        integer,        pointer      :: Ri(:,:)
        integer :: a, b, c, d, ab, cd, i, j, iblk, nblk, mx, nnz, nonzeros, discarded, ithread, nthreads, ierr, d2
        integer :: rbeg, rcur, rend, blocksym, blockcnt, TA_orb(9), CA_orb(9), TM_orb(9), CM_orb(9), iT, iC, iRp, orb
        integer :: iblk_start, iblk_end, no_blocks, i_max_rank, all_iblk_end(nprocs), sendcount_task(nprocs)
        integer :: n_pq, pq, p, q, n_prefetched, pq_type
        logical :: mo_mpi_redistribution, gather_all_AO_integrals, split_within_block

        start_t = omp_get_wtime()

        blockcnt = 0    ! combined CC/CT/TT type of the orbitals with the block indices
        blocksym = 0    ! combined symmetry of the orbitals with the block indices (only used in the second stage)

        ! Get starting (and first continuum) molecular orbitals of individual symmetries
        !
        ! - For example, when the number of target+continuum orbitals per symmetry is 5+4,0+3,1+3,0+0,0+3,0+0,0+0,0+0,
        !   so in total 19 orbitals, then the two helper offset arrays will look like this:
        !
        !      TM_orb = 1  10  13  17  17  20  20  20  20
        !      CM_orb = 6  10  14  17  17  20  20  20  20

        iT = 1 ; TM_orb = this % number_of_functions + 1
        iC = 1 ; CM_orb = this % number_of_functions + 1
        sym_loop: do i = 1, this % no_irr
            orb_loop: do j = 1, this % orbital_data(i) % number_of_functions
                orb = this % relative_to_absolute(j, i)
                ! first orbital of each symmetry (if there is any) will be stored in TM_orb for that symmetry and
                ! also all previous empty symmetries; it will be also stored in CM_orb for all previous empty symmetries
                if (j == 1) then
                    do while (iT <= i) ; TM_orb(iT) = orb ; iT = iT + 1 ; end do
                    do while (iC <  i) ; CM_orb(iC) = orb ; iC = iC + 1 ; end do
                end if
                ! first continuum orbital of each symmetry (if there is any) will be stored in CM_orb
                if (this % is_continuum(orb)) then ; CM_orb(iC) = orb ; iC = iC + 1 ; exit orb_loop ; end if
            end do orb_loop
        end do sym_loop

        ! Now do the same for atomic orbitals.
        !
        ! This is simpler, because the atomic orbitals are not sorted by symmetry; instead
        ! they are stored in one consecutive array with target orbitals at the beginning.

        TA_orb = this % ao_basis % number_of_functions + 1
        CA_orb = this % ao_basis % number_of_functions + 1

        TA_orb(1)  = 1
        do j = 1, this % ao_basis % number_of_functions
            i = this % ao_basis % indices_to_shells(1, j)
            if (this % ao_basis % shell_descriptor(3, i) == 1) CA_orb(1) = min(j, CA_orb(1))
        end do

        if (step == 1) then
            write (level1, '("Transforming the first pair of atomic indices, [AA|AA] -> [AA|MM)")')
        else
            write (level1, '("Transforming the second pair of atomic indices, (MM|AA] -> (MM|MM)")')
        end if
        write (level2, '("Atomic orbital offsets: ")')
        write (level2, '(5x,"T",I6,"   (",I0,")")') TA_orb(1), TA_orb(2)
        write (level2, '(5x,"C",I6)') CA_orb(1)
        write (level2, '("Molecular orbital offsets per symmetry: ")')
        write (level2, '(5x,"T",8I6,"   (",I0,")")') TM_orb(1:8), TM_orb(9)
        write (level2, '(5x,"C",8I6)') CM_orb(1:8)

        ! First, find out the number of blocks so that we can process them concurrently.
        !
        ! - In the first step (ie. [pq|rs] -> [pq|kl)), there will be at most "nthreads" blocks in memory at a single time; we
        !   allocate space for them here. The total number of blocks is known, too: one for each pair of atomic orbitals.
        ! - In the second step (ie. (kl|pq] -> (kl|ij)), the number of blocks is retrieved from the R-tensor itself, and blocks
        !   offsets in the linear storage arrays are stored for later convenience.

        ! Modified, Josh Forer 2026, omp master -> omp masked
        !$omp parallel shared(nthreads)
        !$omp masked
        nthreads = omp_get_num_threads()
        !$omp end masked
        !$omp end parallel

        if (step == 1) then
            nblk = nAO * nAO
            allocate (mo_integrals%Ri(nthreads * nAO * nAO, 1),stat=ierr)  ! This is enough to hold even dense blocks, ...
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 1 failure.', ierr, 1)
            endif
            allocate (mo_integrals%Rv(nthreads * nAO * nAO, 1),stat=ierr)  ! ... so any sparse block should fit in without problems.
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 2 failure.', ierr, 1)
            endif
        else
            allocate (blocks(Nbk * Nbk + 1),stat=ierr)
            if (ierr /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 3 failure.', ierr, 1)
            endif
            nblk = 0 ; mx = -1
            do rcur = 1, size(mo_integrals%Ri, 1)
                ab = mo_integrals%Ri(rcur, 1) / (nAO * nAO)
                if (mx /= ab) then
                    nblk = nblk + 1
                    blocks(nblk) = rcur
                end if
                mx = ab
            end do
            blocks(nblk + 1) = size(mo_integrals%Ri, 1) + 1  ! points after the end of Ri
        end if

        Ri => mo_integrals%Ri
        Rv => mo_integrals%Rv

        ! And here starts the parallel transformation itself. The threads sequentially pick one block at a time and process
        ! it, transforming it from atomic to molecular indices.
        !
        ! - In the first step (ie. [pq|rs] -> [pq|kl)), the blocks to transform are constructed on the fly by the very thread
        !   that is about to do the transformation.
        ! - In the second step (ie. (kl|pq] -> (kl|ij)), the blocks are retrieved from the R-tensor as obtained from the subroutine
        !   arguments.

        discarded = 0
        nonzeros = 0

        iblk_end = 0

        ao_integrals => null()
 
        do !blocks of AO integrals (iblk_start -> iblk_end) which are gathered sequentially from all tasks

           flush(stdout)

           ! The index iblk goes over all pq pairs. In MPI mode we loop over them in chunks of length iblk_end-iblk_start+1.
           ! In case of serial mode (one MPI task) we loop over all pq pairs at once.
           iblk_start = iblk_end + 1

           if (nprocs > 1 .and. step == 1) then

              if (iblk_start == 1) then

                 ! initialize an auxiliary p2d_array_obj with the same number of blocks as this % ao_integral_storage % integral
                 call this % ao_integral_storage % integral % get_array_dim (a, d2, no_blocks)
                 call this % ao_integral_storage % integral % get_column_descriptor (column_descriptor)

                 ! first see if we have enough space to keep all AO integrals from all tasks
                 call mpi_mod_allgather(a, sendcount_task)

                 ! total number of AOs computed by all tasks
                 b = sum(sendcount_task)
                
                 if (max_buffer_size < 0.0_cfp) then
                    ! allocate as much memory as we ideally need
                    i = b
                 else
                    ! Maximum number of elements allowed for the AO integrals.
                    ! Note that this takes into account the array block_offset_all allocated by merge_indexed_array_from_all_tasks.
                    i = int((max_buffer_size*Mib - (nprocs+1)*nblk*storage_unit_int)/(cfp_bytes*1.0_cfp))
   
                    ! I need at least this amount of memory: each task must be able to store one complete block of [pq|rs] integrals with pq-pair fixed.
                    n_pq = nAO * (nAO + 1) / 2
                    i = max(i, n_pq*nprocs)
   
                    ! i == b is the optimal situation
                    i = min(i, b)
                 endif

                 if (i == b) then
                    gather_all_AO_integrals = .true.

                    iblk_end = nblk

                    write(level1,'(/,"All AO integrals will be gathered at once.")')

                 else
                    gather_all_AO_integrals = .false.

                    ! maximum number of elements that each task can insert into the ao_integrals_prefetched%a(:,1) array
                    i_max_rank = i / nprocs

                    ! each pq-block is saved separately
                    no_blocks = nblk

                    column_descriptor = two_el_ints_prefetched

                    write(level1,'(/,"Total number of pq blocks of AO integrals needed: ",i0)') nblk

                 endif

                 write(level1,'(  "Allocating temporary AO storage:",/)')

                 ierr = ao_integrals_prefetched % init (i, d2, no_blocks, column_descriptor)
                 if (ierr /= 0) then
                    call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', &
                                 'ao_integrals_prefetched failed to initialize.', ierr, 1)
                 endif

              endif

              start_prefetch = omp_get_wtime()

              ! AO integrals calculated by myrank
              ao_integrals => this % ao_integral_storage % integral

              ! Copy as many of my AO integrals as we can from ao_integrals into ao_integrals_prefetched
              ao_integrals_prefetched % a = 0.0_cfp
              ao_integrals_prefetched % block_offset = -1

              if (gather_all_AO_integrals) then

                 ! offset for my AO integrals in the combined array
                 rbeg = sum(sendcount_task(1:myrank))

                 ! copy all my AO integrals into the combined array
                 ao_integrals_prefetched % a(rbeg + 1:rbeg + sendcount_task(myrank+1), 1) = &
                            ao_integrals % a(       1:       sendcount_task(myrank+1), 1)

                 ! shift the pointers to the data from ao_integrals % a
                 where (ao_integrals % block_offset /= -1)
                    ao_integrals_prefetched % block_offset = rbeg + ao_integrals % block_offset
                 elsewhere
                    ao_integrals_prefetched % block_offset = -1
                 end where

                 i_max_rank = sendcount_task(myrank+1)
                 split_within_block = .false.

              else !we have to proceed by gathering the AOs in pq blocks

                 split_within_block = .true.

                 ! AO integrals for myrank will be stored in the section of ao_integrals_prefetched % a(rbeg+1:rbeg+i_max_rank)
                 rbeg = myrank * i_max_rank

                 iblk_end = min(iblk_start + (i_max_rank / n_pq) - 1, nblk)

                 n_prefetched = 0
                 iblk = iblk_start

                 do

                    pq = iblk-1
                    p = 1 + pq / nAO
                    q = 1 + mod(pq, nAO)
                    pq = max(p,q) * (max(p,q) - 1) / 2 + min(p,q)

                    ! determine how many rs-indices will be needed for the current pq (iblk) pair (we use the ordered rs-index)
                    nnz = n_pq
                    if (skip2ec) then

                       pq_type = this%ao_basis%get_function_pair_type(p,q)

                       if (pq_type == 1) then     !pq is of TT type
                          nnz = n_pq                         !last rs pair needed is of type CC
                       elseif (pq_type == 2) then !pq is of CT type
                          nnz = this % ao_basis % last_CT_fn !last rs pair needed is of type CT
                       elseif (pq_type == 3) then !pq is of CC type
                          nnz = this % ao_basis % n_TT_pairs !last rs pair needed is of type TT
                       endif

                    endif

                    if (n_prefetched + nnz > i_max_rank) then
                       iblk_end = iblk-1
                       exit
                    else
                       ao_integrals_prefetched % block_offset(iblk) = rbeg
                       n_prefetched = n_prefetched + nnz
                       rbeg = rbeg + nnz
                    endif

                    iblk = iblk + 1
                    if (iblk > nblk) then 
                       iblk_end = iblk-1
                       exit
                    endif

                 enddo

                 !loop over pq (=iblk-1) pairs of indices
                 !$omp parallel &
                 !$omp& shared(iblk_start, iblk_end, ao_integrals, ao_integrals_prefetched, Ri, skip2ec) &
                 !$omp& private(iblk, rbeg)
                 !$omp do
                 do iblk = iblk_start, iblk_end

                    rbeg = ao_integrals_prefetched % block_offset(iblk)

                    ! get AO integrals [pq|rs] for a fixed pq (=iblk-1) pair and all rs pairs excluding the pq-redundant part
                    call this % fetch_atomic_integrals_block(ao_integrals, iblk - 1, rbeg, nnz, &
                                                             Ri, ao_integrals_prefetched % a, skip2ec, index_scheme, &
                                                             ao_integrals_prefetched % block_offset)
                 enddo
                 !$omp end do
                 !$omp end parallel
   
                 write(level2,'("Number of pq blocks of AO integrals prefetched: ",i0)') iblk_end-iblk_start+1
   
                 call mpi_mod_allgather(iblk_end, all_iblk_end)
   
                 ! we have to make sure all processes use the same range of pq-indices
                 iblk_end = minval(all_iblk_end)
   
                 write(level2,'("pq block indices synchronized for all tasks: ",2i10)') iblk_start, iblk_end

              endif !gather_all_AO_integrals

              ! CRUCIAL PART:
              ! Every tasks assembles a single AO integrals array from the integrals gathered from all other tasks.
              ! Note that inside this routine a potentially large temporary array holding the block offsets from all tasks is allocated.
              call ao_integrals_prefetched % merge_indexed_array_from_all_tasks(i_max_rank, 1, split_within_block)

              ! Now all tasks have all AO integrals for the block of pq indices in the range iblk = [iblk_start; iblk_end] assembled in a single block-indexed
              ! array structure so we can proceed with the transformation using the AO integrals from ao_integrals_prefetched as the source.
              ! The MO integrals computed by each task are split below by the 'kl' pair of MO orbitals.
              ao_integrals => ao_integrals_prefetched

              end_prefetch = omp_get_wtime()

              write (level2, '("Prefetch took [s]: ",F25.15,/)') end_prefetch - start_prefetch

           else

              ! in this case the single task has all AO integrals (step == 1) or we're in step == 2.
              iblk_end = nblk

              ao_integrals => this % ao_integral_storage % integral

           endif !(nprocs > 1 .and. step == 1)

           !$omp parallel &
           !$omp& shared  (this, Ri, Rv, Nbk, nblk, blocks, step, nAO, nMO, TA_orb, CA_orb, &
           !$omp&          TM_orb, CM_orb, skip2ec, Cp, Cj, Cv, ao_integrals, iblk_start, iblk_end) &
           !$omp& private (iblk, rbeg, rcur, rend, ab, cd, a, b, c, d, i, j, blocksym, blockcnt, orb, &
           !$omp&          iRp, Rp, Rj, Vp, Vv, Vj, Wq, Wv, Wj, nnz, ithread, ierr, mo_mpi_redistribution) &
           !$omp& reduction (+ : discarded, nonzeros)
   
           ithread = omp_get_thread_num()
   
           ! allocate workspaces (again, use enough space for dense blocks)
           allocate(Rp(nAO + 1), Rj(nAO * nAO, 1), stat=ierr)
           if (ierr /= 0) then
              call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 4 failure.', ierr, 1)
           endif
           allocate(Vp(nMO + 1), Vv(nAO * nMO, 1), Vj(nAO * nMO, 1), stat=ierr)
           if (ierr /= 0) then
              call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 5 failure.', ierr, 1)
           endif
           allocate(Wq(max(nAO, nMO) + 1), Wv(max(nAO, nMO)**2, 1), Wj(max(nAO, nMO)**2, 1), stat=ierr)
           if (ierr /= 0) then
              call xermsg ('molecular_orbital_basis_obj', 'transform_two_indices', 'Memory allocation 6 failure.', ierr, 1)
           endif
   
           !$omp do schedule (dynamic, 1)
   
           ! loop over the current block of the R tensor
           R_loop: do iblk = iblk_start, iblk_end
   
               ! get section of data corresponding to this block
               if (step == 1) then
                   rbeg = nAO * nAO * ithread + 1
                   rend = nAO * nAO * (ithread + 1)
                   ! get AO integrals [pq|rs] for a fixed pq (=iblk-1) pair and all rs pairs in the sparse storage format
                   call this % fetch_atomic_integrals_block(ao_integrals, iblk - 1, rbeg, nnz, Ri, Rv, skip2ec, index_scheme)
                   ! in case of a distributed AO integral calculation some tasks may not have any integrals for certain pq pairs
                   if (nnz == 0) cycle
                   rend = rbeg + nnz - 1
               else
                   rbeg = blocks(iblk)
                   rend = blocks(iblk + 1) - 1
               end if
   
               ! decode the multi-index
               ab = Ri(rbeg, 1) / (nAO * nAO)
               a = 1 + ab / Nbk
               b = 1 + mod(ab, Nbk)

               ! transform only non-redundant half of blocks (automatically satisfied by construction in step == 2)
               if (a < b) cycle
   
               ! get CC/CT/TT (= 2/1/0) type of the two (atomic or molecular) orbitals
               blockcnt = 0
               if (step == 1) then
                   ! type of atomic orbital pair
                   !orb = max(a,b) * (max(a,b) - 1) / 2 + min(a,b)
                   blockcnt = this%ao_basis%get_function_pair_type(a,b) - 1
               else
                   ! type of molecular orbital pair
                   if (this % is_continuum(a)) blockcnt = blockcnt + 1
                   if (this % is_continuum(b)) blockcnt = blockcnt + 1
               end if
   
               ! get combined symmetry of the two (molecular) orbitals (only used in the last step)
               if (step == 2) then
                   blocksym = abel_prod_tab(count(TM_orb <= a), count(TM_orb <= b))
               end if
   
               ! unpack COO sparse triangle (Ri, Rv) to a symmetric CSC sparse matrix (Wq, Wj, Wv)
               call symmetrize_sparse_triangle(nAO, Ri(:, 1), Rv(:, 1), Wq, Wj(:, 1), Wv(:, 1), rbeg, rend)

               mo_mpi_redistribution = .false.
   
               ! multiply C . R -> V, use neither triangular reduction, nor pyramidal reduction
               call transform_one_index(nAO, nMO, Cp, Cj, Cv, &
                                        nAO, nAO, Wq, Wj, 1, Wv, &
                                                  Vp, Vj, Vv, &
                                       .false., -1, blocksym, blockcnt, &
                                       TM_orb, CM_orb, TA_orb, CA_orb, skip2ec, mo_mpi_redistribution)
   
               ! redistribute the 'kl' MO pairs over the MPI tasks in the transformation [pq|rs] -> [pq|kl)
               if (step == 1) then
                  mo_mpi_redistribution = .true.
               else
                  mo_mpi_redistribution = .false.
               endif
   
               ! multiply C . V -> W, use triangular reduction, and (in second step) also pyramidal reduction
               call transform_one_index(nAO, nMO, Cp, Cj, Cv, &
                                        nAO, nMO, Vp, Vj, 1, Vv, &
                                                  Wq, Wj, Wv, &
                                       .true.,  merge(-1, ab, step == 1), blocksym, blockcnt, &
                                       TM_orb, CM_orb, TM_orb, CM_orb, skip2ec, mo_mpi_redistribution)
   
               ! number of non-zero elements (integrals) in the transformed block
               nnz = Wq(nMO + 1) - 1
   
               ! reconstruct rectangular multi-indices from CSC indices and trim integrals by the given threshold
               nnz = 0
               do i = 1, nMO
                   do j = Wq(i), Wq(i + 1) - 1
                       if (abs(Wv(j, 1)) >= tolerance) then
                           nnz = nnz + 1
                           Wj(nnz, 1) = (ab * nMO + (i - 1)) * nMO + (Wj(j, 1) - 1)
                           Wv(nnz, 1) = Wv(j, 1)
                       end if
                   end do
               end do
               discarded = discarded + Wq(nMO + 1) - 1 - nnz
               nonzeros = nonzeros + nnz

               ! copy this thread's data to the auxiliary buffer list
               call mo_integrals % add_block(iblk, nnz, Wj(1:nnz, 1), Wv(1:nnz, 1))

           end do R_loop
   
           !$omp end parallel
  
           ! all pairs of pq indices have been processed 
           if (iblk_end == nblk) exit

        enddo !blocks of AO integrals (iblk_start -> iblk_end)

        write (level2, '("Elements discarded due to threshold: ",I0)') discarded

        ! drop atomic integrals, we will no longer need them (and every bit of free memory is needed below)
        if (step == 1) then
            if (keep_ao_integrals) then
               if (nprocs > 1) then
                  if (ao_integrals % final() /= 0) then
                      write (level2, '("WARNING: Finalization of the temporary atomic integral storage failed!")')
                  end if
               endif
            else
               if (this % ao_integral_storage % integral % final() /= 0) then
                   write (level2, '("WARNING: Finalization of the atomic integral storage failed!")')
               end if
               if (ao_integrals % final() /= 0) then
                  write (level2, '("WARNING: Finalization of the temporary atomic integral storage failed!")')
               end if
            endif
            ! also deallocate the two workspace arrays
            deallocate (mo_integrals%Ri, mo_integrals%Rv)
        end if

        ! overwrite mo_integrals%Ri and mo_integrals%Rv with the output array
        call mo_integrals % combine_blocks

        write (level1, '("Non-zero elements on exit: ",I0)') nonzeros

        end_t = omp_get_wtime()

        write (level1, '("Time spent in transform_two_indices [s]: ",F25.15,/)') end_t - start_t

    end subroutine transform_two_indices


    !> \brief   Multiplication of two sparse matrices (somewhat tweaked)
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Multiply two sparse matrices, producing a sparse matrix as a result.
    !>
    !> The subroutine allows restricting the operation to some elements to allow making use of the index symmetries
    !> of two-particle integrals:
    !>  - If "triangle" is true, then only a triangular part of the resulting matrix will be computed.
    !>  - If "pyramid" is non-negative, then only multi-indices smaller than the given value will be computed.
    !>
    !> If both the above are true, which indicates the second stage of transformation, ie. (kl|pq] -> (kl|ij), then also:
    !>  - Skip combinations of orbitals with incompatible symmetries. The integral can be only nonzero when the overall symmetry
    !>    of the four orbitals is totally symmetric.
    !>  - If "skip2ec" is true, do not calculate molecular integrals of CCCC and CCCT types.
    !>
    !> \param nAr     Number of rows in A (leading dimension).
    !> \param nAc     Number of cols in A (the other dimension).
    !> \param Ap      Positions in Av of the first element of each column of A.
    !> \param Aj      Row indices corresponding to elements in Av.
    !> \param Av      Non-zero elements of the matrix A.
    !> \param nBr     Number of rows in B (leading dimension).
    !> \param nBc     Number of cols in B (the other dimension).
    !> \param Bp      Positions in Bv of the first element of each column of B.
    !> \param Bj      Row indices corresponding to elements in Bv.
    !> \param Bv_beg  Starting index in array Bv.
    !> \param Bv      Non-zero elements of the matrix B.
    !> \param Cp      On return, Positions in Cv of the first element of each column of C.
    !> \param Cj      Row indices corresponding to elements in Cv.
    !> \param Cv      Non-zero elements of the sparse matrix product.
    !> \param Ci      Zero-based multi-indices corresponding to the elements of Cv.
    !> \param triangle  Calculate only a triangular subset of C.
    !> \param pyramid   Further restriction on calculated elements of C (limit on multi-index).
    !> \param blocksym  Combined symmetry of the other pair of orbitals.
    !> \param blockcnt  CC/CT/TT (= 2/1/0) type of the other pair of orbitals.
    !> \param TA_orb    Helper array with the index of the first orbital per symmetry. Corresponds to column index of A.
    !> \param CA_orb    Helper array with the index of the first continuum orbital per symmetry. Corresponds to column index of A.
    !> \param TB_orb    Helper array with the index of the first orbital per symmetry. Corresponds to column index of B.
    !> \param CB_orb    Helper array with the index of the first continuum orbital per symmetry. Corresponds to column index of B.
    !> \param skip2ec   Whether to skip CCCC and CCCT combinations of molecular orbitals.
    !> \param mo_mpi_redistribution Whether to cyclically redistribute pairs of columns 'AB' over MPI tasks.
    !>
    subroutine transform_one_index (nAr, nAc, Ap, Aj, Av, nBr, nBc, Bp, Bj, Bv_beg, Bv, Cp, Cj, Cv, &
                                    triangle, pyramid, blocksym, blockcnt, TA_orb, CA_orb, TB_orb, CB_orb, skip2ec, &
                                    mo_mpi_redistribution)

        use const_gbl, only: abel_prod_tab
        implicit none

        real(kind=cfp), allocatable :: Av(:), Bv(:,:) !intent(in)
        real(kind=cfp), allocatable :: Cv(:,:) !intent(inout)
        integer, allocatable :: Ap(:), Aj(:), Bp(:), Bj(:,:) !intent(in)
        integer, allocatable :: Cp(:), Cj(:,:)               !intent(inout)
        integer, intent(in)  :: nAr, nAc, nBr, nBc, Bv_beg, pyramid, blocksym, blockcnt
        integer, intent(in)  :: TA_orb(9), CA_orb(9), TB_orb(9), CB_orb(9)
        logical, intent(in)  :: triangle, skip2ec, mo_mpi_redistribution

        real(kind=cfp) :: x
        integer :: apos, acur, aend, arow, acol, bpos, bcur, brow, bend, bcol, cpos
        integer :: asym,  bsym    ! symmetry m-value of the current A/B orbital (column index)
        integer :: atype, btype   ! T/C (= 0/1) type of the current A/B orbital (column index)
        integer :: sym            ! combined orbital symmetry of the block and of the A orbital
        integer :: A_T(9), B_T(9) ! zero-based indices of the first A/B orbital in given symmetry
        integer :: A_C(9), B_C(9) ! zero-based indices of the first A/B continuum orbital in given symmetry
        integer :: iA_T, iB_T     ! current positions in arrays A_T, B_T
        integer :: iA_C, iB_C     ! current positions in arrays A_C, B_C
        logical :: last_step      ! transformation of the last index
        integer :: rank           ! rank of the current MPI task

        ! sanity check - the multiplication goes in sync along A and B storage (ie. A has to be transposed)
        if (nAr /= nBr) then
            call xermsg ('molecular_orbital_basis_obj', 'sparse_mmul', 'Non-conformant matrices passed as arguments.', 1, 1)
        end if

        ! determine if this is the last step (and column indices of both A and B correspond to molecular orbitals)
        last_step = triangle .and. pyramid >= 0

        ! scan the matrices to obtain important column pointers to speed up jumping through the columns
        call extract_orb_data(TA_orb, CA_orb, nAr, nAc, Ap, A_T, A_C)
        call extract_orb_data(TB_orb, CB_orb, nBr, nBc, Bp, B_T, B_C)

        ! loop through rows of A
        iA_T = 1; iA_C = 1; asym = 0; atype = 0; cpos = 0
        A_loop: do acol = 1, nAc

            ! starting position of this column in C
            Cp(acol) = cpos + 1

            ! first and last element of the column
            apos = Ap(acol)
            aend = Ap(acol + 1) - 1

            ! exit if at end of matrix
            if (apos == Ap(nAc + 1)) exit A_loop

            ! update current orbital symmetry / type
            do while (A_T(iA_T) <= apos .or. A_C(iA_C) <= apos)
                if (iA_C < iA_T) then
                    atype = 1 ; iA_C = iA_C + 1
                else
                    atype = 0 ; iA_T = iA_T + 1 ; asym = asym + 1
                end if
            end do

            ! skip CCCT and CCCC combinations
            if (skip2ec .and. blockcnt + atype >= 3) cycle A_loop

            ! in the last step calculate the combined [kli] orbital symmetry
            if (last_step) sym = abel_prod_tab(blocksym, asym)

            ! loop through elements of B
            iB_T = 1; iB_C = 1; bsym = 0; btype = 0
            B_loop: do bcol = 1, nBc

                ! first and last element of the column
                bpos = Bp(bcol)
                bend = Bp(bcol + 1) - 1

                ! exit if at end of matrix
                if (bpos == Bp(nBc + 1)) exit B_loop

                ! take into account (i.e. skip) the i <-> j, k <-> l and [ij] <-> [kl] symmetry
                if ((triangle .and. bcol > acol) .or. &
                    (pyramid >= 0 .and. (acol - 1) * nBc + (bcol - 1) > pyramid)) then
                    exit B_loop
                end if

                ! update current orbital symmetry and type
                do while (B_T(iB_T) <= bpos .or. B_C(iB_C) <= bpos)
                    if (iB_C < iB_T) then
                        btype = 1 ; iB_C = iB_C + 1
                    else
                        btype = 0 ; iB_T = iB_T + 1 ; bsym = bsym + 1
                    end if
                end do

                ! skip combinations CTCC
                if (skip2ec .and. blockcnt + atype + btype >= 3) cycle B_loop

                ! during the last step skip (non-contributing) orbitals of symmetries other than sym_i * sym_j * sym_k
                if (last_step .and. sym /= bsym) cycle B_loop

                ! redistribute over the 'kl' pairs of MO indices in the transformation [pq|rs] -> [pq|kl)
                if (mo_mpi_redistribution) then
                   rank = mod(bcol+(acol-1)*nBc, int(nprocs))
                   if (rank /= myrank) cycle B_loop
                endif

                ! calculate dot product of the two sparse columns
                x = 0
                acur = apos
                bcur = bpos
                do while (acur <= aend .and. bcur <= bend)
                    arow = Aj(acur)
                    brow = Bj(bcur, 1)
                    if      (arow < brow) then ; acur = acur + 1
                    else if (arow > brow) then ; bcur = bcur + 1
                    else
                        x = x + Av(acur) * Bv(Bv_beg-1 + bcur, 1)
                        acur = acur + 1
                        bcur = bcur + 1
                    end if
                end do

                ! store the result in the next free cell of C
                if (x /= 0) then
                    cpos = cpos + 1
                    Cj(cpos, 1) = bcol              ! set row index
                    Cv(cpos, 1) = x                 ! set value
                end if

            end do B_loop

        end do A_loop

        ! finalize the C column pointer array
        Cp(nAc + 1) = cpos + 1

    end subroutine transform_one_index


    !> \brief  Get special sparse matrix column pointers
    !> \author Jakub Benda
    !> \date   2018
    !>
    !> Constructs two sets of column pointers: One pointing to the first column corresponding to the first molecular
    !> orbital in every symmetry, the other pointing to the first continuum molecular orbital in every symmetry.
    !>
    !> \param T_orb  Absolute index of the first molecular orbital per symmetry.
    !> \param C_orb  Absolute index of the first continuum molecular orbital per symmetry.
    !> \param nAr    Number of rows in Ai.
    !> \param nAc    Number of columns in Ai (must be equal to the total number of molecular orbitals).
    !> \param An     Number of elements in Ai.
    !> \param A_T    On output, positions in Ai of the column corresponding to the first molecular orbital per symmetry.
    !> \param A_C    On output, positions in Ai of the column corresponding to the first continuum molecular orbital per symmetry.
    !>
    subroutine extract_orb_data(T_orb, C_orb, nAr, nAc, Ap, A_T, A_C)

        integer, allocatable :: Ap(:) !intent(in)
        integer, intent(in)  :: T_orb(9), C_orb(9), nAr, nAc
        integer, intent(out) :: A_T(9), A_C(9)

        integer :: orb, pos, iT, iC

        A_T = Ap(nAc + 1)  ! by default point beyond the storage
        A_C = Ap(nAc + 1)  ! by default point beyond the storage

        iT = 1
        iC = 1

        do orb = 1, nAc
            pos = Ap(orb)
            do while (orb >= T_orb(iT)) ; A_T(iT) = pos ; iT = iT + 1 ; end do
            do while (orb >= C_orb(iC)) ; A_C(iC) = pos ; iC = iC + 1 ; end do
        end do

    end subroutine extract_orb_data


    !> \brief   Twice cycle 4-tensor indices
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Transforms packed multi-indices of a 4-tensor, rotating them twice.
    !> This way, the same element will change its multi-index from [abcd] to [cdab].
    !>
    !> \param Rn   Number of values in Ri.
    !> \param Ri   Zero-based multiindices to transform; originally (((a-1)*nAO + (b-1))*nMO + (c-1))*nMO + (d-1).
    !> \param nAO  Maximal value of the two outer indices (a,b).
    !> \param nMO  Maximal value of the two inner indices (c,d).
    !>
    subroutine invert_indexing(Rn, Ri, nAO, nMO)

        integer, intent(inout) :: Ri(:,:)
        integer, intent(in)    :: Rn, nAO, nMO

        integer :: I, J, p, q, k, l

        write (level1, '("Transposing [AA|MM) -> (MM|AA]",/)')

        !$omp parallel do default(none) private(I, J, k, l, p, q) shared(Rn, Ri) firstprivate(nAO, nMO)
        do I = 1, Rn

            ! decode rectangular index
            J = Ri(I,1) ; l = mod(J, nMO)
            J = J / nMO ; k = mod(J, nMO)
            J = J / nMO ; q = mod(J, nAO)
            J = J / nAO ; p = J

            ! construct zero-based multi-index of a transposed element
            Ri(I,1) = ((k * nMO + l) * nAO + p) * nAO + q

        end do

    end subroutine invert_indexing


    !> \brief   Threaded alternative to heap_sort_int_float
    !> \authors J Benda
    !> \date    2025
    !>
    !> Splits the provided range equally among the available OpenMP threads and sorts each using `heap_sort_int_float`.
    !> If the remainder of integer division of the number of elements by the number of threads is "k", then the first
    !> "k" threads are given chunks one element longer than the other threads.
    !>
    !> Next, the subroutine passes the presorted chunks to `sort_intermediate_integrals`, which uses parallel merges
    !> to finalize the sorting.
    !>
    subroutine sort_final_integrals(Rn, Ri, Rv)

        use omp_lib,  only: omp_get_num_threads, omp_get_thread_num
        use sort_gbl, only: heap_sort_int_float

        integer,   intent(in)    :: Rn
        integer,   intent(inout) :: Ri(:, :)
        real(cfp), intent(inout) :: Rv(:, :)

        integer :: ithread, nthreads, a, b

        !$omp parallel
        !$omp single
        nthreads = min(Rn, int(omp_get_num_threads()))
        !$omp end single
        !$omp end parallel

        !$omp parallel num_threads(nthreads) default(none) shared(nthreads, Rn, Ri, Rv) private(ithread, a, b)
        ithread = omp_get_thread_num()
        a = Rn/nthreads*ithread + 1 + min(ithread, mod(Rn, nthreads))
        b = Rn/nthreads*(ithread + 1) + min(ithread + 1, mod(Rn, nthreads))
        call heap_sort_int_float(b - a + 1, Ri(a:b, 1), Rv(a:b, 1))
        !$omp end parallel

        call sort_intermediate_integrals(Rn, Ri, Rv, nthreads)

    end subroutine sort_final_integrals


    !> \brief   Sort the re-indexed integrals between transformation steps
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Sorts the integrals using a buffered parallel merge sort.
    !>
    !> The advantage of the merge sort is that the input array already consists of concatenated
    !> sorted sub-arrays by construction; there are nAO * nAO sorted segments, roughly of comparable size.
    !>
    !> The parallel merge sort has two stages: In the first stage, adjacent sorted sub-arrays are,
    !> in parallel, periodically collapsed (merged) together to form larger sorted segments, until
    !> the remaining number of (by now fairly large) segments becomes comparable to the number of threads.
    !> In the second stage, only one pair of segments is merged at a time (the one that will produce the smallest
    !> combined segment), but now using a parallel merge. This way, good parallel scaling is achieved
    !> in both stages, at least as far as memory bandwidth permits.
    !>
    !> \param Rn  Size of Ri, Rv.
    !> \param Ri  Integer array to sort.
    !> \param Rv  Real array to sort.
    !> \param n   Number of sorted sequences in Ri.
    !> \param max_memory_opt  Limit (in MiB) on combined memory used by this subroutine.
    !>
    subroutine sort_intermediate_integrals (Rn, Ri, Rv, n, max_memory_opt, scratch_prefix_opt)

        use omp_lib,  only: omp_get_wtime
        use sort_gbl, only: heap_sort_int_float, multiway_merge_float_int

        integer,   intent(in)    :: Rn, n
        integer,   intent(inout) :: Ri(:,:)
        real(cfp), intent(inout) :: Rv(:,:)

        real(cfp),        optional, intent(in) :: max_memory_opt
        character(len=*), optional, intent(in) :: scratch_prefix_opt

        logical, parameter :: force_heap_sort = .false.

        integer  :: i, j, nseg, nchunk
        real(wp) :: start_t, end_t, max_memory_MiB, used_memory_MiB

        integer, allocatable :: sizes(:)

        start_t = omp_get_wtime()

        ! optionally override the algorithm by a serial one (possibly useful for debugging as a reference)
        if (force_heap_sort) then
            write (level1, '("Forced serial heap sort")')
            call heap_sort_int_float(Rn, Ri(:, 1), Rv(:, 1))
            end_t = omp_get_wtime()
            write (level1, '("Time spent in sort_intermediate_integrals [s]: ",F25.15)') end_t - start_t
            return
        endif

        ! get the user-provided memory limit
        max_memory_MiB = 0
        if (present(max_memory_opt)) then
            max_memory_MiB = max_memory_opt
        end if

        ! prepare bookkeeping array for sorted chunk sizes
        allocate (sizes(n))
        sizes = 0

        ! loop over presorted segments and call the sorting routine every time the accumulated memory exceeds the buffer
        i = 1
        used_memory_MiB = 0
        nseg = 0
        nchunk = 0
        do j = 1, Rn
            ! factor 2 is due to the needed duplicity in the merge sort
            used_memory_MiB = used_memory_MiB + 2*(storage_size(Ri)/8 + storage_size(Rv)/8)/1024.**2

            ! avoid sorting within a presorted segment
            if (j < Rn) then
                if (Ri(j, 1) < Ri(j + 1, 1)) cycle
            end if

            ! update the number of presorted segments encountered so far
            nseg = nseg + 1

            ! check if the memory limit has been reached
            if (j == Rn .or. (max_memory_MiB > 0 .and. used_memory_MiB >= max_memory_MiB)) then
                call sort_intermediate_integrals_chunk(j - i + 1, Ri(i:j, :), Rv(i:j, :), nseg)
                i = j + 1
                nchunk = nchunk + 1
                sizes(nchunk) = j
                nseg = 0
                used_memory_MiB = 0
            end if
        end do
        if (nchunk > 1) then
            sizes(2:nchunk) = sizes(2:nchunk) - sizes(1:nchunk-1)
        end if

        ! merge all sorted chunks
        write (level2, '(a,i0)') 'Number of sorted segments to merge: ', nchunk
        write (level3, '(a,*(1x,i0))') 'Segment sizes:', sizes(1:nchunk)
        call multiway_merge_float_int(Rn, Rv(:, 1), Ri(:, 1), sizes(1:nchunk), scratch_prefix_opt)

        ! check that the sort went well
        !$omp parallel do
        do i = 2, Rn
            if (Ri(i - 1, 1) > Ri(i, 1)) then
                write (stdout, '(I0,",",I0,": ",I0," > ",I0)') i - 1, i, Ri(i - 1, 1), Ri(i, 1)
                call xermsg ('molecular_orbital_basis_obj', 'sort_intermediate_integrals', 'Sort failed.', 1, 1)
            end if
        end do

        end_t = omp_get_wtime()

        write (level1, '("Time spent in sort_intermediate_integrals [s]: ",F25.15)') end_t - start_t

    end subroutine sort_intermediate_integrals


    !> \brief   Sort a subset of intermediate integrals
    !> \authors J Benda
    !> \date    2025
    !>
    !> Called from `sort_intermediate_integrals` to sort a subset of the intermediate integrals, as governed
    !> by the user-specified maximal buffer size.
    !>
    subroutine sort_intermediate_integrals_chunk(Rn, Ri, Rv, n)

        use omp_lib, only: omp_get_num_threads, omp_get_wtime

        integer,   intent(in)    :: Rn, n
        integer,   intent(inout) :: Ri(:,:)
        real(cfp), intent(inout) :: Rv(:,:)

        integer, parameter :: k = 10  ! more or less arbitrary small number

        integer  :: ierr0, ierr1, ierr2, nthreads, i, j, nseg, a, b, c, m
        real(wp) :: start_t, end_t

        integer,   allocatable :: Wi(:,:), offsets(:)
        real(cfp), allocatable :: Wv(:,:)

        allocate (offsets(n + 1), stat = ierr0)
        allocate(Wi(size(Ri, 1), size(Ri, 2)), stat = ierr1)
        allocate(Wv(size(Rv, 1), size(Rv, 2)), stat = ierr2)

        ! Modified, Josh Forer 2026, omp master -> omp masked
        !$omp parallel
        !$omp masked
        nthreads = omp_get_num_threads()
        !$omp end masked
        !$omp end parallel

        write (level2, '("Using parallel merge sort with ",I0," threads")') nthreads

        ! get segment offsets
        nseg = 0
        do i = 1, Rn
            if (i == 1) then
                nseg = nseg + 1
                offsets(nseg) = i
            else if (Ri(i - 1, 1) > Ri(i, 1)) then
                nseg = nseg + 1
                offsets(nseg) = i
            end if
        end do
        offsets(nseg + 1) = Rn + 1

        ! reduce the number of segments to a more manageable size
        !$omp parallel default(none) firstprivate(Rn, nthreads) shared(nseg, offsets, Ri, Rv, Wi, Wv) private(i, a, b, c)
        do while (nseg > k * nthreads)
            !$omp do schedule(dynamic, 1)
            do i = 1, nseg / 2
                a = offsets(2*i - 1)
                b = offsets(2*i)
                c = offsets(2*i + 1)
                call parallel_merge_sorted_int_float(1, Ri   , Ri   , Wi   , &
                                                        Rv   , Rv   , Wv   , &
                                                        a,b-1, b,c-1, a,c-1, &
                                                     1)
                Ri(a:c-1, 1) = Wi(a:c-1, 1)
                Rv(a:c-1, 1) = Wv(a:c-1, 1)
            end do
            !$omp barrier     ! <-- should not be necessary (implicit to "end do"), but IFX 2024.1.0 with -O3 hangs if omitted
            !$omp single
            do i = 1, (nseg + 1) / 2
                offsets(i) = offsets(2*i - 1)
            end do
            offsets((nseg + 1) / 2 + 1) = offsets(nseg + 1)
            nseg = (nseg + 1) / 2
            !$omp end single
        end do
        !$omp end parallel

        ! merge segments further until there is just one of them (= the sorted array)
        do while (nseg > 1)
            ! find the smallest double segment
            i = 1
            m = offsets(3) - offsets(1)
            do j = 2, nseg - 1
                if (m > offsets(j + 2) - offsets(j)) then
                    m = offsets(j + 2) - offsets(j)
                    i = j
                end if
            end do
            a = offsets(i)
            b = offsets(i + 1)
            c = offsets(i + 2)
            ! perform the multi-threaded merge
            call parallel_merge_sorted_int_float(nthreads, Ri   , Ri   , Wi   , &
                                                           Rv   , Rv   , Wv   , &
                                                           a,b-1, b,c-1, a,c-1, &
                                                         1)
            ! copy merged segment back to source
            Ri(a:c-1, 1) = Wi(a:c-1, 1)
            Rv(a:c-1, 1) = Wv(a:c-1, 1)
            ! collapse offsets
            offsets(i + 1 : nseg) = offsets(i + 2 : nseg + 1)
            nseg = nseg - 1
        end do

    end subroutine sort_intermediate_integrals_chunk


    !> \brief  Merge sorted indexed arrays
    !> \author Jakub Benda
    !> \date   2018
    !>
    !> Merge two sorted integer arrays, and also corresponding real arrays. Uses the given number of threads.
    !>
    !> \param nthreads  Number of OpenMP threads to use.
    !> \param Ai  First integer array to merge.
    !> \param Bi  Second integer array to merge.
    !> \param Ci  Destination integer array (combined length of Ai, Bi).
    !> \param Af  First real array to merge.
    !> \param Bf  Second real array to merge.
    !> \param Cf  Destination real array (combined length of Af, Bf).
    !> \param As, Ae  Indices defining sections of the arrays Ai,Af to use: Ai(As:Ae), Af(As:Ae).
    !> \param Bs, Be  Indices defining sections of the arrays Bi,Bf to use: Bi(Bs:Be), Bf(Bs:Be).
    !> \param Cs, Ce  Indices defining sections of the arrays Ci,Cf to use: Ci(Cs:Ce), Cf(Cs:Ce).
    !>
    recursive subroutine parallel_merge_sorted_int_float (nthreads, Ai    , Bi    , Ci    , &
                                                                    Af    , Bf    , Cf    , &
                                                                    As, Ae, Bs, Be, Cs, Ce, &
                                                                    icol )

        integer,   intent(in)    :: nthreads, As, Ae, Bs, Be, Cs, Ce, icol
        integer,   intent(in)    :: Ai(:,:), Bi(:,:)
        integer,   intent(inout) :: Ci(:,:)
        real(cfp), intent(in)    :: Af(:,:), Bf(:,:)
        real(cfp), intent(inout) :: Cf(:,:)

        integer :: i, An, Bn, Cn, A_offsets(nthreads + 1), B_offsets(nthreads + 1), C_offsets(nthreads + 1)

        An = Ae-As+1 ! = size(Ai)
        Bn = Be-Bs+1 ! = size(Bi)
        Cn = Ce-Cs+1 ! = size(Ci)

        ! let A be larger
        if (An < Bn) then
            call parallel_merge_sorted_int_float(nthreads, Bi    , Ai    , Ci    , &
                                                           Bf    , Af    , Cf    , &
                                                           Bs, Be, As, Ae, Cs, Ce, &
                                                           icol )
            return
        end if

        ! partition A into "nthreads" segments
        do i = 1, nthreads + 1
            A_offsets(i) = (i - 1) * An / nthreads + 1
        end do

        ! bisect-find corresponding segments in B
        B_offsets(1) = Bs
        do i = 2, nthreads
            B_offsets(i) = lower_bound(Bi, icol, B_offsets(i - 1), Be, Ai(As - 1 + A_offsets(i), icol))
        end do
        B_offsets(nthreads + 1) = Be + 1

        ! here B_offsets must be relative wrt Bs
        B_offsets = B_offsets - Bs + 1

        ! determine segment sizes in C
        do i = 1, nthreads + 1
            C_offsets(i) = A_offsets(i) + B_offsets(i) - 1
        end do

        ! transform relative offsets to absolute offsets
        A_offsets = A_offsets + As - 1
        B_offsets = B_offsets + Bs - 1
        C_offsets = C_offsets + Cs - 1

        ! merge the segments independently
        !$omp parallel do num_threads(nthreads)
        do i = 1, nthreads
            call serial_merge_sorted_int_float (  &
                Ai                           , Bi                           , Ci                           ,  &
                Af                           , Bf                           , Cf                           ,  &
                A_offsets(i),A_offsets(i+1)-1, B_offsets(i),B_offsets(i+1)-1, C_offsets(i),C_offsets(i+1)-1,  &
                icol )
        end do

    end subroutine parallel_merge_sorted_int_float


    !> \brief  Merge sorted indexed arrays
    !> \author J Benda
    !> \date   2018 - 2020
    !>
    !> Merge two sorted integer arrays, and also corresponding real arrays.
    !>
    !> \param Ai  First integer array to merge.
    !> \param Bi  Second integer array to merge.
    !> \param Ci  Destination integer array (combined length of Ai, Bi).
    !> \param Af  First real array to merge.
    !> \param Bf  Second real array to merge.
    !> \param Cf  Destination real array (combined length of Af, Bf).
    !> \param As, Ae  Indices defining sections of the arrays Ai,Af to use: Ai(As:Ae), Af(As:Ae).
    !> \param Bs, Be  Indices defining sections of the arrays Bi,Bf to use: Bi(Bs:Be), Bf(Bs:Be).
    !> \param Cs, Ce  Indices defining sections of the arrays Ci,Cf to use: Ci(Cs:Ce), Cf(Cs:Ce).
    !>
    subroutine serial_merge_sorted_int_float (Ai    , Bi    , Ci    , &
                                              Af    , Bf    , Cf    , &
                                              As, Ae, Bs, Be, Cs, Ce, &
                                              icol)

        integer,   intent(in)    :: As, Ae, Bs, Be, Cs, Ce, icol
        integer,   intent(in)    :: Ai(:,:), Bi(:,:)
        integer,   intent(inout) :: Ci(:,:)
        real(cfp), intent(in)    :: Af(:,:), Bf(:,:)
        real(cfp), intent(inout) :: Cf(:,:)

        integer :: i, j, k, An, Bn, Cn
        logical :: copyA, copyB

        An = Ae-As+1 !size(Ai)
        Bn = Be-Bs+1 !size(Bi)
        Cn = Ce-Cs+1 !size(Ci)  ! == An + Bn

        i = 0
        j = 0

        do k = 0, Cn - 1
            ! are there still any elements to copy from arrays A or B?
            copyA = i < An
            copyB = j < Bn

            ! if both, use the smaller element
            if (copyA .and. copyB) then
                copyA = Ai(As + i, icol) < Bi(Bs + j, icol)
                copyB = .not. copyA
            end if

            ! now perform the copy (also copy the corresponding float) and shift the correct pointer
            if (copyA) then
                Ci(Cs + k, icol) = Ai(As + i, icol)
                Cf(Cs + k, icol) = Af(As + i, icol)
                i = i + 1
            else
                Ci(Cs + k, icol) = Bi(Bs + j, icol)
                Cf(Cs + k, icol) = Bf(Bs + j, icol)
                j = j + 1
            end if
        end do

    end subroutine serial_merge_sorted_int_float


    !> \brief  Index of first element that is not less
    !> \author Jakub Benda
    !> \date   2018
    !>
    !> Funny name of this function comes from the name of the standard C++ function "std::lower_bound".
    !>
    !> Returns index of the first element in a column of sorted array "A" (from "i" to "j") that is not less than
    !> the given value "v", or "j+1" if there is no such element. Takes advantage of the fact that "A"
    !> is sorted and uses interval halving.
    !>
    !> \param A    Sorted (ascending) integer array to bisect-search.
    !> \param icol Column of the array A
    !> \param i    Start of interval to search.
    !> \param j    End of interval to search.
    !> \param v    Value to compare with.
    !>
    integer function lower_bound (A, icol, i, j, v) result(right)

        integer, intent(in) :: icol, i, j, v
        integer, intent(in) :: A(:,:)

        integer :: left, mid

        left  = i
        right = j + 1

        do while (left /= right)
            mid = (left + right) / 2
            if (A(mid,icol) < v) then
                left = mid + 1
            else
                right = mid
            end if
        end do

    end function lower_bound


    !> \brief   Fix indexing of the integral storage
    !> \authors Jakub Benda
    !> \date    2018
    !>
    !> Convert the greedy (rectangular) working indexing scheme used in two_electron_integrals_sparse to the compact (triangular) one.
    !>
    !> \param Rn  Size of the intex array.
    !> \param Ri  Index array.
    !> \param nMO Maximal value of individual sub-indices (rectangular size).
    !>
    subroutine rect_index_to_tri_index(Rn, Ri, nMO, icol)

        integer, intent(in)    :: Rn, nMO, icol
        integer, intent(inout) :: Ri(:,:)

        integer :: i, n, a, b, c, d, u, v

        !$omp parallel do default(none) private(i, n, a, b, c, d, u, v) shared(Rn, Ri, icol) firstprivate(nMO)
        do i = 1, Rn

            ! decode rectangular index
            n = Ri(i, icol)   ; d = 1 + mod(n, nMO)
            n = n / nMO ; c = 1 + mod(n, nMO)
            n = n / nMO ; b = 1 + mod(n, nMO)
            n = n / nMO ; a = 1 + n

            ! construct triangular index
            u     = max(a,b) * (max(a,b) - 1) / 2 + min(a,b)
            v     = max(c,d) * (max(c,d) - 1) / 2 + min(c,d)
            Ri(i, icol) = max(u,v) * (max(u,v) - 1) / 2 + min(u,v)

        end do

    end subroutine rect_index_to_tri_index


    !> \brief   Final stage of calculation of the 2e integrals
    !> \authors Jakub Benda, Zdenek Masin
    !> \date    2018, 2020
    !>
    !> Writes integrals and indices to disk and releases allocated memory. This includes gathering of
    !> the transformed molecular integrals scattered over all MPI tasks. In the sparse integral transform
    !> the parallelization scheme is that both AO and MO integrals are scattered but each task calculates
    !> the full contribution to each assigned MO integral.
    !> Note that this type of finalization only applies to the sparse integral transform.
    !>
    subroutine finalize_two_electron_integrals_sparse(this, Rn, Ri, Rv, integral_storage, integral_options, column_descriptor)

        use const_gbl, only: two_p_sym_ints, ijkl_indices_header, i_min_sparse
        use sort_gbl,  only: sort_int_1d
        use omp_lib,   only: omp_get_wtime

        class(molecular_orbital_basis_obj)         :: this
        class(integral_options_obj), intent(in)    :: integral_options
        class(integral_storage_obj), intent(inout) :: integral_storage

        character(len=line_len), allocatable, intent(in) :: column_descriptor(:)

        real(kind=cfp), intent(inout), target :: Rv(:,:)
        integer,        intent(inout), target :: Ri(:,:)
        integer,        intent(in)            :: Rn

        type(p2d_array_obj), pointer :: mo_integrals
        type(p2d_array_obj), target  :: tmp_mo_integrals
        real(kind=cfp),      pointer :: backup_a(:,:)

        character(len=line_len) :: mo_header, ind_header

        integer :: I, J, k, u, v, a, b, c, d, current_pos, first_record, last_record, lunit, ierr, imaster, from
        integer,                allocatable :: indices(:,:), n_ints_per_task(:)
        integer(kind=mpiint),   allocatable :: displs(:), rcounts(:)
        real(kind=cfp), allocatable :: integrals(:,:)
        real(kind=cfp)              :: tolerance
        integer :: last_record_integrals, last_record_indices, first_record_integrals, first_record_indices
        integer :: current_pos_integrals, current_pos_indices, n_total_integrals, ibeg, iend, n_to_send, stride
        integer :: ind_start, ind_end, current_index, index_range, print_cnt
        integer, parameter :: no_blocks = 0
        real(kind=wp) :: start_t, end_t
        logical :: write_column_headers

!         open (11, form = 'formatted')
!         do I = 1, Rn
!             J = Ri(I,1)
!             u = int(0.5 + sqrt(2.0 * J)) ; v = J - u * (u - 1) / 2
!             a = int(0.5 + sqrt(2.0 * u)) ; b = u - a * (a - 1) / 2
!             c = int(0.5 + sqrt(2.0 * v)) ; d = v - c * (c - 1) / 2
!             write (11, '(I15,1x,I5,1x,I5,1x,I5,1x,I5,1x,E25.15)') J, a, b, c, d, Rv(i,1)
!         end do
!         close (11)

        ! in the distributed mode we have to first reduce the partial contributions to each integral
        if (nprocs > 1) then

           allocate(displs(nprocs),rcounts(nprocs),n_ints_per_task(nprocs),stat=ierr)
           if (ierr /= 0) then
              call mpi_xermsg ('molecular_orbital_basis_obj', 'finalize_two_electron_integrals_sparse', &
                               'Memory allocation 1 failed.', ierr, 1)
           endif
           displs = 0; rcounts = 0; n_ints_per_task = 0

           call mpi_mod_allgather(Rn, n_ints_per_task)

           n_total_integrals = sum(n_ints_per_task)

           do i=1,nprocs
              write(level2,'("Task ",i3," has ",i0," integrals")') i-1, n_ints_per_task(i)
           enddo

           write(level1,'("Total number of integrals ",i0)') n_total_integrals
           flush(stdout)

           if (n_total_integrals == 0) then
              call mpi_xermsg ('molecular_orbital_basis_obj', 'finalize_two_electron_integrals_sparse', &
                               'Number of non-zero integrals is zero...something is wrong.', 1, 1)
           endif

           write (level2, '("Writing headers on the integrals file...")')

           imaster    = master
           lunit      = integral_storage % integral_file % get_unit_no()
           mo_header  = integral_storage % contruct_header_string(this % get_basis_name(), two_p_sym_ints)
           ind_header = integral_storage % contruct_header_string(this % get_basis_name(), ijkl_indices_header)
           
           integral_storage % data_header % name = mo_header
           
           first_record = integral_storage % integral_file % start_record(mo_header)
           call integral_options % write(lunit, first_record, current_pos)

           if (myrank == master) then
              ! write the integral array dimensions
              j = size(Rv, 2)
              write(lunit,pos=current_pos) n_total_integrals, j
              write(lunit) column_descriptor
              inquire(lunit,pos=first_record_integrals)
           endif

           call mpi_mod_bcast(first_record_integrals, master)

           ! calculate the size of the final integrals record
           last_record_integrals = first_record_integrals + n_total_integrals*storage_unit_cfp + 1*storage_unit_int
           call integral_storage % integral_file % close_record(mo_header, first_record, last_record_integrals)

           first_record = integral_storage % integral_file % start_record(ind_header)

           if (myrank == master) then
              ! write the indices array dimensions. Note that the 'column_descriptor' is not part of the ijkl_indices record.
              j = size(Ri, 2)
              write(lunit,pos=first_record) n_total_integrals, j
              inquire(lunit,pos=first_record_indices)
           endif

           call mpi_mod_bcast(first_record_indices, master)

           ! calculate the size of the final indices record
           last_record_indices = first_record_indices + n_total_integrals*storage_unit_int
           call integral_storage % integral_file % close_record(ind_header, first_record, last_record_indices)

           write (level2, '("...done")')
                 
           write (level2, '(/,10X,"Gathering integrals from all tasks onto the master...")')
           flush(stdout)

           if (integral_options%max_ijrs_size < 0.0_cfp) then
              ! allocate as much memory as we ideally need
              I = n_total_integrals
           else
              ! maximum number of elements allowed in the temporary storage
              I = int(max(integral_options%max_ijrs_size*Mib/(max(cfp_bytes,storage_unit_int)*1.0_cfp),1.0_cfp))
   
              ! maximum number of elements in each of the integral and index arrays
              I = min(max(I/2,i_min_sparse), n_total_integrals)
           endif

           if (I.ge.max_data_count) then
               I = max_data_count - 1
               write(stdout,'("max_ijrs has been overridden to : ",i0)') I
           endif

           allocate(indices(I,1),integrals(I,1),stat=ierr)
           if (ierr /= 0) then
              call mpi_xermsg ('molecular_orbital_basis_obj', 'finalize_two_electron_integrals_sparse', &
                               'Memory allocation 2 failed.', ierr, 1)
           endif
           indices = 0; integrals = 0

           ! maximum number of integrals or indices we can receive from each task
           stride = max(I/nprocs,1)

           j = max(nint(n_total_integrals/(stride*1.0_cfp)),1)

           write(level2,'("Estimate of the total number of gather cycles needed: ",i5)') j
           write(level2,'("Maximum number of elements to be gathered in one cycle from each task: ",I0)') stride

           current_pos_integrals = first_record_integrals
           current_pos_indices = first_record_indices

           print_cnt = 0
           write_column_headers = .true.

           ibeg = 1
           iend = 0

           do

              start_t = omp_get_wtime()

              j = min(ibeg + stride-1, Rn)

              ! send to all other tasks the value of my last index that would fit into the buffer
              n_to_send = Ri(j, 1)
              if (j == iend) n_to_send = 0 !this task has sent all of its integrals before

              call mpi_mod_allgather(n_to_send, n_ints_per_task)
              call sort_int_1d(int(nprocs), n_ints_per_task)

              ! check whether all tasks have reached the end of their integral arrays in the last iteration
              if (all(n_ints_per_task == 0)) exit
              
              ! the smallest of the indices sent from the tasks determines the largest index that can be sent from any task
              do k = 1, nprocs
                 if (n_ints_per_task(k) > 0) then
                    current_index = n_ints_per_task(k)
                    exit
                 endif
              enddo

              ! adjust my index range to send so that it doesn't contain indices larger than current_index
              do i = j, ibeg, -1
                 if (Ri(i, 1) <= current_index) then
                    iend = i
                    exit
                 endif
              enddo

              n_to_send = max(iend-ibeg+1,0)
              call mpi_mod_allgather(int(n_to_send,mpiint),rcounts)

              displs = 0
              do j = 2, nprocs
                 displs(j) = displs(j-1) + rcounts(j-1)
              enddo
       
              call mpi_mod_gatherv(Ri(ibeg:iend, 1), indices(:, 1)  , rcounts, displs, int(master,mpiint))
              call mpi_mod_gatherv(Rv(ibeg:iend, 1), integrals(:, 1), rcounts, displs, int(master,mpiint))

              u = sum(rcounts)      
         
              if (myrank == master) then

                 I = 0
                 do j = 1, nprocs
                    if (rcounts(j) > 0) I = I + 1
                 enddo !j
      
                 call sort_intermediate_integrals (u, indices, integrals, I)
            
                 ! if requested print the non-zero integrals
                 if (integral_options % print_integrals) then
                    j = 1
                    if (write_column_headers) then
                       write(level2,'("Column ",i5," descriptor: ",a)') j, column_descriptor(j)
                       write_column_headers = .false.
                    endif
            
                    do i=1,u
                       print_cnt = print_cnt + 1
                       if (integrals(i,j) /= 0.0_cfp) write(level2,'(i0,1x,e25.15)') print_cnt, integrals(i,j)
                    enddo !i
                 endif
            
                 write(lunit,pos=current_pos_integrals) integrals(1:u, 1)
                 inquire(lunit,pos=current_pos_integrals)
         
                 write(lunit,pos=current_pos_indices) indices(1:u, 1)
                 inquire(lunit,pos=current_pos_indices)
         
              endif

              ! prepare for next iteration
              ibeg = iend + 1

              end_t = omp_get_wtime()

              write(level2,'("Gather iteration for ",I0," integrals took ",f25.15," s.")') u, end_t-start_t
              flush(stdout)

           enddo

           if (myrank == master) then

              ! the last piece of information is the number of blocks (= 0)
              write(lunit,pos=current_pos_integrals) no_blocks
              inquire(lunit,pos=current_pos_integrals)

              if (current_pos_integrals /= last_record_integrals) then
                 print *, current_pos_integrals, last_record_integrals
                 call mpi_xermsg ('molecular_orbital_basis_obj', 'finalize_two_electron_integrals_sparse', &
                                  'Unexpected end of integrals record.', 2, 1)
              endif
   
              if (current_pos_indices /= last_record_indices) then
                 print *, current_pos_indices, last_record_indices
                 call mpi_xermsg ('molecular_orbital_basis_obj', 'finalize_two_electron_integrals_sparse', &
                                  'Unexpected end of indices record.', 3, 1)
              endif
           endif

           write (level2, '("...done")')

        else !no MPI: no gather is needed

           mo_integrals => tmp_mo_integrals
   
           ! initialize integral storage
           !todo this allocation is only necessary in serial mode: we're allocating the molecular integrals array twice!!!
           !     in the parallel mode I can write out directly the Rv array.
           if (mo_integrals % init(Rn, 1, 0, column_descriptor) /= 0) then
                   call xermsg ('molecular_basis_mod', 'finalize_two_electron_integrals_sparse', &
                                'Molecular integral storage initialization failed.', 1, 1)
           end if
   
           ! push the sorted integrals and indices to the storage
           backup_a => mo_integrals % a; mo_integrals % a => Rv
           this % ijkl_indices => Ri
           this % ind_ijkl_integral = Rn
   
           flush(stdout)

           ! if requested print the non-zero integrals
           if (integral_options % print_integrals) then
              !call mo_integrals % print(.true.,101)
              call mo_integrals % print(.true.)
           endif

           write (level2, '("Saving integrals to disk...")')
   
           imaster    = master
           lunit      = integral_storage % integral_file % get_unit_no()
           mo_header  = integral_storage % contruct_header_string(this % get_basis_name(), two_p_sym_ints)
           ind_header = integral_storage % contruct_header_string(this % get_basis_name(), ijkl_indices_header)
   
           integral_storage % data_header % name = mo_header
   
           first_record = integral_storage % integral_file % start_record(mo_header)
           call integral_options % write(lunit, first_record, current_pos)
           call mo_integrals % write(lunit, current_pos, last_record, imaster)
           call integral_storage % integral_file % close_record(mo_header, first_record, last_record)
   
           first_record = integral_storage % integral_file % start_record(ind_header)
           call this % write_ijkl_indices(lunit, first_record, last_record)
           call integral_storage % integral_file % close_record(ind_header, first_record, last_record)
           this % ijkl_indices => null()
   
           write (level2, '("...done")')

           ! release molecular integrals
           mo_integrals % a => backup_a
           if (mo_integrals % final() /= 0) then
               call xermsg ('molecular_orbital_basis_obj', 'finalize_two_electron_integrals_sparse', &
                            'Deallocation of the temporary integral array failed.', 5, 1)
           end if

        endif !(nprocs > 1)

    end subroutine finalize_two_electron_integrals_sparse


   subroutine two_electron_integrals(this,integral_storage,integral_options)
      use const_gbl, only: abel_prod_tab, Mib, cache_line_size, two_el_ints, two_p_sym_ints, ijkl_indices_header
      use mpi_gbl
      use omp_lib
      use parallel_arrays_gbl
      use special_functions_gbl, only: ipair, unpack_pq
      use sort_gbl, only: sort_int_float, heap_sort_int_float
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(integral_options_obj), intent(in) :: integral_options
      class(integral_storage_obj), intent(inout) :: integral_storage

      type(integral_storage_obj) :: ao_integrals_disk
      type(integral_options_obj) :: ao_int_opt
      integer :: p, q, pq, r, s, rs, i, j, ij, lunit, no_int, err, first_record, last_record, no_ao, int_type
      integer :: no_pairs, no_tot, my_start, my_end, block, n_threads, thread_id
      integer :: sym_i, sym_j, orb_i, orb_j, orb_j_it, d1, d2, current_pos
      integer :: rs_start, rs_end, rs_block, no_rs_blocks, rs_block_size, last_tgt_ao, n_ints, no_blocks_ao
      integer, allocatable :: rs_ind(:,:), ij_orbital_range(:,:), last_tgt_mo(:), ind_ijkl_integral(:), ijkl_indices(:,:), &
                              ij_offset(:), n_integrals(:)
      integer(kind=1), allocatable :: ij_type(:)
      type(p2d_array_obj), target :: integral_src, integral_tgt !we really need two of these in case disk-to-disk AO-MO run is required
      type(p2d_array_obj), pointer :: ao_integrals, mo_integrals
      real(kind=cfp), allocatable :: cf(:,:), cf_t(:,:), iqrs(:,:), ijrs(:,:), ijks_t(:,:), hlp(:), ijkl_integrals(:,:,:), &
                                     cf_t_non_zero(:,:)
      integer, allocatable :: mo_indices(:,:), n_cf_t_non_zero(:)
      real(kind=wp) :: trmo_start, trmo_end, start_t, end_t, total_trmo, total_ijrs, total_ijkl, val
      logical :: i_is_continuum, j_is_continuum
      logical :: ao_is_local = .true.
      character(len=line_len), allocatable :: column_descriptor(:)
      character(len=line_len) :: ao_header, mo_header, ind_header

      integer, parameter :: no_blocks = 0, padding = 2*cache_line_size/cfp_bytes
      logical, parameter :: mo_is_local = .true. !same meaning as no_blocks = 0

         start_t = omp_get_wtime()
 
         call mpi_mod_barrier(err)
 
         write(level3,'("--------->","molecular_orbital_basis_obj:two_electron_integrals")')
 
         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', 'The basis set has not been initialized.',1, 1)
         end if

         if (this % number_of_functions == 0) then
            call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                         'Number of molecular orbitals on input is zero.', 2, 1)
         end if

         if (.not.associated(this%ao_integral_storage) .or. .not.associated(this%ao_basis)) then
            call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                         'On input at least one of this%ao_integral_storage, this%ao_basis have not been associated.', 4, 1)
         endif

         !Header for the AO integrals that we're looking for
         ao_header = this%ao_integral_storage%contruct_header_string(this%ao_basis%get_basis_name(),two_el_ints)

         !Header for the MO integrals that will be calculated
         mo_header = integral_storage%contruct_header_string(this%get_basis_name(),two_p_sym_ints)

         !Header for the MO integrals that will be calculated
         ind_header = integral_storage%contruct_header_string(this%get_basis_name(),ijkl_indices_header)

         ! In this section we associate ao_integrals which contains the input AO integrals with the appropriate source.
         ! In case the AO integrals are stored in memory then we point directly to the array holding them.
         ! If the AO integrals are on disk then we load them into the local array 'integral' and set the pointer ao_integrals to that.
         ! At the moment 1p integral transform using SHARED input AO integrals is not supported.
         if (this%ao_integral_storage%in_memory()) then
            ao_integrals => this%ao_integral_storage%integral
            if (this%ao_integral_storage%data_header%name .ne. ao_header) then
               call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                            'The AO integrals on input are not compatible with the AO basis set for the MO basis set.', 5, 1)
            endif
            write(level2,'("AO integrals with header: ",a)') ao_header
            ao_is_local = .not.(ao_integrals%have_offsets()) !find out if the AO integrals in the memory are shared or local
         endif

         !load the AO integrals into memory as LOCAL arrays (if they are stored on disk)
         if (this%ao_integral_storage%on_disk()) then
            write(level2,'("Loading AO integrals from the disk...")')

            err = ao_integrals_disk%init(memory=integral_src)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', 'Memory allocation 3 failed.', err, 1)
            end if

            call ao_integrals_disk%read(this%ao_integral_storage,ao_header,ao_int_opt,ao_is_local)
            ao_integrals => ao_integrals_disk%integral !this points to the local array 'integral_src'
            write(level2,'("AO integrals with header: ",a)') ao_header
            write(level2,'("...done")')
         endif

         !BEYOND THIS POINT ao_integrals POINTS TO AN ARRAY CONTAINING THE AO INTEGRALS TO BE TRANSFORMED

         call ao_integrals%get_array_dim(d1,d2,no_blocks_ao) !This gives the dimensions of ao_integrals%a(1:d1,1:d2)
         call ao_integrals%get_column_descriptor(column_descriptor)

         write(level2,'("On input there is ",i0," types of AO integrals")') d2
         write(level2,'("Total number of AO integrals of each type: ",i0)') d1
         if (.not.(ao_is_local)) then
            write(level2,'("AO integrals are distributed over all tasks.")')
         else
            write(level2,'("Every process keeps the full set of AO integrals.")')
         endif

         ! Note that instead of loading the AO integrals we can calculate them now since we have the pointer to the AO basis set
         ! (this%ao_basis) and the AO integral routine (this%ao_integrals)...

         !Assign a pointer to the requested array given by the user.
         if (integral_storage%in_memory()) then
            integral_storage%data_header%name = mo_header

            mo_integrals => integral_storage%integral
         endif

         !If we request the output to be stored on disk then start a new record on the data file that will contain the integrals.
         !Assign a pointer to a temporary array.
         if (integral_storage%on_disk()) then
            !temporary storage for the integrals
            mo_integrals => integral_tgt

            lunit = integral_storage%integral_file%get_unit_no()             !unit that is associated to the file opened
            first_record = integral_storage%integral_file%start_record(mo_header) !position, within the data file, of the first record available for the integral data
         endif

         !BEYOND THIS POINT mo_integrals POINTS TO AN ARRAY CONTAINING THE TRANSFORMED MO INTEGRALS

         ! 2-PARTICLE INTEGRAL TRANSFORM STARTS HERE: the AO integrals are accessed through the pointer ao_integrals;
         ! the MO integrals are accessed through the pointer mo_integrals
!
!--------The orbitals may have changed so recalculate the so2mo,mo2so indices.
!
         call this%determine_auxiliary_indices

         no_ao = this%ao_basis%number_of_functions !total number of AOs

         !todo The ijrs array can be made much smaller once symmetry for the AOs has been implemented
         rs = no_ao*(no_ao+1)/2
         no_pairs = rs

         ij = this%number_of_functions !total number of MOs
         ij = ij*(ij+1)/2

         allocate(cf_t(this%number_of_functions,no_ao),rs_ind(1:2,rs),ij_type(ij),last_tgt_mo(this%no_irr),stat=err)
         if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 4 failed.', err, 1)

         !Copy the orbital coefficients to one array: this relies on the fact that the molecular orbitals are indexed symmetry by symmetry.
         call this%get_orbital_coefficient_matrix(cf)

         !Transposing the MO coefficient matrix and the other matrices speeds up the computation significantly.
         cf_t = transpose(cf)

         !Get index of the last AO representing the target electrons: this procedure assumes all target functions preceed the continuum ones.
         !todo once symmetry for the AOs has been implemented I should get this for each symmetry, i.e. last_tgt_ao will be an array
         last_tgt_ao = 0
         do i=1,no_ao
            j = this%ao_basis%indices_to_shells(1,i) !shell index containing the i-th function
            if (this%ao_basis%shell_descriptor(3,j) .eq. 1) cycle !skip continuum functions
            last_tgt_ao = max(last_tgt_ao,i)
         enddo

         !Prepare the set of (rs)-indices that will be split among all threads
         rs = 0
         do r=1,no_ao
            do s=1,r
               rs = rs + 1
               rs_ind(1:2,rs) = (/r,s/)
            enddo !s
         enddo !r

         !Determine the type of the orbital pairs
         do i=1,this%no_irr
            !Find the last target orbital within each symmetry: this procedure assumes all target orbitals preceed the continuum ones.
            last_tgt_mo(i) = 0
            do j=1,this%orbital_data(i)%number_of_functions
               if (this%is_continuum(this%relative_to_absolute(j,i))) cycle
               last_tgt_mo(i) = max(last_tgt_mo(i),this%relative_to_absolute(j,i))
            enddo !j
         enddo

         ij = 0
         do i=1,this%number_of_functions
            p = this%absolute_to_relative(2,i)
            i_is_continuum = .false.
            if (i > last_tgt_mo(p) .and. i <= this % relative_to_absolute(this % orbital_data(p) % number_of_functions, p)) then
                i_is_continuum = .true.
            end if
            do j=1,i
               q = this%absolute_to_relative(2,j)
               j_is_continuum = .false.
               if (j > last_tgt_mo(q) .and. j <= this % relative_to_absolute(this % orbital_data(q) % number_of_functions, q)) then
                  j_is_continuum = .true.
               end if
               ij = ij + 1
               ij_type(ij) = 1 !TT
               if (i_is_continuum .neqv. j_is_continuum) then
                  ij_type(ij) = 2 !TC or CT
               elseif (i_is_continuum .and. j_is_continuum) then
                  ij_type(ij) = 3 !CC
               endif
            enddo !j
         enddo !i

         !Loop over all types of AO integrals to transform.
         do int_type=1,d2

            trmo_start = omp_get_wtime()

            write(level1,'("Transforming AO integral type ",i0," ...")') int_type

            !We have to do this since ijrs and ij_orbital_range are being allocated below (i.e. inside the int_type loop).
            if (allocated(ijrs)) deallocate(ijrs)
            if (allocated(ij_orbital_range)) deallocate(ij_orbital_range)
            no_int = 0

            ! Note: This parallel section needs to use DEFAULT(SHARED) to allow work with polymorphic objects with gfortran.
            !$OMP PARALLEL &
            !$OMP & DEFAULT(SHARED) &
            !$OMP & PRIVATE(rs,r,s,p,q,ij,i,j,pq,block,n_threads,thread_id,err,rs_block_size,rs_start,rs_end,rs_block, &
            !$OMP &         no_rs_blocks,ij_orbital_range,sym_i,sym_j,orb_i,orb_j,orb_j_it,my_start,my_end,iqrs,ijks_t,hlp, &
            !$OMP &         ij_offset,val,n_ints,no_tot) &
            !$OMP & SHARED(no_pairs,ao_is_local,no_ao,this,int_type,cf,cf_t,ijrs,ao_integrals,rs_ind,integral_options, &
            !$OMP &        mo_integrals,total_trmo,trmo_start,trmo_end,total_ijrs,total_ijkl,level2,last_tgt_ao,ij_type, &
            !$OMP &        ijkl_integrals,ind_ijkl_integral,ijkl_indices,n_integrals,nprocs,cf_t_non_zero, &
            !$OMP &        mo_indices,n_cf_t_non_zero,myrank,column_descriptor,d2) &
            !$OMP & REDUCTION(+:no_int)
            n_threads = omp_get_num_threads()
            thread_id = omp_get_thread_num()

            p = this%no_irr*(this%no_irr+1)/2  !total number of symmetry pairs

            !The ij_orbital_range array is used in the second part of the integral transform to split the (i,j) indices in each symmetry combination among the threads.
            !In the buffers allocated below (iqrs,ijrs,ijks_t) the letters q,r,s stand for atomic orbitals and the letters i,j,k stand for the molecular orbitals.
            !It proved necessary to perform the allocations of the arrays iqrs,ijks_t,hlp explicitly here rather than allocating them before the parallel region and making
            !them PRIVATE. For large calculations (large AO basis) these arrays were causing SEGFAULTs even with unlimited stack set on the compute nodes (the stack can never be
            !unlimited in practice so there is always a maximum that can be breached).
            ij = this%number_of_functions*(this%number_of_functions+1)/2
            allocate(ij_orbital_range(4,p),iqrs(this%number_of_functions,no_ao),&
                        ijks_t(no_ao,this%number_of_functions),hlp(no_ao),ij_offset(ij),stat=err)
            if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 5 failed.', err, 1)
            ij_orbital_range = 0

            !split the unique (i,j) orbital indices in each combination of symmetries among the threads
            do sym_i=1,this%no_irr
               do sym_j=1,sym_i

                  p = sym_i*(sym_i-1)/2 + sym_j !the index of the symmetry block for the pair of symmetries corresponding to the pair (sym_i,sym_j), sym_i .g.e sym_j

                  !calculate the total number of orbital pairs for this combination of symmetries
                  if (sym_i .eq. sym_j) then
                     ij = this%orbital_data(sym_i)%number_of_functions*(this%orbital_data(sym_i)%number_of_functions+1)/2
                  else
                     ij = this%orbital_data(sym_i)%number_of_functions*this%orbital_data(sym_j)%number_of_functions
                  endif

                  ij_orbital_range(1,p) = 1; ij_orbital_range(3,p) = -1
                  if (ij .eq. 0) then 
                     cycle
                  endif

                  if (ij .le. n_threads) then
                     my_start = thread_id+1
                     my_end = thread_id+1
                  else
                     block = ceiling(ij/(n_threads*1.0))
                     my_start = block*thread_id+1 !ij-index for my first orbital pair
                     my_end = min(ij,block*(thread_id+1)) !ij-index for my last orbital pair
                  endif

                  !since my_end is always at most equal to ij we must make sure that if this thread cannot encounter the 'if (ij .eq. my_end)' below 
                  !in case it is not supposed to process any orbitals for this pair of symmetries.
                  if (my_start > my_end) cycle

                  ij = 0
                  do orb_i=1,this%orbital_data(sym_i)%number_of_functions   !over all MOs in symmetry sym_i
                     if (sym_i .eq. sym_j) then
                        orb_j_it = orb_i          !both orbitals are from the same symmetry and therefore the loop over the second orbital must be only from orb_j_it = 1 to orb_j_it = orb_i.
                     else
                        orb_j_it = this%orbital_data(sym_j)%number_of_functions !both orbitals come from different symmetries so the loop over the second loop must be over all orbitals in that symmetry.
                     endif

                     do orb_j=1,orb_j_it
                        ij = ij + 1

                        !Determine the first and the last pair of orbital indices that I'll be transforming, for each pair of symmetries.
                        if (ij .eq. my_start) then
                           ij_orbital_range(1,p) = orb_i
                           ij_orbital_range(2,p) = orb_j
                        endif
                        if (ij .eq. my_end) then
                           ij_orbital_range(3,p) = orb_i
                           ij_orbital_range(4,p) = orb_j
                        endif
                        if (ij > my_end) then
                           exit
                        endif
                     enddo !orb_j
                  enddo !orb_i
               enddo !sym_j
            enddo !sym_i

            ij = this%number_of_functions*(this%number_of_functions+1)/2 !number of (i,j) combinations of orbitals

            if (integral_options%max_ijrs_size > 0.0_cfp) then
               !Limit the size of the ijrs array to the maximum size (in Mib) supplied by the user.
               i = integral_options%max_ijrs_size*Mib/(cfp_bytes*1.0_cfp) !total number of elements allowed in the ijrs array
               rs_block_size = min(i/(ij+padding),no_pairs)
               !$OMP SINGLE
               write(level1,'("Size of the ijrs array limited to (MiB): ",f0.3)') integral_options%max_ijrs_size
               !$OMP END SINGLE
            else
               no_tot = this%block_offset(size(this%block_offset))  !the sum of integrals in each of the unique symmetry blocks
               !Set-up the size of the (rs)-block (i.e. the number of (r,s)-indices in one block) relatively to the size of the array of the transformed integrals.
               rs_block_size = min(no_tot/ij,no_pairs)
            endif

            no_rs_blocks = no_pairs/rs_block_size
            if (no_rs_blocks*rs_block_size < no_pairs) no_rs_blocks = no_rs_blocks + 1

            call generate_ij_offset(this,ij_type,ij_orbital_range,integral_options%two_p_continuum,ij_offset,n_ints)

            !$OMP SINGLE
            allocate(n_integrals(n_threads),stat=err)
            if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 6a failed.', err, 1)
            n_integrals = 0
            !$OMP END SINGLE

            !Total number of ijkl integrals generated by this thread:
            n_integrals(thread_id+1) = n_ints
            !$OMP BARRIER

            !$OMP SINGLE
            write(level1,'(/,"Total number of (rs)-indices: ",i0)') no_pairs
            write(level1,'("Number of (rs)-index blocks: ",i0,", number of (rs)-indices in block: ",i0)') &
                no_rs_blocks, rs_block_size
            write(level1,'(/,"Size of the ijrs buffer to be allocated [Mib]: ",f25.15)') &
                (rs_block_size*(ij+padding)*cfp_bytes*1.0_cfp)/Mib

            !Allocate the output arrays if we request the output to be stored in memory or on the disk.
            !Note that we enforce the mo_integrals array to be LOCAL on each process (no_blocks=0), i.e. every process gets
            !the full integrals array. The present transformation algorithm doesn't 
            !implement splitting of the transformed (ij|kl) integrals among the processes. We allocate space for a non-indexed
            !(that is purely local) array with d2 columns and no_tot rows.
            !The columns in the mo_integrals array correspond to the types of AO integrals we have on input.
            if (int_type .eq. 1) then
               no_tot = sum(n_integrals)
               write(level2,'("Total number of MO integrals of each type that will be obtained: ",i0)') no_tot
               err = mo_integrals%init(no_tot,d2,no_blocks,column_descriptor)
               if (err /= 0) then
                  call xermsg ('molecular_basis_mod', 'two_electron_integrals', &
                               'Array initialization 2 failed; see p2d_array_obj%init.', err, 1)
               end if
               !Allocate the array of indices for the transformed integrals.
               if (associated(this%ijkl_indices)) deallocate(this%ijkl_indices)
               allocate(this%ijkl_indices(no_tot,d2),stat=err)
               if (err /= 0) then
                  call xermsg ('molecular_basis_mod', 'two_electron_integrals', &
                               'Array initialization 3 failed; see p2d_array_obj%init.', err, 1)
               end if
               this%ijkl_indices = 0
               this%ind_ijkl_integral = 0 !index of the last integral stored in mo_integrals%a
            endif

            i = maxval(n_integrals)
            !todo The ijrs array can be made much smaller once symmetry for the AOs has been implemented
            allocate(ijrs(ij+padding,rs_block_size),ijkl_integrals(i,n_threads,int_type:int_type), &
                     ind_ijkl_integral(n_threads),ijkl_indices(i,n_threads),stat=err)
            if (err .ne. 0) call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Memory allocation 6 failed.', err, 1)
            ijrs(:,:) = 0.0_cfp
            total_ijrs = 0.0_cfp
            total_ijkl = 0.0_cfp
            ijkl_integrals = 0.0_cfp
            ind_ijkl_integral = 0
            ijkl_indices = 0
            val = 2*i*n_threads*cfp_bytes/(Mib*1.0_cfp)
            write(level2,'("Memory storage for temporary arrays ijkl_integrals,ijkl_indices (MiB): ",f8.3)') val
            write(level2,'("Memory has been successfuly allocated.",/)')

            call extract_non_zero_cf(cf_t,cf_t_non_zero,mo_indices,n_cf_t_non_zero)
            !$OMP END SINGLE

            !Total number of ijkl integrals generated by this thread:
            ind_ijkl_integral(thread_id+1) = n_integrals(thread_id+1)

            !We divide the whole set of (r,s)-indices into blocks and for each block we transform the first two indices (pq|rs)->(ij|rs) and then calculate their contribution 
            !to the final set of transformed integrals: (ij|rs) -> (ij|kl).
            do rs_block=1,no_rs_blocks

               !Determine the range of rs-indices in the current block
               rs_start = rs_block_size*(rs_block-1)+1
               rs_end = min(rs_block_size*rs_block,no_pairs)

               !$OMP SINGLE
               trmo_start = omp_get_wtime()
               !$OMP END SINGLE

               !I. THE FIRST PART OF THE INTEGRAL TRANSFORM:
               !Transform the first two indices: (pq|rs) -> (ij|rs) for all (i,j)-indices (pairs of orbitals) and the current
               !block of (rs)-indices (rs_start:rs_end). Here we loop over the (rs) indices of the AO integrals (pq|rs) and
               !obtain the partially transformed integrals (ij|rs), where i,j are the molecular orbitals. It is only this part
               !of the integral transformation algorithm where need to accommodate the fact that the AO integrals may be incomplete
               !if they were split among the processes (ao_is_local).
               if (ao_is_local) then !every process keeps all AO integrals
                  call omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local ( &
                            ao_integrals,int_type,rs_start,rs_end,ij_type,iqrs,hlp,this,cf,cf_t,ijrs,no_ao,last_tgt_ao, &
                            no_int,no_pairs,rs_ind,integral_options%two_p_continuum,integral_options%aoints_index_scheme &
                  )
               else !in case the AO integrals are scattered randomly among the processes
                  call omp_two_p_transform_pqrs_block_to_ijrs_AO_is_not_local ( &
                            ao_integrals,int_type,rs_start,rs_end,ij_type,iqrs,hlp,this,cf,cf_t,ijrs,no_ao,last_tgt_ao, &
                            no_int,no_pairs,rs_ind,integral_options%two_p_continuum &
                  )
               endif
               !Remember that we don't need BARRIER here only as long as the loop over rs-indices inside the preceeding routines
               !is split using OMP DO since that implies synchronization.

               !Putting the barrier here ensures accurate timing is reported below: we want the timing for the slowest thread to be reported.
               !$OMP BARRIER

               !$OMP SINGLE
               trmo_end = omp_get_wtime()
               write(level1,'("(rs)-index block: ",i0)') rs_block
               write(level1,'("(pq|rs) -> (ij|rs) for rs-block: ",i0,"->",i0," done in ",f8.3," [s]")') &
                    rs_start,rs_end, trmo_end-trmo_start
               total_ijrs = total_ijrs + trmo_end-trmo_start
               trmo_start = omp_get_wtime()
               !$OMP END SINGLE

               !II. THE SECOND PART OF THE INTEGRAL TRANSFORM:
               !Transform the last two indices: (ij|rs) -> (ij|kl)
               !Note that we use as the ijks array the array iqrs since they have the same size.
               call omp_two_p_transform_ijrs_block_to_ijkl ( &
                        ijrs,iqrs,ijks_t,cf,cf_t_non_zero,mo_indices,n_cf_t_non_zero,no_ao,this,rs_ind,ij_type,ao_is_local, &
                        integral_options%tol,ijkl_integrals,ijkl_indices,int_type,rs_start,rs_end,no_pairs,ij_orbital_range, &
                        integral_options%two_p_continuum,thread_id,ij_offset &
               )

               !Putting the barrier here ensures accurate timing is reported below: we want the timing for the slowest thread to be reported.
               !$OMP BARRIER

               !$OMP SINGLE
               trmo_end = omp_get_wtime()
               write(level1,'("(ij|rs) -> (ij|kl) contribution done in ",f8.3," [s]")') trmo_end-trmo_start
               total_ijkl = total_ijkl + trmo_end-trmo_start
               !$OMP END SINGLE

            enddo !rs_block

            !$OMP SINGLE
            deallocate(ijrs)
            !$OMP END SINGLE

            deallocate(ijks_t,iqrs,hlp)

            !$OMP SINGLE
            !todo put the integrals into blocks using ij_offset but generated
            !for the global array not for individual threads.
            this%ind_ijkl_integral = 0
            do i=1,n_threads
               do j=1,ind_ijkl_integral(i)
                  !At this stage we cannot skip integrals when more processes are used since the integral array must be reduced first.
                  if (nprocs .eq. 1) then
                     if (ijkl_integrals(j,i,int_type) .eq. 0.0_cfp) cycle
                  endif
                  !WARNING: if saving only the non-zero integrals any indexing
                  !making use of ij_offset becomes wrong!
                  this%ind_ijkl_integral = this%ind_ijkl_integral + 1
                  if (this%ind_ijkl_integral > size(mo_integrals%a,1)) then
                     print *,this%ind_ijkl_integral,size(mo_integrals%a,1)
                     stop "error 2:buffer too small"
                  endif
                  mo_integrals%a(this%ind_ijkl_integral,int_type) = ijkl_integrals(j,i,int_type)
                  !WARNING: note that the present implementation may not work correctly if int_type > 1 since the this%ind_ijkl_integral should depend on int_type too!
                  this%ijkl_indices(this%ind_ijkl_integral,int_type) = ijkl_indices(j,i)
!                  print *,'ind',this%ind_ijkl_integral,this%ijkl_indices(this%ind_ijkl_integral,int_type)
               enddo !j
            enddo !i

            deallocate(ind_ijkl_integral,ijkl_integrals,ijkl_indices,ij_offset)
            !$OMP END SINGLE

            !$OMP END PARALLEL

            write(level2,'(/,"Sorting the final array of integrals...")')
            !save the unsorted data so we can recover it later if needed
            !open(newunit=i,file='unsorted_integrals',status='replace',form='unformatted')
            !write(i) this%ind_ijkl_integral
            !write(i) this%ijkl_indices(1:this%ind_ijkl_integral,int_type)
            !write(i) this%ind_ijkl_integral
            !write(i) mo_integrals%a(1:this%ind_ijkl_integral,int_type)
            !close(i)
            trmo_start = omp_get_wtime()
            !Heap sort is generally faster than quicksort for the arrays that
            !we're dealing with. Actually, this may be caused by the current quicksort implementation not being very good.
            call heap_sort_int_float(this%ind_ijkl_integral,this%ijkl_indices(:,int_type),mo_integrals%a(:,int_type))
            !call sort_int_float(this%ind_ijkl_integral,int_type,this%ijkl_indices,mo_integrals%a) !Quick sort
            trmo_end = omp_get_wtime()
            write(level2,'("...done and took [s]: ",f8.3)') trmo_end-trmo_start

!            do i=1,this%ind_ijkl_integral
!               !todo temp debug
!               write(50,'(2i15,e25.15)') i,this%ijkl_indices(i,int_type),mo_integrals%a(i,int_type)
!            enddo

            total_trmo = total_ijrs + total_ijkl
            write(level1,'(/,"Number of AO integrals used: ",i0)') no_int
            write(level1,'("Total time for (pq|rs) -> (ij|rs) ",f0.3," [s]")') total_ijrs
            write(level1,'("Total time for (ij|rs) -> (ij|kl) ",f0.3," [s]")') total_ijkl
            write(level1,'("(pq|rs) -> (ij|kl) done in ",f0.3," [s]")') total_trmo
            write(level1,'("Number of ijkl integrals: ",i15," (",f0.3," MiB)")') &
                this % ind_ijkl_integral, this % ind_ijkl_integral * cfp_bytes / (Mib * 1.0_cfp)

         enddo !int_type

         deallocate(cf_t)

         if (integral_options%keep_ao_integrals) then
            write(level2,'("The AO integrals array will be kept in memory.")')
         else
            !Get rid of the AO integrals now: in case of very large MO integrals array the MPI reduction below may need a lot of memory.
            err = ao_integrals%final()
            if (err /= 0) then
               call xermsg ('molecular_basis_mod', 'two_electron_integrals', 'Finalizing the ao_integrals array failed.', err, 1)
            end if
   
            !If the AO integrals were stored only in memory then we have to inform the user that we've destroyed them now.
            if (this % ao_integral_storage%in_memory()) then
               write(level2,'("WARNING: The AO integrals array supplied on input has been destroyed.")')
            end if
         endif

         !THIS IS AN IMPORTANT STEP: in case the AO integral array was shared (and MO integral array is local) then
         !mo_integrals%a contains only the CONTRIBUTION of each process to the final transformed 
         !integrals. Therefore we need to REUDCE (sum) the contributions to the MO integrals from each process.
         !The results are communicated to all processes so at the end each process keeps its own copy
         !of the full MO integrals array (i.e. this is compatible with the LOCAL attribute of the mo_integrals array).
         if (.not.(ao_is_local) .and. mo_is_local) then
            if (cfp .eq. ep1) then
               call mpi_xermsg ('molecular_basis_mod', 'two_electron_integrals', &
                                'In quad precision alltoall reduction not implemented. &
                                &Only master will have the final MO integrals.', 1, 0)
               call mo_integrals%reduce_a_to_master
            else
               call mo_integrals%reduce_a_local
            endif

            !Get rid of the transformed integrals smaller than the threshold. Note that in the serial run this has been done above.
            !todo this should only be done on the master process if cfp .eq. ep1
            write(level2,'("Deleting integrals smaller than the threshold...")')
            do int_type=1,d2
               !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i) SHARED(mo_integrals,int_type,integral_options)
               !$OMP DO
               do i=1,size(mo_integrals%a,1)
                  if (abs(mo_integrals%a(i,int_type)) < integral_options%tol) mo_integrals%a(i,int_type) = 0.0_cfp
               enddo
               !$OMP END DO
               !$OMP END PARALLEL
            enddo
            write(level2,'("...done")')
         endif

         !If requested print the non-zero integrals
         if (integral_options%print_integrals) then
            if (nprocs > 1) then
               !Only the master process has the final integrals
               if (myrank == master) then
                  call mo_integrals%print(.true.)
               else
                  write(level2,'("The molecular integrals are printed only by the master process, see log_file.0")')
               endif
            else
               call mo_integrals%print(.true.)
            endif
         endif

         !Dump all integrals to disk and close the record.
         if (integral_storage%on_disk()) then

            write(level2,'("Saving integrals to disk...")')

            !The first record are the integral options.
            call integral_options%write(lunit,first_record,current_pos)

            !The second record are the ordered integrals: only master writes them to disk.
            i = master
            call mo_integrals%write(lunit,current_pos,last_record,i)

            !Every process closes the record so that they all keep identical header information.
            call integral_storage%integral_file%close_record(mo_header,first_record,last_record)

            err = mo_integrals%final()

            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'two_electron_integrals', &
                             'Deallocation of the temporary integral array failed.', 5, 1)
            end if

            !Write out the ijkl indices
            first_record = integral_storage%integral_file%start_record(ind_header) !position, within the data file, of the first record available for the integral data

            call this%write_ijkl_indices(lunit,first_record,last_record)

            call integral_storage%integral_file%close_record(ind_header,first_record,last_record)

            write(level2,'("...done")')

         endif
 
         write(level3,'("<---------","molecular_orbital_basis_obj:two_electron_integrals")')
 
         call mpi_mod_barrier(err)

         end_t = omp_get_wtime()

         write(level1,'("Two_electron_integrals took [s]: ",f25.15)') end_t-start_t
 
   end subroutine two_electron_integrals

   function integral_index(this,integral_type,bf_indices,two_p_continuum)
      use const_gbl
      use special_functions_gbl, only: ipair
      implicit none
      class(molecular_orbital_basis_obj) :: this
      character(len=*), intent(in) :: integral_type
      integer, intent(in) :: bf_indices(:,:)
      logical, intent(in) :: two_p_continuum
      integer :: integral_index(size(bf_indices,2))

      integer :: ind,i,j,k,iAB, iCD, a,b,c,d
      logical :: found

         if (size(bf_indices,1) .eq. 2) then !1-electron integral index

            !todo this can be improved by swapping the do-loop with select case for
            !each integral type. Best to pack the indexing computation into an elementary routine.
            do k=1,size(bf_indices,2)
               i = maxval(bf_indices(1:2,k))
               j = minval(bf_indices(1:2,k))

               ind = i*(i-1)/2+j

               !At the moment all 1-electron atomic integrals are index in the
               !same way.
               select case (integral_type)

                  case (overlap_ints)
                     integral_index(k) = ind
                  case (kinetic_ints)
                     integral_index(k) = ind
                  case (property_ints)
                     integral_index(k) = ind
                  case (nuc_rep_att_ints)
                     integral_index(k) = ind
                  case (ecp_ints)
                     integral_index(k) = ind
                  case (one_elham)
                     integral_index(k) = ind
                  case (bbb_ints)
                     integral_index(k) = ind
                  case default
                  call xermsg ('molecular_orbital_basis_obj', 'integral_index', &
                               'Unrecognized one electron molecular integral type on input.', 1, 1)

               end select
            enddo !k

         elseif (size(bf_indices,1) .eq. 4) then !2-electron integral index

            do i=1,size(bf_indices,2)
               a = maxval(bf_indices(1:2,i))
               b = minval(bf_indices(1:2,i))
               c = maxval(bf_indices(3:4,i))
               d = minval(bf_indices(3:4,i))
               if (a >= c) then
                  iAB = ipair(a)+b
                  iCD = ipair(c)+d
               else
                  iAB = ipair(c)+d
                  iCD = ipair(a)+b
               endif
               if (iAB >= iCD) then
                   ind = ipair(iAB)+iCD
               else
                   ind = ipair(iCD)+iAB
               endif

               !Indexing relying on the sorted array of molecular integrals
               call bisect_index(ind,this%ijkl_indices,1,this%ind_ijkl_integral,j,found)
               
               if (found) then
                  integral_index(i) = j
               else
                  integral_index(i) = -1 !this integral is zero
               endif
            enddo !i
         else
            call xermsg ('molecular_orbital_basis_obj', 'integral_index', &
                         'The number of orbital indices on input must be either 2 or 4.', 3, 1)
         endif

   end function integral_index

   function get_basis_name(this)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      character(len=line_len) :: get_basis_name

         get_basis_name = "molecular_orbital_basis_obj"
 
   end function get_basis_name

   function get_shell_name(this,i)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      character(len=line_len) :: get_shell_name
      integer, intent(in) :: i

      type(orbital_data_obj) :: dummy

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_name', 'The basis set has not been initialized.', 1, 1)
         end if

         if (i <= 0 .or. i > this % number_of_shells) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_name', 'On input the value of i was out of range.', 2, 1)
         end if

         get_shell_name = dummy%name()

   end function get_shell_name

   subroutine get_shell_data(this,i,shell_data)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      class(shell_data_obj), intent(out) :: shell_data
      integer, intent(in) :: i

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_data', 'The basis set has not been initialized.', 1, 1)
         end if

         if (i <= 0 .or. i > this % number_of_shells) then
            call xermsg ('molecular_orbital_basis_obj', 'get_shell_data', 'On input the value of i was out of range.', 2, 1)
         end if
         
         select type (shell => shell_data)
            type is (orbital_data_obj)
               shell = this%orbital_data(i)
            class default
               call xermsg ('molecular_orbital_basis_obj', 'get_shell_data', &
                            'On input shell_data must be of orbital_data_obj type.', 3, 1)
         end select

   end subroutine get_shell_data

   subroutine get_all_orbital_sets(this,orbital_sets,number_of_orbital_sets)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(out) :: number_of_orbital_sets
      type(orbital_data_obj), allocatable :: orbital_sets(:)

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_all_orbital_sets', 'The basis set has not been initialized.', 1, 1)
         end if

         number_of_orbital_sets = 0

         !todo this must involve transfering all orbital data from the type-bound arrays into the orbital_sets structure.
         call xermsg ('molecular_orbital_basis_obj', 'get_all_orbital_sets', 'Not implemented yet.',1,1)
 
   end subroutine get_all_orbital_sets

   subroutine get_orbital_coefficient_matrix(this,cf)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      real(kind=cfp), allocatable :: cf(:,:)

      integer :: err, cnt, i

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_coefficient_matrix', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (allocated(cf)) deallocate(cf)
         allocate(cf(this%ao_basis%number_of_functions,this%number_of_functions),stat=err)
         if (err /= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_coefficient_matrix', 'Memory allocation failed.', err, 1)
         end if

         !Place the orbital coefficients into columns of cf symmetry by symmetry.
         cnt = 0
         do i=1,this%no_irr
            if (this%orbital_data(i)%number_of_functions == 0) cycle
            cf(1:this%ao_basis%number_of_functions,cnt+1:cnt+this%orbital_data(i)%number_of_functions) = &
            &this%orbital_data(i)%coefficients(1:this%ao_basis%number_of_functions,1:this%orbital_data(i)%number_of_functions)
            cnt = cnt + this%orbital_data(i)%number_of_functions
         enddo !i

   end subroutine get_orbital_coefficient_matrix

   function is_initialized(this)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      logical :: is_initialized

         is_initialized = this%initialized

   end function is_initialized

   !todo change this work just like the atomic one!!!
   subroutine get_continuum_flags(this,irr,list)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: irr
      logical, allocatable :: list(:)

      integer :: i, err

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_continuum_flag', 'The basis set has not been initialized.', 1, 1)
         end if

         if (irr <= 0 .or. irr > this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'get_continuum_flag', 'On input the value of irr was out of range.', 2, 1)
         end if

         if (allocated(list)) deallocate(list)
         allocate(list(this%orbital_data(irr)%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj', 'get_continuum_flag', 'Memory allocation failed.',err,1)

         list = .false.
         do i=1,this%orbital_data(irr)%number_of_functions
            list(i) = this%is_continuum(this%relative_to_absolute(i,irr))
         enddo !i

   end subroutine get_continuum_flags

   function get_number_of_orbitals(this,irr)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: irr
      integer :: get_number_of_orbitals

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_number_of_orbitals', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (irr <= 0 .or. irr > this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'get_number_of_orbitals', &
                         'On input the value of irr was out of range.', 2, 1)
         end if

         get_number_of_orbitals = this%orbital_data(irr)%number_of_functions

   end function get_number_of_orbitals

   function get_index_within_symmetry(this,absolute_index)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: absolute_index
      integer :: get_index_within_symmetry

         if (.not. this % initialized) then
            call xermsg('molecular_orbital_basis_obj', 'get_index_within_symmetry', &
                        'The basis set has not been initialized.', 1, 1)
         end if

         if (absolute_index <= 0 .or. absolute_index > this % number_of_functions) then
            call xermsg ('molecular_orbital_basis_obj', 'get_index_within_symmetry', &
                         'On input absolute_index was out of range.', 2, 1)
         end if

         get_index_within_symmetry = this%absolute_to_relative(1,absolute_index)

   end function get_index_within_symmetry

   function get_orbital_symmetry(this,absolute_index)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: absolute_index
      integer :: get_orbital_symmetry

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_symmetry', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (absolute_index <= 0 .or. absolute_index > this % number_of_functions) then
            call xermsg ('molecular_orbital_basis_obj', 'get_orbital_symmetry', &
                         'On input absolute_index was out of range.', 2, 1)
         end if

         get_orbital_symmetry = this%absolute_to_relative(2,absolute_index)

   end function get_orbital_symmetry

   function get_absolute_index(this,num,sym)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: num,sym
      integer :: get_absolute_index

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'get_absolute_index', &
                         'The basis set has not been initialized.', 1, 1)
         end if

         if (sym > this % no_irr .or. sym <= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'get_absolute_index', &
                         'On input orbital symmetry was out of range.', 2, 1)
         end if
         if (num > this % orbital_data(sym) % number_of_functions .or. num <= 0) then
            call xermsg ('molecular_orbital_basis_obj', 'get_absolute_index', &
                         'On input orbital index within the given symmetry was out of range.', 3, 1)
         end if

         get_absolute_index = this%relative_to_absolute(num,sym) !the index of the MO within the whole orbital set
  
   end function get_absolute_index 

   subroutine calculate_amplitudes(this,a,normalize_to_a,amplitudes,continuum_channels)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      real(kind=cfp), intent(in) :: a
      logical, intent(in) :: normalize_to_a
      integer, allocatable :: continuum_channels(:,:)
      real(kind=cfp), allocatable :: amplitudes(:,:)

      real(kind=cfp), allocatable :: ao_amplitudes(:,:)
      integer :: n_channels, n_orbitals, err, cnt, i, j, k

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'amplitudes', 'The basis set has not been initialized.', 1, 1)
         end if

         write(level3,'("--------->","molecular_orbital_basis_obj:calculate_amplitudes")')

         if (a .le. 0.0_cfp) call xermsg('molecular_orbital_basis_obj', 'amplitudes', 'On input a .le. 0.0_cfp: invalid input.',2,1)

         !Calculate amplitudes of the atomic functions:
         call this%ao_basis%calculate_amplitudes(a,normalize_to_a,ao_amplitudes,continuum_channels)
         n_channels = size(ao_amplitudes,1)

         !Contract with the molecular orbital coefficients:
         if (allocated(amplitudes)) deallocate(amplitudes)
         allocate(amplitudes(n_channels,this%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj', 'amplitudes', 'Memory allocation failed.',err, 1)

         cnt = 0
         do i=1,this%no_irr

            if (this % orbital_data(i) % number_of_coefficients /= this % ao_basis % number_of_functions) then
                call xermsg ('molecular_orbital_basis_obj', 'amplitudes', 'The AO and MO basis sets are incompatible.', 3, 1)
            end if

            n_orbitals = this%orbital_data(i)%number_of_functions
            if (n_orbitals == 0) cycle

            amplitudes(1:n_channels,cnt+1:cnt+n_orbitals) = matmul(ao_amplitudes,this%orbital_data(i)%coefficients)

            !do j=1,n_orbitals
            !   do k=1,n_channels
            !      if (amplitudes(k,cnt+j) .ne. 0.0_cfp) write(stdout,'(3(i0,1x),e25.15)') i,j,k, amplitudes(k,cnt+j)
            !   enddo !k
            !enddo !j

            cnt = cnt + n_orbitals

         enddo !i

         write(level3,'("<---------","molecular_orbital_basis_obj:calculate_amplitudes")')

   end subroutine calculate_amplitudes

   subroutine orthogonalize(this,overlap_matrix,symmetry,gramm_schmidt,gramm_schmidt_one_by_one,symmetric,sym_ortho_data,&
                            active_start,active_end,passive_start,passive_end,check_overlaps)
      use orthogonalization_gbl, only: GS_ortho_routine, SYM_ortho_routine
      use const_gbl, only: thrs_gs_ortho, thrs_symm_ortho, overlap_ints, int_del_thr, thrs_lin_dep_gs_ortho, thrs_cf_sym_ortho_trans
      use parallel_arrays_gbl, only: p2d_array_obj
      implicit none
      class(molecular_orbital_basis_obj) :: this
      real(kind=cfp), allocatable :: overlap_matrix(:,:)
      type(sym_ortho_io), intent(inout), optional :: sym_ortho_data
      integer, intent(in) :: symmetry
      logical, intent(in), optional :: gramm_schmidt, gramm_schmidt_one_by_one, symmetric, check_overlaps
      integer, intent(in), optional :: active_start, active_end, passive_start, passive_end

      integer :: i, j, k, no_cf, err, a_start, a_end, p_start, p_end, p_start_abs, p_end_abs
      real(kind=cfp), allocatable :: olap_orbs(:,:), thrs_ao(:)
      logical :: do_gs_ortho, do_sym_ortho, all_active, check
      real(kind=cfp) :: thrs, bto_thrs, cgto_thrs

         write(level3,'("--------->","molecular_orbital_basis_obj:orthogonalize")')

         write(level1,'(/,"Orbital orthogonalization")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         do_gs_ortho = .false.
         do_sym_ortho = .false.

         if (present(gramm_schmidt) .and. present(symmetric)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Specify only one of: gramm_schmidt, symmetric but not both.', 2, 1)
         end if
         if (present(gramm_schmidt_one_by_one) .and. (present(gramm_schmidt) .or. present(symmetric))) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'If only gramm_schmidt_one_by_one is given then none of gramm_schmidt and symmetric can be given.', 2, 1)
         end if

         if (present(gramm_schmidt)) do_gs_ortho = gramm_schmidt
         if (present(gramm_schmidt_one_by_one)) do_gs_ortho = gramm_schmidt_one_by_one
         if (present(symmetric)) do_sym_ortho = symmetric

         if (do_gs_ortho .and. present(sym_ortho_data)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'On input sym_ortho_data given but Gramm-Schmidt ortho. requested.', 3, 1)
         end if
         if (do_sym_ortho .and. .not. present(sym_ortho_data)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Symmetric ortho. requested but input sym_ortho_data not given.', 4, 1)
         end if

         if (present(active_start) .neqv. present(active_end)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Only one of active_start,active_end has been specified &
                         &but both are required if one of them is given.', 5, 1)
         endif

         if (present(passive_start) .neqv. present(passive_end)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'Only one of passive_start,passive_end has been specified &
                         &but both are required if one of them is given.', 6, 1)
         endif

         if (symmetry <= 0 .or. symmetry > this % no_irr) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'On input, the value of symmetry was out of range.', 7, 1)
         end if
!
!--------Determine and check absolute indices of the active and passive orbitals. Transfer the values of active_* and passive_* to the variables a_* and p_*.
!
         if (present(active_start)) then

            if (active_end < active_start .or. active_start <= 0 .or. &
                active_start > this % orbital_data(symmetry) % number_of_functions) then
               print *,active_start,active_end,this%orbital_data(symmetry)%number_of_functions
               call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                            'On input, active_start, active_end was incorrect', 8, 1)
            endif

            write(level1,'("Symmetry: ",i1)') symmetry
            write(level1,'("Range of active orbitals within the symmetry: ",2i5)') active_start,active_end

            a_start = active_start
            a_end = active_end

         else
            a_start = 1
            a_end = this%orbital_data(symmetry)%number_of_functions
         endif

         if (present(passive_start)) then

            if (passive_end < passive_start .or. passive_start <= 0 .or. &
                passive_start > this % orbital_data(symmetry) % number_of_functions) then
               print *,passive_start,passive_end
               call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                            'On input, passive_start, passive_end was incorrect', 9, 1)
            endif

            write(level1,'("Symmetry: ",i1)') symmetry
            write(level1,'("Range of passive orbitals within the symmetry: ",2i5)') passive_start,passive_end

            p_start = passive_start
            p_end = passive_end

         else
            p_start = 0
            p_end = 0
         endif

         if (a_start .eq. 1 .and. a_end .eq. this%orbital_data(symmetry)%number_of_functions) then
            all_active = .true.
         else
            all_active = .false.
         endif

         if ((p_end .ge. a_start .and. p_start .le. a_start) .or. (p_end .ge. a_end .and. p_start .le. a_end)) then
            call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                         'The ranges of active and passive orbitals must be disjunct.', 10, 1)
         endif

         !by default we request checking of the overlaps but this can be overrided by the user
         check = .true.
         if (present(check_overlaps)) check = check_overlaps

         no_cf = this%ao_basis%number_of_functions       !the number of AOs in this symmetry
!
!------- Determine the deletion thresholds for the orbital coefficients.
!
         allocate(thrs_ao(no_cf),stat=err)
         if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 3 failed.', err, 1)

         if (do_gs_ortho) then
            bto_thrs = F1MACH(4,cfp_dummy)
            cgto_thrs = thrs_lin_dep_gs_ortho
         elseif (do_sym_ortho) then
            bto_thrs = F1MACH(4,cfp_dummy) !thrs_cf_sym_ortho_trans
            cgto_thrs = thrs_cf_sym_ortho_trans
         endif

         !If the basis contains BTOs then typically there will be many small coefficients 
         !which must not be neglected otherwise the resulting orbitals will not be orthogonal with a sufficient precision.
         if (this%ao_basis%contains_btos()) then
            cgto_thrs = bto_thrs
            write(level2,'(10x,"BTO coefficients with magnitude smaller than: ",e25.15," will be deleted.")') bto_thrs
            write(level2,'(10x,"CGTO coefficients with magnitude smaller than: ",e25.15," will be deleted.")') cgto_thrs
            j = 0
            thrs_ao = cgto_thrs
            do i=1,this%ao_basis%number_of_shells
               k = this%ao_basis%shell_descriptor(5,i) !number of functions in the shell

               !this%ao_basis%shell_descriptor(1,i) == 2 for BTO shells
               if (this%ao_basis%shell_descriptor(1,i) .eq. 2) thrs_ao(j+1:j+k) = bto_thrs

               j = j + k
            enddo !i
         else
            thrs_ao = cgto_thrs
            write(level2,'(10x,"Coefficients with magnitude smaller than: ",e25.15," will be deleted.")') cgto_thrs
         endif
!
!--------Gramm-Schmidt orthogonalization
!
         if (do_gs_ortho) then

            write(level2,'(/,10x,"Gramm-Schmidt orthogonalization requested")')

            !Column indices in the mo2so_range are absolute
            if (p_end .ge. p_start .and. p_end > 0) then
               p_start_abs = this%relative_to_absolute(p_start,symmetry)
               p_end_abs = this%relative_to_absolute(p_end,symmetry)
            else
               p_start_abs = 1
               p_end_abs = 0
               p_start = 1 !setting p_e smaller than p_s guarantees that no orbitals will be treated as fixed in GS_ortho_routine
               p_end = 0
            endif

            if (present(gramm_schmidt_one_by_one)) then
               write(level2,'(10x,"The active orbitals will be orthogonalized one-by-one way.")')
               !Orthogonalize a set of orbitals in the one-by-one fashion (this is most likely to be used for the continuum orthogonalization)
               !$OMP PARALLEL DEFAULT(NONE) &
               !$OMP & PRIVATE(i) &
               !$OMP & SHARED(a_start,a_end,p_start_abs,p_end_abs,no_cf,p_start,p_end,this,overlap_matrix,symmetry,thrs_ao)
               !$OMP DO
               do i=a_start,a_end
                  call GS_ortho_routine (no_cf, p_start, p_end, i, i, this % orbital_data(symmetry) % coefficients, &
                                         overlap_matrix, symmetry, this % mo2so_range(1:2, p_start_abs:max(p_end_abs, 1)), thrs_ao)
               enddo !i
               !$OMP END DO
               !$OMP END PARALLEL
            else
               call GS_ortho_routine (no_cf, p_start, p_end, a_start, a_end, this % orbital_data(symmetry) % coefficients, &
                                      overlap_matrix, symmetry, this % mo2so_range(1:2, p_start_abs:max(p_end_abs,1)), thrs_ao)
            endif

            write(level2,'("...finished.")')

         endif
!
!--------Symmetric orthogonalization
!
         if (do_sym_ortho) then

            err = sym_ortho_data%check(this%no_irr)
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                             'sym_ortho_data%check failed. See the routine for details.', err, 1)
            end if

            if (allocated(sym_ortho_data%to_delete)) deallocate(sym_ortho_data%to_delete)
            allocate(sym_ortho_data%to_delete(1:this%orbital_data(symmetry)%number_of_functions),stat=err)
            if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 1 failed.', err, 1)

            !This array will contain .true. values for those absolute orbital indices which should be deleted from the orbital basis.
            sym_ortho_data%to_delete(:) = .false.

            !Symmetrically orthogonalize a subset of orbitals in the set. Typically, this would be used following G-S orthogonalization of the continuum against the target orbitals 
            !to orthogonalize the subset of continuum orbitals among themselves.
            write(level1,'(/,10x,"Symmetric orthogonalization requested")')

            write(level1,'("Deletion threshold: ",e25.15)') sym_ortho_data%del_thrs(symmetry)

            ! master process will perform the orthogonalization
            if (myrank == master) then
                allocate(olap_orbs(a_start:a_end,a_start:a_end),stat=err)
                if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 2 failed.', err, 1)

                !Form the overlap matrix in the basis of the orbitals to orthogonalize
                !olaps = c**T * S * c
                olap_orbs(a_start:a_end,a_start:a_end) = matmul( &
                    matmul( &
                        transpose(this % orbital_data(symmetry) % coefficients(1:no_cf,a_start:a_end)), &
                        overlap_matrix(1:no_cf,1:no_cf) &
                    ), &
                    this % orbital_data(symmetry) % coefficients(1:no_cf, a_start:a_end) &
                )

                call SYM_ortho_routine (no_cf, a_start, a_end, this % orbital_data(symmetry) % coefficients, &
                                        olap_orbs, sym_ortho_data % del_thrs(symmetry), thrs_ao, sym_ortho_data % to_delete)

                deallocate(olap_orbs)
            end if

            ! master process will share the results, so that all ranks have exactly the same coefficients
            call mpi_mod_bcast(this % orbital_data(symmetry) % coefficients, master)
            call mpi_mod_bcast(sym_ortho_data % to_delete, master)

            k = count(sym_ortho_data%to_delete)
            if (k > 0) then
               write(level1,'(/,"Number of orbitals marked for deletion: ",i0)') k

               write(level2,'("Orbital indices: number.symmetry, absolute index: ")')
               do i=a_start,a_end
                  if (sym_ortho_data % to_delete(i)) then
                     write(level2,'(2i5,i0)') this % absolute_to_relative(1,i), this % absolute_to_relative(2,i), i
                  end if
               enddo
            else
               write(level1,'(/,"There are no orbitals to delete in this symmetry.")')
            endif
         endif
!
!--------The orbitals have changed so recalculate the so2mo,mo2so indices.
!
         call this%determine_auxiliary_indices
!
!--------Optional checking of orbital orthogonality
!
         if (check) then

            if (do_gs_ortho) thrs = thrs_gs_ortho
            if (do_sym_ortho) thrs = thrs_symm_ortho

            if (do_sym_ortho) then !if we performed the symmetric orthogonalization then we must not check orthogonality between the orbitals that are to be deleted
               do i=a_start,a_end
                  if (sym_ortho_data%to_delete(i)) a_end = a_end - 1 !the orbitals to delete are always the last ones in each symmetry so we just lower the index of the last active orbital
               enddo
            endif

            no_cf = this%ao_basis%number_of_functions !the number of AOs in this symmetry

            if (a_end .ge. a_start .and. a_end > 0) then

               write(level1,'(/,10x,"Checking orthogonality of the active orbitals in symmetry: ",i0,/)') symmetry
   
               allocate(olap_orbs(a_start:a_end,a_start:a_end),stat=err)
               if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 3 failed.', err, 1)
   
               !form the overlap matrix in the basis of the orthogonalized orbitals
               !olaps = c**T * S * c
               associate(orbital => this%orbital_data(symmetry))
                  olap_orbs(a_start:a_end,a_start:a_end) = matmul(matmul(transpose(orbital%coefficients(1:no_cf,a_start:a_end)),&
                                                                              &overlap_matrix(1:no_cf,1:no_cf)),&
                                                                              &orbital%coefficients(1:no_cf,a_start:a_end))
               end associate
      
               do i=a_start,a_end
                  write(level3,'(i5,1x,e25.15)') i,olap_orbs(i,i)
                  if (abs(olap_orbs(i,i)-1.0_cfp) .ge. thrs) then
                     write(stdout,'(2e25.15)') thrs,olap_orbs(i,i)
                     call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                                  'Norm of an orbital is not 1 within accuracy given by thrs.', 14, 0)
                  endif
                  !If we orthogonalize in the one-by-one way then we cannot expect that the active orbitals will be orthogonal with each other so this check is skipped.
                  if (.not.(present(gramm_schmidt_one_by_one))) then
                     do j=i+1,a_end
                        if (olap_orbs(j,i) .ge. thrs) then
                           write(stdout,'(i0,1x,i0,2e25.15)') j,i,olap_orbs(j,i),thrs
                           call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                                        'Overlap between two different orbitals is .ge. thrs.', 15, 0)
                        endif
                     enddo
                  endif
               enddo
   
               deallocate(olap_orbs)

               if (p_end .ge. p_start .and. p_end > 0) then
                  write(level1,'(/,10x,"Checking orthogonality of the active orbitals &
                                        &with respect to the fixed orbitals in symmetry: ",i0,/)') symmetry
   
                  allocate(olap_orbs(p_start:p_end,a_start:a_end),stat=err)
                  if (err /= 0) then
                     call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', 'Memory allocation 4 failed.', err, 1)
                  end if
      
                  !form the overlap matrix in the basis of the orthogonalized orbitals
                  !olaps = c**T * S * c
                  associate(orbital => this%orbital_data(symmetry))
                     olap_orbs(p_start:p_end,a_start:a_end) = matmul(matmul(transpose(orbital%coefficients(1:no_cf,p_start:p_end)),&
                                                                                   &overlap_matrix(1:no_cf,1:no_cf)),&
                                                                                   &orbital%coefficients(1:no_cf,a_start:a_end))
                  end associate

                  do i=a_start,a_end
                     do j=p_start,p_end
                        if (olap_orbs(j,i) .ge. thrs) then
                           write(stdout,'(i0,1x,i0,2e25.15)') j,i,olap_orbs(j,i),thrs
                           call xermsg ('molecular_orbital_basis_obj', 'orthogonalize', &
                                        'Overlap between an active and a passive orbital is .ge. thrs.', 17, 0)
                        endif
                     enddo
                  enddo
      
                  deallocate(olap_orbs)
               endif

            endif

            write(level1,'(/,"Orthogonalization complete.")')

         else

            write(level1,'(/,"Checking of orbital overlaps not requested. Orthogonalization complete.")')

         endif

         write(level3,'("<---------","done:molecular_orbital_basis_obj:orthogonalize")')
 
   end subroutine orthogonalize

   subroutine delete_small_coefficients(this)
      use const_gbl, only: thrs_orb_cf
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: symmetry, j, k

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'delete_small_coefficients', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         write(level3,'("--------->","molecular_orbital_basis_obj:delete_small_coefficients")')

         write(level2,'("Removing orbital coefficients with magnitude smaller than: ",e25.15)') thrs_orb_cf

         do symmetry=1,size(this%orbital_data)
            do k=1,size(this%orbital_data(symmetry)%coefficients,2)
               do j=1,this%orbital_data(symmetry)%number_of_coefficients
                  if (abs(this % orbital_data(symmetry) % coefficients(j,k)) < thrs_orb_cf) then
                     this % orbital_data(symmetry) % coefficients(j,k) = 0.0_cfp
                  end if
               enddo !j
            enddo !k
         enddo !i

         write(level3,'("<---------","done:molecular_orbital_basis_obj:delete_small_coefficients")')

   end subroutine delete_small_coefficients

   subroutine delete_orbitals(this,symmetry,to_delete)
      implicit none
      class(molecular_orbital_basis_obj) :: this
      integer, intent(in) :: symmetry
      logical, intent(in) :: to_delete(:)

      type(orbital_data_obj) :: temp_orbitals
      integer :: n_orbitals_to_keep, err, i, cnt

         write(level3,'("--------->","molecular_orbital_basis_obj:delete_orbitals")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'delete_orbitals', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         if (size(to_delete) /= this % orbital_data(symmetry) % number_of_functions) then
            call xermsg ('molecular_orbital_basis_obj', 'delete_orbitals', &
                         'The size of the array on input is not equal to the number of orbitals in the given symmetry.', 2, 1)
         end if

         temp_orbitals = this%orbital_data(symmetry)

         n_orbitals_to_keep = this%orbital_data(symmetry)%number_of_functions - count(to_delete)

         deallocate(this % orbital_data(symmetry) % energy, &
                    this % orbital_data(symmetry) % occup, &
                    this % orbital_data(symmetry) % spin, &
                    this % orbital_data(symmetry) % coefficients)

         allocate(this % orbital_data(symmetry) % energy(n_orbitals_to_keep), &
                  this % orbital_data(symmetry) % occup(n_orbitals_to_keep), &
                  this % orbital_data(symmetry) % spin(n_orbitals_to_keep), &
                  this % orbital_data(symmetry) % coefficients(this % orbital_data(symmetry) % number_of_coefficients, &
                    n_orbitals_to_keep), stat = err)
         if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj', 'delete_orbitals', 'Memory allocation failed.', err, 1)
 
         this%orbital_data(symmetry)%number_of_functions = n_orbitals_to_keep

         !Keep only the orbitals that are not marked for deletion
         cnt = 0
         do i=1,temp_orbitals%number_of_functions
            if (.not.(to_delete(i))) then
               cnt = cnt + 1
               this%orbital_data(symmetry)%energy(cnt) = temp_orbitals%energy(i)
               this%orbital_data(symmetry)%occup(cnt) = temp_orbitals%occup(i)
               this%orbital_data(symmetry)%spin(cnt) = temp_orbitals%spin(i)
               this%orbital_data(symmetry)%coefficients(1:temp_orbitals%number_of_coefficients,cnt) &
                    = temp_orbitals%coefficients(1:temp_orbitals%number_of_coefficients,i)
            else
               write(level1,'("Orbital ",i0," from symmetry ",i1," has been deleted.")') i, symmetry
            endif
         enddo !i

         this%number_of_functions = sum(this%orbital_data(:)%number_of_functions)

         call this%determine_auxiliary_indices

         write(level3,'("<---------","molecular_orbital_basis_obj:delete_orbitals")')

   end subroutine delete_orbitals

   function check_sym_ortho_io(this,last_irr)
      implicit none
      class(sym_ortho_io) :: this
      integer :: check_sym_ortho_io
      integer, intent(in) :: last_irr

      integer :: i

         check_sym_ortho_io = 0

         do i=1,last_irr
            if (this%del_thrs(i) .le. 0.0_cfp) check_sym_ortho_io = i
         enddo
      
   end function check_sym_ortho_io

   subroutine determine_auxiliary_indices(this)
      use const_gbl, only: abel_prod_tab
      implicit none
      class(molecular_orbital_basis_obj) :: this

      integer :: no_tot, no_ao, i, j, k, l, col, err, ij_irr, p, q, kl_irr, kl, ij, no_sym_pairs, n_p, no_pairs
      integer, allocatable :: sym_seq_no(:)
      logical, allocatable :: ao_is_continuum(:)
      logical :: all_zero

         no_tot = this%number_of_functions
         no_ao = this%ao_basis%number_of_functions

         no_pairs = this%number_of_functions*(this%number_of_functions+1)/2
         no_sym_pairs = this%no_irr*(this%no_irr+1)/2   !number of unique pairs of symmetries
         no_sym_pairs = no_sym_pairs*(no_sym_pairs+1)/2 !number of unique pairs of pairs of symmetries 

         n_p = maxval(this%orbital_data(:)%number_of_functions)
         n_p = n_p*(n_p+1)/2

         if (allocated(this%block_offset)) deallocate(this%block_offset)
         if (allocated(this%sym_offset)) deallocate(this%sym_offset)
         if (allocated(this%is_continuum)) deallocate(this%is_continuum)
         if (allocated(this%absolute_to_relative)) deallocate(this%absolute_to_relative)
         if (allocated(this%relative_to_absolute)) deallocate(this%relative_to_absolute)
         if (allocated(this%so2mo_range)) deallocate(this%so2mo_range)
         if (allocated(this%mo2so_range)) deallocate(this%mo2so_range)

         allocate(this % block_offset(no_sym_pairs + 1), &
                  this % sym_offset(n_p), &
                  this % is_continuum(no_tot), &
                  this % absolute_to_relative(2,no_tot), &
                  this % relative_to_absolute(no_ao, this % no_irr), &
                  stat = err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj','add_shell','Memory allocation 2 failed.',err,1)

         this%is_continuum = .false.
         this%block_offset = 0
         this%absolute_to_relative = 0

         allocate(this%so2mo_range(1:2,no_ao),this%mo2so_range(1:2,no_tot),stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj','add_shell','Memory allocation failed.',err,1)

         !determine the so2mo_range, mo2so_range scanning the whole orbital coefficient data and looking for orbitals which have non-zero AO coefficients.
         this%so2mo_range(1,:) = no_tot+1
         this%so2mo_range(2,:) = 0
         this%mo2so_range(1,:) = no_ao+1
         this%mo2so_range(2,:) = 0
         col = 0
         !loop over all symmetries:
         do i=1,this%number_of_shells
            associate(orbital => this%orbital_data(i))
            !Get continuum flags for each AO for this symmetry.
            call this%ao_basis%get_continuum_flags(this%orbital_data(i)%irr,ao_is_continuum)
            !all orbitals within a symmetry:
            do j=1,orbital%number_of_functions
               col = sum(this%orbital_data(1:i-1)%number_of_functions) + j !the MOs are stored sequentially in columns symmetry by symmetry.
               !write(stdout,'("Orbital with num.sym ",i,".",i," has index",i)') orbital%sym_ind, orbital%symmetry, col
               this%absolute_to_relative(1:2,col) = (/j,orbital%irr/) !index within its own symmetry (j) and symmetry of the orbital (orbital%irr)
               this%relative_to_absolute(j,orbital%irr) = col !the index of the MO within the whole orbital set
               all_zero = .true.
               do k=1,no_ao
                  if (orbital%coefficients(k,j) .ne. 0.0_cfp) then
                     all_zero = .false.
                     if (col < this%so2mo_range(1,k)) this%so2mo_range(1,k) = col
                     if (col > this%so2mo_range(2,k)) this%so2mo_range(2,k) = col
                     if (k < this%mo2so_range(1,col)) this%mo2so_range(1,col) = k
                     if (k > this%mo2so_range(2,col)) this%mo2so_range(2,col) = k
                     if (ao_is_continuum(k)) this%is_continuum(col) = .true. !this orbital must represent the continuum since it contains non-zero coefficients for at least one continuum AO.
                  endif
               enddo !k
               ! If the basis includes trivial orbitals then the identification of target and continuum orbitals becomes ambiguous.
               ! This is not necessarily a problem but there is no physical reason why such orbitals should be included in the basis.
               if (all_zero) then
                  print *,i,j
                  print *,orbital%coefficients(:,j)
                  call xermsg ('molecular_orbital_basis_obj', 'add_shell', &
                               'The orbital basis includes at least one orbital with only zero-containing coefficient array. &
                               &This is not allowed.', 7, 1)
               endif
            enddo !j
            end associate
         enddo !i

         !Determine the number of transformed MO 2-particle symmetric integrals for each unique quadruplet of symmetries (ij|kl)
         allocate(sym_seq_no,source=this%block_offset,stat=err)
         if (err .ne. 0) call xermsg('molecular_orbital_basis_obj','add_shell','Memory allocation 2 failed.',err,1)

         this%block_offset = 0
         sym_seq_no = 0
         p = 0
         do i=1,this%no_irr
            do j=1,i
               !How many unique pairs of MOs there are for the symmetry pair (i,j)
               if (i .eq. j) then
                  ij = this%orbital_data(i)%number_of_functions*(this%orbital_data(i)%number_of_functions+1)/2
               else
                  ij = this%orbital_data(i)%number_of_functions*this%orbital_data(j)%number_of_functions
               endif 
               ij_irr = abel_prod_tab(i,j) !IRR of the product of symmetries (i,j)
               p = i*(i-1)/2+j             !index of the unique pair of symmetries (i,j)
               q = 0
               do k=1,i
                  do l=1,k
                     q = k*(k-1)/2+l
                     kl_irr = abel_prod_tab(k,l)
                     if (p .ge. q .and. ij_irr .eq. kl_irr) then !only unique quartets of symmetries; additonally we have the symmetry restriction IRR_ij .eq. IRR_kl
                        if (k .eq. l) then
                           kl = this%orbital_data(k)%number_of_functions*(this%orbital_data(k)%number_of_functions+1)/2
                        else
                           kl = this%orbital_data(k)%number_of_functions*this%orbital_data(l)%number_of_functions
                        endif
                        !below we compute the number of integrals for the quartet of symmetries (p,q), p .ge. q.
                        if (p .eq. q) then  !(i,j) pair = (k,l) pair
                           sym_seq_no(p*(p-1)/2+q) = ij*(ij+1)/2 !only unique quartets of MOs
                        else
                           sym_seq_no(p*(p-1)/2+q) = ij*kl
                        endif
                        !write(stdout,'("block",4i,2i)') max(i,j),min(i,j),max(k,l),min(k,l), p*(p-1)/2+q, sym_seq_no(p*(p-1)/2+q)
                     endif
                  enddo !l
               enddo !k
            enddo !j
         enddo !i

         !Clearly, from the construction below, block_offset(p) contains the total number of integrals that preceed the given symmetry block 'p=s*(s-1)/2+t'. The values s,t are defined as:
         !s=i*(i-1)/2+j, t=k*(k-1)/2+l where i,j,k,l are the IRRs of the MOs in (ab|cd) and we assume i .ge. j, k .ge. l.
         !The 'integral block' is defined here as the set of integrals (ab|cd) where the MOs a,b,c,d come from the symmetry block 'p'.
         !Within each integral block the index of the integral is given depending on the symmetries of the MOs: (II|II), (IJ|IJ), (II|JJ), (IJ|KL).
         !The full index therefore is: block_offset(p) + block_index; the expressions for the block_index for the four cases above can be found in the relevent 2-particle indexing routines.

         j = size(this%block_offset)-1 !total number of unique quartets of symmetries; note that the last value in block_offset is the total number of 2-particle symmetric integrals.
         this%block_offset(j+1) = sum(sym_seq_no)

         do i=2,j
            this%block_offset(i) = sum(sym_seq_no(1:i-1)) !compute the number of integrals preceeding the block i.
         enddo !i
         this%block_offset(1) = 0

         deallocate(sym_seq_no)

    end subroutine determine_auxiliary_indices

!---- ROUTINES FOR 2P TRANSFORMATION:
    !> This routine transposes the 2D matrix on input while determining whether all elements of the array are zero.
    subroutine transpose_2d(batch_in,batch_t,nrow,ncol)
       use const_gbl, only: tile
       implicit none
       integer, intent(in) :: nrow,ncol
       real(kind=cfp), allocatable :: batch_in(:,:) !(nrow,ncol)
       real(kind=cfp), allocatable :: batch_t(:,:) !(ncol,nrow)
       logical :: all_zero
 
       integer :: i, j, ii, jj, iend, jend
 
          all_zero = .true.
 
          do jj = 1,nrow,tile
             jend = min (nrow,jj+tile-1)
             do ii = 1,ncol,tile
                iend = min (ncol,ii+tile-1)
 
                if (.not.(all_zero)) then
                   do j = jj,jend
                      do i = ii,iend
                         batch_t(i,j) = batch_in(j,i)
                      enddo
                   enddo
                else
                   do j = jj,jend
                      do i = ii,iend
                         batch_t(i,j) = batch_in(j,i)
                         if (batch_t(i,j) .ne. 0.0_cfp) all_zero = .false.
                      enddo
                   enddo
                endif
 
             enddo
          enddo
       
    end subroutine transpose_2d
 
    !> \warning Note that the use of the 'allocatable' attribute for the argument arrays is key for performance since this
    !> attribute allows the compiler to assume that the arrays are contiguous in memory.
    !> Alternatively the 'contiguous' attribute can be used but it is a F2008 feature so we omit it here.
    !> It is assumed that the threads have been launched outside of this routine.
    subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local(ao_integrals,int_type,rs_start,rs_end,ij_type,iqrs,&
                                        hlp,mobas,cf,cf_t,ijrs,no_ao,last_tgt_ao,no_int,no_pairs,rs_ind,two_p_continuum,&
                                        index_scheme)
       use omp_lib
       implicit none
       class(molecular_orbital_basis_obj), intent(in) :: mobas
       type(p2d_array_obj), intent(inout) :: ao_integrals
       integer, intent(in) :: int_type,rs_start,rs_end,no_ao,last_tgt_ao,no_pairs,rs_ind(2,no_pairs),index_scheme
       integer(kind=1), allocatable :: ij_type(:)
       integer, intent(inout) :: no_int
       logical, intent(in) :: two_p_continuum
 
       !we use the allocatable attribute so that the compiler can assume that the arrays are unit-stride:
       !this has ~20% impact on the performance of this routine.
       real(kind=cfp), allocatable :: iqrs(:,:), hlp(:), ijrs(:,:), cf(:,:), cf_t(:,:)
 
       integer :: rs, pq, p, q, ij, i, j, max_p, max_q, ijrs_type, pq_mapped, rs_mapped, rs_tmp, iam, r, s
       real(kind=cfp) :: mo_int
       logical :: ao_contains_continuum, dont_two_p_continuum

          if (mobas%ao_basis%n_cont_fns > 0) then
             ao_contains_continuum = .true.
          else
             ao_contains_continuum = .false.
          endif

          dont_two_p_continuum = ao_contains_continuum .and. .not.(two_p_continuum)

          !Split the rs-indices in the block among threads
          !$OMP DO SCHEDULE(DYNAMIC)
          do rs=rs_start,rs_end !rs is the index of the (rs)-pair corresponding to the AO integrals (pq|rs) that we want to process here.
             iam = omp_get_thread_num()
       
             ! Pre-load all (available) AO integrals (pq|O|rs) (== (pq|rs)) with the given (rs) index
             ! into hlp array - this speeds up the 1st step a little bit compared with accessing ao_integrals%a
             ! in the j-loop below. We preload the integrals always only one row at a time so that we don't need
             ! an extra array of size ~no_ao**2.
    
             !1st step - transform the first index: p->i
             ! We loop over all AO integrals (pq|rs) with the given rs index: if integrals for only 1p in the continuum
             ! are required then we skip the AO integrals with 2p in the continuum.

             r = rs_ind(1,rs)
             s = rs_ind(2,rs)
 
             max_p = no_ao
             if (dont_two_p_continuum) then
                !rs pair is CC => we need only (TT|CC) AO integrals
                if (mobas%ao_basis%get_function_pair_type(r,s) .eq. 3) max_p = last_tgt_ao
             endif

             !r = rs_ind(1,rs)
             !s = rs_ind(2,rs)

             !Determine the ordered index of this AO pair
             rs_mapped = mobas%ao_basis%function_pair_ordered_index(r,s)
 
             iqrs = 0.0_cfp
             do p=1,max_p
                max_q = p
                if (dont_two_p_continuum) then
                   !rs pair is CT and p is C => we must load only (TT|CT), (CT|CT) AO integrals
                   if (mobas%ao_basis%get_function_pair_type(r,s) .eq. 2 .and. p > last_tgt_ao) max_q = last_tgt_ao
                endif
                ij = p*(p-1)/2
                !print *,max_p,max_q,ij
                !Preload the integrals (p,q) for the current p into an intermediate buffer hlp.
                do q=1,max_q
                   pq = ij + q

                   !Determine the ordered index of this AO pair
                   pq_mapped = mobas%ao_basis%function_pair_ordered_index(p,q)

                   ! WARNING: note that we have inlined here the indexing function for AOs (i.e. computation of 'i')!
                   ! We are not using the generic index function for performance reasons.
                   pq = max(pq_mapped,rs_mapped)
                   rs_tmp = min(pq_mapped,rs_mapped)

                   i = pq*(pq-1)/2+rs_tmp !standard indexing for TTTT CTCT CTTT classes and for 2p in the continuum

                   if (index_scheme == 2) then
                      i = compact_integral_index(r, s, p, q)
                   else if (dont_two_p_continuum) then !Special indexing for CCTT class in case 1p in the continuum
                      if (pq .le. mobas%ao_basis%n_TT_pairs .or. rs_tmp .le. mobas%ao_basis%n_TT_pairs) then !is pq or rs TT pair?
                         if (pq > mobas%ao_basis%last_CT_fn) then !pq=CC, rs=TT
                            i = mobas%ao_basis%n_prec_ints + rs_tmp + mobas%ao_basis%n_TT_pairs*(pq-mobas%ao_basis%last_CT_fn-1)
                         elseif(rs_tmp > mobas%ao_basis%last_CT_fn) then !pq=TT, rs=CC pair
                            i = mobas%ao_basis%n_prec_ints + pq + mobas%ao_basis%n_TT_pairs*(rs_tmp-mobas%ao_basis%last_CT_fn-1)
                         endif
                      endif
                   endif
                   !if (pq .ge. rs) then
                   !   i = pq*(pq-1)/2+rs
                   !else
                   !   i = rs*(rs-1)/2+pq
                   !endif
                   if (i == 0) then
                      hlp(q) = 0
                   else
                      if (i > size(ao_integrals%a,1)) then
                         print *,p,q,pq,rs_tmp,i
                      endif
                      hlp(q) = ao_integrals%a(i,int_type)
                   end if
                   !write(100+myrank,'(4i4,e25.15)') p,q,r,s,hlp(q)
                enddo !q
 
                !Loop over the preloaded integrals and for each of them calculate its contribution to iqrs.
                do q=1,min(p-1,max_q)
                   if (hlp(q) .ne. 0.0_cfp) then
                      !We assume that the AO hlp(q) is symmetric: (pq|O|rs) = (qp|O|rs) so we calculate at once contributions of (pq|O|rs) to iqrs(i,q) and iqrs(i,p)
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the p-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                      do i = mobas%so2mo_range(1,q), mobas%so2mo_range(2,q) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,p) = iqrs(i,p) + hlp(q)*cf_t(i,q)
                      enddo !i
                   endif
                enddo !q
 
                if (max_q .eq. p) then
                   q = p
                   if (hlp(q) .ne. 0.0_cfp) then
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                   endif
                endif
             enddo !p
       
             no_int = no_int + no_pairs
       
             !2nd step - transform the second index: q->j
             !This step takes a little over half of the compute time for the whole first step (for no symmetry case)
             ij = 0
             p = rs-rs_start+1 !relative rs-index
             do i=1,mobas%number_of_functions
                !Here transpose one row of iqrs: this allows vectorization in the q-loop below and ensures cache locality.
                do q=1,no_ao
                   hlp(q) = iqrs(i,q)
                enddo
                do j=1,i
                   ij = ij + 1
                   !todo implement the same thing as in the not_local equivalent
                   if (dont_two_p_continuum) then
                      !Each pair, (i,j) or (r,s), is assigned a numeric value defining its type: TT=1, TC=2, CC=3.
                      !Within this scheme the pairs of types of [MO,AO] that can occur are: [1,1], [2,1], [1,2], [3,1], [1,3], [2,2].
                      !If we want integrals only for 1p in the continuum then the types that we want to skip have indices: [2,3], [3,2], [3,3].
                      !These pairs can be identified using the sum of their types: 5 and 6. This is what we test below.
                      ijrs_type = ij_type(ij)+mobas%ao_basis%get_function_pair_type(r,s)
                      if (ijrs_type .eq. 5 .or. ijrs_type .eq. 6) cycle
                   endif
                   !todo once AO symmetry has been implemented then I can check here that sym(i*j) = sym(r*s). If not then skip this j.
                   mo_int = 0.0_cfp
                   do q=mobas%mo2so_range(1,j), mobas%mo2so_range(2,j) !over all atomic functions which contribute to the j-th molecular orbital.
                      mo_int = mo_int + hlp(q)*cf(q,j)
                   enddo !q
                   ijrs(ij,p) = mo_int !(ij|O|rs), where i,j are MOs and r,s are AOs.
                enddo !j
             enddo !i
          enddo !rs
          !$OMP END DO
 
    end subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local

    !> \warning Note that the use of the 'allocatable' attribute for the argument arrays is key for performance since this
    !>          attribute allows the compiler to assume that the arrays are contiguous in memory.
    !> Alternatively the 'contiguous' attribute can be used but it is a F2008 feature so we omit it here.
    !> It is assumed that the threads have been launched outside of this routine.
    subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_not_local(ao_integrals,int_type,rs_start,rs_end,ij_type,&
                                iqrs,hlp,mobas,cf,cf_t,ijrs,no_ao,last_tgt_ao,no_int,no_pairs,rs_ind,two_p_continuum)
       use omp_lib
       use gto_routines_gbl, only: find_mapping, index_1p_continuum
       implicit none
       class(molecular_orbital_basis_obj), intent(in) :: mobas
       type(p2d_array_obj), intent(inout) :: ao_integrals
       integer, intent(in) :: int_type,rs_start,rs_end,no_ao,last_tgt_ao,no_pairs
       integer, allocatable :: rs_ind(:,:)
       integer(kind=1), allocatable :: ij_type(:)
       integer, intent(inout) :: no_int
       logical, intent(in) :: two_p_continuum
 
       !we use the allocatable attribute so that the compiler can assume that the arrays are unit-stride: this has ~20% impact on the performance of this routine.
       real(kind=cfp), allocatable :: iqrs(:,:), hlp(:), ijrs(:,:), cf(:,:), cf_t(:,:)
 
       integer :: rs, pq, p, q, ij, i, j, k, l, max_p, max_q, ijrs_type, iam, n_shell_r,n_shell_s,n_shell_p,n_shell_q
       integer :: p_rel,q_rel,r_rel,s_rel,ind_start_r,ind_start_s,ind_start_p,ind_start_q,i1,i2,i3,i4
       integer :: ind_orig(4), n(4), map(4), n_map(3), r_shell,s_shell,p_shell,q_shell,r,s
       real(kind=cfp) :: mo_int
       logical :: ao_contains_continuum, dont_two_p_continuum, rs_is_CC, rs_is_TT, pq_is_CC, pq_is_TT, is_CCTT

          if (mobas%ao_basis%n_cont_fns > 0) then
             ao_contains_continuum = .true.
          else
             ao_contains_continuum = .false.
          endif

          dont_two_p_continuum = ao_contains_continuum .and. .not.(two_p_continuum)

          !Split the rs-indices in the block among threads
          !todo this could be simplified by looping instead over all my quartets
          !of shells in ao_integrals%block_offset(:) and picking out only those
          ![pq|rs] where the pq or rs part falls within rs_start,rs_end
          !$OMP DO SCHEDULE(DYNAMIC)
          do rs=rs_start,rs_end !rs is the index of the (rs)-pair corresponding to the AO integrals (pq|rs) that we want to process here.
             iam = omp_get_thread_num()
       
             ! Pre-load all (available) AO integrals (pq|O|rs) (== (pq|rs)) with the given (rs) index into hlp array - this
             ! speeds up the 1st step a little bit compared with accessing ao_integrals%a in the j-loop below.
             ! We preload the integrals always only one row at a time so that we don't need an extra array of size ~no_ao**2.
    
             !1st step - transform the first index: p->i
             ! We loop over all AO integrals (pq|rs) with the given rs index: if integrals for only 1p in the continuum
             ! are required then we skip the AO integrals with 2p in the continuum.

             r = rs_ind(1,rs)
             s = rs_ind(2,rs)

             rs_is_CC = .false.
             rs_is_TT = .false.
             max_p = no_ao
             if (dont_two_p_continuum) then
                !rs pair is CC => we need only (TT|CC) AO integrals
                if (mobas%ao_basis%get_function_pair_type(r,s) .eq. 3) then
                   max_p = last_tgt_ao
                   rs_is_CC = .true.
                elseif (mobas%ao_basis%get_function_pair_type(r,s) .eq. 1) then
                   rs_is_TT = .true.
                endif
             endif

             !r = rs_ind(1,rs)
             !s = rs_ind(2,rs)

             r_shell = mobas%ao_basis%indices_to_shells(1,r) !index of the shell the r-th function is part of
             s_shell = mobas%ao_basis%indices_to_shells(1,s) !index of the shell the s-th function is part of
             k = max(r_shell,s_shell)
             l = min(r_shell,s_shell)

             r_rel = mobas%ao_basis%indices_to_shells(2,r) !index of the r-th function within the shell to which it belongs
             s_rel = mobas%ao_basis%indices_to_shells(2,s) !index of the s-th function within the shell to which it belongs
             ind_start_r = mobas%ao_basis%shell_descriptor(4,r_shell) !starting index for the shell the r-function is part of
             ind_start_s = mobas%ao_basis%shell_descriptor(4,s_shell) !starting index for the shell the s-function is part of
             n_shell_r = mobas%ao_basis%shell_descriptor(5,r_shell) !number of functions in the shell of which the r-function is part of 
             n_shell_s = mobas%ao_basis%shell_descriptor(5,s_shell) !number of functions in the shell of which the s-function is part of

             iqrs = 0.0_cfp
             do p=1,max_p
                max_q = p

                if (dont_two_p_continuum) then
                   !rs pair is CT and p is C => we must load only (TT|CT), (CT|CT) AO integrals
                   if (mobas%ao_basis%get_function_pair_type(r,s) .eq. 2 .and. p > last_tgt_ao) max_q = last_tgt_ao
                endif

                pq = p*(p-1)/2

                p_shell = mobas%ao_basis%indices_to_shells(1,p) !index of the shell the p-th function is part of
                ind_start_p = mobas%ao_basis%shell_descriptor(4,p_shell) !starting index for the shell the p-function is part of
                n_shell_p = mobas%ao_basis%shell_descriptor(5,p_shell) !number of functions in the shell of which the p-function is part of 
                p_rel = mobas%ao_basis%indices_to_shells(2,p) !index of the p-th function within the shell to which it belongs

                !Preload the integrals [pq|rs] for fixed p,r,s into the intermediate buffer hlp.
                do q=1,max_q
                   pq = pq + 1

                   q_shell = mobas%ao_basis%indices_to_shells(1,q) !index of the shell the q-th function is part of
                   i = max(p_shell,q_shell)
                   j = min(p_shell,q_shell)

                   pq_is_CC = .false.
                   pq_is_TT = .false.
                   is_CCTT = .false.
                   if (dont_two_p_continuum) then
                      if (mobas%ao_basis%get_function_pair_type(p,q) .eq. 1) then
                         pq_is_TT = .true.
                      elseif (mobas%ao_basis%get_function_pair_type(p,q) .eq. 3) then
                         pq_is_CC = .true.
                      endif
                      if (pq_is_TT .and. rs_is_CC) then
                         is_CCTT = .true.
                      elseif (pq_is_CC .and. rs_is_TT) then
                         is_CCTT = .true.
                      endif
                   endif

                   !i = index of the quartet of shells to which the [pq|rs] integral belongs
                   i = index_1p_continuum(mobas % ao_basis % ordered_shell_pairs, &
                                          i, j, k, l, is_CCTT, &
                                          mobas % ao_basis % last_CT_sh, &
                                          mobas % ao_basis % n_prec_sh, &
                                          mobas % ao_basis % n_TT_sh_pairs)

                   if (ao_integrals%block_offset(i) .eq. -1) then !The requested integral is not kept by this task
                      hlp(q) = 0.0_cfp
                      cycle
                   endif

                   !Each function p,q,r,s corresponds to a given shell. In the
                   !section below we permute the p,q,r,s indices to the standard order (determined by the starting indices of the functions in each shell)
                   !which was used to save the integrals (see atomic_orbital_basis_obj%two_electron_integrals).
                   ind_start_q = mobas%ao_basis%shell_descriptor(4,q_shell) !starting index for the shell the q-function is part of
                   n_shell_q = mobas%ao_basis%shell_descriptor(5,q_shell) !number of functions in the shell of which the q-function is part of
                   q_rel = mobas%ao_basis%indices_to_shells(2,q) !index of the q-th function within the shell to which it belongs

                   ind_orig(1:4) = (/ind_start_p,ind_start_q,ind_start_r,ind_start_s/)
                   n(1:4) = (/n_shell_p,n_shell_q,n_shell_r,n_shell_s/)

                   !Map the order of the corresponding p,q,r,s shells to the order in which the integrals within the shell are saved.
                   call find_mapping(ind_orig,n,n_map,map)

                   !Permute the relative indices into the order in which the integrals were saved
                   ind_orig(1:4) = (/p_rel,q_rel,r_rel,s_rel/)
                   i1 = ind_orig(map(1))
                   i2 = ind_orig(map(2))
                   i3 = ind_orig(map(3))
                   i4 = ind_orig(map(4))

                   !Compute the index of the integral [pq|rs] within its own
                   !quartet of shells and add to it the offset for the corresponding quartet of shells.
                   i = ao_integrals%block_offset(i) + i1 + n_map(1)*(i2-1) + n_map(2)*(i3-1) + n_map(3)*(i4-1)

                   if (i > size(ao_integrals%a,1)) then
                      print *,'indexing error:',p,q,r,s,i
                   endif
                   hlp(q) = ao_integrals%a(i,int_type)
                   !if (hlp(q) .ne. 0.0_cfp) write(100+myrank,'(4i4,e25.15)') p,q,r,s,hlp(q)
                enddo !q

                !Loop over the preloaded integrals and for each of them calculate its contribution to iqrs.
                do q=1,min(p-1,max_q)
                   if (hlp(q) .ne. 0.0_cfp) then
                      !We assume that the AO hlp(q) is symmetric: (pq|O|rs) = (qp|O|rs) so we calculate at once contributions of (pq|O|rs) to iqrs(i,q) and iqrs(i,p)
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the p-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                      do i = mobas%so2mo_range(1,q), mobas%so2mo_range(2,q) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,p) = iqrs(i,p) + hlp(q)*cf_t(i,q)
                      enddo !i
                   endif
                enddo !q
 
                if (max_q .eq. p) then
                   q = p
                   if (hlp(q) .ne. 0.0_cfp) then
                      do i = mobas%so2mo_range(1,p), mobas%so2mo_range(2,p) !over all molecular orbitals to which the q-th atomic function contributes.
                         iqrs(i,q) = iqrs(i,q) + hlp(q)*cf_t(i,p)
                      enddo !i
                   endif
                endif
             enddo !p
       
             no_int = no_int + no_pairs

             !2nd step - transform the second index: q->j
             !This step takes a little over half of the compute time for the whole first step (for no symmetry case)
             ij = 0
             p = rs-rs_start+1 !relative rs-index
             do i=1,mobas%number_of_functions
                !Here transpose one row of iqrs: this allows vectorization in the q-loop below and ensures cache locality.
                do q=1,no_ao
                   hlp(q) = iqrs(i,q)
                enddo
                do j=1,i
                   ij = ij + 1
                   !todo implement the same thing as in the not_local equivalent
                   if (dont_two_p_continuum) then
                      !Each pair, (i,j) or (r,s), is assigned a numeric value defining its type: TT=1, TC=2, CC=3.
                      !Within this scheme the pairs of types of [MO,AO] that can occur are: [1,1], [2,1], [1,2], [3,1], [1,3], [2,2].
                      !If we want integrals only for 1p in the continuum then the types that we want to skip have indices: [2,3], [3,2], [3,3].
                      !These pairs can be identified using the sum of their types: 5 and 6. This is what we test below.
                      ijrs_type = ij_type(ij)+mobas%ao_basis%get_function_pair_type(r,s)
                      if (ijrs_type .eq. 5 .or. ijrs_type .eq. 6) cycle
                   endif
                   !todo once AO symmetry has been implemented then I can check here that sym(i*j) = sym(r*s). If not then skip this j.
                   mo_int = 0.0_cfp
                   do q=mobas%mo2so_range(1,j), mobas%mo2so_range(2,j) !over all atomic functions which contribute to the j-th molecular orbital.
                      mo_int = mo_int + hlp(q)*cf(q,j)
                   enddo !q
                   ijrs(ij,p) = mo_int !(ij|O|rs), where i,j are MOs and r,s are AOs.
                enddo !j
             enddo !i

          enddo !rs
          !$OMP END DO
 
    end subroutine omp_two_p_transform_pqrs_block_to_ijrs_AO_is_not_local
 
    subroutine generate_ij_offset(mobas,ij_type,ij_orbital_range,two_p_continuum,ij_offset,n_integrals)
       use const_gbl, only: abel_prod_tab
       use mpi_gbl
       use omp_lib
       use special_functions_gbl, only: ipair
       use sort_gbl, only: sort_int_float
       implicit none
       integer, allocatable :: ij_orbital_range(:,:), ij_offset(:)
       integer, intent(out) :: n_integrals
       integer(kind=1) :: ij_type(:)
       class(molecular_orbital_basis_obj) :: mobas
       logical, intent(in) :: two_p_continuum
 
       integer :: ij, sym_block, sym_i, sym_j, ij_irr, p, orb_i, orb_j, i, j, k, sym_k, sym_l, q, kl_irr
       integer :: orb_k_it, orb_k, orb_l, orb_kl, l, kl, orb_l_it, ijkl_type, offset, orb_ij, orb_j_start, orb_j_end, err
       logical :: case_iiii, case_iijj, case_ijij, case_ijkl 
       integer, allocatable :: tmp(:)

       ij_offset = 0
       ij = 0
       sym_block = 0
       do sym_i=1,mobas%no_irr
          do sym_j=1,sym_i
 
             ij_irr = abel_prod_tab(sym_i,sym_j) !IRR of the product of symmetries of orbitals i,j
             p = sym_i*(sym_i-1)/2 + sym_j       !the index of the symmetry block for the pair of symmetries corresponding to the pair (sym_i,sym_j), sym_i .g.e sym_j
 
             !We loop over those (i,j) orbital indices that have been assigned to this thread for the current combination of symmetries (sym_i,sym_j).
             do orb_i=ij_orbital_range(1,p),ij_orbital_range(3,p)
 
                orb_j_start = 1
                if (sym_i .eq. sym_j) then 
                   orb_j_end = orb_i          !both orbitals are from the same symmetry and therefore the loop over the second orbital must be only from 1 to orb_i.
                else
                   orb_j_end = mobas%orbital_data(sym_j)%number_of_functions !both orbitals come from different symmetries so the loop over the second loop must be over all orbitals in that symmetry.
                endif
 
                if (orb_i .eq. ij_orbital_range(1,p)) orb_j_start = ij_orbital_range(2,p) !my first 'j' orbital index
                if (orb_i .eq. ij_orbital_range(3,p)) orb_j_end = ij_orbital_range(4,p)   !my last 'j' orbital index
                
                do orb_j=orb_j_start,orb_j_end
 
                   i = max(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_i
                   j = min(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_j
                   ij = ipair(i) + j

                   do sym_k=1,sym_i
                      do sym_l=1,sym_k
                         q = sym_k*(sym_k-1)/2 + sym_l
       
                         if (q > p) cycle !symmetry of the integral (ij|O|kl) => we want only unique pairs of pairs of symmetries (p,q) with p=(sym_i,sym_j),q=(sym_k,sym_l)
                         kl_irr = abel_prod_tab(sym_k,sym_l) !IRR of the product of symmetries of orbitals k,l
       
                         !symmetry restriction assuming that the operator O in (ij|O|kl) is totally symmetric; hence the symmetry restriction is given only by the IRRs of the MOs.
                         if (kl_irr .ne. ij_irr) cycle
       
                         sym_block = p*(p-1)/2 + q !index of the unique quartet of symmetries (sym_i sym_j|sym_k sym_l) which defines the symmetry block
                         offset = mobas%block_offset(sym_block) !the total number of integrals preceeding the current symmetry block; this is the base for the index of the quartet of MOs
       
                         case_iiii = .false. 
                         case_iijj = .false.
                         case_ijij = .false.
                         case_ijkl = .false.
       
                         !determine the combination of symmetries we are dealing with; this is used below to ensure only unique combinations of orbitals are produced
                         if (sym_i .eq. sym_j .and. sym_k .eq. sym_l .and. sym_i .eq. sym_k) then
                            case_iiii = .true. !(II|II)
                         else if (sym_i .eq. sym_j .and. sym_k .eq. sym_l) then
                            case_iijj = .true. !(II|JJ)
                         else if (sym_i .eq. sym_k .and. sym_j .eq. sym_l) then
                            case_ijij = .true. !(IJ|IJ)
                         else
                            case_ijkl = .true. !(IJ|KL)
                         endif
       
                         !Compute indices of the (IJ| or (II| bra pairs so that we can compare them against the indices of the ket pairs and make sure orb_ij .ge. orb_kl.
                         !Case (II|: the index of a unique pair of orbitals is a simple triangularization
                         !Case (IJ|: the index of a pair of orbitals (each from a different symmetry) is the same as the sequence number of the element (orb_j,orb_i) in a linearized 2D array.
                         if (case_iiii .or. case_iijj) orb_ij = ipair(orb_i) + orb_j
                         if (case_ijij) orb_ij = orb_j + mobas%orbital_data(sym_j)%number_of_functions*(orb_i-1)
 
                         if (case_iiii .or. case_ijij) orb_k_it = orb_i           !we need (ij|kl) with i .ge. k
                         if (case_iijj .or. case_ijkl) orb_k_it = mobas%orbital_data(sym_k)%number_of_functions !k comes from a symmetry different to i so we have to loop over all orbitals in that symmetry.
 
                         do orb_k=1,orb_k_it
 
                            if (case_iiii .or. case_iijj) orb_l_it = orb_k          !we need (ij|kl) with k .ge. l
                            if (case_ijij .or. case_ijkl) orb_l_it = mobas%orbital_data(sym_l)%number_of_functions !l comes from a symmetry different to k so we have to loop over all orbitals in that symmetry.
 
                            do orb_l=1,orb_l_it
 
                               k = mobas%relative_to_absolute(orb_k,sym_k) !overall index of the orb_k orbital
                               l = mobas%relative_to_absolute(orb_l,sym_l) !overall index of the orb_l orbital
                               kl = ipair(max(k,l)) + min(k,l)
 
                               if (.not.(two_p_continuum)) then
                                  !If integrals for only 1p in the continuum are required then we skip integrals of the type (TC|CC), (CT|CC) and (CC|CC).
                                  ijkl_type = ij_type(ij)+ij_type(kl)
                                  if (ijkl_type .eq. 5 .or. ijkl_type .eq. 6) cycle
                               endif
 
                               !For (II|II) and (IJ|IJ) cases we need to make sure that ij .ge. kl: compute indices of the |IJ) or |II) ket pairs and compare them to the bra indices.
                               if (case_iiii) then
                                  orb_kl = ipair(orb_k) + orb_l
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               if (case_ijij) then 
                                  orb_kl = orb_l + mobas%orbital_data(sym_l)%number_of_functions*(orb_k-1)
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               !Order in which I will loop over the ijkl integrals: total number of integrals for each ij pair.
                               ij_offset(ij) = ij_offset(ij) + 1
 
                            enddo !orb_l
                         enddo !orb_k
 
                      enddo !sym_l
                   enddo !sym_k
                enddo !orb_j
             enddo !orb_i
          enddo !sym_j
       enddo !sym_i

       call move_alloc(ij_offset,tmp)
       allocate(ij_offset(size(tmp)),stat=err)
       if (err .ne. 0) call xermsg ('molecular_basis_mod','generate_ij_offset','Memory allocation failed.',err, 1)

       ij_offset = 0
       do ij=2,size(tmp)
          !Total number of integrals preceeding the current ij pair:
          ij_offset(ij) = ij_offset(ij-1) + tmp(ij-1)
       enddo

       !Total number of integrals generated by this thread
       n_integrals = ij_offset(size(tmp)) + tmp(size(tmp))

    end subroutine generate_ij_offset

    subroutine extract_non_zero_cf(cf_t,cf_t_non_zero,mo_indices,n_non_zero)
       implicit none
       real(kind=cfp), allocatable :: cf_t(:,:), cf_t_non_zero(:,:)
       integer, allocatable :: mo_indices(:,:), n_non_zero(:)

       integer :: i, j, p, n_mo, n_ao, err

          n_mo = size(cf_t,1)
          n_ao = size(cf_t,2)

          if (allocated(cf_t_non_zero)) deallocate(cf_t_non_zero)
          if (allocated(n_non_zero)) deallocate(n_non_zero)
          if (allocated(mo_indices)) deallocate(mo_indices)

          allocate(cf_t_non_zero(n_mo,n_ao),mo_indices(n_mo,n_ao),n_non_zero(n_ao),stat=err)
          if (err .ne. 0) call xermsg ('molecular_basis_mod','extract_non_zero_cf','Memory allocation failed.',err, 1)

          n_non_zero = 0
          cf_t_non_zero = 0.0_cfp
          mo_indices = 0
          do p=1,n_ao
             j = 0
             do i=1,n_mo
                if (cf_t(i,p) .ne. 0.0_cfp) then
                   j = j + 1
                   n_non_zero(p) = n_non_zero(p) + 1
                   cf_t_non_zero(j,p) = cf_t(i,p)
                   mo_indices(j,p) = i
                endif
             enddo !p
          enddo !i

    end subroutine extract_non_zero_cf
 
    !> The routine assumes that the array for output (ijkl_integrals) has been zeroed-out before the algorithm starts and that the threads have been launched outside of this routine.
    !> \warning Note the use of the 'allocatable' attributes for some argument arrays - this strongly affects performance as explained in the comment for omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local.
    subroutine omp_two_p_transform_ijrs_block_to_ijkl(ijrs,ijks,ijks_t,cf,cf_t_non_zero,mo_indices, &
                n_cf_t_non_zero,no_ao,mobas,rs_ind,ij_type,ao_is_local,tol,&
              &ijkl_integrals,ijkl_indices,int_type,rs_start,rs_end,no_pairs,ij_orbital_range,two_p_continuum,thread_id,ij_offset)
       use const_gbl, only: abel_prod_tab
       use mpi_gbl
       use omp_lib
       use special_functions_gbl, only: ipair
       use sort_gbl, only: sort_int_float
       implicit none
       real(kind=cfp), allocatable :: ijrs(:,:), cf(:,:), cf_t_non_zero(:,:), ijks(:,:), ijks_t(:,:), ijkl_integrals(:,:,:)
       real(kind=cfp) :: tol
       integer, allocatable :: rs_ind(:,:), ij_orbital_range(:,:), ijkl_indices(:,:), ij_offset(:), &
                                mo_indices(:,:), n_cf_t_non_zero(:)
       integer, intent(in) :: int_type, rs_start, rs_end, no_ao, no_pairs, thread_id
       integer(kind=1) :: ij_type(:)
       class(molecular_orbital_basis_obj) :: mobas
       logical, intent(in) :: ao_is_local, two_p_continuum
 
       integer :: ij, sym_block, sym_i, sym_j, ij_irr, p, orb_i, orb_j, i, j, rs, r, s, k, sym_k, sym_l, q, &
                    kl_irr, orb_ij, orb_j_start, orb_j_end, a,b,t
       integer :: cnt, orb_k_it, orb_k, orb_l, orb_kl, l, kl, orb_l_it, rs_relative, ijrs_type, ijkl_type
       real(kind=cfp) :: mo_int, threshold
       logical :: case_iiii, case_iijj, case_ijij, case_ijkl, check_small

       !If we will be calculating contributions to the transformed integrals of the last set of rs-indices then we can delete the final integrals smaller than a given threshold.
       !However, this can only be done if all AO integrals are kept by each process - otherwise the accumulated transformed integrals come only from contributions of those AO integrals
       !that are kept by this process. This prevents me from deleting the MO integrals here - this can only be done once all MO integrals have been accumulated from all processes. 
       if (rs_end .eq. no_pairs .and. ao_is_local) then
          check_small = .true.
       else
          check_small = .false.
       endif

       !The value 'threshold' is used to neglect small ijrs integrals. The use of this threshold approximately halves the compute time for the
       !second step (at least for the case: pyrazine 6-311+G_dp/all HF orbitals). Note that we can use the 'tol' value only if this process keeps all
       !AO integrals (ao_is_local .eq. true.)! If the AO integrals are scattered among all processes then ijrs(:,:) contains only PARTIAL contributions to the full
       !integral so we cannot neglect these partial contributions on the 'tol' level here.
       if (ao_is_local) then
          threshold = tol
       else
          threshold = 0.0_cfp
       endif
 
       !II. THE SECOND PART OF THE INTEGRAL TRANSFORM:
       !We loop over all unique symmetry blocks and in each symmetry block we loop over all unique combinations of orbitals. We precompute, for each symmetry block, the ij part of the full integral 
       !index for (ij|kl). Only the kl part is computed in the inner-most loop. Note that the loops over the orbital indices orb_k and orb_l start with 1. These orb_ indices correspond to the internal
       !orbital indices as defined in molecular_orbital_basis_obj. These indices always start with 1 no matter what the actual, i.e. external, indices of the orbitals are. In other words the algorithm works
       !even for cases in which the external indices of the orbitals don't start with 1.
       ij = 0
       sym_block = 0
       do sym_i=1,mobas%no_irr
          do sym_j=1,sym_i
 
             ij_irr = abel_prod_tab(sym_i,sym_j) !IRR of the product of symmetries of orbitals i,j
             p = sym_i*(sym_i-1)/2 + sym_j       !the index of the symmetry block for the pair of symmetries corresponding to the pair (sym_i,sym_j), sym_i .g.e sym_j
 
             !We loop over those (i,j) orbital indices that have been assigned to this thread for the current combination of symmetries (sym_i,sym_j).
             do orb_i=ij_orbital_range(1,p),ij_orbital_range(3,p)
 
                orb_j_start = 1
                if (sym_i .eq. sym_j) then 
                   orb_j_end = orb_i          !both orbitals are from the same symmetry and therefore the loop over the second orbital must be only from 1 to orb_i.
                else
                   orb_j_end = mobas%orbital_data(sym_j)%number_of_functions !both orbitals come from different symmetries so the loop over the second loop must be over all orbitals in that symmetry.
                endif
 
                if (orb_i .eq. ij_orbital_range(1,p)) orb_j_start = ij_orbital_range(2,p) !my first 'j' orbital index
                if (orb_i .eq. ij_orbital_range(3,p)) orb_j_end = ij_orbital_range(4,p)   !my last 'j' orbital index
                
                do orb_j=orb_j_start,orb_j_end
 
                   i = max(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_i
                   j = min(mobas%relative_to_absolute(orb_i,sym_i),mobas%relative_to_absolute(orb_j,sym_j)) !overall index of the orbital orb_j
                   ij = ipair(i) + j
 
                   !3rd step - transform the third index: r->k
                   ijks(:,:) = 0.0_cfp
                   do rs=rs_start,rs_end

                      r = rs_ind(1,rs)
                      s = rs_ind(2,rs)
 
                      if (.not.(two_p_continuum)) then
                         !If integrals for only 1p in the continuum are required then we skip half-transformed integrals of the type (TC|CC], (CT|CC] and (CC|CC].
                         !For explanation of ijrs_type see omp_two_p_transform_pqrs_block_to_ijrs_AO_is_local
                         ijrs_type = ij_type(ij) + mobas%ao_basis%get_function_pair_type(r,s)
                         if (ijrs_type .eq. 5 .or. ijrs_type .eq. 6) cycle
                      endif
 
                      rs_relative = rs-rs_start +1 !the indexing in ijrs is relative to rs_start
                      !r = rs_ind(1,rs)
                      !s = rs_ind(2,rs)
                      mo_int = ijrs(ij,rs_relative) != (ij|O|rs]
                      !write(100+myrank,'(2i10,e25.15)') ij,rs,mo_int
 
                      !The encompassing if statement makes efficient use of the branch predictor so it should use 
                      !only a negligible amount of compute time.
                      if (ao_is_local) then
                         if (abs(mo_int) .le. threshold) cycle
                      endif
 
                      !todo once AO symmetry has been implemented I can check first if a given (ij|rs) integral, ijrs(rs,ij), contributes to the (ij|ks) integrals.
                      if (s .eq. r) then
!                         do k=mobas%so2mo_range(1,r), mobas%so2mo_range(2,r)
!                            ijks(k,s) = ijks(k,s) + mo_int*cf_t(k,r)
!                         enddo !k
                          do k=1,n_cf_t_non_zero(r)
                             a = mo_indices(k,r)
                             ijks(a,s) = ijks(a,s) + mo_int*cf_t_non_zero(k,r)
                          enddo !k
                      else
!                         do k=mobas%so2mo_range(1,r), mobas%so2mo_range(2,r)
!                            ijks(k,s) = ijks(k,s) + mo_int*cf_t(k,r)
!                         enddo !k
                         do k=1,n_cf_t_non_zero(r)
                             a = mo_indices(k,r)
                             ijks(a,s) = ijks(a,s) + mo_int*cf_t_non_zero(k,r)
                         enddo !k
!                         do k=mobas%so2mo_range(1,s), mobas%so2mo_range(2,s)
!                            ijks(k,r) = ijks(k,r) + mo_int*cf_t(k,s)
!                         enddo !k
                         do k=1,n_cf_t_non_zero(s)
                            a = mo_indices(k,s)
                            ijks(a,r) = ijks(a,r) + mo_int*cf_t_non_zero(k,s)
                         enddo !k
                      endif
                   enddo !rs
 
                   !4th step - transform the fourth index: s->l
                   call transpose_2d(ijks,ijks_t,mobas%number_of_functions,no_ao) !ijks_t = ijks**T: this pays off in cache locality (speed up) below
 
                   !Total number of integrals preceeding the current ij pair:
                   cnt = ij_offset(ij)
 
                   do sym_k=1,sym_i
                      do sym_l=1,sym_k
                         q = sym_k*(sym_k-1)/2 + sym_l
       
                         if (q > p) cycle !symmetry of the integral (ij|O|kl) => we want only unique pairs of pairs of symmetries (p,q) with p=(sym_i,sym_j),q=(sym_k,sym_l)
                         kl_irr = abel_prod_tab(sym_k,sym_l) !IRR of the product of symmetries of orbitals k,l
       
                         !symmetry restriction assuming that the operator O in (ij|O|kl) is totally symmetric; hence the symmetry restriction is given only by the IRRs of the MOs.
                         if (kl_irr .ne. ij_irr) cycle
       
                         sym_block = p*(p-1)/2 + q !index of the unique quartet of symmetries (sym_i sym_j|sym_k sym_l) which defines the symmetry block
       
                         case_iiii = .false. 
                         case_iijj = .false.
                         case_ijij = .false.
                         case_ijkl = .false.
       
                         !determine the combination of symmetries we are dealing with; this is used below to ensure only unique combinations of orbitals are produced
                         if (sym_i .eq. sym_j .and. sym_k .eq. sym_l .and. sym_i .eq. sym_k) then
                            case_iiii = .true. !(II|II)
                            orb_ij = ipair(orb_i) + orb_j
                            orb_k_it = orb_i           !we need (ij|kl) with i .ge. k
                         else if (sym_i .eq. sym_j .and. sym_k .eq. sym_l) then
                            case_iijj = .true. !(II|JJ)
                            orb_k_it = mobas%orbital_data(sym_k)%number_of_functions !k comes from a symmetry different to i so we have to loop over all orbitals in that symmetry.
                         else if (sym_i .eq. sym_k .and. sym_j .eq. sym_l) then
                            case_ijij = .true. !(IJ|IJ)
                            orb_ij = orb_j + mobas%orbital_data(sym_j)%number_of_functions*(orb_i-1)
                            orb_k_it = orb_i           !we need (ij|kl) with i .ge. k
                         else
                            case_ijkl = .true. !(IJ|KL)
                            orb_k_it = mobas%orbital_data(sym_k)%number_of_functions !k comes from a symmetry different to i so we have to loop over all orbitals in that symmetry.
                         endif
 
                         do orb_k=1,orb_k_it
 
                            if (case_iiii .or. case_iijj) orb_l_it = orb_k          !we need (ij|kl) with k .ge. l
                            if (case_ijij .or. case_ijkl) orb_l_it = mobas%orbital_data(sym_l)%number_of_functions !l comes from a symmetry different to k so we have to loop over all orbitals in that symmetry.
 
                            do orb_l=1,orb_l_it
 
                               k = mobas%relative_to_absolute(orb_k,sym_k) !overall index of the orb_k orbital
                               l = mobas%relative_to_absolute(orb_l,sym_l) !overall index of the orb_l orbital
                               kl = ipair(max(k,l)) + min(k,l)
 
                               if (.not.(two_p_continuum)) then
                                  !If integrals for only 1p in the continuum are required then we skip integrals of the type (TC|CC), (CT|CC) and (CC|CC).
                                  ijkl_type = ij_type(ij)+ij_type(kl)
                                  if (ijkl_type .eq. 5 .or. ijkl_type .eq. 6) cycle
                               endif
 
                               !For (II|II) and (IJ|IJ) cases we need to make sure that ij .ge. kl: compute indices of the |IJ) or |II) ket pairs and compare them to the bra indices.
                               if (case_iiii) then
                                  orb_kl = ipair(orb_k) + orb_l
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               if (case_ijij) then 
                                  orb_kl = orb_l + mobas%orbital_data(sym_l)%number_of_functions*(orb_k-1)
                                  if (orb_kl > orb_ij) cycle
                               endif
 
                               !In the other cases we just compute indices of the ket pairs.
                               if (case_iijj) orb_kl = ipair(orb_k) + orb_l
                               if (case_ijkl) orb_kl = orb_l+(orb_k-1)*mobas%orbital_data(sym_l)%number_of_functions
 
                               !transform the fourth index: s->l
                               !mo_int = (ij|kl), where i,j,k,l are the molecular orbitals.
                               a = mobas%mo2so_range(1,l); b = mobas%mo2so_range(2,l)
                               mo_int = sum(ijks_t(a:b,k)*cf(a:b,l))
                               !do s=mobas%mo2so_range(1,l), mobas%mo2so_range(2,l)
                               !   mo_int = mo_int + ijks_t(s,k)*cf(s,l)
                               !enddo !s
 
                               !INTEGRAL INDEXING FOLLOWS: this must be consistent with the index function
 
                               a = max(ij,kl)
                               b = min(ij,kl)
                               t = ipair(a) + b

                               cnt = cnt + 1
                               ijkl_indices(cnt,thread_id+1) = t
                               ijkl_integrals(cnt,thread_id+1,int_type) = ijkl_integrals(cnt,thread_id+1,int_type) + mo_int

                               !Delete final MO integrals smaller than threshold. Split the if statement into two to make efficient
                               !use of the branch predictor.
                               if (check_small) then
                                  if (abs(ijkl_integrals(cnt,thread_id+1,int_type)) < tol) then
                                     ijkl_integrals(cnt,thread_id+1,int_type) = 0.0_cfp
                                  end if
                               endif
          
                            enddo !orb_l
                         enddo !orb_k
 
                      enddo !sym_l
                   enddo !sym_k
                enddo !orb_j
             enddo !orb_i
          enddo !sym_j
       enddo !sym_i

    end subroutine omp_two_p_transform_ijrs_block_to_ijkl

    subroutine bisect_index(val,ijkl_indices,col,last_index,ind,found)
       implicit none
       integer, intent(in) :: val, last_index, col
       integer, pointer     :: ijkl_indices(:,:)
       integer, intent(out) :: ind
       logical, intent(out) :: found

       integer :: hlf, A, B, i, j

          found = .false.
          if (last_index .eq. 0) return

          if (last_index .le. 10) then
             do i=1,last_index
                if (ijkl_indices(i,col) .eq. val) then
                   found = .true.
                   ind = i
                   exit
                endif
             enddo !i
             return
          endif
          if (val > ijkl_indices(last_index,col) .or. val < ijkl_indices(1,col)) return

!          write(100+col,*) 'in bisect',val,last_index
!          write(100+col,*) 'max',ijkl_indices(1,col),ijkl_indices(last_index,col) 
!          do i=1,last_index
!             write(100+col,*) 'c',i,ijkl_indices(i,col)
!          enddo
          A = 1
          B = last_index
          i = 0
          do
             i = i + 1
             if (B .le. A+10) then
                do j=A,B
                   if (ijkl_indices(j,col) .eq. val) then
                      found = .true.
                      ind = j
                      exit
                   endif
                enddo !j
                exit
             endif
             hlf = A+(B-A)/2
!             write(100+col,'(3i5,3i10)') A,hlf,B,ijkl_indices(A,col),ijkl_indices(hlf,col),ijkl_indices(B,col)
             if (ijkl_indices(hlf,col) .eq. val) then
                found = .true.
                ind = hlf
                exit
             endif
             if (val < ijkl_indices(hlf,col)) then
                B = hlf-1
             elseif (val > ijkl_indices(hlf,col)) then
                A = hlf+1
             endif
!             if (ijkl_indices(A,col) .eq. val) then
!                found = .true.
!                ind = A
!                exit
!             endif
!             if (ijkl_indices(B,col) .eq. val) then
!                found = .true.
!                ind = B
!                exit
!             endif
!             if (B .eq. A+1) exit !further division is not possible: the index is not there
          enddo

!          write(100+col,*) 'end',found,val,i
                 
    end subroutine bisect_index

    !> Note that whatever change is made here to the format of the ijkl file
    !> must be replicated explicitly in finalize_two_electron_integrals_sparse for the
    !> case of MPI parallel calculation.
    subroutine write_ijkl_indices(this,lunit,record_start,position_after_write)
       use mpi_gbl
       implicit none
       class(molecular_orbital_basis_obj) :: this
       integer, intent(in) :: lunit, record_start
       integer, intent(out) :: position_after_write

       integer :: err
       
         write(level3,'("--------->","molecular_orbital_basis_obj:write_ijkl_indices")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'write_ijkl_indices', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if

         call mpi_mod_barrier(err)

         if (.not. associated(this % ijkl_indices)) then
            call xermsg ('molecular_orbital_basis_obj', 'write_ijkl_indices', &
                         'The ijkl_indices array has not been allocated.', 2, 1)
         end if

         if (myrank .eq. master) then
            write(lunit,pos=record_start,err=10) this%ind_ijkl_integral, size(this%ijkl_indices,2)
            write(lunit,err=10) this%ijkl_indices
            inquire(lunit,pos=position_after_write)
         endif

         !get the position_after_write on all processes
         call mpi_mod_bcast(position_after_write,master)

         write(level3,'("<---------","molecular_orbital_basis_obj:write_ijkl_indices")')

         return

 10      call mpi_xermsg ('molecular_orbital_basis_obj', 'write_ijkl_indices', &
                          'Error executing the write command while writing the array data to the disk.', 2, 1)

    end subroutine write_ijkl_indices

    subroutine read_ijkl_indices(this,lunit,file_name,record_start,position_after_read)
       use iso_c_binding, only: c_associated
       implicit none
       class(molecular_orbital_basis_obj) :: this
       integer, intent(in) :: lunit, record_start
       integer, intent(out) :: position_after_read
       integer :: err, d2, i,j
       integer(mpiint)  :: local_master_array(nprocs)
       character(len=*) :: file_name
       integer(int64) :: offset, length

         write(level3,'("--------->","molecular_orbital_basis_obj:read_ijkl_indices")')

         if (.not. this % initialized) then
            call xermsg ('molecular_orbital_basis_obj', 'read_ijkl_indices', &
                         'The object has not been initialized or not all orbitals have been added.', 1, 1)
         end if
       
         call mpi_mod_barrier(err)

         if (this%shared_window_ijkl_indices /= -1) then
            d2 = size(this%ijkl_indices,2)
            call mpi_memory_deallocate_integer_2dim(this%ijkl_indices,this%ind_ijkl_integral*d2,this%shared_window_ijkl_indices)
            this%shared_window_ijkl_indices = -1
            !this%ijkl_indices => null()
         endif

         if (map .and. c_associated(this % mmap % ptr)) then
            call this % mmap % finalize
         end if

         write(level3,*) record_start,local_rank,lunit

         if (myrank .eq. master) then
            read(lunit,pos=record_start,err=50) this%ind_ijkl_integral, d2
         endif
         call mpi_mod_bcast(d2,master)
         call mpi_mod_bcast(this%ind_ijkl_integral,master)

         write(level3,*) this%ind_ijkl_integral,d2 

         !All local masters allocate memory for the index array: in non-shared mode every process is the local master.
         if (.not. map) then
            this%shared_window_ijkl_indices = mpi_memory_allocate_integer_2dim(this%ijkl_indices,this%ind_ijkl_integral,d2)
         end if
         
         call mpi_memory_synchronize(this%shared_window_ijkl_indices)
         length = this%ind_ijkl_integral * d2 * bit_size(longint)/8

         if (myrank .eq. master) then
            inquire(lunit,pos=position_after_read)
            offset = position_after_read - 1
            if (.not. map) read(lunit,err=50) this%ijkl_indices
            position_after_read = position_after_read + length
         endif

         call mpi_memory_synchronize(this%shared_window_ijkl_indices)

         !get the position_after_read on all processes
         call mpi_mod_bcast(position_after_read,master)
         call mpi_mod_bcast(d2,master)

         !If memory mapping is enabled, let each process map the file independently
         if (map) then
            call mpi_mod_bcast(offset, master)
            call this % mmap % init(file_name, offset, length)
            call c_f_pointer(this % mmap % ptr, this % ijkl_indices, [this%ind_ijkl_integral, d2])
         end if

         if (shared_enabled) then !Every node has only one copy of the indexing array

            local_master_array = 1
            
            call mpi_mod_allgather(local_rank,local_master_array)
            
            call mpi_memory_synchronize(this%shared_window_ijkl_indices)

            if (myrank .eq. master .and. .not. map) then
               do i=1,nprocs
                  if (i-1 .eq. master) cycle
                  if (local_master_array(i) == local_master) then 
                     write(level2,"('Sending array to ',i4)") (i-1)
                     do j=1,d2
                        call mpi_mod_send(int(i-1,mpiint),this%ijkl_indices(:,j),1,this%ind_ijkl_integral)
                     enddo   
                  endif
               enddo
            endif
            
            if (local_rank == local_master .and. myrank /= master .and. .not. map) then
               do j=1,d2
                  write(level2,"('Gathering info my array is size ',i12)") size(this%ijkl_indices,1)
                  call mpi_mod_recv(master,1,this%ijkl_indices(:,j),this%ind_ijkl_integral)
               enddo
            endif
    
            call mpi_memory_synchronize(this%shared_window_ijkl_indices)

         else !every MPI task keeps a copy of the indexing array

            call mpi_mod_barrier(err)

            do i=1,d2
                if (.not. map) call mpi_mod_bcast(this%ijkl_indices(1:this%ind_ijkl_integral,d2),master)
            enddo !i

         endif

         call mpi_mod_barrier(err)
 
         write(level3,'("<---------","molecular_orbital_basis_obj:read_ijkl_indices")')
 
         return
 
 50      call mpi_xermsg ('molecular_orbital_basis_obj', 'read_ijkl_indices', &
                          'Error executing the read command while reading the array data to the disk.', 2, 1)

    end subroutine read_ijkl_indices

    subroutine orbital_radial_charge_density(this,rmat_radius,A,B,delta_r,save_to_disk,charge_densities)
       use cgto_pw_expansions_gbl
       use phys_const_gbl, only: fourpi
       use gto_routines_gbl, only: cms_gto_norm
       use const_gbl, only: line_len, fmat
       implicit none
       class(molecular_orbital_basis_obj) :: this
       real(kind=cfp), intent(in) :: rmat_radius,A,B,delta_r
       logical, intent(in) :: save_to_disk
       real(kind=cfp), allocatable :: charge_densities(:,:)
 
       integer :: i, j, k, err, n, s1, CGTO1_M, s2, CGTO2_M, ao1, ao2, BA_ind, n_shell_pairs, no_orbitals,&
                    lu, number_of_cgto_shells, max_bspline_l, max_prop_l, p
       integer :: bto1_index, bto2_index, BTO_M1, BTO_M2, number_of_bto_shells, lb, lbmb, cnt
       type(CGTO_shell_data_obj), allocatable :: dummy_cgto(:)
       type(BTO_shell_data_obj), allocatable :: dummy_bto(:)
       type(CGTO_shell_pair_pw_expansion_obj), allocatable :: CGTO_shell_pair_pw_expansion(:)
       type(CGTO_shell_pw_expansion_obj), allocatable :: CGTO_shell_pw_expansion(:)
       type(pw_expansion_obj) :: grid
       real(kind=cfp), allocatable :: orb_cf(:,:), bto_amplitude(:,:), jacobian_r(:)
       real(kind=cfp) :: fac, cf, norm, r
       character(len=line_len) :: file_name
       integer, parameter :: der = 0
 
          write(level3,'("--------->","molecular_orbital_basis_obj:radial_charge_density")')
 
          if (.not. this % initialized) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', &
                          'The object has not been initialized or not all orbitals have been added.', 1, 1)
          end if
 
          if (A < 0.0_cfp) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input A < 0.', 2, 1)
          end if
          if (B <= 0.0_cfp) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input B .le. 0.', 3, 1)
          end if
          if (delta_r <= 0.0_cfp) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input delta_r .le. 0.', 4, 1)
          end if
          if (B < A) then
             call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', 'On input B < A.', 5, 1)
          end if

          !Print out distances of atoms from CMS: useful for analyzing the
          !radial charge density plots.
          write(level2,'(/,"Atom distances from CMS")')
          do i=1,this%symmetry_data%no_nuc
             r = sqrt(dot_product(this%symmetry_data%nucleus(i)%center,this%symmetry_data%nucleus(i)%center))
             call this%symmetry_data%nucleus(i)%print
             write(level2,'(5X,"Distance from CMS = ",e25.15," Bohr")') r
          enddo

          !Construct grid for [A,B],delta_r
          call grid%eval_regular_grid(A,B,delta_r)

          n = grid%n_total_points
          allocate(jacobian_r(n),stat=err)
          if (err /= 0) then
             call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 1 error.', err, 1)
          end if

          !write(stdout,'(/,"Radial grid:")')
          do i=1,n
             !write(stdout,'(e25.15)') grid%r_points(i)
             jacobian_r(i) = grid%r_points(i)**2
          enddo !i

          if (this%ao_basis%n_cgto_functions < this%ao_basis%number_of_functions) then
             call this%ao_basis%get_all_BTO_shells(dummy_bto,number_of_bto_shells)
             max_prop_l = 0
             max_bspline_l = maxval(dummy_bto(:)%l)
          else
             number_of_bto_shells = 0
             max_prop_l = 0
             max_bspline_l = -1
          endif

          if (this%ao_basis%n_cgto_functions > 0) then
 
             call this%ao_basis%get_all_CGTO_shells(dummy_cgto,number_of_cgto_shells)
   
             !Calculate radial charge_densities for all pairs of the CGTOs and multiply them in
             !with the orbital coefficients to obtain radial charge_densities of the orbitals.
             n_shell_pairs = number_of_cgto_shells*(number_of_cgto_shells+1)/2
             allocate(CGTO_shell_pair_pw_expansion(n_shell_pairs),CGTO_shell_pw_expansion(number_of_cgto_shells),stat=err)
             if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 3 error.', err, 1)
             end if
    
             do i=1,number_of_cgto_shells
                if (rmat_radius > 0.0_cfp .and. this%ao_basis%shell_descriptor(3,i) .eq. 1) then !Normalize the continuum to the R-matrix sphere
                   fac = cms_gto_norm(rmat_radius, dummy_cgto(i) % l, dummy_cgto(i) % number_of_primitives, &
                                      dummy_cgto(i) % exponents, dummy_cgto(i) % contractions, dummy_cgto(i) % norm, &
                                      dummy_cgto(i) % norms)
                   write(level2,'("Continuum normalization factor",i0,e25.15)') i,fac
                   dummy_cgto(i)%norm = dummy_cgto(i)%norm*fac
                elseif (this%ao_basis%shell_descriptor(3,i) .eq. 1 .and. rmat_radius .le. 0.0_cfp) then
                   write(level2,'("Continuum functions in this shell will not be normalized &
                                    &to the R-matrix sphere since rmat_radius .le. 0.0_cfp.")')
                endif
             enddo !i
    
             write(level2,'("Precalculating radial amplitudes for all pairs of shells of CGTOs...")')
   
             call init_CGTO_pw_expansions_mod(0,maxval(dummy_cgto(:)%l))
   
             cnt = 0
             do i=1,number_of_cgto_shells
                do j=1,i
                   k = i*(i-1)/2+j
                   call CGTO_shell_pair_pw_expansion(k)%init_CGTO_shell_pair_pw_expansion(dummy_cgto(i),i,dummy_cgto(j),j)
                   call CGTO_shell_pair_pw_expansion(k)%assign_grid(grid%r_points,grid%weights)
                   call CGTO_shell_pair_pw_expansion(k)%eval_CGTO_shell_pair_pw_expansion
                   fac = k/real(n_shell_pairs)*100
                   if (mod(nint(fac),5) .eq. 0 .and. nint(fac) .ne. cnt) then
                      write(level2,'(f8.3,"% done")') real(fac)
                      cnt = nint(fac)
                   endif
                enddo !j
             enddo !i
   
             if (number_of_cgto_shells < this%ao_basis%number_of_shells) then
                call init_CGTO_pw_expansions_mod(max_bspline_l,maxval(dummy_cgto(:)%l))
                do i=1,number_of_cgto_shells
                   call CGTO_shell_pw_expansion(i)%init_CGTO_shell_pw_expansion(dummy_cgto(i),i)
                   call CGTO_shell_pw_expansion(i)%assign_grid(grid%r_points,grid%weights)
                   call CGTO_shell_pw_expansion(i) % eval_CGTO_shell_pw_expansion (dummy_bto(1) % bspline_grid % knots, &
                                                                                   max_bspline_l, max_prop_l, 0)
                enddo !i
             endif
   
             write(level2,'("...done")')
          else
             number_of_cgto_shells = 0
          endif

          if (number_of_cgto_shells < this%ao_basis%number_of_shells) then
             !Calculate amplitudes of the radial parts of the BTOs: B(r)/r
             allocate(bto_amplitude(n,dummy_bto(1)%bspline_grid%n),stat=err)
             if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 5 error.', err, 1)
             end if
             bto_amplitude = 0.0_cfp
             norm = 1.0_cfp
             do bto1_index=1,dummy_bto(1)%bspline_grid%n
                do j=1,n
                   bto_amplitude(j,bto1_index) = dummy_bto(1) % bspline_grid % bspline_amplitude(grid % r_points(j), &
                                                                                    norm, bto1_index, der) / grid % r_points(j)
                enddo !j
             enddo !bto1_index
          else
             number_of_bto_shells = 0
          endif

          !Copy the orbital coefficients to one array: this relies on the fact that the molecular orbitals are indexed symmetry by symmetry.
          call this%get_orbital_coefficient_matrix(orb_cf)
    
          !For each orbital calculate the radial charge_densities from the radial amplitudes of the GTOs.
          no_orbitals = size(orb_cf,2)
          if (allocated(charge_densities)) deallocate(charge_densities)
          allocate(charge_densities(n,no_orbitals),stat=err)
          if (err /= 0) then
             call xermsg ('molecular_orbital_basis_obj','orbital_radial_charge_density', 'Memory allocation 5 error.', err, 1)
          end if
    
          charge_densities = 0.0_cfp
          do i=1,no_orbitals
             !Loop over all unique pairs of the CGTO shells
             ao1 = 0
             do s1=1,number_of_cgto_shells
                do CGTO1_M=1,dummy_cgto(s1)%number_of_functions
                   ao1 = ao1 + 1
                   ao2 = 0
                   !CGTO/CGTO contribution
                   do s2=1,s1
                      k = s1*(s1-1)/2+s2
                      do CGTO2_M=1,dummy_cgto(s2)%number_of_functions
                         ao2 = ao2 + 1

                         if (ao2 > ao1) cycle

                         BA_ind = CGTO2_M+dummy_cgto(s2)%number_of_functions*(CGTO1_M-1)
                         if (CGTO_shell_pair_pw_expansion(k)%neglect_m_lm(BA_ind,1)) cycle

                         fac = 2.0_cfp !off diagonal terms contribute twice since we loop only over the unique pairs of CGTOs
                         if (ao2 .eq. ao1) fac = 1.0_cfp !diagonal terms contribute only once
                         !the amplitudes of the CGTO pairs correspond to
                         !projections on X_{0,0} = 1/sqrt(4*pi) so we have to
                         !get rid of the pi-related factor.
                         cf = fac*orb_cf(ao1,i)*orb_cf(ao2,i)*sqrt(fourpi)
                         if (cf .eq. 0.0_cfp) cycle

                         charge_densities(1:n,i) = charge_densities(1:n,i) &
                                                + cf * CGTO_shell_pair_pw_expansion(k) % angular_integrals(1:n,BA_ind,1)

                      enddo !CGTO2_M
                   enddo !s2

                   !CGTO1/BTO contribution
                   fac = 2.0_cfp !off diagonal terms contribute twice since we loop only over the unique pairs of CGTO/BTO
                   ao2 = this%ao_basis%n_cgto_functions
                   do s2=1,number_of_bto_shells
                      lb = dummy_bto(s2)%l
                      bto1_index = dummy_bto(s2)%bspline_index
                      do BTO_M1 = -lb,lb
                         lbmb = lb*lb+lb+BTO_M1+1
                         ao2 = ao2 + 1

                         !ind = CGTO_shell_pw_expansion(s1)%non_neg_indices(CGTO1_M,1,lbmb)
                         !if (ind .eq. 0) cycle

                         !the amplitudes of the CGTOs correspond to projections
                         !on X_{l,m} of the BTO so there is no pi-related factor
                         !to get rid of but we have to multiply-in the BTO norm.
                         cf = fac*orb_cf(ao1,i)*orb_cf(ao2,i)*dummy_bto(s2)%norm
                         if (cf .eq. 0.0_cfp) cycle
                     
                         p = CGTO_shell_pw_expansion(s1)%non_neg_indices_l(CGTO1_M,lbmb)
                         if (p .eq. 0) cycle

                         charge_densities(1:n,i) = charge_densities(1:n,i) &
                                    + cf*CGTO_shell_pw_expansion(s1)%angular_integrals(1:n,p)*bto_amplitude(1:n,bto1_index)

                      enddo !BTO_M1
                   enddo !s2

                enddo !CGTO1_M
             enddo !s1

             !BTO/BTO contribution
             ao1 = this%ao_basis%n_cgto_functions
             do s1=1,number_of_bto_shells
                bto1_index = dummy_bto(s1)%bspline_index
                do BTO_M1=-dummy_bto(s1)%l,dummy_bto(s1)%l
                   ao1 = ao1 + 1
                   ao2 = this%ao_basis%n_cgto_functions
                   do s2=1,s1
                      bto2_index = dummy_bto(s2)%bspline_index
                      do BTO_M2=-dummy_bto(s2)%l,dummy_bto(s2)%l
                         ao2 = ao2 + 1
                         if (ao2 > ao1) cycle

                         !The only angular integrals which are non-zero are
                         !those where the angular parts of the BTOs are the same.
                         if ((dummy_bto(s1)%l .ne. dummy_bto(s2)%l) .or. (BTO_M1 .ne. BTO_M2)) cycle

                         fac = 2.0_cfp !off diagonal terms contribute twice since we loop only over the unique pairs of CGTOs
                         if (ao2 .eq. ao1) fac = 1.0_cfp !diagonal terms contribute only once
                         cf = fac*orb_cf(ao1,i)*orb_cf(ao2,i)*dummy_bto(s1)%norm*dummy_bto(s2)%norm
                         if (cf .eq. 0.0_cfp) cycle

                         charge_densities(1:n,i) = charge_densities(1:n,i) &
                                    + cf*bto_amplitude(1:n,bto1_index)*bto_amplitude(1:n,bto2_index)

                      enddo !BTO_M2
                   enddo !s2
                enddo !BTO_M1
             enddo !s1
 
             if (save_to_disk) then
                !Save the orbital density into a file with name: orb_rad_den_num.sym
                write(file_name,'(i6.6,".",i1)') this%absolute_to_relative(1,i),this%absolute_to_relative(2,i)
                if (this%is_continuum(i)) then
                   file_name = trim("orb_rad_den_continuum_"//adjustl(file_name))
                else
                   file_name = trim("orb_rad_den_target_"//adjustl(file_name))
                endif
                open(file=file_name,newunit=lu,status='replace',form=fmat,iostat=err)
                if (err /= 0) then
                    call xermsg ('molecular_orbital_basis_obj', 'orbital_radial_charge_density', &
                                 'Error opening the file for orbital charge density.', 6, 1)
                end if
 
                write(lu,'("#Radial charge density for orbital: ",3(i0,1x))') &
                            i, this%absolute_to_relative(1,i),this%absolute_to_relative(2,i)
                charge_densities(1:n,i) = charge_densities(1:n,i)*jacobian_r(1:n)
                fac = sum(charge_densities(1:n,i)*grid%weights(1:n))
                write(lu,'("#Integral over the charge density: ",e25.15)') fac
                do j=1,n
                   write(lu,'(2e28.14e4)') grid%r_points(j),charge_densities(j,i)
                enddo !j
                close(lu)
                write(level2,'("Radial charge density for orbital ",i0," has been written to file: ",a)') &
                            i, trim(adjustl(file_name))
                write(level2,'("Integral over the charge density: ",e25.15)') fac
             endif
          enddo !i
 
    end subroutine orbital_radial_charge_density

    subroutine eval_orbital(this,orb_i,r,n_points,orbital_at_r,sign_at_r)
       use cgto_pw_expansions_gbl
       use gto_routines_gbl, only: cms_gto_norm
       use const_gbl, only: line_len, fmat
       implicit none
       class(molecular_orbital_basis_obj) :: this
       integer, intent(in) :: n_points, orb_i
       real(kind=cfp), intent(in) :: r(3,n_points)
       real(kind=cfp), allocatable :: orbital_at_r(:)
       integer, allocatable :: sign_at_r(:)

       integer :: i, j, n, irr, err
       real(kind=cfp), allocatable :: ao_basis_at_r(:,:)

          write(level3,'("--------->","molecular_orbital_basis_obj:eval_orbital")')
 
          if (.not. this % initialized) then
             call xermsg ('molecular_orbital_basis_obj', 'eval_orbital', &
                          'The object has not been initialized or not all orbitals have been added.', 1, 1)
          end if

          call this%ao_basis%eval_basis(r,n_points,ao_basis_at_r)

          if (allocated(orbital_at_r)) deallocate(orbital_at_r)
          if (allocated(sign_at_r)) deallocate(sign_at_r)
          allocate(orbital_at_r(n_points),sign_at_r(n_points),stat=err)
          if (err .ne. 0) call xermsg ('molecular_orbital_basis_obj','eval_orbital', 'Memory allocation 2 failed.', err, 1)

          !Contract the orbital coefficients with the values of the AO functions at the specified points.
          j = this%absolute_to_relative(1,orb_i)
          irr = this%absolute_to_relative(2,orb_i)
          n = this%ao_basis%number_of_functions
          !$omp parallel do
          do i=1,n_points
             orbital_at_r(i) = sum(this%orbital_data(irr)%coefficients(1:n,j)*ao_basis_at_r(1:n,i))
             sign_at_r(i) = nint(sign(1.0_cfp,orbital_at_r(i)))
          enddo !i
          !$omp end parallel do

          write(level3,'("<---------","molecular_orbital_basis_obj:eval_orbital")')

    end subroutine eval_orbital

    subroutine construct_canonical_continuum(this, atomic_1el_integral_storage, atomic_2el_integral_storage, &
                            molecular_integral_storage, integral_options_in, molecular_2el_algorithm, n_cont_bto, n_cont_cgto, &
                            mo_integrals_file_name, can_mo_integrals_file_name)
       use blas_lapack_gbl, only : syev, blasint
       use const_gbl, only: one_p_sym_ints, two_p_sym_ints, one_elham, two_el_ints, ijkl_indices_header, one_electron_ints, &
                            Fock_blocks_header, orbs_line
       implicit none
       class(molecular_orbital_basis_obj) :: this
       class(integral_storage_obj), target, intent(inout) :: atomic_1el_integral_storage, atomic_2el_integral_storage, &
                                                             molecular_integral_storage
       class(integral_options_obj), intent(in) :: integral_options_in
       integer, intent(in) :: molecular_2el_algorithm, n_cont_bto, n_cont_cgto
       character(len=line_len), intent(in) :: mo_integrals_file_name, can_mo_integrals_file_name

       type(integral_storage_obj), target :: tgt_storage
       type(p2d_array_obj), target :: one_electron_integrals, two_electron_integrals, ao_integrals
       type(data_header_obj) :: header
       type(integral_options_obj) :: integral_options, integral_options_dummy
       type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
       
       integer :: err, i, j, irr, lunit, pos, n_target, n_continuum, k, m
       integer :: begin_target_index, end_target_index, begin_continuum_index, end_continuum_index, first_record, last_record
       integer(blasint) :: lwork, info, n, lda
       logical, parameter :: tgt_is_local = .true.
       real(kind=cfp), allocatable :: F_CC_irr_block(:,:), w(:), work(:), cont_MO_wrt_AO(:,:), F_TC_irr_block(:,:), &
                                      F_TT_irr_block(:,:)
       logical, parameter :: VIRTUALS = .false. !testing
       integer :: relative_index !testing
 
          write(level3,'("--------->","molecular_orbital_basis_obj:construct_canonical_continuum")')

          integral_options = integral_options_in

          !1. read-in 1-el and 2-el integrals: ukrmol_interface.f90; also read-in ijkl_indices: ukrmol_interface.f90
          
          write(level1,'("Constructing canonical continuum")')
    
          !set up where to read-in the 1-electron integrals :: into one_electron_integrals
          err = tgt_storage%init(memory=one_electron_integrals)
          tgt_storage%integral_file%identifier = molecular_integral_storage%integral_file%identifier
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error initializing tgt_storage.',err,1)

          !we need to get the full header corresponding to the transformed overlap and kinetic energy integrals.
          err = molecular_integral_storage % integral_file % get_header_containing(header, this % get_basis_name(), &
                                                                                   one_p_sym_ints)
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error locating the transformed 1-electron &
                                      &integrals; see data_header_obj%get_header_containing for details.',err,1)
          !the actual reading-in
          call tgt_storage%read(molecular_integral_storage,header%name,integral_options_dummy,tgt_is_local)
          call tgt_storage%final

          write(level2,'(/,"One electron integrals read-in.")')

          !set up where to read-in the 2-electron integrals :: into two_electron_integrals
          err = tgt_storage%init(memory=two_electron_integrals)
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error initializing tgt_storage.',err,1)

          !we need to get the full header corresponding to the transformed 2-electron integrals.
          err = molecular_integral_storage % integral_file % get_header_containing(header, this % get_basis_name(), &
                                                                                   two_p_sym_ints)
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error locating the transformed 2-electron &
                                      &integrals; see data_header_obj%get_header_containing for details.',err,1)
          !the actual reading-in
          call tgt_storage%read(molecular_integral_storage,header%name,integral_options_dummy,tgt_is_local)
          call tgt_storage%final 

          write(level2,'(/,"Two electron integrals read-in.")')

          !we need to get the full header corresponding to the ijkl_indices
          err = molecular_integral_storage % integral_file % get_header_containing(header, this % get_basis_name(), &
                                                                                   ijkl_indices_header)

          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error locating the ijkl_indices; see &
                                      &data_header_obj%get_header_containing for details.',err,1)

          !we are reading the ijkl_indices into this instance of molecular_orbital_basis_obj
          lunit = molecular_integral_storage%integral_file%get_unit_no()
          call this%read_ijkl_indices(lunit,mo_integrals_file_name,header%first_record,pos)

          !TEST
          call this%print_orbitals()

          !2. construct F_CC of Fock matrix
          do irr = 1,this%no_irr
             write(level1,'(/,"Diagonalizing block of Fock matrix, irr = ",i0,", number of functions = ",i0)') irr, &
                   this%orbital_data(irr)%number_of_functions 
             if (this%orbital_data(irr)%number_of_functions > 0) then  
                
                if (.not. this%is_continuum(this%relative_to_absolute(1,irr))) then
                   begin_target_index=this%relative_to_absolute(1,irr)
                   end_target_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                   do i=1,this%orbital_data(irr)%number_of_functions
                      if (end_target_index == this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr) &
                          .and. this%is_continuum(this%relative_to_absolute(i,irr))) then
                         end_target_index = this%relative_to_absolute(i,irr) - 1
                      endif
                   enddo
                   n_target = end_target_index - begin_target_index + 1 
                   call this%construct_fock_matrix_block(one_electron_integrals,two_electron_integrals,begin_target_index, &
                                                         end_target_index,begin_target_index,end_target_index,F_TT_irr_block)
                    
                   write(level3,'(/,2X,"Target-target block of Fock matrix")')
                   k = begin_target_index - 1
                   do i=1,n_target/orbs_line
                      write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
                      write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1, &
                                                          k+orbs_line)
                      do j=1,n_target
                         write(level3,'(i0,50e25.15)') j+begin_target_index-1,F_TT_irr_block(j,k-begin_target_index+2: &
                                                                                             k+orbs_line-begin_target_index+1)
                      enddo
                      k = k + orbs_line
                   enddo
   
                   m = mod(n_target,orbs_line)
                   if (m > 0) then
                   write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
                   write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1,k+m)
                      do j=1,n_target
                         write(level3,'(i0,50e25.15)') j+begin_target_index-1,F_TT_irr_block(j,k-begin_target_index+2: &
                                                                                             k+m-begin_target_index+1)
                      enddo
                   endif
        
                   if (allocated(F_TT_irr_block)) deallocate(F_TT_irr_block)   
                endif
 
                begin_continuum_index = -1
                end_continuum_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                do i=1,this%orbital_data(irr)%number_of_functions
                   if (begin_continuum_index == -1 .and. this%is_continuum(this%relative_to_absolute(i,irr))) then
                      begin_continuum_index = this%relative_to_absolute(i,irr)
                   endif
                enddo

                if (.not. this%is_continuum(this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr))) then
                   begin_continuum_index = end_continuum_index + 1
                endif

                n_continuum = end_continuum_index - begin_continuum_index + 1

                if (integral_options%canonize_virtuals_instead_of_continuum) then
                   ! continuum -> virtuals
                   begin_continuum_index = -1
                   end_continuum_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                   do i=1,this%orbital_data(irr)%number_of_functions
                      if (begin_continuum_index == -1 .and. this%orbital_data(irr)%occup(i) == 0.0) then
                         begin_continuum_index = this%relative_to_absolute(i,irr)
                      endif
                   enddo

                   if (this%orbital_data(irr)%occup(this%orbital_data(irr)%number_of_functions) > 0.0) then
                      begin_continuum_index = end_continuum_index + 1
                   endif

                   n_continuum = end_continuum_index - begin_continuum_index + 1
                endif

                if (n_continuum > 0) then

                   call this%construct_fock_matrix_block(one_electron_integrals,two_electron_integrals,begin_continuum_index, &
                                                         end_continuum_index,begin_continuum_index,end_continuum_index, &
                                                         F_CC_irr_block)

                   ! if only constructing Fock matrix blocks, then here we load the diagonal elements of the continuum-continuum
                   ! block into the orbital data energies and then we cycle to skip the diagonalisation of the block   
                   ! if canonizing virtuals instead of continuum, it overwrites the real virtuals' energies
                   if (integral_options%only_construct_fock_blocks) then
                      do i=1,n_continuum
                         relative_index = begin_continuum_index - this%relative_to_absolute(1,irr)
                         this%orbital_data(irr)%energy(i+relative_index) = F_CC_irr_block(i,i)
                      enddo  
                      cycle
                   endif

                   write(level3,'(/,2X,"Continuum-continuum block of Fock matrix")')
                   k = begin_continuum_index - 1
                   do i=1,n_continuum/orbs_line
                      write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
                      write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1, &
                                                          k+orbs_line)
                      do j=1,n_continuum
                         write(level3,'(i0,50e25.15)') j+begin_continuum_index-1, F_CC_irr_block(j,k-begin_continuum_index+2: &
                                                                                              k+orbs_line-begin_continuum_index+1)
                      enddo
                      k = k + orbs_line
                   enddo
   
                   m = mod(n_continuum,orbs_line)
                   if (m > 0) then
                   write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
                   write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1,k+m)
                      do j=1,n_continuum
                         write(level3,'(i0,50e25.15)') j+begin_continuum_index-1, F_CC_irr_block(j,k-begin_continuum_index+2: &
                                                                                                 k+m-begin_continuum_index+1)
                      enddo
                   endif
                   
                   !3. Diagonalize F_CC: call syev
                   n = int(n_continuum,blasint)
                   lda = int(n_continuum,blasint)
                   allocate(w(n_continuum))
                   allocate(work(1))
                   lwork=-1
                   call syev('N','U',n,F_CC_irr_block,lda,w,work,lwork,info)
             
                   if (info == 0) then
                      lwork=int(work(1))
                      deallocate(work)
                      allocate(work(lwork))
                      call syev('V','U',n,F_CC_irr_block,lda,w,work,lwork,info)
                   endif
             
                   if (info .ne. 0) call xermsg('molecular_orbital_basis_obj','construct_canonical_continuum', &
                                                'Diagonalisation of the Fock matrix block failed',1,1)
                 
                   !4. obtain continuum MO coefficients in original AO basis.
                   allocate(cont_MO_wrt_AO(this%orbital_data(irr)%number_of_coefficients,n_continuum))
                   cont_MO_wrt_AO = matmul(this%orbital_data(irr)%coefficients(:,this%absolute_to_relative(1,begin_continuum_index &
                                           ):this%absolute_to_relative(1,end_continuum_index)),F_CC_irr_block(:,:))

                   do i=1,n_continuum
                      relative_index = begin_continuum_index - this%relative_to_absolute(1,irr)
                      this%orbital_data(irr)%coefficients(:,i+relative_index) = cont_MO_wrt_AO(:,i)
                      this%orbital_data(irr)%energy(i+relative_index) = w(i)
                   enddo

                   !Making sure that all MPI tasks are working with the same energies and orbitals
                   call mpi_mod_bcast(this%orbital_data(irr)%coefficients,int(master,kind=mpiint))
                   call mpi_mod_bcast(this%orbital_data(irr)%energy,int(master,kind=mpiint))
                   
                   deallocate(w,work,cont_MO_wrt_AO)

                   if (allocated(F_CC_irr_block)) deallocate(F_CC_irr_block)
                endif
 
                !TEST - overlaps
                call this%check_overlaps(one_electron_integrals,this%relative_to_absolute(1,irr),this%relative_to_absolute( &
                                                                    this%orbital_data(irr)%number_of_functions,irr))

             endif
          enddo !irr

          !call this%delete_small_coefficients

          this%ijkl_indices => null()   !to restore initial state (when two_electron_integrals_sparse is used)

          !6. store 2-el integrals in a new file

          call molecular_integral_storage%final() !closing moints, because below we read from it

          call atomic_orbital_basis%read(mo_integrals_file_name) !we have to read this, because the calculation of 2-el ints in scatci_integrals  
                                                      !destroys the data
     
          call atomic_orbital_basis%write(can_mo_integrals_file_name)
          call this%write(can_mo_integrals_file_name)            

          !Describe where the transformed AO->MO integrals will be stored 
          err = molecular_integral_storage%init(disk=can_mo_integrals_file_name)
          if (err /= 0) then
             call mpi_xermsg('molecular_orbital_basis_obj', 'construct_canonical_continuum', 'error initializing the target &
                             &molecular_integral_storage', err, 1)
          end if

          integral_options % keep_ao_integrals = .false.

          if (integral_options_in%keep_ao_integrals) then
             write(level1,'("Using precalculated AO integrals")')
             
             !
             ! TRANSFORM THE 1-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE MOLECULAR ORBITALS:
             !
                     
             this%ao_integral_storage => atomic_1el_integral_storage !point to the storage for the atomic integrals
             call this%one_electron_integrals(molecular_integral_storage,integral_options)

             !
             ! TRANSFORM THE 2-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE MOLECULAR ORBITALS:
             !
 
             this%ao_integral_storage => atomic_2el_integral_storage
             if (molecular_2el_algorithm == 1 .or. (molecular_2el_algorithm == 0 .and. n_cont_bto > 0 .and. n_cont_cgto == 0)) then
                call this % two_electron_integrals_sparse(molecular_integral_storage, integral_options)
             else
                call this % two_electron_integrals(molecular_integral_storage, integral_options)
             end if

          else
             write(level1,'("Recalculating the AO integrals.")')
                             
             !
             ! CALCULATE THE ATOMIC 1-ELECTRON INTEGRALS:
             !
            
             !describe where the AO integrals will be stored 
             err = atomic_1el_integral_storage%init(memory=ao_integrals) 
             if (err /= 0) then
                call mpi_xermsg('molecular_basis', 'construct_canonical_continuum', 'error initializing the target &
                                                                                                &atomic_integral_storage', err, 1)
             endif

             call atomic_orbital_basis%one_electron_integrals(atomic_1el_integral_storage,integral_options)

             !
             ! TRANSFORM THE 1-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE MOLECULAR ORBITALS:
             !
           
             this%ao_integral_storage => atomic_1el_integral_storage !point to the storage for the atomic integrals
 
             call this%one_electron_integrals(molecular_integral_storage,integral_options)

             !
             ! CALCULATE THE ATOMIC 2-ELECTRON INTEGRALS
             !
            
             !describe where the AO integrals will be stored 
             err = atomic_2el_integral_storage%init(memory=ao_integrals) 
             if (err /= 0) then
                call mpi_xermsg('molecular_basis', 'construct_canonical_continuum', 'error initializing the target &
                                                                                                &atomic_integral_storage', err, 1)
             endif

             integral_options%calculate_two_el_ints = .true.
             this%ao_integral_storage => atomic_2el_integral_storage
             call atomic_orbital_basis%two_electron_integrals(atomic_2el_integral_storage,integral_options)

             !
             ! TRANSFORM THE 2-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE CANONICAL MOLECULAR ORBITALS:
             !

             if (molecular_2el_algorithm == 1 .or. (molecular_2el_algorithm == 0 .and. n_cont_bto > 0 .and. n_cont_cgto == 0)) then
                call this % two_electron_integrals_sparse(molecular_integral_storage, integral_options)
             else
                call this % two_electron_integrals(molecular_integral_storage, integral_options)
             end if            
          endif

          write(level1,'("Canonical continuum integrals saved into ",a)') can_mo_integrals_file_name
          
          !  
          !RELOADING THE INTEGRALS BACK TO MEMORY AND CONSTRUCTION OF F_TC BLOCKS
          !

          !set up where to read-in the 1-electron integrals :: into one_electron_integrals
          err = tgt_storage%init(memory=one_electron_integrals)
          tgt_storage%integral_file%identifier = molecular_integral_storage%integral_file%identifier
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error initializing tgt_storage.',err,1)

          !we need to get the full header corresponding to the transformed overlap and kinetic energy integrals.
          err = molecular_integral_storage % integral_file % get_header_containing(header, this % get_basis_name(), &
                                                                                   one_p_sym_ints)
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error locating the transformed 1-electron &
                                      &integrals; see data_header_obj%get_header_containing for details.',err,1)
          !the actual reading-in
          call tgt_storage%read(molecular_integral_storage,header%name,integral_options_dummy,tgt_is_local)
          call tgt_storage%final

          write(level2,'(/,"One electron integrals read-in.")')

          !set up where to read-in the 2-electron integrals :: into two_electron_integrals
          err = tgt_storage%init(memory=two_electron_integrals)
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error initializing tgt_storage.',err,1)

          !we need to get the full header corresponding to the transformed 2-electron integrals.
          err = molecular_integral_storage % integral_file % get_header_containing(header, this % get_basis_name(), &
                                                                                   two_p_sym_ints)
          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error locating the transformed 2-electron &
                                      &integrals; see data_header_obj%get_header_containing for details.',err,1)
          !the actual reading-in
          call tgt_storage%read(molecular_integral_storage,header%name,integral_options_dummy,tgt_is_local)
          call tgt_storage%final 

          write(level2,'(/,"Two electron integrals read-in.")')

          !we need to get the full header corresponding to the ijkl_indices
          err = molecular_integral_storage % integral_file % get_header_containing(header, this % get_basis_name(), &
                                                                                   ijkl_indices_header)

          if (err .ne. 0) call xermsg('molecular_basis','construct_canonical_continuum','Error locating the ijkl_indices; see &
                                      &data_header_obj%get_header_containing for details.',err,1)

          !we are reading the ijkl_indices into this instance of molecular_orbital_basis_obj
          lunit = molecular_integral_storage%integral_file%get_unit_no()
          call this%read_ijkl_indices(lunit,mo_integrals_file_name,header%first_record,pos)

          !call this%Fock_core_test(one_electron_integrals,two_electron_integrals)
          do irr=1,this%no_irr
             write(level3,'(/,"Diagonal block of Fock matrix in basis of canonical continuum functions, irr = ",i0, &
                   &", number of functions = ",i0)') irr,this%orbital_data(irr)%number_of_functions 

             if (this%orbital_data(irr)%number_of_functions > 0) then               
                begin_continuum_index = -1
                end_continuum_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                do i=1,this%orbital_data(irr)%number_of_functions
                   if (begin_continuum_index == -1 .and. this%is_continuum(this%relative_to_absolute(i,irr))) then
                      begin_continuum_index = this%relative_to_absolute(i,irr)
                   endif
                enddo

                if (.not. this%is_continuum(this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr))) then
                   begin_continuum_index = end_continuum_index + 1
                endif

                n_continuum = end_continuum_index - begin_continuum_index + 1
                
                if (n_continuum > 0) then

                   call this%construct_fock_matrix_block(one_electron_integrals,two_electron_integrals,begin_continuum_index, &
                                                     end_continuum_index,begin_continuum_index,end_continuum_index,F_CC_irr_block)
               
                   write(level3,'(/,2X,"Continuum-continuum block of Fock matrix")')
                   k = begin_continuum_index - 1
                   do i=1,n_continuum/orbs_line
                      write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
                      write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1, &
                                                       k+orbs_line)
                      do j=1,n_continuum
                         write(level3,'(i0,50e25.15)') j+begin_continuum_index-1, F_CC_irr_block(j,k-begin_continuum_index+2: &
                                                                                              k+orbs_line-begin_continuum_index+1)
                      enddo
                      k = k + orbs_line
                   enddo
   
                   m = mod(n_continuum,orbs_line)
                   if (m > 0) then
                      write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
                      write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1,k+m)
                      do j=1,n_continuum
                         write(level3,'(i0,50e25.15)') j+begin_continuum_index-1, F_CC_irr_block(j,k-begin_continuum_index+2: &
                                                                                                 k+m-begin_continuum_index+1)
                      enddo
                   endif
                   if (allocated(F_CC_irr_block)) deallocate(F_CC_irr_block)
                endif
       
                if (this%orbital_data(irr)%number_of_functions - n_continuum > 0) then
                   
                   begin_target_index=this%relative_to_absolute(1,irr)
                   end_target_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                   do i=1,this%orbital_data(irr)%number_of_functions
                      if (end_target_index == this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr) &
                          .and. this%is_continuum(this%relative_to_absolute(i,irr))) then
                         end_target_index = this%relative_to_absolute(i,irr) - 1
                      endif
                   enddo
                   n_target = end_target_index - begin_target_index + 1 
     
                   call this%construct_fock_matrix_block(one_electron_integrals,two_electron_integrals,begin_target_index, &
                                                         end_target_index,begin_target_index,end_target_index,F_TT_irr_block)
                   
                   write(level3,'(/,2X,"Target-target block of Fock matrix")')
                   k = begin_target_index - 1
                   do i=1,n_target/orbs_line
                      write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
                      write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1, &
                                                          k+orbs_line)
                      do j=1,n_target
                         write(level3,'(i0,50e25.15)') j+begin_target_index-1,F_TT_irr_block(j,k-begin_target_index+2: &
                                                                                             k+orbs_line-begin_target_index+1)
                      enddo
                      k = k + orbs_line
                   enddo
   
                   m = mod(n_target,orbs_line)
                   if (m > 0) then
                   write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
                   write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1,k+m)
                      do j=1,n_target
                         write(level3,'(i0,50e25.15)') j+begin_target_index-1,F_TT_irr_block(j,k-begin_target_index+2: &
                                                                                             k+m-begin_target_index+1)
                      enddo
                   endif

                   if (allocated(F_TT_irr_block)) deallocate(F_TT_irr_block)
                endif
                call this%check_overlaps(one_electron_integrals,this%relative_to_absolute(1,irr),this%relative_to_absolute( &
                                                                    this%orbital_data(irr)%number_of_functions,irr))
             endif
          enddo

          write(level3,'(/,2X,"Orbitals where continuum is canonical follow")') 
          call this%print_orbitals()

          !at this point, molecular_integral_storage has to be associated with can_mo_integrals_file_name file and open
          lunit = molecular_integral_storage%integral_file%get_unit_no()        !unit that is associated to the file opened
          first_record = molecular_integral_storage%integral_file%start_record(Fock_blocks_header)
          last_record = first_record

          do irr=1,this%no_irr
             if (this%orbital_data(irr)%number_of_functions > 0) then
                ! construct block CT
                begin_target_index = this%relative_to_absolute(1,irr)
                begin_continuum_index = -1
                end_continuum_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                do i=1,this%orbital_data(irr)%number_of_functions
                   if (begin_continuum_index == -1 .and. this%is_continuum(this%relative_to_absolute(i,irr))) then
                      begin_continuum_index = this%relative_to_absolute(i,irr)
                      end_target_index = begin_continuum_index - 1
                   endif
                enddo

                if (.not. this%is_continuum(this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr))) then
                   begin_continuum_index = end_continuum_index + 1
                   end_target_index = end_continuum_index
                endif

                n_target = end_target_index - begin_target_index + 1
                n_continuum = end_continuum_index - begin_continuum_index + 1

                ! the below does not have to be here if we perform canonization (the continuum flags are updated)
                ! but it is neccessary, if we only construct the blocks
                if (integral_options%canonize_virtuals_instead_of_continuum) then
                   begin_target_index = this%relative_to_absolute(1,irr)
                   begin_continuum_index = -1
                   end_continuum_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                   do i=1,this%orbital_data(irr)%number_of_functions
                      if (begin_continuum_index == -1 .and. this%orbital_data(irr)%occup(i) == 0.0) then
                         begin_continuum_index = this%relative_to_absolute(i,irr)
                         end_target_index = begin_continuum_index - 1
                      endif
                   enddo

                   if (this%orbital_data(irr)%occup(this%orbital_data(irr)%number_of_functions) > 0.0) then
                      begin_continuum_index = end_continuum_index + 1
                      end_target_index = end_continuum_index
                   endif

                   n_target = end_target_index - begin_target_index + 1
                   n_continuum = end_continuum_index - begin_continuum_index + 1
                endif

                if (n_target > 0) then
                   call this%construct_fock_matrix_block(one_electron_integrals,two_electron_integrals,begin_target_index, &
                                                         end_target_index, begin_target_index,end_target_index,F_TT_irr_block)
                endif

                if (n_target > 0 .and. n_continuum > 0) then
                   call this%construct_fock_matrix_block(one_electron_integrals,two_electron_integrals,begin_target_index, &
                                                        end_target_index,begin_continuum_index,end_continuum_index,F_TC_irr_block)
                endif

                if (n_continuum > 0) then
                   call this%construct_fock_matrix_block(one_electron_integrals,two_electron_integrals,begin_continuum_index, &
                                                      end_continuum_index,begin_continuum_index,end_continuum_index,F_CC_irr_block)
                endif
             else
                n_target = 0
                n_continuum = 0
             endif

             if (myrank .eq. master) then
                write(lunit,pos=last_record) n_target, n_continuum
                inquire(lunit,pos=last_record)
                if (n_target > 0) then
                   write(lunit,pos=last_record) F_TT_irr_block
                   inquire(lunit,pos=last_record)
                endif
                if (n_target > 0 .and. n_continuum > 0) then
                   write(lunit,pos=last_record) F_TC_irr_block
                   inquire(lunit,pos=last_record)
                endif
                if (n_continuum > 0) then
                   write(lunit,pos=last_record) F_CC_irr_block
                   inquire(lunit,pos=last_record)
                endif
             endif
          enddo !irr

          if (allocated(F_TT_irr_block)) deallocate(F_TT_irr_block)
          if (allocated(F_TC_irr_block)) deallocate(F_TC_irr_block)
          if (allocated(F_CC_irr_block)) deallocate(F_CC_irr_block)

          call mpi_mod_bcast(last_record,int(master,kind=mpiint))

          call molecular_integral_storage%integral_file%close_record(Fock_blocks_header,first_record,last_record)
          write(level1,'("Target-target and target-continuum blocks saved into ",a)') &
                molecular_integral_storage%integral_file%get_file_name()
         
          write(level3,'("<---------","molecular_orbital_basis_obj:construct_canonical_continuum")')
 
    end subroutine construct_canonical_continuum

    subroutine construct_fock_matrix_block(this,one_electron_integrals,two_electron_integrals,begin_index_one,end_index_one,&
                                           begin_index_two,end_index_two,F_block)
       use const_gbl, only : one_elham, two_el_ints
       implicit none
       class(molecular_orbital_basis_obj) :: this
       class(p2d_array_obj),intent(in) :: one_electron_integrals,two_electron_integrals
       integer,intent(in) :: begin_index_one, end_index_one, begin_index_two, end_index_two
       real(kind=cfp), allocatable :: F_block(:,:)

       integer :: i, j, k, ind(1:1), two_ind(1:2,1), four_ind(1:4,1)
       real(kind=cfp) :: integral_1el, coulomb_integral, exchange_integral
       
          if(allocated(F_block)) deallocate(F_block)
          allocate(F_block(end_index_one-begin_index_one+1,end_index_two-begin_index_two+1))

          do i=begin_index_one,end_index_one
             do j=begin_index_two,end_index_two
                !construct F matrix element (i,j)
                !h-matrix
                two_ind(1,1) = i
                two_ind(2,1) = j
                ind(1:1) = this%integral_index(one_elham,two_ind,.false.) !1-el integral for orbitals a,b
                integral_1el = one_electron_integrals%a(ind(1),one_electron_integrals%find_column_matching_name(one_elham))
                      
                coulomb_integral=0.0
                exchange_integral=0.0
                do k=1,this%number_of_functions
                   if (.not. this%is_continuum(k)) then
                
                      !J-matrix
                      four_ind(1,1) = i
                      four_ind(2,1) = j
                      four_ind(3,1) = k
                      four_ind(4,1) = k

                      ind(1:1) = this%integral_index(two_el_ints,four_ind,.false.)
                      !what happens when two_electron_integrals is used instead of two_electron_integrals_sparse ??
     
                      if ((ind(1) > 0) .and.(this%orbital_data(this%absolute_to_relative(2,k))%occup( &
                                             this%absolute_to_relative(1,k)) > 0.0)) then
                         coulomb_integral = coulomb_integral + two_electron_integrals%a(ind(1), &
                                            two_electron_integrals%find_column_matching_name(two_el_ints))
                      endif
 
                      !K-matrix
                      four_ind(1,1) = i
                      four_ind(2,1) = k
                      four_ind(3,1) = j
                      four_ind(4,1) = k

                      ind(1:1) = this%integral_index(two_el_ints,four_ind,.false.)
                            
                      if ((ind(1) > 0) .and. (this%orbital_data(this%absolute_to_relative(2,k))%occup( &
                                              this%absolute_to_relative(1,k)) > 0.0)) then
                         exchange_integral = exchange_integral + two_electron_integrals%a(ind(1), &
                                             two_electron_integrals%find_column_matching_name(two_el_ints))
                         !we are summing for all occupied orbitals, which is valid for closed shell systems
                         !for open shell system it depends on spin, this is not implemented
                      endif
                   endif
                enddo
                      
                F_block(i-begin_index_one+1,j-begin_index_two+1)=integral_1el + 2*coulomb_integral - exchange_integral
                !write(*,'("F_CC element",2i4,e25.15)') a, b, F_CC_irr_block(a,b)
             enddo
          enddo

    end subroutine construct_fock_matrix_block

    subroutine check_overlaps(this,one_electron_integrals,begin_index,end_index)
       use const_gbl, only : overlap_ints, orbs_line
       implicit none
       class(molecular_orbital_basis_obj) :: this
       class(p2d_array_obj),intent(in) :: one_electron_integrals
       integer,intent(in) :: begin_index, end_index

       integer :: i, j, ind(1:1), two_ind(1:2,1), n_indices, k, m
       real(kind=cfp), allocatable :: overlaps(:,:)
       
          n_indices = end_index-begin_index+1
          allocate(overlaps(n_indices,n_indices))
          do i=begin_index,end_index
             do j=begin_index,end_index
                two_ind(1,1) = i
                two_ind(2,1) = j
                ind(1:1) = this%integral_index(overlap_ints,two_ind,.false.) 
                overlaps(i-begin_index+1,j-begin_index+1) = one_electron_integrals%a(ind(1), &
                                                        one_electron_integrals%find_column_matching_name(overlap_ints)) 
             enddo
          enddo

          write(level3,'(/,2X,"Checking molecular orbital overlaps")')
            k = 0
            do i=1,n_indices/orbs_line
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+begin_index,k+orbs_line+begin_index-1)
               do j=1,n_indices
                  write(level3,'(i0,50e25.15)') j+begin_index-1, overlaps(j,k+1:k+orbs_line) !this%orbital_data(symmetry)%coefficients(j,k+1:k+orbs_line)
               enddo
               k = k + orbs_line
            enddo
   
            m = mod(n_indices,orbs_line)
            if (m > 0) then
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+begin_index,k+m+begin_index-1)
               do j=1,n_indices
                  write(level3,'(i0,50e25.15)') j+begin_index-1, overlaps(j,k+1:k+m) !this%orbital_data(symmetry)%coefficients(j,k+1:k+m)
               enddo
            endif

          deallocate(overlaps)

    end subroutine check_overlaps

    subroutine solve_roothan_equations(this,integral_options,orbital_data)
       use blas_lapack_gbl, only : syev, blasint
       use const_gbl, only : orbs_line
       implicit none
       class(molecular_orbital_basis_obj) :: this
       class(integral_options_obj) :: integral_options
       type(orbital_data_obj), allocatable :: orbital_data(:), orbital_data_old(:)

       type(integral_options_obj) :: integral_options_local
       type(integral_storage_obj), target :: atomic_1el_integral_storage, atomic_2el_integral_storage, &
                                             molecular_1el_integral_storage, molecular_2el_integral_storage
       type(p2d_array_obj) :: ao_1el_integrals, ao_2el_integrals, mo_1el_integrals, mo_2el_integrals
       real(kind=cfp), allocatable :: energies_old(:), energies_new(:), F_TT_irr_block(:,:), w(:), work(:), MO_wrt_AO(:,:)
       integer :: err, irr, i, j, begin_index, end_index, iters, n_target, k, m
       integer(blasint) :: lwork, info, n, lda

          write(level3,'("--------->","molecular_orbital_basis_obj:solve_roothan_equations")')

          integral_options_local = integral_options
          integral_options_local%keep_ao_integrals = .true.

          write(level1,'("Solving Roothan equations")')
          write(level2,'("Energies convergence parameter ",e25.15,", maximum number of iterations ",i0)') &
                integral_options_local%HF_convergence, integral_options_local%HF_max_iters
            
          allocate(energies_old(this%number_of_functions))
          allocate(energies_new(this%number_of_functions))

          allocate(orbital_data_old(size(orbital_data)))
          orbital_data_old=this%orbital_data

          do irr=1,this%no_irr
             i = this%orbital_data(irr)%number_of_functions
             if (i == 0) cycle
             begin_index=this%relative_to_absolute(1,irr)
             end_index=this%relative_to_absolute(i,irr)
             energies_new(begin_index:end_index) = this%orbital_data(irr)%energy(1:i)
          enddo
          energies_old = 0.0            

          call this%print_energy_sorted_orbital_table
    
          !
          ! CALCULATE THE ATOMIC 1-ELECTRON INTEGRALS:
          !
            
          !describe where the AO integrals will be stored 
          err = atomic_1el_integral_storage%init(memory=ao_1el_integrals) 
          if (err /= 0) then
             call mpi_xermsg('molecular_basis', 'solve_roothan_equations', 'error initializing the target &
                                                                                             &atomic_integral_storage', err, 1)
          endif

          call this%ao_basis%one_electron_integrals(atomic_1el_integral_storage,integral_options_local)

          !
          ! CALCULATE THE ATOMIC 2-ELECTRON INTEGRALS
          !
            
          !describe where the AO integrals will be stored 
          err = atomic_2el_integral_storage%init(memory=ao_2el_integrals) 
          if (err /= 0) then
             call mpi_xermsg('molecular_basis', 'solve_roothan_equations', 'error initializing the target &
                                                                                             &atomic_integral_storage', err, 1)
          endif

          integral_options_local%calculate_two_el_ints = .true.
          !this%ao_integral_storage => atomic_2el_integral_storage
          call this%ao_basis%two_electron_integrals(atomic_2el_integral_storage,integral_options_local)

          !describe where to store MO integrals
          err = molecular_1el_integral_storage%init(memory=mo_1el_integrals)
          if (err .ne. 0) then
             call mpi_xermsg('molecular_basis', 'solve_roothan_equations', 'error initializing the molecular_integral_storage', &
                             err, 1)
          endif
          err = molecular_2el_integral_storage%init(memory=mo_2el_integrals)
          if (err .ne. 0) then
             call mpi_xermsg('molecular_basis', 'solve_roothan_equations', 'error initializing the molecular_integral_storage', &
                             err, 1)
          endif

          iters = 0
          do while (maxval(abs(energies_old(:)-energies_new(:)))>integral_options_local%HF_convergence .and. &
                    iters < integral_options_local%HF_max_iters)  
             energies_old = energies_new
             iters = iters + 1

             write(level3,'(/,2X,"Iteration ",i0)') iters

             !
             ! TRANSFORM THE 1-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE MOLECULAR ORBITALS:
             !
          
             integral_options_local%calculate_two_el_ints = .false. 
             this%ao_integral_storage => atomic_1el_integral_storage !point to the storage for the atomic integrals
 
             call this%one_electron_integrals(molecular_1el_integral_storage,integral_options_local)

             !
             ! TRANSFORM THE 2-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE CANONICAL MOLECULAR ORBITALS:
             !

             integral_options_local%calculate_two_el_ints = .true.
             this%ao_integral_storage => atomic_2el_integral_storage
             call this % two_electron_integrals(molecular_2el_integral_storage, integral_options_local)

             do irr=1,this%no_irr
                
                if (this%orbital_data(irr)%number_of_functions > 0) then

                   begin_index = this%relative_to_absolute(1,irr)
                   end_index = this%relative_to_absolute(this%orbital_data(irr)%number_of_functions,irr)
                   n_target = this%orbital_data(irr)%number_of_functions
                
                   call this%construct_fock_matrix_block(mo_1el_integrals,mo_2el_integrals,begin_index,end_index,begin_index, &
                                                         end_index,F_TT_irr_block)

                   write(level3,'(/,2X,"Target-target block of Fock matrix, irr = ",i4)') irr
                   k = begin_index - 1
                   do i=1,n_target/orbs_line
                      write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
                      write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1, &
                                                          k+orbs_line)
                      do j=1,n_target
                         write(level3,'(i0,50e25.15)') j+begin_index-1,F_TT_irr_block(j,k-begin_index+2:k+orbs_line-begin_index+1)
                      enddo
                      k = k + orbs_line
                   enddo
   
                   m = mod(n_target,orbs_line)
                   if (m > 0) then
                   write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
                   write(level3,'(1X,50(e25.15))') (this%orbital_data(irr)%energy(this%absolute_to_relative(1,j)),j=k+1,k+m)
                      do j=1,n_target
                         write(level3,'(i0,50e25.15)') j+begin_index-1,F_TT_irr_block(j,k-begin_index+2:k+m-begin_index+1)
                      enddo
                   endif

                   !diagonalize 
                   n = int(n_target,blasint)
                   lda = int(n_target,blasint)
                   allocate(w(n_target))
                   allocate(work(1))
                   lwork=-1
                   call syev('N','U',n,F_TT_irr_block,lda,w,work,lwork,info)
             
                   if (info == 0) then
                      lwork=int(work(1))
                      deallocate(work)
                      allocate(work(lwork))
                      call syev('V','U',n,F_TT_irr_block,lda,w,work,lwork,info)
                   endif
             
                   if (info .ne. 0) call xermsg('molecular_orbital_basis_obj','solve_roothan_equations', &
                                                'Diagonalisation of the Fock matrix block failed',1,1)
               
                   allocate(MO_wrt_AO(this%orbital_data(irr)%number_of_coefficients,n_target))
                   MO_wrt_AO = matmul(this%orbital_data(irr)%coefficients(:,this%absolute_to_relative(1,begin_index &
                                              ):this%absolute_to_relative(1,end_index)),F_TT_irr_block(:,:))
               
                   do i=1,n_target
                      this%orbital_data(irr)%coefficients(:,i) = MO_wrt_AO(:,i)
                      this%orbital_data(irr)%energy(i) = w(i)
                      energies_new(i+begin_index-1) = w(i)
                   enddo

                   !Making sure that all MPI tasks work with the same energies and coefficients
                   call mpi_mod_bcast(this%orbital_data(irr)%coefficients,int(master,kind=mpiint))
                   call mpi_mod_bcast(this%orbital_data(irr)%energy,int(master,kind=mpiint))

                   deallocate(MO_wrt_AO,w,work)
                endif
             enddo !irr

             write(level1,'(/,"Roothan equations convergence: ",i0," ",e25.15)') &
               iters,maxval(abs(energies_old(:)-energies_new(:)))

          enddo

          if (iters==integral_options_local%HF_max_iters) then
             write(level1,'(/,"Roothan equations did not converge to the requested threshold")')
          else
             write(level1,'(/,"Roothan equations converged after ",i0," iterations")') iters
          endif

          orbital_data=this%orbital_data
           
          deallocate(energies_old,energies_new,orbital_data_old)

          call this%print_energy_sorted_orbital_table
        
          write(level3,'("<---------","molecular_orbital_basis_obj:solve_roothan_equations")')

    end subroutine solve_roothan_equations


   !> \brief   Calculate 2-el Coulomb repulsion integrals using Poisson equation method
   !> \authors M Konvalinka, Z Masin, J Benda
   !> \date    2024
   !>
   !> Evaluate the two-particle Coulomb repulsion and exchange integrals. One pair of orbitals acts as a source term
   !> in the Poisson equation with the Bloch operator:
   !> \f[
   !>      (\nabla^2 + L) V = -4 \pi \phi_a \phi_b + L V.
   !> \f]
   !> The product of the two orbitals on the right-hand side is expanded in the molecular orbital basis (reusing
   !> triple overlaps calculated earlier), and the equation is solved by the R-matrix method in the same basis. The
   !> unknown potential in Bloch term on the right-hand side is represented by a known multipole expansion of the
   !> orbital pair. The remaining pair of orbitals is then included by the triple overlaps, again:
   !> \f[
   !>      I = \langle \phi_c \phi_d | V \rangle.
   !> \f]
   !>
   subroutine two_electron_integrals_poisson(this,mo_integral_storage,mo_integral_options)
      use omp_lib, only: omp_get_wtime
      use phys_const_gbl, only: to_ev, pi
      use blas_lapack_gbl, only: blasint, syev
      use const_gbl, only: level1, level2, level3, kinetic_ints, bbb_ints, property_ints, one_p_sym_ints
      use parallel_arrays_gbl, only: p2d_array_obj
      implicit none

      class(molecular_orbital_basis_obj), intent(inout) :: this
      class(integral_storage_obj), intent(inout) :: mo_integral_storage
      class(integral_options_obj), intent(in) :: mo_integral_options

      type(p2d_array_obj), pointer :: mo_integrals
      type(p2d_array_obj), target :: ao_integrals
      type(integral_storage_obj), target :: mo_integral_disk
      type(data_header_obj) :: header
      type(integral_options_obj) :: options
      logical, parameter :: tgt_is_local = .true.

      real(kind=cfp) :: rmatr
      real(kind=wp) :: t

      type(diag_ham_obj) :: diag_ham_all(this%no_irr)
      integer :: err, sym, i, j, k, p, q, add, isym, lam, mu
      integer :: tot_chan, nsym, no_mo, no_ao, d1, d2, no_blocks, max_poisson_l, max_property_l
      integer :: sym_start, sym_end, ind(1:1), col_kei, col_ov3, col_prp
      integer :: l_min, l_max, i_cont, j_cont, p_cont, q_cont
      integer :: n_tgt(this%no_irr), n_cnt(this%no_irr), n_tgt_total
      real(kind=cfp) :: dvout
      real(kind=cfp), allocatable :: cf(:,:), cft(:,:), ov3(:,:,:), dnc(:,:,:), qnc(:,:,:), anc(:,:,:)

      !> chan_lm_all: l,m and index of symmetry
      !> Allocated as chann_lm_all(3,tot_chan)
      !> m is stored in first row (1,:), l in 2nd row (2,:), index of symmetry in the 3rd row
      !> nchannels = (l_max+1)**2-(l-min)**2) = tot_chann ... total number of channels over all symmetries
      integer, allocatable :: chan_lm_all(:,:)
      real(kind=cfp), allocatable :: bamps_all(:,:)
      logical, parameter :: normalize_to_a = .true.

         write(level3,'(/,"--------->molecular_basis_mod:two_electron_integrals_poisson")')

         call mpi_mod_barrier(err)

         if (nprocs > 1) then
            call mpi_xermsg('molecular_basis_gbl', 'two_electron_integrals_poisson', &
                            'Distributed (MPI) mode not yet implemented for Poisson equation method.', 1, 1)
         end if

         if (mo_integral_storage%in_memory()) then
            mo_integrals => mo_integral_storage%integral
         else !on disk
            err = mo_integral_disk%init(memory=ao_integrals)
            if (err /= 0) then
               call mpi_xermsg('molecular_basis_mod', 'two_electron_integrals_poisson', &
                               'Molecular integrals file initialization failed.', err, 1)
            end if
            !we need to get the full header corresponding to the transformed 1-electron integrals.
            header%name = mo_integral_storage%contruct_header_string(this%get_basis_name(),one_p_sym_ints)
            !read the integrals from the disk into mo_integrals
            call mo_integral_disk%read(mo_integral_storage,header%name,options,tgt_is_local)
            mo_integrals => ao_integrals
         endif

         !Get the column number for the mo_integrals%a array where the KEI are stored.
         col_kei = mo_integrals%find_column_matching_name(kinetic_ints)


         !EVALUATE ALL BOUNDARY AMPLITUDESn_tgt(this%no_irr), n_cnt(this%no_irr)
         rmatr = mo_integral_options % a
         call this%calculate_amplitudes(rmatr,normalize_to_a,bamps_all,chan_lm_all)

         tot_chan = size(chan_lm_all,2) !total number of channels

         nsym = this%no_irr
         write(level2,'(/,"Number of symmetries (irr-reps): ",i8)') nsym
         write(level2,'(/,"Total number of channels: ",i8)') tot_chan
         no_mo = this%number_of_functions
         write(level2,'(/,"Total number of molecular orbitals: ",i8)') no_mo
         no_ao = this%ao_basis%number_of_functions
         write(level2,'(/,"Total number of atomic orbitals: ",i8)') no_ao

         n_tgt = 0
         n_cnt = 0
         do i=1,this%number_of_functions
            if (this%is_continuum(i)) then
               n_cnt(this%absolute_to_relative(2,i)) = n_cnt(this%absolute_to_relative(2,i)) + 1
            else
               n_tgt(this%absolute_to_relative(2,i)) = n_tgt(this%absolute_to_relative(2,i)) + 1
            endif
         enddo !i
         write(level2,'("Number of target orbitals: ")')
         write(level2,'(8(i0,1x))') n_tgt(1:this%no_irr)
         n_tgt_total = sum(n_tgt)
         write(level2,'("Total number of target orbitals: ")')
         write(level2,'(i8)') n_tgt_total
         write(level2,'("Number of continuum orbitals: ")')
         write(level2,'(8(i0,1x))') n_cnt(1:this%no_irr)

         call mo_integrals % get_array_dim(d1,d2,no_blocks)
         write(level2,'(/,"Number of columns in storage: ",i8)') d2
         write(level2,'(/,"Number of components in each column(no_mo**2): ",i8)') d1
         write(level2,'(/,"Number of blocks: ",i8)') no_blocks

         max_property_l = mo_integral_options%max_property_l
         write(level2,'(/,"Maximal l of property integrals: ",i8)') max_property_l
         write(level2,'(/"R-matrix boundary: ",f20.10)') rmatr
         write(level2,*)

         call this%ao_basis%get_continuum_l_range(l_min,l_max)
         write(level2,'(/,"Minimal ang. mom. l_min: ",i8)') l_min
         write(level2,'(/,"Maximal ang. mom. l_max: ",i8)') l_max

         max_poisson_l = min(l_max,max_property_l)
         write(level2,'(/,"Maximal ang. mom. used in the Poisson method: ",i8)') max_poisson_l

         do sym=1,nsym
            call kinetic_energy_eigenstates(this,mo_integral_storage,rmatr,sym,chan_lm_all,bamps_all,tot_chan,diag_ham_all(sym))
            write(level3,'(/,"Sym index: ",i8)') sym
            write(level3,'(/,"-- nchan: ",i8)') diag_ham_all(sym)%nchan
            write(level3,'(/,"-- n: ",i8)') diag_ham_all(sym)%n
            write(level3,'(/,"-- l, m:")')
            do i=1,diag_ham_all(sym)%nchan
               write(level3,'(2i8)') diag_ham_all(sym)%chan_lm(2,i), diag_ham_all(sym)%chan_lm(1,i)
            end do
            write(level3,'()')
         end do


         allocate(cf(no_ao,no_mo),cft(no_mo,no_ao))
         call this%get_orbital_coefficient_matrix(cf)
         write(level2,'(/,"Coefficients for transformation between AOs and MOs loaded.")')
         cft = transpose(cf)

         write(level1,'(/,"Calculation of density coefficients started.")')
         t = omp_get_wtime()

         allocate( ov3(no_mo,no_mo,no_mo), dnc(no_mo,no_mo,no_mo) )
         col_ov3 = mo_integrals%find_column_matching_name(bbb_ints)

         !Transform MAM integrals to MMM
         !$omp parallel do collapse(2) private(p, q, i, j, ind)
         do p=1,no_mo
            do q=1,no_mo
               ov3(:,q,p) = 0._cfp
               ind(1:1) = this%integral_index(bbb_ints, reshape([p, q], [2, 1]), mo_integral_options%two_p_continuum)
               do i=1,no_ao
                  if (mo_integrals%a(ind(1),col_ov3+i-1) /= 0._cfp) then
                     ov3(:,q,p) = ov3(:,q,p) + cft(:,i) * mo_integrals%a(ind(1),col_ov3+i-1)
                  end if
               end do
            end do
         end do

         write (level1, '(/,a,f0.1,a)') 'Transformation of triple overlaps took ', omp_get_wtime() - t, ' s'
         t = omp_get_wtime()

         !Transform one index to V_k (Hamiltonain eig-vecs) to obtain dnc (density coefficients) for all pairs of orbitals
         do k=1,no_mo
            add = 0
            do isym = 1, size(diag_ham_all)
               if (k <= add + diag_ham_all(isym)%n) exit
               add = add + diag_ham_all(isym)%n
            end do

            !$omp parallel do collapse(2) private(p,q,j)
            do p=1,no_mo
               do q=1,no_mo
                  dnc(k,q,p) = 0._cfp
                  do j=1,diag_ham_all(isym)%n
                     dnc(k,q,p) = dnc(k,q,p) + diag_ham_all(isym)%ham(j,k-add) * ov3(add+j,q,p)
                  end do
               end do
            end do
         end do

         write (level1, '(/,a,f0.1,a)') 'Calculation of density coefficients took ', omp_get_wtime() - t, ' s'
         t = omp_get_wtime()

         !Calculate coefficients qnc of the multipole expansion
         allocate( anc(no_mo,no_mo,no_mo), qnc((max_poisson_l+1)**2,no_mo,no_mo) )
         col_prp = mo_integrals%find_column_matching_name(property_ints)

         !$omp parallel do collapse(2) private(p,q,ind,lam,mu,j)
         do p=1,no_mo
            do q=1,no_mo
               ind(1:1) = this%integral_index(property_ints, reshape([p, q], [2, 1]), mo_integral_options%two_p_continuum)
               do lam = 0,max_poisson_l
                  do mu=-lam,lam
                     j = (lam+1)*lam + mu + 1
                     qnc(j,q,p)  = qnc(j,q,p) + mo_integrals%a(ind(1),col_prp+j-1) * sqrt(4*pi/(2*lam+1))
                  end do
               end do
            end do
         end do

         write (level1, '(/,a,f0.1,a)') 'Calculation of multipole expansion coefficients took ', omp_get_wtime() - t, ' s'
         t = omp_get_wtime()

         !Compute A_k coefficients (anc) for all pairs of orbitals
         !todo this can be moved to inside the loop 'eval_loop' over p,q indices below
         do k=1,no_mo
            add = 0
            do isym = 1, size(diag_ham_all)
               if (k <= add + diag_ham_all(isym)%n) exit
               add = add + diag_ham_all(isym)%n
            end do

            !$omp parallel do collapse(2) private(p,q,dvout,lam,mu,j,i)
            do q=1,no_mo
               do p=1,no_mo
                  dvout = 0._cfp
                  do lam = 0,max_poisson_l
                     do mu=-lam,lam
                        j = (lam+1)*lam + mu + 1
                        if (chan_lm_all(3,j) == isym) then
                           do i=1,diag_ham_all(isym)%nchan
                              if ((diag_ham_all(isym)%chan_lm(1,i) == mu) .and. (diag_ham_all(isym)%chan_lm(2,i) == lam)) then
                                 dvout = dvout + lam*diag_ham_all(isym)%w(i,k-add)*qnc(j,q,p)/(rmatr**lam)
                              end if
                           end do
                        end if
                     end do
                  end do
                  anc(k,q,p) = anc(k,q,p) + (-4*pi*dnc(k,q,p) + dvout)/diag_ham_all(isym)%e(k-add)
               end do
            end do
         end do

         write (level1, '(/,a,f0.1,a)') 'Calculation of scalar potential expansion coefficients took ', omp_get_wtime() - t, ' s'

         call this%finalize_two_electron_integrals_poisson(no_mo, mo_integral_options, mo_integral_storage, mo_integrals, dnc, anc)

         write (level1, '(/,"Two electron integrals calculated from Poisson equation.")')
         write (level3, '("<---------","molecular_orbital_basis_obj:two_electron_integrals_poisson")')

   end subroutine two_electron_integrals_poisson


   !> \brief   Fill two-electron integral arrays and write them to disk
   !> \authors M Konvalinka, Z Masin, J Benda
   !> \date    2024
   !>
   !> Uses the triple orbital overlaps and the inner region Poisson solutions expanded in molecular orbital basis to evaluate
   !> two-particle Coulomb integrals.
   !>
   subroutine finalize_two_electron_integrals_poisson(this, no_mo, mo_integral_options, mo_integral_storage, mo_integrals, dnc, anc)

      use const_gbl, only: level1, two_p_sym_ints, ijkl_indices_header, two_el_ints
      use omp_lib, only: omp_get_wtime
      use sort_gbl, only: heap_sort_int_float
      use special_functions_gbl, only: ipair
      use utils_gbl, only: reserve_space, append_array

      class(molecular_orbital_basis_obj), intent(inout) :: this
      integer,                            intent(in)    :: no_mo
      type(integral_storage_obj),         intent(inout) :: mo_integral_storage
      type(integral_options_obj),         intent(in)    :: mo_integral_options
      type(p2d_array_obj), pointer,       intent(in)    :: mo_integrals
      real(kind=cfp), allocatable,        intent(in)    :: dnc(:,:,:), anc(:,:,:)

      real(kind=cfp), allocatable, target :: integral_values(:), thread_integral_values(:)
      integer,        allocatable, target :: integral_indices(:), thread_integral_indices(:)

      real(kind=cfp), pointer :: backup_a(:,:)
      character(len=line_len) :: record_name, column_names(1)

      real(kind=cfp) :: aaaa
      real(kind=wp)  :: t
      integer        :: i, i_cont, ij, j, j_cont, n, p, p_cont, pq, q, q_cont, record_begin, record_end, header_end, lunit, err

         write(level1,'(/,"Calculating and storing the two-electron integrals:")')

         t = omp_get_wtime()

         n = 0
         !$omp parallel default(none) reduction(+:n) &
         !$omp& shared(this,no_mo,mo_integral_options,dnc,anc,integral_values,integral_indices) &
         !$omp& private(p,q,i,j,p_cont,q_cont,i_cont,j_cont,pq,ij,aaaa,thread_integral_values,thread_integral_indices)
         allocate (thread_integral_values(100), thread_integral_indices(100))
         !$omp do schedule(dynamic,1)
         do p=1,no_mo
            p_cont = merge(1, 0, this%is_continuum(p))
            do q = 1, p
               q_cont = merge(1, 0, this%is_continuum(q))
               pq = ipair(p) + q
               do i = 1, p
                  i_cont = merge(1, 0, this%is_continuum(i))
                  if (.not. mo_integral_options%two_p_continuum .and. i_cont + p_cont + q_cont == 3) cycle
                  do j = 1, i
                     j_cont = merge(1, 0, this%is_continuum(j))
                     if (.not. mo_integral_options%two_p_continuum .and. i_cont + j_cont + p_cont + q_cont >= 3) cycle
                     ij = ipair(i) + j
                     if (pq >= ij) then
                        aaaa = sum(anc(1:no_mo,q,p)*dnc(1:no_mo,j,i))
                        if (abs(aaaa) > mo_integral_options % tol)  then
                           n = n + 1
                           call reserve_space(thread_integral_values, n)
                           call reserve_space(thread_integral_indices, n)
                           thread_integral_values(n) = aaaa
                           thread_integral_indices(n) = ipair(pq) + ij
                        end if
                     end if
                  end do
               end do
            end do
         end do
         !$omp critical(append_thread_arrays)
         call append_array(integral_values, thread_integral_values(1:n))
         call append_array(integral_indices, thread_integral_indices(1:n))
         !$omp end critical(append_thread_arrays)
         !$omp end parallel

         call heap_sort_int_float(n, integral_indices, integral_values)

         write (level1, '(/,a,f0.1,a)') 'Evaluation of the two-electron integrals took ', omp_get_wtime() - t, ' s'

         !If requested, print the non-zero integrals
         if (mo_integral_options % print_integrals .and. myrank == master) then
            call mo_integrals % print(.true.)
         end if

         !Dump all integrals to disk and close the record.
         if (mo_integral_storage % on_disk()) then

            write (level2, '("Saving integrals to disk...")')

            record_name  = mo_integral_storage % contruct_header_string(this % get_basis_name(), two_p_sym_ints)
            record_begin = mo_integral_storage % integral_file % start_record(record_name)
            lunit        = mo_integral_storage % integral_file % get_unit_no()

            column_names(1) = two_el_ints
            err = mo_integrals % init(n, 1, 0, column_names)
            if (err /= 0) then
               call mpi_xermsg('molecular_basis_mod', 'finalize_two_electron_integrals_poisson', &
                               'Molecular integral storage initialization failed.', err, 1)
            end if

            backup_a => mo_integrals % a
            mo_integrals % a(1:n, 1:1) => integral_values(1:n)
            this % ind_ijkl_integral = n
            this % ijkl_indices(1:n, 1:1) => integral_indices(1:n)

            call mo_integral_options % write(lunit, record_begin, header_end)
            call mo_integrals % write(lunit, header_end, record_end, int(master))
            call mo_integral_storage % integral_file % close_record(record_name, record_begin, record_end)

            mo_integrals % a => backup_a
            err = mo_integrals % final()
            if (err /= 0) then
                call mpi_xermsg('molecular_orbital_basis_obj', 'finalize_two_electron_integrals_poisson', &
                                'Deallocation of the temporary integral array failed.', err, 1)
            end if

            record_name  = mo_integral_storage % contruct_header_string(this % get_basis_name(), ijkl_indices_header)
            record_begin = mo_integral_storage % integral_file % start_record(record_name)

            call this % write_ijkl_indices(lunit, record_begin, record_end)
            call mo_integral_storage % integral_file % close_record(record_name, record_begin, record_end)

            this % ijkl_indices => null()

            write (level2, '("Saving integrals to disk... done")')

         end if

   end subroutine finalize_two_electron_integrals_poisson


   !> \brief   Diagonalize Hamiltonian blocks
   !> \authors M Konvalinka, Z Masin, J Benda
   !> \date    2024
   !>
   !> Obtain eigenstates of the kinetic energy operator separately for each irreducible representation and store the
   !> results in `diag_ham`. Also evaluate inner region boundary amplitudes of the eigenstates in individual channels.
   !>
   subroutine kinetic_energy_eigenstates(this,mo_integral_storage,a,sym,chan_lm_all,bamps_all,tot_chan,diag_ham)
      use phys_const_gbl, only: to_ev
      use blas_lapack_gbl, only: blasint, syev
      use const_gbl, only: level1, level2, level3, kinetic_ints, one_p_sym_ints
      use parallel_arrays_gbl
      implicit none
      class(integral_storage_obj), intent(inout) :: mo_integral_storage
      class(molecular_orbital_basis_obj), intent(inout) :: this
      class(diag_ham_obj), intent(out) :: diag_ham
      integer, intent(in) :: sym, tot_chan, chan_lm_all(3,tot_chan)
      real(kind=cfp), intent(in) :: bamps_all(tot_chan,this%number_of_functions)
      real(kind=cfp), intent(in) :: a

      integer :: err, i, j, k, channel, ind(1:1), sym_start, sym_end
      type(p2d_array_obj), pointer :: mo_integrals
      type(p2d_array_obj), target :: ao_integrals
      type(integral_storage_obj), target :: mo_integral_disk
      type(data_header_obj) :: header
      type(integral_options_obj) :: options
      real(kind=cfp), allocatable :: ham(:,:)
      integer, allocatable :: l_list(:), chan_lm(:,:)
      logical, parameter :: tgt_is_local = .true., two_p_continuum = .false.

      real(kind=cfp), allocatable :: work(:), e(:), w(:,:)
      real(kind=cfp) :: energy, p
      integer(blasint) :: lwork, info, n, nchan
      integer :: col_kei

         write(level3,'(/,"--------->molecular_basis_mod:kinetic_energy_eigenstates")')

         call mpi_mod_barrier(err)

         if (mo_integral_storage%in_memory()) then
            mo_integrals => mo_integral_storage%integral
         else !on disk
            err = mo_integral_disk%init(memory=ao_integrals)
            if (err .ne. 0) call mpi_xermsg('molecular_basis_mod','kinetic_energy_eiegenstates','Memory allocation 1 failed',err,1)
            !we need to get the full header corresponding to the transformed 1-electron integrals.
            header%name = mo_integral_storage%contruct_header_string(this%get_basis_name(),one_p_sym_ints)

            !read the integrals from the disk into mo_integrals
            call mo_integral_disk%read(mo_integral_storage,header%name,options,tgt_is_local)

            mo_integrals => ao_integrals
         end if

         !Get the column number for the mo_integrals%a array where the KEI are stored.
         col_kei = mo_integrals%find_column_matching_name(kinetic_ints)
         !FROM THIS POINT ONWARDS mo_integrals POINTS TO THE STORAGE IN MEMORY OF THE KINETIC ENERGY INTEGRALS

         !BUILD THE HAMILTONIAN = KINETIC ENERGY MATRIX

         write(level2,'(/,"Symmetry: ",i1)') sym

         n = this%get_number_of_orbitals(sym)

         if (n .eq. 0) then
            write(level2,'("No orbitals in this symmetry: skipping")')
            diag_ham%n=n
            diag_ham%nchan=0
            return
         endif

         if (allocated(ham)) deallocate(ham)
         allocate(ham(n,n),stat=err)
         if (err .ne. 0) call mpi_xermsg('molecular_basis_mod','kinetic_energy_eiegenstates','Memory allocation 2 failed',err,1)

         ham = 0.0_cfp

         sym_start = 1
         do i=1,sym-1
            sym_start = sym_start+this%get_number_of_orbitals(i)
         enddo
         sym_end = sym_start-1 + n

         do i=sym_start,sym_end
            do j=sym_start,i
               ind(1:1) = this%integral_index(kinetic_ints,reshape(source=(/j,i/),shape=(/2,1/)),two_p_continuum)
               ham(j-sym_start+1,i-sym_start+1) = mo_integrals%a(ind(1),col_kei)
               ham(i-sym_start+1,j-sym_start+1) = ham(j-sym_start+1,i-sym_start+1)
            enddo
         enddo

         !DIAGONALIZE THE HAMILTONIAN

         if (allocated(work)) deallocate(work)
         if (allocated(e)) deallocate(e)
         allocate(work(1:1),e(1:n),stat=err)
         if (err .ne. 0) call mpi_xermsg('molecular_basis_mod','kinetic_energy_eiegenstates','Memory allocation 3 failed',err,1)
         !determine the size of the auxiliary array WORK needed for diagonalization of H+L
         lwork = -1
         call syev('V','U',n,ham,n,e,work,lwork,info)
         if (info /= 0) then
            call mpi_xermsg ('molecular_basis_mod','kinetic_energy_eiegenstates', &
               'SYEV call 1 failed with an error or you are linking against 32bit integer MKL &
               &while using 64bit default integer.', int(info), 1)
         end if

         lwork=work(1)
         deallocate(work); allocate(work(1:lwork),stat=err)
         if (err .ne. 0) call mpi_xermsg('molecular_basis_mod','kinetic_energy_eiegenstates','Memory allocation 4 failed',err,1)

         !eigenvectors of H+L are normalized; eigenvectors are held in columns of ham
         call syev('V','U',n,ham,n,e,work,lwork,info)
         if (info /= 0) then
            call mpi_xermsg ('molecular_basis_mod','kinetic_energy_eiegenstates', &
               'SYEV call 2 failed with an error or you are linking against 32bit integer MKL &
               &while using 64bit default integer.', int(info), 1)
         end if

         write(level3,'(/,"Eigenvalues of the Kinetic energy matrix:")')
         do i=1,n
            write(level3,'(i0,e25.15)') i,e(i)
            if (e(i) < 0.0_cfp) then
               call mpi_xermsg ('molecular_basis_mod','kinetic_energy_eiegenstates', &
                  'At least one eigenvalue of the kinetic energy matrix is negative: &
                  &severe numerical instability or R-matrix radius is too small.', 2, 1)
            end if
         enddo

         !GET CHANNEL DATA AND CONSTRUCT THE TRUE BOUNDARY AMPLITUDES

         if (allocated(w)) deallocate(w)
         if (allocated(l_list)) deallocate(l_list)

         write(level2,'(/,"Partial waves (m,l) in the current symmetry:")')
         nchan = 0
         do i=1,tot_chan
            if (chan_lm_all(3,i) .eq. sym) then
               nchan = nchan + 1
               write(level2,'(2i5)') chan_lm_all(1:2,i)
            endif
         enddo

         write(level2,'(/,"Number of channels in the current symmetry: ",i0)') nchan

         allocate(w(nchan,n),l_list(nchan),stat=err)
         allocate(chan_lm(2,nchan))
         if (err .ne. 0) call mpi_xermsg('molecular_basis_mod','kinetic_energy_eiegenstates','Memory allocation 5 failed',err,1)

         w(:,:) = 0.0_cfp
         channel = 0
         do i=1,tot_chan
            if (chan_lm_all(3,i) .eq. sym) then
               channel = channel + 1
               l_list(channel) = chan_lm_all(2,i) !save the L-values for each channel
               chan_lm(1,channel) = chan_lm_all(1,i)  ! save the L, M values for diag_ham
               chan_lm(2,channel) = chan_lm_all(2,i)
               !calculate the R-matrix boundary amplitudes in each channel
               do k=1,n
                  do j=1,n !over all coefficients of the R-matrix basis functions
                     w(channel,k) = w(channel,k) + ham(j,k)*bamps_all(i,sym_start-1+j)
                     !write(stdout,'("cf",i10,2e25.15)') j,ham(j,k),bamps_all(i,sym_start-1+j)
                  enddo
                  !write(stdout,'("amplitude",2i10,e25.15)') channel,k,w(channel,k)
               enddo
            endif
         enddo

         !setting the results into the output object diag_ham
         if (allocated(diag_ham%ham)) deallocate(diag_ham%ham)
         if (allocated(diag_ham%w)) deallocate(diag_ham%w)
         if (allocated(diag_ham%e)) deallocate(diag_ham%e)
         if (allocated(diag_ham%chan_lm)) deallocate(diag_ham%chan_lm)

         allocate(diag_ham%ham(n,n))
         allocate(diag_ham%w(nchan,n))
         allocate(diag_ham%e(n))
         allocate(diag_ham%chan_lm(2,nchan))

         diag_ham%nchan=nchan
         diag_ham%n=n
         diag_ham%ham=ham
         diag_ham%w=w
         diag_ham%e=-2*e
         diag_ham%chan_lm = chan_lm

         call mpi_mod_barrier(err)

         write(level3,'("<---------","molecular_orbital_basis_obj:kinetic_energy_eigenstates")')

         err = ao_integrals % final()

   end subroutine kinetic_energy_eigenstates

end module molecular_basis_gbl
