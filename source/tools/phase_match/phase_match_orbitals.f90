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
module phase_match_orbitals_mod
   use precisn_gbl, only: wp, cfp

contains

   !Merge molecular data for two geometries into one and calculate the overlap integrals between the two sets of orbitals.
   subroutine load_orbital_overlaps(integral_options,orbital_olaps,moints_ref,moints_match,ref_indices,match_indices)
      use mpi_gbl, only: mpi_xermsg
      use utils_gbl, only: xermsg
      use const_gbl, only: overlap_ints, level1, level2, C1_id, sym_op_nam_len
      use basis_data_generic_gbl, only: CGTO_shell_data_obj, BTO_shell_data_obj, orbital_data_obj
      use atomic_basis_gbl, only: atomic_orbital_basis_obj
      use molecular_basis_gbl, only: molecular_orbital_basis_obj
      use parallel_arrays_gbl, only: p2d_array_obj
      use integral_storage_gbl, only: integral_storage_obj, integral_options_obj
      use symmetry_gbl, only: geometry_obj
      implicit none
      real(kind=cfp), allocatable :: orbital_olaps(:,:)
      type(integral_options_obj), intent(in) :: integral_options
      character(len=132), intent(in) :: moints_ref, moints_match
      integer, allocatable :: ref_indices(:,:), match_indices(:,:)

      logical, parameter :: force_C1 = .true.  ! for debugging

      integer :: i, j, k, l, ii, t, n_cfs, n_orbs, nnuc_ref, nnuc_match, nnuc_merged, err, cnt, ind(1:1), two_ind(2,1), col, last
      integer :: n_cgto_shells_ref, n_cgto_shells_match, n_cgto_shells_merged, nirr, nirr_ref, nirr_match, pg_ref, pg_match
      integer :: n_bto_shells_ref, n_bto_shells_match, n_bto_shells_merged, pg, nref, nmtc
      integer :: offset_ref_continuum, offset_match_continuum, n_ref_continuum, n_match_continuum, n_ref_target, n_match_target
      logical :: same, have_btos, have_continuum
      type(atomic_orbital_basis_obj), target :: atomic_orbital_basis_ref, atomic_orbital_basis_match, atomic_orbital_basis_merged
      type(molecular_orbital_basis_obj) :: molecular_orbital_basis_ref, molecular_orbital_basis_match,molecular_orbital_basis_merged
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data_ref(:), CGTO_shell_data_match(:)
      type(BTO_shell_data_obj), allocatable :: BTO_shell_data_ref(:), BTO_shell_data_match(:)
      type(geometry_obj) :: geometry_ref, geometry_match, geometry_merged
      type(orbital_data_obj), allocatable :: orbital_data_ref(:), orbital_data_match(:)
      type(orbital_data_obj) :: orbital_data
      type(integral_storage_obj), target :: atomic_integral_storage, molecular_integral_storage
      type(p2d_array_obj), target :: ao_integrals, mo_integrals
      real(kind=cfp), allocatable :: charge_densities(:,:)
      integer, allocatable :: ref_irr_map(:), match_irr_map(:)
      logical, allocatable :: to_delete(:)

         write(*,'(/,5X,"Constructing spin-orbital overlap matrix")')
         write(*,'(  5X,"----------------------------------------")')
         write(*,'(/,5X,"Reference basis sets on file: ",a)') trim(adjustl(moints_ref))
         write(*,'(  5X,"Matching basis sets on file: ",a)') trim(adjustl(moints_match))

         !REFERENCE SET: read-in all basis sets and orbital data
         call atomic_orbital_basis_ref%read(moints_ref)
         molecular_orbital_basis_ref%ao_basis => atomic_orbital_basis_ref
         call molecular_orbital_basis_ref%read(moints_ref)

         call atomic_orbital_basis_ref%get_all_CGTO_shells(CGTO_shell_data_ref,n_cgto_shells_ref)
         if (atomic_orbital_basis_ref%contains_btos()) then
            call atomic_orbital_basis_ref%get_all_BTO_shells(BTO_shell_data_ref,n_bto_shells_ref)
         else
            n_bto_shells_ref = 0
         end if
         call molecular_orbital_basis_ref%symmetry_data%get_geometry(geometry_ref)

         write(*,'(5X,"Reference basis set and orbitals have been read-in")')

         !MATCH SET: read-in all basis sets and orbital data
         call atomic_orbital_basis_match%read(moints_match)
         molecular_orbital_basis_match%ao_basis => atomic_orbital_basis_match
         call molecular_orbital_basis_match%read(moints_match)

         call atomic_orbital_basis_match%get_all_CGTO_shells(CGTO_shell_data_match,n_cgto_shells_match)
         if (atomic_orbital_basis_match%contains_btos()) then
            call atomic_orbital_basis_match%get_all_BTO_shells(BTO_shell_data_match,n_bto_shells_match)
         else
            n_bto_shells_match = 0
         end if
         call molecular_orbital_basis_match%symmetry_data%get_geometry(geometry_match)

         write(*,'(5X,"Matching basis set and orbitals have been read-in")')
         write(*,'(5X,"Number of CGTO shells from the reference basis: ",i10)') n_cgto_shells_ref
         write(*,'(5X,"Number of BTO shells from the reference basis:  ",i10)') n_bto_shells_ref
         write(*,'(5X,"Number of CGTO shells from the basis to match:  ",i10)') n_cgto_shells_match
         write(*,'(5X,"Number of BTO shells from the basis to match:   ",i10)') n_bto_shells_match

         nnuc_ref = geometry_ref%no_nuc
         nnuc_match = geometry_match%no_nuc

         write(*,'(5X,"Number of nuclei from the reference basis:      ",i10)') nnuc_ref
         write(*,'(5X,"Number of nuclei from the basis to match:       ",i10)') nnuc_match

         pg_ref = molecular_orbital_basis_ref % pg
         pg_match = molecular_orbital_basis_match % pg
         nirr_ref = molecular_orbital_basis_ref % no_irr
         nirr_match = molecular_orbital_basis_match % no_irr

         ! ref_irr_map and match_irr_map provide mapping from the original irreducible representations to the irreducible
         ! representations of the subgroup used by the merged basis
         allocate (ref_irr_map(nirr_ref), &
                   match_irr_map(nirr_match), &
                   orbital_data_ref(nirr_ref), &
                   orbital_data_match(nirr_match))

         ! decide on the point group of the merged geometries (original if both the same, C1 otherwise)
         if (pg_ref == pg_match .and. .not. force_C1) then
            print '(/,5X,A)', 'Using the original point groups.'
            pg = pg_ref
            nirr = nirr_ref
            geometry_merged % no_sym_op = geometry_ref % no_sym_op
            geometry_merged % sym_op = geometry_ref % sym_op
            do i = 1, nirr
               ref_irr_map(i) = i
               match_irr_map(i) = i
            end do
         else
            print '(/,5X,A)', 'Downgrading to the C1 point group.'
            pg = C1_id
            nirr = 1
            geometry_merged % no_sym_op = 0
            ref_irr_map = 1
            match_irr_map = 1
         end if

         !JOIN the geometry data, CGTO basis sets and orbital sets.
         nnuc_merged = nnuc_ref+nnuc_match
         allocate(geometry_merged%nucleus(nnuc_merged))
         geometry_merged%no_nuc = nnuc_merged
         geometry_merged%use_symmetry = .false.
  
         have_continuum = .false.
         j = 0
         do i=1,nnuc_ref
            if (.not.(geometry_ref%nucleus(i)%is_continuum())) then
               j = j + 1
               geometry_merged%nucleus(j) = geometry_ref%nucleus(i)
               geometry_merged%nucleus(j)%nuc = j
               call geometry_merged%nucleus(j)%print
            else
               have_continuum = .true.
            endif
         enddo
         do i=1,nnuc_match
            if (.not.(geometry_match%nucleus(i)%is_continuum())) then
               j = j + 1
               geometry_merged%nucleus(j) = geometry_match%nucleus(i)
               geometry_merged%nucleus(j)%nuc = j
               call geometry_merged%nucleus(j)%print
            else
               have_continuum = .true.
            endif
         enddo

         if (have_continuum) then
            !Resize the geometry structure to the appropriate number of nuclei:
            deallocate(geometry_ref%nucleus)
            call move_alloc(geometry_merged%nucleus,geometry_ref%nucleus)
            allocate(geometry_merged%nucleus(j))
            geometry_merged%nucleus(1:j) = geometry_ref%nucleus(1:j)
            geometry_merged%no_nuc = j
 
            call geometry_merged%add_scattering_centre
         endif

         write(*,'(/,5X,"Nuclear data have been merged.")')

         n_cgto_shells_merged = n_cgto_shells_ref+n_cgto_shells_match
         n_bto_shells_merged = n_bto_shells_ref+n_bto_shells_match
         n_cfs = atomic_orbital_basis_ref%number_of_functions + atomic_orbital_basis_match%number_of_functions
         allocate(to_delete(n_cfs))

         !INITIALIZE the joined basis sets
         err = atomic_orbital_basis_merged%init(n_cgto_shells_merged+n_bto_shells_merged,geometry_merged)

         !Target functions: must be first in the basis ("reference" followed by "match")
         offset_ref_continuum = 0
         n_ref_target = 0
         n_match_target = 0
         k = 0
         to_delete = .false.
         do i=1,n_cgto_shells_ref
            if (.not.(CGTO_shell_data_ref(i)%is_continuum())) then
               offset_ref_continuum = offset_ref_continuum + CGTO_shell_data_ref(i)%number_of_functions
               n_ref_target = n_ref_target + CGTO_shell_data_ref(i)%number_of_functions
               if (CGTO_shell_data_ref(i)%l .eq. 0) to_delete(k+1:k+CGTO_shell_data_ref(i)%number_of_functions) = .true.
               k = k + CGTO_shell_data_ref(i)%number_of_functions
               call atomic_orbital_basis_merged%add_shell(CGTO_shell_data_ref(i))
            endif
         enddo
         do i=1,n_cgto_shells_match
            if (.not.(CGTO_shell_data_match(i)%is_continuum())) then
               offset_ref_continuum = offset_ref_continuum + CGTO_shell_data_match(i)%number_of_functions
               n_match_target = n_match_target + CGTO_shell_data_match(i)%number_of_functions
               if (CGTO_shell_data_match(i)%l .eq. 0) to_delete(k+1:k+CGTO_shell_data_match(i)%number_of_functions) = .true.
               k = k + CGTO_shell_data_match(i)%number_of_functions
               call atomic_orbital_basis_merged%add_shell(CGTO_shell_data_match(i))
            endif
         enddo

         !Continuum functions: must be last in the basis ("reference" followed by "match")
         offset_match_continuum = offset_ref_continuum
         n_ref_continuum = 0
         do i=1,n_cgto_shells_ref
            if (CGTO_shell_data_ref(i)%is_continuum()) then
               n_ref_continuum = n_ref_continuum + CGTO_shell_data_ref(i)%number_of_functions
               offset_match_continuum = offset_match_continuum + CGTO_shell_data_ref(i)%number_of_functions
               call atomic_orbital_basis_merged%add_shell(CGTO_shell_data_ref(i))
            endif
         enddo
         do i=1,n_bto_shells_ref
            if (BTO_shell_data_ref(i)%is_continuum()) then
               n_ref_continuum = n_ref_continuum + BTO_shell_data_ref(i)%number_of_functions
               offset_match_continuum = offset_match_continuum + BTO_shell_data_ref(i)%number_of_functions
               call atomic_orbital_basis_merged%add_shell(BTO_shell_data_ref(i))
            endif
         enddo
         n_match_continuum = 0
         do i=1,n_cgto_shells_match
            if (CGTO_shell_data_match(i)%is_continuum()) then
               n_match_continuum = n_match_continuum + CGTO_shell_data_match(i)%number_of_functions
               call atomic_orbital_basis_merged%add_shell(CGTO_shell_data_match(i))
            endif
         enddo
         do i=1,n_bto_shells_match
            if (BTO_shell_data_match(i)%is_continuum()) then
               n_match_continuum = n_match_continuum + BTO_shell_data_match(i)%number_of_functions
               call atomic_orbital_basis_merged%add_shell(BTO_shell_data_match(i))
            endif
         enddo

         write(*,'(5X,"Atomic basis sets have been merged.")')

         molecular_orbital_basis_merged%ao_basis => atomic_orbital_basis_merged
         err = molecular_orbital_basis_merged%init(nirr, geometry_merged)

         orbital_data%number_of_coefficients = n_cfs

         !ref_indices, match_indices map the indices of orbitals within their sets to the merged orbital set.
         if (allocated(ref_indices)) deallocate(ref_indices)
         if (allocated(match_indices)) deallocate(match_indices)
         allocate(ref_indices(molecular_orbital_basis_ref % number_of_functions, nirr_ref), &
                  match_indices(molecular_orbital_basis_match % number_of_functions, nirr_match))
         ref_indices = -1
         match_indices = -1
         last = 0

         do i = 1, nirr
            print '(/,5X,A,I0,A)', 'Building IRR ', i, ' of the merged molecular basis'

            if (allocated(orbital_data%coefficients)) deallocate(orbital_data%coefficients)
            if (allocated(orbital_data%energy)) deallocate(orbital_data%energy)
            if (allocated(orbital_data%spin)) deallocate(orbital_data%spin)
            if (allocated(orbital_data%occup)) deallocate(orbital_data%occup)

            ! initialize orbital list for this irreducible representation of the merged basis
            orbital_data%number_of_functions = 0
            orbital_data%point_group         = pg
            orbital_data%irr                 = i

            ! retrieve all (reference and match) molecular orbitals for this merged irreducible representation
            do j = 1, nirr_ref
               if (ref_irr_map(j) == i) then
                  call molecular_orbital_basis_ref % get_shell_data(j, orbital_data_ref(j))
                  orbital_data % number_of_functions = orbital_data % number_of_functions &
                                                     + orbital_data_ref(j) % number_of_functions
               end if
            end do
            do j = 1, nirr_match
               if (match_irr_map(j) == i) then
                  call molecular_orbital_basis_match % get_shell_data(j, orbital_data_match(j))
                  orbital_data % number_of_functions = orbital_data % number_of_functions &
                                                     + orbital_data_match(j) % number_of_functions
               end if
            end do

            allocate(orbital_data % coefficients(n_cfs, orbital_data%number_of_functions), &
                     orbital_data % energy(orbital_data%number_of_functions), &
                     orbital_data % spin(orbital_data%number_of_functions), &
                     orbital_data % occup(orbital_data%number_of_functions))

            orbital_data%coefficients = 0.0_cfp
            orbital_data%energy = 1.0_cfp
            orbital_data%spin = 0

            ! number of orbitals added from the reference and match basis, respectively
            nref = 0
            nmtc = 0

            ! collect expansion coefficients of reference molecular orbitals in this irreducible representation
            write (level1, '(/,10X,A5,1X,3A5)') 'basis', 'irr ', 'orgI', 'mrgI'
            do j = 1, nirr_ref
               if (ref_irr_map(j) == i) then
                  l = orbital_data_ref(j) % number_of_coefficients
                  do k = 1, orbital_data_ref(j) % number_of_functions
                     nref = nref + 1
                     if (n_ref_target > 0) then
                        orbital_data % coefficients(1:n_ref_target, nref) &
                           = orbital_data_ref(j) % coefficients(1:n_ref_target, k)
                     end if
                     if (n_ref_continuum > 0) then
                        orbital_data % coefficients(offset_ref_continuum+1:offset_ref_continuum+n_ref_continuum, nref) &
                           = orbital_data_ref(j) % coefficients(n_ref_target+1:l, k)
                     end if
                     last = last + 1
                     ref_indices(k, j) = last
                     write (level1, '(10X,A5,1X,3I5)') 'ref', j, k, last
                  end do
               end if
            end do

            ! collect expansion coefficients of match molecular orbitals in this irreducible representation
            do j = 1, nirr_match
               if (match_irr_map(j) == i) then
                  l = orbital_data_match(j) % number_of_coefficients
                  do k = 1, orbital_data_match(j) % number_of_functions
                     nmtc = nmtc + 1
                     if (n_match_target > 0) then
                        orbital_data % coefficients(n_ref_target+1:n_ref_target+n_match_target, nref+nmtc) &
                           = orbital_data_match(j) % coefficients(1:n_match_target, k)
                     end if
                     if (n_match_continuum > 0) then
                        orbital_data % coefficients(offset_match_continuum+1:offset_match_continuum+n_match_continuum, nref+nmtc) &
                        = orbital_data_match(j) % coefficients(n_match_target+1:l, k)
                     end if
                     last = last + 1
                     match_indices(k, j) = last
                     write (level1, '(10X,A5,1X,3I5)') 'match', j, k, last
                  end do
               end if
            end do

            call molecular_orbital_basis_merged%add_shell(orbital_data)
         enddo

         write(*,'(/,5X,"Molecular orbital data have been merged.")')
  
         call molecular_orbital_basis_merged%print_orbitals
!
! CALCULATE THE ATOMIC 1-ELECTRON INTEGRALS
!
         !describe where the AO integrals will be stored 
         err = atomic_integral_storage%init(memory=ao_integrals)
         if (err .ne. 0) then
            print *,err
            call mpi_xermsg('main','main','error initializing the target atomic_integral_storage',1,1)
         endif
   
         call atomic_orbital_basis_merged%one_electron_integrals(atomic_integral_storage,integral_options)

         write(*,'(/,5X,"Atomic overlap integrals have been calculated.")')
!
! TRANSFORM THE 1-ELECTRON ATOMIC INTEGRALS INTO INTEGRALS OVER THE MOLECULAR ORBITALS:
!
         !Describe where the transformed AO->MO integrals will be stored :
         err = molecular_integral_storage%init(memory=mo_integrals)
         if (err .ne. 0) then
            print *,err
            call mpi_xermsg('main','main','error initializing the target molecular_integral_storage',1,1)
         endif
   
         molecular_orbital_basis_merged%ao_integral_storage => atomic_integral_storage !point to the storage for the atomic integrals
         call molecular_orbital_basis_merged%one_electron_integrals(molecular_integral_storage,integral_options)

         write(*,'(5X,"Overlap integrals in the basis of molecular orbitals have been calculated.",/)')

         !Column number in one_electron_integrals%a corresponding to the overlap integrals
         col = mo_integrals%find_column_matching_name(overlap_ints)

         i = molecular_orbital_basis_merged%number_of_functions
         if (allocated(orbital_olaps)) deallocate(orbital_olaps)
         allocate(orbital_olaps(i,i))
         orbital_olaps = 0.0_cfp

         !Print the cross-overlaps between the two orbital sets
         !do i=1,molecular_orbital_basis_ref%number_of_functions
         do i=1,molecular_orbital_basis_merged%number_of_functions
            two_ind(1,1) = i
            !The j-index should correspond to the second orbital set.
            !do j=molecular_orbital_basis_ref%number_of_functions+1,molecular_orbital_basis_merged%number_of_functions
            do j=1,i
               two_ind(2,1) = j
               ind(1:1) = molecular_orbital_basis_merged%integral_index(overlap_ints,two_ind,integral_options%two_p_continuum)

               orbital_olaps(i,j) = mo_integrals%a(ind(1),col)
               orbital_olaps(j,i) = mo_integrals%a(ind(1),col)
               if (i .eq. j) then
                  write(level2,'(5X,"DIAG orbital norm",2i10,e25.15)') i,j,sqrt(orbital_olaps(i,j))
               else
                  write(level2,'(5X,"cross orbital olap",2i10,e25.15)') i,j,orbital_olaps(i,j)
               endif
            enddo
         enddo

         write(*,'(/,5X,"Orbital overlap matrix has been constructed.")')

         !call molecular_orbital_basis_merged%radial_charge_density(10.0_cfp,0.0_cfp,10.0_cfp,0.1_cfp,.true.,charge_densities)

   end subroutine load_orbital_overlaps

   !> Match orbitals between two geometries by matching those with the largest overlap.
   !> Alternatively no matching is done  (if find_matching_orbitals == .false.)
   !> and only the overlaps between specified pairs of orbitals between the two geometries
   !> are reported.
   !> The matching is done by gradually shrinking the pool of orbitals avaiable
   !> for matching by excluding those that have been matched already. This ensures
   !> that always a unique match between the two orbital sets is produced but that
   !> doesn't mean this method 100% reliable: in the future matching based on
   !> comparison of a number of multipoles may be implemented. This should be a
   !> more reliable method since it will effectively compare also the symmetrical
   !> shape of the orbitals.
   subroutine match_orbitals (orbital_olaps, orbital_phases_one_geom, ref_indices, match_indices, n_orbitals_to_match, &
                              orbital_num, orbital_sym, geom, find_matching_orbitals)
      use algorithms, only: findloc
      implicit none
      real(kind=cfp), allocatable :: orbital_olaps(:,:)
      integer, allocatable :: ref_indices(:,:), match_indices(:,:)
      integer, intent(in) :: n_orbitals_to_match, orbital_num(:), orbital_sym(:), geom
      integer, intent(out) :: orbital_phases_one_geom(n_orbitals_to_match)
      logical, intent(in) :: find_matching_orbitals

      integer :: i,j,ii,jj,offset,stride,j_match,ref_orbital,ref_symmetry,match_orbital,match_symmetry
      logical, allocatable :: mask(:)
      character(len=*), parameter :: str_match = "Matching: ", str_compare = "Comparing: "
      character(len=len(str_compare)) :: str

         offset = (geom-2)*n_orbitals_to_match
         stride = n_orbitals_to_match

         ! orbital_olaps keeps track of what orbitals in the combined orbital basis are available for matching
         allocate (mask(size(orbital_olaps, 2)))
         mask = .true.

         ! mask out all reference orbitals from the combined orbital basis (we want to match only to the "match" orbitals)
         do i = 1, size(ref_indices, 1)
            do j = 1, size(ref_indices, 2)
               if (ref_indices(i, j) > 0) then
                  mask(ref_indices(i, j)) = .false.
               end if
            end do
         end do

         do ii=1,n_orbitals_to_match
            ref_orbital = orbital_num(offset + ii)
            ref_symmetry = orbital_sym(offset + ii)
            match_orbital = orbital_num(offset + ii + stride)
            match_symmetry = orbital_sym(offset + ii + stride)

            if (ref_orbital <= 0 .or. ref_symmetry <= 0 .or. match_orbital <= 0 .or. match_symmetry <= 0) then
               print *,'error in input data for geometry', geom
               print *,'orbital_num and/or orbital_sym have not been set (properly)'
               stop 1
            endif

            i = ref_indices(ref_orbital, ref_symmetry)
            !Diagonal orbital index, i.e. as if the orbital has not changed index from geometry to geometry
            j = match_indices(match_orbital, match_symmetry)

            print '(/,5X,3(A,I0),A)', 'Reference orbital ', i, ' (', ref_orbital, '.', ref_symmetry, ')'
            print '(  5X,3(A,I0),A)', 'Diagonal orbital ', j, ' (', match_orbital, '.', match_symmetry, ')'

            if (find_matching_orbitals) then
               !Find the actual matching orbital from the current geometry
               j_match = maxloc(abs(orbital_olaps(i, :)), 1, mask)
               str = str_match
            else
               j_match = j
               str = str_compare
            endif

            orbital_phases_one_geom(ii) = int(sign(1.0_cfp,orbital_olaps(i,j_match)))
            mask(j_match) = .false.  ! the orbital "j_match" from the combined basis will not be matched again
            match_orbital = findloc(match_indices(:, match_symmetry), j_match, 1)

            write(*,'(5X,a10,i5,".",i1," -> ",i5,".",i1)') trim(str), ref_orbital, ref_symmetry, match_orbital, match_symmetry
            if (j_match /= j) then
               write(*,'(5X,"ORBITALS IN THE MATCHED GEOMETRY HAVE SWAPPED!")')
            endif
            write(*,'(5X,"Indices in the merged basis: ",2i5)') i,j
            write(*,'(5X,"Orbital overlap and sign: ",e25.15,",",i2)') orbital_olaps(i,j_match), orbital_phases_one_geom(ii)
         enddo

   end subroutine match_orbitals

end module phase_match_orbitals_mod

!> \brief   Program phase_match_orbitals
!> \author  Z Masin, J Benda
!> \date    2016 - 2022
!>
!> This program reads a given number of GBTOlib data files. It begins by combining the molecular orbitals from the first two
!> files. The molecular orbitals for each irreducible representation are simply concatenated (target orbitals first, then
!> continuum). Overlap integrals for this double molecular orbital basis are calculated, resulting in a set of overlaps
!> that combines overlaps integrals of the first basis set, overlap integrals of the second basis and cross-overlaps of the
!> two basis sets. The cross-overlaps are then used to decide (see \ref phase_match_orbitals_mod::match_orbitals) which
!> orbitals in the two sets are the most similar.
!>
!> If the source files have different point groups, the merged molecular orbital basis is considered to be C1.
!>
!> Sample input:
!>
!> \verbatim
!> &input
!>     n_geom = 2,
!>     moints = 'moints1', 'moints2',
!>
!>     n_orbitals_to_match = 3,
!>     orbital_sym = 1,1,1, 1,1,1,  ! use orbitals from the first irreducible representation
!>     orbital_num = 1,2,3, 1,2,3,  ! match first three orbitals
!>
!>     find_matching_orbitals = .false.,
!>     rmat_radius = 15,
!>
!>     ! the following are only needed for B-spline overlap integrals
!>     mixed_ints_method = 1,
!>     max_l_legendre_1el = 75,
!>     delta_r1 = 0.25,
!>
!>     ! verbosity level (0 = minimal output, 1 = print merged basis indices, 2 = print merged basis overlaps)
!>     ! also controls diagnostic output of GBTOlib
!>     verbosity = 0,
!> /
!> \endverbatim
!>
!> This program cannot be currently used to compute overlaps (and match orbitals) between basis sets that have different
!> B-spline grids.
!>
program phase_match
   use const_gbl, only: set_verbosity_level
   use mpi_gbl, only: mpi_mod_start, mpi_mod_finalize
   use integral_storage_gbl, only: integral_options_obj
   use precisn_gbl, only: wp
   use phase_match_orbitals_mod
   implicit none

   integer, parameter :: max_geom = 800, max_orbs = 200
   logical, parameter :: print_to_stdout = .true.

   integer :: lue
   integer :: i,j,a,b,geom,err
   real(kind=cfp), allocatable :: orbital_olaps(:,:)
   integer, allocatable :: ref_indices(:,:), match_indices(:,:), orbital_phases_one_geom(:), orbital_phases(:,:)

   type(integral_options_obj) :: integral_options

   !Namelist variables
   integer :: n_geom, n_orbitals_to_match, max_l_legendre_1el, mixed_ints_method, verbosity
   integer, allocatable :: orbital_num(:), orbital_sym(:)
   logical :: find_matching_orbitals
   real(kind=cfp) :: rmat_radius, delta_r1

   character(len=132) :: moints(max_geom)
   namelist/INPUT/n_geom,moints,n_orbitals_to_match,orbital_num,orbital_sym,find_matching_orbitals,rmat_radius,max_l_legendre_1el,&
                  mixed_ints_method, delta_r1, verbosity

      call mpi_mod_start(print_to_stdout)

      write(*,'(/,10X,"PHASE-MATCHING OF ORBITALS")')
      write(*,  '(10X,"==========================")')

      allocate(orbital_num(max_orbs*max_geom), orbital_sym(max_orbs*max_geom))

      n_geom = 0
      n_orbitals_to_match = 0
      orbital_num = 0
      orbital_sym =-1
      moints = ''
      find_matching_orbitals = .false.
      rmat_radius = 0
      mixed_ints_method = 0
      max_l_legendre_1el = 0
      delta_r1 = 0
      verbosity = 0

      read(*,nml=input)

      call set_verbosity_level(verbosity)

      write(*,'(/,5X,"Input data")')
      write(*,'(  5X,"----------")')

      write(*,'(/,5X,"Number of geometries: ",i5)') n_geom
      write(*,'(/,5X,"Number of orbitals to phase-match: ",i5)') n_orbitals_to_match

      if (rmat_radius .le. 0) stop "Error: rmat_radius .le. 0"
      if (n_geom .le. 0) stop "Error: n_geom .le. 0"
      if (n_orbitals_to_match .le. 0) stop "Error: n_geom .le. 0"

!      write(*,'(5X,"Indices and symmetries of the orbitals to phase-match:")')
!      do i=1,n_orbitals_to_match
!         write(*,'(5X,3i5)') i,orbital_num(i),orbital_sym(i)
!      enddo

      do i=1,n_geom
         write(*,'(/,5X,"Geometry: ",i5)') i
         write(*,'(5X,"Integrals file: ",a)') trim(adjustl(moints(i)))
      enddo

      write(*,'(/,5X,"End of input data")')
      write(*,'(  5X,"----------")')

      integral_options%a = rmat_radius
      integral_options%max_ijrs_size = 0.0_cfp
      integral_options%calculate_overlap_ints = .true.
      integral_options%calculate_kinetic_energy_ints = .true.
      integral_options%calculate_property_ints = .false.
      integral_options%max_property_l = -1
      integral_options%calculate_nuclear_attraction_ints = .false.
      integral_options%calculate_one_el_hamiltonian_ints = .false.
      integral_options%use_spherical_cgto_alg = .true.
      integral_options%mixed_ints_method = mixed_ints_method
      integral_options%max_l_legendre_1el = max_l_legendre_1el
      integral_options%max_l_legendre_2el = 0
      integral_options%scratch_directory = ''
      integral_options%delta_r1 = delta_r1
      integral_options%two_p_continuum = .false.

      allocate(orbital_phases_one_geom(n_orbitals_to_match),orbital_phases(n_orbitals_to_match,n_geom))

      orbital_phases = 1

      do geom=2,n_geom

         write(*,'(/,5X,"Geometry number: ",i4)') geom

         call load_orbital_overlaps(integral_options,orbital_olaps,moints(geom-1),moints(geom),ref_indices,match_indices)

         !Get the overlaps between the desired orbitals.
         call match_orbitals (orbital_olaps, orbital_phases_one_geom, ref_indices, match_indices, &
                              n_orbitals_to_match, orbital_num, orbital_sym, geom, find_matching_orbitals)

         write(*,'(5X,"Orbital phases for the current pair of geometries:")')
         do i=1,n_orbitals_to_match
            write(*,'(5X,i5,": ",i2)') i,orbital_phases_one_geom(i)
            orbital_phases(i,geom) = orbital_phases(i,geom-1)*orbital_phases_one_geom(i)
         enddo

      enddo !geom

      write(*,'(/,10X,"Phase-correction factors for orbitals:")')

      do geom=1,n_geom
         write(*,'(5X,"Geometry: ",i5)') geom

         do j=1,n_orbitals_to_match
            write(*,'(5X,i5,": ",i2)') j,orbital_phases(j,geom)
         enddo !j

      enddo !geom

!report the final phases for the neutral states: NO2 dyson orbital case only.
!      do geom=1,n_geom
!         write(*,'(2i4)') orbital_phases(1,geom), orbital_phases(3,geom)
!      enddo

      write(*,'(/,10X,"PHASE-MATCHING OF ORBITALS HAS FINISHED")')
      write(*,  '(10X,"=======================================")')

      call mpi_mod_finalize

end program phase_match
