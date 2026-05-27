! Copyright 2025
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

!> \brief   Indexing routines for two-electron integrals over atomic orbitals
!> \authors J Benda
!> \date    2025
!>
!> This module implements routines used for indexing a compact atomic integral storage. The constraints on a good storage scheme
!> of two-electron integrals over atomic orbitals are manyfold:
!>
!>   1. There are three basic permutation symmetries [ab|cd] = [ba|cd] = [ab|dc] = [cd|ab], which reduce the number of integrals
!>      to store by a significant factor.
!>   2. The calculation is likely to assume single-particle continuum, which removes the need to store TCCC and CCCC combinations.
!>      This is a massive reduction of volume. At the same time, two-particle continuum has to be supported as well.
!>   3. Many integrals are zero by symmetry. The product of irreducible representation of all four orbitals has to be the totally
!>      symmetric irreducible representation.
!>   4. As B-splines have compact supports, a pair of B-splines may not overlap at all, further cropping the number of numbers
!>      that have to be stored.
!>
!> If all four sources of sparsity are to be exploited, the layout of the integral storage becomes very complicated. This
!> obviously slows down operations on the storage, but it allows maximal reduction of memory requirements. In the present
!> implementation, the integrals are stored by integral classes:
!>
!>   - TTTT (CGTOs, target orbitals, typically without a definite symmetry)
!>   - CTTT (CGTOs, one continuum orbital)
!>   - CTCT (CGTOs, two continuum orbitals)
!>   - CCTT (CGTOs, two continuum orbitals)
!>   - CCCT (CGTOs, three continuum orbitals; 2-particle continuum only)
!>   - CCCC (CGTOs, four continuum orbitals; 2-particle continuum only)
!>   - BTTT (one BTO continuum orbital, three target CGTOs)
!>   - BTCT (one BTO continuum orbital, two target CGTOs, one continuum CGTO)
!>   - BTCC (one BTO continuum orbital, one target CGTO, two continuum CGTOs; 2-particle continuum only)
!>   - BTBT (two BTO continuum orbitals, two target CGTOs)
!>   - BCBT (two BTO continuum orbitals, one target CGTO, one continuum CGTO; 2-particle continuum only)
!>   - BCBC (two BTO continuum orbitals, two continuum CGTOs; 2-particle continuum only)
!>   - BBTT (two BTO continuum orbitals, two target CGTOs)
!>   - BBCT (two BTO continuum orbitals, one target CGTO, one continuum CGTO)
!>   - BBCC (two BTO continuum orbitals, two continuum CGTOs)
!>   - BBBT (three BTO continuum orbitals, one target CGTOs; 2-particle continuum only)
!>   - BBBC (three BTO continuum orbitals, one continuum CGTOs; 2-particle continuum only)
!>   - BBBB (four BTO continuum orbitals; 2-particle continuum only)
!>
!> Each of the above classes is further binned by combinations of irreducible representations. Assuming that Ia, Ib, Ic, Id are the
!> irreducible representations the orbitals a, b, c, d belong into, then the integral storage contains a block of integrals for
!> each possible combination Ia,Ib,Ic,Id. Blocks, for which the combined irreducible representation Ia×Ib×Ic×Id is not the totally
!> symmetric one, are empty and do not occupy any memory. There is one exception, which is the case of undefined irreducible
!> representation (i.e., one off-centre CGTO in a larger point group). This is labeled as irreducible representation 0 and is not
!> subject to the mentioned selection rule. Each of the symmetry groups (and within each integral class) contains some number of
!> orbitals. If the orbitals a,b have the same type (T/C/B), only groups with Ia >= Ib will be non-empty, taking advantage of the
!> permutation symmetry in the first coordinate. If the orbitals c,d have the same type, only groups with Ic >= Id will be
!> non-empty, taking advantage of the permutation symmetry in the second coordinate. Finally, if the type pair of a,b is the same
!> as the type pair of c,d, then only Ia*(Ia - 1)/2 + Ib >= Ic*(Ic - 1)/2 + Id will be present in the storage, taking advantage
!> of the symmetry in the two integration coordinates.
!>
!> Irreducible representation of an orbital is determined in the following way:
!>
!>   - If the orbital is a CGTO, it is first verified that its centre is symmetric with respect to symmetry operations of
!>     the point group. This is easy, because all symmetry operations consist only of simulatenous reflections along coordinate
!>     axes. The centre is symmetric if the coordinates are zero for all these directions. If the centre is not symmetric,
!>     the orbital is marked as having undefined irreducible representation. This will generally hold for many (if not most)
!>     target orbitals.
!>   - In all other cases (symmetry-compatible CGTO, or a B-spline), the irreducible representation is taken from the real
!>     spherical harmonic associated with the orbital. See \ref Xlm_mgvn.
!>
!> The complete algorithm that assigns a compact index to the integral a,b,c,d is then the following:
!>
!>   1. Sort given orbitals by types (T/C/B) to get integral class and offset in integral array where corresponding integrals start.
!>      This can be done by simply sorting the orbitals by their global index, because B-splines always follow after Gaussians,
!>      and the latter are stored T first, C second (if any).
!>   2. Get irreducible representations of the orbitals and, if there is any symmetry in integral types, sort them further to get
!>      canonically ordered irreducible representations.
!>   3. Update offset in the integral array to the beginning of the current symmetry group. If the group is empty (e.g., due to a
!>      selection rule), return zero.
!>   4. If there is a symmetry both in the orbital types and in the irreducible representations, sort the orbitals further to
!>      a canonical order by their global index. This is not needed if orbitals were sorted by global index already in step 1.
!>   5. Find out the non-redundant orbital pair index from the list of precomputed non-redundant orbital indices for a,b as well
!>      as for c,d. These indices are precomputed at initialization of this module and take into account whether the orbitals
!>      overlap at all. (Such an orbital pair will the have precomputed index within the type-symmetry group equal to 0.)
!>   6. Use these two indices to compute the combined non-redundant index of the orbital quartet. Use triangular index if the
!>      type-symmetry pairs are equal, rectangular index otherwise.
!>   7. This is the position with respect to the beginning of the current type-symmetry block.
!>
submodule (atomic_basis_gbl) integral_indexing_gbl

   implicit none

   integer, allocatable :: orbital_type(:)                  !< for each AO: its shell type (1 = G(T), 2 = G(C), 3 = B)
   integer, allocatable :: orbital_irr(:)                   !< for each AO: its irreducible representation (0 if none)

   integer, allocatable :: type_symm_offset(:, :, :, :, :)  !< how many integrals to skip for a given type-symmetry block
   integer, allocatable :: index_in_type_symm_group(:, :)   !< non-redundant indices of two orbitals per irr. pair
   integer, allocatable :: n_type_symm_indices(:, :, :, :)  !< per 2 irreducible representations and 2 orbital types

   logical :: two_p_continuum_included = .false.

contains

   !> \brief   Set up global module data
   !> \authors J Benda
   !> \date    2025
   !>
   !> Precompute auxiliary data used in the indexing routine \ref integral_index.
   !>
   module subroutine init_compact_integral_indexing(ao_basis, two_p_continuum, n_integrals)

      use const_gbl,    only: abel_prod_tab
      use precisn_gbl,  only: cfp
      use symmetry_gbl, only: is_centre_symmetric, Xlm_mgvn

      type(atomic_orbital_basis_obj), intent(in) :: ao_basis
      type(CGTO_shell_data_obj), allocatable :: cgto_shells(:)
      type(BTO_shell_data_obj), allocatable :: bto_shells(:)
      logical, intent(in) :: two_p_continuum
      integer, intent(out) :: n_integrals
      integer, allocatable :: n_orbitals_per_type_symm(:, :), orbital_symm_type_index(:)
      logical, allocatable :: shells_overlap(:, :)
      integer :: n_ao, n_irrs, n_shells, n_orb_types, n_pair_types, n_quartet_types, n_cgto_shells, n_bto_shells
      integer :: a, b, ia, ib, ic, id, ta, tb, tc, td, sa, sb, ts, ls, tsa, tsb, rsa, rsb, nab, ncd, tab, tcd, tabcd, ba, bb, ra, rb
      logical :: overlaps
      real(cfp) :: ca(3)

      if (allocated(orbital_type)) deallocate (orbital_type)
      if (allocated(orbital_irr)) deallocate (orbital_irr)
      if (allocated(type_symm_offset)) deallocate (type_symm_offset)
      if (allocated(index_in_type_symm_group)) deallocate (index_in_type_symm_group)
      if (allocated(n_type_symm_indices)) deallocate (n_type_symm_indices)

      n_ao = ao_basis % number_of_functions
      n_irrs = ao_basis % symmetry_data % get_no_irrep(ao_basis % symmetry_data % get_pg())
      n_shells = ao_basis % n_target_sh + ao_basis % n_cont_sh
      n_orb_types = 3  ! G(T), G(C), B
      n_pair_types = n_orb_types*(n_orb_types + 1)/2
      n_quartet_types = n_pair_types*(n_pair_types + 1)/2

      allocate (orbital_type(n_ao), orbital_irr(n_ao), orbital_symm_type_index(n_ao), shells_overlap(n_shells, n_shells))
      allocate (type_symm_offset(0:n_irrs, 0:n_irrs, 0:n_irrs, 0:n_irrs, n_quartet_types))
      allocate (index_in_type_symm_group(n_ao, n_ao), n_orbitals_per_type_symm(0:n_irrs, n_orb_types))
      allocate (n_type_symm_indices(0:n_irrs, 0:n_irrs, n_orb_types, n_orb_types))

      type_symm_offset = -1
      n_type_symm_indices = 0
      n_orbitals_per_type_symm = 0
      two_p_continuum_included = two_p_continuum

      call ao_basis % get_all_CGTO_shells(cgto_shells, n_cgto_shells)
      if (any(ao_basis % shell_descriptor(1, :) == 2)) then
         call ao_basis % get_all_BTO_shells(bto_shells, n_bto_shells)
      end if

      ! determine irreducible representations of orbitals and relative indices
      do a = 1, n_ao
         ! get shell index and relative orbital index within the shell
         sa = ao_basis % indices_to_shells(1, a)
         ra = ao_basis % indices_to_shells(2, a)
         ! get type and angular momentum of the shell
         ts = ao_basis % shell_descriptor(1, sa)
         ls = ao_basis % shell_descriptor(6, sa)
         ! get orbital type and centre
         if (ts == 1) then
            ta = merge(1, 2, a <= ao_basis % n_target_fns)
            ca = cgto_shells(sa) % center
         else
            ta = 3
            ca = 0
         end if
         ! to have a definite symmetry, the centre has to be symmetric with respect to the user-specified symmetry elements
         if (is_centre_symmetric(ca, ao_basis % symmetry_data % no_sym_op, ao_basis % symmetry_data % sym_op)) then
            ia = 1 + Xlm_mgvn(ao_basis % symmetry_data % no_sym_op, ao_basis % symmetry_data % sym_op, ls, ra - ls - 1)
         else
            ia = 0
         end if
         ! save orbital information for later
         n_orbitals_per_type_symm(ia, ta) = n_orbitals_per_type_symm(ia, ta) + 1
         orbital_irr(a) = ia
         orbital_type(a) = ta
         orbital_symm_type_index(a) = n_orbitals_per_type_symm(ia, ta)
      end do

      ! set up shell overlap information
      do sa = 1, n_shells
         tsa = ao_basis % shell_descriptor(1, sa)
         rsa = ao_basis % shell_descriptor(2, sa)
         do sb = 1, sa
            tsb = ao_basis % shell_descriptor(1, sb)
            rsb = ao_basis % shell_descriptor(2, sb)
            if (tsa == 1 .or. tsb == 1) then
               ! anything has an overlap with CGTO (in principle)
               overlaps = .true.
            else
               ! two BTOs overlap only if support of the radial B-splines has non-empty intersection
               ba = bto_shells(rsa) % bspline_index
               bb = bto_shells(rsb) % bspline_index
               ! WARNING: we assume that the B-spline grid structure is the same for both shells
               overlaps = abs(ba - bb) < bto_shells(1) % bspline_grid % order
            end if
            shells_overlap(sa, sb) = overlaps
            shells_overlap(sb, sa) = overlaps
         end do
      end do

      ! set up compressed orbital pair indices within a type-symmetry group
      do a = 1, n_ao
         do b = 1, n_ao
            ta = orbital_type(a);  ia = orbital_irr(a);  sa = ao_basis % indices_to_shells(1, a);  ra = orbital_symm_type_index(a)
            tb = orbital_type(b);  ib = orbital_irr(b);  sb = ao_basis % indices_to_shells(1, b);  rb = orbital_symm_type_index(b)
            if (.not. shells_overlap(sa, sb)) then
               ! shells do not overlap
               index_in_type_symm_group(a, b) = 0
            else if (ta /= tb .or. ia /= ib) then
               ! mixed-type and/or mixed-symmetry case (use rectangular index)
               n_type_symm_indices(ia, ib, ta, tb) = n_type_symm_indices(ia, ib, ta, tb) + 1
               index_in_type_symm_group(a, b) = n_type_symm_indices(ia, ib, ta, tb)
            else if (ra >= rb) then
               ! symmetric case (use triangular index, drop non-canonical order)
               n_type_symm_indices(ia, ib, ta, tb) = n_type_symm_indices(ia, ib, ta, tb) + 1
               n_type_symm_indices(ib, ia, tb, ta) = n_type_symm_indices(ia, ib, ta, tb)
               index_in_type_symm_group(a, b) = n_type_symm_indices(ia, ib, ta, tb)
            end if
         end do
      end do

      ! set up global offsets of type-symmetry groups in the integral storage
      n_integrals = 0
      ! loop over all non-redundant canonically ordered integral classes
      do ta = 1, n_orb_types
         do tb = 1, ta
            do tc = 1, ta
               do td = 1, min(ta, tc)
                  tab = tri_index(ta, tb)
                  tcd = tri_index(tc, td)
                  tabcd = tri_index(tab, tcd)
                  ! loop over all symmetry classes
                  do ia = 0, n_irrs
                     do ib = 0, n_irrs
                        do ic = 0, n_irrs
                           do id = 0, n_irrs
                              type_symm_offset(ia, ib, ic, id, tabcd) = n_integrals
                              ! skip this symmetry combination due to non-canonical order of symmetries
                              if (ta == tb .and. ia < ib) cycle
                              if (tc == td .and. ic < id) cycle
                              if (ta == tc .and. tb == td .and. (ia < ic .or. (ia == ic .and. ib < id))) cycle
                              ! skip this symmetry combination if this type combination corresponds to unwanted 2-particle continuum
                              if (count([ta, tb, tc, td] > 1) >= 3 .and. .not. two_p_continuum) cycle
                              ! skip this symmetry combination due to selection rule
                              if (ia /= 0 .and. ib /= 0 .and. ic /= 0 .and. id /= 0) then
                                 if (abel_prod_tab(ia, ib) /= abel_prod_tab(ic, id)) cycle
                              end if
                              ! fill this symmetry combination
                              if (ta == tc .and. tb == td .and. ia == ic .and. ib == id) then
                                 ! symmetric case: triangular index
                                 nab = n_type_symm_indices(ia, ib, ta, tb)
                                 n_integrals = n_integrals + nab*(nab + 1)/2
                              else
                                 ! general case: rectangular index
                                 nab = n_type_symm_indices(ia, ib, ta, tb)
                                 ncd = n_type_symm_indices(ic, id, tc, td)
                                 n_integrals = n_integrals + nab*ncd
                              end if
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

   end subroutine init_compact_integral_indexing


   !> \brief   Interchange two integers
   !> \authors J Benda
   !> \date    2025
   !>
   subroutine swap(a, b)

      integer, intent(inout) :: a, b
      integer :: c

      c = b;  b = a;  a = c

   end subroutine swap


   !> \brief   Calculate triangular index of sorted index pair
   !> \authors J Benda
   !> \date    2025
   !>
   !> Assuming that a >= b > 0, return a single compact ("triangular") index that can be used
   !> to index all such two-integer combinations.
   !>
   integer function tri_index(a, b) result(idx)

      integer, intent(in) :: a, b

      idx = a*(a - 1)/2 + b

   end function tri_index


   !> \brief   Calculate rectangular index of general index pair
   !> \authors J Benda
   !> \date    2025
   !>
   !> The input indices are assumed to be 1-based, the output is also 1-based.
   !>
   integer function rect_index(a, b, nb) result(idx)

      integer, intent(in) :: a, b, nb

      idx = (a - 1)*nb + b

   end function rect_index


   !> \brief   Space-efficient two-electron AO integral ordering
   !> \authors J Benda
   !> \date    2025
   !>
   !> For a given quartet of absolute AO indices, this function returns either zero (in case that the integral is trivially
   !> zero due to selection rules or compact orbital overlaps) or index in a compact one-dimensional integral storage.
   !>
   !> Call \ref init_compact_integral_indexing before using this function.
   !>
   integer module function compact_integral_index(a, b, c, d) result(idx)

      use const_gbl, only: abel_prod_tab

      integer, value :: a, b, c, d
      integer :: ira, irb, irc, ird, ta, tb, tc, td, rab, rcd, tabcd

      ! mark integral as "zero by symmetry" (or due to lack of B-spline overlap)
      idx = 0

      ! reorder orbitals to the canonical non-redundant pyramidal two-electron integral index
      if (a < b) call swap(a, b)
      if (c < d) call swap(c, d)
      if (a < c .or. (a == c .and. b < d)) then
         call swap(a, c)
         call swap(b, d)
      end if

      ! get orbital irreducible representations
      ira = orbital_irr(a);  irb = orbital_irr(b)
      irc = orbital_irr(c);  ird = orbital_irr(d)

      ! combination of the four orbitals must be totally symmetric (at least provided that all have a definite symmetry)
      if (ira /= 0 .and. irb /= 0 .and. irc /= 0 .and. ird /= 0) then
         if (abel_prod_tab(ira, irb) /= abel_prod_tab(irc, ird)) return
      end if

      ! get orbital types (already in canonical order)
      ta = orbital_type(a);  tb = orbital_type(b)
      tc = orbital_type(c);  td = orbital_type(d)

      ! skip this integral if two-particle continuum was not set up
      if (.not. two_p_continuum_included) then
         if (count([ta, tb, tc, td] > 1) >= 3) return
      end if

      ! sort symmetries to canonical order as much as possible while keeping orbital type order unchanged
      if (ta == tb .and. ira < irb) then
         call swap(a, b)
         call swap(ira, irb)
      end if
      if (tc == td .and. irc < ird) then
         call swap(c, d)
         call swap(irc, ird)
      end if
      if (ta == tc .and. tb == td) then
         if (ira < irc .or. (ira == irc .and. irb < ird)) then
            call swap(a, c); call swap(ira, irc)
            call swap(b, d); call swap(irb, ird)
         end if
      end if

      ! sort orbitals to canonical order again, keeping the final type-symmetry order
      if (ta == tc .and. tb == td .and. ira == irc .and. irb == ird) then
         if (a < c .or. (a == c .and. b < d)) then
            call swap(a, c)
            call swap(b, d)
         end if
      end if

      ! get effective orbital pair indices (this takes care of permutation symmetry within pairs, as well as of zero overlaps)
      rab = index_in_type_symm_group(a, b)
      rcd = index_in_type_symm_group(c, d)

      ! skip non-overlapping orbital pairs
      if (rab == 0 .or. rcd == 0) return

      ! get pyramidal identifier for this non-redundant combination of orbital types
      tabcd = tri_index(tri_index(ta, tb), tri_index(tc, td))

      ! compute position in the integral array
      if (ta == tc .and. tb == td .and. ira == irc .and. irb == ird) then
         idx = type_symm_offset(ira, irb, irc, ird, tabcd) + tri_index(rab, rcd)
      else
         idx = type_symm_offset(ira, irb, irc, ird, tabcd) + rect_index(rab, rcd, n_type_symm_indices(irc, ird, tc, td))
      end if

   end function compact_integral_index

end submodule integral_indexing_gbl
