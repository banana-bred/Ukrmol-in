! Copyright 2024
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
module gxg_integrals_gbl
   use precisn_gbl
   use basis_data_generic_gbl
   use cgto_pw_expansions_gbl
   use symmetry_gbl
   use common_obj_gbl, only: nucleus_type
   use utils_gbl

   implicit none

   private

   integer :: inp_max_l_cgto = -1, max_l_pw = -1, inp_max_l_bto = -1
   type(CGTO_shell_data_obj), allocatable :: CGTO_shells_saved(:)
   type(BTO_shell_data_obj), allocatable :: BTO_shells_saved(:)
   type(symmetry_obj) :: symmetry_data_saved
   integer, allocatable :: shell_starting_indices_saved(:)

   type(CGTO_shell_pair_pw_expansion_obj), allocatable :: GG_pair_pw_atoms(:)
   type(CGTO_shell_pair_pw_expansion_obj) :: GG_pair_pw_atoms_aux

   !Used for ECPs
   type(CGTO_shell_pw_expansion_obj), allocatable :: G_shells_pw_atoms(:,:)
   real(kind=cfp), allocatable :: semiloc_ecp_atoms(:,:,:), loc_ecp_atoms(:,:)
   type(pw_expansion_obj) :: radial_grid
   type(nucleus_type), allocatable :: nucleus(:)

   public init_GXG_integrals, GGG_shell_integrals, init_G_ECP_G_integrals, ECP_shell_integrals, GBG_shell_integrals

contains

   subroutine init_GXG_integrals(CGTO_shells, BTO_shells, shell_starting_indices, symmetry_data, a, delta_r1)
      implicit none
      type(CGTO_shell_data_obj), intent(in) :: CGTO_shells(:)
      type(BTO_shell_data_obj), optional, intent(in) :: BTO_shells(:)
      integer, intent(in) :: shell_starting_indices(:)
      type(symmetry_obj), intent(in) :: symmetry_data
      real(kind=cfp), intent(in) :: a, delta_r1

      integer :: err, i

            if (a <= 0.0_cfp) then
               call xermsg ('gxg_integrals_gbl', 'init_GXG_integrals', &
                   'GXG integrals can by calculated only when the R-mat radius has &
                    &been explicitly set.', 1, 1)
            endif

            !Initialize the module calculating the CGTO partial wave expansion on the radial grid.
            inp_max_l_cgto = -1
            do i = 1, size(CGTO_shells)
               inp_max_l_cgto = max(inp_max_l_cgto,CGTO_shells(i)%l)
            enddo

            !we use max() to ensure we don't tamper with the possible setting by another initialization routine
            max_l_pw = max(max_l_pw,inp_max_l_cgto) !the <G_A|G_C|G_B> integral is similar to the property integral with L=L_C

            call init_CGTO_pw_expansions_mod(max_l_pw,inp_max_l_cgto)

            if (allocated(CGTO_shells_saved)) deallocate(CGTO_shells_saved)
            allocate(CGTO_shells_saved,source=CGTO_shells,stat=err)

            if (err /= 0) then
               call xermsg ('gxg_integrals_gbl', 'init_GXG_integrals', &
                   'Memory allocation 1 failed.', err, 1)
            endif

            if (present(BTO_shells)) then
               inp_max_l_bto = -1
               do i = 1, size(BTO_shells)
                  inp_max_l_bto = max(inp_max_l_bto,BTO_shells(i)%l)
               enddo

               if (allocated(BTO_shells_saved)) deallocate(BTO_shells_saved)
               allocate(BTO_shells_saved,source=BTO_shells,stat=err)
   
               if (err /= 0) then
                  call xermsg ('gxg_integrals_gbl', 'init_GXG_integrals', &
                      'Memory allocation 2 failed.', err, 1)
               endif
            endif

            symmetry_data_saved = symmetry_data

            if (allocated(GG_pair_pw_atoms)) deallocate(GG_pair_pw_atoms)
            allocate(GG_pair_pw_atoms(symmetry_data%no_nuc),stat=err)
            if (err /= 0) then
               call xermsg ('gxg_integrals_gbl', 'init_GXG_integrals', &
                   'Memory allocation 2 failed.', err, 1)
            endif

            do i = 1, symmetry_data%no_nuc
               call GG_pair_pw_atoms(i) % eval_regular_grid(0.0_cfp, a, delta_r1)

               GG_pair_pw_atoms(i)%cgto_shell_A_index = -1
               GG_pair_pw_atoms(i)%cgto_shell_B_index = -1
            enddo

            !shell_starting_indices should be equivalent to atomic_orbital_basis_obj%shell_descriptor(4,:)
            !Since CGTOs are always positioned as the first type of basis function in the atomic basis
            !the i-th CGTO shell corresponds to the i-th shell in the whole basis whether it contains BTOs or not.
            !If BTOs are in the basis then the "..._saved" array will be larger than CGTO_shells and is used to store
            !the GBG overlap integrals.
            if (allocated(shell_starting_indices_saved)) deallocate(shell_starting_indices_saved)
            allocate(shell_starting_indices_saved,source=shell_starting_indices,stat=err)
            if (err /= 0) then
               call xermsg ('gxg_integrals_gbl', 'init_GXG_integrals', &
                   'Memory allocation 3 failed.', err, 1)
            endif

   end subroutine init_GXG_integrals

   !> Calculates <G_A|G|G_B> integrals for all G functions in the basis. The radial integrals are calculated numerically.
   !> The radial grid always extends from one of the atoms (or CGTO centers) to r=a measured from that center.
   !> In case at least one of the CGTOs is of target type it is assumed that the function fits inside the R-matrix sphere
   !> centered on the CMS. In case all three CGTOs are of continuum type the method ensures that the integral is calculated
   !> from the CMS to the actual R-matrix radius so no tail subtraction is necessary.
   subroutine GGG_shell_integrals(CGTO_shell_A,CGTO_shell_B,A,B,starting_index_A,starting_index_B,&
                                 &bbb_column, int_index, integrals)
      use const_gbl, only: thrs_center_same
      implicit none
      type(CGTO_shell_data_obj), target, intent(in) :: CGTO_shell_A, CGTO_shell_B
      integer, intent(in) :: A,B, starting_index_A, starting_index_B,bbb_column
      !We assume that these two arrays have been allocated to the appropriate dimensions:
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)

      integer :: i, j, nuc, err, n_pair_m, n_m, offset, i_start, i_end
      real(kind=cfp) :: test, center(3)
      type(CGTO_shell_data_obj) :: shifted_shell_A, shifted_shell_B
      real(kind=cfp), allocatable :: ggg_integrals(:,:)

         n_pair_m = (2*CGTO_shell_A%l+1)*(2*CGTO_shell_B%l+1)
         allocate(ggg_integrals(n_pair_m,(2*inp_max_l_cgto+1)),stat=err)
         if (err /= 0) then
            call xermsg ('gxg_integrals_gbl', 'GGG_integrals', &
                         'Memory allocation failed', err, 1)
         endif

         shifted_shell_A = CGTO_shell_A
         shifted_shell_B = CGTO_shell_B

         !Calculate the PW projections of the G_A*G_B pair wrt all nuclei in the molecule
         do j = 1, symmetry_data_saved%no_nuc
            shifted_shell_A%center = CGTO_shell_A%center - symmetry_data_saved%nucleus(j)%center
            shifted_shell_B%center = CGTO_shell_B%center - symmetry_data_saved%nucleus(j)%center
            call GG_pair_pw_atoms(j)%init_CGTO_shell_pair_pw_expansion(shifted_shell_A,A,shifted_shell_B,B)
            call GG_pair_pw_atoms(j)%eval_CGTO_shell_pair_pw_expansion
         enddo !j

         !Loop over the G_C shells from <G_A|G_C|G_B>
         do i = 1, size(CGTO_shells_saved)

            n_m = 2*CGTO_shells_saved(i)%l+1

            !Identify which center this CGTO is sitting on
            nuc = 0
            do j = 1, symmetry_data_saved%no_nuc
               center = CGTO_shells_saved(i)%center - symmetry_data_saved%nucleus(j)%center
               test = dot_product(center,center)

               if (abs(test) < thrs_center_same) then
                  nuc = j
                  exit
               endif
            enddo !j

            if (nuc == 0) then
               call xermsg ('gxg_integrals_gbl', 'GGG_integrals', &
                   'Matching center not found. shells_A and/or _B are not from the basis?', 1, 1)
            endif

            !Use the PW projections for the matching nucleus to calculate the <G_A|G_C|G_B> integrals
            call GG_pair_pw_atoms(nuc)%eval_GGG_integrals(CGTO_shells_saved(i), &
                                                          starting_index_A,starting_index_B,ggg_integrals,int_index)

            !Transfer the GGG integrals into the appropriate place in the shell integral array
            offset = (bbb_column-1)
            i_start = offset + shell_starting_indices_saved(i)
            i_end   = offset + shell_starting_indices_saved(i) + n_m - 1
            integrals(1:n_pair_m,i_start:i_end) = ggg_integrals(1:n_pair_m,1:n_m)
           
         enddo !i

   end subroutine GGG_shell_integrals

   !> Calculates <G_A|B|G_B> integrals for all BTO functions in the basis. The radial integrals are calculated numerically.
   !> The radial grid always extends from the CMS to r=a.
   !> In case at least one of the CGTOs is of target type it is assumed that the function fits inside the R-matrix sphere
   !> centered on the CMS. No tail subtraction is necessary.
   subroutine GBG_shell_integrals(CGTO_shell_A,CGTO_shell_B,A,B,starting_index_A,starting_index_B,&
                                 &bbb_column, int_index, integrals)
      use const_gbl, only: thrs_center_same
      implicit none
      type(CGTO_shell_data_obj), target, intent(in) :: CGTO_shell_A, CGTO_shell_B
      integer, intent(in) :: A,B, starting_index_A, starting_index_B,bbb_column
      !We assume that these two arrays have been allocated to the appropriate dimensions:
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)

      integer :: i, j, nuc, err, n_pair_m, n_m, offset, i_start, i_end, no_cgto_shells
      real(kind=cfp), allocatable :: gbg_integrals(:,:)

         n_pair_m = (2*CGTO_shell_A%l+1)*(2*CGTO_shell_B%l+1)
         allocate(gbg_integrals(n_pair_m,(2*inp_max_l_bto+1)),stat=err)
         if (err /= 0) then
            call xermsg ('gxg_integrals_gbl', 'GBG_integrals', &
                         'Memory allocation failed', err, 1)
         endif

         !Identify the CMS "atom"
         nuc = 0
         do j = 1, symmetry_data_saved%no_nuc
            if (symmetry_data_saved%nucleus(j)%is_continuum()) then
               nuc = j
               exit
            endif
         enddo

         if (nuc == 0) then
            call xermsg ('gxg_integrals_gbl', 'GBG_integrals', &
                'CMS center not found. Error in the calculation setup.', 1, 1)
         endif

         !Calculate the PW projections of the G_A*G_B pair wrt CMS
         call GG_pair_pw_atoms(nuc)%init_CGTO_shell_pair_pw_expansion(CGTO_shell_A,A,CGTO_shell_B,B)
         call GG_pair_pw_atoms(nuc)%eval_CGTO_shell_pair_pw_expansion

         !Number of CGTO shells in the basis = number of shells preceding BTOs in shell_starting_indices_saved
         no_cgto_shells = size(CGTO_shells_saved)

         offset = (bbb_column-1)

         !Loop over the BTO shells from <G_A|B|G_B>
         do i = 1, size(BTO_shells_saved)

            n_m = 2*BTO_shells_saved(i)%l+1

            !Use the PW projections for CMS to calculate the <G_A|B|G_B> integrals
            call GG_pair_pw_atoms(nuc)%eval_GBG_integrals(BTO_shells_saved(i),&
                                                          starting_index_A,starting_index_B,gbg_integrals,int_index)

            !Transfer the GBG integrals into the appropriate place in the shell integral array. Here we make sure that
            !the integrals for BTOs are appended after the integrals for CGTOs, respecting the order of the shells from
            !the atomic_basis_obj.
            i_start = offset + shell_starting_indices_saved(no_cgto_shells + i)
            i_end   = offset + shell_starting_indices_saved(no_cgto_shells + i) + n_m - 1

            integrals(1:n_pair_m,i_start:i_end) = gbg_integrals(1:n_pair_m,1:n_m)
           
         enddo !i

   end subroutine GBG_shell_integrals

   subroutine init_G_ECP_G_integrals(CGTO_shells, shell_starting_indices, symmetry_data, a, delta_r1)
      use const_gbl, only: level3
      implicit none
      type(CGTO_shell_data_obj), intent(in) :: CGTO_shells(:)
      integer, intent(in) :: shell_starting_indices(:)
      type(symmetry_obj), intent(in) :: symmetry_data
      real(kind=cfp), intent(in) :: a, delta_r1

      integer :: err, i, l, max_bspline_l, max_l_legendre, n, l_loc, nuc
      real(kind=cfp), allocatable :: aux(:), aux_ecp(:,:), aux_ecp_loc(:)
      type(CGTO_shell_data_obj) :: aux_CGTO_shell

            if (a <= 0.0_cfp) then
               call xermsg ('gxg_integrals_gbl', 'init_G_ECP_G_integrals', &
                   'ECP integrals can by calculated only when the R-mat radius has &
                    &been explicitly set.', 1, 1)
            endif

            !Initialize the module calculating the CGTO partial wave expansion on the radial grid.
            inp_max_l_cgto = -1
            do i = 1, size(CGTO_shells)
               inp_max_l_cgto = max(inp_max_l_cgto,CGTO_shells(i)%l)
            enddo

            do i = 1, symmetry_data%no_nuc
               if (symmetry_data%nucleus(i)%ecp%n_core > 0) then
                  l = max(symmetry_data%nucleus(i)%ecp%lmax,symmetry_data%nucleus(i)%ecp%lmaxp)
                  !we use max() to ensure we don't tamper with the possible setting by another initialization routine
                  max_l_pw = max(max_l_pw,l) !the <G_A|ECP_C|G_B> integral is similar to the property integral with L=L_C
               endif
            enddo

            write(level3,'(/," Maximum L in ECP potentials: ",i0)') l

            call init_CGTO_pw_expansions_mod(max_l_pw,inp_max_l_cgto)

            allocate(nucleus,source=symmetry_data%nucleus,stat=err)
            if (err /= 0) then
               call xermsg ('gxg_integrals_gbl', 'init_G_ECP_G_integrals', &
                   'Memory allocation 1 failed.', err, 1)
            endif

            write(level3,'(/," Delta_r for CGTO/CGTO ECP quadratures: ",e25.15)') delta_r1
            write(level3,'(  " Radial range of the quadratures: ",e25.15)') a

            call radial_grid % eval_regular_grid(0.0_cfp, a, delta_r1)
            n = radial_grid % n_total_points

            if (allocated(G_shells_pw_atoms)) deallocate(G_shells_pw_atoms)
            if (allocated(semiloc_ecp_atoms)) deallocate(semiloc_ecp_atoms)
            if (allocated(loc_ecp_atoms)) deallocate(loc_ecp_atoms)

            i = size(CGTO_shells)
            allocate(G_shells_pw_atoms(symmetry_data%no_nuc,i),aux(1),semiloc_ecp_atoms(n,0:l,symmetry_data%no_nuc),&
                     loc_ecp_atoms(n,symmetry_data%no_nuc),stat=err)
            if (err /= 0) then
               call xermsg ('gxg_integrals_gbl', 'init_G_ECP_G_integrals', &
                   'Memory allocation 2 failed.', err, 1)
            endif

            !Precalculate CGTO pw expansion for all functions in the basis with respect to all different ECP centers 
            !and the maximum L from the ECP.

            aux = 0.0_cfp
            !set the other angular momenta to zero: these are in fact added to 'l' by eval_CGTO_shell_pw_expansion.
            max_bspline_l = 0
            max_l_legendre = 0

            do nuc = 1, symmetry_data%no_nuc
               if (symmetry_data%nucleus(nuc)%ecp%n_core > 0) then

                  if (symmetry_data%nucleus(nuc)%ecp%have_ecp) then

                     call symmetry_data%nucleus(nuc)%ecp%eval_local_ecp(radial_grid % r_points, aux_ecp_loc)
                     loc_ecp_atoms(1:n,nuc) = aux_ecp_loc(1:n)

                     if (symmetry_data%nucleus(nuc)%ecp%have_semi_local_ecp) then
                        call symmetry_data%nucleus(nuc)%ecp%eval_semi_local_ecp(radial_grid % r_points, aux_ecp)
                        l_loc = symmetry_data%nucleus(nuc)%ecp%lmax-1
                        semiloc_ecp_atoms(1:n,0:l_loc,nuc) = aux_ecp(1:n,0:l_loc)
                     endif

                  endif

                  if (symmetry_data%nucleus(nuc)%ecp%have_so_ecp) then
                     call xermsg ('gxg_integrals_gbl', 'init_G_ECP_G_integrals', &
                                  'SO ECPs not implemented yet.', 3, 1)
                  endif

                  do i = 1, size(CGTO_shells)

                     !SHIFT THE CENTER OF THE CGTO TO BE WITH RESPECT TO THE I-TH ATOM
                     aux_CGTO_shell = CGTO_shells(i)
                     aux_CGTO_shell % center = aux_CGTO_shell % center - symmetry_data%nucleus(nuc)%center
                     
                     call G_shells_pw_atoms(nuc, i) % init_CGTO_shell_pw_expansion(aux_CGTO_shell, i)
                     call G_shells_pw_atoms(nuc, i) % assign_grid(radial_grid % r_points, radial_grid % weights)

                     !Evaluate the pw expansion with respect to the i-th atom
                     call G_shells_pw_atoms(nuc, i) % eval_CGTO_shell_pw_expansion(aux, &
                                                                                  max_bspline_l, &
                                                                                  l, &
                                                                                  max_l_legendre)
                  enddo !i

               endif
            enddo !nuc

            !shell_starting_indices should be equivalent to atomic_orbital_basis_obj%shell_descriptor(4,:)
            !Since CGTOs are always positioned as the first type of basis function in the atomic basis
            !the i-th CGTO shell corresponds to the i-th shell in the whole basis whether it contains BTOs or not.
            if (allocated(shell_starting_indices_saved)) deallocate(shell_starting_indices_saved)
            allocate(shell_starting_indices_saved,source=shell_starting_indices,stat=err)
            if (err /= 0) then
               call xermsg ('gxg_integrals_gbl', 'init_G_ECP_G_integrals', &
                   'Memory allocation 3 failed.', err, 1)
            endif

            call GG_pair_pw_atoms_aux % eval_regular_grid(0.0_cfp, a, delta_r1)

   end subroutine init_G_ECP_G_integrals

   !> Calculates <G_A|ECP|G_B> integrals for a pair of CGTO functions. The radial integrals are calculated numerically.
   !> The radial grid always extends from the atoms on which the ECP is centered to r=a measured from that center.
   !> In case at least one of the CGTOs is of target type it is assumed that the function fits inside the R-matrix sphere
   !> centered on the CMS. In case both CGTOs are of continuum type we assume naturally that the ECPs are fully contained
   !> inside the sphere ensuring that the integrand vanished before the boundary.
   subroutine ECP_shell_integrals(CGTO_shell_A,CGTO_shell_B,A,B,starting_index_A,starting_index_B,&
                                 &ecp_column, integrals, int_index)
      use const_gbl, only: thrs_center_same, int_rel_prec
      use phys_const_gbl, only: fourpi
      use cgto_hgp_gbl, only: index_1el
      use eri_sph_coord_gbl, only: atomic_ecp_GG_int
      implicit none
      type(CGTO_shell_data_obj), target, intent(in) :: CGTO_shell_A, CGTO_shell_B
      integer, intent(in) :: A,B, starting_index_A, starting_index_B,ecp_column
      !We assume that these two arrays have been allocated to the appropriate dimensions:
      integer, allocatable :: int_index(:,:)
      real(kind=cfp), allocatable :: integrals(:,:)

      integer :: i, j, lA, lB, mA, mB, jA, jB, terms, l_loc, n, n_mA, n_mB, ind, l, m, lm, nuc, pA, pB, ind_AB
      real(kind=cfp) :: cAB, testA(3), testB(3)
      logical, parameter :: calculate_indices = .false.
      real(kind=cfp) :: rad_ecp_int
      real(kind=cfp), allocatable :: cA(:), cB(:)
      type(CGTO_shell_data_obj) :: shifted_shell_A, shifted_shell_B

         lA = CGTO_shell_A%l
         n_mA = (2*lA+1)
         lB = CGTO_shell_B%l
         n_mB = (2*lB+1)
         terms = n_mA * n_mB

         integrals(1:terms,ecp_column) = 0.0_cfp

         !check indices for the integrals
         call index_1el(CGTO_shell_A%l,CGTO_shell_B%l,starting_index_A,starting_index_B,1,int_index,calculate_indices)

         n = radial_grid % n_total_points

         allocate(cA(CGTO_shell_A%number_of_primitives),cB(CGTO_shell_B%number_of_primitives))

         do jA=1,CGTO_shell_A%number_of_primitives
            cA(jA) = sqrt(fourpi/(2*lA+1))*CGTO_shell_A%norm*CGTO_shell_A%contractions(jA)*CGTO_shell_A%norms(jA)
         enddo

         do jB=1,CGTO_shell_B%number_of_primitives
            cB(jB) = sqrt(fourpi/(2*lB+1))*CGTO_shell_B%norm*CGTO_shell_B%contractions(jB)*CGTO_shell_B%norms(jB)
         enddo

         !Cycle over the nuclei which have an ECP and accumulate the ECP integral for the given pair of shells
         do nuc = 1, size(nucleus)

            if (nucleus(nuc)%ecp%have_ecp) then

               if (nucleus(nuc)%ecp%have_so_ecp) then
                  call xermsg ('gxg_integrals_gbl', 'ECP_shell_integrals', &
                               'SO ECPs not implemented yet.', 1, 1)
               endif

               !SINGLE-CENTER (ATOMIC, i.e. ANALYTIC) CASE
               testA = CGTO_shell_A%center-nucleus(nuc)%center
               testB = CGTO_shell_B%center-nucleus(nuc)%center

               if (dot_product(testA,testA) < thrs_center_same .and. dot_product(testB,testB) < thrs_center_same) then

                  !Selection rule for the atomic ECP integrals
                  if (CGTO_shell_A%l /= CGTO_shell_B%l) cycle

                  rad_ecp_int = 0.0_cfp
                  do jA=1,CGTO_shell_A%number_of_primitives
                     do jB=1,CGTO_shell_B%number_of_primitives
         
                        cAB = cA(jA)*cB(jB)
         
                        !V_loc
                        do j=1,nucleus(nuc)%ecp%nlmax
                           rad_ecp_int = rad_ecp_int + cAB*nucleus(nuc)%ecp%cloc(j)&
                                         *atomic_ecp_GG_int(lA,lB,nucleus(nuc)%ecp%mloc(j),CGTO_shell_A%exponents(jA),&
                                                                  CGTO_shell_B%exponents(jB),nucleus(nuc)%ecp%gloc(j))
                        enddo !j
         
                        !The semi-local term contributes only if there is a matching projector
                        if (nucleus(nuc)%ecp%have_semi_local_ecp .and. (lA <= nucleus(nuc)%ecp%lmax-1)) then
                           do j=1,nucleus(nuc)%ecp%nl(lA)
                              rad_ecp_int = rad_ecp_int + cAB*nucleus(nuc)%ecp%c(j,lA)*&
                                            atomic_ecp_GG_int(lA,lB,nucleus(nuc)%ecp%m(j,lA),CGTO_shell_A%exponents(jA),&
                                                                   CGTO_shell_B%exponents(jB),nucleus(nuc)%ecp%g(j,lA))
                           enddo !j
                        endif
            
                     enddo !jB
                  enddo !jA
         
                  i = 0
                  do mB=-lB,lB
                    do mA=-lA,lA
                       i = i + 1
                       if (mA == mB) then
                          integrals(i,ecp_column) = integrals(i,ecp_column) + rad_ecp_int
                       endif
                    enddo
                  enddo

                  cycle !go to the next atom
               endif

               !MULTI-CENTER CASE

               !FIRST EVALUATE THE CONTRIBUTION OF THE LOCAL POTENTIAL: CALCULATE SQRT(4*PI)*<X_{00}|G_{A}G_{B}>_{R}
               if (any(abs(loc_ecp_atoms(1:n,nuc)) > int_rel_prec)) then

                  !1. Shift the center of coordinates to the center of the ECP.
                  shifted_shell_A = CGTO_shell_A
                  shifted_shell_B = CGTO_shell_B
                  shifted_shell_A%center = CGTO_shell_A%center - nucleus(nuc)%center
                  shifted_shell_B%center = CGTO_shell_B%center - nucleus(nuc)%center
   
                  !2. Initialize the pw expanstion of the GG pair.
                  GG_pair_pw_atoms_aux%cgto_shell_A_index = -1
                  GG_pair_pw_atoms_aux%cgto_shell_B_index = -1
                  call GG_pair_pw_atoms_aux%init_CGTO_shell_pair_pw_expansion(shifted_shell_A,A,shifted_shell_B,B)
   
                  !3. Calculate <X_{00}|G_{A}G_{B}>_{r} on the radial grid.
                  call GG_pair_pw_atoms_aux%eval_CGTO_shell_pair_pw_expansion
   
                  !4. Integrate with the ECP
   
                  lm = 1 !index of the X_{00} projector
                  ind = 0
   
                  !For all m-values of the CGTO B
                  do mB = -lB,lB
                     !For all m-values of the CGTO A
                     do mA = -lA,lA
         
                        ind = ind + 1
                        ind_AB = mB+lB+1 + n_mB*(mA+lA)
   
                        !The radial integal extends from the center of the current atom to distance r=a away from it.
                        integrals(ind,ecp_column) = integrals(ind,ecp_column) + &
                                sqrt(fourpi) * sum(loc_ecp_atoms(1:n,nuc) * GG_pair_pw_atoms_aux%angular_integrals(1:n,ind_AB,lm) *&
                                                   radial_grid % r_points(1:n)**2 * radial_grid % weights(1:n))
   
                     enddo !mA
                  enddo !mB
               endif

               if (nucleus(nuc)%ecp%have_semi_local_ecp) then
                  l_loc = nucleus(nuc)%ecp%lmax-1
   
                  !Contribution of the semi-local potential
                  do l = 0,l_loc
                     do m = -l,l
      
                        !index of the L,M projector from the ECP
                        lm = l*l+l+m+1
      
                        ind = 0
                        !For all m-values of the CGTO B
                        do mB = 1,n_mB
                           pB = G_shells_pw_atoms(nuc, B)%non_neg_indices_l(mB,lm)
            
                           !For all m-values of the CGTO A
                           do mA = 1,n_mA
                              pA = G_shells_pw_atoms(nuc, A)%non_neg_indices_l(mA,lm)
      
                              ind = ind + 1
            
                              if (pA /= 0 .and. pB /= 0) then
                                 
                                 !<A|ECP_{LM}^{semi-local}|B> = \int dr r^2 <A|LM>ECP^{semi-local}(R)<LM|B>, where |LM> are real spherical harmonics.
                                 !The radial integal extends from the center of the current atom to distance r=a away from it.
                                 integrals(ind,ecp_column) = integrals(ind,ecp_column) + &
                                          sum(semiloc_ecp_atoms(1:n,l,nuc) *&
                                          G_shells_pw_atoms(nuc, B) % angular_integrals(1:n,pB) * &
                                          G_shells_pw_atoms(nuc, A) % angular_integrals(1:n,pA) * &
                                          radial_grid % r_points(1:n)**2 * radial_grid % weights(1:n))
                              endif                        
                               
                           enddo !mA
                        enddo !mB
                     enddo !m
                  enddo !l

               endif !have_semi_local_ecp

            endif !have_ecp

          enddo !nuc

   end subroutine ECP_shell_integrals

end module gxg_integrals_gbl
