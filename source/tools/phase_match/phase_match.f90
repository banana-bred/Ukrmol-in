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
module phase_match_mod

   use precisn_gbl,   only: wp, cfp
   !use cdenprop_io
   use cdenprop_defs, only: CIvect, CSFheader, CSFbody
   use cdenprop_aux, only: makemg
   use omp_lib

contains

   !Use the overlap integrals between two sets of wavefunctions to match the wavefunctions between the two sets, phase correct them and put them into the corresponding columns of CI vectors set.
   subroutine match_wavefunctions(wfs_olap,ci_vec_match,CV,ei)
      implicit none
      type(CIvect), intent(inout) :: ci_vec_match
      real(kind=wp), allocatable :: wfs_olap(:,:), CV(:,:), ei(:)

      real(kind=wp), allocatable :: wfs_olap_tmp(:,:)
      integer, allocatable :: match_ind(:)
      integer, allocatable :: final_map(:,:)
      integer :: nstat, a, b, phase, i, j, nocsf
      real(kind=wp) :: s

         nocsf = size(ci_vec_match%CV,1)
         nstat = size(ci_vec_match%CV,2)
         if (nocsf .le. 0 .or. nstat .le. 0) then
            print *,nocsf,nstat
            stop "match_wavefunctions: bad CI vector not allocated"
         endif

         if (allocated(CV)) deallocate(CV)
         if (allocated(ei)) deallocate(ei)
         allocate(wfs_olap_tmp(nstat,nstat),CV(nocsf,nstat),ei(nstat),match_ind(nstat),final_map(nstat,nstat))

         write(*,'(/,5X,"Constructing table of mappings between the states for the two geometries...")')
         match_ind = 0
         wfs_olap_tmp = wfs_olap
         do a=1,nstat
            s = maxval(abs(wfs_olap_tmp(1:nstat,a)))
            do b=1,nstat
               if (abs(wfs_olap_tmp(b,a)) .eq. s) then
                  match_ind(a) = b
                  wfs_olap_tmp(b,:) = 0.0_wp !this ensures this b-state never gets matched by two different a-states.
                  exit
               endif
            enddo !b
         enddo !a

         write(*,'(/,5X,"Matching of the target wavefunctions a -> b, relative phase, overlap")')
         final_map = 0
         do a=1,nstat
            b = match_ind(a)
            phase = nint(sign(1.0_wp,wfs_olap(b,a)))
            write(*,'(5X,i5," -> ",i5,",",i2,",",e25.15)') a,b,phase,wfs_olap(b,a)
            !Apply the phase-correction to the corresponding target vector and store it in the column corresponding to the a-state.
            !note that in case the matching for two states should not be done we still multiply a phase factor but that is harmeless...

            !Flip the phase of the CI vectors and permute their order to the canonical order, i.e. order given by the 1st geometry.
            CV(:,a) = phase*ci_vec_match%CV(:,b)
            ei(a) = ci_vec_match%ei(b)

            !Keep the order of the CI vectors the same and only flip the phase of the CI vectors.
            ci_vec_match%CV(:,b) = phase*ci_vec_match%CV(:,b)

            final_map(a,b) = final_map(a,b) + 1
         enddo !a

         do a=1,nstat
            i = sum(final_map(a,:))
            j = sum(final_map(:,a))
            if (i .ne. 1 .or. j .ne. 1) then
               print *,a,i,j
               stop "Error in WF mapping"
            endif
         enddo !a

   end subroutine match_wavefunctions

   !Use the expanded determinants, orbital overlaps and CI vectors to calculate overlaps between the two sets of CI vectors.
   subroutine calculate_wavefunction_overlaps (ndtr_ref, csf_head_ref, csf_body_ref, ndtr_match, csf_head_match, csf_body_match, &
                                               ci_vec_ref, ci_vec_match, spin_orb_olaps, wfs_olap, nstat)
      implicit none
      type(CSFheader), intent(in) :: csf_head_ref, csf_head_match 
      type(CSFbody), intent(in)   :: csf_body_ref, csf_body_match
      type(CIvect), intent(in) :: ci_vec_ref,ci_vec_match
      integer, allocatable :: ndtr_ref(:,:), ndtr_match(:,:)
      real(kind=wp), allocatable :: spin_orb_olaps(:,:)
      real(kind=wp), allocatable :: wfs_olap(:,:)
      integer, intent(out) :: nstat

      integer :: num_threads,iam,i,j,k,a,b,det_i,det_j,no_spin_orbitals_ref,n_dets_ref,no_spin_orbitals_match,n_dets_match
      real(kind=wp) :: olap,det
      real(kind=wp), allocatable :: wfs_olap_thread(:,:,:)
      logical :: error1, error2

         n_dets_ref = size(ndtr_ref,2)
         no_spin_orbitals_ref=2*sum(csf_head_ref%nob)
         n_dets_match = size(ndtr_match,2)
         no_spin_orbitals_match=2*sum(csf_head_match%nob)
         nstat = size(ci_vec_ref%CV,2)

         error1 =  (n_dets_ref .le. 0 .or. no_spin_orbitals_ref .ne. size(spin_orb_olaps,1))
         error2 =  (n_dets_match <= 0 .or. no_spin_orbitals_match /= size(spin_orb_olaps,2) .or. nstat /= size(ci_vec_match%CV,2))

         if (error1 .or. error2) then
            print *,error1,error2
            stop "calculate_wavefunction_overlaps: error in input data"
         endif

         error1 = csf_head_ref%NELT .ne. csf_head_match%NELT
         if (error1) then
            print *,csf_head_ref%NELT,csf_head_match%NELT
            stop "calculate_wavefunction_overlaps: number of electrons differs between the two sets of CSFs"
         endif

         if (allocated(wfs_olap)) deallocate(wfs_olap)
         allocate(wfs_olap(nstat,nstat))

         write(*,'(/,5X,"Calculating wavefunction overlaps...")')
         !$OMP PARALLEL DEFAULT(NONE) PRIVATE(num_threads,iam,i,j,k,olap,a,det_i,b,det_j,det) &
         !$OMP& SHARED(nstat,csf_head_ref,csf_body_ref,ndtr_ref,csf_head_match,csf_body_match,ndtr_match, &
         !$OMP&        n_dets_ref,n_dets_match,spin_orb_olaps,no_spin_orbitals_ref,no_spin_orbitals_match,&
         !$OMP&        ci_vec_ref,ci_vec_match,wfs_olap,wfs_olap_thread)
         num_threads = omp_get_num_threads()
         iam = omp_get_thread_num()
         !$OMP SINGLE
         allocate(wfs_olap_thread(nstat,nstat,num_threads))
         wfs_olap_thread = 0.0_wp
         !$OMP END SINGLE
         do i=1,csf_head_ref%nocsf
            do j=1,csf_head_match%nocsf
               k = j+csf_head_match%nocsf*(i-1)
               if (mod(k,num_threads) .ne. iam) cycle !Work redistribution
               olap = 0.0_wp
               !print *,'n',csf_head%nodo(i),csf_head%nodo(j)
               do a=1,csf_head_ref%nodo(i)
                  det_i = csf_body_ref%icdo(i)+a-1
                  !write(*,'("NDTR I",i10,100i3)') det_i,ndtr(1:csf_head%NELT,det_i)
                  do b=1,csf_head_match%nodo(j)
                     det_j = csf_body_match%icdo(j)+b-1
                     !write(*,'("NDTR J",i10,100i3)') det_j,ndtr(1:csf_head%NELT,det_j)
                     call overlap_determinant_pair (ndtr_ref, ndtr_match, det_i, det_j, csf_head_ref % NELT, n_dets_ref, &
                                                    n_dets_match, spin_orb_olaps, no_spin_orbitals_ref, no_spin_orbitals_match, det)
                     if (det .ne. 0.0_wp) then
                        !write(*,'("det",2i4,3e25.15)') a,b,csf_body%cdo(det_i),csf_body%cdo(det_j),det
                        olap = olap + det*csf_body_ref%cdo(det_i)*csf_body_match%cdo(det_j)
                     endif
                  enddo
               enddo
               !if (j .eq. i) write(*,'(5X,"diag csf olap",2i6,e25.15)') j,i,olap
               do a=1,nstat
                  do b=1,nstat
                     wfs_olap_thread(b,a,iam+1) = wfs_olap_thread(b,a,iam+1) + olap*ci_vec_ref%CV(i,a)*ci_vec_match%CV(j,b)
                  enddo
               enddo
            enddo
         enddo
         !$OMP BARRIER
         !$OMP SINGLE
         wfs_olap = 0.0_wp
         do i=1,num_threads
            wfs_olap(1:nstat,1:nstat) = wfs_olap(1:nstat,1:nstat) + wfs_olap_thread(1:nstat,1:nstat,i)
         enddo
         deallocate(wfs_olap_thread)
         !$OMP END SINGLE
         !$OMP END PARALLEL

   end subroutine calculate_wavefunction_overlaps

   !Take the determinants in the packed form as output from CONGEN and expand them into full determinants.
   subroutine expand_determinants(csf_head,csf_body,ndtr)
      implicit none
      type(CSFheader), intent(in) :: csf_head
      type(CSFbody), intent(in)   :: csf_body
      integer, allocatable :: ndtr(:,:)

      integer, allocatable :: ndtrf_point_back(:)
      integer :: no_spin_orbitals, n_closed, i, j, k, ii, n_ex, index_ndo_start_i, index_cdo_start_i, ir, n_dets
      
      !MAKEMG constructs an array that indexes the reference determinant
      !array. Given an orbital number it holds the position of that
      !orbital in the reference determinant.

      no_spin_orbitals=2*sum(csf_head%nob)
      allocate(ndtrf_point_back(no_spin_orbitals))

      call MAKEMG(ndtrf_point_back,no_spin_orbitals,csf_head%NELT,csf_head%ndtrf)

      write(*,'(5X,"Reference determinant: ",200i3)') (csf_head%ndtrf(i),i=1,csf_head%NELT)

      n_dets = sum(csf_head%nodo)
      if (allocated(ndtr)) deallocate(ndtr)
      allocate(ndtr(csf_head%NELT,n_dets))

      write(*,'(/,5X,"Number of determinants: ",i15)') n_dets
      write(*,'(  5X,"Expanding all determinants in full...")')
      !Process the CSFs: expand each determinant in full and save it in the ndtr array.
      n_closed = csf_head%NELT
      do i=1,csf_head%nocsf
         index_ndo_start_i=csf_body%indo(i) !starting index in array NDO: number of excitations for each determinant
         index_cdo_start_i=(csf_body%icdo(i)) !starting index in array CDO: coefficients for each determinant
         j = index_ndo_start_i
         do k=1,csf_head%nodo(i) !over all determinants for CSF i
            n_ex = csf_body%ndo(j)
            ii=j
            !Replace the spin-orbitals in the reference orbital with the spin-orbitals building the current determinant:
            ndtr(:,index_cdo_start_i+k-1) = csf_head%ndtrf
            do j=ii+1,ii+n_ex
               !Spin-orbital csf_body%ndo(j) from the reference determinant has been replaced with spin-orbital csf_body%ndo(j+n_ex)
               ir = ndtrf_point_back(csf_body%ndo(j))
               ndtr(ir,index_cdo_start_i+k-1) = csf_body%ndo(j+n_ex)
            enddo
            !ndtr lists all spin-orbitals comprising the determinant
            do j=1,csf_head%NELT
               if (ndtr(j,index_cdo_start_i+k-1) .ne. csf_head%ndtrf(j)) then
                  n_closed = min(j,n_closed)
               endif
            enddo
            j = ii+2*n_ex+1
         enddo
      enddo
      write(*,'(5X,"Number of closed spin-orbitals: ",i5)') n_closed

   end subroutine expand_determinants

   subroutine map_spin_orbitals(csf_head,orb_data)
      implicit none
      type (CSFheader), intent(in) :: csf_head
      integer, allocatable :: orb_data(:,:)

      integer :: i,j,k,no_spin_orbitals,ind,cnt

         no_spin_orbitals=2*sum(csf_head%nob)
         if (allocated(orb_data)) deallocate(orb_data)
         allocate(orb_data(3,no_spin_orbitals))

         write(*,'(/,5X,"Spin-orbital table: spin-orbital index, orbital number within symmetry, spin, orbital symmetry")')
         write(*,'(  5X,"----------------------------------------------------------------------------------------------")')
         ind = 0
         cnt = 0
         do i=1,csf_head%nsym
            do j=1,csf_head%nob(i)
               cnt = cnt + 1
               do k=0,1
                  ind = ind + 1
                  orb_data(1:3,ind) = (/j,k,i/)
                  write(*,'(5X,i8,4i5)') ind,orb_data(1:3,ind)
               enddo
            enddo
         enddo
   
   end subroutine map_spin_orbitals

   !Merge molecular data for two geometries into one and calculate the overlap integrals between the two sets of orbitals.
   subroutine load_spin_orbital_overlaps(orb_data_ref,orb_data_match,spin_orb_olaps,moints_ref,moints_match)
      use mpi_gbl, only: mpi_xermsg
      use utils_gbl, only: xermsg
      use const_gbl, only: overlap_ints
      use basis_data_generic_gbl, only: CGTO_shell_data_obj, orbital_data_obj
      use atomic_basis_gbl, only: atomic_orbital_basis_obj
      use molecular_basis_gbl, only: molecular_orbital_basis_obj
      use parallel_arrays_gbl, only: p2d_array_obj
      use integral_storage_gbl, only: integral_options_obj, integral_storage_obj
      use symmetry_gbl, only: geometry_obj
      implicit none
      integer, intent(in) :: orb_data_ref(:,:), orb_data_match(:,:)
      real(kind=wp), allocatable :: spin_orb_olaps(:,:)
      character(len=132), intent(in) :: moints_ref, moints_match

      integer :: i, j, k, no_spin_orbitals_ref, no_spin_orbitals_match, nshells_ref, nshells_match, n_cfs, n_orbs, nshells_merged, &
                 nnuc_ref, nnuc_match, nnuc_merged, err, cnt, ind(1:1), two_ind(2,1), n_ref(8), col
      logical :: same, have_btos
      type(atomic_orbital_basis_obj), target :: atomic_orbital_basis_ref, atomic_orbital_basis_match, atomic_orbital_basis_merged
      type(molecular_orbital_basis_obj) :: molecular_orbital_basis_ref, molecular_orbital_basis_match,molecular_orbital_basis_merged
      type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data_ref(:), CGTO_shell_data_match(:), CGTO_shell_data_merged(:)
      type(geometry_obj) :: geometry_ref, geometry_match, geometry_merged
      type(orbital_data_obj) :: orbital_data, orbital_data_tmp
      type(integral_options_obj) :: integral_options
      type(integral_storage_obj), target :: atomic_integral_storage, molecular_integral_storage
      type(p2d_array_obj), target :: ao_integrals, mo_integrals

         write(*,'(/,5X,"Constructing spin-orbital overlap matrix")')
         write(*,'(  5X,"----------------------------------------")')
         write(*,'(/,5X,"Reference basis sets on file: ",a)') trim(adjustl(moints_ref))
         write(*,'(  5X,"Matching basis sets on file: ",a)') trim(adjustl(moints_match))

         !REFERENCE SET: read-in all basis sets and orbital data
         call atomic_orbital_basis_ref%read(moints_ref)
         molecular_orbital_basis_ref%ao_basis => atomic_orbital_basis_ref
         call molecular_orbital_basis_ref%read(moints_ref)

         call atomic_orbital_basis_ref%get_all_CGTO_shells(CGTO_shell_data_ref,nshells_ref)
         call molecular_orbital_basis_ref%symmetry_data%get_geometry(geometry_ref)

         write(*,'(5X,"Reference basis set and orbitals have been read-in")')

         !MATCH SET: read-in all basis sets and orbital data
         call atomic_orbital_basis_match%read(moints_match)
         molecular_orbital_basis_match%ao_basis => atomic_orbital_basis_match
         call molecular_orbital_basis_match%read(moints_match)

         call atomic_orbital_basis_match%get_all_CGTO_shells(CGTO_shell_data_match,nshells_match)
         call molecular_orbital_basis_match%symmetry_data%get_geometry(geometry_match)

         write(*,'(5X,"Matching basis set and orbitals have been read-in")')
         write(*,'(5X,"Number of CGTO shells from the reference basis: ",i10)') nshells_ref
         write(*,'(5X,"Number of CGTO shells from the basis to match: ",i10)') nshells_match

         have_btos = atomic_orbital_basis_ref%contains_btos() .or. atomic_orbital_basis_match%contains_btos()
         if (have_btos) then
            print *,'not implemented for btos but can work the same way'
            stop "error in load_spin_orbital_overlaps"
         endif

         nnuc_ref = geometry_ref%no_nuc
         nnuc_match = geometry_ref%no_nuc

         write(*,'(5X,"Number of nuclei from the reference basis: ",i10)') nnuc_ref
         write(*,'(5X,"Number of nuclei from the basis to match: ",i10)') nnuc_match

         !JOIN the geometry data, CGTO basis sets and orbital sets.
         nnuc_merged = nnuc_ref+nnuc_match
         allocate(geometry_merged%nucleus(nnuc_merged))
         geometry_merged%no_nuc = nnuc_merged
         geometry_merged%no_sym_op = geometry_ref%no_sym_op
         geometry_merged%sym_op = geometry_ref%sym_op
         geometry_merged%use_symmetry = geometry_ref%use_symmetry
  
         j = 0
         do i=1,nnuc_ref
            j = j + 1
            geometry_merged%nucleus(j) = geometry_ref%nucleus(i)
            geometry_merged%nucleus(j)%nuc = j
            call geometry_merged%nucleus(j)%print
         enddo
         do i=1,nnuc_match
            j = j + 1
            geometry_merged%nucleus(j) = geometry_match%nucleus(i)
            geometry_merged%nucleus(j)%nuc = j
            call geometry_merged%nucleus(j)%print
         enddo

         write(*,'(/,5X,"Nuclear data have been merged.")')

         nshells_merged = nshells_ref+nshells_match
         n_cfs = atomic_orbital_basis_ref%number_of_functions + atomic_orbital_basis_match%number_of_functions
         n_orbs = molecular_orbital_basis_ref%number_of_functions + molecular_orbital_basis_match%number_of_functions
         allocate(CGTO_shell_data_merged(nshells_merged))

         j = 0
         do i=1,nshells_ref
            j = j + 1
            CGTO_shell_data_merged(j) = CGTO_shell_data_ref(i)
         enddo
         do i=1,nshells_match
            j = j + 1
            CGTO_shell_data_merged(j) = CGTO_shell_data_match(i)
         enddo

         !INITIALIZE the joined basis sets
         err = atomic_orbital_basis_merged%init(nshells_merged,geometry_merged)

         !add the target CGTO shells
         do i=1,nshells_merged
            call atomic_orbital_basis_merged%add_shell(CGTO_shell_data_merged(i))
         enddo

         write(*,'(5X,"Atomic basis sets have been merged.")')

         if (molecular_orbital_basis_ref%no_irr .ne. molecular_orbital_basis_match%no_irr) then
            print *,'number of IRRs differs for the two geometries'
            stop "error in load_spin_orbital_overlaps"
         endif

         molecular_orbital_basis_merged%ao_basis => atomic_orbital_basis_merged
         err = molecular_orbital_basis_merged%init(molecular_orbital_basis_ref%no_irr,geometry_merged)

         orbital_data%number_of_coefficients = n_cfs

         n_ref = 0
         do i=1,molecular_orbital_basis_ref%no_irr
            call molecular_orbital_basis_ref%get_shell_data(i,orbital_data_tmp)

            if (allocated(orbital_data%coefficients)) deallocate(orbital_data%coefficients)
            if (allocated(orbital_data%energy)) deallocate(orbital_data%energy)
            if (allocated(orbital_data%spin)) deallocate(orbital_data%spin)
            if (allocated(orbital_data%occup)) deallocate(orbital_data%occup)
            n_orbs = orbital_data_tmp%number_of_functions*2
            allocate(orbital_data % coefficients(n_cfs, n_orbs), &
                     orbital_data % energy(n_orbs), &
                     orbital_data % spin(n_orbs), &
                     orbital_data % occup(n_orbs))
            orbital_data%number_of_functions = n_orbs
            orbital_data%coefficients = 0.0_cfp
            orbital_data%energy = 1.0_cfp
            orbital_data%spin = 0
            orbital_data%coefficients = 0.0_cfp
            orbital_data%point_group = orbital_data_tmp%point_group
            orbital_data%irr = orbital_data_tmp%irr

            j = orbital_data_tmp%number_of_coefficients
            k = orbital_data_tmp%number_of_functions
            orbital_data%coefficients(1:j,1:k) = orbital_data_tmp%coefficients(1:j,1:k)
            n_ref(orbital_data%irr) = k

            cnt = k

            call molecular_orbital_basis_match%get_shell_data(i,orbital_data_tmp)
            k = orbital_data_tmp%number_of_functions
            orbital_data % coefficients(j + 1 : j + orbital_data_tmp % number_of_coefficients, cnt + 1 : cnt + k) = &
                orbital_data_tmp % coefficients(1 : orbital_data_tmp % number_of_coefficients, 1 : k)

            call molecular_orbital_basis_merged%add_shell(orbital_data)
         enddo

         write(*,'(5X,"Molecular orbital data have been merged.")')
  
         call molecular_orbital_basis_merged%print_orbitals
!
!GENERATE THE OVERLAP INTEGRALS OVER BOTH ORBITAL SETS INCLUDING THE CROSS OVERLAPS
!
         integral_options%a = -1.0_cfp
         integral_options%max_ijrs_size = 0.0_cfp
         integral_options%calculate_overlap_ints = .true.
         integral_options%calculate_kinetic_energy_ints = .true.
         integral_options%calculate_property_ints = .false.
         integral_options%max_property_l = -1
         integral_options%calculate_nuclear_attraction_ints = .false.
         integral_options%calculate_one_el_hamiltonian_ints = .false.
         integral_options%use_spherical_cgto_alg = .true.
         integral_options%mixed_ints_method = 0
         integral_options%max_l_legendre_1el = 0
         integral_options%max_l_legendre_2el = 0
         integral_options%scratch_directory = ''
         integral_options%delta_r1 = 0.0_cfp
         integral_options%two_p_continuum = .false.
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

         no_spin_orbitals_ref=size(orb_data_ref,2)
         no_spin_orbitals_match=size(orb_data_match,2)
         if (allocated(spin_orb_olaps)) deallocate(spin_orb_olaps)
         allocate(spin_orb_olaps(no_spin_orbitals_ref,no_spin_orbitals_match))
         spin_orb_olaps = 0.0_wp

         do i=1,no_spin_orbitals_ref
            two_ind(1,1) = molecular_orbital_basis_merged%get_absolute_index(orb_data_ref(1,i),orb_data_ref(3,i))
            !The j-index should correspond to the second orbital set.
            do j=1,no_spin_orbitals_match
               k = orb_data_match(1,j) + n_ref(orb_data_match(3,j))
               two_ind(2,1) = molecular_orbital_basis_merged%get_absolute_index(k,orb_data_match(3,j))
               ind(1:1) = molecular_orbital_basis_merged%integral_index(overlap_ints,two_ind,integral_options%two_p_continuum)

               !Spin selection rule:
               if (orb_data_ref(2,i) .eq. orb_data_match(2,j)) spin_orb_olaps(i,j) = mo_integrals%a(ind(1),col)

               same = .true.
               do k=1,3
                  if (orb_data_ref(k,i) .ne. orb_data_match(k,j)) same = .false.
               enddo
               if (same) write(*,'(5X,"DIAG spin-orb olap",2i10,e25.15)') i,j,spin_orb_olaps(i,j)
            enddo
         enddo

         write(*,'(/,5X,"Spin-orbital overlap matrix has been constructed.")')

   end subroutine load_spin_orbital_overlaps

   !Calculate overlap between a pair of determinants.
   subroutine overlap_determinant_pair (ndtr_ref, ndtr_match, det_i, det_j, nelt, n_dets_ref, n_dets_match, spin_orb_olaps, &
                                        no_spin_orbitals_ref, no_spin_orbitals_match, det)
      implicit none
      integer, intent(in) :: det_i, det_j, nelt, n_dets_ref, n_dets_match, no_spin_orbitals_ref, no_spin_orbitals_match
      integer, intent(in) :: ndtr_ref(nelt,n_dets_ref), ndtr_match(nelt,n_dets_match)
      real(kind=wp), intent(in) :: spin_orb_olaps(no_spin_orbitals_ref,no_spin_orbitals_match)
      real(kind=wp), intent(out) :: det

      integer :: i,j,info,piv(nelt)
      real(kind=wp) :: a(nelt,nelt)

         a = 0.0_wp
         do j=1,nelt
            do i=1,nelt
               !ndtr(j,det_j) is the orbital from the second orbital set.
               a(i,j) = spin_orb_olaps(ndtr_ref(i,det_i),ndtr_match(j,det_j))
            enddo
         enddo

         !Determinant of the overlap matrix
         call dgetrf(nelt,nelt,a,nelt,piv,info)
         det = 0.0_wp
         if (info.ne.0) then !one of the U(i,i) elements is zero: the determinant is zero too
            !print *,'det=0'
            return
         endif
         det = 1.0_wp
         do i=1,nelt
            if (piv(i).ne.i) then
               det = -det*a(i,i)
            else
               det = det*a(i,i)
            endif
         enddo

   end subroutine overlap_determinant_pair

!*==writcip.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
!modified by ZM not to include reading of geometry from the integrals file.
      SUBROUTINE WRITCIP(nftw,nset,NAME,e0,nocsf,nstat,mgvn,s,sz,nelt,&
     & EIG,VEC,iphase,dgem,ntgsym,mcont,notgs,nnuc,cname,xnuc,ynuc,znuc,charge,&
     & NFT,npflg,npcvc)
      USE params, ONLY : ccidata, c8stars, cblank, cpoly
      USE scatci_data, ONLY : MEIG
      USE SCATCI_ROUTINES, ONLY: movep, search, civio
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      REAL(KIND=wp) :: E0, S, SZ
      INTEGER :: MGVN, NELT, NFT, NFTI, NFTW, NOCSF, NPCVC, NPFLG, NSET, NSTAT, NTGSYM, NNUC
      character(len=120) :: NAME
      REAL(KIND=wp), DIMENSION(*) :: DGEM, EIG, VEC
      INTEGER, DIMENSION(*) :: IPHASE
      INTEGER, DIMENSION(ntgsym) :: MCONT, NOTGS
      REAL(KIND=wp), intent(in) :: XNUC(NNUC),YNUC(NNUC),ZNUC(NNUC),CHARGE(NNUC)
      character(len=8), intent(in) :: CNAME(NNUC)
      INTENT (IN) E0, MCONT, MGVN, NAME, NELT, NOTGS, NPFLG, NTGSYM, S, SZ, NNUC
      INTENT (INOUT) NSET
!
! Local variables
!
      INTEGER :: I, IFAIL, II, NALM, NEIG, NREC, NTH
!
!*** End of declarations rewritten by SPAG
!
!
!**********************************************************************
!
!     WRITCIP writes CI data to unit NFTW in format used for
!     polyatomic targets.  Data is always appended to the end of
!     existing files
!
!**********************************************************************
!
!.... Position file for next set
      NTH=NSET
!
      CALL MOVEP(NFTW,NTH,NALM,NPCVC,NFT)
      IF(NALM.NE.0)THEN
         WRITE(NFT,2940)
 2940    FORMAT(' ERROR POSITIONING FOR OUTPUT OF CI COEFFICIENTS',//)
         STOP
      END IF
!
!
!.... Write header
      WRITE(nftw)c8stars, cblank, cblank, ccidata
      nset=nth
      nrec=nnuc+nstat+1
      WRITE(NFT,3822)Nset
 3822 FORMAT(/,' CI data stored as set number',I3)
      WRITE(NFTw)Nset, nrec, NAME, nnuc, nocsf, nstat, mgvn, s, sz, &
     &           nelt, e0, ntgsym, (mcont(i),notgs(i),i=1,ntgsym)
      DO i=1, nnuc
         WRITE(nftw)cname(i), xnuc(i), ynuc(i), znuc(i), charge(i)
      END DO
!
!     Write CI eigenvalues and eigenvectors
!
      CALL CIVIO(NFTW,0,NOCSF,nstat,EIG,VEC,NALM,iphase,dgem)
      IF(NALM.NE.0)GO TO 2900
!
!...  Print summary of output data
      WRITE(NFT,100)NSET, name
      WRITE(nft,103)mgvn, s, sz, nelt, nnuc
      WRITE(NFT,120)(cname(i),xnuc(i),ynuc(i),znuc(i),charge(i),i=1,nnuc)
      WRITE(NFT,101)nocsf, nstat, e0
      IF(ntgsym.GT.0)WRITE(nft,104)ntgsym, (i,mcont(i),notgs(i),i=1,ntgsym)
      neig=nstat
      IF(npflg.LE.0 .AND. nstat.GT.meig)neig=meig
      WRITE(NFT,102)(EIG(I)+E0,I=1,NEIG)
!
      RETURN
!
 2900 WRITE(NFT,2960)
 2960 FORMAT(' ERROR WRITING CI OUTPUT FILE',//)
      STOP
!
 100  FORMAT(/' SET',I4,4X,A)
 101  FORMAT(/' NOCSF=',I5,4X,'NSTAT=',I5,4X,'EN   =',F20.10)
 103  FORMAT(/' MGVN =',I2,4x,'s =',f6.1,4x,'sz =',f6.1,4x,'NELT =',I5,4x,'NNUC =',I3)
 104  FORMAT(/' NTGSYM=',i4/'  I   MCONT  NOTGT'/(i3,2I7))
 102  FORMAT(/' EIGEN-ENERGIES',/(16X,5F20.10))
 120  FORMAT(/' Nuclear data     X         Y         Z       Charge'/(3x,a8,2x,4F10.6))
!
      END SUBROUTINE WRITCIP

end module phase_match_mod
!
program phase_match
   use mpi_gbl, only: mpi_mod_start, mpi_mod_finalize
   use precisn_gbl, only: wp
   use cdenprop_io
   use cdenprop_defs
   use scatci_data, only : ntgtmx
   use phase_match_mod
   implicit none

   integer, parameter :: max_geom = 800

   type(CSFheader) :: csf_head_ref,csf_head_match
   type(CSFbody)   :: csf_body_ref,csf_body_match
   type(CIvect) :: ci_vec_ref,ci_vec_match
   type(CIvect), allocatable :: ci_vec_temp(:)

   integer :: lue
   integer :: i,a,b,geom,nstat,nset,total_sets,err,set
   integer :: notgs(ntgtmx)
   integer, allocatable :: ndtr_ref(:,:), ndtr_match(:,:), orb_data_ref(:,:), orb_data_match(:,:)
   real(kind=wp), allocatable :: spin_orb_olaps(:,:), wfs_olap(:,:), dgem(:), CV(:,:), ei(:)
   logical :: opn

   !Namelist variables
   integer :: lucivec(max_geom),nciset(max_geom), lucsf(max_geom),n_geom
   logical :: replace_with_phase_corrected_vectors
   character(len=132) :: moints(max_geom), energies_file
   character(len=2) :: str_set
   namelist/INPUT/lucivec,nciset,lucsf,n_geom,moints,replace_with_phase_corrected_vectors,energies_file

      call mpi_mod_start

      write(*,'(/,10X,"PHASE-MATCHING OF TARGET WAVEFUNCTIONS")')
      write(*,  '(10X,"======================================")')

      lucivec = 0
      nciset = 0
      lucsf = 0
      n_geom = 0
      moints = ''
      replace_with_phase_corrected_vectors = .false.
      energies_file = 'matched_energies'

      read(5,nml=input)

      write(*,'(/,5X,"Input data")')
      write(*,'(  5X,"----------")')

      write(*,'(/,5X,"Number of geometries: ",i5)') n_geom

      if (replace_with_phase_corrected_vectors) then
         write(*,'(5X,"The phase corrected target vectors will replace the corresponding CI vectors on the input file.")')
      else
         write(*,'(5X,"The phase corrected target vectors will not be written to a file.")')
      endif

      if (n_geom .le. 0) stop "Error: n_geom .le. 0"

      do i=1,n_geom
         write(*,'(/,5X,"Geometry: ",i5)') i
         write(*,'(5X,"CI vectors to match on file unit and set number: ",i10,1X,i5)') lucivec(i), nciset(i)
         write(*,'(5X,"CSFs on file unit: ",i10)') lucsf(i)
         write(*,'(5X,"Integrals file: ",a)') trim(adjustl(moints(i)))
         if (lucivec(i) .le. 0 .or. nciset(i) .le. 0) stop "Input error: at least one of lucivec(i),nciset(i) is .le. 0"
         if (lucsf(i) .le. 0) stop "Input error: lucsf.le. 0"
      enddo

      write(*,'(/,5X,"End of input data")')
      write(*,'(  5X,"----------")')

      write(*,'(/,5X,"Reading CSF file...")')

      call cwbopn(lucsf(1))
      call read_csf_head(lucsf(1),csf_head_ref,6)
      call read_csf_body(lucsf(1),csf_head_ref,csf_body_ref,6)
      close(lucsf(1))
      write(*,'(5X,"CSFs read in")')

      write(*,'(/,5X,"Number of configurations: ",i15)') csf_head_ref%nocsf

      if (csf_head_ref%ntgsym > 0) then
         stop "CSF file contains target*continuum expansion: not implemented yet"
      endif

      call expand_determinants(csf_head_ref,csf_body_ref,ndtr_ref)

      call map_spin_orbitals(csf_head_ref,orb_data_ref)

      allocate(dgem(csf_head_ref%nocsf))
      notgs = 0 !Only valid for target runs
      dgem = 0.0_wp !Diagonal elements of the Hamiltonian are skipped during the CI vector read in read_ci_vector

      write(*,'(/,5X,"Reading the starting CI vectors...")')

      call cwbopn(lucivec(1))
      call read_ci_vector(lucivec(1),nciset(1),ci_vec_ref,0,6)
      close(lucivec(1))
      nstat = size(ci_vec_ref%CV,2)

      write(*,'(5X,"Number of states: ",i15)') nstat

      if (csf_head_ref%nocsf .ne. size(ci_vec_ref%CV,1)) then
         print *,'number of CSFs on the CSF file and on the CI file dont match',csf_head_ref%nocsf,size(ci_vec_ref%CV,1)
         stop
      endif
      if (csf_head_ref%mgvn .ne. ci_vec_ref%mgvn) then
         print *,'spatial symmetries of the CSF wfs and those from the CI file dont match',ci_vec_ref%mgvn,csf_head_ref%mgvn
         stop
      endif
      if (csf_head_ref%S .ne. ci_vec_ref%S) then
         print *,'spin symmetries of the CSF wfs and those from the CI file dont match',csf_head_ref%S,ci_vec_ref%S
         stop
      endif

      !Write the starting set of target energies
      open(newunit=lue,file=energies_file,status='replace',form='formatted')
      geom = 1
      write(lue,'(i5,100e25.15)') geom, (ci_vec_ref%ei(i)+ci_vec_ref%e0,i=1,nstat)

      do geom=2,n_geom

         write(*,'(/,5X,"Geometry number: ",i4)') geom

         write(*,'(/,5X,"Reading CSF file...")')
         call cwbopn(lucsf(geom))
         call read_csf_head(lucsf(geom),csf_head_match,6)
         call read_csf_body(lucsf(geom),csf_head_match,csf_body_match,6)
         close(lucsf(geom))
         write(*,'(5X,"CSFs read in")')
   
         write(*,'(/,5X,"Number of configurations: ",i15)') csf_head_match%nocsf
   
         if (csf_head_match%ntgsym > 0) then
            stop "CSF file contains target*continuum expansion: not implemented yet"
         endif
   
         call expand_determinants(csf_head_match,csf_body_match,ndtr_match)
   
         call map_spin_orbitals(csf_head_match,orb_data_match)
  
         if (allocated(dgem)) deallocate(dgem)
         allocate(dgem(csf_head_match%nocsf))
         notgs = 0 !Only valid for target runs
         dgem = 0.0_wp !Diagonal elements of the Hamiltonian are skipped during the CI vector read in read_ci_vector

         !Read CI vectors for the geometry to match.
         write(*,'(/,5X,"Reading CI vectors file...")')
         call cwbopn(lucivec(geom))
         call read_ci_vector(lucivec(geom),nciset(geom),ci_vec_match,0,6)
         close(lucivec(geom))
         write(*,'(5X,"CI vectors read in")')

         !Check compatibility with the CI vector
         if (ci_vec_match%mgvn .ne. ci_vec_match%mgvn) then
            print *,'spatial symmetries of the wfs dont match',ci_vec_match%mgvn,ci_vec_match%mgvn
            stop
         endif
         if (ci_vec_match%S .ne. ci_vec_match%S) then
            print *,'spin symmetries of the wfs dont match',ci_vec_match%S,ci_vec_match%S
            stop
         endif
         if (ci_vec_match%nocsf .ne. ci_vec_match%nocsf) then
            print *,'number of CSFs of the wfs dont match',ci_vec_match%nocsf,ci_vec_match%nocsf
            stop
         endif
         if (size(ci_vec_match%CV,2) .ne. nstat) then
            print *,'number of vectors doesnt match',size(ci_vec_match%CV,2),nstat
            stop
         endif

         call load_spin_orbital_overlaps(orb_data_ref,orb_data_match,spin_orb_olaps,moints(geom-1),moints(geom))

         call calculate_wavefunction_overlaps (ndtr_ref, csf_head_ref, csf_body_ref, ndtr_match, csf_head_match, csf_body_match, &
                                               ci_vec_ref, ci_vec_match, spin_orb_olaps, wfs_olap, nstat)
   
         write(*,'(/,5X,"Overlap matrix S(b,a) of target wavefunctions: a = reference WF, b = WF for the current geometry")')
         do a=1,nstat
            do b=1,nstat
               write(*,'(5X,"S(",i4,",",i4,") = ",e25.15)') b,a,wfs_olap(b,a)
            enddo
         enddo

         !Match the wfs and phase-correct them. On output: CV,ei are the phase-corrected wfs and energies in the energy order given by the 1st geometry,
         !ci_vec_match are the phase-corrected wfs and energies in the energy order for the present geometry.
         call match_wavefunctions(wfs_olap,ci_vec_match,CV,ei)

         !Write the matched target energies for each geometry
         write(lue,'(i5,100e25.15)') geom, (ei(i)+ci_vec_match%e0,i=1,nstat)

         if (geom < n_geom) then
            !Prepare for the next geometry: make the reference CI vectors the vectors that we've just matched; the same for CSFs
            !ci_vec_ref = ci_vec_match
            ci_vec_ref%CV = CV
            ci_vec_ref%ei = ei
            csf_head_ref = csf_head_match
            csf_body_ref = csf_body_match
            deallocate(orb_data_ref,ndtr_ref)
            call move_alloc(orb_data_match,orb_data_ref)
            call move_alloc(ndtr_match,ndtr_ref)
            write(*,'(/,5X,"The CSFs and the phase-corrected vectors have been moved to the reference set.")')
         endif

         if (replace_with_phase_corrected_vectors) then
            write(str_set,'(i2)') nciset(geom)
            ci_vec_match%name = 'Phase-matched target vectors for the original SET '//str_set
   
            !nset = 0 !append to the CI vectors set for the current geometry
            nset = nciset(geom) !replace the CI vectors in the same set
            write(*,'(/,5X,"Output of the diagonal Hamiltonian elements to the CI vector set will be omitted.")')
            write(*,'(5X,"The following SET will be replaced with its phase-matched CI vectors: ",i4,/)') nciset(geom)

            !Rewrite the set that has been phase-matched and write again the sets that have not been touched.
            total_sets = 0
            CALL MOVEP(lucivec(geom),total_sets,err,0,6)
            total_sets = total_sets-1 !index of the last set stored on the file

            allocate(ci_vec_temp(total_sets),stat=err)
            if (err .ne. 0) then
               stop "error allocating ci_vec_temp"
            endif

            do i=1,total_sets
               set = i
               call read_ci_vector(lucivec(geom),set,ci_vec_temp(i),0,6)
            enddo

            do i=1,total_sets
               set = i
               if (set .eq. nset) then
                  call writcip (lucivec(geom), set, ci_vec_match % name, ci_vec_match % e0, ci_vec_match % nocsf, &
                                ci_vec_match % nstat, ci_vec_match % mgvn, ci_vec_match % s, ci_vec_match % sz ,&
                                ci_vec_match % nelt, ci_vec_match % ei, ci_vec_match % CV, ci_vec_match % iphz, dgem, &
                                csf_head_match % ntgsym, csf_head_match % mcont, notgs, ci_vec_match % nnuc, ci_vec_match % cname, &
                                ci_vec_match % xnuc, ci_vec_match % ynuc, ci_vec_match % znuc, ci_vec_match % charge, 6, 0, 0)
               else
                  call writcip (lucivec(geom), set, ci_vec_temp(i) % name, ci_vec_temp(i) % e0, ci_vec_temp(i) % nocsf, &
                                ci_vec_temp(i) % nstat, ci_vec_temp(i) % mgvn, ci_vec_temp(i) % s, ci_vec_temp(i) % sz, &
                                ci_vec_temp(i) % nelt, ci_vec_temp(i) % ei, ci_vec_temp(i) % CV, ci_vec_temp(i) % iphz, dgem, &
                                csf_head_match % ntgsym, csf_head_match % mcont, notgs, ci_vec_temp(i) % nnuc, &
                                ci_vec_temp(i) % cname, ci_vec_temp(i) % xnuc, ci_vec_temp(i) % ynuc, ci_vec_temp(i) % znuc, &
                                ci_vec_temp(i) % charge, 6, 0, 0)
               endif
            enddo

            close(lucivec(geom))
            deallocate(ci_vec_temp)
         endif

      enddo !geom

      close(lue)

      write(*,'(/,10X,"PHASE-MATCHING OF TARGET WAVEFUNCTIONS HAS FINISHED")')
      write(*,  '(10X,"===================================================")')

      call mpi_mod_finalize

end program phase_match
