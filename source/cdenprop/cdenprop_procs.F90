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

module cdenprop_procs

      use cdenprop_defs,   only: idp, tsym_max
      use iso_fortran_env, only: output_unit

      implicit none

      integer :: iwrite = output_unit

!     constants
      integer, parameter, private :: BraState = 1
      integer, parameter, private :: KetState = 2

!     namelist variables
      integer :: lucivec(2), lutdip, ntset, luneutralci, ludyson, ludyson_gcsf, lucivec_l2contract, lutarg_phase, &
     &           min_multipole, max_multipole, iuse_bound, ibra_or_ket_bound, lutargci,luprop,lupropw,isw,nbset,  &
     &           gflag, lucont_orbs, verbosity
      character(len=11) :: bform
      real(kind=idp) :: threshold, sparsity
      integer, dimension(tsym_max) :: lucsf, nciset, nstat, numtgt, notgt, lunp1targ, no_l2_virtuals, nciset_l2contract, &
                                      lutransmat_l2contract, itarget_spin
      logical :: ukrmolp_ints,target_phase_correction

!     Namelists
      namelist /DENINP/ nciset, lucsf, lucivec, luneutralci, numtgt, luprop, lupropw, lutargci,    &
     &                  nstat, lutdip, ntset, threshold, lunp1targ, no_l2_virtuals, lucont_orbs,   &
     &                  ludyson, ludyson_gcsf, itarget_spin, lutarg_phase, max_multipole,          &
     &                  iuse_bound, ibra_or_ket_bound,ukrmolp_ints, isw, target_phase_correction,  &
     &                  sparsity, nbset, bform, gflag, verbosity

      namelist /L2CONTRACT/ nciset_l2contract, lucivec_l2contract,lutransmat_l2contract

contains

      !> The main driver routine for CDENPROP
      subroutine cdenprop_drv(ci_vec_i,ci_vec_j,cden_namelists,mol_prop_data)
      use maths, only: maths_matrix_multiply_blas95
      use cdenprop_defs
      use cdenprop_io
      use cdenprop_aux
      use class_COOSparseMatrix_integer
      use class_COOSparseMatrix_real
      use ukrmol_interface_gbl, only: write_molden_dyson_orbitals, molecular_orbital_basis, destroy_ukrmolp
      use const_gbl, only: sym_op_nam_len, level1, level2, set_verbosity_level
      use mpi_gbl, only: mpi_mod_start, mpi_mod_finalize, mpi_started, mpi_running, myrank, master, mpi_xermsg
      use class_namelists
      use class_molecular_properties_data

      use class_namelists
      use class_molecular_properties_data

      implicit none
!     Arguments
      type(cdenprop_namelists), optional :: cden_namelists
      type(molecular_properties_data), optional, target :: mol_prop_data
      type(molecular_properties_data), target :: mol_prop_data_loc
      type(molecular_properties_data), pointer :: mol_prop_data_ptr
      type (CIvect), intent(inout) :: ci_vec_i,ci_vec_j !The N+1 CI vectors
!     Local
      integer :: i, j, k, l, m, icontains_continuum_orbitals, istate, jstate, icount, cdeninp, nspinsym, verbosity
      integer, allocatable, dimension(:) :: itarget_overall_phase_i, itarget_overall_phase_j, id_map_targ_to_congen, idtarg, &
                                            full_idtarg, idtarg_test, itarget_irrep
      integer, allocatable, dimension(:,:) :: itarget_symmetry_order,icontinuum_orbital_table_i,icontinuum_orbital_table_j
      real(kind=idp), allocatable, dimension(:):: properties_cc,  properties_cl2, properties_l2c, properties_l2l2, target_spin
      REAL(KIND=idp), allocatable, dimension(:,:) ::  density_matrix_coefficients_cl2
      INTEGER, allocatable, DIMENSION(:,:,:) ::  density_matrix_orbital_pairs_cl2
      real(kind=idp), allocatable, DIMENSION(:,:,:) ::  dyson_orbitals_cl2, dyson_orbitals_l2c, properties_for_write, &
                                                        single_property(:,:)

!     Sparse Matrixes for the L2-L2 block
      type(COOMatrix_real) :: density_matrix_coefficients_l2l2
      type(COOMatrix_integer) :: density_matrix_orbital_pairs_l2l2_1,density_matrix_orbital_pairs_l2l2_2

!     The target CI vector
      integer :: itarg, iphase_np1_targ_ci, iphase_targ_ci, iskip,all_ntgsym
      type (CIvect), allocatable, dimension(:) :: ci_vec_target, all_ci_vec_target
      type (CIvect) :: ci_vec_target_temp
      real(kind=idp) :: start, finish,spin_temp
      real(kind=idp), allocatable, dimension(:) :: np1_calc_target_vectors
      integer, allocatable ::  all_numtgt(:)

!     The N+1 CI vectors
      type (CIvect) :: ci_vec_contracted_i,ci_vec_contracted_j,states_from_bound
      real(kind=idp), allocatable, dimension(:,:) :: bound_vec_i, bound_vec_j
      character(len=1) :: TRANSI, TRANSJ, TRANS

!     The N+1 CSF file.
      type (CSFheader) :: csf_head_i,csf_head_j
      type (CSFbody)   :: csf_body_i,  csf_body_j

!     The propert integral table
      type (property_integrals), pointer :: pintegrals => null()
      integer, allocatable :: lm_to_l_m(:,:)

!     For reading target properties
      integer :: no_target_states,  ismax, ifail, nnuc
      real(kind=idp), allocatable, dimension(:,:,:) ::target_properties
      integer, allocatable, dimension(:) ::prop_order

!     Variables for printing properties table
      integer :: non_zero_properties
      character(len=sym_op_nam_len) :: nam
      logical :: skip

!     Variables for new matrix scheme
      type (CIvect) :: properties_times_civec_j,prototype_properties_matrix,prototype_properties_matrix_temp,properties_temp
      type (CIvect) :: temp_ci_vector_i,temp_ci_vector_j
      integer :: inumber_of_gcsf_j,inumber_of_gcsf_i

!     Dyson orbitals
      real(kind=idp) :: dyson_orbital_magnitude
      REAL(KIND=idp), allocatable, dimension(:,:,:) :: dyson_orbitals
      real(kind=idp), allocatable, dimension(:,:) :: dyson_orbital_norms
      integer, allocatable, dimension(:) :: idyson_orbital_irrep,idyson_itarget_irrep,idyson_itarget_spin
      integer,allocatable,dimension(:,:)::orbs_diag_l2l2
      logical :: scattering_calculation, mpi_i_write

!     Variables for the L2 contraction
      integer :: noriginal_l2_basis_i, ncontracted_l2_basis_i, noriginal_l2_basis_j, ncontracted_l2_basis_j
      integer :: ifirst_l2_csf_i, ifirst_l2_csf_j
      REAL(KIND=idp), allocatable, dimension(:,:) ::l2block_contract_i, l2block_contract_j, &
                                                    l2block_contract_transpose_i, l2block_contract_transpose_j

!     SCALAPACK variables
      logical :: use_scalapack_mat_ops
      integer :: icontext

!     Namelist default values
!---------------------------------------------------------------------

!     Defaults for DENINP
      isw=0
      luprop=17                ! File unit for input of property integrals.
      ludyson=123              ! File unit for dyson orbital output.
      ludyson_gcsf=9978        ! Output information needed to contract the L2 configurations.
      lupropw=667              ! File unit for dipole output.
      lutargci=26              ! Target CI vectors from target run
      lucont_orbs=222          ! File unit for continuum orbital table.
      lutarg_phase=333         ! File unit for scattering run target vectors.
      nstat=0                  ! Number of states to calculate dipoles for defaults to all (0,0).
      numtgt=0                 ! Number of TARGET states per symmetry (for reconstruction of bound+continuum configurations).
      lutdip=24                ! Fileunit for bound state transition dipoles.
      gflag=1                  ! flag for CSFs that are grouped by symmetry in congen (1=on [default])
      ntset=1                  ! Set for bound state transtion dipoles.
      threshold =1.E-10_idp    ! Printing threshold for dipoles.
      no_l2_virtuals=0         ! Don't change.
      itarget_spin=2           ! Array of dim no_target_symmetries: defaults to doublet.
      min_multipole=1          ! Omit l = 0 properties by default.
      max_multipole=1          ! 1: Dipole transitions to the continuum.
                               ! 2: Quad transitions (Not implemented).
      iuse_bound=0             ! Use bound states from BOUND.
      ibra_or_ket_bound=0      ! Set which are the bound states - only needed for Dyson orbital
                               ! output when iuse_bound=0 and both bra and ket contain continuum
                               ! orbitals.
      ukrmolp_ints = .true.    ! Are we using property integrals generated using the new integral code?
      !ZM: switch to enable/disable matching of phases of the target
      !vectors from the target run with those of the N+1 run. This
      !should be set to .false. ONLY IF the target vectors used in N+1
      !SCATCI run were forced from the target run (using the NFTG flag).
      target_phase_correction = .true.
      !
      !ZM: the units for both I and J vectors must be specified. In case
      !of a pure target run lucivec(1) is used for the target vectors.
      !Taking the CI vector data directly from the argument ci_vec_i is specified setting lucivec(1) to a value .le. 0
      !Taking the CI vector data directly from the argument ci_vec_j is specified setting lucivec(2) to a value .le. 0
      lucivec(1) = -1 !unit for the I states
      lucivec(2) = -1 !unit for the J states
      !
      !ZM: unit number for output of BOUND
      luneutralci = -1
      states_from_bound%CV_is_scalapack = .false. !so far we only allow bound states from the disk without distribution in SCALAPACK format.
      nbset = 1
      bform = 'UNFORMATTED'
      !
      lunp1targ(1)=2001
      lunp1targ(2)=2002
      !
      !ZM: sparsity of the matrix. Default is 1%
      sparsity =  0.01d0
      verbosity = 1

!     Defaults for L2CONTRACT
      nciset_l2contract=0
      lucivec_l2contract=0
      lutransmat_l2contract(1)=3001
      lutransmat_l2contract(2)=3002

!     Read input file
      if (.not. present(cden_namelists) ) then
         if (mpi_running) then !each MPI task reads the namelist input from the file den.inp
            open(newunit=cdeninp,file='./den.inp',form='formatted',status='unknown')
            read(cdeninp,nml=DENINP)
            read(cdeninp,nml=L2CONTRACT)
            close(cdeninp)
         else
            read(*,DENINP)
            read(*,L2CONTRACT)
         endif
         scattering_calculation=.true.
         call set_verbosity_level(verbosity)
      else
         call cden_namelists % return_namelist_variables (lucsf, lucivec(1), nciset, nstat, numtgt, lutdip, luprop, &
                                                          lupropw, max_multipole, ukrmolp_ints, isw, gflag, itarget_spin)
         lucivec(2)=lucivec(1)
         scattering_calculation=.false.
         target_phase_correction=.false.
      end if

      if (ukrmolp_ints) then
         if (.not. mpi_started) call mpi_mod_start(.true.)
         iwrite = level1
         mpi_i_write = .false.
         if (myrank .eq. master) mpi_i_write = .true. !in parallel calculation only master writes
      else
         mpi_i_write = .true. !serial calculation
      endif

      call density_matrix_coefficients_l2l2%set_iwrite(level2)
      call density_matrix_orbital_pairs_l2l2_1%set_iwrite(level2)
      call density_matrix_orbital_pairs_l2l2_2%set_iwrite(level2)

      ! number of target spin-symmetries with non-zero number of target states
      nspinsym = count(numtgt > 0)

      write(iwrite,'(/, " INPUT FILES")')
      write(iwrite,'(" -----------")')
      write(iwrite,'(" Configuration state functions to be read from files",5x,30i5)') (lucsf(i),i=1,2) ! change for more than 2 syms
      if (lucivec(1) .le. 0) then
         write(iwrite,'(" CI vectors for I states to be taken from ci_vec_i  ",5x)')
         if (ci_vec_i%CV_is_scalapack)  write(iwrite,'(" CI vectors ci_vec_i are in SCALAPACK format        ",5x)')
      else
         write(iwrite,'(" CI vectors for I states to be read from file       ",5x,i5)') lucivec(1)
      endif
      if (lucivec(2) .le. 0) then
         write(iwrite,'(" CI vectors for J states to be taken from ci_vec_j  ",5x)')
         if (ci_vec_j%CV_is_scalapack)  write(iwrite,'(" CI vectors ci_vec_j are in SCALAPACK format        ",5x)')
      else
         write(iwrite,'(" CI vectors for J states to be read from file       ",5x,i5)') lucivec(2)
      endif
      write(iwrite,'(" CI vectors to be read from sets                    ",5x,30i5)') (nciset(i),i=1,2) !change for more than 2 syms
      write(iwrite,'(" Number of states                                   ",5x,2i5)') (nstat(i),i=1,2)
      write(iwrite,'(" Group configurations by states (gflag)             ",5x,i5)') gflag
      write(iwrite,'(" Use UKRmol+ integrals                              ",5x,l1)') ukrmolp_ints
      write(iwrite,'(" Verbosity level                                    ",5x,i5)') verbosity
      if (gflag == 0) then
          write(iwrite,'(" Number of target states per spin-symmetry ")')
          write(iwrite,'(10i5)') numtgt(1:nspinsym)
          write(iwrite,'(" Target spin per spin-symmetry ")')
          write(iwrite,'(10i5)') itarget_spin(1:nspinsym)
      else
          write(iwrite,'(" Number of target states per spin-symmetry ",*(i5))') numtgt(1:nspinsym)
          write(iwrite,'(" Target spin per spin-symmetry             ",*(i5))') itarget_spin(1:nspinsym)
      end if
      write(iwrite,'(" Which of the states is bound (ibra_or_ket_bound)   ",5x,i5)') ibra_or_ket_bound
      if (iuse_bound .ne. 0) then
         write(iwrite,'(" Bound states to be read from file                  ",5x,30i5)') luneutralci
         write(iwrite,'(" Dataset to read bound state from                   ",5x,i5)') nbset
         write(iwrite,'(" Format of bound state file                         ",5x,i5)') bform
      end if
      if (ukrmolp_ints) then
         write(iwrite,'(" Assuming UKRmol+ property integrals on input")')
         if (.not. mpi_started) call mpi_mod_start
      else
         write(iwrite,'(" Assuming SWEDEN property integrals on input")')
      endif
      write(iwrite,'(" Target phase correction                            ",5X,l1)') target_phase_correction
      write(iwrite,'(" Property integrals to be read from                 ",5x,30i5)') luprop
      write(iwrite,'(" Target properties to be read from                  ",5x,30i5)') lutdip
      write(iwrite,'(" Dataset to read target properties from             ",5x,i5)') ntset
      write(iwrite,'(" Target run: target vectors                         ",5x,30i5)') lutargci
      if (target_phase_correction) then
         write(iwrite,'(" Inner region run: target vectors                   ",5x,30i5)') (lunp1targ(i),i=1,2) !change for more than 2 syms
      endif
      write(iwrite,'(" L2 transformation matrices                         ",5x,30i5)') (lutransmat_l2contract(i),i=1,2)
      write(iwrite,'(" Sparsity of the density matrices                   ",5x,f10.6,"%")') sparsity*100
      write(iwrite,'(" Max multipole to calculate (max_multipole)         ",5x,i5)') max_multipole
      write(iwrite,'(" Nuclear multipole flag (isw)                       ",5x,i5)') isw
      write(iwrite,'(" Property drop threshold (threshold)                ",5x,e25.15)') threshold

      write(iwrite,'(/, " OUTPUT FILES")')
      write(iwrite,'(   " ------------")')
      write(iwrite,'(" Inner region properties                            ",5x,30i5)') lupropw
      write(iwrite,'(" Data for Lsquared contraction                      ",5x,30i5)') ludyson_gcsf
      write(iwrite,'(" The Dyson orbitals                                 ",5x,30i5)') ludyson
      write(iwrite,'(" Continuum orbital table                            ",5x,30i5)') lucont_orbs
      write(iwrite,'(" Target phases                                      ",5x,30i5)') lutarg_phase
      write(iwrite,'("")')
!
!     ZM error check on lucivec(1) which is always used (the check on lucivec(2) is done below if it is used).
      if (lucivec(1) .le. 0) then
         if (allocated(ci_vec_i%CV)) then
            write(iwrite,'(/," Dimensions of the CI vector matrix in ci_vec_i: ",2i10)') size(ci_vec_i%CV,1), size(ci_vec_i%CV,2)
         else
            write(iwrite,'(" On input the matrix of CI vectors from the argument ci_vec_i was not allocated")')
            stop
         endif
      endif

!     ZM: IMPORTANT - we only allow to use ci_vec_i, ci_vec_j if they are both in the same format.
!
      if (ci_vec_i % CV_is_scalapack .neqv. ci_vec_j % CV_is_scalapack) then
         print *,'Only one of ci_vec_i,ci_vec_j is in the SCALAPACK format. This situation is not supported yet.'
         stop
      endif

!     ZM: if the ci_vec_i, ci_vec_j CI vectors are in SCALAPACK format then all the other arrays which
!     are coupled to them must use the SCALAPACK format too.
!
      icontext = 0
      if (ci_vec_i%CV_is_scalapack .and. ci_vec_j%CV_is_scalapack) then
         use_scalapack_mat_ops = .true.
         !
         if (ci_vec_i%CV_is_scalapack) icontext = ci_vec_i%blacs_context
         if (ci_vec_j%CV_is_scalapack) icontext = ci_vec_j%blacs_context
         !
         properties_times_civec_j%CV_is_scalapack = .true.
         properties_times_civec_j%blacs_context = icontext

         prototype_properties_matrix%CV_is_scalapack = .true.
         prototype_properties_matrix%blacs_context = icontext

         prototype_properties_matrix_temp%CV_is_scalapack = .true.
         prototype_properties_matrix_temp%blacs_context = icontext

         properties_temp%CV_is_scalapack = .true.
         properties_temp%blacs_context = icontext

         temp_ci_vector_i%CV_is_scalapack = .true.
         temp_ci_vector_i%blacs_context = icontext

         temp_ci_vector_j%CV_is_scalapack = .true.
         temp_ci_vector_j%blacs_context = icontext
      else
         use_scalapack_mat_ops = .false.
         !
         properties_times_civec_j%CV_is_scalapack = .false.
         prototype_properties_matrix%CV_is_scalapack = .false.
         prototype_properties_matrix_temp%CV_is_scalapack = .false.
         properties_temp%CV_is_scalapack = .false.
         temp_ci_vector_i%CV_is_scalapack = .false.
         temp_ci_vector_j%CV_is_scalapack = .false.
      endif

      if (use_scalapack_mat_ops .and. (lucivec_l2contract .ne. 0)) then
         print *,'L2 contraction option with SCALAPACK matrices not supported'
         stop
      endif
!
!     Read in the CSF files for the two symmetries
!     ------------------------------------------------------------------
      write(iwrite,'(/, " READ CSF FILES")')
      write(iwrite,'(   " --------------")')

!     Read  CSF file I
      call cwbopn(lucsf(1))
      call read_csf_head(lucsf(1),csf_head_i,iwrite)
      call read_csf_body(lucsf(1),csf_head_i,csf_body_i,iwrite)

!     Continuum orbitals are involved we need to figure out some
!     details of the contraction: The target state symmetries
!     the number of target states per symmetry, the number
!     of continuum orbitals per symmetry.

      if (csf_head_i%ntgsym .ne. -1) then
         csf_head_i%no_l2_virtuals(1:csf_head_i%nsym) = no_l2_virtuals(1:csf_head_i%nsym)
         notgt(1 : csf_head_i % nsym) = csf_head_i % nob(1 : csf_head_i % nsym) &
                                      - csf_head_i % nob0(1 : csf_head_i % nsym) &
                                      - csf_head_i % no_l2_virtuals(1 : csf_head_i % nsym)
         allocate(csf_head_i%numtgt(csf_head_i%ntgsym))

         do i=1,csf_head_i%ntgsym
            csf_head_i%notgt(i)=notgt(csf_head_i%mcont(i)+1)

         end do

         csf_head_i%numtgt=numtgt(1:csf_head_i%ntgsym)

!        Calculate target state symmetries from the continuum orbitals
!        symmetry

         allocate( itarget_symmetry_order(csf_head_i%ntgsym,2))
         itarget_symmetry_order=0

         do i=1, csf_head_i%ntgsym
            itarget_symmetry_order(i,1)=IPD2H(csf_head_i%mgvn+1,csf_head_i%mcont(i)+1) !Target irrep
            itarget_symmetry_order(i,2)=itarget_spin(i)     !Target spin

         end do

         allocate(csf_head_i%itarget_symmetry_order(csf_head_i%ntgsym,2))
         csf_head_i%itarget_symmetry_order=itarget_symmetry_order

!-------

      end if

!     Read  CSF file J

      call cwbopn(lucsf(2))
      call read_csf_head(lucsf(2),csf_head_j,iwrite)
      call read_csf_body(lucsf(2),csf_head_j,csf_body_j,iwrite)


      if (csf_head_j%ntgsym .ne. -1) then
         csf_head_j%no_l2_virtuals(1:csf_head_i%nsym) = no_l2_virtuals(1:csf_head_i%nsym)
         notgt(1 : csf_head_j % nsym) = csf_head_j % nob(1 : csf_head_j % nsym) - &
                                        csf_head_j % nob0(1 : csf_head_j % nsym) - &
                                        csf_head_j % no_l2_virtuals(1 : csf_head_j % nsym)
         allocate(csf_head_j%numtgt(csf_head_j%ntgsym))

         do i=1,csf_head_j%ntgsym
            csf_head_j%notgt(i)=notgt(csf_head_j%mcont(i)+1)

         end do

         csf_head_j%numtgt=numtgt(1:csf_head_j%ntgsym)

!        Check that target state symmetries match in both CSF (FELIPE: ONLY IF I HAS CONTINUUM ORBITALS
         if (csf_head_i%ntgsym .ne. -1) then
            do i=1, csf_head_j%ntgsym
               if ( itarget_symmetry_order(i,1) .ne. IPD2H(csf_head_j%mgvn+1,csf_head_j%mcont(i)+1) ) then
                  print '(1x,a)', 'TARGET STATE SYMMETRIES IN I AND J DIFFER'
                  print '(4x,*(a,i0))', 'I: total ', csf_head_i%mgvn+1, ', continuum ', csf_head_i%mcont(i)+1, &
                                        ', deduced target ', itarget_symmetry_order(i,1)
                  print '(4x,*(a,i0))', 'J: total ', csf_head_j%mgvn+1, ', continuum ', csf_head_j%mcont(i)+1, &
                                        ', deduced target ', IPD2H(csf_head_j%mgvn+1,csf_head_j%mcont(i)+1)
                  stop 1
               end if

            end do

         else
            allocate( itarget_symmetry_order(csf_head_j%ntgsym,2))
            do i=1, csf_head_j%ntgsym
               itarget_symmetry_order(i,1)=IPD2H(csf_head_j%mgvn+1,csf_head_j%mcont(i)+1)
               itarget_symmetry_order(i,2)=itarget_spin(i)

            end do

         end if

         allocate(csf_head_j%itarget_symmetry_order(csf_head_j%ntgsym,2))
         csf_head_j%itarget_symmetry_order=itarget_symmetry_order

      end if

!     Sanity check orbital set is the same for both symmetries and
!     determine which if any of the input CSF files contain continuum
!     orbitals (ntgsym=-1 implies target run).

      if ((csf_head_i%ntgsym .ne. -1) .and. (csf_head_j%ntgsym .ne. -1))  then
         if ((sum(csf_head_i%nob).ne.sum(csf_head_j%nob)).or.(csf_head_i%nsym.ne.csf_head_j%nsym)) then
            stop ' PROBLEM WITH READING THE CSFS: MISMATCH IN NOB OR NSYM 1'
         endif
         write(iwrite,'(/," Both I and J contain continuum orbitals",/)')
         icontains_continuum_orbitals=3

      else if  ((csf_head_i%ntgsym .ne. -1) .and. (csf_head_j%ntgsym .eq. -1))  then
         write(iwrite,'(/," Only I contains continuum orbitals",/)')
         icontains_continuum_orbitals=1

      else if  ((csf_head_i%ntgsym .eq. -1) .and. (csf_head_j%ntgsym .ne. -1))  then
         write(iwrite,'(/," Only J contains continuum orbitals",/)')
         icontains_continuum_orbitals=2

      else if  ((csf_head_i%ntgsym .eq. -1) .and. (csf_head_j%ntgsym .eq. -1))  then
         write(iwrite,'(/," No continuum orbitals. This is a target run.",/)')
         icontains_continuum_orbitals=0

     end if

      if (present(mol_prop_data)) then
         mol_prop_data % contains_continuum = icontains_continuum_orbitals > 0
      end if

      if (icontains_continuum_orbitals > 0 .and. isw == 0) then
         call mpi_xermsg('cdenprop_procs', 'cdenprop_drv', &
                         'You are calculating continuum properties and including nuclear multipole (isw = 0).' // new_line(' ') // &
                         ' *  Nuclear multipole should be off for photoionization, both at DENPROP and CDENPROP stage!', 0, 0)
      end if

!     Determine the number of L^2 CSFs and the index of the last
!     CSF involving a continuum orbital.

!     For target runs all CSFs are L2

      csf_head_i%l2nocsf=csf_head_i%nocsf
      csf_head_j%l2nocsf=csf_head_j%nocsf
      csf_head_i%last_continuum_csf=0
      csf_head_j%last_continuum_csf=0

      select case(icontains_continuum_orbitals)
      case(0)
!        Target run: nothing done
         ibra_or_ket_bound=0
      case(1)
!        Only I contains continuum orbitals
         csf_head_i%l2nocsf=csf_head_i%nocsf-2*csf_head_i%nctarg
         csf_head_i%last_continuum_csf=csf_head_i%nocsf-csf_head_i%l2nocsf

!        J is the bound state
         ibra_or_ket_bound=2
      case(2)
!         Only J contains continuum orbitals
         csf_head_j%l2nocsf=csf_head_j%nocsf-2*csf_head_j%nctarg
         csf_head_j%last_continuum_csf=csf_head_j%nocsf-csf_head_j%l2nocsf

!        I is the bound state
         ibra_or_ket_bound=1

      case(3)
!         Both I and J contain continuum orbitals
         csf_head_i%l2nocsf=csf_head_i%nocsf-2*csf_head_i%nctarg
         csf_head_i%last_continuum_csf=csf_head_i%nocsf-csf_head_i%l2nocsf

         csf_head_j%l2nocsf=csf_head_j%nocsf-2*csf_head_j%nctarg
         csf_head_j%last_continuum_csf=csf_head_j%nocsf-csf_head_j%l2nocsf

!        Either I or J could be the bound state and contains continuum orbitals
         select case(iuse_bound)
         case(1)
!           I is the bound state and contains continuum orbitals (BOUND used)
            ibra_or_ket_bound=1
            if (nstat(1) /= 0) then
               call mpi_xermsg('cdenprop_procs', 'cdenprop_drv', &
                               'Bra BOUND coefficients provided, but nstat(1) is not 0. Result may be incomplete!', 0, 0)
            end if
         case(2)
!           J is the bound state and contains continuum orbitals (BOUND used)
            ibra_or_ket_bound=2
            if (nstat(2) /= 0) then
               call mpi_xermsg('cdenprop_procs', 'cdenprop_drv', &
                               'Ket BOUND coefficients provided, but nstat(2) is not 0. Result may be incomplete!', 0, 0)
            end if
         case(0)
!           BOUND not used I or J could be the bound state
            select case(ibra_or_ket_bound)
            case(0)
!              We don't know or dont care which is the bound state
               write(iwrite,*) " Warning: ibra_or_ket_bound must be set if iuse_bound=0 &
                                &and both bra and ket contain continuum orbitals."
               write(iwrite,*) "          otherwise dyson orbitals and target phases won't be written"
            case(1)
!              I is the bound state and contains continuum orbitals
            case(2)
!              J is the bound state and contains continuum orbitals
            end select

         end select


      end select


!     Read in the CI vectors
!     ZM: the CI vectors for the I and J states are taken either from the
!     disk or directly from the cdenprop_drv argument structures ci_vec_i, ci_vec_j.
!     ------------------------------------------------------------------
      call read_I_J_states (icontains_continuum_orbitals, iuse_bound, lucivec, nciset, luneutralci, &
                            nbset, bform, nstat, ci_vec_i, ci_vec_j, states_from_bound)
!
!

!
!     ZM: error checking on compatibility of the symmetries of the CI vectors and CSFs.
!
      if (ci_vec_i%mgvn .ne. csf_head_i%mgvn .or. ci_vec_i%s .ne. csf_head_i%s .or. ci_vec_i%sz .ne. csf_head_i%sz) then
         print *,ci_vec_i%mgvn,csf_head_i%mgvn,ci_vec_i%s,csf_head_i%s,ci_vec_i%sz,csf_head_i%sz
         print *,"Symmetries of the CSFs and of the CI vectors for the I-states are not compatible"
         stop
      endif
      if (ci_vec_j%mgvn .ne. csf_head_j%mgvn .or. ci_vec_j%s .ne. csf_head_j%s .or. ci_vec_j%sz .ne. csf_head_j%sz) then
         print *,ci_vec_j%mgvn,csf_head_j%mgvn,ci_vec_j%s,csf_head_j%s,ci_vec_j%sz,csf_head_j%sz
         print *,"Symmetries of the CSFs and of the CI vectors for the J-states are not compatible"
         stop
      endif
      if (iuse_bound .eq. 1) then
         if (states_from_bound%mgvn .ne. csf_head_i%mgvn) then
            print *,iuse_bound,states_from_bound%mgvn,csf_head_i%mgvn
            print *,"Symmetry of the I CSFs and of the bound CI vector are not compatible"
            stop
         endif
      elseif (iuse_bound .eq. 2) then
         if (states_from_bound%mgvn .ne. csf_head_j%mgvn) then
            print *,iuse_bound,states_from_bound%mgvn,csf_head_j%mgvn
            print *,"Symmetry of the J CSFs and of the bound CI vector are not compatible"
            stop
         endif
      endif

      if (nstat(1) .eq. 0) then
         nstat(1) = ci_vec_i%nstat
      end if

      if (nstat(2) .eq. 0) then
         nstat(2) = ci_vec_j%nstat
      end if

      if ((nstat(1) .le. 0) .or. (nstat(2) .le. 0)) then
         print *,'nstat assignment error',nstat(1:2),ci_vec_i%nstat,ci_vec_j%nstat
         stop
      endif

      if (present(mol_prop_data)) then
         if (associated(mol_prop_data % pintegrals)) then
            pintegrals => mol_prop_data % pintegrals
         end if
      end if

      if (.not. associated(pintegrals)) then
!        Read in the property integral table
!        ------------------------------------------------------------------
         write(iwrite,'(/, " ALLOCATE PROPERTY INTEGRAL")')
         allocate(pintegrals)
         write(iwrite,'(/, " READ PROPERTY INTEGRAL")')
         write(iwrite,'(   " ----------------------")')
         call read_property_integrals2(luprop,pintegrals,iwrite,ukrmolp_ints)

!        Set up pointer array from symmetry block to position in the
!        integral table.
         call symmetry_block_to_integral_table(pintegrals%istart,csf_head_i%nob,csf_head_i%nsym)
      end if

!     Construct the density matrix
!     ------------------------------------------------------------------
      write(iwrite,'(/, " DENSITY MATRIX CONSTRUCTION")')
      write(iwrite,'(   " ----------------------------------")')

!     Blocks containing continuum configurations if necessary

      select case(icontains_continuum_orbitals)
      case(1)
! Read the N+1 CI vectors correctly for j

!        Only I contains continuum orbitals
!        Read in the target CI vectors

         call read_all_ci_vectors(lutargci,csf_head_i%ntgsym,numtgt,itarget_symmetry_order,&
     &                                all_ntgsym,all_numtgt,all_ci_vec_target,iwrite,gflag)

         call prop_order_ci_vec(prop_order, all_ntgsym, sum(all_numtgt), all_ci_vec_target)

         call congen_to_target_energy_order_map(all_ci_vec_target,csf_head_i%ntgsym,numtgt,all_ntgsym, all_numtgt,&
     &                                          idtarg,ifail)

         allocate(ci_vec_target(csf_head_i%ntgsym))
         do itarg=1,csf_head_i%ntgsym
            ci_vec_target(itarg)=all_ci_vec_target(itarg)
            call all_ci_vec_target(itarg) % dealloc

         end do

! DO FIRST: PHASES

!        figure out the phase
         no_target_states=sum(numtgt(1:csf_head_i%ntgsym)) !ZM: sum only over the defined elements
         allocate(csf_head_i%itarget_overall_phase(no_target_states),csf_head_j%itarget_overall_phase(no_target_states))
         csf_head_i%itarget_overall_phase=0
         csf_head_j%itarget_overall_phase=0
!ZM
         if (target_phase_correction) then

            !Determine overall target state phase introduced by the rediagonalisation of the target
            !Hamiltonian during the N+1 inner region Hamiltonian diagonalisation. Note, this is
            !done in scatci so as to deal with phase issues coming from 'dictionary' orbital ordering
            !and construction of spin orbitals as excitations from reference determinant between
            !the target and scattering CONGEN runs. This phase is also recorded in iphase, but
            !unfortunately an overall phase coming from the diagonalisation is not.
            icount=0
            itarg=0
            do i=1, csf_head_i%ntgsym
               do j=1,numtgt(i)
                  itarg=itarg+ci_vec_target(i)%nocsf
                  do k=1,csf_head_i%nctgt(i)
                     icount=icount+1

                  end do

               end do

            end do

            allocate(np1_calc_target_vectors(icount))
            np1_calc_target_vectors=0.0_idp
            write(iwrite,*) lunp1targ(1)
            call CWBOPN(lunp1targ(1))
            read(lunp1targ(1)) np1_calc_target_vectors
            close(lunp1targ(1))

            call determine_target_phase(csf_head_i%iphz,np1_calc_target_vectors,&
     &                ci_vec_target,csf_head_i%ntgsym,numtgt,&
     &                csf_head_i%itarget_overall_phase, ifail)
         else
            csf_head_i%itarget_overall_phase=1
         endif

         call phases_in_target_energy_order_map(ci_vec_target,csf_head_i%ntgsym,numtgt,&
     &             id_map_targ_to_congen,idtarg_test,itarget_irrep,target_spin,ifail)


         call write_target_phases(csf_head_i%ntgsym,csf_head_i%numtgt,csf_head_i%itarget_overall_phase)

         if (mpi_i_write) then
            write(lutarg_phase,*) 'Target phases'
            write(lutarg_phase,*) '-------------'
            write(lutarg_phase,*) 'State No.    I.Rep.     Spin       Phase'
            do i=1,no_target_states
               write(lutarg_phase,'(i7,8x,i2,8x,f4.1,8x,i2)') i,itarget_irrep(id_map_targ_to_congen(i)),&
                                                              target_spin(id_map_targ_to_congen(i)),  &
                                                              csf_head_i%itarget_overall_phase(id_map_targ_to_congen(i))
            end do
         end if

! DO FIRST: CREATE icontinuum_orbital_table_i (code inside continuum_continuum_drv)
         call read_transdip2 (iwrite, lutdip, ntset, no_target_states, numtgt(1 : csf_head_i % ntgsym), &
                              itarget_symmetry_order(1 : csf_head_i % ntgsym, :), ci_vec_i % nnuc, ismax, target_properties, &
                              gflag, prop_order, ifail)

         call continuum_continuum_only_one_drv(csf_head_i,csf_body_i, icontinuum_orbital_table_i)

      case(2)
!      Only J contains continuum orbitals
!        Read in the target CI vectors
         call read_all_ci_vectors(lutargci,csf_head_j%ntgsym,numtgt,itarget_symmetry_order,&
     &                                all_ntgsym,all_numtgt,all_ci_vec_target,iwrite,gflag)

         call prop_order_ci_vec(prop_order, all_ntgsym, sum(all_numtgt), all_ci_vec_target)

         call congen_to_target_energy_order_map(all_ci_vec_target,csf_head_j%ntgsym,numtgt,all_ntgsym, all_numtgt,&
     &                                          idtarg,ifail)

         allocate(ci_vec_target(csf_head_j%ntgsym))
         do itarg=1,csf_head_j%ntgsym
            ci_vec_target(itarg)=all_ci_vec_target(itarg)
            call all_ci_vec_target(itarg) % dealloc

         end do

         no_target_states=sum(numtgt(1:csf_head_j%ntgsym)) !ZM: sum only over the defined values
         allocate(csf_head_i%itarget_overall_phase(no_target_states),csf_head_j%itarget_overall_phase(no_target_states))
         csf_head_i%itarget_overall_phase=0
         csf_head_j%itarget_overall_phase=0
!ZM
         if (target_phase_correction) then

            icount=0
            itarg=0
            do i=1, csf_head_j%ntgsym
               do j=1,numtgt(i)
                  itarg=itarg+ci_vec_target(i)%nocsf
                  do k=1,csf_head_j%nctgt(i)
                     icount=icount+1

                  end do

               end do

            end do

            allocate(np1_calc_target_vectors(icount))
            np1_calc_target_vectors=0.0_idp
            write(iwrite,*) lunp1targ(2)
            call CWBOPN(lunp1targ(2))
            read(lunp1targ(2)) np1_calc_target_vectors
            close(lunp1targ(2))

            call determine_target_phase(csf_head_j%iphz,np1_calc_target_vectors,&
     &                ci_vec_target,csf_head_j%ntgsym,numtgt,&
     &                csf_head_j%itarget_overall_phase, ifail)
         else
            csf_head_j%itarget_overall_phase=1
         endif


         call phases_in_target_energy_order_map(ci_vec_target,csf_head_j%ntgsym,numtgt,&
     &             id_map_targ_to_congen,idtarg_test,itarget_irrep,target_spin,ifail)

         call write_target_phases(csf_head_j%ntgsym,csf_head_j%numtgt,csf_head_j%itarget_overall_phase)

         if (mpi_i_write) then
            write(lutarg_phase,*) 'Target phases'
            write(lutarg_phase,*) '-------------'
            write(lutarg_phase,*) 'State No.    I.Rep.     Spin       Phase'
            do i=1,no_target_states
               write(lutarg_phase,'(i7,8x,i2,8x,f4.1,8x,i2)') i,itarget_irrep(id_map_targ_to_congen(i)),&
                                                              target_spin(id_map_targ_to_congen(i)),  &
                                                              csf_head_j%itarget_overall_phase(id_map_targ_to_congen(i))
            end do
         end if

         if (csf_head_i%ntgsym .ne. csf_head_j%ntgsym) then
            !ZM in this case the I and J wavefunctions come from different N+1 calculations: the CDENPROP assumption is that both come from the same calculation.
            write(iwrite,'("Number of target states in I and J wavefunctions is not equal")')
            print *,"Error: this calculation is not allowed; see code for details"
            stop
         endif

         call read_transdip2 (iwrite, lutdip, ntset, no_target_states, numtgt(1 : csf_head_i % ntgsym), &
                              itarget_symmetry_order(1 : csf_head_i % ntgsym, :), ci_vec_j % nnuc, ismax, target_properties, &
                              gflag, prop_order, ifail)

         call continuum_continuum_only_one_drv(csf_head_j,csf_body_j, icontinuum_orbital_table_j)

      case(3)
!        Both I and J contain continuum orbitals
!        Read in the target CI vectors, we ensure that
!        only target symmetries included in the inner region
!        target state expansion are read in here.
         call read_all_ci_vectors(lutargci,csf_head_i%ntgsym,numtgt,itarget_symmetry_order,&
     &                                all_ntgsym,all_numtgt,all_ci_vec_target,iwrite,gflag)

         call prop_order_ci_vec(prop_order, all_ntgsym, sum(all_numtgt), all_ci_vec_target)

         call congen_to_target_energy_order_map(all_ci_vec_target,csf_head_i%ntgsym,numtgt,all_ntgsym, all_numtgt,&
     &                                          idtarg,ifail)

         allocate(ci_vec_target(csf_head_i%ntgsym))

         do itarg=1,csf_head_i%ntgsym
            ci_vec_target(itarg)=all_ci_vec_target(itarg)
            call all_ci_vec_target(itarg) % dealloc
         end do

!        Figure out the phase
         no_target_states=sum(numtgt(1:csf_head_i%ntgsym)) !ZM: sum only over the elements that have been defined
         allocate(csf_head_i%itarget_overall_phase(no_target_states),csf_head_j%itarget_overall_phase(no_target_states))
         csf_head_i%itarget_overall_phase=0
         csf_head_j%itarget_overall_phase=0
!ZM
         if (target_phase_correction) then

           !Determine overall target state phase introduced by the rediagonalisation of the target
           !Hamiltonian during the N+1 inner region Hamiltonian diagonalisation. Note, this is
           !done in scatci so as to deal with phase issues coming from 'dictionary' orbital ordering
           !and construction of spin orbitals as excitations from reference determinant between
           !the target and scattering CONFIG runs. This phase is also recorded in iphase, but
           !unfortunately an overall phase coming from the diagonalisation is not.

            icount=0
            itarg=0
            do i=1, csf_head_i%ntgsym
               do j=1,numtgt(i)
                  itarg=itarg+ci_vec_target(i)%nocsf
                  do k=1,csf_head_i%nctgt(i)
                     icount=icount+1

                  end do

               end do

            end do

            allocate(np1_calc_target_vectors(icount))
            np1_calc_target_vectors=0.0_idp
            call CWBOPN(lunp1targ(1))
            read(lunp1targ(1)) np1_calc_target_vectors
            close(lunp1targ(1))

            call determine_target_phase(csf_head_i%iphz,np1_calc_target_vectors,&
     &                ci_vec_target,csf_head_i%ntgsym,numtgt,&
     &                csf_head_i%itarget_overall_phase, ifail)
         else
            csf_head_i%itarget_overall_phase=1
         endif

         call phases_in_target_energy_order_map(ci_vec_target,csf_head_i%ntgsym,numtgt,&
     &             id_map_targ_to_congen,idtarg_test,itarget_irrep,target_spin,ifail)


         if(ibra_or_ket_bound.eq.2) call write_target_phases(csf_head_i%ntgsym,csf_head_i%numtgt,csf_head_i%itarget_overall_phase)

         allocate(csf_head_i%idtarg(no_target_states),csf_head_j%idtarg(no_target_states))
         csf_head_i%idtarg=idtarg
         csf_head_j%idtarg=idtarg_test

         if (mpi_i_write) then
            write(lutarg_phase,*) 'Target phases in state I'
            write(lutarg_phase,*) '------------------------'
            write(lutarg_phase,*) 'State No.    I.Rep.     Spin       Phase'
            do i=1,no_target_states
               write(lutarg_phase,'(i7,8x,i2,8x,f4.1,8x,i2)') i,itarget_irrep(id_map_targ_to_congen(i)),&
                                                              target_spin(id_map_targ_to_congen(i)),  &
                                                              csf_head_i%itarget_overall_phase(id_map_targ_to_congen(i))
            end do
         end if

         if (target_phase_correction) then

            np1_calc_target_vectors=0.0_idp
            call CWBOPN(lunp1targ(2))

            read(lunp1targ(2)) np1_calc_target_vectors
            close(lunp1targ(2))

            call determine_target_phase(csf_head_j%iphz,np1_calc_target_vectors,&
     &                ci_vec_target,csf_head_j%ntgsym,numtgt,&
     &                csf_head_j%itarget_overall_phase, ifail)
         else
            csf_head_j%itarget_overall_phase=1
         endif


         if(ibra_or_ket_bound.eq.1) call write_target_phases(csf_head_j%ntgsym,csf_head_j%numtgt,csf_head_j%itarget_overall_phase)

         if (mpi_i_write) then
            write(lutarg_phase,*) ''
            write(lutarg_phase,*) 'Target phases in state J'
            write(lutarg_phase,*) '------------------------'
            write(lutarg_phase,*) 'State No.    I.Rep.     Spin       Phase'
            do i=1,no_target_states
               write(lutarg_phase,'(i7,8x,i2,8x,f4.1,8x,i2)') i,itarget_irrep(id_map_targ_to_congen(i)),&
                                                              target_spin(id_map_targ_to_congen(i)),  &
                                                              csf_head_j%itarget_overall_phase(id_map_targ_to_congen(i))
            end do
         end if

!        Read in target properties
         no_target_states=sum(numtgt(1:csf_head_i%ntgsym))
         ismax=2
         ifail=0

         call read_transdip2 (iwrite, lutdip, ntset, no_target_states, numtgt(1 : csf_head_i % ntgsym), &
                              itarget_symmetry_order(1 : csf_head_i % ntgsym, :), ci_vec_i % nnuc, ismax, target_properties, &
                              gflag, prop_order, ifail)

!        Construct the density matrix blocks involving continuum orbitals
         call continuum_continuum_drv (csf_head_i, csf_head_j, csf_body_i, csf_body_j, &
                                       icontinuum_orbital_table_i, icontinuum_orbital_table_j, mpi_i_write, lucont_orbs)

      end select

!     Calculate 1 electron properties
!     ------------------------------------------------------------------

!     NEW MATRIX MULTIPLICATION SCHEME
      select case(icontains_continuum_orbitals)
      case(0)
         inumber_of_gcsf_i = csf_head_i%l2nocsf
         inumber_of_gcsf_j = csf_head_j%l2nocsf

      case(1)
         inumber_of_gcsf_i = size(icontinuum_orbital_table_i,1) + csf_head_i%l2nocsf
         ifirst_l2_csf_i=size(icontinuum_orbital_table_i,1) + 1
         inumber_of_gcsf_j = csf_head_j%l2nocsf

      case(2)
         inumber_of_gcsf_i = csf_head_i%l2nocsf
         inumber_of_gcsf_j = size(icontinuum_orbital_table_j,1) + csf_head_j%l2nocsf
         ifirst_l2_csf_j=size(icontinuum_orbital_table_j,1) + 1

      case(3)
         inumber_of_gcsf_i = size(icontinuum_orbital_table_i,1) + csf_head_i%l2nocsf
         ifirst_l2_csf_i=size(icontinuum_orbital_table_i,1) + 1
         inumber_of_gcsf_j = size(icontinuum_orbital_table_j,1) + csf_head_j%l2nocsf
         ifirst_l2_csf_j=size(icontinuum_orbital_table_j,1) + 1

      end select

!     L2 target contraction for removal of pseudoresonances if requested
!     ---------------------------------------------------------------
!     The transformed N+1 CI vectors (these must be stacked in sets in a single file)

      if (lucivec_l2contract .ne. 0) then
         select case(icontains_continuum_orbitals)
         case(0)
         case(1)
            call read_ci_vector(lucivec_l2contract,nciset_l2contract(1),ci_vec_contracted_i,0,iwrite)
            inumber_of_gcsf_i=ci_vec_contracted_i%nocsf
            call read_transformation_matrix(lutransmat_l2contract(1),l2block_contract_i)
            allocate(l2block_contract_transpose_i(size(l2block_contract_i,2),size(l2block_contract_i,1)))
            l2block_contract_transpose_i=transpose(l2block_contract_i)
            nstat(1) = ci_vec_contracted_i%nstat

         case(2)
            call read_ci_vector(lucivec_l2contract,nciset_l2contract(2),ci_vec_contracted_j,0,iwrite)
            inumber_of_gcsf_j=ci_vec_contracted_j%nocsf
            call read_transformation_matrix(lutransmat_l2contract(2),l2block_contract_j)
            nstat(2) = ci_vec_contracted_j%nstat

         case(3)
            call read_ci_vector(lucivec_l2contract,nciset_l2contract(1),ci_vec_contracted_i,0,iwrite)
            call read_ci_vector(lucivec_l2contract,nciset_l2contract(2),ci_vec_contracted_j,0,iwrite)
            inumber_of_gcsf_i=ci_vec_contracted_i%nocsf
            inumber_of_gcsf_j=ci_vec_contracted_j%nocsf
            call read_transformation_matrix(lutransmat_l2contract(1),l2block_contract_i)
            call read_transformation_matrix(lutransmat_l2contract(2),l2block_contract_j)
            allocate(l2block_contract_transpose_i(size(l2block_contract_i,2),size(l2block_contract_i,1)))
            l2block_contract_transpose_i=transpose(l2block_contract_i)
! TODO FM (NEEDS FIXING SO NUMBER OF STATES CAN BE DEFINED)
            nstat(1) = ci_vec_contracted_i%nstat
            nstat(2) = ci_vec_contracted_j%nstat

         end select

      end if

      if (nstat(1) < nstat(2)) then
         !ZM in this case we change the order of multiplications below so we need to allocate differently
         call properties_times_civec_j%init_CV(nstat(1),inumber_of_gcsf_j)
      else
         call properties_times_civec_j%init_CV(inumber_of_gcsf_i ,nstat(2))
      endif

      ! set up the distributed generalized CSF overlap matrix
      call allocate_prototype_property_matrix(prototype_properties_matrix, inumber_of_gcsf_i, inumber_of_gcsf_j, iwrite)

      ! The L2-L2 Block of the density matrix is calculated. For target runs this is
      ! The only part needed as no continuum orbitals are involved.
      call l2l2_drv(csf_head_i, csf_head_j, csf_body_i, csf_body_j, density_matrix_orbital_pairs_l2l2_1, &
                    density_matrix_orbital_pairs_l2l2_2, density_matrix_coefficients_l2l2, orbs_diag_l2l2, &
                    prototype_properties_matrix, inumber_of_gcsf_i, inumber_of_gcsf_j, sparsity)

      ! Now the continuum-L2 blocks (evaluate Dyson orbitals)
      if (csf_head_j % l2nocsf /= 0 .and. iand(BraState, icontains_continuum_orbitals) /= 0) then
         call continuum_l2_drv(ci_vec_target, csf_head_i, csf_head_j, csf_body_i, csf_body_j, &
                               icontinuum_orbital_table_i, dyson_orbitals_cl2, sparsity, &
                               prototype_properties_matrix, inumber_of_gcsf_j, BraState)
      end if
      if (csf_head_i % l2nocsf /= 0 .and. iand(KetState, icontains_continuum_orbitals) /= 0) then
         call continuum_l2_drv(ci_vec_target, csf_head_j, csf_head_i, csf_body_j, csf_body_i, &
                               icontinuum_orbital_table_j, dyson_orbitals_l2c, sparsity, &
                               prototype_properties_matrix, inumber_of_gcsf_i, KetState)
      end if

      ! Write dyson-like orbitals to file for input to dyson code (which transforms back to the
      ! the GTO basis). These are target states overlapped with the close coupling basis
      ! configurations (Generalised CSF i.e. gcsf).
      if (icontains_continuum_orbitals /= 0 .and. mpi_i_write) then
         ! Avoid calling the I/O routine in the distributed case, because then no process has complete arrays dyson_orbitals_xxx
         if (prototype_properties_matrix % CV_is_scalapack) then
            write (iwrite, '(A,I0,A)') '  Warning: Not writing Dyson orbitals in GCSF format (unit ', ludyson_gcsf, ') in parallel.'
         else
            ! Send to the subroutine only the L2-part of the dyson_orbitals_xxx arrays
            call write_dyson_gcsf_orbitals(ludyson_gcsf, icontains_continuum_orbitals, csf_head_i, csf_head_j, &
                                           icontinuum_orbital_table_i, icontinuum_orbital_table_j, &
                                           dyson_orbitals_cl2(:, ifirst_l2_csf_j:, :), &
                                           dyson_orbitals_l2c(:, ifirst_l2_csf_i:, :))
         end if
      end if

!     If we take states from bound we need to multiply the appropriate
!     symmetrys CI vectors by the bound coefficients.

      select case(iuse_bound)
      case(1)
         nstat(1)=states_from_bound%nstat
      case(2)
         nstat(2)=states_from_bound%nstat
      end select

!ZM TRANSI, TRANSJ are flags that tell me if the ci_vec_i%CV and ci_vec_j%CV arrays should be transposed when performing the matrix operations below.
!   From now on the arrays ci_vec_*%CV should always be used together with the flags TRANSI, TRANSJ.
      TRANSI = 'N'
      TRANSJ = 'N'
      if  (lucivec_l2contract .ne. 0) then
         allocate(ci_vec_j%CV,source=ci_vec_contracted_j%CV)
         TRANSI = 'T' ! Neutral: This is case 2
         TRANSJ = 'N'
      else

   !     We need to multiply in the bound state if given
         select case(iuse_bound)
         case(0)
            TRANSI = 'T'
            TRANSJ = 'N'
         case(1)
            !ZM: use the properties of matrix transposition to rearrange the multiplications so the additional transposition is not necessary just like in the cas(2) below.
            call temp_ci_vector_i%init_CV(inumber_of_gcsf_i,nstat(1))
            call temp_ci_vector_i%A_B_matmul(ci_vec_i,states_from_bound,'N','N')
            call ci_vec_i%final_CV
            call ci_vec_i%init_CV(inumber_of_gcsf_i,nstat(1))
            ci_vec_i%CV = temp_ci_vector_i%CV
            ci_vec_i%ei(1:nstat(1)) = states_from_bound%ei(1:nstat(1))
            ci_vec_i%e0 = 0.0_wp !the bound state energies from BOUND contain the nuclear repulsion energy
            write(iwrite,'("Bound state energies: ",100e25.15)') ci_vec_i%ei(1:nstat(1))
            call temp_ci_vector_i%final_CV
            TRANSI = 'T'
            TRANSJ = 'N'
         case(2)
            call temp_ci_vector_j%init_CV(inumber_of_gcsf_j,nstat(2))
            call temp_ci_vector_j%A_B_matmul(ci_vec_j,states_from_bound,'N','N')
            call ci_vec_j%final_CV
            call ci_vec_j%init_CV(inumber_of_gcsf_j,nstat(2))
            ci_vec_j%CV = temp_ci_vector_j%CV
            ci_vec_j%ei(1:nstat(2)) = states_from_bound%ei(1:nstat(2))
            ci_vec_j%e0 = 0.0_wp !the bound state energies from BOUND contain the nuclear repulsion energy
            write(iwrite,'("Bound state energies: ",100e25.15)') ci_vec_j%ei(1:nstat(2))
            call temp_ci_vector_j%final_CV
            TRANSI = 'T'
            TRANSJ = 'N'
         end select

   !     Construct and write out the Dyson orbitals
   !     ZM: the Dyson orbitals and their symmetries are returned in the arrays dyson_orbitals,idyson_orbital_irrep,idyson_itarget_irrep,idyson_itarget_spin
         if (icontains_continuum_orbitals .ne. 0) then
            !
            if (allocated(bound_vec_i)) deallocate(bound_vec_i)
            if (allocated(bound_vec_j)) deallocate(bound_vec_j)
            !
            select case(ibra_or_ket_bound)
               case (1)
                  call ci_vec_i%gather_vectors(bound_vec_i,nstat(1),-1,-1)
               case (2)
                  call ci_vec_j%gather_vectors(bound_vec_j,nstat(2),-1,-1)
            end select
            call write_dyson_orbitals (ludyson, icontains_continuum_orbitals, iuse_bound, ibra_or_ket_bound, csf_head_i, &
                                       csf_head_j, icontinuum_orbital_table_i, icontinuum_orbital_table_j, dyson_orbitals_cl2, &
                                       dyson_orbitals_l2c, bound_vec_i, nstat(1), TRANSI, bound_vec_j, nstat(2), TRANSJ, &
                                       dyson_orbitals, idyson_orbital_irrep, idyson_itarget_irrep, idyson_itarget_spin, &
                                       dyson_orbital_norms, prototype_properties_matrix, mpi_i_write, iwrite)
            if (ukrmolp_ints) then
               if (scattering_calculation) then
                  !Use the UKRmol interface module to write out the Dyson orbitals into the Molden file using the integral code's objects.
                  call write_molden_dyson_orbitals(dyson_orbitals, idyson_orbital_irrep, idyson_itarget_irrep, &
                                                   idyson_itarget_spin, dyson_orbital_norms)
               end if
            else
               write(iwrite,'(10X,"Using SWEDEN integrals and orbitals: Dyson orbitals will NOT be saved in the Molden format.")')
            endif
         end if

      end if

      !ZM moved allocation from inside the k-loop: nstat(:) is always the same.
      call properties_temp%init_CV(nstat(1),nstat(2))

      !ZM construct table mapping property index to l,m values
      allocate(lm_to_l_m(2,(max_multipole+1)**2))
      k = 0
      do l=0,max_multipole
         do m=-l,l
            k = k + 1
            lm_to_l_m(1,k) = l
            lm_to_l_m(2,k) = m
         enddo !m
      enddo !l

      ! decide which molecular_properties_data object to use (given or custom)
      if (present(mol_prop_data)) then
         mol_prop_data_ptr => mol_prop_data
      else
         call mol_prop_data_loc % preallocate_property_blocks((max_multipole + 1)**2)
         mol_prop_data_ptr => mol_prop_data_loc
      end if

      ! store the requested nuclear multipole flag within the properties class
      mol_prop_data_ptr % nuclear_prop_flag = isw

      ! loop over all requested properties
      do k = min_multipole**2 + 1, (max_multipole + 1)**2

         !ZM find out if this multipole moment can be non-zero:
         if (ukrmolp_ints) then
            l = lm_to_l_m(1,k)
            m = lm_to_l_m(2,k)
            !IRR of the property:
            i = molecular_orbital_basis%symmetry_data%sph_harm_pg_sym(l,m,nam)
            !Wigner-Eckart test using D2h product table:
            if (IPD2H(i,csf_head_i%MGVN+1) .ne. csf_head_j%MGVN+1) then
               write(iwrite,*) 'Skipping property',l,m,'for this pair of symmetries'
               cycle
            endif
         endif

!        Create the prototype properties matrix
!        ---------------------------------------------------------------

         call create_prototype_properties_matrix (csf_head_i, csf_head_j, icontinuum_orbital_table_i, icontinuum_orbital_table_j, &
                target_properties, dyson_orbitals_cl2, dyson_orbitals_l2c, density_matrix_orbital_pairs_l2l2_1, &
                density_matrix_orbital_pairs_l2l2_2, density_matrix_coefficients_l2l2, orbs_diag_l2l2, pintegrals, k, &
                icontains_continuum_orbitals, prototype_properties_matrix, ukrmolp_ints)


!        L2 target contraction for removal of pseudoresonances if requested
!        ---------------------------------------------------------------

         if (lucivec_l2contract .ne. 0) then
!        Transform prototype_properties matrix

            if (allocated(prototype_properties_matrix_temp%CV)) deallocate(prototype_properties_matrix_temp%CV)
            allocate(prototype_properties_matrix_temp%CV(inumber_of_gcsf_i,inumber_of_gcsf_j) )
            prototype_properties_matrix_temp%CV=0._idp

            select case(icontains_continuum_orbitals)
            case(1)
               prototype_properties_matrix_temp%CV=prototype_properties_matrix%CV(1:csf_head_i%last_continuum_csf,:)
               call maths_matrix_multiply_blas95 (l2block_contract_transpose_i, &
                                                  prototype_properties_matrix % CV(csf_head_i % last_continuum_csf + 1 :, :), &
                                                  prototype_properties_matrix_temp % CV(csf_head_i % last_continuum_csf + 1 :, :))

            case(2)
               prototype_properties_matrix_temp%CV=prototype_properties_matrix%CV(:,1:ifirst_l2_csf_j-1)
               call maths_matrix_multiply_blas95 (prototype_properties_matrix % CV(:, ifirst_l2_csf_j :), &
                                                  l2block_contract_j, &
                                                  prototype_properties_matrix_temp % CV(:, ifirst_l2_csf_j :))

            case(3)
               prototype_properties_matrix_temp % CV = &
               prototype_properties_matrix % CV(1 : size(icontinuum_orbital_table_i,1), 1 : size(icontinuum_orbital_table_j,1))

            end select

            deallocate(prototype_properties_matrix%CV)
            allocate(prototype_properties_matrix%CV(inumber_of_gcsf_i,inumber_of_gcsf_j))
            prototype_properties_matrix%CV=prototype_properties_matrix_temp%CV

         end if

         write(iwrite,'(/,1x,a8," prototype property matrix constructed.")')  pintegrals%property_name(k)

!        Multiply in the CI vectors to give the dipoles between inner region states
!        ---------------------------------------------------------------
         properties_times_civec_j%CV=0.0_idp
         properties_temp%CV=0.0_idp

         TRANS = 'N'
         !ZM this if statement must be consistent with the one above where properties_times_civec_j%CV is allocated.
         if (nstat(1) < nstat(2)) then
            !ZM multiply-in the I-CI vector first so we have matrix x vector in both cases
            call properties_times_civec_j%A_B_matmul(ci_vec_i,prototype_properties_matrix,TRANSI,TRANS) !properties_times_civec_j = matmul(ci_vec_i,prototype_properties_matrix)
            call properties_temp%A_B_matmul(properties_times_civec_j,ci_vec_j,TRANS,TRANSJ)             !properties_temp = matmul(properties_times_civec_j,ci_vec_j)
         else
            !ZM multiply-in the J-CI vector first
            call properties_times_civec_j%A_B_matmul(prototype_properties_matrix,ci_vec_j,TRANS,TRANSJ) !properties_times_civec_j = matmul(prototype_properties_matrix,ci_vec_j)
            call properties_temp%A_B_matmul(ci_vec_i,properties_times_civec_j,TRANSI,TRANS)             !properties_temp = matmul(ci_vec_i,properties_times_civec_j)
         endif

         ! Apply the spin part overlap factor: At the moment, the property operator has no spin dependence,
         ! so the spin part overlap is either zero or one.
         if (ci_vec_i%s /= ci_vec_j%s .or. ci_vec_i%sz /= ci_vec_j%sz) then
            properties_temp%CV = 0
         end if

         ! Addition of the nuclear dipole/quadrupole term. The nuclear dipole/quadrupole term is non-zero
         ! only when <psi_i |r | psi_j> .ne. 0 and  psi_i=psi_j.
         if (isw == 0 .and. ci_vec_i%mgvn == ci_vec_j%mgvn .and. ci_vec_i%s == ci_vec_j%s .and. ci_vec_i%sz == ci_vec_j%sz) then
            call add_nuclear(ci_vec_i, ci_vec_j, properties_temp, k, ukrmolp_ints)
         end if

         ! Store property matrix for this l,m
         call mol_prop_data_ptr % insert_property_block(csf_head_i, ci_vec_i, nstat(1), &
                                                        csf_head_j, ci_vec_j, nstat(2), &
                                                        threshold, pintegrals % lp(k), pintegrals % mp(k) * pintegrals % qp(k), &
                                                        properties_temp)

         write(iwrite,'(/,1x,a8," for all state pairs calculated.", 3i5)') &
            pintegrals % property_name(k), k, pintegrals % lp(k), pintegrals % mp(k) * pintegrals % qp(k)

      end do !k

!     Write to file or to mol_prop_data
      if ( .not. present(mol_prop_data) ) then
         !ZM todo: it would be more memory-efficient if the properties were written one-by-one (i.e. in the loop over k)
         !         rather than all properties (dipole components, quadrupoles, etc.) at once.
         call mol_prop_data_ptr % write_properties(lupropw, ukrmolp_ints)

         !we can finish now since if cdenprop was supposed to be called multiple
         !times then it should be using the mol_prop_data on input and therefore
         !the other branch of the outer if statement.
         if (ukrmolp_ints) then
            call destroy_ukrmolp
            call mpi_mod_finalize
         end if
      end if

      if (.not.ukrmolp_ints) close(luprop) !close SWEDEN property integrals file

      deallocate(lm_to_l_m)
      call prototype_properties_matrix%final_CV
      call properties_temp%final_CV
      call properties_times_civec_j%final_CV

      if (present(mol_prop_data)) mol_prop_data % pintegrals => pintegrals
      nullify (pintegrals)

      end subroutine cdenprop_drv

      subroutine read_I_J_states (icontains_continuum_orbitals, iuse_bound, lucivec, nciset, &
                                  luneutralci, nbset, bform, nstat, ci_vec_i, ci_vec_j, states_from_bound)
         use cdenprop_defs, only: CIvect
         use cdenprop_io, only: read_ci_vector, read_states_from_bound
         implicit none
         integer, intent(in) :: lucivec(2), luneutralci, nbset, iuse_bound, icontains_continuum_orbitals
         character(len=11) :: bform
         integer, dimension(tsym_max), intent(in) :: nciset, nstat
         type (CIvect), intent(inout) :: ci_vec_i,ci_vec_j,states_from_bound !The N+1 CI vectors

         integer :: i

         select case(icontains_continuum_orbitals)
         case(0)
            write(iwrite,'(/, " READ TARGET CI FILES")')
            write(iwrite,'(   " --------------------")')
         case(1:)
            write(iwrite,'(/, " READ N+1 CI FILES")')
            write(iwrite,'(   " -----------------")')
         end select

         !ZM error check on lucivec(1:2):
         if (lucivec(1) .le. 0) then
            if (allocated(ci_vec_i%CV)) then
               write(iwrite,'(/," Dimensions of the CI vector matrix in ci_vec_i: ",2i10)') size(ci_vec_i%CV,1), size(ci_vec_i%CV,2)
            else
               write(iwrite,'(" On input the matrix of CI vectors from the argument ci_vec_i was not allocated")')
               stop
            endif
         endif
         if (lucivec(2) .le. 0) then
            if (allocated(ci_vec_j%CV)) then
               write(iwrite,'(/," Dimensions of the CI vector matrix in ci_vec_j: ",2i10)') size(ci_vec_j%CV,1), size(ci_vec_j%CV,2)
            else
               write(iwrite,'(" On input the matrix of CI vectors from the argument ci_vec_j was not allocated")')
               stop
            endif
         endif

         if (lucivec(1) > 0) then
            call read_ci_vector(lucivec(1),nciset(1),ci_vec_i,nstat(1),iwrite)
         else
            if (ci_vec_i%CV_is_scalapack) then
               i = ci_vec_i%mat_dimen
            else
               i = size(ci_vec_i%cv,2)
            endif
            if (i < nstat(1)) then
               write(iwrite,'(" Case 0: the number of states in ci_vec_i%cv is smaller than requested: ",2i10)') i,nstat(1)
               stop
            endif
         endif

         if (lucivec(2) > 0) then
            call read_ci_vector(lucivec(2),nciset(2),ci_vec_j,nstat(2),iwrite)
         else
            if (ci_vec_j%CV_is_scalapack) then
               i = ci_vec_j%mat_dimen
            else
               i = size(ci_vec_j%cv,2)
            endif
            if (i < nstat(2)) then
               write(iwrite,'(" Case 0: the number of states in ci_vec_j%cv is smaller than requested: ",2i10)') i,nstat(2)
               stop
            endif
         endif

        !Bound states from bound on a separate unit luneutralci
        if (iuse_bound .ne. 0) then
           write(iwrite,'(/, " READ STATES FROM BOUND")')
           write(iwrite,'(   " ----------------------")')
           call read_states_from_bound(luneutralci,nbset,bform,states_from_bound,iwrite)
        end if

      end subroutine read_I_J_states


    !> \brief   Add nuclear multipole to properties
    !> \authors A Harvey, J Benda
    !> \date    2011 - 2019
    !>
    !> Add nuclear multipole contribution to the diagonal properties. This subroutine only operates on the subset
    !> of the distributed property matrix available to the current process.
    !>
    !> \param[in]    ci_vec_i  Bra vector; only used to get nuclear information.
    !> \param[in]    ci_vec_j  Ket vector; not used at all.
    !> \param[inout] props     Distributed dense property state-to-state matrix.
    !> \param[in]    iprop     Label of the Ylm property (1 = Y00, 2 = Y1-1, 3 = Y10, ...)
    !> \param[in]    ukrmolp_ints  Quadrupole convention flag (SWEDEN vs. GBTOlib).
    !>
    subroutine add_nuclear (ci_vec_i, ci_vec_j, props, iprop, ukrmolp_ints)

        use cdenprop_defs, only: CIvect
        use ukrmol_interface_gbl, only: GET_ECP_CORE_CHARGES

        type(CIvect), intent(in)    :: ci_vec_i, ci_vec_j
        type(CIvect), intent(inout) :: props
        integer,      intent(in)    :: iprop
        logical,      intent(in)    :: ukrmolp_ints

        integer   :: nnuc, i, inuc, q
        real(idp) :: x, y, z, Ylm
        integer, allocatable :: core_charges(:)

        if (iprop < 2) return

        nnuc = ci_vec_i % nnuc

        if (ukrmolp_ints) then
           call GET_ECP_CORE_CHARGES(core_charges)
        else
           allocate(core_charges(nnuc))
           core_charges = 0
        endif

        do inuc = 1, nnuc
            x = ci_vec_i % xnuc(inuc)
            y = ci_vec_i % ynuc(inuc)
            z = ci_vec_i % znuc(inuc)
            q = ci_vec_i % charge(inuc) - core_charges(inuc)
            Ylm = multipole_nuclear_term(iprop, ukrmolp_ints, x, y, z)
            do i = 1, min(props % mat_dimen_r, props % mat_dimen_c)
                call props % add_to_CV_element(-q * Ylm, i, i)
            end do
        end do

    end subroutine add_nuclear


      subroutine  cdenprop_orb_find(indo,nodo,ndo,lastcont,icdo,cdo)
      use cdenprop_defs
      implicit none
      integer :: lastcont
      integer, dimension(:) :: indo,ndo,nodo, orbs(2)
      integer, optional,dimension(:) :: icdo
      real(kind =idp), optional,dimension(:) :: cdo !testing
      integer :: nocsf
      integer :: i,j,k,m,n, nexa,nexb, jj,kk

!      INPUT
!      -----
!      INDO: Gives the array index for 1st element of ith CSF. (NOCSF+1) last values is lndof
!      NODO:   no. of determinants per CSF. I(1:nocsf)
!      NDO:  For each CSF: 1st element(indexed by INDO) gives number of excitations from the reference determinant
!      e.g N, this is followed by the N orbitals excited from, then by the N orbitals excited to. This is repeated
!      for each determinant contained in the CSF, so the total number of array elements for the ith CSF is given by
!      (2*NDO(INDO(i))+1)*NODO(i).  I(1: )
!      Last continuum CSF

      nocsf=size(nodo)

!      write configurations

      do i=1,nocsf
         write(lastcont,'(/" CSF NO. ",I7, " INVOLVING ",i3," DETERMINANTS")') i, nodo(i)
         k=indo(i)
         write(lastcont,'(/"  Space-Spin")')
         if(present(cdo)) then
            n=icdo(i)

         end if

         do j=1,nodo(i)
            if(present(cdo)) then
               write(lastcont,'(I3,D20.5,20(I3,1x))') ndo(k),cdo(n),(ndo(k+m), m=1,2*(ndo(k)))
               k=k+2*(ndo(k))+1
               n=n+1

            else
               write(lastcont,'(20(I3,1x))') ndo(k),(ndo(k+m), m=1,2*(ndo(k)))
               k=k+2*(ndo(k))+1

            end if

         end do

         k=indo(i)
         write(lastcont,'(/"  Space only")')
         do j=1,nodo(i)
            write(lastcont,'(20(I3,1x))') ndo(k),((ndo(k+m)+1)/2, m=1,2*(ndo(k)))
            k=k+2*(ndo(k))+1

         end do

      end do

      end subroutine cdenprop_orb_find



      subroutine continuum_continuum_drv (csf_head_i, csf_head_j, csf_body_i, csf_body_j, &
                                          icontinuum_orbital_table_i, icontinuum_orbital_table_j, mpi_i_write, lucont_orbs)
      use cdenprop_defs
      use cdenprop_io
      implicit none
!     ******************************************************************
!
!     Constructs the density matrix block corresponding to continuum
!     orbitals in both configuration sets.
!     Many of the routines are taken from DENPROP
!     Note: In contrast to DENPROP the density matrix between a state
!     and itself is not constructed.
!
!     ******************************************************************


      integer, dimension(tsym_max) :: nciset,nstat
      integer :: lucivec, lucont_orbs
      logical :: mpi_i_write

!     The CSF files.
      type (CSFheader) :: csf_head_i,csf_head_j
      type (CSFbody)   :: csf_body_i,csf_body_j

      integer :: i,j,k,iorb,jorb,no_target_symmetries, itarget_symmetry,  no_continuum_orbitals, &
     &           no_CSFS_in_contracted_basis_i,no_CSFS_in_contracted_basis_j, &
     &           no_CSF_in_continuum_part_i, no_CSF_in_continuum_part_j,icontinuum_orbital_index,&
     &           icsf_index,itarget_index

      integer, dimension(tsym_max) ::  icontinuum_orbital_symmetry, &
     &                         no_targ_states_per_symmetry, no_cont_orbitals_per_symmetry, &
     &                         no_bound_orbitals_per_symmetry, ibound_orbital_start_indices, &
     &                         icontinuum_orbital_start_indices
      integer, allocatable, dimension(:,:) :: icontinuum_orbital_table_i,icontinuum_orbital_table_j
      INTEGER, allocatable, dimension(:,:,:) :: density_matrix_orbital_pairs_cc


!     First we make a continnum orbital table that links each GCSF
!     (generalised configuration state function) to a
!     continuum orbital. We work entirely in spatial orbitals as
!     spin does not need to be considered for the overlap of
!     two GCSFs as the overlap of the spin coupled targets gives
!     one or zero simply leaving the the continuum orbitals as the
!     unmatched pair.

!     number of states in the CI contracted basis

      no_CSF_in_continuum_part_i=dot_product(csf_head_i%numtgt,csf_head_i%notgt)
      no_CSF_in_continuum_part_j=dot_product(csf_head_j%numtgt,csf_head_j%notgt)

      no_CSFS_in_contracted_basis_i=csf_head_i%l2nocsf +  no_CSF_in_continuum_part_i
      no_CSFS_in_contracted_basis_j=csf_head_j%l2nocsf +  no_CSF_in_continuum_part_j

      write(iwrite,'(/," Size of contracted basis for I is ",I6," and for J is",I6,/ )') &
          no_CSFS_in_contracted_basis_i,no_CSFS_in_contracted_basis_j
      write(iwrite,'(  " Size of continuum part of the basis for I is ",I6," and for J is",I6,/ )')  &
          no_CSF_in_continuum_part_i,no_CSF_in_continuum_part_j
!     we need the target symmetry order


!     We assume that the rest of the arrays that denote
!     a quantity per symmetry are ordered in the
!     target symmetry order (I think this is how they
!     are stored in the files - double check)

!     we need the number of targets states per symmetry
      no_targ_states_per_symmetry=0

!     then we need the symmetry of the continuum orbitals
!     for each target symmetry
      icontinuum_orbital_symmetry=0

!     number of continuum orbitals per symmetry
      no_cont_orbitals_per_symmetry=0
      no_continuum_orbitals=sum(no_cont_orbitals_per_symmetry)

!     number of bound orbitals per symmetry
      no_bound_orbitals_per_symmetry=0

!     determine the start points for the bound  and continuum orbitals

      call calc_orbital_start_indices(csf_head_i,icontinuum_orbital_start_indices, ibound_orbital_start_indices)

!     now we want a continuum orbital table
      allocate(icontinuum_orbital_table_i(no_CSF_in_continuum_part_i,6),icontinuum_orbital_table_j(no_CSF_in_continuum_part_j,6))
      icontinuum_orbital_table_i=0
      icontinuum_orbital_table_j=0


      icsf_index=1
      itarget_index=1
      do i=1,csf_head_i%ntgsym
         do j=1,csf_head_i%numtgt(i)
            icontinuum_orbital_index=icontinuum_orbital_start_indices(csf_head_i%mcont(i)+1)
            do k=1,csf_head_i%notgt(i)
               icontinuum_orbital_table_i(icsf_index,1)=itarget_index
               icontinuum_orbital_table_i(icsf_index,2)=j !Target number within symmetry
               icontinuum_orbital_table_i(icsf_index,3)=csf_head_i%itarget_symmetry_order(i,1) !Target irrep.
               icontinuum_orbital_table_i(icsf_index,4)=csf_head_i%mcont(i) ! Continuum smmetry
               icontinuum_orbital_table_i(icsf_index,5)=icontinuum_orbital_index
               icontinuum_orbital_table_i(icsf_index,6)=i !Symmetry index
               icontinuum_orbital_index=icontinuum_orbital_index+1
               icsf_index=icsf_index+1

            end do
            itarget_index=itarget_index+1

         end do

      end do

      if (mpi_i_write) then
         write(lucont_orbs,*) "Continuum orbital table - State I"
         do i = 1, no_CSF_in_continuum_part_i
            write(lucont_orbs, '(i5,i5,i5,i5,i5,i5)') (icontinuum_orbital_table_i(i,j), j=1,6)
         end do
         write(lucont_orbs,*) ' '
      end if

      icsf_index=1
      itarget_index=1

      do i=1,csf_head_j%ntgsym
         do j=1,csf_head_j%numtgt(i)
            icontinuum_orbital_index=icontinuum_orbital_start_indices(csf_head_j%mcont(i)+1)
            do k=1,csf_head_j%notgt(i)
               icontinuum_orbital_table_j(icsf_index,1)=itarget_index
               icontinuum_orbital_table_j(icsf_index,2)=j !Target number within symmetry
               icontinuum_orbital_table_j(icsf_index,3)=csf_head_j%itarget_symmetry_order(i,1) !Target irrep
               icontinuum_orbital_table_j(icsf_index,4)=csf_head_j%mcont(i) ! Continuum smmetry
               icontinuum_orbital_table_j(icsf_index,5)=icontinuum_orbital_index
               icontinuum_orbital_table_j(icsf_index,6)=i !Symmetry index
               icontinuum_orbital_index=icontinuum_orbital_index+1
               icsf_index=icsf_index+1

            end do
            itarget_index=itarget_index+1

         end do

      end do

      if (mpi_i_write) then
         write(lucont_orbs,*) "Continuum orbital table - State J"
         do i =1, no_CSF_in_continuum_part_j
            write(lucont_orbs, '(i5,i5,i5,i5,i5,i5)') (icontinuum_orbital_table_j(i,j), j=1,6)
         end do
      end if

      end subroutine continuum_continuum_drv



      subroutine continuum_continuum_only_one_drv(csf_head_i,csf_body_i, icontinuum_orbital_table_i)
      use cdenprop_defs
      use cdenprop_io
      implicit none
!     ******************************************************************
!
!     Constructs the density matrix block corresponding to continuum
!     orbitals in both configuration sets.
!     Many of the routines are taken from DENPROP
!     Note: In contrast to DENPROP the density matrix between a state
!     and itself is not constructed.
!
!     ******************************************************************


      integer, dimension(tsym_max) :: nciset,nstat
      integer :: lucivec

!     The CSF files.
      type (CSFheader) :: csf_head_i,csf_head_j
      type (CSFbody)   :: csf_body_i,csf_body_j

      integer :: i,j,k,iorb,jorb,no_target_symmetries, itarget_symmetry,  no_continuum_orbitals, &
     &           no_CSFS_in_contracted_basis_i,no_CSFS_in_contracted_basis_j, &
     &           no_CSF_in_continuum_part_i, no_CSF_in_continuum_part_j,icontinuum_orbital_index,&
     &           icsf_index,itarget_index

      integer, dimension(tsym_max) :: icontinuum_orbital_symmetry, &
     &                         no_targ_states_per_symmetry, no_cont_orbitals_per_symmetry, &
     &                         no_bound_orbitals_per_symmetry, ibound_orbital_start_indices, &
     &                         icontinuum_orbital_start_indices
      integer, allocatable, dimension(:,:) :: icontinuum_orbital_table_i,icontinuum_orbital_table_j
      INTEGER, allocatable, dimension(:,:,:) :: density_matrix_orbital_pairs_cc


!     First we make a continnum orbital table that links each GCSF
!     (generalised configuration state function) to a
!     continuum orbital. We work entirely in spatial orbitals as
!     spin does not need to be considered for the overlap of
!     two GCSFs as the overlap of the spin coupled targets gives
!     one or zero simply leaving the the continuum orbitals as the
!     unmatched pair.

!     number of states in the CI contracted basis

      no_CSF_in_continuum_part_i=dot_product(csf_head_i%numtgt,csf_head_i%notgt)

      no_CSFS_in_contracted_basis_i=csf_head_i%l2nocsf +  no_CSF_in_continuum_part_i

      write(iwrite,'(/," Size of contracted basis is ",I6/ )') no_CSFS_in_contracted_basis_i!,no_CSFS_in_contracted_basis_j
      write(iwrite,'(  " Size of continuum part of the basis is ",I6/ )')  no_CSF_in_continuum_part_i!,no_CSF_in_continuum_part_j
!     we need the target symmetry order


!     We assume that the rest of the arrays that denote
!     a quantity per symmetry are ordered in the
!     target symmetry order (I think this is how they
!     are stored in the files - double check)

!     we need the number of targets states per symmetry
      no_targ_states_per_symmetry=0

!     then we need the symmetry of the continuum orbitals
!     for each target symmetry
      icontinuum_orbital_symmetry=0

!     number of continuum orbitals per symmetry
      no_cont_orbitals_per_symmetry=0
      no_continuum_orbitals=sum(no_cont_orbitals_per_symmetry)

!     number of bound orbitals per symmetry
      no_bound_orbitals_per_symmetry=0

!     determine the start points for the bound  and continuum orbitals

      call calc_orbital_start_indices(csf_head_i,icontinuum_orbital_start_indices, ibound_orbital_start_indices)

!     now we want a continuum orbital table
      allocate(icontinuum_orbital_table_i(no_CSF_in_continuum_part_i,6))!,icontinuum_orbital_table_j(no_CSF_in_continuum_part_j,5))
      icontinuum_orbital_table_i=0


      icsf_index=1
      itarget_index=1
      do i=1,csf_head_i%ntgsym
         do j=1,csf_head_i%numtgt(i)
            icontinuum_orbital_index=icontinuum_orbital_start_indices(csf_head_i%mcont(i)+1)
            do k=1,csf_head_i%notgt(i)
               icontinuum_orbital_table_i(icsf_index,1)=itarget_index
               icontinuum_orbital_table_i(icsf_index,2)=j !Target number within symmetry
               icontinuum_orbital_table_i(icsf_index,3)=csf_head_i%itarget_symmetry_order(i,1) !Target irrep
               icontinuum_orbital_table_i(icsf_index,4)=csf_head_i%mcont(i) ! Continuum smmetry
               icontinuum_orbital_table_i(icsf_index,5)=icontinuum_orbital_index
               icontinuum_orbital_table_i(icsf_index,6)=i !Symmetry index
               icontinuum_orbital_index=icontinuum_orbital_index+1
               icsf_index=icsf_index+1

            end do
            itarget_index=itarget_index+1

         end do

      end do
! write(99,*) 'mcont', csf_head_i%mcont

      do i =1, no_CSF_in_continuum_part_i
         write(222, '(i5,i5,i5,i5,i5)') (icontinuum_orbital_table_i(i,j), j=1,6)

      end do
      write(222,*) ' '

      end subroutine continuum_continuum_only_one_drv



      subroutine calc_orbital_start_indices(csf_head,icontinuum_orbital_start_indices, ibound_orbital_start_indices)
      use cdenprop_defs
      type (CSFheader) :: csf_head
      integer, dimension(:) :: icontinuum_orbital_start_indices,ibound_orbital_start_indices
      integer :: i, num_ir

      num_ir=csf_head%nsym

!     Determine the index for the first bound spatial orbital in each symmetry
       ibound_orbital_start_indices(1)=1
       do i=2,num_ir
          ibound_orbital_start_indices(i)=ibound_orbital_start_indices(i-1)+csf_head%nob(i-1)

       end do

       do i=1,num_ir
          icontinuum_orbital_start_indices(i)= ibound_orbital_start_indices(i)+(csf_head%nob0(i))+csf_head%no_l2_virtuals(i)

       end do

      end subroutine calc_orbital_start_indices



!     Diagonal matrix indices functions
      function ij_to_index_lower_rowwise(i,j)
      implicit none
      integer ::ij_to_index_lower_rowwise,i,j,itemp, jtemp
      itemp=i
      jtemp=j
      if (i .lt. j)then
!       if (i .gt. j)then
         itemp=j
         jtemp=i

      end if

      ij_to_index_lower_rowwise=((itemp-1)*(itemp))/2 +jtemp
!       ij_to_index_lower_rowwise=((jtemp-1)*(jtemp))/2 +itemp


      end function ij_to_index_lower_rowwise


      function ij_to_index_full_colwise(i,j,size_matrix)
      implicit none
      integer ::ij_to_index_full_colwise,i,j,size_matrix
      intent(in):: i,j,size_matrix

      ij_to_index_full_colwise=(j-1)*size_matrix +i

      end function ij_to_index_full_colwise


      subroutine symmetry_block_to_integral_table(istart,nob,nsym)
      implicit none
      integer :: nsym
      integer,dimension(:) :: nob
      integer,allocatable,dimension(:,:) :: istart

      integer :: isym,jsym,last,m,n,k

!.... Set up pointers to start of each block of transformed integrals
!
!     Note last element of istart points to the last element of the
!     the integral table

      allocate(istart(nsym,nsym))
      istart=0

!     Diagonal symmetry blocks come first (main diagonal)
!     they are stored in column-wise, upper triangle form
!     ( or row-wise, lower triangle)
!     Example for 4 symmetries
!     Symmetry block order is
!     1,1, 2,2, 3,3, 4,4,
! write(99,*) 'nsym in symmetry_block_to_integral_table routine', nsym
! write(99,*) 'nob in symmetry_block_to_integral_table routine', nob

      last=0
      k=1
      do m=1, nsym
         istart(m,m)=last+1
         last=last+nob(m)*(nob(m)+1)/2

      end do

!     Now come the off diagonal elements
!     We loop over the diagonals
!     Example for 4 symmetries
!     Symmetry blocks order is
!     1,2, 2,3, 3,4, 1,3, 2,4, 1,4,

      do k=1, nsym-1
         do m=1,nsym-k
            istart(m,m+k)=last+1
            last=last+nob(m)*nob(m+k)

         end do

      end do

!     Note: This ordering can be obtained by cyclic permutation of j
!           and dropping i>j pairs, in the following manner
!
!     i=1,2,3,4
!     j=1,2,3,4

!       1,1,1,1, MDEL+1

!     i=1,2,3,  4
!     j=2,3,4,  1

!       2,4,2,   MDEL+1

!     i=1,2,  3,4
!     j=3,4,  1,2
!
!       3,3,     MDEL+1

!     i=1,  2,3,4
!     j=4,  1,2,3

!       4,       MDEL+1

      end subroutine symmetry_block_to_integral_table


      subroutine l2l2_drv (csf_head_i, csf_head_j, csf_body_i, csf_body_j, density_matrix_orbital_pairs_l2l2_1, &
                           density_matrix_orbital_pairs_l2l2_2, density_matrix_coefficients_l2l2, orbs_diag_l2l2, &
                           pmat, inumber_of_gcsf_i, inumber_of_gcsf_j, sparsity)

      use blas_lapack_gbl, only: blasint
      use cdenprop_defs,   only: CSFheader, CSFbody, CIvect
      use cdenprop_aux,    only: mkorbs, change_refdet, dryrun, makemg, tmtma_sparse

      use class_COOSparseMatrix_real,    only: COOMatrix_real
      use class_COOSparseMatrix_integer, only: COOMatrix_integer

      implicit none
!     ******************************************************************
!
!     Constructs the density matrix for configurations that do not
!     contain any continuum orbitals
!     Many of the routines are taken from DENPROP
!     Note: In contrast to DENPROP the density matrix between a state
!     and itself is not constructed.
!
!     ******************************************************************

      integer :: l2nocsf, lastcont, no_spin_orbitals,maxso,npflg,isymtype, n_cont_gcsf_i, n_cont_gcsf_j
      integer :: i, j, k, m, n, li, lj, gi, gj, index_ndo_start_i, index_ndo_start_j, index_cdo_start_i,&
                 index_cdo_start_j, index_ndo_end_i, index_ndo_end_j, inumber_of_gcsf_i, inumber_of_gcsf_j, nestimated, nproc
      integer, allocatable, dimension(:) :: list_space,list_ug, list_quant_no, list_spin
      integer, allocatable, dimension(:) :: ndtrf_point_back
      integer, dimension(tsym_max) :: lucsf ! Input CSF file units
      logical :: evalue
      REAL(KIND=idp), allocatable, dimension(:,:) :: density_matrix_sorted
      REAL(KIND=idp), allocatable, dimension(:) :: properties
      REAL(KIND=idp) :: sparsity

      type(CIvect), intent(in) :: pmat
      type(COOMatrix_real) :: density_matrix_coefficients_l2l2
      type(COOMatrix_integer) :: density_matrix_orbital_pairs_l2l2_1 , density_matrix_orbital_pairs_l2l2_2


      integer, dimension(tsym_max) :: nciset,nstat
      integer :: lucivec

!     The CSF files.
      type (CSFheader) :: csf_head_i,csf_head_j
      type (CSFbody)   :: csf_body_i,csf_body_j

!     The property integrals
      integer :: luprop,no_integrals,nlmq,index_for_integral_table, index_for_symmetry_block,index_for_integral_table2
      integer, allocatable, dimension(:) :: nilmq,lp,mp,qp, mob,mpob
      integer, allocatable, dimension(:,:) :: indexv, istart
      CHARACTER(LEN=8), allocatable, DIMENSION(:) :: pname
      REAL(KIND=idp), allocatable, DIMENSION(:,:) :: XBUF
      integer, allocatable, dimension(:,:) :: inverted_indexv,orbs_diag_l2l2

      write(iwrite,'(/," L^2-L^2 Block: Density Matrix Construction",/)')


!     Construct the spin orbital table
!     ------------------------------------------------------------------

      no_spin_orbitals=2*sum(csf_head_i%nob)
      allocate(list_space(no_spin_orbitals), list_ug(no_spin_orbitals), &
               list_quant_no(no_spin_orbitals), list_spin(no_spin_orbitals))

      list_space=0
      list_ug=0
      list_quant_no=0
      list_spin=0
      npflg=0
      isymtype=2

      CALL MKORBS(isymtype,csf_head_i%nob,csf_head_i%nsym, list_space,list_ug, list_quant_no, list_spin,maxso,npflg)

!     Transform the configurations in J so that they correspond to the reference determinant in I
!     subsequent routines depend on both sets of configurations being defined relative to the same
!     reference determinant.
!     ------------------------------------------------------------------

      call change_refdet(csf_head_j, csf_body_j, csf_head_i % ndtrf, no_spin_orbitals)

      n_cont_gcsf_i = inumber_of_gcsf_i - csf_head_i % l2nocsf
      n_cont_gcsf_j = inumber_of_gcsf_j - csf_head_j % l2nocsf

!     Now we identify CSF pairs that contribute to the density matrix
!     by applying the Slater-Condon rules (only those configurations that
!     differ by a single unmatched pair of spatial orbitals).
!     ------------------------------------------------------------------

!     Allocate space for the CSF density matrix (no CI vectors
!     multiplied in yet.
!     The coefficient is stored in density_matrix_coefficients_l2l2
!     The orbital pairs are stored in density_matrix_orbital_pairs_l2l2
!     in the last dimension, 1 is the Ith orbital and 2 the Jth orbital

      nestimated = csf_head_i%nocsf * csf_head_j%nocsf * sparsity
      nproc = merge(pmat%nprow * pmat%npcol, 1_blasint, pmat%CV_is_scalapack)

      call density_matrix_coefficients_l2l2    % init_matrix(csf_head_i%nocsf, csf_head_j%nocsf, nestimated/nproc)
      call density_matrix_orbital_pairs_l2l2_1 % init_matrix(csf_head_i%nocsf, csf_head_j%nocsf, nestimated/nproc)
      call density_matrix_orbital_pairs_l2l2_2 % init_matrix(csf_head_i%nocsf, csf_head_j%nocsf, nestimated/nproc)

!     MAKEMG constructs an array that indexes the reference determinant
!     array. Given an orbital number it holds the position of that
!     orbital in the reference determinant.

      allocate(ndtrf_point_back(no_spin_orbitals))

      call MAKEMG(ndtrf_point_back,no_spin_orbitals,csf_head_i%NELT,csf_head_i%ndtrf)

!     Apply Slater-Condon rules first for spatial orbitals in DRYRUN
!     then inbetween each spin orbital configuration within the CSFs
!     involved in the density matrix
!     Note: Each CSF corresponds to one spatial configuration, but
!     there are a number of ways of choosing the spin orbitals to
!     give the correct total spin quantum numbers for the configuration
      allocate(orbs_diag_l2l2(csf_head_i%NELT,csf_head_i%nocsf))
      orbs_diag_l2l2=0
      !ZM note that the parallelization scheme is different here and in the continuum_l2_drv.
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP & PRIVATE(i,j,li,lj,gi,gj,evalue,index_ndo_start_i,index_cdo_start_i,index_ndo_start_j,index_cdo_start_j) &
      !$OMP & SHARED(csf_head_i,csf_head_j,csf_body_i,csf_body_j,maxso,no_spin_orbitals,list_space,list_spin,ndtrf_point_back, &
      !$OMP &        density_matrix_coefficients_l2l2,density_matrix_orbital_pairs_l2l2_1,density_matrix_orbital_pairs_l2l2_2, &
      !$OMP &        orbs_diag_l2l2,pmat,n_cont_gcsf_i,n_cont_gcsf_j)
      !$OMP DO
      do li = 1, pmat % local_row_dimen
         do lj = 1, pmat % local_col_dimen

            ! convert local index in the prototype properties matrix to a global one
            call pmat % local_to_global(li, lj, gi, gj)
            if (gi <= n_cont_gcsf_i) cycle
            if (gj <= n_cont_gcsf_j) cycle

            ! convert the global index to a (congen) CSF index
            i = csf_head_i % last_continuum_csf + gi - n_cont_gcsf_i
            j = csf_head_j % last_continuum_csf + gj - n_cont_gcsf_j
            if (i > csf_head_i % nocsf) cycle
            if (j > csf_head_j % nocsf) cycle

            call DRYRUN(csf_body_i%indo,csf_body_j%indo,i,j,EVALUE,list_space,csf_body_i%ndo,csf_body_j%ndo,maxso)
            if (evalue) then
               ! Extract array segments corresponding to the ith and jth CSF
               index_ndo_start_i=csf_body_i%indo(i)
               index_cdo_start_i=(csf_body_i%icdo(i))
               index_ndo_start_j=csf_body_j%indo(j)
               index_cdo_start_j=(csf_body_j%icdo(j))

               call TMTMA_SPARSE(csf_head_i % nodo(i), csf_body_i % cdo(index_cdo_start_i:), csf_body_i % ndo(index_ndo_start_i:), &
                                 csf_head_j % nodo(j), csf_body_j % cdo(index_cdo_start_j:), csf_body_j % ndo(index_ndo_start_j:), &
                                 no_spin_orbitals, csf_head_i % NELT, csf_head_i % ndtrf, list_space, list_spin, &
                                 ndtrf_point_back, I, J, density_matrix_coefficients_l2l2, density_matrix_orbital_pairs_l2l2_1, &
                                 density_matrix_orbital_pairs_l2l2_2, orbs_diag_l2l2)

            end if

         end do

      end do

      !$OMP END PARALLEL

      call density_matrix_coefficients_l2l2%end_matrix_construction()
      call density_matrix_orbital_pairs_l2l2_1%end_matrix_construction()
      call density_matrix_orbital_pairs_l2l2_2%end_matrix_construction()

      deallocate(list_space,list_ug, list_quant_no, list_spin)
      deallocate(ndtrf_point_back)

      end  subroutine l2l2_drv


      subroutine continuum_l2_drv (ci_vec_target, csf_head_i, csf_head_j,csf_body_i, csf_body_j, &
                                   icontinuum_orbital_table_i, dyson_orbitals_cl2, sparsity, pmat, no_gcsf_j, cont_side)

      use blas_lapack_gbl, only: blasint
      use cdenprop_defs,   only: CSFheader, CSFbody, CIvect
      use cdenprop_aux,    only: mkorbs, change_refdet, dryrun, tmtma_sparse, makemg
      use const_gbl,       only: level2
      use mpi_gbl,         only: mpi_xermsg

      use class_COOSparseMatrix_real,    only: COOMatrix_real
      use class_COOSparseMatrix_integer, only: COOMatrix_integer

      implicit none

#if defined(mpi) && defined(scalapack)
      external dgsum2d
#endif
!     ******************************************************************
!
!     Constructs the density matrix for case where one configuration
!     contains a continuum orbital and the other doesn't.
!     Many of the routines are taken from DENPROP
!     Note: In contrast to DENPROP the density matrix between a state
!     and itself is not constructed.
!
!     ******************************************************************

      integer :: l2nocsf, lastcont, no_spin_orbitals, maxso, npflg, isymtype, itargsym, itarg, itarg_base, icsf, jcsf, iorb, &
                 icount, itarg_start_csf, no_gcsf_j, cont_side, lj, gj, no_csfs_local
      integer :: i, j, k, index_ndo_start_i, index_ndo_start_j, index_cdo_start_i, &
                 index_cdo_start_j, index_ndo_end_i, index_ndo_end_j, proc_i, proc_j, no_procs_i, no_procs_j
      integer, allocatable, dimension(:) :: list_space,list_ug, list_quant_no, list_spin,ibound_or_continuum,map_to_relative_index
      integer, allocatable, dimension(:) :: ndtrf_point_back
      integer, dimension(tsym_max) :: lucsf ! Input CSF file units
      logical :: evalue
      REAL(KIND=idp) :: sparsity
      REAL(KIND=idp), allocatable, dimension(:,:) :: density_matrix_sorted
      real(kind=idp), allocatable, DIMENSION(:,:,:) ::  dyson_orbitals_cl2
      REAL(KIND=idp), allocatable, dimension(:) :: properties
      integer, dimension(:,:) :: icontinuum_orbital_table_i
      type(CIvect), intent(in) :: pmat

      type(COOMatrix_real) :: density_matrix_coefficients_cl2
      type(COOMatrix_integer) :: density_matrix_orbital_pairs_cl2_1 , density_matrix_orbital_pairs_cl2_2

!     The  CI vectors
      type (CIvect), dimension(:) :: ci_vec_target
      integer, dimension(tsym_max) :: nciset,nstat
      integer :: lucivec
      integer, allocatable, dimension(:) :: np12target_order_map


!     The CSF files.
      type (CSFheader) :: csf_head_i,csf_head_j
      type (CSFbody)   :: csf_body_i,csf_body_j

!     The property integrals
      integer :: luprop,no_integrals,nlmq,index_for_integral_table, index_for_symmetry_block,index_for_integral_table2
      integer, allocatable, dimension(:) :: nilmq,lp,mp,qp, mob,mpob
      integer, allocatable, dimension(:,:) :: indexv, istart
      CHARACTER(LEN=8), allocatable, DIMENSION(:) :: pname
      REAL(KIND=idp), allocatable, DIMENSION(:,:) :: XBUF
      integer, allocatable, dimension(:,:) :: inverted_indexv,orbs_diag_cl2

      integer :: itarg_set,l2csf_index
      integer(blasint) :: lngth, zero = 0, one = 1
      character(len=1) :: scope
      character(len=100) :: msg

!     Local variables for Sparse matrix
      integer :: index_out, orbital_pair_tmp_1, orbital_pair_tmp_2, icsf_temp, jcsf_temp, no_elements_per_sym_cumulative, &
                 index_sparse_initial, index_sparse, nestimated
      real(kind=idp) :: value

      write(iwrite,'(/," Continuum Block-L^2: Density Matrix Construction",/)')

      ! local portion of CSFs in the J symmetry
      no_csfs_local = merge(pmat % local_col_dimen, pmat % local_row_dimen, cont_side == BraState)

      ! parallel distribution of the I and J symmetries
      proc_i     = merge(merge(pmat % myrow, pmat % mycol, cont_side == BraState), zero, pmat % CV_is_scalapack)
      proc_j     = merge(merge(pmat % mycol, pmat % myrow, cont_side == BraState), zero, pmat % CV_is_scalapack)
      no_procs_i = merge(merge(pmat % nprow, pmat % npcol, cont_side == BraState), one,  pmat % CV_is_scalapack)
      no_procs_j = merge(merge(pmat % npcol, pmat % nprow, cont_side == BraState), one,  pmat % CV_is_scalapack)

!     The first section is almost identical to l2-l2 block except that in
!     one of the wavefunctions we loop over the continuum configurations
!     (odd numbered CSFs only due to prototype CSF structure)


!     Construct the spin orbital table
!     ------------------------------------------------------------------

      no_spin_orbitals=2*sum(csf_head_i%nob)

      allocate(list_space(no_spin_orbitals), list_ug(no_spin_orbitals), &
               list_quant_no(no_spin_orbitals), list_spin(no_spin_orbitals))
      list_space=0
      list_ug=0
      list_quant_no=0
      list_spin=0
      npflg=0
      isymtype=2

      CALL MKORBS(isymtype,csf_head_i%NOB,csf_head_i%NSYM, list_space,list_ug, list_quant_no, list_spin,maxso,npflg)

!     Transform the configurations in J so that they correspond to the reference determinant in I
!     subsequent routines depend on both sets of configurations being defined relative to the same
!     reference determinant.
!     ------------------------------------------------------------------

      call change_refdet(csf_head_j, csf_body_j, csf_head_i % ndtrf, no_spin_orbitals)

!     Now we identify CSF pairs that contribute to the density matrix
!     by applying the Slater-Condon rules (only those configurations that
!     differ by a single unmatched pair of spatial orbitals).
!     ------------------------------------------------------------------

!     Allocate space for the CSF density matrix (no CI vectors
!     multiplied in yet.
!     The coefficient is stored in density_matrix_coefficients_l2l2
!     The orbital pairs are stored in density_matrix_orbital_pairs_l2l2
!     in the last dimension, 1 is the Ith orbital and 2 the Jth orbital
      call density_matrix_coefficients_cl2%set_iwrite(level2)
      call density_matrix_orbital_pairs_cl2_1%set_iwrite(level2)
      call density_matrix_orbital_pairs_cl2_2%set_iwrite(level2)

      ! based on the loop limits below and a typical worst case
      nestimated = min(csf_head_i%last_continuum_csf/2/no_procs_i * no_csfs_local * sparsity, &
                       csf_head_i % nocsf * csf_head_j % nocsf * sparsity)

      call density_matrix_coefficients_cl2    % init_matrix(csf_head_i % nocsf, csf_head_j % nocsf, nestimated)
      call density_matrix_orbital_pairs_cl2_1 % init_matrix(csf_head_i % nocsf, csf_head_j % nocsf, nestimated)
      call density_matrix_orbital_pairs_cl2_2 % init_matrix(csf_head_i % nocsf, csf_head_j % nocsf, nestimated)

!     MAKEMG constructs an array that indexes the reference determinant
!     array. Given an orbital number it holds the position of that
!     orbital in the reference determinant.

      allocate(ndtrf_point_back(no_spin_orbitals))
      call MAKEMG(ndtrf_point_back,no_spin_orbitals,csf_head_i%NELT,csf_head_i%ndtrf)

!     Apply Slater-Condon rules first for spatial orbitals in DRYRUN
!     then inbetween each spin orbital configuration within the CSFs
!     involved in the density matrix
!     Note: Each CSF corresponds to one spatial configuration, but
!     there are a number of ways of choosing the spin orbitals to
!     give the correct total spin quantum numbers for the configuration

      allocate(orbs_diag_cl2(csf_head_i%NELT,csf_head_i%nocsf))
      orbs_diag_cl2=0
      !ZM we have to parallelize the inner loop and not the outer loop (as in l2l2_drv) since the construction of dyson_orbitals_cl2 below
      !   relies on the ORDERED first index in the list of the density matrix values.
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP & PRIVATE(i,j,evalue,index_ndo_start_i,index_cdo_start_i,index_ndo_start_j,index_cdo_start_j,lj,gj) &
      !$OMP & FIRSTPRIVATE(proc_i,proc_j,no_procs_i,no_procs_j,no_csfs_local,cont_side,no_gcsf_j) &
      !$OMP & SHARED(csf_head_i,csf_head_j,csf_body_i,csf_body_j,list_space,maxso,no_spin_orbitals,list_spin,ndtrf_point_back,pmat,&
      !$OMP & density_matrix_coefficients_cl2,density_matrix_orbital_pairs_cl2_1,density_matrix_orbital_pairs_cl2_2,orbs_diag_cl2)
      do i =1, csf_head_i%last_continuum_csf,2

         ! Distribute continuum configurations over I dimension of the BLACS context. "I" is the row dimension when `cont_side` is
         ! BraState, or column dimension when `cont_side` is KetState; "J" denotes the other dimension. When distributed in
         ! this way, individual processes along the I dimension will obtain incomplete Dyson orbitals. These are all-reduced later
         ! at the end of the subroutine, so that all processes within the same J position have the same data.
         if (mod(i/2, no_procs_i) /= proc_i) cycle

         ! In the following, we loop over all L2 CSFs; however, in a distributed manner over J dimension of the BLACS contenxt.
         ! Each process loops over the appropriate local dimension (row or column, depending on `cont_side`) of the prototype
         ! properties matrix, whose row/column indices correspond to individual contracted CSFs (or "GCSFs" in CDENPROP) forming
         ! the bra/ket states. Indices are translated to CONGEN (prototype) CSF indices. Continuum CSFs are skipped.

         !$OMP DO
         do lj = 1, no_csfs_local

            ! find out the GCSF index of local element 'lj' ('j' serves as dummy here)
            if (cont_side == BraState) then
               call pmat % local_to_global(1, lj, j, gj)
            else
               call pmat % local_to_global(lj, 1, gj, j)
            end if

            ! convert to CONGEN (prototype CSF) indexing, skip continuum CSFs
            j = gj - (no_gcsf_j - csf_head_j % l2nocsf)
            if (j < 1) cycle
            j = j + csf_head_j % last_continuum_csf

            call DRYRUN(csf_body_i%indo,csf_body_j%indo,i,j,EVALUE,list_space,csf_body_i%ndo,csf_body_j%ndo,maxso)
            if (evalue) then
!              Extract array segments corresponding to the ith and jth CSF
               index_ndo_start_i=csf_body_i%indo(i)
               index_cdo_start_i=(csf_body_i%icdo(i))
               index_ndo_start_j=csf_body_j%indo(j)
               index_cdo_start_j=(csf_body_j%icdo(j))

               call TMTMA_SPARSE(csf_head_i % nodo(i), csf_body_i % cdo(index_cdo_start_i:), csf_body_i % ndo(index_ndo_start_i:), &
                                 csf_head_j % nodo(j), csf_body_j % cdo(index_cdo_start_j:), csf_body_j % ndo(index_ndo_start_j:), &
                                 no_spin_orbitals, csf_head_i % NELT, csf_head_i % ndtrf, list_space, list_spin, &
                                 ndtrf_point_back, I, J, density_matrix_coefficients_cl2, density_matrix_orbital_pairs_cl2_1, &
                                 density_matrix_orbital_pairs_cl2_2, orbs_diag_cl2)

            end if

         end do
         !$OMP END DO

      end do
      !$OMP END PARALLEL

      deallocate(orbs_diag_cl2)

      call density_matrix_coefficients_cl2%end_matrix_construction()
      call density_matrix_orbital_pairs_cl2_1%end_matrix_construction()
      call density_matrix_orbital_pairs_cl2_2%end_matrix_construction()

!     Now we make dyson-like orbitals between target states and each L^2 CSF
!     They are stored in dyson_orbitals_cl2. For each target-L^2 csf pair
!     the set of bound orbitals are of a single symmetry given by the target
!     symmetry X the L^2 function symmetry. A coefficient for each orbital
!     is constructed out of the spin-coupling coefficients, the  target CI
!     coefficients and the phase vector iphz.
!     ------------------------------------------------------------------

!     need function that returns the bound orbital of a pair (add array to orbital table)
!     give orbitals a relative index within each symmetry

      allocate(ibound_or_continuum(no_spin_orbitals/2), map_to_relative_index(no_spin_orbitals/2))
      ibound_or_continuum=0
      map_to_relative_index=0

      icount=1
      do i =1,csf_head_i%nsym
!        Bound orbitals
         do j=1, csf_head_i%nob0(i)+csf_head_i%no_l2_virtuals(i)
            ibound_or_continuum(icount)=1
            icount=icount+1

         end do
!        Continuum or virtual orbitals
         do j=1, csf_head_i%nob(i)-csf_head_i%nob0(i)-csf_head_i%no_l2_virtuals(i)
             ibound_or_continuum(icount)=0
             icount=icount+1

         end do

      end do

!     Construct array giving relative index within symmetry
      icount=1
      do i =1,csf_head_i%nsym
         do j=1, csf_head_i%nob(i)
            map_to_relative_index(icount)=j
            icount=icount+1

         end do

      end do

!     Construct the N+1 congen order to target mapping
!     Could potentially be removed as target vectors are now read in
!     in the N+1 congen target expansion order.

      allocate(np12target_order_map(csf_head_i%ntgsym))
      np12target_order_map=0

      icount=0
      do i=1,csf_head_i%ntgsym
          do j=1,csf_head_i%ntgsym
             if (((ci_vec_target(i) % mgvn + 1) == csf_head_i % itarget_symmetry_order(j,1)) .and. &
                 ((2 * ci_vec_target(i) % S + 1) == real(csf_head_i % itarget_symmetry_order(j,2)))) then
                 np12target_order_map(j) = i
                 icount=icount+1

             end if

          end do

      end do

      if (gflag == 0) then
          icount=0
          do i=1,csf_head_i%ntgsym
              np12target_order_map(i) = i
              icount=icount+1
          end do
      endif
!--------------------------------

      if (icount .ne. csf_head_i%ntgsym) then
         write (msg, '(a,2(1x,i0))') 'Target symmetry mismatch between target and N+1 CI vectors', icount, csf_head_i%ntgsym
         call mpi_xermsg('cdenprop_procs', 'continuum_l2_drv', msg, 1, 1)
      end if

      allocate(dyson_orbitals_cl2(maxval(csf_head_i%nob0+csf_head_i%no_l2_virtuals), no_csfs_local, sum(csf_head_i%numtgt)))
      dyson_orbitals_cl2=0.0_idp

!     loop over L^2 csfs first

      no_elements_per_sym_cumulative=0
      index_sparse_initial=1
      itarg_base=0
      index_sparse=1
      do itargsym=1,csf_head_i%ntgsym
         no_elements_per_sym_cumulative =no_elements_per_sym_cumulative+csf_head_i%nctgt(itargsym)*2

         index_sparse= index_sparse_initial
         if(index_sparse .le. density_matrix_orbital_pairs_cl2_2%nnz) then
            call density_matrix_orbital_pairs_cl2_2%get_indexed_element(index_sparse,icsf_temp,jcsf_temp,orbital_pair_tmp_2)
         else
            exit
         end if

         do while (icsf_temp / no_elements_per_sym_cumulative == 0 .and. index_sparse <= density_matrix_orbital_pairs_cl2_2 % nnz)

            call density_matrix_orbital_pairs_cl2_1%get_indexed_element(index_sparse,icsf_temp,jcsf_temp,orbital_pair_tmp_1)
            call density_matrix_orbital_pairs_cl2_2%get_indexed_element(index_sparse,icsf_temp,jcsf_temp,orbital_pair_tmp_2)
            call density_matrix_coefficients_cl2%get_indexed_element(index_sparse,icsf_temp,jcsf_temp,value)

            ! Pick out the bound orbital from the orbital pair

            if ((orbital_pair_tmp_1 .ne. 0) .and. (orbital_pair_tmp_2 .ne. 0)) then
               if( ibound_or_continuum(orbital_pair_tmp_1) .eq. 1) then
                  iorb=orbital_pair_tmp_1 !density_matrix_orbital_pairs_cl2(icsf,jcsf,1)
               else if ( ibound_or_continuum(orbital_pair_tmp_2) .eq. 1) then
                  iorb=orbital_pair_tmp_2 !density_matrix_orbital_pairs_cl2(icsf,jcsf,2)
               else
                  write (msg, '(a,2(1x,i0),a,2(1x,i0))') 'No bound orbital in orbital pair', &
                     orbital_pair_tmp_1, orbital_pair_tmp_2, ' for CSFs', icsf_temp, jcsf_temp
                  call mpi_xermsg('cdenprop_procs', 'continuum_l2_drv', msg, 1, 1)
               end if

               iorb=map_to_relative_index(iorb)

               itarg_set= np12target_order_map(itargsym)

               ! covert CONGEN (contracted) CSF index to (global) generalized CSF index
               i = (jcsf_temp - csf_head_j%last_continuum_csf) + (no_gcsf_j - csf_head_j%l2nocsf)

               ! convert global index of GCSF to local index
               if (cont_side == BraState) then
                  call pmat % global_to_local(1, i, j, l2csf_index)  ! j is dummy here
               else
                  call pmat % global_to_local(i, 1, l2csf_index, j)  ! j is dummy here
               end if

!              multiply in spin-coupling, target vector and iphz coefficents store sum in dyson orbital vector
!              use relative orbital indices for the bound orbitals.

               j=(icsf_temp+1)/2 - (no_elements_per_sym_cumulative/2-csf_head_i%nctgt(itargsym)) !+ 1

               itarg=itarg_base
               do i=1,csf_head_i%numtgt(itargsym)
                  itarg=itarg+1
                  !ZM the iphz phase is applied to the configuration containing the continuum orbital. We're performing the target contraction here:
                  !the target configurations come from target congen run
                  !which use a different ordering of the orbitals
                  !resulting in phases of the target configurations wrt configurations
                  !generated in the scattering run. The iphz phase is applied
                  !to the N+1 configuration to eliminate this phase difference.
                  !where the target configuration is generated in the target congen run and therefore is not necessarily in the same dictionary order.
                  dyson_orbitals_cl2(iorb,l2csf_index,itarg) = &
                  dyson_orbitals_cl2(iorb,l2csf_index,itarg) + value &
                    * ci_vec_target(itarg_set) % CV(j,i) &
                    * csf_head_i % iphz((icsf_temp+1)/2) &
                    * csf_head_i % itarget_overall_phase(itarg)
               enddo

            else
               write(iwrite,*) 'Non zero coeff but zero orbitals, why?'
               write(iwrite,'(2i4,d20.5)') orbital_pair_tmp_1,orbital_pair_tmp_2,value
            end if

            index_sparse =index_sparse+1
            if(index_sparse .le. density_matrix_orbital_pairs_cl2_2%nnz) then
               call density_matrix_orbital_pairs_cl2_2%get_indexed_element(index_sparse,icsf_temp,jcsf_temp,orbital_pair_tmp_2)
            else
               exit
            end if

         end do !while
         index_sparse_initial=index_sparse
         itarg_base = itarg_base+csf_head_i%numtgt(itargsym)
      end do !itargsym

!       if (icsf .ne. csf_head_i%last_continuum_csf+1) then
!          stop ' ERROR: Target configurations selection problem'
!       end if

#if defined(usempi) && defined(scalapack)
      ! now combine the partial results obtained by individual members of this BLACS context
      if (pmat % CV_is_scalapack) then
         lngth = size(dyson_orbitals_cl2)
         scope = merge('C', 'R', cont_side == BraState)
         call dgsum2d(pmat % blacs_context, scope, ' ', lngth, one, dyson_orbitals_cl2, lngth, -one, -one)
      end if
#endif

      end subroutine continuum_l2_drv



      function get_property_integral(pintegrals,iprop,iorb1,jorb1)
      use cdenprop_defs
      use ukrmol_interface_gbl, only: get_ukrmolp_property_integral => get_property_integral
      use mpi_gbl, only: mpi_xermsg
      implicit none

!     Argument variables
      type (property_integrals) :: pintegrals
      integer :: iorb1,jorb1,itemp_orb, iprop
      real(kind=idp) :: get_property_integral

!     Local variables
      character(len=200) :: msg
      integer ::  index_for_integral_table,index_for_integral_table2, index_for_symmetry_block, iorb,jorb
      integer :: i,j ! TESTING ONLY to be deleted

      get_property_integral=0.0_idp

      if (iorb1 < 1 .or. jorb1 < 1) then
         write (msg, '(a,i0,a,i0,a)') 'Invalid orbital indices ', iorb1, ' ', jorb1, &
                                      ' on entry. Are there duplicate CSFs in the CONGEN file?'
         call mpi_xermsg('cdenprop_procs', 'get_property_integral', trim(msg), 1, 1)
      end if

      if (pintegrals % ukrmolp_ints) then
         get_property_integral = get_ukrmolp_property_integral(iprop,iorb1,jorb1)
         return
      end if

!     Property integrals are stored upper triangle by symmetry block. The following condition
!     enforces this for orbital pairs.
      if((pintegrals%mob(iorb1) .gt. pintegrals%mob(jorb1))) then
         iorb=jorb1
         jorb=iorb1

      else
         iorb=iorb1
         jorb=jorb1

      end if

      index_for_integral_table=0
      index_for_symmetry_block= pintegrals%istart(pintegrals%mob(iorb)+1,pintegrals%mob(jorb)+1)
! write(iwrite,*) 'pintegrals%istart'
! do i=1,8
! write(iwrite,'(20i6)')  (pintegrals%istart(i,j),j=1,8)
! end do
! write(iwrite,*) iprop, index_for_symmetry_block
!     Mapping between density matrix index and indexing scheme of the integral table
      if ((pintegrals%mob(iorb)+1) .eq. (pintegrals%mob(jorb)+1)) then
        ! Lower triangle rowwise=upper triangle columnwise
         index_for_integral_table = index_for_symmetry_block + &
            ij_to_index_lower_rowwise(pintegrals % mpob(jorb), pintegrals % mpob(iorb)) - 1

! print *, "orig", index_for_integral_table
! index_for_integral_table = index_for_symmetry_block +ij_to_index_full_colwise(pintegrals%mpob(jorb),pintegrals%mpob(iorb),pintegrals%nob(pintegrals%mob(iorb)+1) )-1
! print *, "new", index_for_integral_table

         if (pintegrals%inverted_indexv(index_for_integral_table,iprop) .ne. 0) then

            index_for_integral_table2 = pintegrals%inverted_indexv(index_for_integral_table,iprop)
            get_property_integral = pintegrals%xintegrals(index_for_integral_table2,iprop)
! write(iwrite,'(3i5,d20.5)') index_for_integral_table2,iorb,jorb,get_property_integral
         end if

      else if ((pintegrals%mob(iorb)+1) .lt. (pintegrals%mob(jorb)+1)) then
         index_for_integral_table = index_for_symmetry_block + &
            ij_to_index_full_colwise(pintegrals % mpob(iorb), &
                                     pintegrals % mpob(jorb), &
                                     pintegrals % nob(pintegrals % mob(iorb) + 1)) - 1

         if (pintegrals%inverted_indexv(index_for_integral_table,iprop) .ne. 0) then

            index_for_integral_table2 = pintegrals%inverted_indexv(index_for_integral_table,iprop)
            get_property_integral = pintegrals%xintegrals(index_for_integral_table2,iprop)

         end if

      else
         get_property_integral = 0.0_idp

      end if

      return
      end function get_property_integral


    !> \brief    (Re)allocate the prototype property matrix
    !> \authors  A Harvey, Z Masin, J Benda
    !> \date     2011 - 2019
    !>
    !> Check if the prototype property matrix is allocated and has the correct size for the current
    !> combination of bra and ket CSFs. If not, reallocate it.
    !>
    subroutine allocate_prototype_property_matrix (pmat, m, n, iwrite)

        use cdenprop_defs, only: CIvect

        type(CIvect), intent(inout) :: pmat
        integer,      intent(in)    :: m, n, iwrite

        ! allocate only if necessary: for large problems this allocation takes a while
        if (allocated(pmat % CV)) then
            if (pmat % mat_dimen_r /= m .or. pmat % mat_dimen_c /= n) then
                call pmat % final_CV
            end if
        end if
        if (.not. allocated(pmat % CV)) then
            write (iwrite, '(/," ALLOCATING prototype_properties_matrix",2I10)') m, n
            call pmat % init_CV(m, n)
        end if

    end subroutine allocate_prototype_property_matrix


      subroutine create_prototype_properties_matrix (csf_head_i, csf_head_j, icontinuum_orbital_table_i, &
            icontinuum_orbital_table_j, target_properties, dyson_orbitals_cl2, dyson_orbitals_l2c, &
            density_matrix_orbital_pairs_l2l2_1, density_matrix_orbital_pairs_l2l2_2, density_matrix_coefficients_l2l2, &
            orbs_diag_l2l2, pintegrals, iprop, icontains_continuum_orbitals, prototype_properties_matrix, ukrmolp_ints)
      use cdenprop_defs
      use class_COOSparseMatrix_integer
      use class_COOSparseMatrix_real
      implicit none

!     Argument variables
      type (CSFheader) :: csf_head_i,csf_head_j
      integer,        allocatable, dimension(:,:)   :: icontinuum_orbital_table_i, icontinuum_orbital_table_j
      real(kind=idp), allocatable, dimension(:,:,:) :: dyson_orbitals_cl2, dyson_orbitals_l2c, target_properties
      type (property_integrals) :: pintegrals
      integer :: iprop,icontains_continuum_orbitals
      type (CIvect) :: prototype_properties_matrix
      integer, dimension(csf_head_i%NELT,*):: orbs_diag_l2l2

      type(COOMatrix_real) :: density_matrix_coefficients_l2l2
      type(COOMatrix_integer) :: density_matrix_orbital_pairs_l2l2_1,density_matrix_orbital_pairs_l2l2_2
      logical :: ukrmolp_ints
!     Local variables
      integer :: igcsf,jgcsf,icsf,jcsf, number_of_gcsf_i, number_of_gcsf_j, no_cont_gcsf_i, no_cont_gcsf_j,&
     &           iorb,jorb, idyson_orb,jdyson_orb, max_dyson_orb, i,j,k, icount, itarget_index, jtarget_index,&
     &           last_continuum_gcsf_i, last_continuum_gcsf_j, itarg, jtarg, err, li, lj
      integer, allocatable, dimension(:) :: idyson_orbital_symmetry
      integer, allocatable, dimension(:,:) :: map_relative_to_absolute_index

      real(kind=idp) :: value, element

      select case(icontains_continuum_orbitals)
      case(0)
         no_cont_gcsf_i = 0
         no_cont_gcsf_j = 0

      case(1)
         no_cont_gcsf_i = size(icontinuum_orbital_table_i,1)
         no_cont_gcsf_j = 0

      case(2)
         no_cont_gcsf_i = 0
         no_cont_gcsf_j = size(icontinuum_orbital_table_j,1)

      case(3)
         no_cont_gcsf_i = size(icontinuum_orbital_table_i,1)
         no_cont_gcsf_j = size(icontinuum_orbital_table_j,1)

      end select

      number_of_gcsf_i = no_cont_gcsf_i + csf_head_i%l2nocsf
      number_of_gcsf_j = no_cont_gcsf_j + csf_head_j%l2nocsf

      prototype_properties_matrix % CV = 0

!     Diagonal Blocks
!     ------------------------------------------------------------------

!     L2-L2 Block
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP & PRIVATE(i,igcsf,jgcsf,value,icsf,jcsf,iorb,jorb,element) &
      !$OMP & SHARED(density_matrix_coefficients_l2l2,density_matrix_orbital_pairs_l2l2_1,density_matrix_orbital_pairs_l2l2_2,&
      !$OMP &        csf_head_i,csf_head_j,prototype_properties_matrix,pintegrals,iprop,no_cont_gcsf_i,no_cont_gcsf_j)
      !$OMP DO
      do i = 1,density_matrix_coefficients_l2l2%nnz
         call density_matrix_coefficients_l2l2%get_indexed_element(i,icsf,jcsf,value)

         igcsf=icsf-csf_head_i%last_continuum_csf+no_cont_gcsf_i
         jgcsf=jcsf-csf_head_j%last_continuum_csf+no_cont_gcsf_j

         call density_matrix_orbital_pairs_l2l2_1%get_indexed_element(i,icsf,jcsf,iorb)
         call density_matrix_orbital_pairs_l2l2_2%get_indexed_element(i,icsf,jcsf,jorb)

         element = value*get_property_integral(pintegrals,iprop,iorb,jorb)
!         prototype_properties_matrix%CV(igcsf,jgcsf)= value*get_property_integral(pintegrals,iprop,iorb,jorb)
         call prototype_properties_matrix%set_CV_element(element,igcsf,jgcsf)

      end do
      !$OMP END PARALLEL


      do iorb=csf_head_i%last_continuum_csf+1,csf_head_i%last_continuum_csf+csf_head_i%l2nocsf
         igcsf=iorb-csf_head_i%last_continuum_csf+no_cont_gcsf_i
         do jorb=1,csf_head_i%NELT
            if (orbs_diag_l2l2(jorb,iorb).ne.0) then
               element = get_property_integral(pintegrals,iprop,orbs_diag_l2l2(jorb,iorb),orbs_diag_l2l2(jorb,iorb))
               call prototype_properties_matrix%add_to_CV_element(element,igcsf,igcsf)
!                 prototype_properties_matrix%CV(igcsf,igcsf)=prototype_properties_matrix%CV(igcsf,igcsf)+&
!     &              get_property_integral(pintegrals,iprop,orbs_diag_l2l2(jorb,iorb),orbs_diag_l2l2(jorb,iorb))

            end if

         end do

      end do

      select case(icontains_continuum_orbitals)
      case(1)
!     Off-diagonal Blocks
!     ------------------------------------------------------------------

         if (allocated(dyson_orbitals_cl2)) then
            max_dyson_orb = size(dyson_orbitals_cl2, 1)
         else
            max_dyson_orb = 0
         end if

!     determine the dyson orbital symmetry for each target state
!     target symmetry X L^2 symmetry

         allocate(idyson_orbital_symmetry(sum(csf_head_i%numtgt)))
         idyson_orbital_symmetry=0

         k=0
         do i=1,size(csf_head_i%numtgt)
            do j=1,csf_head_i%numtgt(i)
               k=k+1
               idyson_orbital_symmetry(k)=IPD2H(csf_head_i%itarget_symmetry_order(i,1),csf_head_j%mgvn+1)

            end do

         end do

!     we need the mapping from relative orbital index to absolute index
         allocate(map_relative_to_absolute_index(maxval(csf_head_i%nob0+csf_head_i%no_l2_virtuals),csf_head_i%nsym))
         map_relative_to_absolute_index=0

         icount=1
         do i =1,csf_head_i%nsym
!        Bound orbitals
            do j=1, csf_head_i%nob0(i)+csf_head_i%no_l2_virtuals(i)
               map_relative_to_absolute_index(j,i)=icount
               icount=icount+1

            end do
!        Skip continuum and virtual orbitals
            do j=1, csf_head_i%nob(i)-csf_head_i%nob0(i)-csf_head_i%no_l2_virtuals(i)
               icount=icount+1

            end do

         end do


!     Continuum-L2 block
         !$OMP PARALLEL DEFAULT(SHARED) &
         !$OMP & PRIVATE(igcsf,itarget_index,jgcsf,jcsf,iorb,idyson_orb,jorb,element,lj) &
         !$OMP & SHARED(no_cont_gcsf_i,icontinuum_orbital_table_i,no_cont_gcsf_j,number_of_gcsf_j,max_dyson_orb, &
         !$OMP &        map_relative_to_absolute_index,idyson_orbital_symmetry,prototype_properties_matrix, &
         !$OMP &        dyson_orbitals_cl2,pintegrals,iprop)
         !$OMP DO
         do igcsf=1,no_cont_gcsf_i
            itarget_index=icontinuum_orbital_table_i(igcsf,1)
            iorb=icontinuum_orbital_table_i(igcsf,5)
            do lj = 1, prototype_properties_matrix % local_col_dimen
               call prototype_properties_matrix % local_to_global(1, lj, jcsf, jgcsf)  ! jcsf is a dummy here
               if (jgcsf <= no_cont_gcsf_j .or. number_of_gcsf_j < jgcsf) cycle  ! only loop over L2 section of J
               do idyson_orb=1,max_dyson_orb
                  jorb=map_relative_to_absolute_index(idyson_orb,idyson_orbital_symmetry(itarget_index))
                  if (iorb /= 0 .and. jorb /= 0 .and. dyson_orbitals_cl2(idyson_orb,lj,itarget_index) /= 0) then
                     element = dyson_orbitals_cl2(idyson_orb,lj,itarget_index) * get_property_integral(pintegrals,iprop,iorb,jorb)
                     call prototype_properties_matrix%add_to_CV_element(element,igcsf,jgcsf)
                  end if
               end do
            end do
         end do
         !$OMP END PARALLEL

      case(2)
!     Off-diagonal Blocks
!     ------------------------------------------------------------------

         if (allocated(dyson_orbitals_l2c)) then
            max_dyson_orb = size(dyson_orbitals_l2c, 1)
         else
            max_dyson_orb = 0
         end if

!     determine the dyson orbital symmetry for each target state
!     target symmetry X L^2 symmetry

         allocate(idyson_orbital_symmetry(sum(csf_head_j%numtgt)))
         idyson_orbital_symmetry=0

         k=0
         do i=1,size(csf_head_j%numtgt)
            do j=1,csf_head_j%numtgt(i)
               k=k+1
               idyson_orbital_symmetry(k)=IPD2H(csf_head_j%itarget_symmetry_order(i,1),csf_head_i%mgvn+1)

            end do

         end do

!     we need the mapping from relative orbital index to absolute index
         allocate(map_relative_to_absolute_index(maxval(csf_head_j%nob0+csf_head_j%no_l2_virtuals),csf_head_j%nsym))
         map_relative_to_absolute_index=0

         icount=1
         do i =1,csf_head_j%nsym
!        Bound orbitals
            do j=1, csf_head_j%nob0(i)+csf_head_j%no_l2_virtuals(i)
               map_relative_to_absolute_index(j,i)=icount
               icount=icount+1

            end do
!        Skip continuum and virtual orbitals
            do j=1, csf_head_j%nob(i)-csf_head_j%nob0(i)-csf_head_j%no_l2_virtuals(i)
               icount=icount+1

            end do

         end do


!     Continuum-L2 block

         !$OMP PARALLEL DEFAULT(SHARED) &
         !$OMP & PRIVATE(jgcsf,jtarget_index,jorb,igcsf,icsf,jdyson_orb,iorb,element,li) &
         !$OMP & SHARED(no_cont_gcsf_j,icontinuum_orbital_table_j,no_cont_gcsf_i,number_of_gcsf_i,max_dyson_orb, &
         !$OMP &        map_relative_to_absolute_index,idyson_orbital_symmetry,dyson_orbitals_l2c, &
         !$OMP &        prototype_properties_matrix,pintegrals,iprop)
         !$OMP DO
         do jgcsf=1,no_cont_gcsf_j ! loop over continuum part of J
            jtarget_index=icontinuum_orbital_table_j(jgcsf,1)
            jorb=icontinuum_orbital_table_j(jgcsf,5)
            do li = 1, prototype_properties_matrix % local_row_dimen
               call prototype_properties_matrix % local_to_global(li, 1, igcsf, icsf)  ! icsf is a dummy here
               if (igcsf <= no_cont_gcsf_i .or. number_of_gcsf_i < igcsf) cycle  ! only loop over L2 section of I
               do jdyson_orb=1,max_dyson_orb
                  iorb=map_relative_to_absolute_index(jdyson_orb,idyson_orbital_symmetry(jtarget_index))
                  if (iorb /= 0 .and. jorb /= 0 .and. dyson_orbitals_l2c(jdyson_orb,li,jtarget_index) /= 0) then
                     element = dyson_orbitals_l2c(jdyson_orb,li,jtarget_index) * get_property_integral(pintegrals,iprop,iorb,jorb)
                     call prototype_properties_matrix%add_to_CV_element(element,igcsf,jgcsf)
                  end if
               end do
            end do
         end do
         !$OMP END PARALLEL

      case(3)
!     Continuum-Continuum Block
         if (iprop .eq. 1) then
!        Multiply by number of electrons for the overlaps
            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP & PRIVATE(igcsf,jgcsf,iorb,jorb,element) &
            !$OMP & SHARED(no_cont_gcsf_i,no_cont_gcsf_j,icontinuum_orbital_table_i,icontinuum_orbital_table_j, &
            !$OMP & prototype_properties_matrix,pintegrals,iprop,csf_head_i)
            !$OMP DO
            do igcsf=1,no_cont_gcsf_i
               do jgcsf=1,no_cont_gcsf_j
                  if (icontinuum_orbital_table_i(igcsf,1) .eq. icontinuum_orbital_table_j(jgcsf,1)) then
                     iorb=icontinuum_orbital_table_i(igcsf,5)
                     jorb=icontinuum_orbital_table_j(jgcsf,5)
                     element = get_property_integral(pintegrals,iprop,iorb,jorb)*(csf_head_i%nelt)
                     call prototype_properties_matrix%set_CV_element(element,igcsf,jgcsf)
!                     prototype_properties_matrix%CV(igcsf,jgcsf)= get_property_integral(pintegrals,iprop,iorb,jorb)*(csf_head_i%nelt)
                  end if

               end do

            end do
            !$OMP END PARALLEL

          else
             !$OMP PARALLEL DEFAULT(SHARED) &
             !$OMP & PRIVATE(igcsf,jgcsf,itarg,jtarg,iorb,jorb,element) &
             !$OMP & SHARED(no_cont_gcsf_i,no_cont_gcsf_j,icontinuum_orbital_table_i,icontinuum_orbital_table_j,&
             !$OMP &        prototype_properties_matrix,pintegrals,iprop,target_properties,csf_head_i,csf_head_j,ukrmolp_ints)
             !$OMP DO
             do igcsf=1,no_cont_gcsf_i
                do jgcsf=1,no_cont_gcsf_j

                   itarg=icontinuum_orbital_table_i(igcsf,1)
                   jtarg= icontinuum_orbital_table_j(jgcsf,1)
                   if (itarg .eq. jtarg) then
                      iorb=icontinuum_orbital_table_i(igcsf,5)
                      jorb=icontinuum_orbital_table_j(jgcsf,5)
                      element = get_property_integral(pintegrals,iprop,iorb,jorb) &
                              * csf_head_i % itarget_overall_phase(itarg) &
                              * csf_head_j % itarget_overall_phase(jtarg)
                      call prototype_properties_matrix%set_CV_element(element,igcsf,jgcsf)
!                      prototype_properties_matrix%CV(igcsf,jgcsf)= get_property_integral(pintegrals,iprop,iorb,jorb) *csf_head_i%itarget_overall_phase(itarg)*csf_head_j%itarget_overall_phase(jtarg)

                   end if
                   !target states not equal and continuum orbitals equal, get relevant target state dipole
                   if (icontinuum_orbital_table_i(igcsf,4) == icontinuum_orbital_table_j(jgcsf,4) .and. &
                       icontinuum_orbital_table_i(igcsf,5) == icontinuum_orbital_table_j(jgcsf,5))  then
                      element = target_properties(itarg, jtarg, dipole_component_map(iprop,ukrmolp_ints)) &
                              * csf_head_i % itarget_overall_phase(itarg) &
                              * csf_head_j % itarget_overall_phase(jtarg)
                      call prototype_properties_matrix%add_to_CV_element(element,igcsf,jgcsf)
!                       prototype_properties_matrix%CV(igcsf,jgcsf)=prototype_properties_matrix%CV(igcsf,jgcsf)+target_properties(itarg,jtarg,dipole_component_map(iprop,ukrmolp_ints)) *csf_head_i%itarget_overall_phase(itarg)*csf_head_j%itarget_overall_phase(jtarg)

                   end if

                end do

             end do
             !$OMP END PARALLEL

          end if
!
!     Off-diagonal Blocks
!     ------------------------------------------------------------------

         if (allocated(dyson_orbitals_cl2)) then
            max_dyson_orb = size(dyson_orbitals_cl2, 1)
         else
            max_dyson_orb = 0
         end if

!     determine the dyson orbital symmetry for each target state
!     target symmetry X L^2 symmetry

         allocate(idyson_orbital_symmetry(sum(csf_head_i%numtgt)))
         idyson_orbital_symmetry=0

         k=0
         do i=1,size(csf_head_i%numtgt)
            do j=1,csf_head_i%numtgt(i)
               k=k+1
               idyson_orbital_symmetry(k)=IPD2H(csf_head_i%itarget_symmetry_order(i,1),csf_head_j%mgvn+1)

            end do

         end do

!     we need the mapping from relative orbital index to absolute index
         allocate(map_relative_to_absolute_index(maxval(csf_head_i%nob0+csf_head_i%no_l2_virtuals),csf_head_i%nsym))
         map_relative_to_absolute_index=0

         icount=1
         do i =1,csf_head_i%nsym
!        Bound orbitals
            do j=1, csf_head_i%nob0(i)+csf_head_i%no_l2_virtuals(i)
               map_relative_to_absolute_index(j,i)=icount
               icount=icount+1

            end do
!        Skip continuum and virtual orbitals
            do j=1, csf_head_i%nob(i)-csf_head_i%nob0(i)-csf_head_i%no_l2_virtuals(i)
               icount=icount+1

            end do

         end do


!     Continuum-L2 block
         !$OMP PARALLEL DEFAULT(SHARED) &
         !$OMP & PRIVATE(igcsf,itarget_index,jgcsf,jcsf,iorb,idyson_orb,jorb,element,lj) &
         !$OMP & SHARED(no_cont_gcsf_i,no_cont_gcsf_j,number_of_gcsf_j,icontinuum_orbital_table_i,max_dyson_orb, &
         !$OMP &        map_relative_to_absolute_index,idyson_orbital_symmetry,prototype_properties_matrix, &
         !$OMP &        dyson_orbitals_cl2,pintegrals,iprop)
         !$OMP DO
         do igcsf=1,no_cont_gcsf_i
            itarget_index=icontinuum_orbital_table_i(igcsf,1)
            iorb=icontinuum_orbital_table_i(igcsf,5)
            do lj = 1, prototype_properties_matrix % local_col_dimen
               call prototype_properties_matrix % local_to_global(1, lj, jcsf, jgcsf)  ! jcsf is a dummy here
               if (jgcsf <= no_cont_gcsf_j .or. number_of_gcsf_j < jgcsf) cycle  ! only loop over L2 section of J
               do idyson_orb=1,max_dyson_orb
                  jorb=map_relative_to_absolute_index(idyson_orb,idyson_orbital_symmetry(itarget_index))
                  if (iorb /= 0 .and. jorb /= 0 .and. dyson_orbitals_cl2(idyson_orb,lj,itarget_index) /= 0) then
                     element = dyson_orbitals_cl2(idyson_orb,lj,itarget_index) * get_property_integral(pintegrals,iprop,iorb,jorb)
                     call prototype_properties_matrix%add_to_CV_element(element,igcsf,jgcsf)
                  end if
               end do
            end do
         end do
         !$OMP END PARALLEL
!~ end if
!       write(iwrite,*) 'After C-L2 Block'

!     L2-Continuum Block
         if (allocated(dyson_orbitals_l2c)) then
            max_dyson_orb = size(dyson_orbitals_l2c, 1)
         else
            max_dyson_orb = 0
         end if
!     determine the dyson orbital symmetry for each target state
!     target symmetry X L^2 symmetry
         idyson_orbital_symmetry=0

         k=0
         do i=1,size(csf_head_j%numtgt)
            do j=1,csf_head_j%numtgt(i)
               k=k+1
               idyson_orbital_symmetry(k)=IPD2H(csf_head_j%itarget_symmetry_order(i,1),csf_head_i%mgvn+1)

            end do

         end do

         !$OMP PARALLEL DEFAULT(SHARED) &
         !$OMP & PRIVATE(jgcsf,jtarget_index,jorb,igcsf,icsf,jdyson_orb,iorb,element,li) &
         !$OMP & SHARED(no_cont_gcsf_j,icontinuum_orbital_table_j,no_cont_gcsf_i,number_of_gcsf_i,max_dyson_orb,&
         !$OMP &        map_relative_to_absolute_index,idyson_orbital_symmetry,dyson_orbitals_l2c,prototype_properties_matrix, &
         !$OMP &        pintegrals,iprop)
         !$OMP DO
         do jgcsf=1,no_cont_gcsf_j ! loop over continuum part of J
            jtarget_index=icontinuum_orbital_table_j(jgcsf,1)
            jorb=icontinuum_orbital_table_j(jgcsf,5)
            do li = 1, prototype_properties_matrix % local_row_dimen
               call prototype_properties_matrix % local_to_global(li, 1, igcsf, icsf)  ! icsf is a dummy here
               if (igcsf <= no_cont_gcsf_i .or. number_of_gcsf_i < igcsf) cycle  ! only loop over L2 section of I
               do jdyson_orb=1,max_dyson_orb
                  iorb=map_relative_to_absolute_index(jdyson_orb,idyson_orbital_symmetry(jtarget_index))
                  if (iorb /= 0 .and. jorb /= 0 .and. dyson_orbitals_l2c(jdyson_orb,li,jtarget_index) /= 0) then
                     element = dyson_orbitals_l2c(jdyson_orb,li,jtarget_index) * get_property_integral(pintegrals,iprop,iorb,jorb)
                     call prototype_properties_matrix%add_to_CV_element(element,igcsf,jgcsf)
                  end if
               end do
            end do
         end do
         !$OMP END PARALLEL

      end select

      end subroutine  create_prototype_properties_matrix

      real(kind=idp) function multipole_nuclear_term(iprop,ukrmolp_ints,x,y,z)
      implicit none
      real(kind=idp) :: x,y,z
      integer :: iprop
      logical :: ukrmolp_ints
      if (ukrmolp_ints) then

         select case(iprop)
         ! dipoles
         case(2)
            multipole_nuclear_term=y

         case(3)
            multipole_nuclear_term=z

         case(4)
            multipole_nuclear_term=x

         case(5)
            multipole_nuclear_term=sqrt(3._idp)*x*y

         case(6)
            multipole_nuclear_term=(sqrt(3._idp)*y*z)

         case(7)
            multipole_nuclear_term=(0.5_idp*(2._idp*z*z-x*x-y*y))


         case(8)
            multipole_nuclear_term=(sqrt(3._idp)*x*z)

         case(9)
            multipole_nuclear_term=(0.5_idp*sqrt(3._idp)*(x*x-y*y))

         end select

      else
         select case(iprop)
         ! dipoles
         case(2)
            multipole_nuclear_term=(x)

         case(3)
            multipole_nuclear_term=(y)

         case(4)
            multipole_nuclear_term=(z)

         case(5)
            multipole_nuclear_term=(0.5_idp*(2._idp*z*z-x*x-y*y))

         case(6)
            multipole_nuclear_term=(sqrt(3._idp)*x*y)

         case(7)
            multipole_nuclear_term=(sqrt(3._idp)*x*z)


         case(8)
            multipole_nuclear_term=(0.5_idp*sqrt(3._idp)*(x*x-y*y))

         case(9)
            multipole_nuclear_term=(sqrt(3._idp)*y*z)

         end select
      end if

      end function



      integer function dipole_component_map(iprop,ukrmolp_ints)
      implicit none
      integer :: iprop
      logical :: ukrmolp_ints
! Overlap
      dipole_component_map=0
      if (ukrmolp_ints) then
      ! The property integrals are stored in the same order as the target properties are.
      ! i.e. (00, 1-1, 10, 11 ... )
      !      (olap, d_y, d_z, d_x, ... )
            dipole_component_map=iprop

      else
      ! The property integrals are stored in a different order to the target properties.
      ! i.e. (00,   11,  1-1, 10 ...)
      !      (olap, d_x, d_y, d_z, ...
         select case(iprop)
         ! dipoles
         case(2)
            dipole_component_map=4

         case(3)
            dipole_component_map=2

         case(4)
            dipole_component_map=3

         case(5)
            dipole_component_map=7

         case(6)
            dipole_component_map=5

         case(7)
            dipole_component_map=8

         case(8)
            dipole_component_map=9

         case(9)
            dipole_component_map=6
       ! quadrupoles
!~       case(5)
!~          dipole_component_map=2
!~
!~       case(6)
!~          dipole_component_map=3
!~
!~       case(7)
!~          dipole_component_map=1
!~
!~       case(8)
!~          dipole_component_map=2
!~
!~       case(9)
!~          dipole_component_map=2
         end select
      end if

      end function dipole_component_map

      subroutine determine_target_phase(iphz,target_vectors_np1_run,&
     &    ci_vec_target_run,ntgsym,numtgt,&
     &    itarget_overall_phase, ifail)
!     ------------------------------------------------------------------
!
!     INPUT.
!
!     iphz: phase map vector from congen
!     target_vectors_np1_run: Target CI vectors created during the N+1
!       calculation
!     target_vectors_target_run: Target CI vectors created in the target
!       run
!     numtgt: number of target states per symmetry (N+1 congen order)
!     nctgt: number of target configurations per symmetry (N+1 congen
!       order)
!     ntgsym: number of symmetries involved in the target calculation
!
!     OUTPUT
!.
!     itarget_overall_phase: the phase between target states created in
!       the N and N+1 calculations
!
!     Note: I assume that the read in target CI vectors are sorted into
!     the N+1 congen order (if necessary) and held in a one dimensional
!     array. If not the order mapping needs to be determined.
!
!     Alex Harvey. Oct 2011
!     ------------------------------------------------------------------
      use cdenprop_defs
      implicit none
      integer, parameter :: thresh=1.0E-5_idp

!     Argument variables
      real(kind=idp), dimension(:) :: target_vectors_np1_run
      integer, dimension(:) :: iphz,numtgt,itarget_overall_phase
      integer :: ntgsym,ifail, iphase_tmp
      type (CIvect), dimension(:) :: ci_vec_target_run


!     Local variables
      integer :: no_target_states, itarg, icsf,iskip,&
     &           iphase_np1_run,iphase_targ_run, i,j,k, icsf2, istart, iend
      real(kind=idp), allocatable :: packed_ci_vec_target_run(:)
      ifail=0

!     Phase determination by calculating overlap between target CI vectors

      allocate( packed_ci_vec_target_run( size(target_vectors_np1_run) ) )
      packed_ci_vec_target_run=0._idp

      icsf=0
      itarg=0
      icsf2=0
      do i=1, ntgsym
!        Loop over target state in each symmetry
         do j=1,numtgt(i)
            itarg=itarg+1
!           Loop over configurations for each target state
            do k=1,ci_vec_target_run(i)%nocsf
               icsf=icsf+1
               icsf2=icsf2+1
               packed_ci_vec_target_run(icsf2)=ci_vec_target_run(i)%CV(k,j) * iphz(icsf)
            end do
            if (j .ne. numtgt(i)) icsf=icsf-ci_vec_target_run(i)%nocsf
         end do
      end do

 !     Loop over symmetries
      do i=1, ntgsym
!        Loop over target state in each symmetry
         do j=1,numtgt(i)
            itarg=itarg+1

          end do
      end do

      itarget_overall_phase=0
      itarg = 0; istart=1; iend=1


      do i=1, ntgsym
         do j=1,numtgt(i)
            itarg=itarg+1
            iend=istart+ci_vec_target_run(i)%nocsf -1

!~             write(55555,'(3i5,d20.5)') i,j,itarg, dot_product(packed_ci_vec_target_run(istart:iend),target_vectors_np1_run(istart:iend))
            itarget_overall_phase(itarg) = int(sign(1.0_idp, &
                dot_product(packed_ci_vec_target_run(istart:iend),target_vectors_np1_run(istart:iend)) ))

            istart=iend+1

         end do

      end do

      end subroutine determine_target_phase

      subroutine phases_in_target_energy_order_map(ci_vec_target_run,ntgsym,numtgt,&
     &    id_map_targ_to_congen,idtarg,itarget_irrep,target_spin,ifail)
!     ------------------------------------------------------------------
!
!     INPUT.
!

!     numtgt: number of target states per symmetry (N+1 congen order)
!     ntgsym: number of symmetries involved in the target calculation
!     ci_vec_target_run: ci vectors and energies from the target run
!
!     OUTPUT
!     id_map_targ_to_congen:  id_map_targ_to_congen(I)=K where
!     I is theindex of the target state in energy order and
!     K is the index of the state in the N+1 congen run
!
!     Note: I assume that the read in target CI vectors are sorted into
!     the N+1 congen order (if necessary) and held in a one dimensional
!     array. If not the order mapping needs to be determined.
!
!     Alex Harvey. Oct 2011
!     ------------------------------------------------------------------
      use cdenprop_defs
      use cdenprop_aux
      implicit none
      integer, parameter :: thresh=1.0E-10_idp

!     Argument variables
      integer, dimension(:) :: numtgt
      integer,allocatable, dimension(:) :: id_map_targ_to_congen,idtarg
!       integer, dimension(:,:) :: itarget_symmetry_order
      integer :: ntgsym,ifail
      type (CIvect), dimension(:) :: ci_vec_target_run

!     Local variables
      integer :: no_target_states, isym,istate,itarg
      integer,allocatable, dimension(:) ::itarget_irrep
      real(kind=idp), allocatable,dimension(:) :: target_energies,target_spin, real_map_targ_to_congen

!     Construct single array holding target energie in N+1 congen run order.
!     Might be easier to just write/read them from scatci like the vectors

      no_target_states=0
      do isym=1,ntgsym
         do istate=1,numtgt(isym)
            no_target_states= no_target_states+1
         end do
      end do
      allocate( target_energies(no_target_states),itarget_irrep(no_target_states),target_spin(no_target_states))

      itarg=0
      do isym=1,ntgsym
         do istate=1,numtgt(isym)
            itarg=itarg+1
            target_energies(itarg)=ci_vec_target_run(isym)%ei(istate)
            itarget_irrep(itarg)=ci_vec_target_run(isym)%mgvn
            target_spin(itarg)=ci_vec_target_run(isym)%S
         end do

      end do

      allocate(id_map_targ_to_congen(no_target_states),idtarg(no_target_states))
      id_map_targ_to_congen=0
      idtarg=0

      call INDEXX(no_target_states,target_energies,id_map_targ_to_congen)
      allocate(real_map_targ_to_congen(no_target_states))
      real_map_targ_to_congen=real(id_map_targ_to_congen)
      call INDEXX(no_target_states,real_map_targ_to_congen,idtarg)
!       write(332, *) idtarg

      end subroutine phases_in_target_energy_order_map

      subroutine congen_to_target_energy_order_map(ci_vec_target_run,ntgsym,numtgt,all_ntgsym, all_numtgt,&
     &    idtarg,ifail)
!     ------------------------------------------------------------------
!
!     INPUT.
!

!     numtgt: number of target states per symmetry (N+1 congen order)
!     ntgsym: number of symmetries involved in the target calculation
!     ci_vec_target_run: ci vectors and energies from the target run
!
!     OUTPUT
!     id_map_targ_to_congen:  id_map_targ_to_congen(I)=K where
!     I is theindex of the target state in energy order and
!     K is the index of the state in the N+1 congen run
!
!     Note: I assume that the read in target CI vectors are sorted into
!     the N+1 congen order (if necessary) and held in a one dimensional
!     array. If not the order mapping needs to be determined.
!
!     Alex Harvey. Oct 2011
!     ------------------------------------------------------------------
      use cdenprop_defs
      use cdenprop_aux
      implicit none
      integer, parameter :: thresh=1.0E-10_idp

!     Argument variables
      integer, dimension(:) :: numtgt, all_numtgt
      integer,allocatable, dimension(:) :: id_map_targ_to_congen,idtarg
!       integer, dimension(:,:) :: itarget_symmetry_order
      integer :: ntgsym,ifail,all_ntgsym
      type (CIvect), dimension(:) :: ci_vec_target_run

!     Local variables
      integer :: no_target_states, isym,istate,itarg,all_no_target_states, irestricted_targ, icongen_order
      integer,allocatable, dimension(:) ::itarget_irrep, ire_index
      real(kind=idp), allocatable,dimension(:) :: target_energies,target_spin, real_map_targ_to_congen,all_target_energies
      integer,allocatable, dimension(:) :: all_id_map_targ_to_congen,all_idtarg

!     Construct single array holding target energies in N+1 congen run order.

      all_no_target_states=sum(all_numtgt)
      no_target_states=sum(numtgt(1:ntgsym)) !ZM: this sum must be only over the elements that have been defined

      allocate(all_target_energies(all_no_target_states),ire_index(all_no_target_states))
      ire_index=0
      itarg=0
      irestricted_targ=0

      do isym=1,all_ntgsym
         do istate=1,all_numtgt(isym)
            itarg=itarg+1
            all_target_energies(itarg)=ci_vec_target_run(isym)%ei(istate)
            if ((isym .le. ntgsym)) then
               if ((istate .le. numtgt(isym))) then
                  irestricted_targ=irestricted_targ+1
                  ire_index(itarg)=irestricted_targ

               end if

            end if

         end do

      end do

      allocate(all_id_map_targ_to_congen(all_no_target_states),all_idtarg(all_no_target_states))
      all_id_map_targ_to_congen=0
      all_idtarg=0


      call INDEXX(all_no_target_states,all_target_energies,all_id_map_targ_to_congen)

      allocate(idtarg(no_target_states))
      do itarg=1, all_no_target_states
         irestricted_targ=ire_index(all_id_map_targ_to_congen(itarg))
         if (irestricted_targ .ne. 0) then
            idtarg(irestricted_targ)=itarg

         end if

      end do

      end subroutine congen_to_target_energy_order_map

      !Works out the correct mapping between the order of target states in the CI vectors unit "lutargci" and the target states in
      !the properties unit "lutdip". cdenprop currently assumes that target states in "lutargci" are in energy order for a given
      !symmetry. For all normal calculations this generally true. However, for PCCHF model A, the eigenstates are computed
      !separately and appeneded one-at-a-time to unit "lutargci" in the order they are generated which is not necessarily in energy
      !order.
      subroutine prop_order_ci_vec(prop_order, ncsfs, no_states, all_ci_vec_target)
          use cdenprop_defs, only: CIvect
          use algorithms, only: indexx

          integer, allocatable, intent(out) :: prop_order(:)
          type (CIvect), intent(in)         :: all_ci_vec_target(:)
          integer, intent(in)               :: ncsfs, no_states

          integer :: i, nstat
          integer :: idx
          real(kind=idp), allocatable ::  eigenvalues(:)

          allocate(eigenvalues(no_states))
          allocate(prop_order(no_states))

          idx = 1
          prop_order = 0
          eigenvalues = 0._idp

          do i = 1, ncsfs
              nstat = all_ci_vec_target (i) % nstat
              eigenvalues(idx:idx+nstat-1) = all_ci_vec_target (i) % ei
              idx = idx + nstat
          end do

          call indexx(no_states, eigenvalues, prop_order)

      end subroutine prop_order_ci_vec

end module cdenprop_procs

