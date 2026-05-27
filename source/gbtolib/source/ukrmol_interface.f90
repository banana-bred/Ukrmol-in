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
module ukrmol_interface_gbl
   use const_gbl,   only: line_len, data_file_obj_id, sym_op_nam_len, one_electron_ints, overlap_ints, kinetic_ints, stdout, &
                          property_ints, nuc_rep_att_ints, one_elham, two_el_ints, one_p_sym_ints, two_p_sym_ints, &
                          ijkl_indices_header,abel_prod_tab
   use mpi_gbl,     only: master, mpi_mod_bcast, mpi_mod_start, mpi_mod_finalize
   use precisn_gbl, only: wp
   use utils_gbl,   only: xermsg

   use atomic_basis_gbl
   use molecular_basis_gbl
   use parallel_arrays_gbl
   use integral_storage_gbl
   use mpi_memory_gbl

   private

   !> These two variables are used in SCATCI_ROUTINES to identify the format of the integral file on input.
   public data_file_obj_id, line_len

   !> The following routines are used by UKRmol-in an UKRmol-out programs (SCATCI, DENPROP, CDENPROP) and routines (SWINTERF).
   public READ_UKRMOLP_BASIS, GET_GEOM, GET_NAME_SYM, READ_UKRMOLP_INTS, READ_UKRMOLP_PROPERTY_INTS, TMG_UKPLUS, EVAL_AMPLITUDES
   public UKP_PREAMP, UKP_READAMP, GET_INTEGRAL, GET_KINETIC_ENERGY_INTEGRAL, GET_PROPERTY_INTEGRAL, GET_ECP_CORE_CHARGES
   public GET_ORB_ENERGY, CONSTRUCT_PINTEGRALS, WRITE_MOLDEN_DYSON_ORBITALS,DESTROY_UKRMOLP, START_MPI, FINALIZE_MPI

   !> ADC-RELATED ROUTINES AND VARIABLES -- START
   !> Routines to read basis information and two-electron integrals as needed in ADC codes
   public READ_ADC_INTS, GET_ENERSYMOCC
   !
   !> Returns two-electron integrals -- simplifies GET_INTEGRAL, compatible with
   !> integer(2) indices in ADC codes
   public LPQRS
   !
   !> Fock matrix
   real(kind=cfp),allocatable,public :: fock(:,:)
   !> ADC-RELATED ROUTINES AND VARIABLES -- END


   !> These arrays are visible but not changeable from outside this module.
   type(p2d_array_obj), target :: one_electron_integrals, one_positron_integrals, two_electron_integrals

   !> Integral options as read-in by READ_UKRMOLP_PROPERTY_INTS or READ_UKRMOLP_INTS.
   type(integral_options_obj) :: options

   !> Numbers mapping the type of integral with its column number in one_electron_integrals%a(:,:) and one_positron_integrals%a(:,:)
   integer :: one_elham_column = -1, kei_column = -1, nai_column = -1, property_column = -1, one_posham_column = -1
   
   !> Number of different (l,m) properties stored in one_electron_integrals.
   integer :: number_of_properties = -1

   !> Basis set atomic orbitals. Read-in by READ_UKRMOLP_BASIS.
   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis

   !> Basis set of molecular orbitals. Read-in by READ_UKRMOLP_BASIS.
   type(molecular_orbital_basis_obj), public :: molecular_orbital_basis

   !> Number of molecular orbitals in each symmetry. Determined by READ_UKRMOLP_BASIS.
   integer, allocatable :: num_orbs_sym(:)

   !> Array of dimensions (1:3,1:lmq) where lmq is the total number of (l,m) properties for which the integrals have been calculated. It is determined by READ_UKRMOLP_PROPERTY_INTS.
   integer, allocatable :: msym(:,:)

   !> Channel information (l,m,IRR) and the number of orbitals in each IRR. Determined by EVAL_AMPLITUDES.
   integer, allocatable :: continuum_channels_m_l_irr(:,:), ncont(:)

   !> Amplitudes at r=a for each orbital. Determined by EVAL_AMPLITUDES.
   real(kind=cfp), allocatable :: amplitudes(:,:)

contains

   !> This is used to start MPI in SWINTERF and HAMDIAG. In the serial integral code this is just a dummy routine that does nothing.
   subroutine START_MPI
      use mpi_gbl
      implicit none

         if (.not. mpi_started) call mpi_mod_start

   end subroutine START_MPI

   !> This is used to start MPI in SWINTERF and HAMDIAG. In the serial integral code this is just a dummy routine that does nothing.
   subroutine FINALIZE_MPI
      use mpi_gbl
      implicit none

         call mpi_mod_finalize

   end subroutine FINALIZE_MPI

   !> Converts the integer number XXX corresponding to the fort.XXX unit into actual string 'fort.XXX'.
   subroutine unit_number_2_file_name(lu,file_name)
      implicit none
      integer, intent(in) :: lu
      character(len=line_len), intent(out) :: file_name

      character(len=line_len) :: cnum

         if (lu.gt.0) then
           file_name = 'fort.'
           write(cnum,'(i0)') lu 
           file_name = trim(file_name)//trim(adjustl(cnum))
         else if (lu.eq.-1) then
           file_name = "moints_canonical"
         else
           file_name = "moints"
         end if

   end subroutine unit_number_2_file_name

   subroutine DESTROY_UKRMOLP()
     implicit none
     integer :: dum

        dum = molecular_orbital_basis%final()
        dum = two_electron_integrals%final()
        dum = one_electron_integrals%final()
        dum = one_positron_integrals%final()

   end subroutine DESTROY_UKRMOLP

   !> Reads-in the AO and MO basis sets from the fort file with a given number.
   subroutine READ_UKRMOLP_BASIS(NFTINT)
      implicit none
      integer, intent(in) :: NFTINT

      type(integral_storage_obj) :: mo_ijkl_storage
      type(data_header_obj) :: header
      character(len=line_len) :: file_name
      integer :: irr, err, pos, lunit
      logical :: op

         !we assume that the integrals are saved on the unit number 'nfti' so we need to construct the appropriate file name.
         call unit_number_2_file_name(NFTINT,file_name)

         !make sure the integral file is closed before we attempt to open it: this can be
         !an issue when calling swinterf repeatedly in one program.
         inquire(file=file_name,opened=op)
         if (op) close(NFTINT)

         !read-in all basis sets and orbital data
         call atomic_orbital_basis%read(file_name)
         molecular_orbital_basis%ao_basis => atomic_orbital_basis
         call molecular_orbital_basis%read(file_name)

         err = mo_ijkl_storage%init(disk=file_name)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'READ_UKRMOLP_BASIS', 'Error initializing the integrals storage object.', err, 1)
         end if

         !we need to get the full header corresponding to the ijkl_indices.
         err = mo_ijkl_storage % integral_file % get_header_containing (header, &
                                        molecular_orbital_basis % get_basis_name(), ijkl_indices_header)
         if (err == 0) then
            call mpi_mod_bcast(header%first_record, master)
            lunit = mo_ijkl_storage%integral_file%get_unit_no()
            call mpi_mod_bcast(lunit, master)
            call molecular_orbital_basis%read_ijkl_indices(lunit, file_name, header%first_record, pos)
         end if

         call mo_ijkl_storage%final

         if (allocated(num_orbs_sym)) deallocate(num_orbs_sym)
         allocate(num_orbs_sym(molecular_orbital_basis%no_irr),stat=err)
          if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_BASIS','Memory allocation failed.',err,1)

         do irr=1,molecular_orbital_basis%no_irr
            num_orbs_sym(irr) = molecular_orbital_basis%get_number_of_orbitals(irr)
         enddo !i

   end subroutine READ_UKRMOLP_BASIS

   !> Returns the nuclear data for the problem as contained in the molecular_orbital_basis object.
   !> \warning This routine ALWAYS returns double precision values due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   subroutine GET_GEOM(nnuc,cname,xnuc,ynuc,znuc,charge)
      implicit none
      integer, intent(out) :: nnuc
      real(kind=wp), allocatable, dimension(:), intent(out) :: charge, xnuc, ynuc, znuc
      character(len=8), allocatable, dimension(:), intent(out) :: cname

      integer :: i, j, k, err

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'GET_GEOM', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         nnuc = molecular_orbital_basis%symmetry_data%no_nuc
         allocate(charge(nnuc), xnuc(nnuc), ynuc(nnuc), znuc(nnuc), cname(nnuc),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','GET_GEOM','Memory allocation failed.',err,1)

         associate(nucleus => molecular_orbital_basis%symmetry_data%nucleus, nuclear_data => molecular_orbital_basis%symmetry_data)
            k = 0
            do i=1,nuclear_data%no_sym_nuc       !over all symmetrically non-redundant nuclei
               k = k + 1
               !data for the non-redundant nucleus itself
               write(cname(k),'(a2,i2)') nuclear_data%non_red_nuc(i)%name, i
               xnuc(k) = nuclear_data%non_red_nuc(i)%center(1)
               ynuc(k) = nuclear_data%non_red_nuc(i)%center(2)
               znuc(k) = nuclear_data%non_red_nuc(i)%center(3)
               charge(k) = nuclear_data%non_red_nuc(i)%charge
               do j=1,nuclear_data%no_eqv_nuc(i) !over all symmetrical partners
                  k = k + 1
                  !data for the symmetrical partner
                  write(cname(k),'(a2,i2)') nuclear_data%eqv_nuc(i,j)%name, i
                  xnuc(k) = nuclear_data%eqv_nuc(i,j)%center(1)
                  ynuc(k) = nuclear_data%eqv_nuc(i,j)%center(2)
                  znuc(k) = nuclear_data%eqv_nuc(i,j)%center(3)
                  charge(k) = nuclear_data%eqv_nuc(i,j)%charge
               enddo !j
            enddo !i
         end associate

   end subroutine GET_GEOM

   !> Returns the number of symmetries and some other basic parameters as obtained from the initialized object molecular_orbital_basis and property_integrals.
   subroutine GET_NAME_SYM(name,no_sym,nob,nlmq)
      implicit none
      character(len=132), intent(out) :: name
      integer, intent(out) :: no_sym, nob(:), nlmq

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'GET_NO_SYM', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         no_sym = molecular_orbital_basis%no_irr
         name = property_ints

         if (size(nob) < no_sym) call xermsg('ukrmol_interface','GET_NO_SYM','The array NOB on input is too small.',2,1)
         nob(1:no_sym) = num_orbs_sym(1:no_sym)

         if (number_of_properties > 0) then
            nlmq = number_of_properties !how many (l,m) property integrals we have (with l=0,1,...).
         else
            call xermsg('ukrmol_interface','GET_NO_SYM','Cannot supply nlmq: property have not been read-in.',3,1)
         endif

   end subroutine GET_NAME_SYM

   !> This routine reads-in all one electron and two electron integrals needed for SCATCI calculations. It is an analogoue of the INTSIN routine from SCATCI.
   !> \warning This routine ALWAYS requires double precision value of scalem on input due to the need for compatibility with UKRmol-in/out which is always compiled using double precision. Note that 
   !> this results in reducing the values in ke_integrals%a into double precision.
   !> NEW FEATURE: WORKS WITH SHARED MEMORY OF THE MPI-3.0 STANDARD. This means that all integrals are read-in the shared mode so in the MPI regime there is only one copy of the integral arrays PER NODE.
   !> If we do not have the shared memory capability then the code resorts to the local mode, i.e. each MPI task keeps its own copy of all integral arrays.
   subroutine READ_UKRMOLP_INTS(nfte,nfti,lembf,nint1e,nint2e,nocsf,nfta,isymtp,nsym1,nob1,iposit,scalem,name,nalm,qmoln)
      use symmetry_gbl
      implicit none
      integer, intent(in) :: nfte,nfti,lembf,nocsf,nfta,isymtp,nsym1,nob1(nsym1),iposit
      integer, intent(inout) :: nint1e, nint2e
      real(kind=wp), intent(in) :: scalem
      character(len=120), intent(in) :: name
      logical, intent(in) :: qmoln
      integer, intent(inout) :: nalm

      !local variables
      character(len=line_len) :: file_name
      type(integral_storage_obj) :: mo_integral_storage
      type(integral_storage_obj), target :: tgt_storage
      type(data_header_obj) :: header
      integer :: i, j, err, nnuc, ecp_i, ecp_j
      real(kind=cfp) :: inv_scale
      real(kind=wp) :: nuc_repulsion, distance
      real(kind=wp), dimension(41) :: dtnuc
      integer, dimension(20) :: nhe
      integer, parameter :: number_of_blocks = 0
      logical, parameter :: tgt_is_local = .true.
      character(len=line_len), allocatable :: column_descriptor(:)

         !error reporting flag used in UKRmol but not here: we report errors immediately
         nalm = 0
         ! BC - down grade this to a warning rather than a hard exit
         if (qmoln) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Use of the QMOLN flag has not been implemented.',1,0)

         !we assume that the integrals are saved on the unit number 'nfti'
         call READ_UKRMOLP_BASIS(NFTI)

         !check for data input errors
         if (molecular_orbital_basis%no_irr < nsym1) then !nsym1 is the number of the last IRR where we have some orbitals
            print *,molecular_orbital_basis%no_irr, nsym1
            call xermsg ('ukrmol_interface', 'READ_UKRMOLP_INTS', &
                         'The number of IRRs obtained from the orbital basis data is not compatible &
                         &with the one given by SCATCI.', 2, 1)
         endif

         do i=1,nsym1
            if (num_orbs_sym(i) .ne. nob1(i)) then
               write(nfta,451)
               write(nfta,452) (num_orbs_sym(j),j=1,nsym1)
               write(nfta,452) (nob1(j),j=1,nsym1)
 451           format(/' MISMATCH IN NOB BETWEEN INTEGRALS AND FORMULAE')
 452           format(' NOB=',20i5)
               call xermsg ('ukrmol_interface', 'READ_UKRMOLP_INTS', &
                            'The number of orbitals contained in the orbital basis data is not equal &
                            &to the one expected by SCATCI.', 3, 1)
            endif
         enddo

         !open the data file and analyse the headers manually
         call unit_number_2_file_name(NFTI,file_name)
         err = mo_integral_storage%init(disk=file_name)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'READ_UKRMOLP_INTS', 'Error initializing the integrals storage object.', err, 1)
         end if

         write(nfta,1600) nint1e, nint2e
 1600    format(/,' On input SCATCI expects:',/,' 1-electron integrals, NINT1e =',i0,/,' 2-electron integrals, NINT2e =',i0)

         !read-in the 1-electron integrals into one_electron_integrals
         err = tgt_storage%init(memory=one_electron_integrals)
         tgt_storage%integral_file%identifier = mo_integral_storage%integral_file%identifier
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error initializing tgt_storage.',err,1)

         !we need to get the full header corresponding to the transformed overlap and kinetic energy integrals.
         err = mo_integral_storage % integral_file % get_header_containing(header, &
                                                                           molecular_orbital_basis % get_basis_name(), &
                                                                           one_p_sym_ints)
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error locating the transformed 1-electron integrals; &
                                     &see data_header_obj%get_header_containing for details.',err,1)

         call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
         call tgt_storage%final

         write(nfta,'(/,"One electron integrals read-in.")')

         !Column number in one_electron_integrals%a corresponding to the 1-electron Hamiltonian integrals: these integrals must ALWAYS be calculated.
         one_elham_column = one_electron_integrals%find_column_matching_name(one_elham)

         if (scalem .ne. 1.0_cfp) then
            write(nfta,'("Kinetic energy integrals will be rescaled (divided) by: ",e25.15)') scalem

            !Column number in one_electron_integrals%a corresponding to the kinetic energy integrals.
            kei_column = one_electron_integrals%find_column_matching_name(kinetic_ints)
   
            !Column number in one_electron_integrals%a corresponding to the nuclear attraction integrals.
            nai_column = one_electron_integrals%find_column_matching_name(nuc_rep_att_ints)

            !Rescale the kinetic energy integrals and recompute the 1-electron Hamiltonian integrals for the electron case.
            inv_scale = 1.0_cfp/scalem

            call one_electron_integrals%multiply_column(kei_column,inv_scale)
            !Direct summation below implies that the indexing of the kinetic and nuclear integrals is the same
            call one_electron_integrals%combine_columns(kei_column,'+',nai_column,one_elham_column)
            
         endif

         !Construct the one-positron integrals: K - V
         if (iposit .ne. 0) then
            write(nfta,'("Recomputing the 1-particle Hamiltonian integrals for the positron case...")')

            !Column number in one_electron_integrals%a corresponding to the kinetic energy integrals.
            kei_column = one_electron_integrals%find_column_matching_name(kinetic_ints)
   
            !Column number in one_electron_integrals%a corresponding to the nuclear attraction integrals.
            nai_column = one_electron_integrals%find_column_matching_name(nuc_rep_att_ints)

            call one_electron_integrals%get_column_descriptor(column_descriptor)

            !Allocate space for the 1p Hamiltonian integrals
            err = one_positron_integrals % init(size(one_electron_integrals % a, 1), &
                                                size(one_electron_integrals % a, 2), &
                                                number_of_blocks, &
                                                column_descriptor) 
            if (err /= 0) then
                call xermsg ('molecular_orbital_basis_obj', 'READ_UKRMOLP_INTS', &
                             'Initialization of one_positron_integrals failed; see p2d_array%init.', err, 1)
            end if

            one_posham_column = one_electron_integrals%find_column_matching_name(one_elham)

            one_positron_integrals%a = one_electron_integrals%a

            !Direct subtraction below implies that the indexing of the kinetic and nuclear integrals is the same
            call one_positron_integrals%combine_columns(kei_column,'-',nai_column,one_posham_column)
            
            write(nfta,'("...done")')
            
         endif

         !read-in the 2p integrals
         err = tgt_storage%init(memory=two_electron_integrals)
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error initializing tgt_storage.',err,1)

         !we need to get the full header corresponding to the transformed 2-electron integrals.
         err = mo_integral_storage % integral_file % get_header_containing(header, &
                                                                           molecular_orbital_basis % get_basis_name(), &
                                                                           two_p_sym_ints)
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Error locating the transformed 2-electron integrals; &
                                     &see data_header_obj%get_header_containing for details.',err,1)

         call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
         call tgt_storage%final

         write(nfta,'(/,"Two electron integrals read-in.")')

         !Now update the integral counters by the actual number of integrals read-in
         nint1e = size(one_electron_integrals%a,1)
         nint2e = size(two_electron_integrals%a,1)
         write(nfta,1700) nint1e, nint2e
 1700    format(/,' UKRMol+ integrals read successfully:',&
                /,' 1-electron integrals, NINT1e =',i0,&
                /,' 2-electron integrals, NINT2e =',i0)

         !Compute the nuclear repulsion energy
         associate(nucleus => molecular_orbital_basis%symmetry_data%nucleus, nuclear_data => molecular_orbital_basis%symmetry_data)
            nuc_repulsion = 0.0_cfp
            nnuc = nuclear_data%no_nuc
            do j=2,nnuc
               do i=1, j-1
                  if (nucleus(j)%charge*nucleus(i)%charge .eq. 0.0_cfp) cycle
                  ecp_i = nucleus(i)%ecp%n_core !ECP charges on nucleus i
                  ecp_j = nucleus(j)%ecp%n_core !ECP charges on nucleus j
                  distance = sqrt(dot_product(nucleus(j)%center-nucleus(i)%center,nucleus(j)%center-nucleus(i)%center))
                  if (distance .eq. 0.0_cfp) call xermsg('ukrmol_interface','READ_UKRMOLP_INTS','Nuclei too close.',10,1)
                  nuc_repulsion = nuc_repulsion + (nucleus(j)%charge-ecp_j)*(nucleus(i)%charge-ecp_i)/distance
               end do
            end do
            !The assignments to dtnuc effectively replicate the equivalence statements in INTSIN in scatci_build. 
            !However we drop the assignment to dtnuc(22) which is set only when Alchemy integrals are used and seems completely redundant.
            dtnuc(:) = 0.0_cfp
            dtnuc(1) = nuc_repulsion
            dtnuc(2) = nucleus(1)%charge !is this redundant?
         end associate

         write(nfta,1510) nuc_repulsion
 1510    format(/,10x,'Nuclear repulsion energy = ',f15.7)
 
         !The assignment to nhe effectively replicates the equivalence statement in INTSIN in scatci_build
         nhe(:) = 0
         nhe(1) = num_orbs_sym(1) !is this redundant?

         !Write header on energy matrix file  
         open(unit=nfte,form='unformatted')
         write(nfte) nocsf, lembf, 0, nocsf, 0, nsym1, 0, 0, 0, 0, nnuc, 0, name, nhe, dtnuc

         call mo_integral_storage%final

         write(nfta,'(/,10x,"Finished reading UKRMol+ integrals.")')

   end subroutine READ_UKRMOLP_INTS

   !> Reads-in the basis sets and the property integrals.
   subroutine READ_UKRMOLP_PROPERTY_INTS(NFTINT,IWRITE)
      use symmetry_gbl
      implicit none
      integer, intent(in) :: NFTINT, IWRITE

      !local variables
      character(len=line_len) :: file_name
      type(integral_storage_obj) :: mo_integral_storage
      type(integral_storage_obj), target :: tgt_storage
      type(data_header_obj) :: header
      integer :: i, err, l, m
      character(len=sym_op_nam_len) :: nam
      integer, parameter :: no_blocks = 0
      logical, parameter :: tgt_is_local = .true.
      character(len=line_len), allocatable :: column_descriptor(:)

         !we assume that the integrals are saved on the unit number 'nfti'
         call READ_UKRMOLP_BASIS(NFTINT)

         !open the data file and analyse the headers manually
         call unit_number_2_file_name(NFTINT,file_name)
         err = mo_integral_storage%init(disk=file_name)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'READ_UKRMOLP_PROPERTY_INTS', &
                         'Error initializing the integrals storage object.', err, 1)
         end if

         !we need to get the full header corresponding to the transformed property integrals.
         err = mo_integral_storage % integral_file % get_header_containing(header, &
                                                                           molecular_orbital_basis % get_basis_name(), &
                                                                           one_p_sym_ints)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'READ_UKRMOLP_PROPERTY_INTS', &
                         'Error locating the transformed 1-electron integrals; &
                         &see data_header_obj%get_header_containing for details.', err, 1)
         end if

         !read-in all 1-electron integrals
         err = tgt_storage%init(memory=one_electron_integrals)
         tgt_storage%integral_file%identifier = mo_integral_storage%integral_file%identifier
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_PROPERTY_INTS','Error initializing tgt_storage.',err,1)

         call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
         call tgt_storage%final

         write(IWRITE,'(/,"1-electron integrals read-in.")')

         property_column = one_electron_integrals%find_column_matching_name(property_ints)

         call one_electron_integrals%get_column_descriptor(column_descriptor)
         
         number_of_properties = 0
         do i=1,size(column_descriptor)
            if (column_descriptor(i) .eq. property_ints) number_of_properties = number_of_properties + 1
         enddo
         
         write(IWRITE,'(/,"Number of types of property integrals: ",i4)') number_of_properties

         if (allocated(msym)) deallocate(msym)
         allocate(msym(1:3,number_of_properties),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_UKRMOLP_PROPERTY_INTS','Memory allocation 2 error',err,1)

         !build the (l,m) property numbers up to the last property.
         i = 1
         l = 0
         m = 0 
         do
            msym(1,i) = m
            msym(2,i) = l

            !Find out the symmetry of the property operator (l,m) which is the same as the IRR of the spherical harmonic in the P-G of the molecule.
            msym(3,i) = molecular_orbital_basis%symmetry_data%sph_harm_pg_sym(l,m,nam)

            write(level2,'(/,10x,"Property (l,m): (",i4,",",i4,")")') l,m
            write(level2,'(10x,"Property symmetry: ",i2)') msym(3,i)

            if (i .eq. number_of_properties) exit

            if (m .eq. l) then
               l = l + 1
               m = -l
            else
               m = m + 1
            endif

            i = i + 1
         enddo

         !Make sure we always close the integrals file to ensure we don't get errors in CDENPROP when the properties are being read-in multiple times.
         call mo_integral_storage%final

   end subroutine READ_UKRMOLP_PROPERTY_INTS

   !> Returns a one electron/positron integral (K+-V) or two-electron integral given the basis function indices and the positron flag. If c,d == 0 then the one electron/positron integral is returned.
   !> If positron_flag > 0 then the positronic one-particle integral (K-V) is returned.
   !> \warning This routine ALWAYS returns double precision value due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   !> NEW FEATURE: WORKS WITH SHARED MEMORY OF THE MPI-3.0 STANDARD
   function GET_INTEGRAL(a,b,c,d,positron_flag)
      implicit none
      integer, intent(in) :: a,b,c,d, positron_flag
      real(kind=wp) :: GET_INTEGRAL

      integer :: ind(1), two_ind(1:2,1), four_ind(1:4,1)

         GET_INTEGRAL = 0.0_wp
         if (c .le. 0 .or. d .le. 0) then !one particle case
            two_ind(1,1) = a
            two_ind(2,1) = b
            ind(1:1) = molecular_orbital_basis%integral_index(one_elham,two_ind,options%two_p_continuum)
            if (positron_flag .ne. 0) then !positron case: K-V, where V is calculated for electrons
               GET_INTEGRAL = one_positron_integrals%a(ind(1),one_posham_column)
            else !electron case: K+V, where V is calculated for electrons
               GET_INTEGRAL = one_electron_integrals%a(ind(1),one_elham_column)
            endif
            !write(*,'("1el",2i3,i,e25.15)') a,b,ind(1),GET_INTEGRAL
            !write(100,'("1el",2i3,e25.15)') a,b,GET_INTEGRAL
         else !two electron case
            four_ind(1,1) = a
            four_ind(2,1) = b
            four_ind(3,1) = c
            four_ind(4,1) = d
            ind(1:1) = molecular_orbital_basis%integral_index(two_el_ints,four_ind,options%two_p_continuum)
            if (ind(1) > 0) GET_INTEGRAL = two_electron_integrals%a(ind(1),1)
            !write(*,'("2el",4i3,i,e25.15)') a,b,c,d,ind(1),GET_INTEGRAL
            !write(101,'("2el",4i3,e25.15)') a,b,c,d,GET_INTEGRAL
         endif  

   end function GET_INTEGRAL

   function GET_ORB_ENERGY(a)
      implicit none
      integer, intent(in) :: a
      real(kind=wp) :: GET_ORB_ENERGY
      integer :: i, j, num
      type(orbital_data_obj) :: orbital_data

         GET_ORB_ENERGY = 0.0

         do j = 1,molecular_orbital_basis%no_irr
           call molecular_orbital_basis%get_shell_data(j,orbital_data)
           do i = 1, orbital_data%number_of_functions
              num = molecular_orbital_basis%get_absolute_index(i,j)
              if (num .eq. a) then
                GET_ORB_ENERGY = orbital_data%energy(i)
                return
              endif
           enddo
         enddo
         write(*,"(A,I3)") "Orbital index ", a
         call xermsg ('ukrmol_interface', 'GET_ORB_ENERGY', 'Orbital energy not found.', 1, 1)

   end function GET_ORB_ENERGY

   !> \brief   Function to return the kinetic energy integral as these integrals
   !>          are used in the calculation of BEB cross sections.
   !> \authors Bridgette Cooper
   !> \date    October 2019
   !> \warning This routine ALWAYS returns double precision value due to the need for compatibility
   !>          with UKRmol-in/out which is always compiled using double precision.
   
   function GET_KINETIC_ENERGY_INTEGRAL(a,b)
      implicit none
      integer, intent(in) :: a,b 
      real(kind=wp) :: GET_KINETIC_ENERGY_INTEGRAL
      integer :: ind(1), two_ind(1:2,1)
    
         kei_column = one_electron_integrals%find_column_matching_name(kinetic_ints)
         GET_KINETIC_ENERGY_INTEGRAL = 0.0_wp
         two_ind(1,1) = a
         two_ind(2,1) = b
         ind(1:1) = molecular_orbital_basis%integral_index(kinetic_ints, two_ind, options%two_p_continuum)
         GET_KINETIC_ENERGY_INTEGRAL = one_electron_integrals%a(ind(1),kei_column)

   end function GET_KINETIC_ENERGY_INTEGRAL


   !> \brief   Retrieve value of a one-particle property integral
   !> \authors J Benda
   !> \date    2023
   !>
   !> Assuming that the subroutine READ_UKRMOLP_PROPERTY_INTS has been called before, retrieve the one-electron property
   !> integral. This function is currently used in CDENPROP.
   !>
   !> \param[in] iprop   Index of the property (1 = s0, 2 = py, 3 = pz, 4 = px, etc)
   !> \param[in] a       First molecular orbital index (continuous index in concatenated symmetries)
   !> \param[in] b       Second molecular orbital index (continuous index in concatenated symmetries)
   !> \param[in] iposit  Return electronic property (if absent or zero) or positronic one (otherwise)
   !>
   !> \return Value of the integal as a double precision real number.
   !>
   real(wp) function GET_PROPERTY_INTEGRAL(iprop, a, b, iposit) result(intg)
      implicit none
      integer, intent(in) :: iprop, a, b
      integer, optional   :: iposit
      integer :: indx(1)

         indx = molecular_orbital_basis % integral_index(property_ints, reshape([a,b], [2,1]), options%two_p_continuum)

         if (iprop > number_of_properties) then
            call xermsg('ukrmol_interface', 'GET_PROPERTY_INTEGRAL', 'Requested property out of range', number_of_properties, 1)
         end if

         if (indx(1) > size(one_electron_integrals % a, 1)) then
            call xermsg('ukrmol_interface', 'GET_PROPERTY_INTEGRAL', 'Requested property integral out of range', 1, 1)
         end if

         if (present(iposit)) then
            if (iposit /= 0) then
               intg = one_positron_integrals % a(indx(1), property_column - 1 + iprop)
               return
            end if
         end if

         intg = one_electron_integrals % a(indx(1), property_column - 1 + iprop)

   end function GET_PROPERTY_INTEGRAL


   !> Multiplies the density matrix elements with the set of property integrals with (l,m) given by the linear index lmq = l*l+l+m+1. This routine may only be called following a call to 
   !> READ_UKRMOLP_PROPERTY_INTS. We assume that this routine is called with the first call having lmq=1. This is needed for the index function initialization. The input values are: the wavefunction symmetry
   !> difference (mdel), the number of elements in the density matrix (lden), the property l,m,q sequence number (lmq) and the block data constructed in DENPROP/TMGP. The output are: the l,m,q values 
   !> corresponding to lmq and the value for the electronic contribution to the property (prlmq).
   !> \warning This routine ALWAYS uses double precision values of the property integrals due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   subroutine TMG_UKPLUS(mdel,block_data,no_blocks,nob,den,lden,prlmq,lmq,l,m,q)
      implicit none
      integer, intent(in) :: mdel, lden, lmq, block_data(:,:), nob(:), no_blocks
      integer, intent(out) :: l,m,q
      real(kind=wp), dimension(lden) :: den
      real(kind=wp) :: prlmq
 
      integer :: isym, jsym, nsym, i, orb_i, orb_j, last, o_sym_i_start, o_sym_j_start, two_ind(1:2,1:1), ind(1:1), block

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'TMG_UKPLUS', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         l = msym(2,lmq)
         m = abs(msym(1,lmq))
         if (m .eq. 0) then 
            q = 0
         else
            q = sign(1,msym(1,lmq))
         endif

         nsym = molecular_orbital_basis%no_irr

         prlmq = 0.0_wp

         !Here we loop over the pairs of symmetries of orbitals that contribute to the density matrix. The order of the blocks of symmetries is important since this order corresponds to the order of the
         !blocks (and elements) of the density matrix (indexed by 'i').
         i = 0
         do block = 1,no_blocks
            if (block_data(1,block) .eq. mdel) then

               isym = block_data(2,block)
               jsym = block_data(3,block)

               !offsets for the orbital indices from the symmetries isym,jsym
               o_sym_i_start = sum(num_orbs_sym(1:isym-1))
               o_sym_j_start = sum(num_orbs_sym(1:jsym-1))
   
               !Note that the loops over orbitals are loops only over the orbitals that are present in the CONGEN input, i.e. we don't loop over all orbitals in the basis set: that number might be different
               !to the one used to construct the wavefunctions. The order of the loops orb_j, orb_i is also important:
               !we loop over the property integrals in the order in which the density matrix elements DEN(1:lden) were constructed (see DENPROP routines SETMBA, RWDIJ, INDEX1).
               last = nob(isym)
               do orb_j = 1,nob(jsym)
                  if (isym .eq. jsym) last = orb_j !only unique combinations of orbitals are used in the density matrix expressions
                  do orb_i = 1,last
                     i = i + 1
                     if (i .le. lden .and. DEN(i) .ne. 0.0_wp) then !there may be a non-zero contribution to the property
                        two_ind(1,1) = o_sym_i_start + orb_i !the index of the orbital in the set of orbitals for all symmetries, sim. for the second orbital.
                        two_ind(2,1) = o_sym_j_start + orb_j
                        ind(1:1) = molecular_orbital_basis%integral_index(property_ints,two_ind,options%two_p_continuum) !index of the property integral <orb_i|(l,m)|orb_j> for the orbitals from the symmetries isym, jsym.
                        prlmq = prlmq + DEN(i)*one_electron_integrals%a(ind(1),property_column-1+lmq)
                     endif
                  enddo !orb_j
               enddo !orb_i

            endif
         enddo

   end subroutine TMG_UKPLUS

   !> This routine interfaces with CDENPROP. It transfers the necessary orbital data into the argument variables. We use
   !> READ_UKRMOLP_PROPERTY_INTS to read the integrals from disk and they are ready for retrieval using GET_PROPERTY_INTEGRAL.
   subroutine construct_pintegrals (nftprop,iwrite,no_of_integrals,no_of_properties,nilmq,lp,mp,qp,property_name,nob,mob,mpob)
      implicit none
      integer, intent(in) :: nftprop, iwrite
      integer, intent(out) :: no_of_integrals,no_of_properties
      integer, allocatable :: nilmq(:), lp(:), mp(:), qp(:), nob(:), mob(:), mpob(:)
      character(len=8), allocatable :: property_name(:)

      integer :: nobt, lmq, i, j, two_ind(1:2,1:1), ind(1:1), last, k, m, integral_index, sym_i, sym_j, iorb, jorb, &
                 relative_index, err
      integer, allocatable :: istart(:,:)

         !Read all basis sets and integrals into memory
         call READ_UKRMOLP_PROPERTY_INTS(nftprop,iwrite)

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'construct_pintegrals', &
                         'The orbital basis set data has not been initialized.', 1, 1)
         end if

         nobt = sum(num_orbs_sym(:))
         allocate(nob(molecular_orbital_basis%no_irr),mob(nobt),mpob(nobt),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','construct_pintegrals','Memory allocation 1 failed.',err,1)

         !get dimensions of the integral array:
         no_of_properties = number_of_properties
         no_of_integrals = size(one_electron_integrals%a,1)

         allocate(lp(no_of_properties), mp(no_of_properties), nilmq(no_of_properties), &
                  qp(no_of_properties),property_name(no_of_properties), stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','construct_pintegrals','Memory allocation 2 failed.',err,1)
         nilmq(:) = no_of_integrals

         nob(1:molecular_orbital_basis%no_irr) = num_orbs_sym(1:molecular_orbital_basis%no_irr)
         do i=1,nobt
            mpob(i) = molecular_orbital_basis%get_index_within_symmetry(i)  !orbital index within its symmetry
            mob(i) = molecular_orbital_basis%get_orbital_symmetry(i)-1 != symmetry of the orbital - 1
         enddo

         do lmq=1,no_of_properties

            lp(lmq) = msym(2,lmq) 
            mp(lmq) = abs(msym(1,lmq))
            if (mp(lmq) .eq. 0) then
               qp(lmq) = 0
            else 
               qp(lmq) = sign(1,msym(1,lmq))
            endif
            if (lp(lmq) .eq. 0) property_name(lmq) = 'overlap'
            if (lp(lmq) .eq. 1) property_name(lmq) = 'dipole'
            if (lp(lmq) .eq. 2) property_name(lmq) = 'quadpole'
            if (lp(lmq) > 2) property_name(lmq) = '>2 mpole'

         enddo !lmq

   end subroutine construct_pintegrals

   !> This routine takes the Dyson orbitals produced by CDENPROP and writes them out in the Molden file format along with the GTO basis set. The CDENPROP Dyson's are expressed as
   !> a linear combination of the target and continuum MOs. Additionally to producing the Molden file this routine always saves the complete set of Dyson orbitals (along with the AO
   !> basis set) into an object of orbital_data_obj type and then writes this to disk. This ensures that we always retain the full information about the Dyson orbitals (in case
   !> the AO basis contains high L that Molden cannot handle).
   !> The value of the R-matrix radius stored in 'options' is used to normalize the continuum functions to the R-matrix radius. These normalization coefficients are multiplied in with the orbital coefficients
   !> to produce the Molden file. The orbital coefficients saved in the UKRmol+ format are modified, too. Therefore when calculating the radial densities of the Dysons from the UKRmol+ file using the
   !> routine orbital_basis_data_obj%radial_charge_density the input value of rmat_radius for this routine must be set to a value .le. 0.0_cfp. See also orbital_basis_data_obj%radial_charge_density.
   !>
   !> By default, the Molden file will use Cartesian basis. This can be overriden by setting the optional parameter
   !> `molden_format_opt` to 2, which will instead produce a spherical Molden file.
   !>
   !> \warning This routine multiplies in the spin-related CG coefficient (see commit 879 to cdenprop_io) but it is hard coded to the value 1/sqrt(2). This must be changed for your particular case.
   !> \warning This routine ALWAYS REQUIRES double precision values of the coefficients on input due to the need for compatibility with CDENPROP which is always compiled 
   !> using double precision.
   subroutine write_molden_dyson_orbitals(dyson_orbitals,idyson_orbital_irrep,itarget_irrep,itarget_spin,dyson_orbital_norms,&
                                          molden_format_opt)
      use molden_gbl
      use molden_const_gbl, only: cart_basis_fmt
      use symmetry_gbl
      use mpi_gbl, only: myrank
      use gto_routines_gbl, only: cms_gto_norm
      use utils_gbl, only: delete_file
      implicit none
      real(kind=wp), intent(in) :: dyson_orbitals(:,:,:)
      integer, intent(in) :: idyson_orbital_irrep(:),itarget_irrep(:),itarget_spin(:)
      real(kind=wp), intent(in) :: dyson_orbital_norms(:,:)
      integer, optional, intent(in) :: molden_format_opt

      integer :: no_orbitals, no_target_states, no_bound_states, molden_format
      type(geometry_obj) :: geometry
      type(molecular_orbital_basis_obj) :: Dyson_orbital_basis
      type(orbital_data_obj), allocatable :: orbital_data(:), orbital_data_no_btos(:)
      type(molden_input_obj) :: molden_file
      integer :: i, j, k, n_cf, ao, err, irr, number_of_cgto_shells
      character(len=line_len) :: molden_file_name, cnum, dyson_basis_set_file
      integer, allocatable :: sym_ind(:), n_dyson_orbitals(:)
      type(CGTO_shell_data_obj), allocatable :: CGTO_shells(:)
      real(kind=cfp), allocatable :: cf(:,:), normalization_factor(:)
      real(kind=cfp) :: rmat_radius, spin_cg, fac

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'The orbital basis set data has not been initialized.', 1, 1)
         end if

         no_orbitals = size(dyson_orbitals,1)
         no_target_states = size(dyson_orbitals,2)
         no_bound_states = size(dyson_orbitals,3)

         if (size(dyson_orbital_norms,1) /= no_target_states .or. &
             size(dyson_orbital_norms,2) /= no_bound_states) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'On input the array dyson_orbital_norms does not have the expected dimensions.', 2, 1)
         end if

         if (sum(num_orbs_sym(:)) /= no_orbitals) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'Number of orbitals in CDENPROP is different to the number of orbitals in the orbital basis set.', 3, 1)
         end if

         allocate(orbital_data(1:molecular_orbital_basis%no_irr),&
                  sym_ind(1:molecular_orbital_basis%no_irr),&
                  n_dyson_orbitals(1:molecular_orbital_basis%no_irr),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Memory allocation 1 failed',err,1)

         !Count the number of Dyson orbitals in each symmetry (this replicates the double loop from below)
         n_dyson_orbitals = 0
         do k=1,no_bound_states
            do j=1,no_target_states
               !the Dyson orbital symmetries given by CDENPROP start with 0
               n_dyson_orbitals(idyson_orbital_irrep(j)+1) = n_dyson_orbitals(idyson_orbital_irrep(j)+1) + 1
            enddo !j
         enddo !k

         !Note that this must be amended once symmetry for the AOs is used!
         n_cf = molecular_orbital_basis%ao_basis%number_of_functions
         do i=1,molecular_orbital_basis%no_irr

            allocate(orbital_data(i) % coefficients(n_cf,n_dyson_orbitals(i)), &
                     orbital_data(i) % energy(n_dyson_orbitals(i)), &
                     orbital_data(i) % spin(n_dyson_orbitals(i)), &
                     orbital_data(i) % occup(n_dyson_orbitals(i)), stat = err)
            if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Memory allocation 2 failed',err,1)

            orbital_data(i)%number_of_coefficients = n_cf
            orbital_data(i)%irr = i
            orbital_data(i)%number_of_functions = n_dyson_orbitals(i)
            orbital_data(i)%point_group = molecular_orbital_basis%pg
            orbital_data(i)%occup = 1.0_cfp
            orbital_data(i)%spin = 1
            orbital_data(i)%energy = 0.0_cfp
         enddo

         call molecular_orbital_basis%ao_basis%get_all_CGTO_shells(CGTO_shells,number_of_cgto_shells)

         !Calculate normalization factors for the continuum GTOs and multiply them in with the orbital coefficients.
         allocate(normalization_factor(molecular_orbital_basis%ao_basis%number_of_functions),stat=err)
         if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Memory allocation -1 failed.',err,1)
         normalization_factor = 1.0_cfp

         rmat_radius = options%a
         ao = 0
         do i=1,number_of_cgto_shells
            if (rmat_radius > 0.0_cfp .and. molecular_orbital_basis%ao_basis%shell_descriptor(3,i) .eq. 1) then !Normalize the continuum to the R-matrix sphere
               fac = cms_gto_norm (rmat_radius, &
                                   CGTO_shells(i) % l, &
                                   CGTO_shells(i) % number_of_primitives, &
                                   CGTO_shells(i) % exponents, &
                                   CGTO_shells(i) % contractions, &
                                   CGTO_shells(i) % norm, &
                                   CGTO_shells(i) % norms)
               write(level2,'("Continuum normalization factor",i0,e25.15)') i,fac
               normalization_factor(ao+1:ao+2*CGTO_shells(i)%l+1) = fac
            elseif (molecular_orbital_basis%ao_basis%shell_descriptor(3,i) .eq. 1 .and. rmat_radius .le. 0.0_cfp) then
               write(level2,'("Continuum functions in this shell will not be normalized to &
                                &the R-matrix sphere since rmat_radius .le. 0.0_cfp.")')
            endif
            ao = ao + 2*CGTO_shells(i)%l+1
         enddo !i

         call molecular_orbital_basis%get_orbital_coefficient_matrix(cf)
         !
         !TRANSFORM THE DYSON ORBITALS FROM THE MOLECULAR ORBITAL BASIS INTO THE ATOMIC ORBITAL BASIS:
         !
         spin_cg = 1.0_cfp !1.0_cfp/sqrt(2.0_cfp) !Spin related Clebsch-Gordan coefficient (see commit 879 to trunk:cdenprop_io).
         write(level2,'("The value of the spin related Clebsch-Gordan coefficient is: ",e25.15)') spin_cg

         write(level2,'("Dyson orbital index to target state mapping: bound state,index within sym.,symmetry,target &
            &state sym.,target state spin")')
         sym_ind(:) = 0
         do k=1,no_bound_states
            do j=1,no_target_states
               irr = idyson_orbital_irrep(j)+1 !the Dyson orbital and target state symmetries given by CDENPROP start with 0
               sym_ind(irr) = sym_ind(irr) + 1
               write(level2,'("Mapping: ",i3,1X,i5,1X,i1,1X,i1,1X,i1)') k,sym_ind(irr),irr, itarget_irrep(j)+1, itarget_spin(j)
               orbital_data(irr)%coefficients(:,sym_ind(irr)) = 0.0_wp
               orbital_data(irr)%occup(sym_ind(irr)) = dyson_orbital_norms(j,k)

               !For this Dyson orbital accumulate the contributions to the Dyson orbital's AO coefficients from all MO orbitals.
               do i=1,no_orbitals
                  if (dyson_orbitals(i,j,k) .eq. 0.0_cfp) cycle
                  do ao=1,n_cf
                     orbital_data(irr)%coefficients(ao,sym_ind(irr)) = orbital_data(irr)%coefficients(ao,sym_ind(irr)) &
                                                                        + cf(ao,i)*dyson_orbitals(i,j,k)
                  enddo !ao
               enddo !i

               !Include normalization of the continuum functions (needed for Molden) and the spin-related CG coefficient.
               do ao=1,n_cf
                  orbital_data(irr)%coefficients(ao,sym_ind(irr)) = orbital_data(irr)%coefficients(ao,sym_ind(irr)) &
                                                                    *normalization_factor(ao)*spin_cg
               enddo !ao

            enddo !j
         enddo !k
         !
         !OUTPUT THE DYSONS TO DISK:
         !
         dyson_basis_set_file = './dyson_orbitals.ukrmolp'
         molden_file_name = './dyson_orbitals.molden'
         write(cnum,'(i0)') myrank
         molden_file_name = trim(molden_file_name)//trim(adjustl(cnum))
         dyson_basis_set_file = trim(dyson_basis_set_file)//trim(adjustl(cnum))

         !Get the geometry and symmetry information from the GTO basis set
         call molecular_orbital_basis%ao_basis%symmetry_data%get_geometry(geometry)

         !Delete any already existing file with the target file name
         call delete_file(molden_file_name)

         !Write the CGTO basis and all Dyson orbitals into the Molden file
         call molden_file%init(molden_file_name,2)

         if (present(molden_format_opt)) then
            molden_format = molden_format_opt
         else
            molden_format = cart_basis_fmt
         end if

         if (molecular_orbital_basis%ao_basis%contains_btos()) then

            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
              'The Molden file will be written without BTOs since the basis contains BTOs which are not supported by Molden.', 3, 0)

            allocate(orbital_data_no_btos(size(orbital_data)),stat=err)
            if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Memory allocation 3 failed.',err,1)

            orbital_data_no_btos = orbital_data
            orbital_data_no_btos(:)%number_of_coefficients = molecular_orbital_basis%ao_basis%n_cgto_functions

            i = molecular_orbital_basis%ao_basis%n_cgto_functions
            do irr=1,molecular_orbital_basis%no_irr
               deallocate(orbital_data_no_btos(irr)%coefficients)
               allocate(orbital_data_no_btos(irr)%coefficients(i,size(orbital_data(irr)%coefficients,2)))
               orbital_data_no_btos(irr)%coefficients = 0.0_cfp
               orbital_data_no_btos(irr)%coefficients(1:i,:) = orbital_data(irr)%coefficients(1:i,:)
            enddo

            call molden_file%write(geometry%nucleus,CGTO_shells,orbital_data_no_btos,molden_format)
         else
            call molden_file%write(geometry%nucleus,CGTO_shells,orbital_data,molden_format)
         endif

         call molden_file%final()

         !Initialize the basis set of Dyson orbitals - this will store all Dyson orbital coefficients and the full GTO basis. We are doing this so that
         !we always keep all information about the Dyson orbitals (the Molden file may not be able to store information for all L in the basis).
         Dyson_orbital_basis%ao_basis => molecular_orbital_basis%ao_basis
         err = Dyson_orbital_basis%init(molecular_orbital_basis%no_irr,geometry)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'write_molden_dyson_orbitals', &
                         'Error initializing the Dyson orbital basis data structure.', 3, 1)
         end if

         do irr=1,molecular_orbital_basis%no_irr
            !Remove the normalization factor for the continuum functions from the orbital coefficients: that's only needed for Molden since it does not normalize to the R-matrix radius.
            do ao=1,n_cf
               orbital_data(irr)%coefficients(ao,:) = orbital_data(irr)%coefficients(ao,:)/normalization_factor(ao)
            enddo !ao
            !Add the Dyson orbitals for each symmetry to the full set of Dyson orbitals
            call Dyson_orbital_basis%add_shell(orbital_data(irr))
         enddo

         !Delete any already existing file with the same name as dyson_basis_set_file
         call delete_file(dyson_basis_set_file)

         !Finally, save the full set of Dysons and the AO basis set to a separate file that can be used for later processing
         call molecular_orbital_basis%ao_basis%write(dyson_basis_set_file)
         call Dyson_orbital_basis%write(dyson_basis_set_file)

         err = Dyson_orbital_basis%final()
         if (err .ne. 0) call xermsg('ukrmol_interface','write_molden_dyson_orbitals','Error finalizing Dyson_orbital_basis.',4,1)

   end subroutine write_molden_dyson_orbitals

!OUTER REGION (SWINTERF) INTERFACING ROUTINES:

   !> Calculates the orbital amplitudes and determines the channel information. The result is stored in the module private variables amplitudes,continuum_channels_m_l_irr for later use by UKP_PREAMP, UKP_READAMP.
   !> Call to this routine must be preceeded by calling e.g. READ_UKRMOLP_BASIS in order to initialize the molecular_orbital_basis object which performs the amplitude evaluation.
   !> \warning This routine ALWAYS requires double precision value of the R-matrix radius on input due to the need for compatibility with UKRmol-out which is always compiled using double precision.
   subroutine EVAL_AMPLITUDES(a,normalize_to_a)
      implicit none
      real(kind=wp), intent(in) :: a
      logical, intent(in) :: normalize_to_a
      real(kind=cfp) :: a_conv

      logical, allocatable :: is_continuum(:)
      integer :: err, irrcont

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'eval_amplitudes', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

!         if (a .eq. real(options%a,kind=wp)) then
!            normalize_to_a = .true.
!         else
!            normalize_to_a = .false.
!         endif

         a_conv = a
         call molecular_orbital_basis%calculate_amplitudes(a_conv,normalize_to_a,amplitudes,continuum_channels_m_l_irr)

         if (.not.allocated(ncont)) then
            allocate(ncont(molecular_orbital_basis%no_irr),stat=err)
            if (err .ne. 0) call xermsg('ukrmol_interface','EVAL_AMPLITUDES','Memory allocation failed.',err,1)
            ncont = 0
         endif

         do irrcont=1,molecular_orbital_basis%no_irr
            call molecular_orbital_basis%get_continuum_flags(irrcont,is_continuum)
            ncont(irrcont) = count(is_continuum)

            write(level2,'(/,10X,"IRR = ",i2," number of continuum orbitals = ",i0)') irrcont, ncont(irrcont)
         enddo !irrcont

   end subroutine EVAL_AMPLITUDES

   !> Transfers the channel angular and IRR numbers to the channel variables from SWINTERF. It is the equivalent of the PREAMP routine from SWINTERF. The call to this routine must be preceeded by call to
   !> EVAL_AMPLITUDES.
   subroutine UKP_PREAMP(IRRCONT,IX,LCHL,MCHL,QCHL,NCH,MAXNMO)
      implicit none
      integer, intent(in) :: IRRCONT,IX
      integer, intent(out) :: LCHL(*),MCHL(*),QCHL(*),NCH
      integer, intent(inout) :: MAXNMO

      integer :: lm, nch_lm, irr

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'UKP_PREAMP', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         nch_lm = size(continuum_channels_m_l_irr,2)

         nch = 0
         do lm=1,nch_lm
            irr = continuum_channels_m_l_irr(3,lm)
            if (irr .eq. irrcont) then
               nch = nch + 1
               LCHL(IX+nch-1) = continuum_channels_m_l_irr(2,lm)
               MCHL(IX+nch-1) = abs(continuum_channels_m_l_irr(1,lm))
               QCHL(IX+nch-1) = sign(1,continuum_channels_m_l_irr(1,lm))
               if (MCHL(IX+nch-1) .eq. 0) QCHL(IX+nch-1) = 0
            endif
         enddo

         maxnmo = max(maxnmo,ncont(irrcont))

   end subroutine UKP_PREAMP

   !> This routine transfers into WAMPS the R-matrix amplitudes for the orbitals. It is equivalent of the READAMP routine. Note that the orbital amplitudes are divided by sqrt(2). 
   !> See the routine outerio/READRM routine and the equation for the R-matrix there. It implies that the 1/2 factor from the R-matrix formula is absorbed into the product of the two amplitudes 
   !> in the summation over the R-matrix poles. Hence the factor 1/sqrt(2) by which the amplitudes must be multiplied.
   !> The call to this routine must be preceeded by call to UKP_PREAMP.
   !> \warning This routine ALWAYS outputs double precision values of the amplitudes in WAMPS due to the need for compatibility with UKRmol-in/out which is always compiled using double precision.
   subroutine UKP_READAMP(WAMPS,NCHAN,IRRCHL,LCHL,MCHL,QCHL,NCONTCSF,MCONT,IPRINT)
      implicit none
      integer :: IRRCHL(*), LCHL(*), MCHL(*), QCHL(*), NCONTCSF(8), MCONT(:), NCHAN, IPRINT
      real(kind=wp) :: WAMPS(NCHAN,*)

      integer :: i, l, lm, nch_lm, orbs_start, orbs_end, nfound, nunique, ntgtsym
      logical :: match

         if (.not. molecular_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'UKP_READAMP', 'The orbital basis set data has not been initialized.', 1, 1)
         end if

         if (iprint > 0) WRITE(stdout,1000)
         if (iprint > 1) WRITE(stdout,1010) NCHAN,(I,IRRCHL(I),I=1,NCHAN)

         nch_lm = size(continuum_channels_m_l_irr,2) !number of different (l,m) channel numbers
         ntgtsym = size(mcont)    !number of symmetries of the target states

         if (iprint > 0) then
            do i=1,molecular_orbital_basis%no_irr
               nunique = 0
               do lm=1,nch_lm
                  if (continuum_channels_m_l_irr(3,lm) .eq. i) nunique = nunique + 1
               enddo
               WRITE(stdout,2010) i,molecular_orbital_basis%ao_basis%number_of_functions,num_orbs_sym(i),nunique
            enddo
         endif

         !loop over all channels in the problem
         nfound = 0
         do L=1,NCHAN

            if (IRRCHL(L) > molecular_orbital_basis%no_irr .or. IRRCHL(L) <= 0) then
                call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                             'IRRCHL value is out of range and is incompatible with the orbital basis data.', 2, 1)
            end if

            !make sure that for each target state symmetry there really is a continuum channel with symmetry that matches the contiuum symmetry given by MCONT (from SCATCI).
            match = .false.
            do i=1,ntgtsym
               if (IRRCHL(L) .eq. mcont(i)+1) match = .true.
            enddo
            if (.not. match) then
                call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                             'Mismatch between orbital basis and CI vector data: &
                             &symmetry of at least one continuum channel lies outside of the set of continuum symmetries &
                             &given by the MCONT values obtained from CI (SCATCI) data.', 3, 1)
            end if

            !loop over all channels from the boundary amplitude data and look for the matching channel data
            do lm=1,nch_lm
               if (IRRCHL(L) == continuum_channels_m_l_irr(3,lm) .and. &
                   MCHL(L)*QCHL(L) == continuum_channels_m_l_irr(1,lm) .and. &
                   LCHL(L) == continuum_channels_m_l_irr(2,lm)) then
                  nfound = nfound + 1
                  orbs_end = sum(num_orbs_sym(1:IRRCHL(L)))
                  orbs_start = orbs_end-ncont(IRRCHL(L))+1 !index of the first continuum orbital in symmetry IRRCHL(L).
                  i = sum( num_orbs_sym( 1:(IRRCHL(L)-1) ) )+1 !index of the first orbital in symmetry IRRCHL(L).
                  if (orbs_start < i) then
                     print *,orbs_start,i
                     call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                                  'Wrong number of continuum orbitals or an error in orbital data.', 3, 1)
                  endif
                  WAMPS(L,1:ncont(IRRCHL(L))) = amplitudes(lm,orbs_start:orbs_end)/sqrt(2.0_cfp)
                  if (iprint > 0) then
                      write(stdout,2810) L, continuum_channels_m_l_irr(3,lm), continuum_channels_m_l_irr(2,lm), &
                                            continuum_channels_m_l_irr(1,lm), (WAMPS(L,I),I=1,ncont(IRRCHL(L)))
                  end if
                  exit
               endif
            enddo !lm

         enddo !L

         ncontcsf = 0
         ncontcsf(1:molecular_orbital_basis%no_irr) = ncont(1:molecular_orbital_basis%no_irr)

         if (nfound /= nchan) then
            call xermsg ('ukrmol_interface', 'UKP_READAMP', &
                         'The number of channels for which boundary amplitudes have been transfered &
                         &does not match the actual number of channels for which the amplitudes are required.', 4, 1)
         end if

 1000 FORMAT(/,10X,'====> UKP_READAMP - READ BOUNDARY AMPS <====',/)
 1010 FORMAT(/,10X,'Searching for data on ',I3,' channels ',//,&
     &         10X,'Seq No.  IRRCHL',/,&
     &         10X,'-------  ------',//,&
     &        (10X,I6,2X,I4))
 2010 FORMAT(/,10X,'Symmetry number = ',I3,//,10X,'No. of basis functions = ',I3,/,&
               10X,'No. of orbitals        = ',I3,/,10X,'No. of angular parts   = ',I3,/)
 2810 FORMAT(/1X,'Chan No. ',I2,' Symmetry = ',I2,' Angular behaviour (l,m)',' = ',2I5,/,(1X,5(F12.8,1X)))

   end subroutine UKP_READAMP

   subroutine GET_ECP_CORE_CHARGES(core_charges)
      implicit none
      integer, allocatable :: core_charges(:)

      integer :: i

         if (.not. atomic_orbital_basis % is_initialized()) then
            call xermsg ('ukrmol_interface', 'GET_ECP_CORE_CHARGES', 'The atomic basis set data has not been initialized.', 1, 1)
         end if

         allocate(core_charges(atomic_orbital_basis%symmetry_data%no_nuc))
         core_charges = 0
 
         do i = 1, atomic_orbital_basis%symmetry_data%no_nuc
            core_charges(i) = atomic_orbital_basis%symmetry_data%nucleus(i)%ecp%n_core
            if (core_charges(i) /= 0) then
               write(stdout,'(" ECP found on atom ",i2,1X,a2)') atomic_orbital_basis%symmetry_data%nucleus(i)%nuc, &
                                                                atomic_orbital_basis%symmetry_data%nucleus(i)%name
            endif
         enddo

   end subroutine GET_ECP_CORE_CHARGES


   !>****************************************************************
   !>                                                               *
   !>    ROUTINES PROVIDING SIMPLIFIED INTERFACE FOR ADC CODES      *
   !>                                                               *
   !>****************************************************************

   !> This routine reads-in all two electron integrals needed for ADC codes. It is derived from the READ_UKRMOLP_INTS -- parameter interface is significantly simplified and redundant/UKRmol-specific functionality is removed.
   !> no positrons are expected, corresponding switch parameter (iposit) removed
   !> \warning This routine ALWAYS requires double precision value of scalem on input due to the need for compatibility with UKRmol-in/out which is always compiled using double precision. Note that
   !> this results in reducing the values in ke_integrals%a into double precision.
   !> NEW FEATURE: WORKS WITH SHARED MEMORY OF THE MPI-3.0 STANDARD. This means that all integrals are read-in the shared mode so in the MPI regime there is only one copy of the integral arrays PER NODE.
   !> If we do not have the shared memory capability then the code resorts to the local mode, i.e. each MPI task keeps its own copy of all integral arrays.
   subroutine READ_ADC_INTS(nfti,nfta,nint1e,nint2e,nsym)
      use symmetry_gbl
      use const_gbl, only: set_verbosity_level, Fock_blocks_header, &
        level1, orbs_line
      implicit none
      integer, intent(in) :: nfti,nfta
      integer(kind=shortint), intent(inout) :: nsym
      integer, intent(inout) :: nint1e, nint2e

      !local variables
      character(len=line_len) :: file_name
      type(integral_storage_obj) :: mo_integral_storage
      type(integral_storage_obj), target :: tgt_storage
      type(data_header_obj) :: header
      type(data_file_obj) :: data_file
      integer :: i, j, err, nnuc, lunit, first_record, last_record
      logical, parameter :: tgt_is_local = .true.

      integer irr,n_target,n_continuum,k,m,nbas,imin,imax,jmin,jmax
      type(orbital_data_obj) :: orbital_data
      logical, allocatable :: is_continuum(:)
      real(kind=cfp), allocatable :: TT_block(:,:), TC_block(:,:), CC_block(:,:)

         write(nfta,'(/,"READ_ADC_INTS called.")')
         !we assume that the integrals are saved on the unit number 'nfti'
         call READ_UKRMOLP_BASIS(NFTI)

         !check for data input errors
         if (molecular_orbital_basis%no_irr.ne.nsym) then !nsym is the number of the last IRR where we have some orbitals
            write(nfta,'(/,"WARNING: The number of IRRs obtained from the MO")')
            write(nfta,'("    basis data not compatible with the input.")')
            write(nfta,*) molecular_orbital_basis%no_irr, nsym
            nsym = molecular_orbital_basis%no_irr
            write(nfta,'("    Using the value from actual MO basis:",I2)') nsym
!             call xermsg ('ukrmol_interface', 'READ_ADC_INTS', &
!                          'The number of IRRs obtained from the orbital basis data is not compatible &
!                          &with the one given by ADC calling program.', 2, 1)
         endif

         !store the MO basis information


         !open the data file and analyse the headers manually
         call unit_number_2_file_name(NFTI,file_name)
         err = mo_integral_storage%init(disk=file_name)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'READ_ADC_INTS', 'Error initializing the integrals storage object.', err, 1)
         end if

         write(nfta,1600) nint1e, nint2e
         write(nfta,'(/," ... 1-e integrals will not be loaded")')
 1600    format(/,' On input ADC expects:',/,' 1-electron integrals, NINT1e = ',i0,/,' 2-electron integrals, NINT2e = ',i0)

          !read-in the 1-electron integrals into one_electron_integrals
!          err = tgt_storage%init(memory=one_electron_integrals)
!          tgt_storage%integral_file%identifier = mo_integral_storage%integral_file%identifier
!          if (err .ne. 0) call xermsg('ukrmol_interface','READ_ADC_INTS','Error initializing tgt_storage.',err,1)
!
!          !we need to get the full header corresponding to the transformed overlap and kinetic energy integrals.
!          err = mo_integral_storage % integral_file % get_header_containing(header, &
!                                                                            molecular_orbital_basis % get_basis_name(), &
!                                                                            one_p_sym_ints)
!          if (err .ne. 0) call xermsg('ukrmol_interface','READ_ADC_INTS','Error locating the transformed 1-electron integrals; &
!                                      &see data_header_obj%get_header_containing for details.',err,1)
!
!          call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
!          call tgt_storage%final
!
!          write(nfta,'(/,"One electron integrals read-in.")')
!
!          !Column number in one_electron_integrals%a corresponding to the 1-electron Hamiltonian integrals: these integrals must ALWAYS be calculated.
!          one_elham_column = one_electron_integrals%find_column_matching_name(one_elham)

         !read-in the 2p integrals
         err = tgt_storage%init(memory=two_electron_integrals)
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_ADC_INTS','Error initializing tgt_storage.',err,1)

         !we need to get the full header corresponding to the transformed 2-electron integrals.
         err = mo_integral_storage % integral_file % get_header_containing(header, &
                                                                           molecular_orbital_basis % get_basis_name(), &
                                                                           two_p_sym_ints)
         if (err .ne. 0) call xermsg('ukrmol_interface','READ_ADC_INTS','Error locating the transformed 2-electron integrals; &
                                     &see data_header_obj%get_header_containing for details.',err,1)

         call tgt_storage%read(mo_integral_storage,header%name,options,tgt_is_local)
         call tgt_storage%final

         write(nfta,'(/,"Two electron integrals read-in.")')

         !Now update the integral counters by the actual number of integrals read-in
         nint1e = size(one_electron_integrals%a,1)
         nint2e = size(two_electron_integrals%a,1)
         write(nfta,1700) nint1e, nint2e
 1700    format(/,' UKRMol+ integrals read successfully:',&
                /,' 1-electron integrals, NINT1e =',i0,&
                /,' 2-electron integrals, NINT2e =',i0)

         call mo_integral_storage%final

         write(nfta,'(/,10x,"Finished reading GBTOlib integrals.")')

     ! READING AND PRINTING THE ORBITAL ENERGIES AND THE FOCK MATRIX BLOCKS
         call set_verbosity_level(3)

         call data_file%open(file_name)
         lunit = data_file%get_unit_no()

         err = data_file%find_header(Fock_blocks_header,first_record,last_record)
         if (err /= 0) then
            call xermsg ('ukrmol_interface', 'READ_ADC_INTS', &
                         'Error locating the header corresponding to the Fock matrix blocks.', err, 1)
         end if

         if (first_record <= 0 .or. last_record <= 0) then
            call xermsg ('ukrmol_interface', 'READ_ADC_INTS', &
                         'The requested data record is missing or not complete.', 2, 1)
         end if
         last_record = first_record

         ! ALLOCATE THE FOCK MATRIX ARRAY
         nbas = 0
         do irr = 1,molecular_orbital_basis%no_irr
           nbas = nbas+molecular_orbital_basis%get_number_of_orbitals(irr)
         end do
         allocate(fock(nbas,nbas))
         fock = 0.0
         write(level1,'(/,3X,"Number of MOs:",i0,/)') nbas

         do irr = 1,molecular_orbital_basis%no_irr

           write(level1,'(/,3X,"Symmetry: ",i2,/)') irr
           call molecular_orbital_basis%get_shell_data(irr,orbital_data)
           allocate(is_continuum(orbital_data%number_of_functions))
           call molecular_orbital_basis%get_continuum_flags(irr,is_continuum)
           do i = 1, orbital_data%number_of_functions
              k = molecular_orbital_basis%get_absolute_index(i,irr)
              if (is_continuum(i)) then
                 write(level1,'(5X,i3,5X,i3,2X," Continuum, energy: ",e25.15)') i, k, orbital_data%energy(i)
              else
                 write(level1,'(5X,i3,5X,i3,2X," Target,    energy: ",e25.15)') i, k, orbital_data%energy(i)
              endif
           enddo !i

           read(lunit,pos=last_record,err=10) n_target, n_continuum
           write(level1,*)
           write(level1,*) "n_targ/n_cont:",n_target,n_continuum
           inquire(lunit,pos=last_record)

           if (n_target > 0) then
             allocate(TT_block(n_target,n_target))
             read(lunit,pos=last_record,err=10) TT_block
             inquire(lunit,pos=last_record)

             imin = molecular_orbital_basis%get_absolute_index(1,irr)
             imax = molecular_orbital_basis%get_absolute_index(n_target,irr)

             fock(imin:imax,imin:imax) = TT_block(1:n_target,1:n_target)

             ! Writing the block into the log file
             write(level3,'(/,2X,"Target-target block of Fock matrix")')
             write(level3,*) imin,imax
             k = 0
             do i=1,n_target/orbs_line
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
               do j=1,n_target
                  write(level3,'(i0,50e25.15)') j,TT_block(j,k+1:k+orbs_line)
               enddo !j
               k = k + orbs_line
             enddo !i

             m = mod(n_target,orbs_line)
             if (m > 0) then
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
               do j=1,n_target
                 write(level3,'(i0,50e25.15)') j,TT_block(j,k+1:k+m)
               enddo !j
             endif
             ! End of writing the block into the log file

           endif !n_target>0

           if (n_target > 0 .and. n_continuum > 0) then
             allocate(TC_block(n_target,n_continuum))
             read(lunit,pos=last_record,err=10) TC_block
             inquire(lunit,pos=last_record)

             imin = molecular_orbital_basis%get_absolute_index(1,irr)
             imax = molecular_orbital_basis%get_absolute_index(n_target,irr)

             jmin = molecular_orbital_basis%get_absolute_index(n_target+1,irr)
             jmax = molecular_orbital_basis%get_absolute_index(n_target+n_continuum,irr)

             fock(imin:imax,jmin:jmax) = TC_block(1:n_target,1:n_continuum)
             fock(jmin:jmax,imin:imax) = transpose(TC_block(1:n_target,1:n_continuum))

             ! Writing the block into the log file
             write(level3,'(/,2X,"Continuum-target block of Fock matrix (a transpose of the saved block)")')
             write(level3,*) imin,imax
             write(level3,*) jmin,jmax
             k = 0
             do i=1,n_target/orbs_line
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
               do j=1,n_continuum
                  write(level3,'(i0,50e25.15)') j+n_target,TC_block(k+1:k+orbs_line,j)
               enddo !j
               k = k + orbs_line
             enddo !i

             m = mod(n_target,orbs_line)
             if (m > 0) then
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
               do j=1,n_continuum
                  write(level3,'(i0,50e25.15)') j+n_target,TC_block(k+1:k+m,j)
               enddo !j
             endif
            ! End of writing the block into the log file
           end if !n_target > 0 .and. n_continuum > 0

           if (n_continuum > 0) then
              allocate(CC_block(n_continuum,n_continuum))
              read(lunit,pos=last_record,err=10) CC_block
              inquire(lunit,pos=last_record)

              jmin = molecular_orbital_basis%get_absolute_index(n_target+1,irr)
              jmax = molecular_orbital_basis%get_absolute_index(n_target+n_continuum,irr)

              fock(jmin:jmax,jmin:jmax) = CC_block(1:n_continuum,1:n_continuum)

              ! Writing the block into the log file
              write(level3,'(/,2X,"Continuum-continuum block of Fock matrix")')
              write(level3,*) jmin,jmax
              k = 0
              do i=1,n_continuum/orbs_line
                 write(level3,'(/,10X,50(i0,2X))') (j+n_target,j=k+1,k+orbs_line)
                 do j=1,n_continuum
                    write(level3,'(i0,50e25.15)') j+n_target,CC_block(k+1:k+orbs_line,j)
                 enddo !j
                 k = k + orbs_line
              enddo !i

              m = mod(n_continuum,orbs_line)
              if (m > 0) then
                 write(level3,'(/,10X,50(i0,2X))') (j+n_target,j=k+1,k+m)
                 do j=1,n_continuum
                    write(level3,'(i0,50e25.15)') j+n_target,CC_block(k+1:k+m,j)
                 enddo !j
              endif

              write(level3,*)
              k = 0
              do i = 1,(jmax-jmin+1)/orbs_line
                write(level3,'(/,10X,50(i0,2X))') (j+jmin-1,j=k+1,k+orbs_line)
                do j = jmin,jmax
                  write(level3,'(i0,50e25.15)') j,fock(j,k+jmin:k+jmin-1+orbs_line)
                end do
                k = k+orbs_line
              end do
              ! End of writing the block into the log file
           end if

           deallocate(is_continuum)
           if (allocated(TT_block)) deallocate(TT_block)
           if (allocated(TC_block)) deallocate(TC_block)
           if (allocated(CC_block)) deallocate(CC_block)

         end do !irr

         call data_file%close()

         return

10       call xermsg ('ukrmol_interface', 'READ_ADC_INTS', 'Error reading Fock matrix.', 4, 2)

   end subroutine READ_ADC_INTS

   !> Returns molecular orbital data needed for the construction of ADC
   !> Hamiltonian
   subroutine GET_ENERSYMOCC(pgroup,nsym,nbas,nocc,ener,sym,occ,iscont,mt)
     integer(kind=shortint),intent(in) :: nsym
     integer(kind=shortint),intent(out) :: pgroup,nbas,nocc
     real(kind=wp),allocatable,intent(out) :: ener(:),occ(:)
     integer(kind=shortint),allocatable,intent(out) :: sym(:)
     logical,allocatable,intent(out) :: iscont(:)
     integer(kind=shortint),dimension(8,8),intent(out) :: mt

     type(orbital_data_obj),allocatable :: orbital_data(:)
     logical,allocatable :: contflag(:)
     integer i,j,ind

     allocate(orbital_data(nsym))

     nbas = 0
     do i = 1,nsym
       nbas = nbas+molecular_orbital_basis%get_number_of_orbitals(i)
       call molecular_orbital_basis%get_shell_data(i,orbital_data(i))
     end do
     pgroup = int(orbital_data(1)%point_group,kind=shortint)

     allocate(sym(nbas),ener(nbas),occ(nbas),iscont(nbas))
     nocc = 0
     do i = 1,nsym
       allocate(contflag(molecular_orbital_basis%get_number_of_orbitals(i)))
       call molecular_orbital_basis%get_continuum_flags(i,contflag)
       do j = 1,molecular_orbital_basis%get_number_of_orbitals(i)
         ind = molecular_orbital_basis%get_absolute_index(j,i)
         ener(ind) = real(orbital_data(i)%energy(j),wp)
         occ(ind) = real(orbital_data(i)%occup(j),wp)
         if (occ(ind).eq.2) nocc = nocc+1
         iscont(ind) = contflag(j)
         sym(ind) = i
       end do
       deallocate(contflag)
     end do
     deallocate(orbital_data)

     mt(1:8,1:8) = int(abel_prod_tab(1:8,1:8),shortint)

   end subroutine GET_ENERSYMOCC

   !> Returns double-precision two-electron integral given the four INTEGER(2)
   !> indices p,q,r,s -- required by the ADC codes; one-electron integrals or
   !> positron flag not supported
   function LPQRS(p,q,r,s)
     implicit none
     integer(kind=shortint2),intent(in) :: p,q,r,s
     real(kind=wp) :: LPQRS

     integer :: ind(1), four_ind(1:4,1)

        four_ind(1,1) = p
        four_ind(2,1) = q
        four_ind(3,1) = r
        four_ind(4,1) = s
        ind(1:1) = molecular_orbital_basis%integral_index(two_el_ints,four_ind,options%two_p_continuum)
        if (ind(1).gt.0) then
          LPQRS = two_electron_integrals%a(ind(1),1)
        else
          LPQRS = 0.0_wp
        end if
  end function LPQRS

end module ukrmol_interface_gbl
