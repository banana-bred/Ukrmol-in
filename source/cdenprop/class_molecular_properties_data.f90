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

!> \brief   Property storage
!> \authors A Harvey, Z Masin, J Benda
!> \date    2011 - 2019
!>
!> This module contains the type \ref molecular_properties_data, which is used by cdenprop_target
!> and mpi_scatci for storage and operation with properties extracted by cdenprop.
!>
module class_molecular_properties_data

    use cdenprop_defs, only: CIvect, idp, CSFheader, property_integrals, maxnuc, tsym_max, ir_max, &
                             amass, asymb, cc2v_symb, cd2h_symb, spin_symb, &
                             small_swintf_no_states, large_swintf_prop_states_nuclei, &
                             large_swintf_no_nuclei, large_swintf_no_states, large_swintf_non_zero_properties
    use precisn,       only: wp

    implicit none

    private

    !> \brief   Property storage
    !> \authors A Harvey, Z Masin, J Benda
    !> \date    2011 - 2019
    !>
    !> This data structure is (optionally) used for storing properties extracted by cdenprop; currently, it is used
    !> in cdenprop_target and mpi-scatci. It offers two independent way of storing the properties:
    !>  1. As a set of dense matrices, one per property and bra and ket IRR.
    !>  2. As a single list of non-zero properties, where all properties, bra and ket states are merged.
    !>
    type, public :: molecular_properties_data

        integer :: no_nuclei = 0            !< Number of nuclei forming the molecule
        integer :: no_symmetries = 0        !< Total number of distinct spin-symmetries present in the storage
        integer :: no_states = 0            !< Total number of states in all spin-symmetries
        integer :: no_blocks = 0            !< Number of dense blocks in the dense property storage
        integer :: non_zero_properties = 0  !< Number of non-zero properties in the sparse storage
        integer :: nuclear_prop_flag = 0    !< Is nuclear part of the property operator considered (= 0) or not (= 1)?

        logical :: swintf_format = .false.  !< Controls whether to use SWINTERF-compatible property file format
        
        logical :: contains_continuum = .false.  !< Defines whether the property csf(s) contain any continuum 
                                                 !< wavefunctions (can be used to differentiate N and N+1 calculations).

        type (property_integrals),pointer :: pintegrals => null()

        character(len=8), dimension(:) :: cname(maxnuc)
        real(kind=idp),   dimension(:) :: charge(maxnuc), mass(maxnuc), xnuc(maxnuc), ynuc(maxnuc), znuc(maxnuc)

        !> List of data for all spin-symmetries present in the property set
        !>  1. irreducible representation
        !>  2. spin multiplicity
        !>  3. number of states
        integer :: symmetries_index(3, tsym_max) = 0

        !> Metadata for states:
        !>  1. index of the state within the spin-symmetry
        !>  2. index of the spin-symmetry
        integer, allocatable :: energies_index(:,:)

        !> Energies of states, in the same order as in \ref energies_index.
        real(kind=idp), allocatable :: energies(:)

        !> Set of dense property matrices \f$ \langle I | Y[\ell,m] | J \rangle \f$.
        type(CIvect), allocatable :: dense_blocks(:)

        !> Metadata for all dense property matrices
        !>  1. bra spin-symmetry index
        !>  2. ket spin-symmetry index
        !>  3. angular momentum
        !>  4. its projection multiplied by projectile charge
        integer, allocatable :: dense_index(:,:)

        !> Property data:
        !>  1. bra state
        !>  2. ket state
        !>  3. angular momentum
        !>  4. its projection multiplied by projectile charge
        integer, allocatable :: properties_index(:,:)

        !> Property values, in the same order as in \ref properties_index.
        real(kind=idp), allocatable :: properties(:)

    contains

        ! subroutines operating on sparse storage
        procedure, public  :: read_properties
        procedure, public  :: decompress
        procedure, public  :: calculate_polarisability

        ! subroutines operating on dense storage
        procedure, public  :: preallocate_property_blocks
        procedure, public  :: insert_property_block

        ! subroutines operating on both
        procedure, public  :: sort_by_energy
        procedure, public  :: write_properties
        procedure, private :: append_states
        procedure, public  :: clean

        ! Finalisation routine
        final :: clean_final

    end type molecular_properties_data

contains

    !> \brief   Release memory used by the object
    !> \authors J Benda
    !> \date    2019
    !>
    !> Deallocates all allocatable attributes, resets all counters to zero.
    !>
    !> \param[inout] this     Property object to update.
    !>
    subroutine clean (this)

        class(molecular_properties_data), intent(inout) :: this

        this % no_nuclei = 0
        this % no_blocks = 0
        this % no_states = 0
        this % no_symmetries = 0
        this % non_zero_properties = 0
        this % symmetries_index = 0

        if (allocated(this % energies))         deallocate (this % energies)
        if (allocated(this % energies_index))   deallocate (this % energies_index)
        if (allocated(this % dense_blocks))     deallocate (this % dense_blocks)
        if (allocated(this % dense_index))      deallocate (this % dense_index)
        if (allocated(this % properties_index)) deallocate (this % properties_index)
        if (allocated(this % properties))       deallocate (this % properties)

        if (associated(this % pintegrals)) then
           deallocate (this % pintegrals)
           nullify (this % pintegrals)
        end if

    end subroutine clean


    subroutine clean_final (this)

       type(molecular_properties_data), intent(inout) :: this

        if (allocated(this % energies))         deallocate (this % energies)
        if (allocated(this % energies_index))   deallocate (this % energies_index)
        if (allocated(this % dense_blocks))     deallocate (this % dense_blocks)
        if (allocated(this % dense_index))      deallocate (this % dense_index)
        if (allocated(this % properties_index)) deallocate (this % properties_index)
        if (allocated(this % properties))       deallocate (this % properties)

        if (associated(this % pintegrals)) then
           deallocate (this % pintegrals)
           nullify (this % pintegrals)
        end if

    end subroutine clean_final


    !> \brief   Prepare the dense storage
    !> \authors J Benda
    !> \date    2019
    !>
    !> To allow seamless addition of dense property blocks, the array of blocks is pre-allocated in this subroutine
    !> to the requested length. This call destroys any previously stored dense blocks.
    !>
    !> \param[inout] this     Property object to update.
    !> \param[in]    nblocks  How many blocks to allocate.
    !>
    subroutine preallocate_property_blocks (this, nblocks)

        class(molecular_properties_data), intent(inout) :: this
        integer, intent(in) :: nblocks

        if (allocated(this % dense_blocks)) deallocate (this % dense_blocks)
        if (allocated(this % dense_index))  deallocate (this % dense_index)

        this % no_blocks = 0

        allocate (this % dense_blocks(nblocks), this % dense_index(6, nblocks))

    end subroutine preallocate_property_blocks


    !> \brief   Append a new dense property block to the dense property storage
    !> \authors J Benda
    !> \date    2019
    !>
    !> This subroutine checks whether there are any non-zero properties in the given (potentially distributed) property block
    !> and if yes, assuming that there are enough pre-allocated dense property block slots, it appends the block to the
    !> preallocated storage.
    !>
    !> New states are added to the index if not already present. This subroutine requires that the states are groupped by their
    !> spin symmetry. Calling this subroutine after sorting the states by energiy will not work properly!
    !>
    !> This subroutine does not affect the sparse storage.
    !>
    !> \param[inout] this        Property object to update.
    !> \param[in]    csf_head_i  Bra eigenvector metadata (mgvn, spin).
    !> \param[in]    ci_vec_i    Bra eigenvectors and eigenvalues.
    !> \param[in]    nstat_i     Number of bra eigenvectors.
    !> \param[in]    csf_head_j  Ket eigenvector metadata (mgvn, spin).
    !> \param[in]    ci_vec_j    Ket eigenvectors and eigenvalues.
    !> \param[in]    nstat_j     Number of ket eigenvectors.
    !> \param[in]    threshold   Minimal magnitude of a property to be considered non-zero.
    !> \param[in]    l           Property angular momentum.
    !> \param[in]    mq          Property angular momentum projections.
    !> \param[in]    props       Distributed dense matrix of properties.
    !>
    subroutine insert_property_block (this, csf_head_i, ci_vec_i, nstat_i, csf_head_j, ci_vec_j, nstat_j, threshold, l, mq, props)

        use mpi_gbl, only: nprocs, mpiint, mpi_mod_allgather

        class(molecular_properties_data), intent(inout) :: this

        type (CSFheader), intent(in) :: csf_head_i, csf_head_j
        type (CIvect),    intent(in) :: ci_vec_i, ci_vec_j, props

        integer,  intent(in) :: nstat_i, nstat_j, l, mq
        real(wp), intent(in) :: threshold

        integer(mpiint) :: nnz, nprops(nprocs)
        integer         :: i, mgvn_i, mgvn_j, spin_i, spin_j, spinsym_i, spinsym_j

        ! find out the number of non-zero properties
        nnz = count(abs(props % CV) >= threshold)
        call mpi_mod_allgather(nnz, nprops)
        if (csf_head_i%ntgsym /= -1 .or. csf_head_j%ntgsym /= -1) then
            if (sum(int(nprops)) == 0) return
        end if

        this % no_nuclei = ci_vec_i % nnuc
        this % non_zero_properties = 0

        this % cname  = ci_vec_i % cname
        this % charge = ci_vec_i % charge
        this % xnuc   = ci_vec_i % xnuc
        this % ynuc   = ci_vec_i % ynuc
        this % znuc   = ci_vec_i % znuc
        this % mass   = 0

        mgvn_i = csf_head_i % mgvn
        mgvn_j = csf_head_j % mgvn

        spin_i = nint(2 * csf_head_i % S + 1)
        spin_j = nint(2 * csf_head_j % S + 1)

        ! get atomic masses from the CDENPROP database
        do i = 1, this % no_nuclei
            if (1 <= nint(ci_vec_i % charge(i)) .and. nint(ci_vec_i % charge(i)) <= size(AMASS)) then
                this % mass(i) = AMASS(nint(ci_vec_i % charge(i)))
            end if
        end do

        ! add the bra/ket spin-symmetry states (if not already present) and get offsets in energies_index
        spinsym_i = this % append_states(mgvn_i, spin_i, nstat_i, ci_vec_i % ei, ci_vec_i % e0)
        spinsym_j = this % append_states(mgvn_j, spin_j, nstat_j, ci_vec_j % ei, ci_vec_j % e0)

        ! store the properties
        this % no_blocks = this % no_blocks + 1
        this % dense_blocks(this % no_blocks) = props
        this % dense_index(1, this % no_blocks) = spinsym_i
        this % dense_index(2, this % no_blocks) = spinsym_j
        this % dense_index(3, this % no_blocks) = l
        this % dense_index(4, this % no_blocks) = mq

        ! erase small properties
        where (abs(this % dense_blocks(this % no_blocks) % CV) < threshold)
            this % dense_blocks(this % no_blocks) % CV = 0
        end where

    end subroutine insert_property_block


    !> \brief   Add states to index
    !> \authors J Benda
    !> \date    2019
    !>
    !> For the given spin-symmetry, add given number of states to the energies index. When the spin-symmetry is already
    !> fully present in the index, do nothing. If this spin-symmetry is partially present (this happens in CDENPROP, where
    !> there is one bra state, but many ket states), all states are added (with unique indices) even if they duplicate some
    !> already existing ones.
    !>
    !> \param[in] mgvn   Irreducible representation (zero-based index).
    !> \param[in] spin   Spin multiplicity.
    !> \param[in] nstat  Number of states.
    !> \param[in] ei     Energies of the states (relative to e0).
    !> \param[in] e0     Base energy.
    !>
    integer function append_states (this, mgvn, spin, nstat, ei, e0) result (isym)

        use mpi_gbl, only: mpi_xermsg

        class(molecular_properties_data), intent(inout) :: this
        real(wp), allocatable,            intent(in)    :: ei(:)

        integer,  intent(in) :: mgvn, spin, nstat
        real(wp), intent(in) :: e0

        real(wp), allocatable :: energies(:)
        integer,  allocatable :: energies_index(:,:)

        integer :: istate

        ! exit immediately if this spin-symmetry is already present
        do isym = 1, this % no_symmetries
            if (this % symmetries_index(1, isym) == mgvn .and. &
                this % symmetries_index(2, isym) == spin .and. &
                this % symmetries_index(3, isym) == nstat) return
        end do

        ! expanding symmetry index is not implemented
        if (this % no_symmetries == tsym_max) then
            call mpi_xermsg('class_molecular_properties_data', 'append_states', &
                            'Cannot add next symmetry - out of space.', 1, 1)
        end if

        ! expand the state index
        call move_alloc(this % energies, energies)
        call move_alloc(this % energies_index, energies_index)
        allocate (this % energies_index(2, this % no_states + nstat), &
                  this % energies(this % no_states + nstat))
        if (this % no_states > 0) then
            this % energies_index(:, 1:this % no_states) = energies_index(:, 1:this % no_states)
            this % energies(1:this % no_states) = energies(1:this % no_states)
        end if

        ! append new states
        do istate = 1, nstat
            this % energies_index(1, this % no_states + istate) = istate
            this % energies_index(2, this % no_states + istate) = isym
            this % energies(this % no_states + istate) = e0 + ei(istate)
        end do

        ! update metadata
        this % no_states = this % no_states + nstat
        this % no_symmetries = this % no_symmetries + 1
        this % symmetries_index(1, this % no_symmetries) = mgvn
        this % symmetries_index(2, this % no_symmetries) = spin
        this % symmetries_index(3, this % no_symmetries) = nstat

    end function append_states


    !> \brief   Read property file
    !> \authors J Benda
    !> \date    2019 - 2021
    !>
    !> Reads atomic, targets and properties from the selected set in given file unit.
    !> Assumes FORMATTED file. Erases all stored data before reading the file.
    !> Adapted from SWTARG (UKRmol-out).
    !>
    !> \note This subroutine reads the data as they are, i.e. if only upper triangle of each istate-jstate matrix is stored,
    !>       this objects will contain only those. This also means that if some property occurs multiple times, this subroutine
    !>       will neither discover nor correct it. To symmetrize the properties and check possible duplicates, use the subroutine
    !>       \ref symmetrize_properties.
    !>
    !> \param[inout] this    Properties object to update.
    !> \param[in]    lutarg  File unit to read from.
    !> \param[in]    pset    Selected data set in the file.
    !>
    subroutine read_properties (this, lutarg, pset)

        use mpi_gbl, only: mpi_xermsg

        class(molecular_properties_data), intent(inout) :: this
        integer,                          intent(in)    :: pset, lutarg

        integer  :: i, j, keyh, key, iset, nrecs, nnuc, ntarg, nmom, isw, iseq, icharg, inx(8), ierr, symm, spin, spinsym
        real(wp) :: rmass, x, y, z, dnx, inertia(3)

        character(len=3)   :: catom
        character(len=26)  :: head
        character(len=256) :: line

        call this % clean

        rewind lutarg

        iset = 0

        do while (iset /= pset)

            ! read the first number (file type)
            read (lutarg, *, end = 10) keyh
            backspace lutarg

            if (keyh /= 6) then
                call mpi_xermsg('class_molecular_properties_data', 'read_properties', &
                                'Invalid property file key (expected "6").', 1, 1)
            end if

            ! read the first line (file format)
            read (lutarg, '(A)') line
            backspace lutarg

            ! read the set header
            select case (len_trim(line))
                case (80)  ! old DENPROP format
                    read (lutarg, '(I1,6I3,1X,3D20.12)') keyh, iset, nrecs, nnuc, ntarg, nmom, isw, inertia
                case (87)  ! current DENPROP format
                    read (lutarg, '(I1,I3,I6,I3,I4,I6,I3,1X,3D20.12)') keyh, iset, nrecs, nnuc, ntarg, nmom, isw, inertia
                case (104)  ! CDENPROP format
                    read (lutarg, '(I1,I3,I9,I9,I9,I9,I3,1X,3D20.12)') keyh, iset, nrecs, nnuc, ntarg, nmom, isw, inertia
                case default
                    call mpi_xermsg('class_molecular_properties_data', 'read_properties', &
                                    'Failed to detect version of the property file (80/87/104)', 1, 1)
            end select

            if (iset == pset) then
                this % no_nuclei = nnuc
                this % no_states = ntarg
                this % non_zero_properties = nmom
                allocate (this % energies_index(2, ntarg), this % energies(ntarg), &
                          this % properties_index(4, nmom), this % properties(nmom), stat = ierr)
                if (ierr /= 0) then
                    call mpi_xermsg('class_molecular_properties_data', 'read_properties', 'Memory initialization failed.', 1, 1)
                end if
            end if

            this % nuclear_prop_flag = isw

            ! read nuclei information
            do i = 1, min(nnuc, maxnuc)
                read (lutarg, *, end = 11) key, iseq, catom, icharg, rmass, x, y, z
                if (key /= 8) then
                    call mpi_xermsg('class_molecular_properties_data', 'read_properties', &
                                    'Invalid nuclear data key (expected 8).', 1, 1)
                end if
                if (iset == pset) then
                    ! store nucleus data
                    this % cname(i)  = catom
                    this % charge(i) = icharg
                    this % mass(i) = rmass
                    this % xnuc(i) = x
                    this % ynuc(i) = y
                    this % znuc(i) = z
                end if
            end do

            ! read target information
            do i = 1, ntarg
                read (lutarg, *, end = 11) (inx(j),j=1,8), dnx, head
                if (inx(1) /= 5) then
                    call mpi_xermsg('class_molecular_properties_data', 'read_properties', &
                                    'Invalid target data key (expected 5). Or too many nuclei (maxnuc)?', 1, 1)
                end if
                if (iset == pset) then
                    ! find this spin-symmetry in symmetries index (or end with spinsym = this % no_symmetries + 1 if not found)
                    do spinsym = 1, this % no_symmetries
                        if (this % symmetries_index(1, spinsym) == inx(5) .and. &
                            this % symmetries_index(2, spinsym) == inx(6)) exit
                    end do

                    ! update symmetry count
                    this % no_symmetries = max(spinsym, this % no_symmetries)
                    if (this % no_symmetries > tsym_max) then
                        call mpi_xermsg('class_molecular_properties_data', 'read_properties', &
                                        'Too many symmetries in property file for tsym_max limit.', tsym_max, 1)
                    end if

                    ! update data for this symmetry
                    this % symmetries_index(1, spinsym) = inx(5)
                    this % symmetries_index(2, spinsym) = inx(6)
                    this % symmetries_index(3, spinsym) = this % symmetries_index(3, spinsym) + 1

                    ! store target data
                    this % energies_index(1, i) = 1 + count(this % energies_index(2, 1:i-1) == spinsym)  ! target index
                    this % energies_index(2, i) = spinsym ! spin-symmetry index
                    this % energies(i) = dnx
                end if
            end do

            ! read property information
            do i = 1, nmom
                read (lutarg, *, end = 11) (inx(j),j=1,8), dnx
                if (inx(1) /= 1) then
                    call mpi_xermsg('class_molecular_properties_data', 'read_properties', &
                                    'Invalid property data key (expected 1).', 1, 1)
                end if
                if (iset == pset) then
                    ! store property data
                    this % properties_index(1, i) = inx(2)  ! index of the first state
                    this % properties_index(2, i) = inx(4)  ! index of the second state
                    this % properties_index(3, i) = inx(7)  ! angular momentum of the property
                    this % properties_index(4, i) = inx(8)  ! angular momentum projection multiplied by projectile charge
                    this % properties(i) = dnx
                end if
            end do

        end do

    10 return
    11 call mpi_xermsg('class_molecular_properties_data', 'read_properties', 'Unexpected end of file.', 1, 1)

    end subroutine read_properties


    !> \brief   Symmetrize the sparse property matrices
    !> \authors J Benda
    !> \date    2019
    !>
    !> Dump selected sparse property matrices to dense matrices, correctly duplicating off-diagonal
    !> elements to both positions.
    !>
    !> \param[in]    this  Properties object to read.
    !> \param[in]    lmin  Lower limit on included angular momentum transfer.
    !> \param[in]    lmax  Upper limit on included angular momentum transfer.
    !> \param[inout] prop  Output property matrix (will be re/allocated if needed).
    !>
    subroutine decompress (this, lmin, lmax, prop)

        use mpi_gbl, only: mpi_xermsg

        class(molecular_properties_data), intent(in)    :: this
        real(wp), allocatable,            intent(inout) :: prop(:,:,:)
        integer,                          intent(in)    :: lmin, lmax

        integer :: istate, jstate, i, l, m, ip, ierr, nstates, npropts

        nstates = this % no_states
        npropts = max(0, (lmax + 1)**2 - lmin**2)

        ! deallocate the output matrix if too small
        if (allocated(prop)) then
            if (size(prop, 1) /= nstates .or. &
                size(prop, 2) /= nstates .or. &
                size(prop, 3) <  npropts) then
                deallocate (prop)
            end if
        end if

        ! allocate the output matrix
        if (.not. allocated(prop)) then
            allocate (prop(nstates, nstates, npropts), stat = ierr)
            if (ierr /= 0) then
                call mpi_xermsg('class_molecular_properties_data', 'decompress', 'Failed to allocate properties matrix.', 1, 1)
            end if
        end if

        ! reset output matrix
        prop(:,:,:) = 0

        ! populate the matrix (duplicate off-diagonal elements)
        do i = 1, this % non_zero_properties
            istate = this % properties_index(1, i)
            jstate = this % properties_index(2, i)
            l      = this % properties_index(3, i)
            m      = this % properties_index(4, i)
            if (lmin <= l .and. l <= lmax) then
                ip = l*l + l + m
                prop(istate, jstate, ip - lmin**2 + 1) = this % properties(i)
                prop(jstate, istate, ip - lmin**2 + 1) = this % properties(i)
            end if
        end do

    end subroutine decompress


    !> \brief   Write the DENPROP-like property file
    !> \authors A Harvey, Z Masin, J Benda
    !> \date    2011 - 2019
    !>
    !> This subroutine will write the DENPROP-like "fort.24" property file. It automatically adjusts the output format
    !> when there are too many states for fit into old restrictions. If possible, it uses the most restrictive format to
    !> achieve maximal backward compatibility.
    !>
    !> In case that the property storage is distributed, the master process writes the header and then all processes, one
    !> by one, write properties they have access to. Note that currently all communication happens on MPI world communicator.
    !>
    !> \param[in] this          Property storage to dump to disk.
    !> \param[in] lupropw       Output file unit number.
    !> \param[in] ukrmolp_ints  Format flag.
    !>
    subroutine write_properties (this, lupropw, ukrmolp_ints)

        use ukrmol_interface_gbl, only: molecular_orbital_basis
        use const_gbl,            only: stdout
        use mpi_gbl,              only: myrank, master, nprocs, mpiint, mpi_mod_allgather, mpi_mod_barrier

        ! arguments
        class(molecular_properties_data), intent(in) :: this
        integer, intent(in) :: lupropw
        logical, intent(in) :: ukrmolp_ints

        ! local
        character(len=36) :: str_description
        character(len=5)  :: str_istate
        character(len=:), allocatable :: f_header, f_states, f_propty, irrname
        real(wp) :: p
        integer, allocatable :: to_absolute(:,:)
        integer(mpiint) :: nnz_dense, nnz_sparse, nprops(nprocs)
        logical, parameter :: triangle = .false.
        integer  :: istate, i, j, ielement, iprop, ierr, irank, non_zero_properties, iblock, state_i, state_j, &
                    mgvn_i, mgvn_j, r, c, loc_r, loc_c, l, mq, spin_i, spin_j, nstat_i, nstat_j, spinsym_i, spinsym_j

        ! get total number of non-zero properties in sparse storage
        nnz_sparse = this % non_zero_properties
        call mpi_mod_allgather(nnz_sparse, nprops)
        nnz_sparse = sum(nprops)

        ! get total number of non-zero, non-redundant properties in dense storage
        nnz_dense = 0
        do iblock = 1, this % no_blocks
            spinsym_i = this % dense_index(1, iblock)
            spinsym_j = this % dense_index(2, iblock)
            mgvn_i = this % symmetries_index(1, spinsym_i)
            spin_i = this % symmetries_index(2, spinsym_i)
            mgvn_j = this % symmetries_index(1, spinsym_j)
            spin_j = this % symmetries_index(2, spinsym_j)
            nstat_i = count(this % energies_index(2,:) == spinsym_i)
            nstat_j = count(this % energies_index(2,:) == spinsym_j)
            associate (D => this % dense_blocks(iblock))
                do loc_c = 1, size(D % CV, 2)
                    do loc_r = 1, size(D % CV, 1)
                        if (D % CV(loc_r, loc_c) == 0) cycle
                        call D % local_to_global(loc_r, loc_c, r, c)
                        if (r < 1 .or. nstat_i < r) cycle
                        if (c < 1 .or. nstat_j < c) cycle
                        if (triangle .and. mgvn_i == mgvn_j .and. spin_i == spin_j .and. r < c) cycle
                        nnz_dense = nnz_dense + 1
                    end do
                end do
            end associate
        end do
        call mpi_mod_allgather(nnz_dense, nprops)
        nnz_dense = sum(nprops)

        ! pick one of them
        non_zero_properties = merge(nnz_sparse, nnz_dense, nnz_sparse > 0)

        ! determine write formats
        if (this % swintf_format .and. this % no_states < small_swintf_no_states) then
            ! small, swinterf-compatible
            f_header = '(I1,I3,I6,I3,I4,I6,I3,1x,3D20.12)'
            f_states = '(I1,7I3,D20.12,2X,A36)'
            f_propty = '(I1,7I3,D20.12,2X,A36)'
        else if (this % swintf_format .and. &
                 non_zero_properties + this % no_states + this % no_nuclei < large_swintf_prop_states_nuclei .and. &
                 this % no_nuclei < large_swintf_no_nuclei .and. this % no_states < large_swintf_no_states .and. &
                  non_zero_properties < large_swintf_non_zero_properties) then
            ! large, but still swinterf-compatible
            f_header = '(I1,I3,I6,I3,I4,I6,I3,1x,3D20.12)'
            f_states = '(I1,I8,6I3,D20.12,2X,A36)'
            f_propty = '(I1,I8,I3,I8,4I3,D20.12,2X,A35)'
        else
            ! extra large, only used with photo_outerio and rmt_interface
            f_header = '(I1,I3,I9,I9,I9,I9,I3,1x,3D20.12)'
            f_states = '(I1,I8,I3,I8,4I3,D20.12,2X,A36)'
            f_propty = '(I1,I8,I3,I8,4I3,D20.12)'
        end if

        !ZM correct labels are guaranteed only with the new integral code.
        if (.not. ukrmolp_ints) then
            write (*, '("WARNING: State labels on fort.24 are correct only for C2v and D2h point groups &
                        &and in their canonical orientations.")')
        end if

        ! master writes header information to file
        if (myrank == master) then
            if (lupropw /= stdout) open (lupropw, form = 'formatted')
            write (lupropw, f_header) 6, 1, non_zero_properties + this % no_states + this % no_nuclei, &
                this % no_nuclei, this % no_states, non_zero_properties, this % nuclear_prop_flag, 0.d0, 0.d0, 0.d0
            do i = 1, this % no_nuclei
                ! determine the element index
                ielement = 0
                do j = 1, 103
                    if (this % cname(i)(1:2) == ASYMB(j)) then
                        ielement = j
                    end if
                end do

                ! write information about the i-th nucleus
                if (this % cname(i)(1:4) == 'Scat' .or. this % cname(i)(1:4) == 'Pseu' .or. ielement == 0) then
                    write (lupropw, '(I1,I3,a3,i3,f10.4,3F20.10)') 8, i, this % cname(i)(1:2), nint(this % charge(i)), 0.0_idp, &
                        this % xnuc(i), this % ynuc(i), this % znuc(i)
                else
                    write (lupropw, '(I1,I3,a3,i3,f10.4,3F20.10,I0)') 8, i, this % cname(i)(1:2), nint(this % charge(i)), &
                        AMASS(ielement), this % xnuc(i), this % ynuc(i), this % znuc(i), ielement
                end if
            end do

            do istate = 1, this % no_states
                mgvn_i = this % symmetries_index(1, this % energies_index(2, istate))
                spin_i = this % symmetries_index(2, this % energies_index(2, istate))
                write (str_istate, '(i5)') this % energies_index(1,istate)
                if (ukrmolp_ints) then
                    irrname = ' ' // molecular_orbital_basis % symmetry_data % get_irr_name(mgvn_i + 1)
                else if (this % no_symmetries <= 4) then
                    irrname = cc2v_symb(mgvn_i + 1)
                else if (this % no_symmetries == 8) then
                    irrname = cd2h_symb(mgvn_i + 1)
                end if
                str_description = 'State No. ' // str_istate // ' ' // spin_symb(spin_i) // irrname
                write (lupropw, f_states) 5, istate, 0, 0, mgvn_i, spin_i, 0, 0, this % energies(istate), str_description
            end do
            flush (lupropw)
            if (lupropw /= stdout) close (lupropw)
        end if

        ! only proceed further if some non-zero properties exist
        if (non_zero_properties == 0) return

        ! construct relative index (within spin-symmetry) to absolute index (in energy sorted order) map
        allocate (to_absolute(maxval(this % energies_index(1, :)), this % no_symmetries))
        to_absolute = 0
        do istate = 1, this % no_states
            state_i    = this % energies_index(1, istate)
            spinsym_i  = this % energies_index(2, istate)
            to_absolute(state_i, spinsym_i) = istate
        end do

        ! when no non-zero properties are found in the sparse storage, write from dense storage
        do irank = 0, nprocs - 1
            if (myrank == irank) then
                if (lupropw /= stdout) open (lupropw, form = 'formatted', action = 'write', position = 'append')
                if (nnz_sparse == 0) then
                    do iblock = 1, this % no_blocks
                        spinsym_i = this % dense_index(1, iblock)
                        spinsym_j = this % dense_index(2, iblock)
                        mgvn_i = this % symmetries_index(1, spinsym_i)
                        mgvn_j = this % symmetries_index(1, spinsym_j)
                        spin_i = this % symmetries_index(2, spinsym_i)
                        spin_j = this % symmetries_index(2, spinsym_j)
                        l      = this % dense_index(3, iblock)
                        mq     = this % dense_index(4, iblock)
                        nstat_i = count(this % energies_index(2,:) == spinsym_i)
                        nstat_j = count(this % energies_index(2,:) == spinsym_j)
                        associate (D => this % dense_blocks(iblock))
                            do loc_c = 1, size(D % CV, 2)
                                do loc_r = 1, size(D % CV, 1)
                                    p = D % CV(loc_r, loc_c)
                                    if (p == 0) cycle
                                    call D % local_to_global(loc_r, loc_c, r, c)
                                    if (r < 1 .or. nstat_i < r) cycle
                                    if (c < 1 .or. nstat_j < c) cycle
                                    if (triangle .and. mgvn_i == mgvn_j .and. spin_i == spin_j .and. r < c) cycle
                                    state_i = to_absolute(r, spinsym_i)
                                    state_j = to_absolute(c, spinsym_j)
                                    if (state_i == 0 .or. state_j == 0) cycle
                                    write (lupropw, f_propty) 1, state_i, mgvn_i, state_j, mgvn_j, this % no_nuclei + 1, l, mq, p
                                end do
                            end do
                        end associate
                    end do
                else
                    do iprop = 1, this % non_zero_properties
                        state_i = this % properties_index(1, iprop)
                        state_j = this % properties_index(2, iprop)
                        mgvn_i  = this % symmetries_index(1, this % energies_index(2, state_i))
                        mgvn_j  = this % symmetries_index(1, this % energies_index(2, state_j))
                        l       = this % properties_index(3, iprop)
                        mq      = this % properties_index(4, iprop)
                        p       = this % properties(iprop)
                        write (lupropw, f_propty) 1, state_i, mgvn_i, state_j, mgvn_j, this % no_nuclei + 1, l, mq, p
                    end do
                end if
                flush (lupropw)
                if (lupropw /= stdout) close (lupropw)
            end if
            call mpi_mod_barrier(ierr)
        end do

    end subroutine write_properties


    !> \brief   Sort states by energy
    !> \authors A Harvey, Z Masin, J Benda
    !> \date    2011 - 2019
    !>
    !> Reorder states in the index by their energy and renumber their indices in the property index accordingly.
    !> Does not affect the dense storage in any way.
    !>
    !> \param[in] this  Property storage to update.
    !>
    subroutine sort_by_energy (this)

        use algorithms, only: indexx_stable, rank

        class(molecular_properties_data), intent(inout) :: this

        integer :: istate, iprop, stat_i, symm_i, spin_i
        integer,   allocatable :: energies_index(:,:), sort_index(:), rank_index(:), map_index(:,:,:)
        real(idp), allocatable :: energies(:)

        allocate(energies(this % no_states), energies_index(2, this % no_states))

        ! first sort energies
        allocate (sort_index(this % no_states))
        sort_index = 0
        call indexx_stable(this % no_states, this % energies, sort_index)

        ! then create rank index for re-indexing states
        allocate(rank_index(this % no_states))
        rank_index = 0
        call rank(this % no_states, sort_index, rank_index)

        ! finally apply the re-indexing
        do istate = 1, this % no_states
            energies(istate) = this % energies(sort_index(istate))
            energies_index(:, istate) = this % energies_index(:, sort_index(istate))
        end do

        do iprop = 1, this % non_zero_properties
            this % properties_index(1,iprop) = rank_index(this % properties_index(1,iprop))
            this % properties_index(2,iprop) = rank_index(this % properties_index(2,iprop))
        end do

        call move_alloc(energies, this % energies)
        call move_alloc(energies_index, this % energies_index)

    end subroutine sort_by_energy

    !> \brief   Calculate molecular polarisability from sparse storage.
    !> \authors D Darby-Lewis
    !> \date    2019
    !>
    !> This subroutine calculates molecular polarisabilities for each electronic state.
    !>
    !> It calculates a XX, a YY and a ZZ component by mapping M, the projection
    !> of angular momentum L as follows, M = -1 to Y, 0 to Z and 1 to X.
    !>
    !> It also calculates a mean spherical polarisability as the mean of the
    !> three above components, mean = (XX + YY + ZZ) / 3.
    !>
    !> This subroutine uses but does not affect the sparse storage.
    !>
    !> \param[in]    this        Property object to containing data used.
    !> \param[in]    iwrite      Unit number for text output.
    !> \param[in]    qmoln       Switch for printing molecule.polarisability file. (mgvn, spin).
    !>
    subroutine calculate_polarisability(this, iwrite, qmoln)

        use mpi_gbl, only: mpi_xermsg

        class(molecular_properties_data), intent(in) :: this
        integer,                          intent(in) :: iwrite
        logical,                          intent(in) :: qmoln

        real(wp), allocatable    :: polar(:,:), prop(:,:,:), diff(:,:), gs_sos(:,:)
        real(wp)                 :: p
        integer                  :: i, j, k, l, m, spin_i, spin_j, luqm

        allocate ( prop(this%no_states,this%no_states,-1:1), diff(this%no_states,2), &
                   gs_sos(this%no_states,-1:1), polar(this%no_states,-1:1) )
        polar = 0.0_wp
        prop = 0.0_wp
        diff = 0.0_wp
        gs_sos = 0.0_wp

        write(iwrite, '("Calculating polarisabilities with M = -1 (Y), 0 (Z), 1 (X)")')

!        !$OMP parallel do default(none) shared(this,prop) private(k,i,j,l,m,spin_i,spin_j) schedule(guided)
        do k = 1, this%non_zero_properties
            i = this%properties_index(1,k)
            j = this%properties_index(2,k)
            l = this%properties_index(3,k)
            if ( l /= 1 .or. i == j ) cycle
            spin_i = this%symmetries_index(2,this%energies_index(2,i))
            spin_j = this%symmetries_index(2,this%energies_index(2,j))
            if ( spin_i /= spin_j ) cycle
            m = this%properties_index(4,k)
            p = this%properties(k)
            if ( abs(this%energies(j)-this%energies(i)) < 10E-10_wp ) p = 0.0_wp 
            prop(i,j,m) = p
            prop(j,i,m) = p
        end do
        prop = prop ** 2.0

        !$OMP parallel do default(none) shared(gs_sos,prop,this) private(j,diff) reduction(+:polar) schedule(guided)
        do j = 1, this%no_states
            if ( any( prop(j,:,:) /= prop(:,j,:) ) ) &
                 call mpi_xermsg('class_molecular_properties_data', &
                 'calculate_polarisability', 'Dipoles are not symmetrical.', 1, 1)

            diff(:,1) = this%energies(j)-this%energies(:)
            diff(1:j,2) = this%energies(1:j)-this%energies(1)
            where ( abs(diff) < 10E-10_wp ) diff = 1.0_wp
            do m = -1,1,1
                polar(:, m) = polar(:,m) + prop(:,j,m) / diff(:,1)
                gs_sos(j,m) = sum( prop(1:j,1,m) / diff(1:j,2) )
            end do !m
        end do !j
        polar = polar * 2.0_wp
        gs_sos = gs_sos * 2.0_wp
    
        do i = 1, this%no_states
            write(iwrite, '("State ", i5," alpha_XX  ", D20.12)') i,     polar(i, 1)
            write(iwrite, '("State ", i5," alpha_YY  ", D20.12)') i,     polar(i,-1)
            write(iwrite, '("State ", i5," alpha_ZZ  ", D20.12)') i,     polar(i, 0)
            write(iwrite, '("State ", i5," alpha_MS  ", D20.12)') i, sum(polar(i, :))/3.0
        end do ! i

        if ( qmoln ) then
            open (newunit = luqm, file = 'molecule.polarisability', status = 'unknown')
            write(luqm, *)'********************************************************************'
            write(luqm, *)'                Target GS Polarisability Calculation                '
            write(luqm, *)'********************************************************************'
            write(luqm, *)'                                                                    '
            write(luqm, *)'Polarisability computed in a.u. using 1st order perturbation theory.'
            write(luqm, *)'                                                                    '
            write(luqm, '("Ground State alpha_XX = ", E25.15)') polar(1, 1)
            write(luqm, '("Ground State alpha_YY = ", E25.15)') polar(1,-1)
            write(luqm, '("Ground State alpha_ZZ = ", E25.15)') polar(1, 0)
            write(luqm, *)'                                                                    '
            write(luqm, *)'--------------------------------------------------------------------'
            write(luqm, *)'                                                                    '
            write(luqm, '("Mean Spherical Polarisability alpha_MS = ", E25.15)') sum(polar(1, :))/3.0
            close(luqm)
        end if ! qmoln

    end subroutine calculate_polarisability

end module class_molecular_properties_data
