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

!> \brief   SWEDEN integral module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module SWEDEN_module

    use const_gbl,              only: stdout
    use consts_mpi_ci,          only: NIDX
    use global_utils,           only: INDFUNC, mprod
    use mpi_memory_gbl,         only: mpi_memory_allocate_real, mpi_memory_deallocate_real, mpi_memory_synchronize, local_master
    use precisn,                only: longint, wp
    use BaseIntegral_module,    only: BaseIntegral
    use MemoryManager_module,   only: master_memory
    use Options_module,         only: Options
    use Parallelization_module, only: grid => process_grid
    use Timing_module,          only: master_timer
    use Utility_module,         only: pack_ints, unpack_ints

    implicit none

    public SWEDENIntegral

    private

    type, extends(BaseIntegral) :: SWEDENIntegral

        !Our integrals
#ifdef mpithree
        real(wp), pointer      :: one_electron_integral(:)
        real(wp), pointer      :: two_electron_integral(:)
#else
        real(wp), allocatable  :: one_electron_integral(:)
        real(wp), allocatable  :: two_electron_integral(:)
#endif

        integer                :: one_electron_window
        integer                :: two_electron_window

        integer                :: num_one_electron_integrals
        integer                :: num_two_electron_integrals

        real(wp), allocatable  :: xnuc(:),ynuc(:),znuc(:),charge(:)
        character(len=8), allocatable :: CNAME(:)

        integer                :: num_unique_pairs          !> How many IJKL pairs we have
        integer(longint), allocatable :: pair_labels(:,:)   !> The list of unique labels
        integer, allocatable   :: num_orbitals_sym(:)       !> The number of labels per symmetry

        integer                :: max_number_pair_sets
        integer                :: num_two_electron_blocks
        integer, allocatable   :: one_electron_pointer(:)
        integer, allocatable   :: two_electron_pointer(:)

        integer                :: num_PQ, num_RS
        integer                :: num_pair_idx

        !>the pair id
        integer, allocatable   :: pair_idx(:)
        integer, allocatable   :: orbital_idx(:)
        integer, allocatable   :: symmetry_idx(:)

    contains

        procedure, public  :: initialize_self => initialize_SWEDEN
        procedure, public  :: finalize_self => finalize_SWEDEN
        procedure, public  :: load_integrals => load_integrals_SWEDEN
        procedure, public  :: get_integral_ijklm => get_integral_SWEDEN
        procedure, public  :: write_geometries => write_geometries_SWEDEN
        procedure, public  :: destroy_integrals => destroy_integral_SWEDEN
        procedure, private :: count_num_pairs
        procedure, private :: generate_pairs
        procedure, private :: generate_pointer_table
        procedure, private :: generate_pair_index
        procedure, private :: generate_orbital_index
       !procedure, private :: read_one_electron_integrals
       !procedure, private :: read_two_electron_integrals
        procedure, private :: get_one_electron_index
        procedure, private :: get_two_electron_index

    end type SWEDENIntegral

contains

    integer function count_num_pairs (this)
        use global_utils, only: mprod

        class(SWEDENIntegral)  :: this
        integer                :: i, j, k, l
        integer                :: jend, lend
        integer                :: ij, kl, ijkltest

        count_num_pairs = 0

        do i = 1, this % num_symmetries
            count_num_pairs = count_num_pairs + 1
        end do

        if (this % num_symmetries == 1) return

        !Generate IJIJ pairs
        do i = 2, this % num_symmetries
            jend = i - 1
            do j = 1, jend
                lend = j - 1
                count_num_pairs = count_num_pairs + 1
            end do
        end do

        !Generate IIJJ pairs
        do i = 2, this % num_symmetries
            jend = i - 1
            do j = 1, jend
                lend = j - 1
                count_num_pairs = count_num_pairs + 1
            end do
        end do

        if (this % num_symmetries < 4) return

        do i = 2, this % num_symmetries
            jend = i - 1
            do j = 1, jend
                do k = 2, i
                    lend = k - 1
                    if (i == k) lend = j
                    do l = 1, lend
                        if (i == k .and. j == l) cycle
                        ij = mprod(i, j, 0, stdout)
                        kl = mprod(k, l, 0, stdout)
                        ijkltest = mprod(ij, kl, 0, stdout)
                        if (ijkltest /= 1) cycle
                        count_num_pairs = count_num_pairs + 1
                    end do
                end do
            end do
        end do

    end function count_num_pairs


    subroutine generate_pairs (this)
        class(SWEDENIntegral)  :: this
        integer                :: i, j, k, l
        integer                :: jend, lend
        integer                :: ij, kl, ijkltest
        integer                :: pair_number

        pair_number = 0

        write (stdout, *) 'All IIII blocks of integrals'

        do i = 1, this % num_symmetries
            j = i - 1
            pair_number = pair_number + 1
            call pack_ints(j, j, j, j, 0, 0, 0, 0, this % pair_labels(:, pair_number))
            write (stdout, *) pair_number, j, j, j, j
        end do

        if (this % num_symmetries == 1) return

        !Generate IJIJ pairs
        write (stdout, *) 'All IJIJ blocks'
        do i = 2, this % num_symmetries
            jend = i - 1
            do j = 1, jend
                lend = j - 1
                pair_number = pair_number + 1
                call pack_ints(jend, lend, jend, lend, 0, 0, 0, 0, this % pair_labels(:, pair_number))
                write (stdout, *) pair_number, jend, lend, jend, lend
            end do
        end do

        !Generate IIJJ pairs
        write (stdout, *) 'All IIJJ blocks'
        do i = 2, this % num_symmetries
            jend = i - 1
            do j = 1, jend
                lend = j - 1
                pair_number = pair_number + 1
                call pack_ints(jend, jend, lend, lend, 0, 0, 0, 0, this % pair_labels(:, pair_number))
                write (stdout, *) pair_number, jend, jend, lend, lend
            end do
        end do

        if (this % num_symmetries < 4) return
        write (stdout, *)'All IJKL blocks'
        do i = 2, this % num_symmetries
            jend = i - 1
            do j = 1, jend
                do k = 2, i
                    lend = k - 1
                    if (i == k) lend = j
                    do l = 1, lend
                        if (i == k .and. j == l) cycle
                        ij = mprod(i, j, 0, stdout)
                        kl = mprod(k, l, 0, stdout)
                        ijkltest = mprod(ij, kl, 0, stdout)
                        if (ijkltest /=1 ) cycle
                        pair_number = pair_number + 1
                        call pack_ints(i - 1, j - 1, k - 1, l - 1, 0, 0, 0, 0, this % pair_labels(:, pair_number))
                        write (stdout, *) pair_number, i - 1, j - 1, k - 1, l - 1
                    end do
                end do
            end do
        end do

    end subroutine generate_pairs


    subroutine generate_pointer_table (this)
        class(SWEDENIntegral) :: this

        integer :: label(8), num_orbs(4)
        integer :: ido, jdo, nam1, nam2, num_pq, num_rs, num_pr
        integer :: block_number, err

        this % num_two_electron_blocks = INDFUNC(this % num_symmetries + 1, 0)
        this % num_two_electron_blocks = INDFUNC(this % num_two_electron_blocks + 1, 0)
        this % num_two_electron_integrals = 1

        allocate(this % one_electron_pointer(this % num_symmetries + 1), stat = err)

        call master_memory % track_memory(storage_size(this % one_electron_pointer)/8, &
                                          size(this % one_electron_pointer), err, 'SWEDEN::One electron pointer')

        this % one_electron_pointer(1) = 1

        do ido = 1, this % num_symmetries
            this % one_electron_pointer(ido + 1) = this % one_electron_pointer(ido) + INDFUNC(this % num_orbitals_sym(ido) + 1, 0)
        end do

        this % num_one_electron_integrals = this % one_electron_pointer(this % num_symmetries + 1) - 1

        allocate(this % two_electron_pointer(this % num_two_electron_blocks), stat = err)

        call master_memory % track_memory(storage_size(this % two_electron_pointer)/8, &
                                          size(this % two_electron_pointer), err, 'SWEDEN::Two electron pointer')

        this % two_electron_pointer(:) = 0

        write (stdout, "(//,5X,'Integral Storage Table follows : ',/)")

        do ido = 1, this % num_unique_pairs
            !> Get our lsabels
            call unpack_ints(this % pair_labels(:, ido), label)

            do jdo = 1, 4
                label(jdo) = label(jdo) + 1
                num_orbs(jdo) = this % num_orbitals_sym(label(jdo))
            end do

            ! if IIJJ
            if (label(1) == label(2)) then
                num_PQ = INDFUNC(num_orbs(1) + 1, 0)
            else
                num_PQ = num_orbs(1) * num_orbs(2)
            end if
            NAM1 = INDFUNC(max(label(1), label(2)) + 1, 0) - abs(label(1) - label(2))

            ! if JJ
            if (label(3) == label(4)) then
                num_RS = INDFUNC(num_orbs(3) + 1, 0)
            else
                num_RS = num_orbs(3) * num_orbs(4)
            end if

            NAM2 = INDFUNC(max(label(3), label(4)) + 1, 0) - abs(label(3)-label(4))
            block_number = INDFUNC(max(NAM1, NAM2) + 1, 0) - abs(NAM1 - NAM2)
            this % two_electron_pointer(block_number) = this % num_two_electron_integrals

            write (stdout, "(3X,4I3,1X,I5,1X,I10)") label(1), label(2), label(3), label(4), &
                                                    block_number, this % num_two_electron_integrals

            if (NAM1 == NAM2) then
                num_PR = INDFUNC(num_PQ + 1, 0)
            else
                num_PR = num_PQ * num_RS
            end if

            this % max_number_pair_sets = max(this % max_number_pair_sets, num_PQ, num_RS)
            if (num_PR > 0) this % num_two_electron_integrals = this % num_two_electron_integrals + num_PR + 1
        end do

        !Move to the end of the block
        this % num_two_electron_integrals = this % num_two_electron_integrals - 1

        write (stdout, "('No of 1 electron integrals = ',i8)") this % num_one_electron_integrals
        write (stdout, "('No of 2 electron integrals = ',i8)") this % num_two_electron_integrals

    end subroutine generate_pointer_table


    subroutine generate_orbital_index (this)
        class(SWEDENIntegral) :: this
        integer :: total_num_orbitals, orbital, iorb, isym, max_orbitals, ifail

        total_num_orbitals = sum(this % num_orbitals_sym(:))

        allocate(this % orbital_idx(total_num_orbitals), this % symmetry_idx(total_num_orbitals), stat = ifail)

        call master_memory % track_memory(storage_size(this%orbital_idx)/8,  size(this % orbital_idx),  ifail, 'SWEDEN::pair_idx')
        call master_memory % track_memory(storage_size(this%symmetry_idx)/8, size(this % symmetry_idx), ifail, &
                                          'SWEDEN::symmetry_idx')

        this % orbital_idx(:) = 0
        this % symmetry_idx(:) = 0

        max_orbitals = 0
        orbital = 0

        do isym = 1, this % num_symmetries
            max_orbitals = max(max_orbitals, this % num_orbitals_sym(isym)**2)
            do iorb = 1, this % num_orbitals_sym(isym)
                orbital = orbital + 1
                this % orbital_idx(orbital) = iorb
                this % symmetry_idx(orbital) = isym - 1
            end do
        end do

        this % num_pair_idx = max((this % num_symmetries**2 + this % num_symmetries) / 2 + 1, &
                                  this % max_number_pair_sets, &
                                  300, &
                                  max_orbitals)

    end subroutine generate_orbital_index


    subroutine generate_pair_index (this)
        class(SWEDENIntegral) :: this
        integer :: ido, idx, ifail

        allocate(this % pair_idx(0:this%num_pair_idx), stat = ifail)
        call master_memory % track_memory(storage_size(this % pair_idx)/8, size(this % pair_idx), ifail, 'SWEDEN::pair_idx')

        idx = 0
        this % pair_idx(0) = 0
        do ido = 1, this % num_pair_idx
            this % pair_idx(ido) = idx
            idx = idx + ido
        end do

    end subroutine generate_pair_index


    subroutine initialize_SWEDEN (this, option)
        class(SWEDENIntegral)      :: this
        class(Options), intent(in) :: option
        integer :: ido, jdo, ifail

        write (stdout, "('Using SWEDEN integral')")

        this % num_unique_pairs = 0
        this % num_symmetries = option % num_syms

        allocate(this % num_orbitals_sym(this % num_symmetries), stat = ifail)
        call master_memory % track_memory(storage_size(this % num_orbitals_sym)/8, &
                                          size(this % num_orbitals_sym), ifail, 'SWEDEN::num_orbitals_sym')

        this % num_orbitals_sym(:) = option % num_electronic_orbitals_sym(:)
        this % num_unique_pairs = this % count_num_pairs()

        allocate(this % pair_labels(NIDX, this % num_unique_pairs), stat = ifail)
        call master_memory % track_memory(storage_size(this % pair_labels)/8, &
                                          size(this % pair_labels), ifail, 'SWEDEN::Pair labels')

        write (stdout, *) 'SWEDEN Integral -- Generating pairs'

        !>Generate indexes
        call this % generate_pairs
        call this % generate_orbital_index
        call this % generate_pair_index

        write (stdout, "('1',/,10x,'D2h Two Electron Integral Box ','Information')")
        write (stdout, "(    /,10x,'No. of Boxes = ',i4,/)") this % num_unique_pairs
        write (stdout, "(      10x,'Box Descripti ons: ',/)")
        write (stdout, *) 'SWEDEN Integral -- Generating block pointers'

        call this % generate_pointer_table

    end subroutine initialize_SWEDEN


    subroutine finalize_SWEDEN (this)
        class(SWEDENIntegral) :: this
    end subroutine finalize_SWEDEN


    !>This is just a copy from scatci_routines with only the relevant SWEDEN parts
    subroutine load_integrals_SWEDEN (this, iounit)
        use params,          only: ctrans1, cpoly
        use global_utils,    only: INTAPE
        use scatci_routines, only: SEARCH, CHN2E

        class(SWEDENIntegral) :: this
        integer, intent(in)   ::  iounit

        !     Sweden header arrays
        character(len=4), dimension(4)  :: LABEL
        character(len=4), dimension(33) :: NAM1
        integer(longint), dimension(8)  :: NAO, NCORE
        integer(longint), dimension(20) :: NHE, NOB
        real(wp), allocatable           :: xtempe(:), xtempp(:)

        integer(longint) :: I, NSYM, IONEIN
        integer          :: ifail, k, l, ind, ii
        real(wp)         :: EN

        call master_timer % start_timer('SWEDEN load')

        write (stdout, *) ' Loading SWEDEN Integrals into core'

        !allocate the integral space, this will be replaced by a caching system
        !One electron integral
        !allocate(this%one_electron_integral(2*this%num_one_electron_integrals),stat=ifail)
        this % one_electron_window = mpi_memory_allocate_real(this % one_electron_integral, 2 * this % num_one_electron_integrals)

        !if(ifail /=0) stop "could not allocate 1 electron integral"

        !Two electron integral
        this % two_electron_window = mpi_memory_allocate_real(this % two_electron_integral, this % num_two_electron_integrals)

        call mpi_memory_synchronize(this % one_electron_window)
        call mpi_memory_synchronize(this % two_electron_window)

        !if( (this%two_electron_window == -1 .and. this%one_electron_window == -1)    .or. local_rank == local_master) then
            !if(ifail /=0) stop "could not allocate 2 electron integral"

        read (iounit) (LABEL(I), I = 1, 4)

        write (stdout, fmt='(/,5X,''Label on Sweden tape = '',4A)') (LABEL(I), I = 1, 4)

        read (iounit) NSYM, EN, (NAO(I), NOB(I), NCORE(I), I = 1, NSYM)

        write (stdout, "(' NSYM =',I5,3X,'NOB  =',20I5)") NSYM, (NOB(I), I = 1, NSYM)
        write (stdout, "(/,10x,'Nuclear repulsion energy = ',f15.7)") EN
        write (stdout, "(/10x,'  NAO  NMO  NCORE',/,(10x,2I5,i7))") (NAO(I), NOB(I), NCORE(I), I = 1, NSYM)

        this % core_energy = EN

        do I = 1, NSYM
            if (NCORE(I) /= 0) then
                write (stdout, "(10X,'Non-zero CORE in Sweden. Sym no. = ',i5,' Ncore = ',I5,/)") I, NCORE(I)
                write (stdout, "(10X,'Core Energy = ',F15.7,' Hartrees ',/)") EN
            end if
        end do

        !Do Erro checkin here will ignore for now

        IONEIN = 0
        do I = 1, NSYM
            IONEIN = IONEIN + (NOB(I) * (NOB(I) + 1)) / 2
        end do

        if ((this % two_electron_window == -1 .and. this % one_electron_window == -1) .or. grid % lrank == local_master) then
            call SEARCH(iounit, ctrans1, ifail)
            read (iounit)
            call INTAPE(iounit, this % one_electron_integral, this % num_one_electron_integrals)

            !Positron case
            if (this % positron_flag /= 0) then
                allocate(xtempe(this % num_one_electron_integrals), stat = ifail)
                call master_memory % track_memory(storage_size(xtempe)/8, size(xtempe), ifail, 'SWEDEN::xtempe')

                allocate(xtempp(this % num_one_electron_integrals), stat = ifail)
                call master_memory % track_memory(storage_size(xtempp)/8, size(xtempp), ifail, 'SWEDEN::xtempp')

                xtempe(1:this % num_one_electron_integrals) = this % one_electron_integral(1:this % num_one_electron_integrals)
                call INTAPE(iounit, xtempp, this % num_one_electron_integrals)

                !Are we using Quantemol-N?
                if (this % quantamoln_flag) then
                    ind = 0
                    do k = 1, NSYM
                        DO l = 1, NOB(k)
                            ind = ind + l
                            write (413, "(i6,3x,i6,3x,i6,5x,f12.7)") k - 1, l, l, XTEMPP(ind)
                        end do
                    end do
                end if

                do I = 1, this % num_one_electron_integrals
                    this % one_electron_integral(I + this % num_one_electron_integrals) = 2.0_wp * XTEMPP(I) - XTEMPE(I)
                end do

                call master_memory % free_memory(storage_size(xtempp)/8, size(xtempp))
                deallocate(xtempp)
                call master_memory % free_memory(storage_size(xtempe)/8, size(xtempe))
                deallocate(xtempe)
            end if

            call CHN2E(iounit, stdout, this % two_electron_integral, this % num_two_electron_integrals)
        end if

        1600 format(' Integrals read successfully:',       /, &
                    ' 1-electron integrals, NINT1e =',i10, /, &
                    ' 2-electron integrals, NINT2e =',i10)

        write (stdout, 1600) this % num_one_electron_integrals, this % num_two_electron_integrals

        this % nhe(:) = nhe(:)
        this % nnuc = 0
        this % dtnuc(:) = 0
        this % dtnuc(1) = en
        this % nhe(:) = nob(:)
        this % dtnuc(22:41) = 0
        this % dtnuc(2:21) = 0

        !Reading Geometries
        rewind iounit
        call SEARCH(iounit, cpoly, ifail)
        read (iounit) this % nnuc

        allocate(this % xnuc(this % nnuc), this % ynuc(this % nnuc), this % znuc(this % nnuc), &
                 this % charge(this % nnuc), this % CNAME(this % nnuc))

        do i = 1, this % nnuc
            read (iounit) this % cname(i), ii, this % xnuc(i), this % ynuc(i), this % znuc(i), this % charge(i)
        end do

        call master_timer % stop_timer('SWEDEN load')

        !endif

        call mpi_memory_synchronize(this % one_electron_window)
        call mpi_memory_synchronize(this % two_electron_window)

    end subroutine load_integrals_SWEDEN


    function get_integral_SWEDEN (this, i, j, k, l, m) result(coeff)
        class(SWEDENIntegral) :: this
        integer, intent(in)   :: i, j, k, l, m
        real(wp) :: coeff
        integer  :: ia, ib, integral_idx

        coeff = 0.0

        !One electron case
        if (i == 0) then
            integral_idx = this % get_one_electron_index(j, l, m)
            coeff = this % one_electron_integral(integral_idx)
        else
            integral_idx = this % get_two_electron_index(i, j, k, l, m)
            coeff = this % two_electron_integral(integral_idx)
        end if

    end function get_integral_SWEDEN


    integer function get_one_electron_index (this, i, j, pos)
        class(SWEDENIntegral), intent(in) :: this
        integer,               intent(in) :: i, j, pos
        integer :: m, p, q, ii, jj

        ii = this % orbital_mapping(i)
        jj = this % orbital_mapping(j)

        p = this % orbital_idx(ii)
        q = this % orbital_idx(jj)

        m = this % symmetry_idx(ii)

        if (p < q) then
            get_one_electron_index = this % pair_idx(q) + p + this % one_electron_pointer(m + 1) - 1
        else
            get_one_electron_index = this % pair_idx(p) + q + this % one_electron_pointer(m + 1) - 1
        end if

        if (pos /= 0) get_one_electron_index = get_one_electron_index + this % num_one_electron_integrals
        !write(stdout,"(10i8)") 0,p,0,q,0,get_one_electron_index,m,this%one_electron_pointer(m+1),this%pair_idx(p)

    end function get_one_electron_index


    integer function get_two_electron_index (this, i, j, k, l, m)
        class(SWEDENIntegral), intent(in) :: this
        integer, intent(in)    :: i, j, k, l, m
        integer                :: symmetry
        integer, dimension(4)  :: orb_num_sym, orb_sym, num_orb_sym
        integer                :: mapped(4)
        integer                :: ido
        integer                :: mpq, mrs, mpr, aaa, block_pointer, NOBRS, NOBPQ, IPQRS

        mapped(1) = this % orbital_mapping(i)
        mapped(2) = this % orbital_mapping(j)
        mapped(3) = this % orbital_mapping(k)
        mapped(4) = this % orbital_mapping(l)

        do ido = 1, 4
            orb_num_sym(ido) = this % orbital_idx(mapped(ido))
            symmetry = this % symmetry_idx(mapped(ido)) + 1
            orb_sym(ido) = symmetry
            num_orb_sym(ido) = this % num_orbitals_sym(symmetry)
        end do

        !write(stdout,"(4i4,2x,4i4,2x,4i4)") orb_num_sym(1:4),orb_sym(1:4),num_orb_sym(1:4)

        mpq = this % pair_idx(orb_sym(1)) + orb_sym(2)
        mrs = this % pair_idx(orb_sym(3)) + orb_sym(4)
        mpr = this % pair_idx(mpq) + mrs

        if (MPQ == MRS .and. &
            ( orb_num_sym(1) <  orb_num_sym(3)                                        .or. &
             (orb_num_sym(1) == orb_num_sym(3) .and. orb_num_sym(2) < orb_num_sym(4)) .or. &
             (orb_num_sym(1) <  orb_num_sym(3) .and. orb_num_sym(2) < orb_num_sym(4)) )) then
            aaa = orb_num_sym(1) ; orb_num_sym(1) = orb_num_sym(3) ; orb_num_sym(3) = aaa
            aaa = orb_sym(1)     ; orb_sym(1)     = orb_sym(3)     ; orb_sym(3)     = aaa
            aaa = orb_num_sym(2) ; orb_num_sym(2) = orb_num_sym(4) ; orb_num_sym(4) = aaa
            aaa = orb_sym(2)     ; orb_sym(2)     = orb_sym(4)     ; orb_sym(4)     = aaa
        end if

        block_pointer = this % two_electron_pointer(mpr) - 1

        if (orb_sym(1) == orb_sym(2) .and. orb_sym(3) == orb_sym(4) .and. orb_sym(1) == orb_sym(3)) then
            IPQRS = this % pair_idx(orb_num_sym(1)) + orb_num_sym(2)
            !IF(ipqrs.GT.nobtest)WRITE(6,910)ipqrs, nobtest
            IPQRS = this % pair_idx(IPQRS - 1) + IPQRS - 1 + this % pair_idx(orb_num_sym(3)) + orb_num_sym(4)
        else if (orb_sym(1) == orb_sym(3) .and. orb_sym(2) == orb_sym(4)) then
            IPQRS = (orb_num_sym(1) - 1) * num_orb_sym(2) + orb_num_sym(2) - 1
            ! IF(ipqrs.GT.nobtest)WRITE(6,910)ipqrs, nobtest
            IPQRS = this % pair_idx(IPQRS) + IPQRS + (orb_num_sym(3) - 1) * num_orb_sym(4) + orb_num_sym(4)
        else if (orb_sym(1) == orb_sym(2) .and. orb_sym(3) == orb_sym(4)) then
            NOBRS = this % pair_idx(num_orb_sym(3)) + num_orb_sym(3)
            IPQRS = (this % pair_idx(orb_num_sym(1)) + orb_num_sym(2) - 1)*NOBRS + this % pair_idx(orb_num_sym(3)) + orb_num_sym(4)
        else
            NOBRS = num_orb_sym(3) * num_orb_sym(4)
            IPQRS = ((orb_num_sym(1) - 1) * num_orb_sym(2) + orb_num_sym(2) - 1) * NOBRS + (orb_num_sym(3) - 1) * num_orb_sym(4) &
                    + orb_num_sym(4)
        end if

        get_two_electron_index = IPQRS + block_pointer + 1

        !write(stdout,"(6i8)") i,j,k,l,m,get_two_electron_index

    end function get_two_electron_index


    subroutine write_geometries_SWEDEN (this, iounit)
        use params,          only: cpoly
        use scatci_routines, only: SEARCH

        class(SWEDENIntegral), intent(in) :: this
        integer,               intent(in) :: iounit
        integer :: ido, ifail, i

        do ido = 1, this % nnuc
            write (iounit) this % cname(ido), this % xnuc(ido), this % ynuc(ido), this % znuc(ido), this % charge(ido)
        end do

        write (stdout, "(/,' Nuclear data     X         Y         Z       Charge',/,(3x,a8,2x,4F10.6))") &
            (this % cname(i), this % xnuc(i), this % ynuc(i), this % znuc(i), this % charge(i), i = 1, this % nnuc)

    end subroutine write_geometries_SWEDEN


    subroutine destroy_integral_SWEDEN (this)
        class(SWEDENIntegral) :: this

        !if(allocated(this%one_electron_integral)) then
            !deallocate(this%one_electron_integral)
        !endif

        call mpi_memory_deallocate_real(this % one_electron_integral, &
                                        size(this % one_electron_integral), this % one_electron_window)
        call mpi_memory_deallocate_real(this % two_electron_integral, &
                                        size(this % two_electron_integral), this % two_electron_window)

        !if(allocated(this%two_electron_integral)) then
            !deallocate(this%two_electron_integral)
        !endif

        if (allocated(this % pair_labels)) then
            call master_memory % free_memory(storage_size(this % pair_labels)/8, size(this % pair_labels))
            deallocate(this % pair_labels)
        end if

        !> The number of labels per symmetry
        if (allocated(this % num_orbitals_sym)) then
            call master_memory % free_memory(storage_size(this % num_orbitals_sym)/8, size(this % num_orbitals_sym))
            deallocate(this % num_orbitals_sym)
        end if

        if (allocated(this % one_electron_pointer)) then
            call master_memory % free_memory(storage_size(this % one_electron_pointer)/8, size(this % one_electron_pointer))
            deallocate(this % one_electron_pointer)
        end if

        if (allocated(this % two_electron_pointer)) then
            call master_memory % free_memory(storage_size(this % two_electron_pointer)/8, size(this % two_electron_pointer))
            deallocate(this % two_electron_pointer)
        end if

        if (allocated(this % pair_idx)) then
            call master_memory % free_memory(storage_size(this % pair_idx)/8, size(this % pair_idx))
            deallocate(this % pair_idx)
        end if

        if (allocated(this%orbital_idx)) then
            call master_memory % free_memory(storage_size(this % orbital_idx)/8, size(this % orbital_idx))
            deallocate(this %orbital_idx)
        end if

        if (allocated(this % symmetry_idx)) then
            call master_memory % free_memory(storage_size(this % symmetry_idx)/8, size(this % symmetry_idx))
            deallocate(this % symmetry_idx)
        end if

    end subroutine destroy_integral_SWEDEN

end module SWEDEN_module
