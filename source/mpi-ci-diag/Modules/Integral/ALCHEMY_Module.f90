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

!> \brief   ALCHEMY integral module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module ALCHEMY_module

    use precisn,              only: longint, wp
    use const_gbl,            only: stdout
    use consts_mpi_ci,        only: NIDX
    use BaseIntegral_module,  only: BaseIntegral
    use MemoryManager_module, only: master_memory
    use Options_module,       only: Options
    use global_utils,         only: INDFUNC
    use scatci_routines,      only: RDINT, RDINTS
    use Utility_module,       only: pack_ints, unpack_ints

    implicit none

    public ALCHEMYIntegral

    private

    type, extends(BaseIntegral) :: ALCHEMYIntegral

        !Our integrals
        real(wp), pointer :: one_electron_integral(:)
        real(wp), pointer :: two_electron_integral(:)
        integer           :: num_one_electron_integrals
        integer           :: num_two_electron_integrals
        integer           :: integral_ordering
        integer           :: use_SCF = 0

        !> How many IJKL pairs we have
        integer              :: num_unique_pairs
        !> The list of unique labels
        integer(longint), allocatable :: pair_labels(:,:)
        !> The number of labels per symmetry
        integer, allocatable :: num_orbitals_sym(:)

        integer              :: max_number_pair_sets
        integer, allocatable :: num_two_electron_blocks(:)
        integer, allocatable :: num_one_electron_blocks(:)
        integer, allocatable :: one_electron_pointer(:)
        integer, allocatable :: two_electron_pointer(:)

        integer :: num_PQ, num_RS
        integer :: num_pair_idx

        !>the pair id
        integer, allocatable :: pair_idx(:)
        integer, allocatable :: orbital_idx(:)
        integer, allocatable :: symmetry_idx(:)

    contains

        procedure, public  :: initialize_self    => initialize_ALCHEMY
        procedure, public  :: finalize_self      => finalize_ALCHEMY
        procedure, public  :: load_integrals     => load_integrals_ALCHEMY
        procedure, public  :: get_integral_ijklm => get_integral_ALCHEMY
        procedure, public  :: destroy_integrals  => destroy_integral_ALCHEMY
        procedure, public  :: write_geometries   => write_geometries_ALCHEMY
        procedure, private :: count_num_pairs
        procedure, private :: generate_pairs
        procedure, private :: generate_pointer_table
        procedure, private :: generate_pair_index
        procedure, private :: generate_orbital_index
       !procedure, private :: read_one_electron_integrals
       !procedure, private :: read_two_electron_integrals
        procedure, private :: get_one_electron_index
        procedure, private :: get_two_electron_index

    end type ALCHEMYIntegral

contains

    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer function count_num_pairs (this)
        class(ALCHEMYIntegral) :: this
        integer :: pair_number, mdo, ndo, ido, jdo, nn(5), mSym, begin, pass, istart, num_sym

        num_sym = this % num_symmetries
        mSym = num_sym + num_sym - 1
        count_num_pairs = 0

        if (this % use_SCF == 0) then
            !Standard Alchemy list
            do mdo = 1, mSym
                begin = mdo / 2 + 1
                nn(1) = mdo - 1
                do ndo = begin, num_sym
                    nn(2) = ndo - 1
                    nn(3) = abs(ndo - mdo)
                    do ido = begin, ndo
                        nn(4) = ido - 1
                        nn(5) = abs(ido - mdo)
                        count_num_pairs = count_num_pairs + 1
                    end do
                end do
            end do
        else
            do pass = 1, 2
                do mdo = pass, mSym
                    begin = mdo / 2 + 1
                    nn(1) = mdo - 1
                    do ndo = mdo / 2 + pass, num_sym
                        nn(2) = ndo - 1
                        nn(3) = abs(ndo - mdo)
                        istart = ndo + 1 - pass
                        if (mdo /= 1 .and. pass == 1) begin = istart
                        do ido = begin, istart
                            nn(4) = ido - 1
                            nn(5) = abs(ido - mdo)
                            count_num_pairs = count_num_pairs + 1
                        end do
                    end do
                end do
                mSym = mSym - 2
            end do
        end if

    end function count_num_pairs


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine generate_pairs (this)
        class(ALCHEMYIntegral) :: this
        integer :: pair_number, npqrs, nn(5), mdo, ndo, ido, jdo, I, IBGN, IF, IPASS, J, K, M, MSYM, N, NBGN, begin, pass, num_sym
        integer(longint) :: packed_label(NIDX)
        integer,allocatable :: lpqrs(:,:)

        num_sym = this % num_symmetries
        mSym = num_sym + num_sym - 1

        allocate(lpqrs(5, this % num_unique_pairs))

        pair_number = 0

        if (this % use_SCF == 0) then
            !Standard Alchemy list
            do mdo = 1, mSym
                begin = mdo / 2 + 1
                nn(1) = mdo - 1
                do ndo = begin, num_sym
                    nn(2) = ndo - 1
                    nn(3) = abs(ndo - mdo)
                    do ido = begin, ndo
                        nn(4) = ido - 1
                        nn(5) = abs(ido - mdo)
                        pair_number = pair_number + 1
                        do jdo = 1, 5
                            lpqrs(jdo,pair_number) = NN(jdo)
                        end do
                    end do
                end do
            end do
        else
            do IPASS = 1, 2
                do M = IPASS, MSYM
                    IBGN = M / 2 + 1
                    NN(1) = M - 1
                    do N = M / 2 + IPASS, num_sym
                        NN(2) = N - 1
                        NN(3) = ABS(N - M)
                        IF = N + 1 - IPASS
                        if (M /= 1 .and. IPASS == 1) IBGN = IF
                        do I = IBGN, IF
                            NN(4) = I - 1
                            NN(5) = ABS(I - M)
                            pair_number = pair_number + 1
                            do J = 1, 5
                                LPQRS(J,pair_number) = NN(J)
                            end do
                        end do
                     end do
                 end do
                 mSym = mSym - 2
            end do
        end if

        if (pair_number /= this % num_unique_pairs) stop "Error in pair count"

        do ido = 1, this % num_unique_pairs
            call pack_ints(LPQRS(2,ido), LPQRS(3,ido), LPQRS(4,ido), LPQRS(5,ido), 0, 0, 0, 0, packed_label)
            if (LPQRS(2,ido) + LPQRS(3,ido) /= LPQRS(1,ido)) packed_label(1) = -packed_label(1)
            this % pair_labels(:,ido) = packed_label
        end do

        write (stdout, "(' NMPRS=',I5//' LPQRS'/(2X,5(2X,I3)))") this % num_unique_pairs, &
            ((LPQRS(ido,jdo), ido = 1, 5), jdo = 1, this % num_unique_pairs)

        deallocate(lpqrs)

    end subroutine generate_pairs


    subroutine generate_pointer_table (this)
        class(ALCHEMYIntegral) :: this
        integer                :: label(8), num_orbs(4), size_two_pointer
        integer(longint)       :: packed_label(NIDX)
        integer                :: ido, jdo, nam1, nam2, num_pq, num_rs, num_pr
        integer                :: block_number, mvl, md, mm, m, mpa, mpb, mp

        !Lets allocate our one elctron pointers

        size_two_pointer = this % num_symmetries * (this % num_symmetries + 1) / 2
        do ido = 2, this % num_symmetries
            size_two_pointer = size_two_pointer + (this % num_symmetries - ido + 2) * (this % num_symmetries - ido + 1)
        end do

        allocate(this % num_two_electron_blocks(2 * this % num_symmetries + 1))

        this % num_two_electron_blocks(1) = 0
        this % num_two_electron_blocks(2) = INDFUNC(this % num_symmetries + 1, 0)

        mvl = 2 * this% num_symmetries - 1

        do ido = 2, mvl
            md = (ido - 2) / 2
            mm = this % num_symmetries - md - 1
            this% num_two_electron_blocks(ido + 1) = this % num_two_electron_blocks(ido) + INDFUNC(mm + 1, 0)
        end do

        allocate(this % one_electron_pointer(this % num_symmetries + 1))

        this % one_electron_pointer(1) = 1

        do ido = 1, this % num_symmetries
            this % one_electron_pointer(ido + 1) = this % one_electron_pointer(ido) + INDFUNC(this % num_orbitals_sym(ido) + 1, 0)
        end do

        allocate(this % two_electron_pointer(this % num_two_electron_blocks(MVL + 1)))

        do ido =1, this % num_two_electron_blocks(MVL + 1)
            this % two_electron_pointer(ido) = 0
        end do

        this % num_one_electron_integrals = this % one_electron_pointer(this % num_symmetries + 1) - 1
        this % num_two_electron_integrals = 1

        write (stdout, "(//,5X,'Integral Storage Table follows : ',/)")

        do ido = 1, this % num_unique_pairs
            !> Get our lsabels
            packed_label(:) = this % pair_labels(:,ido)

            if (packed_label(1) > 0) then
                md = 0
            else
                md = 1
                packed_label(1) = -packed_label(1)
            end if

            call unpack_ints(packed_label, label)

            do jdo = 1, 4
                num_orbs(jdo) = this % num_orbitals_sym(label(jdo) + 1)
            end do

            if (md == 0) then
                m = label(1) + label(2)
            else
                m = abs(label(1) - label(2))
            end if

            if (m == 0) then
                md = -1
            else
                md = (m - 1) / 2
            end if
            mpa = label(1) - md
            mpb = label(3) - md
            if (MPA < MPB) then
                !MP=(MPB*(MPB-1))/2+MPA+this%num_two_electron_blocks(M+1)
                MP = (MPB * (MPB - 1)) / 2 + MPA + this % num_two_electron_blocks(M + 1)
            else
                !MP=(MPA*(MPA-1))/2+MPB+this%num_two_electron_blocks(M+1)
                MP = (MPA * (MPA - 1)) / 2 + MPB + this % num_two_electron_blocks(M + 1)
            end if

            this % two_electron_pointer(mp) = this % num_two_electron_integrals
            !write(stdout,"('ELEC ',3i8)") mp,this%num_two_electron_integrals,this%two_electron_pointer(mp)
            if (label(1) == label(2)) then
                num_PQ = INDFUNC(num_orbs(1) + 1, 0)
            else
                num_PQ = num_orbs(1) * num_orbs(2)
            end if

            if (label(3) == label(4)) then
                num_RS = INDFUNC(num_orbs(3) + 1, 0)
            else
                num_RS = num_orbs(3) * num_orbs(4)
            end if
            this % max_number_pair_sets = max(this % max_number_pair_sets, num_PQ, num_RS)
            if (label(1) == label(3) .and. label(2) == label(4)) then
                this % num_two_electron_integrals = this % num_two_electron_integrals + INDFUNC(num_PQ + 1, 0)
            else
                this % num_two_electron_integrals = this % num_two_electron_integrals + num_PQ * num_RS
            end if
        end do
        !Move to the end of the block
        this % num_two_electron_integrals = this % num_two_electron_integrals - 1

        write (stdout, "('No of 1 electron integrals = ',i8)") this % num_one_electron_integrals
        write (stdout, "('No of 2 electron integrals = ',i8)") this % num_two_electron_integrals

    end subroutine generate_pointer_table


    subroutine generate_orbital_index (this)
        class(ALCHEMYIntegral) :: this
        integer :: total_num_orbitals, orbital, iorb, isym, max_orbitals, ifail

        total_num_orbitals = sum(this % num_orbitals_sym(:))

        allocate(this % orbital_idx(total_num_orbitals), &
                 this % symmetry_idx(total_num_orbitals), stat = ifail)

        call master_memory % track_memory(storage_size(this%orbital_idx)/8,  size(this % orbital_idx),  ifail, 'SWEDEN::pair_idx')
        call master_memory % track_memory(storage_size(this%symmetry_idx)/8, size(this % symmetry_idx), ifail, &
                                          'SWEDEN::symmetry_idx')

        this % orbital_idx(:) = 0
        this % symmetry_idx(:) = 0

        max_orbitals = 0
        orbital = 0

        do isym = 1,this % num_symmetries
            max_orbitals = max(max_orbitals, this % num_orbitals_sym(isym)**2)
            do iorb = 1, this % num_orbitals_sym(isym)
                orbital = orbital + 1
                this % orbital_idx(orbital) = iorb
                this % symmetry_idx(orbital) = isym - 1
            end do
        end do

        this % num_pair_idx = max((this% num_symmetries**2 + this % num_symmetries) / 2 + 1, &
                                  this % max_number_pair_sets, &
                                  300, &
                                  max_orbitals)

    end subroutine generate_orbital_index


    subroutine generate_pair_index (this)
        class(ALCHEMYIntegral) :: this
        integer                :: ido, idx, ifail

        allocate(this % pair_idx(0:this % num_pair_idx), stat = ifail)

        call master_memory % track_memory(storage_size(this % pair_idx)/8, size(this % pair_idx), ifail, 'SWEDEN::pair_idx')

        idx = 0
        this % pair_idx(0) = 0
        do ido = 1, this % num_pair_idx
            this % pair_idx(ido) = idx
            idx = idx + ido
        end do

    end subroutine generate_pair_index


    subroutine initialize_ALCHEMY (this, option)
        class(ALCHEMYIntegral)     :: this
        class(Options), intent(in) :: option
        integer                    :: ido, jdo

        write (stdout, "('Using ALCHEMY integral')")

        this % num_unique_pairs = 0
        this % num_symmetries = option % num_syms
        this % integral_ordering = option % integral_ordering
        this % use_SCF = option % use_SCF

        !Lets begin counting

        allocate(this % num_orbitals_sym(this % num_symmetries))

        this % num_orbitals_sym(:) = option % num_electronic_orbitals_sym(:)
        this % num_unique_pairs = this % count_num_pairs()

        write (stdout, *) 'Number of pairs =',this % num_unique_pairs

        allocate(this % pair_labels(2, this % num_unique_pairs))

        write (stdout,*) 'ALCHEMY Integral -- Generating pairs'

        !>Generate indexes
        call this % generate_pairs
        call this % generate_orbital_index
        call this % generate_pair_index

        write (stdout, "('1',/,10x,'D2h Two Electron Integral Box ','Information')")
        write (stdout, "(    /,10x,'No. of Boxes = ',i4,/)") this % num_unique_pairs
        write (stdout, "(      10x,'Box Descriptions: ',/)")
        write (stdout, *) 'ALCHEMY Integral -- Generating block pointers'

        call this % generate_pointer_table

    end subroutine initialize_ALCHEMY


    subroutine finalize_ALCHEMY (this)
        class(ALCHEMYIntegral) :: this
    end subroutine finalize_ALCHEMY


    !> This is just a copy from scatci_routines with only the relevant ALCHEMY parts.
    !> This will be replaced when the prototype is completed with a 'caching' system.
    subroutine load_integrals_ALCHEMY (this, iounit)
        use params,          only: ctrans1, cpoly
        use global_utils,    only: INTAPE
        use scatci_routines, only: SEARCH, CHN2E

        class(ALCHEMYIntegral)    :: this
        integer, intent(in)       :: iounit
        !     ALCHEMY header arrays
        character(len=4), dimension(4)  :: LABEL
        character(len=4), dimension(33) :: NAM1
        integer, dimension(8)  :: NAO, NCORE
        integer, dimension(20) :: NHE, NOB
        real(wp), allocatable           :: xtempe(:), xtempp(:)

        integer  :: I, NSYM, IONEIN, NT, ND, LTRI, NNUC, N, IA, IAT, NALM
        integer  :: ifail, k, l, ind, IMAX, IMIN, NEND, IMAXP
        real(wp) :: EN, SIGN, x1e
        real(wp) :: vsmall = 1.E-10_wp
        real(kind=wp), dimension(20) :: CHARG, XNUC

        write (stdout, *) ' Loading ALCHEMY Integrals into core'

        !allocate the integral space, this will be replaced by a caching system
        !One electron integral
        allocate(this % one_electron_integral(2 * this % num_one_electron_integrals), stat = ifail)
        if (ifail /= 0) stop "could not allocate 1 electron integral"

        !Two electron integral
        allocate(this % two_electron_integral(this % num_two_electron_integrals), stat = ifail)
        if (ifail /= 0) stop "could not allocate 2 electron integral"

        read (iounit) (NAM1(I), I = 1, 33), NSYM, NT, NNUC, ND, LTRI, (NOB(I), I = 1, NSYM), (ND, I = 1, NT), (ND, I = 1, NT), &
                      (CHARG(I), I = 1, NNUC), (XNUC(I), I = 1, NNUC)

        write (stdout, "(/' Transformed integrals read:',/5X,30A4)") (NAM1(I), I = 1, 30)

        EN = 0.0_wp
        do N = 2, NNUC
            if (abs(charg(n)) < vsmall) cycle
            do I = 1, N - 1
               if (abs(charg(i)) < vsmall) cycle
               EN = EN + CHARG(N) * CHARG(I) / abs(XNUC(N) - XNUC(I))
            end do
        end do

        this % core_energy = EN

        write (stdout, "(' NSYM =',I5,3X,'NOB  =',20I5)") NSYM, (NOB(I), I = 1, NSYM)
        write (stdout, "(/,10x,'Nuclear repulsion energy = ',f15.7)") EN
        write (stdout, "(' NNUC =',I5,3X,'CHARG=',10F10.0/21X,10F10.0)") NNUC, (CHARG(I), I = 1, NNUC)
        write (stdout, "(15X,'XNUC =',10F10.5/21X,10F10.5)") (XNUC(I), I = 1, NNUC)

        IA = 1
        call RDINTS(iounit, NSYM, NOB, LTRI, this % num_one_electron_integrals + 1, IA, this % one_electron_integral, NALM)
        ! print*,this%one_electron_integral(1)
        if (NALM /= 0) stop "UNABLE TO READ ALCHEMY INTEGRALS"

        IA = 1
        call RDINTS(iounit, NSYM, NOB, LTRI, this % num_one_electron_integrals + 1, IA, this % one_electron_integral, NALM)
        ! print*,this%one_electron_integral(1)
        if (NALM /= 0) stop "UNABLE TO READ ALCHEMY INTEGRALS"

        IAT = IA - 1
        ia = 1
        allocate(xtempp(this % num_one_electron_integrals))
        call RDINTS(iounit, NSYM, NOB, LTRI, this % num_one_electron_integrals + 1, IA, xtempp, NALM)
        ! print*,xtempp(1)

        if (NALM == 0) then
            if (IAT /= this % num_one_electron_integrals) stop "INCORRECT NUMBER OF INTEGRALS"
            SIGN = 1.0
            if(abs(this % positron_flag) == 1) SIGN = -1.0
            do I = 1, IAT
                x1e = this % one_electron_integral(I) + xtempp(I)
                this % one_electron_integral(I) = x1e
            end do
            !this%one_electron_integral(1:) = this%one_electron_integral(:) + xtempp(:)
        end if

        !print *,this%one_electron_integral(1)

        IMAX = 0
        IMIN = 0
        NEND = -1
        IMAXP = 0
        IA = 1
        call RDINT(IMIN, IMAX, NEND, imaxp, iounit, IA, this % num_two_electron_integrals + 1, this % two_electron_integral)
        if (NEND /= 1) stop "UNABLE TO READ ALCHEMY INTEGRALS"

        deallocate(xtempp)

        !endif

        1600 FORMAT(' Integrals read successfully:',       /, &
                    ' 1-electron integrals, NINT1e =',i10, /, &
                    ' 2-electron integrals, NINT2e =',i10)
        write (stdout, 1600) this % num_one_electron_integrals, this % num_two_electron_integrals

        this % nhe(:) = nhe(:)
        this % nnuc = NNUC
        this % dtnuc(:) = 0
        this % dtnuc(1) = en
        this % nhe(:) = nob(:)
        this % dtnuc(22:21+NNUC) = XNUC(1:NNUC)
        this % dtnuc(2:1+NNUC) = CHARG(1:NNUC)

    end subroutine load_integrals_ALCHEMY


    function get_integral_ALCHEMY (this, i, j, k, l, m) result(coeff)
        class(ALCHEMYIntegral) :: this
        integer, intent(in)    :: i,j,k,l,m
        real(wp) :: coeff
        integer  :: ia, ib, integral_idx

        coeff = 0.0

        !One electron case
        if (i == 0) then
            integral_idx = this % get_one_electron_index(j, l, m)
            coeff = this % one_electron_integral(integral_idx)
            !print *,coeff
        else
            integral_idx = this % get_two_electron_index(i, j, k, l, m)
            coeff = this % two_electron_integral(integral_idx)
        end if

        !WRITE(stdout,*)i, j, k, l, m
        !write(stdout,*)coeff,integral_idx
        !write(stdout,*) i,j,k,l,m,coeff,integral_idx

    end function get_integral_ALCHEMY


    integer function get_one_electron_index (this, i, j, pos)
        class(ALCHEMYIntegral), intent(in) :: this
        integer,                intent(in) :: i, j, pos
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
        class(ALCHEMYIntegral), intent(in) :: this
        integer,                intent(in) :: i, j, k, l, m
        integer :: symmetry, MPP(4), NBP(4), NPP(4), mapped(4), npq(4), ido, block_pointer, ipqrs, mv, md, mpq, mrs, mpr

        mapped(1) = this % orbital_mapping(i)
        mapped(2) = this % orbital_mapping(j)
        mapped(3) = this % orbital_mapping(k)
        mapped(4) = this % orbital_mapping(l)

        do ido = 1, 4
            NPP(ido) = this % orbital_idx(mapped(ido))
            symmetry = this % symmetry_idx(mapped(ido))
            MPP(ido) = symmetry
            NBP(ido) = this % num_orbitals_sym(symmetry + 1)
        end do

        if (m /= 0) then
            MV = MPP(1) + MPP(2)
        else
            MV = abs(MPP(1) - MPP(2))
        end if

        if (MV /= 0) then
            MD = (MV - 1) / 2
        else
            MD = -1
        end if

        MPQ = MPP(1) - MD
        MRS = MPP(3) - MD
        MPR = this % pair_idx(MPQ) + MRS + this % num_two_electron_blocks(MV + 1)

        block_pointer = this % two_electron_pointer(mpr) - 1

        do ido = 1, 4, 2
            if (MPP(ido) /= MPP(ido + 1)) then
                NPQ(ido)     = NBP(ido) * (NPP(ido + 1) - 1) + NPP(ido)
                NPQ(ido + 1) = NBP(ido) *  NBP(ido + 1)
            else
                if (NPP(ido) < NPP(ido + 1)) then
                    NPQ(ido) = this % pair_idx(NPP(ido + 1)) + NPP(ido)
                else
                    NPQ(ido) = this % pair_idx(NPP(ido))     + NPP(ido + 1)
                end if
                NPQ(ido + 1) = this % pair_idx(NBP(ido) + 1)
            end if
        end do

        if (MPP(1) == MPP(3) .and. MPP(2) == MPP(4)) then
            if (NPQ(1) < NPQ(3)) then
                IPQRS = this % pair_idx(NPQ(3)) + NPQ(1)
            else
                IPQRS = this % pair_idx(NPQ(1)) + NPQ(3)
            end if
        else
            if (this % integral_ordering == 0) then
                IPQRS = NPQ(2) * (NPQ(3) - 1) + NPQ(1)
            else
                IPQRS = NPQ(4) * (NPQ(1) - 1) + NPQ(3)
            end if
        end if

        get_two_electron_index = IPQRS + block_pointer
        !print *,'ALC ',block_pointer,     IPQRS
        !write(stdout,"(6i8)") i,j,k,l,m,get_two_electron_index

    end function get_two_electron_index


    subroutine write_geometries_ALCHEMY (this, iounit)
        use params,          only : cpoly
        use scatci_routines, only : SEARCH

        class(ALCHEMYIntegral), intent(in) :: this
        integer,                intent(in) :: iounit
        integer :: ido, ifail

        !REWIND iounit
        !CALL SEARCH(iounit,cpoly,ifail)
        !READ(iounit)nnuc
        ! ALLOCATE(xnuc(nnuc),ynuc(nnuc),znuc(nnuc),charge(nnuc),CNAME(nnuc))
        !DO i=1, nnuc
        !         READ(iounit)cname(i), ii, xnuc(i), ynuc(i), znuc(i), charge(i)
        !END DO

    end subroutine write_geometries_ALCHEMY


    subroutine destroy_integral_ALCHEMY (this)
        class(ALCHEMYIntegral) :: this

        if (associated(this % one_electron_integral)) deallocate(this % one_electron_integral)
        if (associated(this % two_electron_integral)) deallocate(this % two_electron_integral)

        if (allocated(this % pair_labels)) deallocate(this % pair_labels)

        !> The number of labels per symmetry
        if (allocated(this % num_orbitals_sym)) deallocate(this % num_orbitals_sym)

        if (allocated(this % one_electron_pointer)) deallocate(this % one_electron_pointer)
        if (allocated(this % two_electron_pointer)) deallocate(this % two_electron_pointer)

        if (allocated(this % pair_idx)) deallocate(this % pair_idx)
        if (allocated(this % orbital_idx)) deallocate(this % orbital_idx)
        if (allocated(this % symmetry_idx)) deallocate(this % symmetry_idx)

        if (allocated(this % num_two_electron_blocks)) deallocate(this % num_two_electron_blocks)
        if (allocated(this % num_one_electron_blocks)) deallocate(this % num_one_electron_blocks)

    end subroutine destroy_integral_ALCHEMY

end module ALCHEMY_module
