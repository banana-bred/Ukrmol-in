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
!> \brief Distribute electrons to available orbitals
!>
module congen_distribution

    implicit none

    ! module entry point
    public distrb

    ! internal support routines
    private assign
    private cplea, cplem
    private cgcoef
    private getsa, getsm, getso
    private packdet
    private print1, print2, print3, print4, print5
    private state
    private wfn

contains

    !> \brief Assign quantum numbers to real shells.
    !>
    !> \param nshl    Number of shells.
    !> \param ndist   Number of distributions generated from a given set of shells.
    !> \param nused   Number of slots used for finally assigned shells.
    !> \param refcon  ?
    !> \param excon   ?
    !> \param qnstor  Finally assigned quantum numbers.
    !> \param nx      Dimension of x overlayed with \c qnstor.
    !> \param nftw    Text output file unit.
    !>
    subroutine assign (nshl, ndist, nused, refcon, excon, qnstor, nx, nftw)

        ! Reminder of some used global variables:
        !  pqn(1,nshl)   0 for pseudoshell, sequence number for real shell
        !  pqn(2,nshl)   starting index for pseudoshell
        !  pqn(3,nshl)   ending   index for pseudoshell
        !  occshl(nshl)  occupation
        !  sshl(nshl)    sym-values per shell
        !  nst(m,n)      pointer to shell index and order count
        !  nshsym(nsym)  counter for number of shells in s symmetry
        !  qnshlr(nshl)  work area stores index of shells

        use congen_data, only : ntcon, test, nrcon, nobt, nst, nshsym, nobi, ift, occshl, pqn => pqnst, sshl => sshlst, qnshlr

        integer :: ndist, nftw, nshl, nused, nx
        integer, dimension(nobt) :: excon
        integer, dimension(*) :: qnstor
        integer, dimension(nobt,max(ntcon,1)) :: refcon ! JMC changing the 2nd dimension from 2
        intent (in) nftw, nshl, nx, refcon
        intent (out) nused, qnstor
        intent (inout) excon, ndist

        integer, save :: allow, i, ic, it, ita, itest, j, ks, kss, mrun, mt, nav, navm1, ndrop, nrep, nshrun, nt

        ! workspace limits
        if (ift /= 0) then
            nav = nx
            navm1 = nav - 1
            ift = 0
        end if

        ! categorize shells by symmetry (into NST array)

        mrun = 0                    ! 1-based index of current symmetry being processed (-> number of symmetries processed)
        nshrun = 0                  ! 1-based index of current shell being processed (-> number of shells proccessed)

        do while (nshrun /= nshl)   ! process all NSHL shells
            mrun = mrun + 1
            ita = 0                 ! index of shell within symmetry (-> number of shells within symmetry)
            do i = 1, nshl
                if (mrun == sshl(i)) then
                    ita = ita + 1
                    nst(mrun,ita) = i
                end if
            end do
            nshrun = nshrun + ita
        end do

        ! initializations

        nshsym(1:mrun) = 0          ! number of shells processed so far per symmetry
        ndist = 0                   ! number of distributions generated (since the start of the program)
        ndrop = 0                   ! number of distributed shells stored (since the start of the program)
        ks = 0                      ! 1-based index of the current shell

        ! ???

    100 ks = ks + 1                 ! next shell, please!
        mt = sshl(ks)               ! get symmetry number of the current shell
        it = nshsym(mt)             ! get number of shells with the same symmetry processed before the current one
        nshsym(mt) = it + 1         ! update number of processed shells with that symmetry
        qnshlr(ks) = pqn(2,ks)      ! starting orbital of the set of orbitals this shell belongs to

        ! skip to 200 if any of the previous orbitals for this symmetry violate rules below
    120 do i = 1, it
            kss = nst(mt,i)
            if (qnshlr(kss) == qnshlr(ks)) go to 200
            if (occshl(kss) == occshl(ks) .and. &
                qnshlr(kss) >= qnshlr(ks) .and. &
                qnshlr(ks)  >= pqn(2,kss) .and. &
                qnshlr(kss) <= pqn(3,ks)) go to 200
        end do

        if (ks < nshl) go to 100

        ! take into account possible constraints
        if (ntcon /= 0) then
            excon(1:nobt) = 0
            do i = 1, nshl
                j = nobi(sshl(i)) + qnshlr(i)
                excon(j) = occshl(i)
            end do

            allow = 0
            itest = 0

            do ic = 1, ntcon
                nrep = 0
                do i = 1, nobt
                    nt = excon(i) - refcon(i,ic)
                    nrep = abs(nt) + nrep
                end do
                if (test(ic) /= 1) then
                    if (nrep <= nrcon(ic) + nrcon(ic)) go to 200
                else
                    if (nrep <= nrcon(ic) + nrcon(ic)) allow = 1
                end if
                itest = itest + test(ic)
            end do

            if (allow == 0 .and. itest /= 0) go to 200
        end if

        ! allowed assignment is to be stored - check that there is enough space
        if (ndrop + nshl >= navm1) then
            write(nftw,'("1",31("*"),/," ",31("*"),"STORAGE OVERFLOW IN ASSIGN:",I8, " WORDS AVAILABLE")') nx
            stop 70
        end if

        ! store the assignment
        ndist = ndist + 1
        do i = 1, nshl
            ndrop = ndrop + 1
            qnstor(ndrop) = qnshlr(i)
        end do

        ! ascend in assignment loops
    200 do while (ks > 0)
            mt = sshl(ks)           ! get symmetry number of the current shell
            it = nshsym(mt) - 1     ! get number of shells with the same symmetry processed before the current one
            qnshlr(ks) = qnshlr(ks) + 1
            if (qnshlr(ks) <= pqn(3,ks)) go to 120
            nshsym(mt) = it
            ks = ks - 1
        end do

        nused = ndrop

    end subroutine assign


    !> \brief Clebsch-Gordan coefficients.
    !>
    !> Evaluates all Clebsch-Gordan coefficients for given \f$ j_1, j_2, j_2 \f$ and \f$ m_3 \f$.
    !>
    !> \param j1      Angular momentum 1.
    !> \param j2      Angular momentum 2.
    !> \param j3      Angular momentum 3.
    !> \param m3      Angular momentum projection 3.
    !> \param n       On return, number of (\c m1, \c m2) pairs giving non-zero Clebsch-Gordan coefficients.
    !> \param ms      On return, pairs of projections \c m1, \c m2 for angular momenta j1 and j2.
    !> \param c       On return, array of resulting \c n Clebsch-Gordan coefficients.
    !> \param intpfg  Requested text output intensity.
    !>
    subroutine cgcoef (j1, j2, j3, m3, n, ms, c, intpfg)

        use precisn,     only : wp
        use consts,      only : xzero, xone
        use congen_data, only : binom, ind, nftw, thresh1

        integer :: intpfg, j1, j2, j3, m3, n
        real(kind=wp), dimension(*) :: c
        integer, dimension(2,*) :: ms
        intent (in) intpfg, j1, j2, j3, m3
        intent (inout) c, ms, n

        real(kind=wp) :: a, b, t
        integer :: i, i1, i2, ii, jj, js, lb, lb1, lb2, lbh, lbl, m, m1, m2
        integer, dimension(3) :: j, k, l

        js = (j1 + j2 + j3 - 1) / 2
        j(1) = j1
        j(2) = j2
        j(3) = j3
        k(1) = js - j2
        k(2) = js - j3
        k(3) = js - j1
        n = 0
        m = m3
        if (j3 - 1 < abs(m - 1)) return
        if (any(k < 0)) return
        a = xone / (binom(ind(js+1)+k(2)) * binom(ind(j3)+k(1)))
        l(3) = (j3 + m3 - 2) / 2
        m1 = j1

        do jj = 1, j1
            m2 = m3 - m1 + 1
            if (abs(m2 - 1) <= j2 - 1) then
                n = n + 1
                ms(1,n) = m1
                ms(2,n) = m2
                l(1) = (j1 - m1) / 2
                l(2) = (j2 + m2 - 2) / 2
                b = a
                do ii = 1, 3
                    b = binom(ind(j(ii)) + k(ii)) / binom(ind(j(ii)) + l(ii)) * b
                end do
                b = sqrt(b)
                i1 = max(l(1) - k(1), l(2) - k(3), 0)
                i2 = min(l(1), l(2), k(2))
                t = xzero
                if (i2 >= i1) then
                    lbl = ind(k(2) + 1) + i1
                    lbh = ind(k(2) + 1) + i2
                    lb1 = ind(k(1) + 1) + l(1) - i1
                    lb2 = ind(k(3) + 1) + l(2) - i1
                    do lb = lbl, lbh
                        t = binom(lb) * binom(lb1) * binom(lb2) - t
                        lb1 = lb1 - 1
                        lb2 = lb2 - 1
                    end do
                end if
                c(n) = b * t * (-xone)**i2
                if (abs(c(n)) <= thresh1) n = n - 1
            end if
            m1 = m1 - 2
        end do

        if (intpfg /= 0) then
            write(nftw,'(" CGCOEF : CLEBSCH-GORDAN COEFFICIENTS FOR")')
            write(nftw,'(" J1 =",I4," J2 =",I4," J3 =",I4," M3 =",I4,/)') j1, j2, j3, m3
            write(nftw, '(/(E25.15,2I5))') (c(i), ms(1,i), ms(2,i), i=1,n)
        end if

    end subroutine cgcoef


    !> \brief Form \c ne electrons in a shell coupled to (l,is).
    !>
    !> This version is for linear molecules (Alchemy).
    !>
    !> \param ne   Number of electrons.
    !> \param l    Lambda of shell.
    !> \param is   Spin of coupled shell (is=2*s+1).
    !> \param isz  2*sz+1
    !> \param m    Projection of lambda of coupled shell.
    !> \param nc   Number of determinants required for shell.
    !> \param c    Coefficient of determinant,
    !> \param iso  Spin-orbitals for determinants according to the following table
    !>             \verbatim
    !>                  0    sa        0    l+a
    !>                  1    sb        1    l+b
    !>                                 2    l-a
    !>                                 3    l-b
    !>             \endverbatim
    !>
    subroutine getsa (ne, l, is, isz, m, nc, c, iso)

        use precisn,     only : wp
        use consts,      only : one => xone
        use congen_data, only : root2

        integer :: is, isz, l, m, nc, ne
        real(kind=wp), dimension(*) :: c
        integer, dimension(*) :: iso
        intent (in) is, isz, l, m, ne
        intent (out) c, nc
        intent (inout) iso

        nc = 1
        c(1) = one
        if (ne == 0) return
        if (ne /= 4) then
            if (ne < 2) then
                iso(1) = 1 - isz / 2
                if (l == 0) return
                if (m < 0) iso(1) = iso(1) + 2
                return
            else if (ne == 2) then
                if (l == 0) then
                    iso(1) = 0
                    iso(2) = 1
                    return
                end if
                if (isz + m /= 1) then
                    if (m /= 0) then
                        iso(1) = 1 - (l + l) / m
                        iso(2) = iso(1) + 1
                        return
                    end if
                    iso(1) = (3 - isz) / 4
                    iso(2) = iso(1) + 2
                    return
                end if
                nc = 2
                iso(1) = 0
                iso(2) = 3
                iso(3) = 1
                iso(4) = 2
                c(1) = root2
                c(2) = root2
                if (is == 1) c(2) = -root2
                return
            else
                if (m < 0) then
                iso(1) = 1 - isz / 2
                iso(2) = 2
                iso(3) = 3
                return
                end if
                iso(1) = 0
                iso(2) = 1
                iso(3) = 3 - isz / 2
                return
            end if
        end if

        iso(3) = 2
        iso(4) = 3
        iso(1) = 0
        iso(2) = 1

    end subroutine getsa


    !> \brief Form \c ne electrons in a shell coupled to (l,is).
    !>
    !> This version is for non-linear molecules (Molecule).
    !>
    !> \param ne   Number of electrons
    !> \param isz  2*sz+1
    !> \param nc   Number of determinants required for shell.
    !> \param c    Coefficient of determinant.
    !> \param iso  Spin-orbitals for determinants according to the following table:
    !>             \verbatim
    !>                  0    sa        0    l+a
    !>                  1    sb        1    l+b
    !>                                 2    l-a
    !>                                 3    l-b
    !>             \endverbatim
    !>
    subroutine getsm (ne, isz, nc, c, iso)

        use precisn, only : wp
        use consts,  only : one => xone

        integer :: isz, nc, ne
        real(kind=wp), dimension(*) :: c
        integer, dimension(*) :: iso
        intent (in) isz, ne
        intent (out) c, iso, nc

        nc = 1
        c(1) = one
        iso(1) = 0
        iso(2) = 1
        if (ne == 1) iso(1) = 1 - isz / 2

        ! So the results will be:
        !
        !   isz    ne    iso(1)    iso(2)
        ! ---------------------------------
        !   any    /= 1     0       1
        !   1      == 1     1       1
        !   2      == 1     1       1
        !   3      == 1     0       1
        !   4      == 1    -1       1

    end subroutine getsm


    !> \brief ?
    !>
    !> \param ns      ?
    !> \param intpfg  Print level flag.
    !> \param nti     ?
    !> \param iqns    ?
    !> \param ci      ?
    !> \param nd      ?
    !> \param id      ?
    !> \param cd      ?
    !> \param last    ?
    !>
    subroutine getso (ns, intpfg, nti, iqns, ci, nd, id, cd, last)

        use precisn,     only : wp
        use congen_data, only : nes => occshl, ms => mshl, iqn => qnshl, ne => nnlecg, symtyp, nftw

        integer :: intpfg, last, nd, ns, nti
        real(kind=wp), dimension(*) :: cd, ci
        integer, dimension(*) :: id
        integer, dimension(2,ns,*) :: iqns
        intent (in) ci, intpfg, iqns, last, ns, nti
        intent (inout) cd, id, nd

        real(kind=wp), dimension(100) :: c ! JMC change dimension to ns*(max nc=2) ??? (an overestimate, for safety).
        real(kind=wp), dimension(ns) :: cs ! JMC changing the dimension from 50
        integer :: i, ie1, ie2, is, isz, iti, kc, ke, kso, ld, ld1, ld2, ml, nc
        integer, dimension(200) :: iso ! JMC change dimension to (max nc=2)*sum(nes(i),i=1,ns) ??? (an overestimate, for safety).
        integer, dimension(150) :: jso ! JMC change dimension to max(ne, sum(nes(i),i=1,ns)) ???
        integer, dimension(ns+1) :: lc, lso ! JMC changing the dimension from 51.  Think LSO could be (NS) not (NS+1)...
        integer, dimension(ns) :: lcs, lsos ! JMC changing the dimension from 50
        real(kind=wp) :: t

        if (intpfg /= 0) then
            write(nftw,'(" GETSO : NTI =",I10,/," IQNS :",/)') nti
            do iti = 1, nti
                write(nftw,'(20I5)') (iqns(1,is,iti), is=1,ns)
                write(nftw,'(20I5)') (iqns(2,is,iti), is=1,ns)
            end do
        end if

        nd = 0
        ld = 1

        do iti = 1, nti

            kc = 1
            kso = 1

            do is = 1, ns
                isz = iqns(1,is,iti)
                if (symtyp <= 1) then
                    ml = iqns(2,is,iti)
                    call getsa (nes(is), ms(is), iqn(1,is), isz, ml, nc, c(kc), iso(kso))
                else
                    call getsm (nes(is), isz, nc, c(kc), iso(kso))
                end if
                lc(is) = kc
                lso(is) = kso
                kc = kc + nc
                kso = kso + nc * nes(is)
            end do

            lc(ns + 1) = kc
            is = 1
            t = ci(iti)
            ie2 = 0

    300     cs(is) = t
            lsos(is) = lso(is)
            lcs(is) = lc(is)

            if (nes(is) == 0) go to 415

    400     ie1 = ie2 + 1
            ie2 = ie2 + nes(is)
            kso = lsos(is)
            do ke = ie1, ie2
                jso(ke) = iso(kso)
                kso = kso + 1
            end do

    415     lsos(is) = kso
            t = cs(is) * c(lcs(is))
            lcs(is) = lcs(is) + 1
            is = is + 1
            if (is <= ns) go to 300

            nd = nd + 1
            if (nd > last) then
                write(nftw,'("0storage overflow")')
                nd = 0
                return
            end if

            do ke = 1, ne
                id(ld) = jso(ke)
                ld = ld + 1
            end do

            cd(nd) = t

            do while (is > 1)
                is = is - 1
                ie2 = ie2 - nes(is)
                if (lcs(is) < lc(is + 1)) go to 400
            end do

        end do

        if (intpfg /= 0) then
            write(nftw,'(" GETSO : ND =",I6,/," CD, ID :",/)') nd
            ld2 = 0
            do i = 1, nd
                ld1 = ld2 + 1
                ld2 = ld2 + ne
                write(nftw,'(E25.15,20I5)') cd(i), (id(ld), ld=ld1,ld2)
            end do
        end if

    end subroutine getso


    !> \brief Loop through (and fill) all allowed couplings for a given electron distribution into shells.
    !>
    !> \verbatim
    !>    CUPSET IS THE ENTRY TO SET LOCATIONS OF ARRAYS
    !>        MSHL(NSHL)     SYMMETRY NUMBER FROM ZERO TO N-1 (MVALUE)
    !>        QNSHL(3,2*NSHL-1)
    !>                       1 -- MULT / 2 -- SYMMETRY / 3 -- +- (NOT USED)
    !>        CUP(3,NSHL-1)
    !>        QNTOT(3)       TOTAL QN'S
    !>        SPNMIN(NSHL-1) TEMP STORAGE FOR LOWEST SPIN COUPLING
    !>        X(NX)          WORK AREA (R*8)
    !>        NSHL           NUMBER OF TRUE SHELLS OCCUPIED
    !>        NSTATE         NUMBER OF COUPLINGS (COMPUTED)
    !>        NTYPE          PROTO-TYPE NUMBER (INPUT)
    !>        NDIST          NUMBER OF PQN ASSIGNMENTS
    !>        NCSF           RUNNING CSF NUMBER
    !>        SYMTYP         GE 2 FOR MOLECULE
    !>        CONPF          PRINT FLAG
    !>
    !>        GUSHELL(NSHL)  GU VALUE FOR EACH SHELL
    !>        SHLMX(NSHL)    MAX OCCUPATION FOR A SHELL
    !>        OCCSHL(NSHL)   COMPRESSED SHELL OCC'S / ALL ZEROS DELETED
    !>                       AND PSEUDO SHELLS EXPANDED
    !>        QNSHL(3,NSHL)  FIRST INDEX = 1 MULTIPLICITY FOR EACH COUPLING
    !>                       FIRST INDEX = 2 M-VALUE FOR EACH COUPLING
    !>                       FIRST INDEX = 3 +1 FOR (S+) -1 FOR (S-)
    !>                                       OTHERWIZE ZERO
    !>        CUP(3,NSHL)    COUPLING SCEME IN ORDER STATE A TO STATE B
    !>                       GIVES STATE C
    !>        QNTOT(3)       INPUT STATE TO BE SEARCHED FOR. ORDER: MULTI-
    !>                       PLICITY,ANG.MOM., PLUS(+1) OR MINUS(-1)
    !>        GUTOT          G(+1) OR U(-1) FOR INPUT STATE
    !>        KSLIM(2,NSHL)  FIRST INDEX = 1  IS NUMBER OF STATES
    !>                       ALREDY LOOPED OVER IN A SHELL
    !>                       FIRST INDEX = 2 IS THE NUMBER OF STATES FOR
    !>                       A SHELL
    !>        MCLIM(2,NSHL)  FIRST INDEX = 1 IS THE NUMBER OF STATES
    !>                       ALREADY LOOPED OVER IN A COUPLING
    !>                       FIRST INDEX = 2 IS THE MAXIMUM NUMBER OF STATES
    !>                       IN A COUPLING
    !>        QNTMP(2,3,NSHL) FIRST INDEX POINTS ON + OR
    !>        QNTMP(2,3,NSHL) FIRST INDEX POINTS ON +- AND M(2 AND 1 RESP.)
    !>                       SECOND INDEX POINTS ON POSSIBLE COUPLINGS
    !>                       (M+M,M-M,2*M,S+,S-) SAVE AREA DOWN THE LOOPS
    !>        SPNMIN(NSHL)   MIN MULTIPLICITY FOR A COUPLING
    !> \endverbatim
    !>
    subroutine cplea (nncsf, nadel, x, nftw)

        use iso_c_binding, only : c_loc, c_f_pointer
        use precisn,       only : wp
        use congen_data,   only : lg, next, nx, icdi, iexcon, indi, inodi, iqnstr, irfcon, ndel, confpf, ndist, nshl, &
                                  nstate, occshl, mshl, gushl, cup, qnshl, spnmin, kslim, mclim, shlmx1, qntot,       &
                                  qntar, gutot, nndel, cdimx, ndimx, nodimx

        integer :: nadel, nftw, nncsf
        real(wp), dimension(*), target  :: x
        integer,  dimension(:), pointer :: int_ptr, irfcon_ptr, iexcon_ptr, iqnstr_ptr, nodi_ptr, ndi_ptr, ndel_ptr, nnext_ptr
        integer,  dimension(2,3,lg)     :: qntmp

        integer :: gutry, i, idumm, iidis1, iidist, iocc, kc, kclim, kcs, kcsp, ks, &
                   m, mtry, nav, nnext, nused, shl1, shl2, spntry, z1=1

        real(wp), pointer :: x_ptr

        nstate = 0  ! number of states which may be formed with current distribution into shells
        iidis1 = 0  ! number of the current shell distribution

        ! perform gross checks on symmetries allowed for dist

        if (gutot /= 0) then
            gutry = 1
            do i = 1, nshl
                if (gushl(i) > 0) cycle
                if (iand(occshl(i), z1) /= 0) gutry = -gutry
            end do
            if (gutry /= gutot) return
        end if
        spntry = 1
        mtry = 0
        do i = 1, nshl
            m = mshl(i) + 1
            if (gutot /= 0 ) m = m + m
            iocc = occshl(i)
            if (iocc > shlmx1(m) / 2) iocc = shlmx1(m) - iocc
            spntry = spntry + iocc
            mtry = mtry + iocc * mshl(i)
        end do
        if (spntry > qntot(1)) return
        if (mtry - qntot(2) > 0 .or. iand(mtry - qntot(2), z1) /= 0) return

        ! initialize third shell coupling

        kclim = nshl - 1
        do i = 1, kclim
            qntmp(1,3,i) = 0
            qntmp(2,3,i) = -1
        end do

        ! initialize shell quantum numbers and loop limits for shell coupling

        spntry = 1
        mtry = 0
        do i = 1, nshl
            m = mshl(i) + 1
            if (gutot /= 0) m = m + m
            qnshl(1,i) = 1
            qnshl(2,i) = 0
            qnshl(3,i) = 1
            kslim(2,i) = 1
            if (occshl(i) == 0 .or. occshl(i) == shlmx1(m)) cycle
            spntry = spntry + 1
            mtry = mtry + mshl(i)
            qnshl(1,i) = 2
            qnshl(2,i) = mshl(i)
            qnshl(3,i) = 0
            if (mshl(i) == 0) qnshl(3,i) = 1
            if (occshl(i) /= 2) cycle
            kslim(2,i) = 3
            spntry = spntry - 1
            mtry = mtry - mshl(i)
        end do

        ! begin to decend into loops

        ks = 1
    100 kslim(1,ks) = 1
        if (kslim(2,ks) == 1) go to 120
        qnshl(1,ks) = 3
        qnshl(2,ks) = 0
        qnshl(3,ks) = -1
        spntry = spntry + 2
    120 ks = ks + 1
        if (ks < nshl) go to 100

        if (spntry < qntot(1) .or. mtry < qntot(2)) go to 600
        kc = 1
        kcs = nshl + 1
        if (kclim /= 0) go to 300
        if (qnshl(1,1) /= qntot(1) .or. qnshl(2,1) /= qntot(2)) go to 500
        if (qnshl(2,1) == 0 .and. qnshl(3,1) /= qntot(3)) go to 500
        go to 420
    300 shl1 = cup(1,kc)
        shl2 = cup(2,kc)
        qnshl(1,kcs) = qnshl(1,shl1) + qnshl(1,shl2) - 1
        spnmin(kc) = abs(qnshl(1,shl1) - qnshl(1,shl2)) + 1
        qntmp(1,1,kc) = qnshl(2,shl1) + qnshl(2,shl2)
        qnshl(2,kcs) = qnshl(2,shl1) + qnshl(2,shl2)
        qntmp(1,2,kc) = abs(qnshl(2,shl1) - qnshl(2,shl2))
        qntmp(2,1,kc) = qnshl(3,shl1) * qnshl(3,shl2)
        qnshl(3,kcs) = qnshl(3,shl1) * qnshl(3,shl2)
        qntmp(2,2,kc) = 0
        mclim(1,kc) = 1
        mclim(2,kc) = 1
        if (qntmp(1,1,kc) == qntmp(1,2,kc)) go to 320
        mclim(2,kc) = 2
        if (qntmp(1,2,kc) /= 0) go to 320
        mclim(2,kc) = 3
        qntmp(2,2,kc) = 1
    320 kc = kc + 1
        kcs = kcs + 1
        if (kc < kclim) go to 300

        ! test if the final coupling is in the range permitted

        kc = kclim
        kcs = kcs - 1
        if (qnshl(1,kcs) < qntot(1) .or. spnmin(kc) > qntot(1)) go to 500
        qnshl(1,kcs) = qntot(1)
        if (qnshl(2,kcs) /= qntot(2) .and. qntmp(1,2,kc) /= qntot(2)) go to 500
        qnshl(2,kcs) = qntot(2)
        qnshl(3,kcs) = qntot(3)
        if (qntmp(1,1,kc) /= 0) go to 420
        qnshl(3,kcs) = qntmp(2,1,kc)
        if (qnshl(3,kcs) /= qntot(3)) go to 500

        ! for R-matrix calculations: reject couplings which do not preserve the target quantum numbers

    420 if (qntar(1) >= 0) then
            kcsp = kcs - 1
            if (kcs == 3) kcsp = 1
            if (qnshl(1,kcsp) /= qntar(1)) go to 500
            if (qnshl(2,kcsp) /= qntar(2)) go to 500
            if (qnshl(3,kcsp) /= qntar(3)) go to 500
        end if
        nstate = nstate + 1
        if (nstate /= 1) go to 400

        ! convert some real pointers to integer pointers for use in ASSIGN
        x_ptr => x(irfcon) ; call c_f_pointer (c_loc(x_ptr), irfcon_ptr, (/1/))
        x_ptr => x(iexcon) ; call c_f_pointer (c_loc(x_ptr), iexcon_ptr, (/1/))
        x_ptr => x(iqnstr) ; call c_f_pointer (c_loc(x_ptr), iqnstr_ptr, (/1/))

        ! Assign pqn value and print type and distrib data for allowed state
        call assign (nshl, ndist, nused, irfcon_ptr, iexcon_ptr, iqnstr_ptr, nx, nftw)

        if (ndist == 0) then
            nstate = 0
            return
        end if

        nav = nx - nused
        nnext = next + nused
        if (nndel /= 0) go to 405
        call print1 (iqnstr_ptr)
    400 if (confpf >= 10 .and. nndel == 0) call print2 (0)
    405 continue

        ! convert some real pointers to integer pointers for use in WFN
        x_ptr => x(inodi)  ; call c_f_pointer (c_loc(x_ptr), nodi_ptr,   (/1/))
        x_ptr => x(indi)   ; call c_f_pointer (c_loc(x_ptr), ndi_ptr,    (/1/))
        x_ptr => x(ndel)   ; call c_f_pointer (c_loc(x_ptr), ndel_ptr,   (/1/))
        x_ptr => x(iqnstr) ; call c_f_pointer (c_loc(x_ptr), iqnstr_ptr, (/1/))
        x_ptr => x(nnext)  ; call c_f_pointer (c_loc(x_ptr), nnext_ptr,  (/1/))

        call wfn (nncsf, nadel, iidist, iidis1, nodi_ptr, ndi_ptr, x(icdi), ndel_ptr, iqnstr_ptr, x(nnext), nnext_ptr, nav)

        if (nndel > 0 .and. nadel > nndel .and. iidis1 == 0) return
        if (nndel == 0) go to 500
        if (iidist == 0) go to 500
        if (nstate == 1) call print1 (iqnstr_ptr)
        if (confpf >= 10 .and. iidis1 /= 0) call print2 (iidis1)

        ! ascend in coupling tree - ascend in the shell to shell coupling loops

    500 kc = kc - 1
        kcs = kcs - 1
        if (kc == 0) go to 600
        mclim(1,kc) = mclim(1,kc) + 1
        if (mclim(1,kc) <= mclim(2,kc)) go to 520
        qnshl(1,kcs) = qnshl(1,kcs) - 2
        if (qnshl(1,kcs) < spnmin(kc)) go to 500
        mclim(1,kc) = 1
    520 qnshl(2,kcs) = qntmp(1,mclim(1,kc),kc)
        qnshl(3,kcs) = qntmp(2,mclim(1,kc),kc)
        go to 320

        ! ascend in the loops which couple shells to themselves

    600 ks = ks - 1
        if (ks == 0) go to 900
        kslim(1,ks) = kslim(1,ks) + 1
        if (kslim(1,ks) > kslim(2,ks)) go to 600
        if (kslim(1,ks) == 3) go to 625

        ! couple to 1(2*l) and adjust spntry and mtry

        qnshl(1,ks) = 1
        spntry = spntry - 2
        qnshl(2,ks) = mshl(ks) + mshl(ks)
        mtry = mtry + mshl(ks) + mshl(ks)
        qnshl(3,ks) = 0
        go to 120

        ! couple to 1(s+) and adjust spntry and mtry

    625 mtry = mtry - qnshl(2,ks)
        qnshl(2,ks) = 0
        qnshl(3,ks) = 1
        go to 120

    900 if (nndel /= 0 .and. iidis1 /= 0) call print3 (iidis1, 1)
        if (nndel == 0 .and. nstate /= 0) call print3 (idumm, 2)

    end subroutine cplea


    !> \brief ?
    !>
    !> \param nncsf  Number of CSFs (total).
    !> \param nadel  Number of processed configurations.
    !> \param x      Real workspace (or the free part of it).
    !> \param nftw   Text output unit.
    !>
    subroutine cplem (nncsf, nadel, x, nftw, noex)

        use iso_c_binding, only : c_loc, c_f_pointer
        use precisn,       only : wp
        use global_utils,  only : mprod
        use congen_data,   only : lg, next, nx, icdi, iexcon, indi, inodi, iqnstr, irfcon, ndel, confpf, ndist, nshl, &
                                  nstate, occshl, mshl, cup, qnshl, spnmin, shlmx1, qntot, qntar, nndel

        ! Reminder of some used global variables:
        !
        !  occshl(nshl)         compressed shell occ's / all zeros deleted and pseudo shells expanded
        !  mshl(nshl)           symmetry number from zero to n-1 (mvalue)
        !  qnshl(3,2*nshl-1)    1 -- mult / 2 -- symmetry / 3 -- +- (not used)
        !  cup(3,nshl-1)
        !  qntot(3)             total qn's
        !  qntmp(2*nshl-1)      temp storage for true molecule sym values -- zeros passed to bowen in qnshl
        !  spnmin(nshl-1)       temp storage for lowest spin coupling
        !  x(nx)                real work area
        !  nshl                 number of true shells occupied
        !  nstate               number of couplings (computed)
        !  ntype                proto-type number (input)
        !  ndist                number of pqn assignments
        !  ncsf                 running csf number
        !  symtyp               >= 2 for molecule
        !  conpf                print flag

        integer :: nadel, nftw, nncsf, noex
        real(wp), dimension(*), target :: x
        integer,  dimension(:), pointer :: int_ptr, nodi_ptr, ndi_ptr, ndel_ptr, irfcon_ptr, iexcon_ptr, iqnstr_ptr, nnext_ptr
        real(wp), pointer :: x_ptr

        integer :: i, idumm, iidis1, iidist, kc, kclim, kcs, kcsp, mtry, mu1, mu2, nav, nnext, nused, spntry
        integer, dimension(2*lg) :: qntmp ! JMC increasing the dimension from LG to match 2nd dimension of QNSHL

        ! initialize

        nstate = 0                                                  ! ?
        mtry = 0                                                    ! combination of IRRs of all shells
        spntry = 1                                                  ! maximal combination of spins of all shels
        iidis1 = 0                                                  ! ?

        do i = 1, nshl                                              ! loop over all shells in the current set of orbitals
            qnshl(1,i) = 1                                          ! default to singlet (2S+1)
            qnshl(2,i) = 0                                          ! default to maximal symmetry (IRR = 1)
            qnshl(3,i) = 1                                          ! default to +1 inversion symmetry
            qntmp(i) = 0
            if (occshl(i) == shlmx1(mshl(i) + 1)) cycle             ! fully occupied shells do not contribute anything valuable -> skip
            qnshl(1,i) = 2                                          ! but shells not fully occupied are doublets
            spntry = spntry + 1                                     ! compose maximal possible total spin
            qntmp(i) = mshl(i)                                      ! the symmetry number of the shell
            mtry = mprod(mtry + 1, qntmp(i) + 1, 0, nftw) - 1       ! combine IRR of all shells so far
        end do

        if (spntry < qntot(1)) return                               ! jump out if impossible to compose the wanted total spin
        if (mtry /= qntot(2)) return                                ! jump out if the combined IRR of all shells is not the wanted one
        if (nshl == 1 .and. qnshl(1,1) /= qntot(1)) return          ! jump out if the *only* shell violates spin requirement

        ! complete M-values for coupling of shells to each other

        kc = 0
        kcs = nshl + 1

        if (nshl == 1 .and. qnshl(1,1) == qntot(1)) go to 360       ! trivial combination -- one contributing shell

        do i = 1, nshl - 1
            qntmp(kcs) = mprod(qntmp(cup(1,i)) + 1, qntmp(cup(2,i)) + 1, 0, nftw) - 1
            qnshl(2,kcs) = 0
            qnshl(3,kcs) = 1
            kcs = kcs + 1
        end do

        ! begin to descend into shell to shell coupling tree

        kcs = nshl

    300 do while (kc < nshl - 1)
            kc = kc + 1
            kcs = kcs + 1
            mu1 = qnshl(1,cup(1,kc))
            mu2 = qnshl(1,cup(2,kc))
            qnshl(1,kcs) = mu1 + mu2 - 1
            spnmin(kc) = abs(mu1 - mu2) + 1
        end do

        ! test for allowed state

        if (qnshl(1,kcs) < qntot(1) .or. spnmin(kc) > qntot(1)) go to 500
        qnshl(1,kcs) = qntot(1)

        ! for R-matrix calculations: reject couplings which do not preserve the target quantum numbers

    360 if (qntar(1) >= 0) then
            kcsp = kcs - 1
            if (kcs == 3) kcsp = 1
            if (qnshl(1,kcsp) /= qntar(1)) go to 500
        end if

        nstate = nstate + 1

        ! first iteration: assign pqn value and print type and distrib data for allowed state

        if (nstate == 1) then
            iexcon = 1

            ! convert some real pointers to integer pointers for use in ASSIGN 
            x_ptr => x(irfcon) ; call c_f_pointer (c_loc(x_ptr), irfcon_ptr, (/1/))
            x_ptr => x(iexcon) ; call c_f_pointer (c_loc(x_ptr), iexcon_ptr, (/1/))
            x_ptr => x(iqnstr) ; call c_f_pointer (c_loc(x_ptr), iqnstr_ptr, (/1/))

            call assign (nshl, ndist, nused, irfcon_ptr, iexcon_ptr, iqnstr_ptr, nx, nftw)

            if (ndist == 0) then
                nstate = 0
                return
            end if

            nav = nx - nused
            nnext = next + nused
            if (nndel == 0) call print1 (iqnstr_ptr)
        end if

        if (confpf >= 10 .and. nndel == 0) call print2 (0)

        ! convert some real pointers to integer pointers for use in WFN
        x_ptr => x(inodi)  ; call c_f_pointer (c_loc(x_ptr), nodi_ptr,   (/1/))
        x_ptr => x(indi)   ; call c_f_pointer (c_loc(x_ptr), ndi_ptr,    (/1/))
        x_ptr => x(ndel)   ; call c_f_pointer (c_loc(x_ptr), ndel_ptr,   (/1/))
        x_ptr => x(iqnstr) ; call c_f_pointer (c_loc(x_ptr), iqnstr_ptr, (/1/))
        x_ptr => x(nnext)  ; call c_f_pointer (c_loc(x_ptr), nnext_ptr,  (/1/))

        call wfn (nncsf, nadel, iidist, iidis1, nodi_ptr, ndi_ptr, x(icdi), ndel_ptr, iqnstr_ptr, x(nnext), nnext_ptr, nav)

        if (nndel > 0 .and. nadel > nndel .and. iidis1 == 0) return

        if (nndel /= 0 .and. iidist /= 0) then
            if (nstate == 1) call print1 (iqnstr_ptr)
            if (confpf >= 10 .and. iidis1 /= 0) call print2 (iidis1)
        end if

        ! ascend in coupling tree

    500 do while (kc > 1)
            kc = kc - 1
            kcs = kcs - 1
            qnshl(1,kcs) = qnshl(1,kcs) - 2
            if (qnshl(1,kcs) >= spnmin(kc)) go to 300
        end do

        if (nndel /= 0 .and. iidis1 /= 0) call print3 (iidis1, 1)
        if (nndel == 0) call print3 (idumm, 2)

    end subroutine cplem


    !> \brief Distribute electrons to spin-orbitals.
    !>
    !> This subroutine will generate all possible configuration state functions subject to setup and constraints set
    !> by the input namelist. All generated configurations are normally stored in the global arrays \c nodi,
    !> \c ndi and \c cdi, but if the storage is about to overflow (which is very likely for
    !> large molecules), it is flushed to the output file sooner (see \ref congen_distribution::wfn), providing space for
    !> further configurations.
    !>
    !> \param nelecp        Number of electrons per set of orbitals.
    !> \param nshlp         Number of shells in a set of orbitals.
    !> \param shlmx         Max occupation of a shell or a pseudoshell.
    !> \param occshl        On return, occupation of shells / pseudoshells.
    !> \param nslsv         Save area for index, per shell.
    !> \param kdssv         ?
    !> \param loopf         Pointer for zero occup or pseudoshell, per shell.
    !> \param ksss          Save area for index, per shell.
    !> \param pqn           Shell data (triplet per shell): 1st == 0 (pseudoshell) or /= 0 (real shell index),
    !>                      2nd start index for pseudoshell, 3rd end index for pseudoshell.
    !> \param occst         ?
    !> \param shlmx1        ?
    !> \param pqnst         ?
    !> \param mshl          Quantum number of a shell or of a pseudoshell, per shell.
    !> \param mshlst        ?
    !> \param gushl         Gerade/ungerade flag, per shell.
    !> \param gushst        ?
    !> \param cup           Coupling scheme.
    !> \param cupst         Expanded shells with zeroes deleted.
    !> \param ndprod        Number of sets of orbitals.
    !> \param symtyp        Symmetry group flag (< 2 diatomic, >= 2 abelian).
    !> \param confpf        ?
    !> \param sshl          ?
    !> \param sshlst        ?
    !> \param ncsf          ?
    !> \param nncsf         ?
    !> \param nadel         ?
    !> \param nndel         ?
    !> \param x             Real workspace.
    !> \param nftw          Unit for text output of the program.
    !>
    subroutine distrb (nelecp, nshlp, shlmx, occshl, nslsv, kdssv, loopf, ksss, pqn, occst, shlmx1, pqnst, mshl, mshlst,    &
                      gushl, gushst, cup, cupst, ndprod, symtyp, confpf, sshl, sshlst, ncsf, nncsf, nadel, nndel, x, nftw, &
                       noex)

        use precisn,     only : wp
        use congen_data, only : ntyp, ndist, nstate, ncsft => ncsf, confpt => confpf, nshrun => nshl

        integer :: confpf, nadel, ncsf, ndprod, nftw, nncsf, nndel, symtyp
        integer, dimension(3,*) :: cup, cupst, pqn, pqnst
        integer, dimension(*) :: gushl, gushst, kdssv, ksss, loopf, mshl, mshlst, nelecp, nshlp, nslsv, occshl, &
                                occst, shlmx, shlmx1, sshl, sshlst
        real(kind=wp), dimension(*) :: x
        intent (in) confpf, cup, gushl, mshl, ncsf, ndprod, nelecp, nndel, nshlp, pqn, shlmx, shlmx1, sshl, symtyp
        intent (out) gushst, mshlst, pqnst, sshlst
        intent (inout) cupst, kdssv, ksss, loopf, nslsv, occshl, occst

        integer :: i, ibias, ic1, id, idd, it, it1, it2, ita, j, kds, kdsb, kdst, kprod, krun, ksi, kss, nadd, nadd2, ncrun, &
                   nela, neleft, ninitx, nshlw, nslots, noex

        confpt = confpf                 ! Print flag.
        ntyp   = 0                      ! Prototype number ...?
        nstate = 0                      ! Global number of couplings computed (later) in "cplea" or "cplem".
        ndist  = 0                      ! Number of distributions generated from set of shells (set in "cplea"/"cplem" via "assign").
        ibias  = 0                      ! Offset in the list of shells.
        kprod  = 0                      ! Set of orbitals currently being processed.
        ninitx = sum(nshlp(1:ndprod))   ! Total number of shells in all sets of orbitals.

        orbital_set_loop: do

            kprod = kprod + 1                               ! Move on to the next set of orbitals.
            kds = 0                                         ! Shell index within the set of orbitals.
            neleft = nelecp(kprod)                          ! (Remaining) number of movable electrons in the current set of orbitals.
            nshlw = nshlp(kprod)                            ! Number of shells in current set of orbitals.
            nslots = sum(shlmx(ibias + 1 : ibias + nshlw))  ! Total number of spin-orbitals in the current set of orbitals.
            occshl(ibias + 1 : ibias + nshlw) = 0           ! All shells in this set are un-occupied at the beginning.

            shell_loop: do

                kds = kds + 1                               ! Move on to the next shell (orbital) within the current set of orbitals.
                kdsb = kds + ibias                          ! Global index of the current shell.
                occshl(kdsb) = min(neleft, shlmx(kdsb))     ! Put all available electrons in that shell, if possible.

                neleft = neleft - occshl(kdsb)              ! Electrons that still did not fit into any shell.
                nslots = nslots - shlmx(kdsb)               ! Number of free spin-orbitals for redistribution of remaining electrons.

                if (neleft /= 0) then
                    if (nslots >= neleft) cycle shell_loop  ! Enough room to redistribute the remaining electrons.
                    do while (nslots <= neleft .and. kds /= 0)  ! Repeat until there is either enough space, or all shells are empty.
                        neleft = neleft + occshl(kdsb)      ! Get back electrons from the last populated shell.
                        occshl(kdsb) = 0                    ! Its occupancy is now zero.
                        nslots = nslots + shlmx(kdsb)       ! The number of slots accordingly rises.
                        kds = kds - 1                       ! Return before the last shell.
                        kdsb = kdsb - 1                     ! Return before the last shell (global index).
                    end do
                    go to 130
                else if (kprod /= ndprod) then              ! Not all sets of orbitals processed?
                    nslsv(kprod) = nslots                   ! Store the number of free spin-orbitals (of the current set of orbitals).
                    kdssv(kprod) = kds                      ! Store the last processed shell (of the current set of orbitals).
                    ibias = ibias + nshlp(kprod)            ! Move the global shell index offset to the beginning of the next set of orbitals.
                    cycle orbital_set_loop                  ! Next set of orbitals, please!
                end if

                ! No electrons left and all sets of orbitals processed here...

                where (occshl(1:ninitx) == 0 .or. pqn(1,1:ninitx) /= 0)
                    loopf(1:ninitx) = 1                     ! Un-occupied or movable-to real shells.
                elsewhere
                    loopf(1:ninitx) = 0                     ! Occupied and non-movable shells.
                end where

                ksss(1:ninitx) = 0                          ! 
                kdst = 0                                    ! 
                ksi = 0                                     ! Orbital set index.

                ! ----------------------------------------------------------------------------------------------------------
                ! expansion of pseudoshells and delete of zeroes follow the "st" arrays but for cupst are set up

                ksi_loop: do

                    if (ksi == ninitx) then

                        ! expand and compress the coupling scheme
                        nshrun = ninitx
                        ncrun = nshrun - 1
                        cupst(1:3,1:ncrun) = cup(1:3,1:ncrun)
                        krun = 0

                        all_shell_loop: do i = 1, ninitx
                            krun = krun + 1
                            if (ksss(i) < 1) then
                                ! delete section
                                do id = 1, ncrun
                                    if (cupst(1,id) == krun) then ; it1 = cupst(2,id) ; it2 = cupst(3,id) ; exit ; end if
                                    if (cupst(2,id) == krun) then ; it1 = cupst(1,id) ; it2 = cupst(3,id) ; exit ; end if
                                end do

                                ncrun = ncrun - 1

                                do idd = id, ncrun
                                    cupst(1:3,idd) = cupst(1:3,idd+1)
                                end do

                                do idd = 1, ncrun
                                    if (cupst(1,idd) == it2) then ; cupst(1,idd) = it1 ; exit ; end if
                                    if (cupst(2,idd) == it2) then ; cupst(2,idd) = it1 ; exit ; end if
                                end do

                                do idd = 1, ncrun
                                    do j = 1, 3
                                        if (cupst(j,idd) <= krun) cycle
                                        cupst(j,idd) = cupst(j,idd) - 1
                                        if (cupst(j,idd) >= it2) cupst(j,idd) = cupst(j,idd) - 1
                                    end do
                                end do

                                nshrun = nshrun - 1
                                krun = krun - 1
                            else if (ksss(i) > 1) then
                                ! add section
                                nadd = ksss(i) - 1
                                nadd2 = nadd + nadd
                                do idd = 1, ncrun
                                    cupst(3,idd) = cupst(3,idd) + nadd2
                                    do j = 1, 2
                                        if (cupst(j,idd) == krun) then
                                            cupst(j,idd) = nshrun + nadd2
                                        else if (cupst(j,idd) > krun .and. cupst(j,idd) <= nshrun) then
                                            cupst(j,idd) = cupst(j,idd) + nadd
                                        else if (cupst(j,idd) > krun .and. cupst(j,idd) > nshrun) then
                                            cupst(j,idd) = cupst(j,idd) + nadd2
                                        end if
                                    end do
                                end do
                                ic1 = krun
                                nshrun = nshrun + nadd
                                do idd = 1, nadd
                                    cupst(1,ncrun+idd) = ic1
                                    cupst(2,ncrun+idd) = krun + idd
                                    cupst(3,ncrun+idd) = nshrun + idd
                                    ic1 = nshrun + idd
                                end do
                                ncrun = ncrun + nadd
                                krun = krun + nadd
                            end if
                        end do all_shell_loop ! end of add/delete loop

                        ! now sort the expanded coupling array
                        do i = 1, ncrun - 1
                            do while (i /= cupst(3,i) - kdst)
                                it = cupst(3,i) - kdst
                                do j = 1, 3
                                    ita = cupst(j,it)
                                    cupst(j,it) = cupst(j,i)
                                    cupst(j,i) = ita
                                end do
                            end do
                        end do

                        ! generate a new wave function prototpe
                        ntyp = ntyp + 1
                        if (symtyp < 2) then
                            call cplea (nncsf, nadel, x, nftw)
                        else
                            call cplem (nncsf, nadel, x, nftw, noex)
                        end if
                        if (nstate == 0) ntyp = ntyp - 1

                        ! end of coupling scheme section
                        if (loopf(ksi) == 0) then
                            go to 320
                        else
                            kdst = kdst - ksss(ksi)
                            go to 360
                        end if
                    end if

                    ksi = ksi + 1                                       ! Next set of orbitals, please.
                    kss = 1                                             ! Shell index.
                    if (occshl(ksi) == 0) cycle ksi_loop                ! Do not process un-occupied shells.
                    if (occshl(ksi) /= 0) continue

                    nela = occshl(ksi)                                  ! Occupation of current (ksi) shell.
                    kss = 0                                             ! ... ???

                410 kss_loop: do while (kss <= pqn(3,ksi) - pqn(2,ksi)) ! Loop over prescribed orbitals in the current set.
                        kss = kss + 1                                   ! Next orbital within ksi, please.
                        kdst = kdst + 1                                 ! ... ???
                        occst(kdst) = min(nela, shlmx1(sshl(ksi)))      ! 
                        if (kss /= 1 .and. kdst > 1) occst(kdst) = min(nela, occst(kdst-1))
                        pqnst(1:3,kdst) = pqn(1:3,ksi)
                        mshlst(kdst) = mshl(ksi)
                        gushst(kdst) = gushl(ksi)
                        sshlst(kdst) = sshl(ksi)
                        nela = nela - occst(kdst)
                        if (nela /= 0) then
                            cycle kss_loop
                        else
                            ksss(ksi) = kss
                            cycle ksi_loop
                        end if
                    end do kss_loop

                    ! second part of the expansion of pseudoshells and deletion of zeros
                320 unkss_loop: do
                        occst(kdst) = occst(kdst) - 1
                        nela = nela + 1
                        if (occst(kdst) /= 0) go to 410
                        kss = kss - 1
                        kdst = kdst - 1
                        if (kss >= 1) then
                            cycle unkss_loop
                        else
                            nela = 0
                            go to 360
                        end if
                    end do unkss_loop

                360 unksi_loop: do
                        ksi = ksi - 1
                        if (ksi <= 0) exit ksi_loop  ! sic
                        if (loopf(ksi) == 0) then
                            kss = ksss(ksi)
                            go to 320
                        end if
                        kdst = kdst - ksss(ksi)
                    end do unksi_loop

                    ! end of pseudoshell expansion loops

                end do ksi_loop

                ! ----------------------------------------------------------------------------------------------------------

                ! last part of outer loop follows

            120 if (occshl(kdsb) == 0) go to 140
                occshl(kdsb) = occshl(kdsb) - 1                     ! Remove one electron from the current shell.
                neleft = 1                                          ! Remember that we have it
                if (nndel /= 0 .and. nadel > nndel) exit orbital_set_loop
                if (kds < nshlp(kprod)) cycle shell_loop            ! If we can push it into the NEXT shell, do it.
                neleft = neleft + occshl(kdsb)                      ! Otherwise, get back all electrons from the last shell.
                occshl(kdsb) = 0                                    ! Clear its population.
                nslots = nslots + shlmx(kdsb)                       ! Update available number of slots.
                kds = kds - 1                                       ! Move back one shell.
                kdsb = kdsb - 1                                     ! Global index, too.

            130 unkds_loop: do while (kds /= 0)                     ! Reset preceding shells one by one...
                    if (occshl(kdsb) /= 0) then                     ! If shell occupied:
                        occshl(kdsb) = occshl(kdsb) - 1             !  - remove one electron
                        neleft = neleft + 1                         !  - put among the unused electrons
                        if (nndel /= 0 .and. nadel > nndel) exit orbital_set_loop
                        cycle shell_loop                            !  - try to add the electrons to shells
                    end if
                    nslots = nslots + shlmx(kdsb)                   ! If not occupied, add these shells to list of known slots.
                    kds = kds - 1                                   ! Recede by one more shell...
                    kdsb = kdsb - 1                                 ! Global index too.
                end do unkds_loop

            140 kprod = kprod - 1                                   ! Recede by one group of orbitals.
                if (kprod == 0) exit orbital_set_loop               ! If at beginning, exit subroutine.
                ibias = ibias - nshlp(kprod)                        ! Decrease global shell index offset.
                nslots = nslsv(kprod)                               ! Get (previously stored) number of free slots.
                kds = kdssv(kprod)                                  ! Get (previously stored) index of the last shell in that set.
                kdsb = kds + ibias                                  ! Get its global index.
                go to 120

            end do shell_loop

        end do orbital_set_loop

        ncsft = ncsf

    end subroutine distrb


    !> \brief Copy determinant data.
    !>
    !> Whatever was this function originally designed for, it just straghtforwardly copies arrays of data.
    !>
    !> \param iqn  Source data for \c jqn.
    !> \param ni   Second rank of \c iqn.
    !> \param jqn  Destination for data from \c iqn.
    !> \param nj   Second rank of \c jqn.
    !> \param ci   Source data for \c cj.
    !> \param cj   Destination for data from \c ci.
    !> \param nt   Minimal dimension of of \c ci and \c cj and of the third rank of \c iqn and \c jqn.
    !>
    subroutine packdet (iqn, ni, jqn, nj, ci, cj, nt)

        use precisn, only : wp

        integer :: ni, nj, nt
        real(kind=wp), dimension(*) :: ci, cj
        integer, dimension(2,ni,*) :: iqn
        integer, dimension(2,nj,*) :: jqn
        intent (in) ci, iqn, ni, nj, nt
        intent (out) cj, jqn

        jqn(1:2,1:nj,1:nt) = iqn(1:2,1:nj,1:nt)
        cj(1:nt) = ci(1:nt)

    end subroutine packdet


    !> \brief ?
    !>
    subroutine print1 (qnstor)

        use global_utils, only : getin
        use congen_data , only : confpf, ndist, ne, nshl, ntyp, occshl, pqnr => pqnst, mshl, gushl, qnshl, &
                                 nnelcg => nnlecg, symtyp, nftw, blnk43, header, lp, star, nitem
        use congen_pagecontrol, only : addl, newpg, space, taddl, taddl1

        integer, dimension(*), intent(in) :: qnstor


        integer :: i, ii, imax, j, jj, k, kf, ki, klabel, lt, lta, nlex, nstar
        character(len=4), dimension(3,4) :: label = reshape((/ ' NTY', 'P=XX', 'X   ', &
                                                               ' NDI', 'ST=Y', 'YYY ', &
                                                               ' NST', 'ATE=', 'ZZZ ', &
                                                               '    ', '    ', '    '  /) , (/ 3, 4/))
        character(len=16), parameter :: frmt = '(3A4,A8,I8,8I12)'
        character(len=1), dimension(12,4) :: labell

        equivalence(label(1,1), labell(1,1))

        ! print type(confpf.ge.1) and distribution data (confpf.ge.20)
        ne = nnelcg
        if (confpf < 1) return
        nlex = 0
        if (symtyp == 1) nlex = 1
        nstar = 20 + min(nitem, nshl) * 12

        klabel = 1
        call getin (ntyp, 3, labell(7,1), 0)

        if (confpf <= 1 .and. ntyp /= 1) then
            call taddl1 (5 + (6 + nlex) * ((nshl + nitem - 1) / nitem), lta)
            if (lta /= 0) then
                call space (2)
            else
                call newpg
            end if
        else
            call newpg
        end if

        do i = 1, nshl, nitem
            imax = min(i + nitem - 1, nshl)
            call addl (5 + nlex)
            write(nftw,frmt) (label(j,klabel), j=1,3), header(1), (j, j=i,imax)
            write(nftw,frmt) blnk43, header(2), (occshl(j), j=i,imax)
            write(nftw,frmt) blnk43, header(3), (mshl(j), j=i,imax)
            if (symtyp  == 1) write(nftw,frmt) blnk43, header(4), (gushl(j), j=i,imax)

            ! print occshl and sym data

            klabel = 4
            write(nftw,'(3A4,A8,9(A3,2(I3,","),I3,")"))') blnk43, header(5), (lp, (pqnr(jj,j), jj=1,3), j=i,imax)
            if (symtyp < 2) then
                write(nftw,'(3A4,A8,9(A3,2(I2,","),I2,")"))') blnk43, header(6), (lp, (qnshl(jj,j), jj=1,3), j=i,imax)
            else
                write(nftw,frmt) blnk43, header(6), (qnshl(1,j), j=i,imax)
            end if
            call space (1)
        end do

        ! bypass printing of distributions if not required

        if (confpf >= 20) then

            ! print row of star separator

            call taddl (1, lt)
            if (lt > 0) write(nftw,'(" ",132A1)') (star, j=1,nstar)
            call space (1)
            kf = 0
            do ii = 1, ndist
                klabel = 2
                call getin (ii, 4, labell(8,2), 0)
                do i = 1, nshl, nitem
                    imax = min(i + nitem - 1, nshl)
                    call addl (4 + nlex)
                    write(nftw,frmt) (label(j,klabel), j=1,3), header(1), (j, j=i,imax)
                    write(nftw,frmt) blnk43, header(2), (occshl(j), j=i,imax)
                    write(nftw,frmt) blnk43, header(3), (mshl(j), j=i,imax)
                    if (symtyp == 1) write(nftw,frmt) blnk43, header(4), (gushl(j), j=i,imax)
                    klabel = 4
                    ki = kf + 1
                    kf = ki + (imax - i)
                    write(nftw,frmt) blnk43, header(5), (qnstor(k), k=ki,kf)
                    call space (1)
                end do
            end do

        end if

        ! merge confpf >= 1 and confpf >= 20 paths

        call addl (1)
        write(nftw,'(" ",19("*"),5X,"TOTAL NUMBER OF DISTRIBUTIONS FOR NTYP =",I3," IS",I5)') ntyp, ndist
        if (confpf < 10) return

        ! prepare for state printing

        call space (1)
        call taddl (1, lt)
        if (lt > 0) write(nftw,'(" ",132A1)') (star, j=1,nstar)
        call space (1)

    end subroutine print1


    !> \brief ?
    !>
    subroutine print2 (iidis1)

        use global_utils, only : getin
        use congen_data,  only : confpf, ncsf, ndist, nshl, nstate, cup, qnshl, symtyp, nftw, blnk43, header, lp, nitem
        use congen_pagecontrol, only : addl, space

        integer :: iidis1
        intent (in) iidis1

        integer :: i, imax, j, jj, kf, ki, klabel, ncsff, ncsfi, nshlm1
        character(len=4), dimension(3,4) :: label = reshape((/' NTY', 'P=XX', 'X   ', &
                                                              ' NDI', 'ST=Y', 'YYY ', &
                                                              ' NST', 'ATE=', 'ZZZ ', &
                                                              '    ', '    ', '    '/) , (/ 3, 4/))
        character(len=1), dimension(12,4) :: labell

        equivalence(label(1,1), labell(1,1))

        ! print state data

        nshlm1 = nshl - 1
        if (confpf < 10) return
        call getin (nstate, 3, labell(9,3), 0)
        if (nshlm1 <= 0) then
            call addl (1)
            write(nftw,'(3A4,A8,I8,8I12)') (label(j,3), j=1,3)
            return
        end if
        kf = nshl
        klabel = 3
        do i = 1, nshlm1, nitem
            imax = min(i + nitem - 1, nshlm1)
            call addl (2) ! JMC argument should probably be 3 as there are 3 writes below???
            write(nftw,'(3A4,A8,9(A3,2(I2,","),I2,")"))') (label(j,klabel), j=1,3), header(7), (lp, (cup(jj,j), jj=1,3), j=i,imax)
            klabel = 4
            ki = kf + 1
            kf = ki + (imax - i)
            write(nftw,'(3A4,A8,9(A3,2(I2,","),I2,")"))') blnk43, header(6), (lp, (qnshl(jj,j), jj=1,3), j=ki,kf)
            if (symtyp < 2) then
                write(nftw,'(3A4,A8,9(A3,2(I2,","),I2,")"))') blnk43, header(6), (lp, (qnshl(jj,j), jj=1,3), j=ki,kf)
            else
                write(nftw,'(3A4,A8,I8,8I12)') blnk43, header(6), (qnshl(1,j), j=ki,kf)
            end if
            call space (1)
        end do

        ! print CSF numbers for this state

        if (iidis1 /= 0) then
            ncsfi = ncsf - iidis1 + 1
            ncsff = ncsf
        else
            ncsfi = ncsf + 1
            ncsff = ncsf + ndist
        end if

        call addl (1)
        write(nftw,'(" ",19("*"),5X,"CSF NUMBERS",I6," TO",I6," GENERATED FOR NSTATE=",I3)') ncsfi, ncsff, nstate
        call space (1)

    end subroutine print2


    !> \brief ?
    !>
    subroutine print3 (iidis1, i13)

        use congen_data, only : confpf, ncsf, ndist, nshl, nstate, ntyp, nftw, star, nitem
        use congen_pagecontrol, only : addl, space, taddl

        integer, intent(in) :: i13, iidis1

        integer :: i, j, lt, ncsfi, nstar, nstot

        ! print summary data
        if (i13 /= 2) then
            nstot = iidis1
        else
            if (confpf < 1) return
            nstot = nstate * ndist
        end if
        ncsfi = ncsf + 1 - nstot
        if (confpf >= 40) call space (1)
        call addl (2)

        write(nftw,'(" ",19("*"),5X,"TOTAL NUMBER OF STATES FOR NTYP=",I3," IS",I4)') ntyp, nstate
        write(nftw,'(20X,"CSF NUMBERS",I10," TO",I10," (",I9," CSFS) GENERATED")') ncsfi, ncsf, nstot

        nstar = 20 + min(nitem, nshl) * 12 ! JMC see PRINT1 to understand where this hardwiring has come from
        do i = 1, 2
            call taddl (1, lt)
            if (lt > 0) write(nftw,'(" ",132A1)') (star, j=1,nstar)
        end do

    end subroutine print3


    !> \brief ?
    !>
    subroutine print4 (nd,x,ix)

        use precisn,     only : wp
        use congen_data, only : ne, nftw
        use congen_pagecontrol, only : addl, space

        integer :: nd
        integer, dimension(*) :: ix
        real(kind=wp), dimension(*) :: x
        intent (in) ix, nd, x

        integer :: ie, ind, ip, ipi, it, j, nl, nitem

        ! print determinant and spin-orbital data

        ipi = 0
        it = 1
        nitem = 20
        nl = (ne + nitem - 1) / nitem

        do ind = 1, nd
            call addl (nl + it)
            if (it /= 0) then
                write(nftw,'(30X,"NUMBER OF DET IS",I3)') nd
            it = 0
            end if
            do ie = 1, ne, nitem
                ip = ipi + 1
                ipi = ip + min(nitem - 1, ne - ie)
                if (ie /= 1) then
                    write(nftw,'(45X,20I4)') (ix(j), j=ip,ipi)
                else
                    write(nftw,'(I25,5X,E15.8,20I4)') ind, x(ind), (ix(j), j=ip,ipi)
                end if
            end do
        end do

        call space (1)

    end subroutine print4


    !> \brief ?
    !>
    subroutine print5 (nd, ndi)

        use congen_data, only : ncsf, nftw
        use congen_pagecontrol, only : addl

        integer :: nd
        integer, dimension(*) :: ndi
        intent (in) nd, ndi

        integer :: ind, ip, ipi, j, nrep

        ip = 1
        do ind = 1, nd
            nrep = ndi(ip)
            if (nrep == 0) then
                call addl (1)
                write(nftw, '(35X,2I5)') ind, nrep
                if (ind == 1) then
                    write(nftw,'("+",29X,I5)') ncsf
                end if
            else
                call addl (2)
                ipi = ip + 1
                ip = ip + nrep

                if (ind /= 1) then
                    write(nftw,'(35X,2I5,21I4/(45X,21I4))') ind, nrep, (ndi(j), j=ipi,ip)
                else
                    write(nftw,'(30X,3I5,21I4/(45X,21I4))') ncsf, ind, nrep, (ndi(j), j=ipi,ip)
                end if

                ipi = ip + 1
                ip = ip + nrep

                write(nftw,'(45X,21I4)') (ndi(j), j = ipi,ip)
            end if
            ip = ip + 1
        end do

    end subroutine print5


    !> \brief ?
    !>
    !> \param ns      Actually, this is just reference to the global \ref congen_data::nshl (number of shells in the set of orbitals
    !>                currently being processed.)
    !> \param x       Real workspace (or the currently free part of it).
    !> \param last    Number of available elements in the workspace.
    !> \param nd      ?
    !> \param confpf  Print level flag.
    !>
    subroutine state (ns, x, last, nd, confpf)

        use iso_c_binding, only : c_loc, c_f_pointer
        use precisn,       only : wp
        use congen_data,   only : lratio, cup, iqn => qnshl, ne => nnlecg, isz => nisz

        integer :: confpf, last, nd, ns
        real(kind=wp), dimension(*), target :: x
        intent (in) confpf, last
        intent (inout) nd, x

        integer :: ic, id, intpfg, iqns, jqns, lc, ld, ldp, ldq, nam, ndmx, nt, ntmx

        real(wp), pointer :: x_ptr
        integer,  pointer, dimension(:) :: iqns_ptr, jqns_ptr, ld_ptr

        x_ptr => x(1)

        intpfg = merge(1, 0, confpf > 40)

        nd     = 0
        nam    = ns + ns - 1
        ntmx   = last / (nam + nam + 1)
        ic     = 1
        iqns   = ic + ntmx

        ! convert real pointer to integer pointer for use in WFCPLE
        x_ptr => x(iqns) ; call c_f_pointer (c_loc(x_ptr), iqns_ptr, (/1/))

        call wfcple (nam, iqn, isz, cup, iqns_ptr, x(ic), ntmx, nt, intpfg)

        if (nt == 0) return

        jqns = last - 2 * ns * nt + 1
        lc = jqns - nt

        ! convert real pointer to integer pointer for use in PACK
        x_ptr => x(jqns) ; call c_f_pointer (c_loc(x_ptr), jqns_ptr, (/1/))

        call packdet (iqns_ptr, nam, jqns_ptr, ns, x(1), x(lc), nt)

        ndmx = (lc - 1) / (ne + 1)
        ld = 1 + ndmx

        ! convert real pointer to integer pointer for use in GETSO
        x_ptr => x(ld)   ; call c_f_pointer (c_loc(x_ptr), ld_ptr, (/1/))

        call getso (ns, intpfg, nt, jqns_ptr, x(lc), nd, ld_ptr, x(ic), ndmx)

        if (nd == 0) return

        ! Finalize

        ldp = nd + 1
        ldq = nd + (ne * nd + lratio - 1) / lratio
        do id = ldp, ldq
            x(id) = x(ld)
            ld = ld + 1
        end do

    end subroutine state


    !> \brief ?
    !>
    subroutine wfcple (nam, iqn, isz, icup, iqns, c, last, lc2, intpfg)

        use precisn,     only : wp
        use consts,      only : one => xone
        use congen_data, only : root2, nftw, lg

        integer :: intpfg, isz, last, lc2, nam
        real(kind=wp), dimension(*) :: c
        integer, dimension(3,lg) :: icup
        integer, dimension(3,*) :: iqn
        integer, dimension(2,nam,*) :: iqns
        intent (in) icup, isz, last, nam
        intent (inout) c, iqns, lc2

        integer :: i, iam, ic, init, j, l, lc, lc1, lc3, m, mp, ms, n, n1, n2, n3, nc, niam, niam1
        logical, dimension(nam) :: ind ! JMC changing the dimension from 150
        integer, dimension(2,200) :: iszt ! JMC not sure how to dimension this...
        real(kind=wp) :: sign

        iqns(1,nam,1) = isz
        iqns(2,nam,1) = iqn(2,nam)
        c(1) = one
        lc2 = 1
        if (nam == 1) return
        do i = 1, nam
            ind(i) = .true.
        end do
        niam = (nam + 1) / 2
    100 n3 = nam
    110 if (ind(n3)) then
            if (n3 <= niam) go to 299
            ind(n3) = .false.
            init = last - lc2 + 1
            l = last
            lc = lc2
            do while (lc > 0)
                iqns(1:2,1:nam,l) = iqns(1:2,1:nam,lc)
                c(l) = c(lc)
                lc = lc - 1
                l = l - 1
            end do
            n = 1
    300     if (icup(3,n) == n3) then
                n1 = icup(1,n)
                n2 = icup(2,n)
                if (n1 > n3 .or. n2 > n3) go to 299
                if (.not.ind(n1) .or. .not.ind(n2)) go to 299
                if (n1 <= niam) ind(n1) = .false.
                if (n2 <= niam) ind(n2) = .false.
                lc2 = 0
                do l = init, last
                    iqns(2,n1,l) = iqn(2,n1)
                    iqns(2,n2,l) = iqn(2,n2)
                    m = iqns(2,n3,l)
                    if (abs(m) /= iqn(2,n1) + iqn(2,n2)) then
                        mp = iqn(2,n1) - iqn(2,n2)
                        if (abs(m) /= abs(mp)) go to 511
                        n = n2
                        if (m /= mp) n = n1
                        iqns(2,n,l) = -iqns(2,n,l)
                    else if (iqns(2,n3,l) < 0) then
                        iqns(2,n1,l) = -iqns(2,n1,l)
                        iqns(2,n2,l) = -iqns(2,n2,l)
                    end if
                    lc1 = lc2 + 1
                    ms = iqns(1,n3,l)
                    call cgcoef (iqn(1,n1), iqn(1,n2), iqn(1,n3), ms, nc, iszt, c(lc1), intpfg)
                    if (nc > 0) then
                        lc2 = lc2 + nc
                        if (lc2 >= l) go to 899
                        ic = 1
                        do lc = lc1, lc2
                            iqns(1:2,1:nam,lc) = iqns(1:2,1:nam,l)
                            iqns(1,n1,lc) = iszt(1,ic)
                            iqns(1,n2,lc) = iszt(2,ic)
                            ic = ic + 1
                            c(lc) = c(lc) * c(l)
                        end do
                        if (iqns(2,n3,l) < 0) then
                            if (iqn(3,n1) >= 0 .and. iqn(3,n2) >= 0) cycle
                            c(lc1:lc2) = -c(lc1:lc2)
                        else if (iqns(2,n3,l) == 0) then
                            if (iqn(2,n1) /= 0) then
                                if (lc2 + nc >= l) go to 899
                                lc3 = lc2
                                sign = one
                                if (iqn(3,n3) < 0) sign = -one
                                do lc = lc1, lc3
                                    lc2 = lc2 + 1
                                    iqns(1:2,1:nam,lc2) = iqns(1:2,1:nam,lc)
                                    iqns(2,n1,lc2) = -iqns(2,n1,lc2)
                                    iqns(2,n2,lc2) = -iqns(2,n2,lc2)
                                    c(lc) = c(lc) * root2
                                    c(lc2) = c(lc) * sign
                                end do
                                cycle
                            end if
                            if (iqn(3,n1) * iqn(3,n2) /= iqn(3,n3)) go to 511
                        end if
                        cycle
                    end if
                    go to 511
                end do
                go to 100
            end if
            n = n + 1
            if (n < niam) go to 300
            go to 299
        end if
        n3 = n3 - 1
        if (n3 > 0) go to 110
        return

    299 niam1 = niam - 1
        write(nftw,'("0ERROR IN COUPLING TREE"//(3I5))') ((icup(i,j), i=1,3), j=1,niam1)
        lc2 = 0
        return

    511 write(nftw,'("0COUPLING IMPOSSIBLE  MULTIPLICITIES FOLLOW"//(3I5))') (iqn(i,n1), iqn(i,n2), iqn(i,n3), i=1,3)
        lc2 = 0
        return

    899 write(nftw,'("0STORAGE OVERFLOW IN VECTOR COUPLING")')
        lc2 = 0

    end subroutine wfcple


    !> \brief TODO...
    !>
    !> One of interesting features of this subroutine is that it can detect if the storage for determinants is about to
    !> overflow. In such a case it will automatically flush the data to the output file (not waiting for the call to
    !> \ref congen_driver::csfout), providing space for further calculation.
    !>
    !> \param nncsf   Number of CSFs (total).
    !> \param nadel   Number of processed configurations.
    !> \param iidist  ?
    !> \param iidis3  ?
    !> \param nodi    Number of determinants per CSF (1 ... \c noi).
    !> \param ndi     Compressed determinants (1 ... \c nid).
    !> \param cdi     Determinant coefficients (1 ... \c nid).
    !> \param ndel    Configurations read from input.
    !> \param pqnshl  ?
    !> \param x       Real workspace, far enough after the chunk corresponding to \c pqnshl.
    !> \param ix      Actually, the same address as \c x, so the same workspace (just cast to integer).
    !> \param nx      Available elements in the real workspace.
    !>
    subroutine wfn (nncsf, nadel, iidist, iidis3, nodi, ndi, cdi, ndel, pqnshl, x, ix, nx)

        ! noi # of states (nodi)
        ! nid # of determinants (cdi)
        ! ni  # of replacements (ndi)

        use precisn,     only : wp
        use congen_data, only : lratio, iidis2, lcdi, lndi, ni, nid, noi, exdet, exref, norep, nonew, ntso, &
                                noimx => nodimx, nidmx => cdimx, jmx => ndimx, megul, nsoi, ncall, nftw,    &
                                confpf, ncsf, ndist, nshl, occshl, sshl => sshlst, shlmx1, nndel, nelec => nnlecg, ndimx

        integer :: iidis3, iidist, nadel, nncsf, nx
        real(kind=wp), dimension(*) :: cdi, x
        integer, dimension(*) :: ix, ndel, ndi, nodi, pqnshl
        intent (in) ndel, pqnshl
        intent (out) iidis3
        intent (inout) cdi, iidist, nadel, ndi, nncsf, nodi

        integer :: det, i, ia, ib, id, idist, iidis1, irep, j, k, k1, kshl, kshlst, lf, li, nd, nii

        save iidis1 ! JMC adding this in consultation with Jonathan Tennyson because the variable was found to be used but not set in tests...
                    ! The other variables initialized in the 1st call to this routine are (former) common block variables now in congen_data.

        ! reset counters for a new wave function group ("ncall = 1" is set in CSFGEN before calling DISTRB)
        if (ncall == 1) then
            lcdi   = 0
            lndi   = 0
            noi    = 0
            nid    = 0
            iidis1 = 0
            iidis2 = 0
            ni     = 0
            ncall  = 0
        end if

        ! here most work is done
        call state (nshl, x, nx, nd, confpf)

        if (nd <= 0) then
            write(nftw,'("1","*******   ERROR IN WFN ND=0"//)')
            stop 70
        end if

        ! no distributions generated in ASSIGN
        if (ndist == 0) return

        if (nndel /= 0 .and. nadel > nndel) return

        iidist = 0
        kshlst = 0

        ! for all distributions generated by ASSIGN
        do idist = 1, ndist

            nncsf = nncsf + 1

            if (nndel /= 0 .and. ndel(nadel) /= nncsf) then
                kshlst = kshlst + nshl
                cycle
            end if

            ncsf = ncsf + 1
            if (nndel /= 0) then
                nadel = nadel + 1
                iidist = iidist + 1
            end if

            k  = nd * lratio
            ia = nd * (2 * nelec + 1) + ni
            ib = nd + nid
            k1 = noi + 1

            if (ia > jmx .or. ib > nidmx .or. k1 > noimx) then
                if (nndel == 0 .or. iidis2 /= 0) then
                    write(megul) noi, (nodi(i), i=1,noi)
                    write(megul) nid, (cdi(i),  i=1,nid)
                    write(megul) ni,  (ndi(i),  i=1,ni)
                    lcdi = lcdi + nid
                    lndi = lndi + ni
                end if
                noi = 0
                ni = 0
                nid = 0
            end if

            nii = ni + 1
            do id = 1, nd
                exdet(1:ntso) = 0
                kshl = kshlst
                lf = 0
                do i = 1, nshl
                    li = lf + 1
                    lf = lf + occshl(i)
                    kshl = kshl + 1
                    do j = li, lf
                        k = k + 1
                        det = ix(k) + nsoi(sshl(i)) + (pqnshl(kshl) - 1) * shlmx1(sshl(i))
                        exdet(det) = 1
                    end do
                end do
                ia = 1
                ib = 1
                do i = 1, ntso
                    if (exdet(i) < exref(i)) then
                        norep(ia) = i
                        ia = ia + 1
                    else if (exdet(i) > exref(i)) then
                        nonew(ib) = i
                        ib = ib + 1
                    end if
                end do
                irep = ib - 1
                ni = ni + 1

                ! ZM check that NI does not overflow NDIMX
                !    If NDI was allocatable we could resize it here
                if (ni + 2 * irep > ndimx) then
                    write(nftw,'("*******   ERROR IN WFN, NDIMX TOO SMALL",2I10)') ni + 2 * irep, ndimx
                    stop 81
                end if

                ndi(ni) = irep
                do i = 1, irep
                    ni = ni + 1
                    ndi(ni) = norep(i)
                    ndi(ni + irep) = nonew(i)
                end do
                ni = ni + irep

            end do

            if (confpf >= 40) call print5 (nd, ndi(nii))
            noi = noi + 1

            ! ZM check that NOI does not overflow NODIMX
            !    If NODI was allocatable we could resize it here
            if (noi > noimx) then
                write(nftw,'("*******   ERROR IN WFN, NODIMX TOO SMALL",2I10)') noi, noimx
                stop 82
            end if

            nodi(noi) = nd

            ! ZM check that NID+ND does not overflow CDIMX
            !    If CDI was allocatable we could resize it here
            if (nid + nd > nidmx) then
                write(nftw,'("*******   ERROR IN WFN, CDIMX TOO SMALL",2I10)') nid + nd, nidmx
                stop 83
            end if

            do i = 1, nd
                nid = nid + 1
                cdi(nid) = x(i)
            end do

            kshlst = kshl

        end do

        iidis1 = iidis1 + iidist
        iidis2 = iidist + iidis2
        if (confpf >= 30 .and. iidist > 0) call print4 (nd, x, ix(nd*lratio+1))
        iidis3 = iidis1

    end subroutine wfn

end module congen_distribution
