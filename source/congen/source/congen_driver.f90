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
!> \brief Main CONGEN subroutines.
!>
module congen_driver

    implicit none

    ! module entry point
    public csfgen

    ! internal support routines
    private csfout
    private getcon
    private getcup
    private getref
    private icgcf
    private stwrit
    private subdel
    private wfnin
    private wfnin0
    private wrnmlt
    private wvwrit

contains

    !> \brief Central CONGEN subroutine
    !>
    !> Main driver routine: this is a transportable (Cray/IBM) version of CONGEN with the D2h symmetry for MOLECULE re-instated.
    !>
    subroutine csfgen

        use precisn,     only : wp
        use consts,      only : zero => xzero, two => xtwo
        use congen_data, only : mxtarg, nx => nu, ny, nz, jx, jy, jz, ky, kz, ntcon, test, nrcon, nobt,             &
                                nsym, nobi, ift, next, navail => nx, icdi, iexcon, indi, inodi, iqnstr,             &
                                irfcon, ndel, ndel1, ndel2, exref, ntso, nodimx, cdimx, ndimx, megul, nsoi, nncsf,  &
                                ncall, ncsf, confpf, occst => occshl, pqnst, mshlst => mshl, gushst => gushl,       &
                                cupst => cup, sshlst, shlmx1, nqntot => qntot, nqntar => qntar, ngutot => gutot,    &
                                nnndel => nndel, nnlecg, nisz, nsmtyp => symtyp, nftw, nftr, exdet, exref

        use iso_c_binding,       only : c_loc, c_f_pointer
        use congen_distribution, only : distrb
        use congen_pagecontrol,  only : addl, newpg, space
        use congen_projection,   only : projec

        real(wp), allocatable, dimension(:), target :: bmx
        real(wp), pointer :: bmx_ptr
        integer,  dimension(:), pointer :: int_ptr, ndel_ptr, ndel1_ptr, ndel2_ptr, irfcon_ptr
        integer :: byproj, cdimn, conmes, d, defltc, guold, gutot, i, idiag, iposit, is, iscat, iso, isz, itu, j, lc, lcdi,     &
               lcdo, lcdt, ln, lndi, lndo, lndt, lpp, lsquare, ltri, maxtgsym, megu, mflag, mold, mt, n, nadel, nbmx,           &
               nconmx, ncsf0, ncupp, ndimn, ndpmax, ndpp, ndprod, negmax, negr, nelecg, nelect, nelp, nelr, nemax, nerfg,       &
               nerfs, nextk, nfto, nndel, nobep, nodimn, npcupf, npmult, nrefo, nrefog, nrefop, nrerun, nrfgmx, nrfgoe,         &
               nrfoe, nrfomx, ns, nshgmx, nshlp0, nshlpt, nspf, nspi, nsymmx, nsymp, ntconp, err, gflag, NOEX
        integer, dimension(3,jx) :: cup, pqn
        logical :: ene, ener, enob, enrefo, epsno, eqnt, erefd, error, errorg, espace, esymt, qmoln
        logical, dimension(11) :: erfg
        logical, dimension(9) :: erfs
        character(len=80) :: gname, sname
        character(len=256) :: errmsg
        integer, dimension(mxtarg) :: gucont, mcont, mdegen, mrkorb, nctgt, notgt
        integer, dimension(jx) :: gushl, ksss, loopf, mshl, occshl, shlmx, sshl
        integer, dimension(jy) :: kdssv, nelecp, nshlp, nslsv
        integer, dimension(nx) :: nob, nob0, nob0l, nobl
        integer, dimension(nx) :: nobe, nobp, nobv
        integer, dimension(6) :: npflg
        integer, dimension(jz) :: nshcon
        integer :: ntgsmx, ntgsym, nwfngp, pqn2, symtyp
        real(kind=wp) :: pin, r, s, sz, thres
        integer, dimension(3) :: qntar, qntar1, qntot
        integer, dimension(nz) :: refdet        ! Reference determinant (assignment of electrons to spin-orbitals).
        integer, dimension(kz) :: refdtg        ! Wfngrp reference determinant (assignment of electrons to spin-orbitals).
        integer, dimension(ny) :: refgu         ! Gerade/ungerade quantum numbers for the reference determinant (per quintet).
        integer, dimension(ky) :: refgug        ! Gerade/ungerade quantum numbers for the wfngrp reference determinant (per quintet).
        integer, dimension(5,ny) :: reforb      ! Reference determinant setup (quintets from &state).
        integer, dimension(5,ky) :: reforg      ! Movable electrons setup (quintets from &wfngrp).
        integer, dimension(3,jx,jz) :: tcon     ! Constraints. ???

        equivalence(reforg(1,1), reforb(1,1))
        equivalence(refgu(1), refgug(1))

        ! named error flags
        equivalence(erfs(1), esymt)
        equivalence(erfs(2), enob)
        equivalence(erfs(3), epsno)
        equivalence(erfs(4), ene)
        equivalence(erfs(5), enrefo)
        equivalence(erfs(6), ener)
        equivalence(erfs(7), erefd)
        equivalence(erfs(8), eqnt)
        equivalence(erfs(9), espace)

        ! input namelists
        namelist /state / megul, nrerun, lcdt, lndt, nfto, ltri, idiag, thres, megu, npflg, nodimx, ndimx, cdimx,   &
                          byproj, lcdo, lndo, nftw, iscat, ntgsym, sname, lpp, confpf, symtyp, qntot, gutot, isz,   &
                          npmult, nob, reforb, refgu, nrefo, nelect, nndel, qmoln, gflag,                           &
                          iposit, nob0, nbmx, nobe, nobp, nobv, NOEX  ! this last row: positron control data

        namelist /wfngrp/ nelecg, ndprod, nelecp, nshlp, gname, reforg, refgug, nrefog, mshl, gushl, pqn, cup, defltc, &
                          npcupf, test, nrcon, nshcon, tcon, ntcon, qntar, lsquare

        !
        !     Hardcoded limits
        !

        conmes = 400    ! Quantenmol-specific output file unit with time estimate.
        nrfoe  = 30     ! initial value for "nrefop" (see below).
        nrfgoe = 10     ! ... ?? ... Used in GETCON.
        nerfs  = 9      ! Number of error condition flags in erfs.
        nerfg  = 11     ! Number of error condition flags in erfg.
        nsymmx = nx     ! Maximal number of symmetries.
        nrfomx = ny     ! Maximal number of reference orbital quintets.
        nemax  = nz     ! Maximal number of electrons.
        nrfgmx = ky     ! Maximal number of quintets used in single wave function group.
        negmax = kz     ! Maximal number of movable electrons.
        nshgmx = jx     ! Maximal number of pqn triplets per wave function group.
        ndpmax = jy     ! Maximal number of sets of orbitals per wave function group.
        nconmx = jz     ! Maximal number of constraints.

        nobl(:)   = 0
        nob0l(:)  = 0
        mdegen(:) = 0

        !
        !     Default input data
        !

        qmoln  = .false.        ! Indicates CONGEN running from Quantenmol.
        gflag  = 1              ! Flag to disable CSF grouping based on symmetry 0=OFF 1=ON (default)
        iscat  = 0              ! Output format: 0 = SPEEDY, 1 = SCATCI, 2 = SCATCI + information for phase correcting target CI wavefunctions.
        ntgsym = mxtarg         ! Number of target state symmetries
        nndel  = 0              ! Whether to read CSFs from input stream.
        megul  = 13             ! Congen output file with generated configurations.
        nrerun = 0              ! Restart flag for SPEEDY (not used in SCATCI).
        lcdt   = 500            ! Not used.
        lndt   = 5000           ! Not used.
        lcdo   = 500            ! Limit on number of determinants.
        lndo   = 5000           ! Limit on number of integers used in determinant description.
        nfto   = 15             ! Fortran data set number for output unit of SPEEDY.
        ltri   = 300            ! Not used.
        idiag  = -1             ! Matrix element evaluation flag passed to SPEEDY or SCATCI.
        nbmx   = 2000000        ! Workspace: Number of real numbers.
        byproj = 1              ! Project wave functions in 1 = CONGEN, 0 = SPEEDY, -1 = nowhere (should be avoided).
        megu   = 14             ! Output data set with the namelist input to be passed to SPEEDY.
        thres  = 1.e-10_wp      ! Threshold for storing coefficients of matrix elements that contribute to Hij.

        nobe(:) = 0             ! Number of electronic orbitals.
        nobp(:) = 0             ! Number of positronic orbitals.
        nobv(:) = 0             ! Not used, set to number of target (non-continuum) orbitals.

        npflg(1:6) = 0          ! Print flags for SPEEDY / SCATCI.

        cdimx  = 400            ! Workspace size for determinant (real) multiplication factors.
        ndimx  = 4000           ! Workspace size for packed determinant (integer) descriptions.
        nodimx = 100            ! Workspace size for number of packed determinants per CSF.
        cdimn  = 100            ! Default value for user namelist override of cdimx.
        ndimn  = 1000           ! Default value for user namelist override of ndimx.
        nodimn = 25             ! Default value for user namelist override of nodimx.

        isz    = 0              ! 2*Sz + 1
        lpp    = 0              ! Length of line printer page. Used to determine where to put page banner.
        confpf = 1              ! Print flag for the amount of print given of configurations and prototypes.
        symtyp = -1             ! Symmetry type, 0 = C_infv, 1 = D_infh, 2 = {D_2h, C_2v, C_s, E}.
        gutot  = 0              ! Inversion symmetry, used for D_infh only (+1 = gerade, -1 = ungerade).
        npmult = 0              ! Print flag for the D2h multiplication table (0 = no, 1 = yes).
        nelect = 0              ! Number of electrons.
        iposit = 0              ! Exotic particle flag: 0 = electrons only, +/-1 = positive non-electron, +/-2 = negative non-electron.
        noex = 0                ! Flag to switch off exchange between projectile electron and target electrons: 0 = no exchange

        qntot(1:3) = -2         ! Total quantum numbers (spin multiplicity, irreducible representation, inversion symmetry).

        nob(1:nsymmx)  = 0      ! Number of orbitals for configuration generation per symmetry, including continuum orbitals.
        nob0(1:nsymmx) = 0      ! As above, excluding continuum orbitals.
        nsoi(1:nsymmx) = 0      ! For each symmetry, index of the first spin-orbital (index runs across all symmetries).
        nobi(1:nsymmx) = 0      ! For each symmetry, index of the first orbital (index runs across all symmetries).

        nrefo = 0               ! How many quintets are needed to define the reference determinant.

        refdet(1:nemax)      = 0        ! Reference determinant (assignment of electrons to spin-orbitals).
        refgu(1:nrfomx)      = -2       ! Gerade (= +1) or underage (= -1) symmetry for the M value of each reference quintet.
        reforb(1:5,1:nrfomx) = -2       ! Reference determinant setup (quintets).
        erfs(1:nerfs)        = .false.  ! Error condition indicators.

        nsym   = nsymmx         ! Number of symmetries: initialized to maximal number of symmetries.
        nsymp  = nsymmx         ! Number of symmetries, whose info to print in stwrit.
        nrefop = nrfoe          !
        nelp   = 0              ! Auxiliary variable typically holding the (calculated) number of electrons.

        ! read the input namelist
        read(nftr, state, iostat = i, iomsg = errmsg)

        ! abort if namelist read failed
        if (i /= 0) then
            write(nftw,'("1*****   FAILED TO READ NAMELIST &STATE:",/,9X,A)') trim(errmsg)
            return
        end if

        ! choose default value for idiag if not set by input
        if (idiag < 0) then
            if (iscat <= 0) idiag = 0   ! SPEEDY format
            if (iscat >  0) idiag = 1   ! SCATCI format
        end if

        ! Check the NOB-values
        !
        !  NOBE(i) number of electronic orbitals, default nobe(i)=nob(i)
        !  NOBP(i) number of positronic orbitals, default nobp(i)=0
        !  NOBV(i) at the moment not used, but printed on output, default nobv(i)=nob0(i)

        do i = 1, nsym
            if (nobp(i) == 0) then
                nobe(i) = nob(i)
            end if
            nobep = nobe(i) + nobp(i)
            if (nobep /= nob(i)) then
                write(nftw,*) 'ERROR on input:'
                write(nftw,*) 'not: NOB(i)=NOBE(i)+NOBP(i)'
                write(nftw,*) 'i=', i
                write(nftw,*) 'NOB(i)=', nob(i)
                write(nftw,*) 'NOBE(i)=', nobe(i)
                write(nftw,*) 'NOBP(i)=', nobp(i)
            end if
            if (nobv(i) == 0) then
                nobv(i) = nob0(i)
            end if
            if (nob0(i) > nob(i)) then
                write(nftw,*) 'ERROR on input:'
                write(nftw,*) 'NOB0(i) > NOB(i)'
                write(nftw,*) 'i=', i
                write(nftw,*) 'NOB(i)=', nob(i)
                write(nftw,*) 'NOB0(i)=', nob0(i)
            end if
        end do
 
        ! get index (-> nsym) of last symmetry with at least one orbital
        nsym = 0
        do i = 1, nsymmx
            if (nob(i) /= 0) nsym = i
        end do
        nsymp = nsym

        ! set limit on shell occupancy
        !   symtyp = 0 (C_infv) : 2,4,4,4,...
        !   symtyp = 1 (D_infh) : 2,2,4,4,...
        !   symtyp other        : 2,2,2,2,...
        shlmx1(1:nsymmx) = 2
        if (symtyp == 0) shlmx1(2:nsymmx) = 4
        if (symtyp == 1) shlmx1(3:nsymmx) = 4

        ! prepare global indices of orbitals and spin-orbitals across symmetries
        nobt = 0
        ntso = 0
        do i = 1, nsym
            nsoi(i) = ntso + 1                  ! 1-based index of spin-orbital across symmetries
            nobi(i) = nobt                      ! 0-based index of orbital across symmetries
            nob(i)  = abs(nob(i))
            nobt = nobt + nob(i)
            ntso = ntso + shlmx1(i) * nob(i)
        end do

        allocate(exdet(ntso), exref(ntso), stat = err)

        ! set up error flags
        enob  = (nsym == 0)
        epsno = (err /= 0)
        ene   = (nelect <= 0 .or. nelect > nemax)
        esymt = (symtyp < 0)

        ! continue only if no error condition is raised
        if (.not. any(erfs)) then

            ! assemble reference list of spin-orbitals (target "reference determinant")
            call getref (reforb, refgu, nrefo, nelect, refdet, nelr, nsoi, nob, shlmx1, nsym, symtyp, nrfomx, enrefo, ener, erefd)

            ! more initializations
            nelp = nelect
            if (.not. enrefo) nrefop = nrefo
            confpf = max(confpf, 1)
            if (isz == 0) isz = qntot(1)

            ! check validity of total quantum numbers
            eqnt = (qntot(1) <= 0 .or. qntot(2) < 0)                                &
                   .or. (qntot(2) == 0 .and. symtyp <= 1 .and. abs(qntot(3)) /= 1)  &
                   .or. (abs(isz - 1) > qntot(1) - 1)                               &
                   .or. (symtyp == 1 .and. abs(gutot) /= 1)

            ! update array bounds based on namelist input
            cdimx  = max(cdimn,  cdimx)
            ndimx  = max(ndimn,  ndimx)
            nodimx = max(nodimn, nodimx)

            ! allocate the wp-real (!) workspace
            allocate(bmx(nbmx))

            ! calculate positions of data in the workspace
            icdi   = 1
            indi   = icdi  + cdimx
            inodi  = indi  + ndimx
            ndel   = inodi + nodimx
            ndel1  = ndel  + nndel
            ndel2  = ndel1 + nndel
            next   = ndel2 + nndel
            navail = nbmx  - next + 1
            nextk  = next  / 1024

            ! Workspace is composed of the following consecutive chunks of real numbers
            !    icdi  : icdi  + cdimx   ... determinant coefficients (real number per determinant)
            !    indi  : indi  + ndimx   ... packed determinants (size + replacements per determinant)
            !    inodi : inodi + nodimx  ... states ()
            !    ndel  : ndel  + nndel   ... space for configurations read from input (unit nftr)
            !    ndel1 : ndel1 + nndel   ... workspace for reading configurations from input (unit nftr)
            !    ndel2 : ndel2 + nndel   ... workspace for reading configurations from input (unit nftr)
            !
            ! + further chunks per every wave function group:
            !    irfcon : irfcon + nobt * ntcon ... orbital constraints
            !    iexcon : iexcon + nobt         ... more orbital contrains?
            !    iqnstr : ?                     ... only temporarily, overwritten by next wave function group

            write(nftw,'(" NBMX =",I9," WORDS")') nbmx
            write(nftw,'(" ***** REGION USED FOR INPUT DATA ",I7," WORDS ",I5," K")') next, nextk
            write(nftw,'(" ***** LEFT ",I9,"  WORDS")') navail

            ! are we out of space already?
            espace = (navail <= 0)

        end if

        ! print current setup information
        call stwrit (nelect, confpf, qntot, cdimx, icdi, ntso, symtyp, ndimx, indi, nrefo, nodimx, inodi, nsym, gutot, nbmx,    &
                     isz, navail, idiag, megu, thres, lcdt, megul, lndt, nfto, nrerun, ltri, npflg, nndel, nob, nsoi, nsymp,    &
                     refdet, nerfs, erfs, nrefop, reforb, refgu, nelp, lpp, sname, error, byproj, lndo, lcdo, iposit, nob0,     &
                     npmult, ntgsym, mxtarg, nobe, nobp, nobv)

        ! if requested, read CSFs directly from input stream
        if (nndel /= 0) then
            bmx_ptr => bmx(ndel)  ; call c_f_pointer (c_loc(bmx_ptr),  ndel_ptr, (/1/))
            bmx_ptr => bmx(ndel1) ; call c_f_pointer (c_loc(bmx_ptr), ndel1_ptr, (/1/))
            bmx_ptr => bmx(ndel2) ; call c_f_pointer (c_loc(bmx_ptr), ndel2_ptr, (/1/))
            call subdel (ndel_ptr, ndel1_ptr, ndel2_ptr, nndel)
        end if

        ! general default for wfngrp arrays
        call wfnin (nwfngp, nadel, nncsf, ncsf, lcdi, lndi, nelecg, ndprod, nrefog, npcupf, negmax, refdtg, nrfgmx, refgug, &
                    ntcon, reforg, nshgmx, nsymmx, mshl, gushl, pqn, cup, ndpmax, nshlp, nconmx, test, nrcon, nshcon, tcon, errorg)

        ncsf0     = 0           ! number of CSFs generated up to the previous wfn group
        ntgsmx    = ntgsym      ! 
        ntgsym    = 0           ! number of target symmetries (assuming that all wfn groups are "target × continuum")
        maxtgsym  = 0           ! number of target symmetries (actually taking into account the value of lsquare)
        mold      = -1          ! symmetry (M-value) of the last shell (pqn triplet) in the previous wfn group
        guold     = 0           ! g/u symmetry of the last shell in the previous wfn group
        qntar1(1) = -1          ! target quantum numbers in the previous wfn group
        pqn2      = 0           ! orbital sequence start of the last shell in the previous wfn group
        lsquare   = 0           ! do not assume L2 configuration by default

        ! loop over wfngrp input sets
        wfngrp_loop: do while (nndel == 0 .or. nadel < nndel)

            ! set default values to some parameters (others are carried over from the previous namelist, if any)
            call wfnin0 (nelecp, defltc, nerfg, erfg, gname, qntar, errorg, ndpmax)

            ! read the user settings from the next &wfngrp namelist
            read(nftr, wfngrp, iostat = i)

            ! check success
            if (i /= 0) then
                call space (2)
                call addl (1)
                if (nwfngp == 0) then
                    write(nftw,'(" *****  NO WFNGRP INPUT FOUND")')
                    return
                end if
                write(nftw,'(" *****  END OF FILE ON INPUT")')
                exit wfngrp_loop
            end if

            call wfngrp_sanity_check(nftw, nob, ndprod, nelecp, nshlp, pqn, mshl)

            ! total number of pqn triplets for this &wfngrp
            nshlpt = sum(abs(nshlp(1:ndprod)))

            ! assemble information for phase correction data
            if (iscat < 2 .or. mold < -1) then

                ! no phase correction data required (or all that was needed is already done, do nothing)

            else if (mold == mshl(nshlpt) .and. (symtyp /= 1 .or. guold == gushl(nshlpt)) .and. &
                     all(qntar(1:3) == qntar1(1:3)) .and. pqn2 == pqn(2,nshlpt) .and. lsquare /= 1 .and. gflag /= 0) then

                ! the same target state? - check continuum orbitals consistent!
                if (notgt(ntgsym) /= pqn(3,nshlpt) - pqn(2,nshlpt) + 1) then
                    write(nftw,'(//," Attempt to perform CI target calculation with different length continua for same target:")')
                    write(nftw,'(   " Number of continua, last WFNGRP",I4)') notgt(ntgsym)
                    write(nftw,'(   " Number of continua, this WFNGRP",I4)') pqn(3,nshlpt) - pqn(2,nshlpt) + 1
                    write(nftw,'(   " STOP")')
                    stop
                end if

            else

                ! new target state detected: first save data about old one
                if (ntgsym >= 1) then
                    nctgt(ntgsym) = (ncsf - ncsf0) / notgt(ntgsym)
                    mrkorb(ntgsym) = nspi
                end if

                if (qntar(1) <= 0 ) then

                    ! all target states are finished
                    mold = -2

                else

                    ncsf0 = ncsf
                    ntgsym = ntgsym + 1
                    if (ntgsym > ntgsmx) then
                        write(nftw,'(" Number of target states [",I5,"] has exceeded mxtarg [",I5,"].")') ntgsym, ntgsmx
                        write(nftw,'(" Increase mxtarg and recompile.")')
                        stop 1
                    end if
                    if (lsquare == 0) maxtgsym = maxtgsym + 1
                    write(nftw,'(/," Target state number",I3) ') ntgsym
                    write(nftw,'(" TARGET MULTIPLICITY =",I5,5X,"TARGET SYMMETRY =",I5,5X,"TARGET INVERSION SYMMETRY =",I5)') qntar
                    write(nftw,'(" Coupling to continuum with M =",I3)') mshl(nshlpt)
                    if (symtyp == 1) write(nftw,'(27x,"GU =",I3)') gushl(nshlpt)

                    ! save continuum electron data for future use
                    notgt(ntgsym) = pqn(3,nshlpt) - pqn(2,nshlpt) + 1
                    mcont(ntgsym) = mshl(nshlpt)
                    if (symtyp == 1) gucont(ntgsym) = gushl(nshlpt)
                    pqn2 = pqn(2,nshlpt)

                    ! for degenerate symmetries/degenerate target states, need extra
                    ! information to sort out coupling of possible second continuum
                    if (symtyp <= 1 .and. qntot(2) > 0 .and. qntar(2) > 0) then
                        mdegen(ntgsym) = max(qntot(2), qntar(2)) - mcont(ntgsym)
                    end if

                    qntar1(1:3) = qntar(1:3)
                    mold = mshl(nshlpt)
                    if (symtyp == 1) guold = gushl(nshlpt)

                end if

            end if

            ndpp    = ndpmax
            nelp    = 0
            ncupp   = 0
            ntconp  = 0
            nwfngp  = nwfngp + 1
            erfg(1) = .true.
            if (nelecg <= 0 .or. nelecg > negmax) go to 150
            erfg(1) = .false.
            ndprod  = max(1, ndprod)
            erfg(2) = .true.
            if (ndprod > ndpmax) go to 150
            erfg(2) = .false.
            ndpp    = ndprod
            nshlpt  = 0

            nelecp(1:ndprod) = abs(nelecp(1:ndprod))
            nshlp(1:ndprod)  = abs(nshlp(1:ndprod))
            negr   = sum(nelecp(1:ndprod))
            nshlpt = sum(nshlp(1:ndprod))

            erfg(3) = (any(nshlp(1:ndprod) == 0))

            nshlp0 = nshlpt - nshlp(ndprod)
            if (nshlpt > nshgmx) erfg(3) = .true.
            if (negr /= nelecg) erfg(4) = .true.
            if (erfg(3) .or. erfg(4)) go to 150

            ! reference determinant (refdtg) for the wave function group
            call getref (reforg, refgug, nrefog, nelecg, refdtg, nelr, nsoi, nob, shlmx1, &
                        nsym, symtyp, nshgmx, erfg(5), erfg(6), erfg(7))

            if (.not. erfg(5)) nrefop = nrefog
            if (.not. erfg(6)) nelp   = nelecg
            if (erfg(5) .or. erfg(6) .or. erfg(7)) go to 150

            ! Set "exref" to contain only the non-movable subset of the global reference determinant.
            exref(1:ntso) = 0                                       ! start with unpopulated spin-orbitals
            exref(refdet(1:nelect)) = 1                             ! populate spin-orbitals according to the reference determinant
            if (any(exref(refdtg(1:nelecg)) == 0)) erfg(7) = .true. ! non-movable spin-orbitals must be a subset of all reference spin-orbitals
            exref(refdtg(1:nelecg)) = 0                             ! un-populate spin-orbitals with movable electrons
            if (erfg(7)) go to 150

            ! check shell data...
            erfg(8) = .true.
            if (NOEX .eq. 0) then
              if (iposit /= 0 .and. nelecp(ndprod) /= 1) go to 150  ! ... but skip the check if for non-electrons ...?
            endif

            do i = 1, nshlpt                                        ! scan all sets of orbitals
                mshl(i) = abs(mshl(i))                              ! use absolute value of the M-value (symmetry)
                mt = mshl(i) + 1                                    ! change M-value from 0-based index to 1-based index
                if (symtyp == 1) then                               ! special section for degenerate point groups
                    if (abs(gushl(i)) /= 1) go to 150               ! gerage (+1) or ungerade (-1) only possible
                    itu = 1                                         ! let's call this a kind of parity,...
                    if (mod(mt,2) == 0) itu = -1                    ! ...as it changes its sign across the symmetry eigenspaces
                    mt = mt + mt - abs((gushl(i) + itu) / 2)        ! expanded symmetry index combining the orig symmetry and g/u
                end if
                sshl(i) = mt                                        ! 1-based index of symmetry (expanded, in case of degenerate groups)
                if (mt > nsym) go to 150                            ! out of bounds symmetry index!
                if (pqn(1,i) /= 0) then                             ! transform pqn=x00 to xxx for unified processing later
                    pqn(2,i) = pqn(1,i)
                    pqn(3,i) = pqn(2,i)
                end if
                d = pqn(3,i) - pqn(2,i)                             ! distance between orbitals (= orbital count - 1)
                if (d < 0) go to 150                                ! wrong order of orbitals in triplet!
                shlmx(i) = shlmx1(mt) * (d + 1)                     ! maximal occupation of i-th set of orbitals
                nspi = nsoi(mt) + (pqn(2,i) - 1) * shlmx1(mt)       ! index of the first spin-orbital in the current set orbitals
                nspf = nspi + shlmx(i) - 1                          ! index of the last spin-orbital in the current set of orbitals
                if (iposit /= 0 .and. i >= nshlp0) cycle            ! ... some extra for non-electrons ...?
                if (any(exref(nspi:nspf) /= 0)) go to 150           ! the set of orbitals mustn't contain non-movable ones
                if (pqn(3,i) > nob(mt)) then                        ! out of bounds orbital index!
                    write(nftw,'(//," Error: PQN number",i3," accesses orbital number",i3)') i, pqn(3,i)
                    write(nftw,'(   " symmetry",i2," only contains",i3," orbitals")') mt, nob(mt)
                    go to 150
                end if
            end do
            erfg(8) = .false.

            ! ???
            ncupp = nshlpt - 1
            call getcup (nshlpt, defltc, ndprod, nshlp, cup(1,1), erfg(9))
            if (erfg(9)) npcupf = 1
            ntcon = min(abs(ntcon), nconmx)

            ! calculate positions of new data in the workspace
            irfcon = next
            if (ntcon == 0) go to 150
            iexcon = irfcon + nobt * ntcon
            next = iexcon + nobt
            navail = nbmx - next + 1
            nextk = next / 1024
            write(nftw,'(" ***** REGION USED FOR DETERMINANTS ",I7," WORDS ",I5," K")') next, nextk
            write(nftw,'(" **** LEFT ",I7,"  WORDS")') navail
            if (navail <= 0) erfg(11) = .true.

            ! set up constraints
            bmx_ptr => bmx(irfcon) ; call c_f_pointer (c_loc(bmx_ptr), irfcon_ptr, (/1/))
            call getcon (ntcon, nshcon, nrcon, nshgmx, nrfgoe, ntconp, nelecg, nsym,  &
                         nobt, nob, nobi, nsoi, shlmx1, exref, tcon, irfcon_ptr,      &
                         erfg(10))
        
            ! position of "qnstr"(?) data in the workspace
        150 iqnstr = next

            ! write information about wave function group to output...
            bmx_ptr => bmx(irfcon) ; call c_f_pointer (c_loc(bmx_ptr), irfcon_ptr, (/1/))
            call wvwrit (nwfngp, gname, nelecg, defltc, irfcon, ndprod, symtyp, ntcon,    &
                         navail, nrefog, nelecp, nshlp, qntar, nshlpt, nshgmx, mshl,      &
                         gushl, pqn, cup, ncupp, npcupf, refdtg, nelp, ntconp, nshcon,    &
                         nobt, reforb, refgu, test, nrcon, tcon, irfcon_ptr, nerfg,       &
                         erfg, ndpp, nrefop, lsquare, errorg)

            ! ... and do nothing more if there are errors
            if (errorg) return

            ! push data to global arrays
            nqntot(1:3) = qntot(1:3)
            nqntar(1:3) = qntar(1:3)
            nisz   = isz
            ngutot = gutot
            nnlecg = nelecg
            nsmtyp = symtyp
            nnndel = nndel

            exref(1:ntso)           = 0     ! start with un-populated spin-orbitals
            exref(refdtg(1:nelecg)) = 1     ! populate spin-orbitals with movable electrons only

            ncall = 1                       ! this will make "wfn" routine initialize some arrays
            ift = 1                         ! this will make "assign" routine initialize some arrays

            call icgcf                      ! precompute needed binomial coefficients

            call distrb (nelecp, nshlp, shlmx, occshl, nslsv, kdssv, loopf, ksss, pqn,     &
                         occst, shlmx1, pqnst, mshl, mshlst, gushl, gushst, cup, cupst,    &
                         ndprod, symtyp, confpf, sshl, sshlst, ncsf, nncsf, nadel,         &
                         nndel, bmx, nftw, noex)

            ! write constructed (non-yet-projected) CSFs to output file
            bmx_ptr => bmx(1) ; call c_f_pointer (c_loc(bmx_ptr), int_ptr, (/1/))
            call csfout (lc, ln, megul, nndel, bmx, int_ptr)

            ! update total number of stored coefficients and integers forming packed determinants
            lcdi = lcdi + lc
            lndi = lndi + ln

        end do wfngrp_loop

        ! check that data for final target state has been saved
        mflag = 0
        if (ntgsym >= 1) then
            if (mold > -2) then
                nctgt(ntgsym) = (ncsf - ncsf0) / notgt(ntgsym)
                mrkorb(ntgsym) = nspi
            end if

            ! for degenerate symmetries/degenerate target states, need extra
            ! information to sort out coupling of possible second continuum
            ! check this for errors and missed couplings
            if (symtyp <= 1 .and. qntot(2) > 0) then
                if (mdegen(ntgsym) > 0) mdegen(ntgsym) = 0
                do n = 2, ntgsym
                    if (mdegen(n) > 0) then
                        if (mdegen(n-1) <= 0) then
                            write(nftw,'(/," WARNING: for target state number",I3)') ntgsym
                            write(nftw,'(  " Coupling to upper continuum only detected")')
                            write(nftw,'(  " Calculation may give target phase problems")')
                            mdegen(n) = 0
                        else if (nctgt(n) /= nctgt(n-1)) then
                            write(nftw,'(/," Target states",I3," and",I3)') n - 1, n
                            write(nftw,'(  " analysed for degenerate coupling to the continuum")')
                            write(nftw,'(  " But number of CSFs differ:",I6," and",I6," respectively: STOP")') nctgt(n-1), nctgt(n)
                            stop
                        end if
                    else if (mdegen(n) > 0) then
                        if (mdegen(n-1) > 0) mdegen(n-1) = 0
                    end if
                    if (mdegen(n) /= 0) mflag = max(mflag, nctgt(n))
                end do
            end if
        end if

        call space (2)
        call addl (1)

        write(nftw,'(" ********** TOTAL NUMBER OF CSF''S GENERATED IS ",I9)') ncsf

        !
        ! Quantenmol: open a file fort.400 to inform the user about estimated time of the run
        !

        if (qmoln .and. iscat < 2 .and. megul == 70) then

            open(unit = conmes, status = 'unknown')
            write(conmes,'(/," *** TOTAL NUMBER OF GENERATED CSF''S FOR THE GROUND STATE IS ",I6)') ncsf

            if (ncsf <= 600) then
                write(conmes,'("*** This target calculation should not take long !",/)')
            else if (ncsf > 600 .and. ncsf < 12000) then
                write(conmes,'("*** This target calculation will take a few hours !")')
                write(conmes,'("*** You can have a cup of tea and come back later.",/)')
            else if (ncsf > 12000 .and. ncsf <= 22000) then
                write(conmes,'("*** Oops, This target calculation is very big !")')
                write(conmes,'("*** You can come back tomorrow.",/)')
            else if (ncsf > 22000 .and. ncsf <= 80000) then
                write(conmes,'("*** Oops, This target calculation is very big ")')
                write(conmes,'("***    and can take several days to run !",/)')
            else if (ncsf > 80000) then
                write(conmes,'("*** Oops, This target calculation is too big ")')
                write(conmes,'("***    to be computationally possible !")')
                write(conmes,'("*** Rerun with a smaller basis or contact ")')
                write(conmes,'("*** technical support: support@quantemol.com")')
            end if

            close(conmes)

        end if

        s  = real(qntot(1) - 1,kind = wp) / two
        sz = real(isz - 1, kind = wp) / two
        r  = zero
        if (qntot(2) == 0 .and. symtyp <= 1) r = real(qntot(3), kind = wp)
        pin = zero

        ! save contents of nob and nob0 prior to any symmetry conversion
        nobl(1:nsym) = nob(1:nsym)
        nob0l(1:nsym) = nob0(1:nsym)
        if (symtyp == 1) then
            pin = real(gutot, kind = wp)
            j = 0
            do i = 1, nsym, 2
                j = j + 1
                nob(j) = nob(i) + nob(i+1)
                nob0(j) = nob0(i) + nob0(i+1)
            end do
            nsym = j
        end if

        if (allocated(bmx)) deallocate(bmx)

        call newpg

        ! stage 2: read CSFs back and project
        rewind megul

        ! SCATCI data output
        if (byproj > 0 .or. iscat > 0) then
            call projec (sname, megul, symtyp, qntot(2), s, sz, r, pin, ncsf, byproj, idiag, npflg, thres, nelect, nsym,    &
                         nob, refdet, nftw, iposit, nob0, nobl, nob0l, iscat, ntgsym, notgt, nctgt, mcont, gucont, mrkorb,  &
                         mdegen, mflag, nobe, nobp, nobv, maxtgsym)
        end if

        ! SPEEDY data output
        if (iscat <= 0) then
            call wrnmlt (megu, sname, nrerun, megul, symtyp, qntot(2), s, sz, r, pin, ncsf, byproj, lcdi, lndi, lcdo, lndo, &
                         lcdt, lndt, nfto, ltri, idiag, npflg, thres, nelect, nsym, nob, refdet, nftw, iposit, nob0, nobl,  &
                         nob0l, nx, nobe, nobp, nobv)
            rewind megu
        end if

        ! final text information output
        if (confpf >= 1) then
            call wrnmlt (nftw, sname, nrerun, megul, symtyp, qntot(2), s, sz, r, pin, ncsf, byproj, lcdi, lndi, lcdo, lndo, &
                         lcdt, lndt, nfto, ltri, idiag, npflg, thres, nelect, nsym, nob, refdet, nftw, iposit, nob0, nobl,  &
                         nob0l, nx, nobe, nobp, nobv)
        end if

    end subroutine csfgen


    !> \brief   Check electron distribution
    !> \authors J Benda
    !> \date    2024 - 2025
    !>
    !> Perform some basic checks on the electron distribution parameters to filter out trivial errors.
    !>
    !> \param nftw      Output unit.
    !> \param ndprod    Number of sets of shells.
    !> \param nelecp    Number of electrons in each set.
    !> \param nshlp     Number of shells (triplets) in each set.
    !> \param pqn       Orbitals in each shell.
    !> \param mshl      Symmetry of shells above.
    !>
    subroutine wfngrp_sanity_check (nftw, nob, ndprod, nelecp, nshlp, pqn, mshl)

        use congen_data, only: nu, shlmx1

        integer, intent(in) :: nftw, nob(:), ndprod, nelecp(:), nshlp(:), pqn(:, :), mshl(:)

        integer :: iset, jset, ishl, jshl, ipqn, jpqn, iorb, jorb, irr, npqn, nso

        character(len=*), parameter :: err_few_shells = &
            '(/,"ERROR: Set of shells ",i0," has to accommodate ",i0," particles but has only ",i0," spinorbitals!")'
        character(len=*), parameter :: err_missing_pqn = &
            '(/,"ERROR: PQN triplet ",i0," unspecified, but used by shell ",i0," in set ",i0,"!")'
        character(len=*), parameter :: err_invalid_mshl = &
            '(/,"ERROR: Invalid MSHL value ",i0," for PQN triplet ",i0,"!")'
        character(len=*), parameter :: err_overlapping_pqns = &
            '(/,"ERROR: Orbital ",i0," with MSHL ",i0," is used in two PQN triplets ",i0," and ",i0,"!")'
        character(len=*), parameter :: err_out_of_orbitals = &
            '(/,"ERROR: PQN triplet ",i0," references orbitals beyond NOB (",i0,")!")'

        ! check that there are enough PQN triplets for the given NDPROP and NSHLP
        ipqn = 0
        do iset = 1, ndprod
            nso = 0
            do ishl = 1, nshlp(iset)
                ipqn = ipqn + 1
                if (pqn(1, ipqn) == -1) then
                    write (nftw, err_missing_pqn) ipqn, ishl, iset
                    stop 1
                end if
                irr = mshl(ipqn)
                if (irr < 0 .or. irr >= nu) then
                    write (nftw, err_invalid_mshl) irr, ipqn
                    stop 1
                end if
                if (any(pqn(:, ipqn) > nob(irr + 1))) then
                    write (nftw, err_out_of_orbitals) ipqn, nob(irr + 1)
                    stop 1
                end if
                if (pqn(1, ipqn) /= 0) then
                    nso = nso + shlmx1(irr + 1)  ! direct orbital reference mode (pqn = K,0,0)
                else
                    nso = nso + (pqn(3, ipqn) - pqn(2, ipqn) + 1)*shlmx1(irr + 1)  ! sequence mode (pqn = 0,I,J)
                end if
            end do
            if (nelecp(iset) > nso) then
                write (nftw, err_few_shells) iset, nelecp(iset), nso
                stop 1
            end if
        end do

        ! check non-overlapping orbital ranges (provided that all non-empty PQNs use the sequence mode, i.e. 0,I,J)
        ipqn = 0
        npqn = sum(nshlp(1:ndprod))
        if (all(pqn(1, 1:npqn) == 0)) then
            ! scan all PQN triplets
            do iset = 1, ndprod
                do ishl = 1, nshlp(iset)
                    ipqn = ipqn + 1
                    ! check only non-empty ones
                    if (nelecp(iset) > 0) then
                        do iorb = pqn(2, ipqn), pqn(3, ipqn)
                            ! scan all previous PQN triplets
                            jpqn = 0
                            do jset = 1, ndprod
                                do jshl = 1, nshlp(jset)
                                    jpqn = jpqn + 1
                                    ! check only non-empty ones with the same symmetry
                                    if (jpqn < ipqn .and. nelecp(jset) > 0 .and. mshl(jpqn) == mshl(ipqn)) then
                                        do jorb = pqn(2, jpqn), pqn(3, jpqn)
                                            if (iorb == jorb) then
                                                write (nftw, err_overlapping_pqns) iorb, mshl(jpqn), ipqn, jpqn
                                                stop 1
                                            end if
                                        end do
                                    end if
                                end do
                            end do
                        end do
                    end if
                end do
            end do
        end if

    end subroutine wfngrp_sanity_check


    !> \brief Store the wave function to file and reset arrays
    !>
    !> Writes all determinants present in the CONGEN global arrays (i.e. for the current CSFs) to the output file.
    !> The data written are as follows:
    !> \verbatim
    !>     1. NOI           ... number of states
    !>     2. NODI(1:NOI)   ... number of determinants per state
    !>     3. NID           ... number of all determinants
    !>     4. CR(1:NID)     ... determinant multiplicative factors
    !>     5. NI            ... length of integer array defining all determinants
    !>     6. NR(1:NI)      ... integer array defining all determinants as differences w.r.t. reference determinant
    !> \endverbatim
    !>
    !> Note that this is not the only moment when these data are written. The subroutine \ref congen_distribution::wfn
    !> does writing, too, whenever it detects that any further addition of determinants would exceed the storage.
    !>
    !> Once the writing is done, this subroutine erases all determinant data.
    !>
    !> \param ia     On output, accumulated number of written determinants (value of \ref congen_data::lcdi).
    !> \param ib     On output, accumulated number of written integers describing the determinants (value of \ref congen_data::lndi).
    !> \param megul  Output file unit number.
    !> \param nndel  How many CSFs were read from standard input. If non-zero, nothing will be written to the file and the
    !>               variables \c ia and \c ib will be assigned just the current values of \c lcdi and \c lcdt.
    !> \param cr     Real workspace array; determinant coefficients are expected to be located at index \c icdi and onwards.
    !> \param nr     Real workspace array; determinant array pointers are expected to be located at position \c inodi and onwards,
    !>               while packed determinant descriptions are expected at position \c indi and onwards. More specifically, this is
    !>               true if integer has the same byte size as the real data type. Otherwise the offsets are adjusted so that
    !>               they allow for \c inodi - 1 or \c indi - 1 real numbers before start of the determinant integer data.
    !>
    subroutine csfout (ia, ib, megul, nndel, cr, nr)

        use precisn,     only : wp
        use congen_data, only : lratio, iidis2, lcdi, lndi, ni, nid, noi, icdi, indi, inodi

        integer :: ia, ib, megul, nndel
        real(kind=wp), dimension(*) :: cr
        integer, dimension(*) :: nr
        intent (in) cr, megul, nndel, nr
        intent (out) ia, ib

        integer :: i

        if (noi /= 0 .or. ni /= 0 .or. nid /= 0) then

            if (nndel == 0 .or. iidis2 /= 0) then
                write(megul) noi, (nr(i+(inodi-1)*lratio), i=1,noi)
                write(megul) nid, (cr(i+icdi-1),           i=1,nid)
                write(megul) ni,  (nr(i+(indi-1)*lratio),  i=1,ni)
                lcdi = lcdi + nid
                lndi = lndi + ni
            end if

            noi = 0
            nid = 0
            ni  = 0

        end if

        ia = lcdi
        ib = lndi

    end subroutine csfout


    !> \brief Check tcon data and form refcon array
    !>
    subroutine getcon (ntcon, nshcon, nrcon, nshgmx, npmax, nc, nelecg, nsym, nobt, nob, nobi, nsoi, shlmx1, exref, tcon,   &
                      refcon, error)

        logical :: error
        integer :: nc, nelecg, nobt, npmax, nshgmx, nsym, ntcon
        integer, dimension(*) :: exref, nob, nobi, nrcon, nshcon, nsoi, shlmx1
        integer, dimension(nobt,*) :: refcon
        integer, dimension(3,nshgmx,*) :: tcon ! JMC the dimensions could be changed to (3,JX,JZ) if desired.
        intent (in) exref, nelecg, nob, nobi, nobt, npmax, nrcon, nshgmx, nsoi, nsym, ntcon, shlmx1, tcon
        intent (out) error
        intent (inout) nc, nshcon, refcon

        integer :: ic, j, k, nel, nesym, net, netc, noc, ns, nshcr, nspf, nspi, pqnt, pqntm, symt

        error = .false.
        nc = 0
        if (ntcon == 0) return
        error = .true.
        do ic = 1, ntcon
            nc = nc + 1
            do j = 1, nobt
                refcon(j,ic) = 0
            end do
            nshcr = nshcon(ic)
            if (nshcr <= nshgmx) then
            else
                nshcon(ic) = npmax
                return
            end if
            if (nrcon(ic) < 0 .or. nrcon(ic) > nelecg) return
            netc = 0
            do j = 1, nshcr
                symt = tcon(1,j,ic) + 1
                pqnt = tcon(2,j,ic)
                net = tcon(3,j,ic)
                if (symt <= 0 .or. symt > nsym) return
                if (net <= 0) return
                netc = netc + net
                nesym = shlmx1(symt)
                pqntm = pqnt + (net - 1) / nesym
                if (pqnt <= 0 .or. pqntm > nob(symt)) return
                nspi = nsoi(symt) + (pqnt - 1) * nesym
                nspf = nspi + ((net + 1) / nesym) * nesym - 1
                do ns = nspi, nspf
                if (exref(ns) /= 0) return
                end do
                noc = nobi(symt) + pqnt - 1
                nel = net
                do k = 1, nel, nesym
                noc = noc + 1
                if (refcon(noc,ic) /= 0) return
                refcon(noc,ic) = min(nesym, net)
                net = net - nesym
                end do
            end do
            if (netc /= nelecg) return
        end do
        error = .false.

    end subroutine getcon


    !> \brief Set up coupling scheme
    !>
    subroutine getcup (nshlt, def, nd, nshlp, cup, error)

        integer :: def, nd, nshlt
        logical :: error
        integer, dimension(3,nshlt) :: cup
        integer, dimension(nd) :: nshlp
        intent (in) def, nd, nshlp, nshlt
        intent (out) error
        intent (inout) cup

        integer :: i, i1cup, i2cup, ifc, ifcup, ii, iicup, itest, j, nc, ncup, ncup2, ndm1, ns1, ns2, nsc1, nsc2

        error = .false.
        if (nshlt == 1) return

        if (def == 0) then

            ifcup = nshlt
            ifc = 0
            i1cup = 1
            do i = 1, nd
                iicup = nshlp(i) - 1
                if (iicup /= 0) then
                    i2cup = i1cup
                    do ii = 1, iicup
                        ifc = ifc + 1
                        ifcup = ifcup + 1
                        i2cup = i2cup + 1
                        cup(1,ifc) = i1cup
                        cup(2,ifc) = i2cup
                        cup(3,ifc) = ifcup
                        i1cup = ifcup
                    end do
                    i1cup = i2cup
                end if
                i1cup = i1cup + 1
            end do
            ndm1 = nd - 1

            ! complete shell to shell couplings
            ns1 = 1
            nsc1 = nshlp(1) - 1
            do i = 1, ndm1
                ns2 = ns1 + nshlp(i)
                nsc2 = nsc1 + (nshlp(i + 1) - 1)
                ifc = ifc + 1
                ifcup = ifcup + 1
                i1cup = ns1
                i2cup = ns2
                if (nshlp(i) /= 1) i1cup = cup(3,nsc1)
                if (nshlp(i+1) /= 1) i2cup = cup(3,nsc2)
                if (i > 1) i1cup = cup(3,ifc-1)
                cup(1,ifc) = i1cup
                cup(2,ifc) = i2cup
                cup(3,ifc) = ifcup
                ns1 = ns2
                nsc1 = nsc2
            end do

        end if

        ! check cup array for allowed values
        error = .true.
        ncup = nshlt - 1
        itest = nshlt
        do i = 1, ncup
            itest = itest + 1
            if (cup(3,i) /= itest) return
            if (cup(1,i) >= itest .or. cup(2,i) >= itest) return
        end do

        ncup2 = ncup + ncup
        do i = 1, ncup2
            nc = 0
            do j = 1, ncup
                if (cup(1,j) == i) nc = nc + 1
                if (cup(2,j) == i) nc = nc + 1
            end do
            if (nc /= 1) return
        end do
        error = .false.
        
    end subroutine getcup


    !> \brief Form reference list of spin-orbital numbers
    !>
    !> This subroutine construct the "reference determinant", which is nothing else than a list of unique numbers, one for each
    !> electron, where each number corresponds to one of molecular spin-orbitals defined by the \c &state namelist of CONGEN.
    !> The spin-orbitals are labeled sequentially across all symmetries, where all spin-orbitals sharing an orbital (but differing
    !> by spin) are kept together with consecutive indices. The order of spins within groups corresponding to a single orbital
    !> is "alpha spin, then beta spin". On successful return, the array \c refdet contains the reference list with indices sorted
    !> from smallest to largest.
    !>
    !> This subroutine can run into several error conditions that are indicated by logical flags \c e1, \c e2, \c e3 (.true.
    !> means that the error occured). They are all related to inconsistencies in the input file.
    !>
    !> \param reforb   Reference determinant quintets specification (1 = M-value of shell, 2 = \c pqn, 3 = number of electrons in shell,
    !>                 4 = spin orb 1, 5 = spin orb 2)
    !> \param refgu    Inversion symmetry numbers for orbitals (+1 = gerade, -1 = ungerade). Only used for \c symtyp = 1 (D_infh).
    !> \param nrefo    Number of quintets in \c reforb.
    !> \param nelec    Input value of number of electrons.
    !> \param refdet   List of spin-orbital numbers (one for each electron, in increasing order) in reference determinant.
    !> \param nelr     Computed number of electrons from \c reforb(3,i).
    !> \param nsoi     Initial value of spin orb 1 for all symmetries (index runs through all symmetries).
    !> \param nob      ?
    !> \param symtyp   Point group (0 = C_infv, 1 = D_infh, 2 = D_2h,C_2v,C_s,E).
    !> \param shlmx    Shell occupancy limit for each symmetry.
    !> \param nsym     Number of available irreducible representations.
    !> \param nrfomx   Max number of orbitals for ref state in arrays \c reforb and \c refgu.
    !> \param e1       Output error flag: \c nrefo out of range.
    !> \param e2       Output error flag: \c nelr does not match \c nelec.
    !> \param e3       Output error flag: Errors in \c reforb and \c refgu arrays.
    !>
    subroutine getref (reforb, refgu, nrefo, nelec, refdet, nelr, nsoi, nob, shlmx, nsym, symtyp, nrfomx, e1, e2, e3)

        logical :: e1, e2, e3
        integer :: nelec, nelr, nrefo, nrfomx, nsym, symtyp
        integer, dimension(*) :: nob, nsoi, refdet, shlmx
        integer, dimension(*) :: refgu
        integer, dimension(5,*) :: reforb
        intent (in) nelec, nob, nrefo, nrfomx, nsoi, nsym, refgu, reforb, shlmx, symtyp
        intent (out) e1, e2, e3
        intent (inout) nelr, refdet

        integer :: i, isot, j, jj, ne, neb, nelrm1, neo, ner, nesr, pqnr, pqnrm, symr, t1, tgu

        ! so far no errors
        e1 = .false.
        e2 = .false.
        e3 = .false.

        ! make sure that 'nrefo' is in interval 1 .. 'nrfomx'
        if (nrefo <= 0 .or. nrefo > nrfomx) then
            e1 = .true.
            return
        end if

        ! make sure that there are really 'nelec' electrons used in quintets
        nelr = sum(abs(reforb(3,1:nrefo)))
        if (nelr /= nelec) then
            e2 = .true.
            return
        end if

        ! loop over all quintets
        nelr = 0
        do i = 1, nrefo

            ! unpack quintet data
            symr = reforb(1,i) + 1
            pqnr = reforb(2,i)
            ne   = reforb(3,i)

            ! D_infh-specific code
            if (symtyp == 1) then
                ! check if the mirror symmetry q-number is valid
                tgu = refgu(i)
                if (abs(tgu) /= 1) then
                    e3 = .true.
                    return
                end if

                ! remap symr according to this table:
                !
                !     |     original symr 
                ! tgu |  1   2   3   4   5   ...
                ! ----|-----------------------
                !  +1 |  1   4   5   8   9   ...
                !  -1 |  2   3   6   7  10   ...
                !
                if (mod(symr,2) /= 0) tgu = -tgu
                symr = 2 * symr - (1 - tgu) / 2
            end if

            ! abort on reference to non-existent symmetry
            if (symr <= 0 .or. symr > nsym .or. ne <= 0) then
                e3 = .true.
                return
            end if

            ! make sure that electrons in this quintet actually fit in provided orbitals
            nesr = shlmx(symr)
            pqnrm = pqnr + (ne - 1) / nesr
            if (pqnr <= 0 .or. pqnrm > nob(symr)) then
                e3 = .true.
                return
            end if

            ! spin-orbital indices
            isot = nsoi(symr) + (pqnr - 1) * nesr - 1   ! global index of the spin-orbital just before the starting one for this quintet
            neo = mod(ne, nesr)                         ! number of electrons in this quintet that are in open shell
            ner = nelr + ne - neo                       ! number of all electrons in closed shells so far
            neb = min(neo, nesr - neo)                  ! number of electrons or holes in the shell, whichever is smaller

            ! label electrons in this quintet by the next free spin-orbitals
            do j = 1, ne
                nelr = nelr + 1
                isot = isot + 1
                refdet(nelr) = isot
            end do

            ! And that's all for closed shell!
            ! ... also cycle if there are no data for treatment of open shells
            if (neo == 0) cycle
            if (all(reforb(3+1:3+neb,i) == -1)) cycle

            ! revisit the first electron that is not in closed shell
            isot = isot - neo + 1

            ! if the shell is not filled over half...
            if (neo == neb) then
                do j = 1, neb
                    if (reforb(3+j,i) < 0 .or. reforb(3+j,i) > nesr) then
                        e3 = .true.
                        return
                    end if
                    ner = ner + 1
                    ! ... just correct the electron spin-orbital indices by offset ("spin") value given by the user from reforb(4:,i)
                    refdet(ner) = isot + reforb(3+j,i)
                end do
                cycle
            end if

            ! when the shell is filled over half ...
            nesr_loop: do j = 1, nesr
                do jj = 1, neb
                    if (reforb(3+jj,i) < 0 .or. reforb(3+jj,i) > nesr) then
                        e3 = .true.
                        return
                    end if
                    if (reforb(3+jj,i) < j - 1) cycle nesr_loop
                end do
                ner = ner + 1
                ! ... do some other magic ???
                refdet(ner) = isot + j - 1
            end do nesr_loop

        end do

        ! order spin-orbital numbers in refdet (if more than one electron)
        nelrm1 = nelr - 1
        if (nelrm1 > 0) then
            do i = 1, nelrm1
                t1 = refdet(i)
                j = i + 1
                do jj = 1, ner
                    if (t1 == refdet(jj)) return
                    if (t1 <= refdet(jj)) cycle
                    t1 = refdet(jj)
                    refdet(jj) = refdet(i)
                    refdet(i) = t1
                end do
            end do
        end if

    end subroutine getref


    !> \brief Precomputes needed binomial coefficients
    !>
    !> Calculates the Pascal triangle of a given order. The binomial coefficients are stored
    !> in the array \c binom as concatenated lines of the Pascal triangle, where the consecutive
    !> lines share the common value of 1. The array \c ind contains pointers to line beginnings.
    !> The leading elements of the two arrays are:
    !> \verbatim
    !>    ind = 1, 1, 2, 4, 7, ...
    !>  binom = 1, 1, 2, 1, 3, 3, 1, 4, 6, 4, 1, ...
    !> \endverbatim
    !>
    subroutine icgcf

        use precisn,     only : wp
        use consts,      only : one => xone
        use congen_data, only : binom, ind, jsmax

        integer :: i, jj, js, lb, lb1

        ind(1) = 1
        ind(2) = 1
        binom(1) = one
        binom(2) = one
        lb = 3
        lb1 = 1
        js = jsmax + 1
        do i = 2, js
            do jj = 2, i
                binom(lb) = binom(lb1) + binom(lb1 + 1)
                lb = lb + 1
                lb1 = lb1 + 1
            end do
            ind(i + 1) = lb1
            binom(lb) = one
            lb = lb + 1
        end do

    end subroutine icgcf


    !> \brief Write information obtained from &state
    !>
    !> This subroutine is called by \ref csfgen after the \c &state namelist is read, workspace allocated and reference
    !> determinant formed. Its purpose is to summarize the so far read input data on standard output (or other unit specified
    !> by \c nftw) and also to print information about the reference determinant.
    !>
    !> This subroutine does not alter any of its arguments except for \c error, which is set to .true. if an error condition
    !> occured and to .false. otherwise.
    !>
    subroutine stwrit (nelect, confpf, qntot, cdimx, icdi, ntso, symtyp, ndimx, indi, nrefo, nodimx, inodi, nsym, gutot, nbmx,  &
                       isz, navail, idiag, megu, thres, lcdt, megul, lndt, nfto, nrerun, ltri, npflg, nndel, nob, nsoi, nsymp,  &
                       refdet, nerfs, erfs, nrefop, reforb, refgu, nelp, lpp, sname, error, byproj, lndo, lcdo, iposit, nob0,   &
                       npmult, ntgsym, mxtarg, nobe, nobp, nobv)

        use precisn,      only : wp
        use global_utils, only : mprod
        use congen_data,  only : nftw, nu, rhead
        use congen_pagecontrol, only : addl, ctlpg1, newpg, space

        integer :: byproj, cdimx, confpf, gutot, icdi, idiag, indi, inodi, &
                   iposit, isz, lcdo, lcdt, lndo, lndt, lpp, ltri, megu,   &
                   megul, mxtarg, navail, nbmx, ndimx, nelect, nelp,       &
                   nerfs, nfto, nndel, nodimx, npmult, nrefo, nrefop,      &
                   nrerun, nsym, nsymp, ntgsym, ntso, symtyp
        logical :: error
        character(len=80) :: sname
        real(kind=wp) :: thres
        logical, dimension(9)   :: erfs
        integer, dimension(nu)  :: nob, nob0, nobe, nobp, nobv, nsoi ! JMC changing the dimension from 10
        integer, dimension(6)   :: npflg
        integer, dimension(3)   :: qntot
        integer, dimension(*)   :: refdet, refgu
        integer, dimension(5,*) :: reforb
        intent (in) byproj, cdimx, confpf, erfs, gutot, icdi, idiag, indi, &
                    inodi, iposit, isz, lcdo, lcdt, lndo, lndt, ltri,      &
                    megu, megul, mxtarg, navail, nbmx, ndimx, nelect,      &
                    nelp, nerfs, nfto, nndel, nob, nob0, nobe, nobp, nobv, &
                    nodimx, npflg, nrefo, nrefop, nrerun, nsoi, nsym,      &
                    nsymp, ntgsym, ntso, qntot, refdet, refgu, reforb,     &
                    symtyp, thres
        intent (inout) error

        character(len=32), dimension(9) :: ersnts
        character(len=30) :: head = 'CONGEN 1.0  IBM SAN JOSE      '
        integer :: i, ii, imax, ip, it, junk, lsn = 64, nitem = 30

        data ersnts/'SYMMETRY TYPE OUT OF RANGE      ', &
                    'NO ORBITALS GIVEN               ', &
                    'EXDET, EXREF ALLOCATION FAILED  ', &
                    'NELECT OUT OF RANGE             ', &
                    'NREFO OUT OF RANGE              ', &
                    'SUM NELEC IN REF ORBS NE NELECT ', &
                    'ERROR IN REFORB DATA            ', &
                    'ERROR IN TOTAL QN DATA          ', &
                    'NO CORE FOR CDI, NDI, AND NODI  '/

        call ctlpg1 (lpp, head, len(head), sname, lsn)
        call newpg
        call addl (15)

        write(nftw,'(T2,"NELECT",T8,I4,T15,"CONFPF",I3,T27,"MULT  ",I2,T38,"CDIMX ",I5,T52,"ICDI ",I6)') &
            nelect, confpf, qntot(1), cdimx, icdi
        write(nftw,'(" NTSO  ",I4,T15,"SYMTYP",I3,T27,"MVAL  ",I2,T38,"NIDMX ",I5,T52,"INDI ",I6)') &
            ntso, symtyp, qntot(2), ndimx, indi
        write(nftw,'(" NREFO ",I4,T27,"REFLC",I3,T38,"NODIMX",I5,T52,"INODI",I6)') nrefo, qntot(3), nodimx, inodi
        write(nftw,'(" NSYM  ",I4,T27,"GUTOT",I3,T38,"NCORE",I10)') nsym, gutot, nbmx
        write(nftw,'(T27,"ISZ  ",I3,T38,"NAV  ",I10)') isz, navail
        
        write(nftw,'(//,T14," DATA FOR SPEEDY INPUT")')
        write(nftw,'(/, T14,"IDIAG ",I4,T27,"MEGU ",I3,T38,"THRES ",1PE12.4)') idiag, megu, thres
        write(nftw,'(T14,"LCDT  ",I4,T27,"MEGUL",I3)') lcdt, megul
        write(nftw,'(T14,"LNDT  ",I4,T27,"NFTO ",I3)') lndt, nfto
        write(nftw,'(T14,"NRERUN",I4,T27,"LTRI ",I3)') nrerun, ltri
        write(nftw,'(/,T14,"NPFLG =",6I3,2X,"NNDEL =  ",I5)') npflg, nndel

        call addl (1)

        write(nftw,'(T14,"BYPROJ",I2,3x,"LNDO",I10,3x,"LCDO",I10)') byproj, lndo, lcdo

        if (ntgsym < mxtarg) then
            call addl (1)
            write(nftw,'(" ntgsym",I4)') ntgsym
        end if

        call space (2)
        call addl (3)

        write(nftw,'(" NSYM",30I5)') (ip, ip=1,nsymp)
        write(nftw,'(" NOB ",30I5)') (nob(ip), ip=1,nsymp)
        if (iposit /= 0) then
            call addl (4) ! JMC adding this line
            write(nftw,'(" NOB0",30I5)') (nob0(ip), ip=1,nsymp)
            write(nftw,'(" NOBE",30I5)') (nobe(ip), ip=1,nsymp)
            write(nftw,'(" NOBP",30I5)') (nobp(ip), ip=1,nsymp)
            write(nftw,'(" NOBV",30I5)') (nobv(ip), ip=1,nsymp)
        end if
        write(nftw, '(" NSOI",30I5)') (nsoi(ip), ip=1,nsymp)
        call space (1)
        if (iposit /= 0) then
            call addl (1) ! JMC adding this line
            write(nftw,'(5X,"POSITRON SCATTERING CASE: IPOSIT =",I3)') iposit
            call space (1)
        end if

        it = 1
        if (symtyp == 1) it = 2
        do i = 1, nrefop, nitem
            imax = min(i + nitem - 1, nrefop)
            call addl (6) ! JMC the argument will be too small in some cases (i=1 and symtyp=1)???
            if (i == 1) then
                write(nftw,'(" REFERENCE DETERMINANT INPUT DATA")')
                it = it - 1
            end if
            write(nftw,'(1X,A4,I5,29I4)') rhead(1), (ip,ip=i,imax)
            do ii = 1, 5
                write(nftw, '(1X,A4,I5,29I4)') rhead(ii+1), (reforb(ii,ip), ip=i,imax)
            end do
            if (symtyp == 1) write(nftw, '(1X,A4,I5,29I4)') rhead(7), (refgu(ip), ip=i,imax)
            call space (1)
        end do
        call space (1)
        it = (nelp + nitem - 1) / nitem
        if (mod(nelp,nitem) == 0) it = it + 1

        if (nelp /= 0) then
            call addl (it)
            write(nftw,'(" REFDET =",30(I3,",")/(9X,30(I3,",")))') (refdet(ip), ip=1,nelp)
        end if

        ! print D_2h multiplication table
        if (symtyp >= 2 .and. npmult /= 0) then
            call addl (25)
            junk = mprod(1, 1, npmult, nftw)
        end if

        ! process &state errors
        error = .false.
        do i = 1, nerfs
            if (erfs(i)) then
                if (.not. error) then
                    call space (2)
                    call addl (2)
                    write(nftw,'(" **** ERROR DATA FOR &STATE FOLLOWS (&WFNGRP NOT PROCESSED)"/12X,A32)') ersnts(i)
                    error = .true.
                else
                    call addl (1)
                    write(nftw,'(12X,A32)') ersnts(i)
                end if
            end if
        end do

    end subroutine stwrit


    !> \brief Read CSFs from input stream.
    !>
    subroutine subdel (ndel, ndel1, ndel2, nndel)

        use congen_data, only : nftw, nftr

        integer :: nndel
        integer, dimension(*) :: ndel, ndel1, ndel2 ! JMC changing dimension from 2.
        intent (inout) ndel, ndel1, ndel2, nndel

        integer :: i, j, k, m, nndel1

        read(nftr,'(16I5)') nndel
        read(nftr,'(16I5)') (ndel(i), i=1,nndel)

        read(nftr,'(16I5)') nndel1

        do while (nndel1 /= 0)

            read(nftr,'(16I5)') (ndel1(j), j=1,nndel1)

            i = 1
            j = 1
            k = 0

            read_loop: do

                if (ndel(i) <= ndel1(j)) then
                    k = k + 1
                    if (ndel(i) == ndel1(j)) j = j + 1
                    ndel2(k) = ndel(i)
                    i = i + 1
                    if (i > nndel) then
                        ndel(1:k) = ndel2(1:k)
                        do m = j, nndel1
                            k = k + 1
                            ndel(k) = ndel1(m)
                        end do
                        exit read_loop
                    end if
                else
                    k = k + 1
                    ndel2(k) = ndel1(j)
                    j = j + 1
                end if

                if (j > nndel1) then
                    do m = i, nndel
                        k = k + 1
                        ndel2(k) = ndel(m)
                    end do
                    ndel(1:k) = ndel2(1:k)
                    exit read_loop
                end if

            end do read_loop

            nndel = k
            read(nftr,'(16I5)') nndel1

        end do

        write(nftw,'(" ***** NUMBER OF CHOSEN CONFIGURATIONS ***",I10)') nndel
        write(nftw,'(/,24I5)')(ndel(i),i=1,nndel)
        write(7,'(16I5)') (ndel(i), i=1,nndel) ! JMC writing to unit 7 (hardwired)???

    end subroutine subdel


    !> \brief Set defaults for wave function group parameters
    !>
    subroutine wfnin (nwfngp, nadel, nncsf, ncsf, lcdi, lndi, nelecg, ndprod, nrefog, npcupf, negmax, refdtg, nrfgmx, refgug, &
                      ntcon, reforg, nshgmx, nsymmx, mshl, gushl, pqn, cup, ndpmax, nshlp, nconmx, test, nrcon, nshcon, tcon, &
                      errorg)

        logical :: errorg
        integer :: lcdi, lndi, nadel, nconmx, ncsf, ndpmax, ndprod, negmax, nelecg, nncsf, npcupf, nrefog, nrfgmx, nshgmx, &
                   nsymmx, ntcon, nwfngp
        integer, dimension(3,*) :: cup, pqn
        integer, dimension(*)   :: gushl, mshl, nrcon, nshcon, nshlp, refdtg, refgug, test
        integer, dimension(5,*) :: reforg
        integer, dimension(3,nshgmx,*) :: tcon ! JMC the dimensions could be changed to (3,JX,JZ) if desired.
        intent (in)  nconmx, ndpmax, negmax, nrfgmx, nshgmx, nsymmx
        intent (out) cup, errorg, gushl, lcdi, lndi, mshl, nadel, ncsf, ndprod, nelecg, nncsf, npcupf, nrcon, nrefog, nshcon, &
                     nshlp, ntcon, nwfngp, pqn, refdtg, refgug, reforg, tcon, test

        integer :: i, j, k

        ! general default for wfngrp arrays

        nwfngp = 0
        nadel  = 1
        nncsf  = 0
        ncsf   = 0
        lcdi   = 0
        lndi   = 0
        ntcon  = 0

        nelecg = 0
        ndprod = 0
        nrefog = 0
        npcupf = 0

        refdtg(1:negmax)     =  0
        refgug(1:nrfgmx)     = -2
        reforg(1:5,1:nrfgmx) = -2

        mshl(1:nshgmx)  = nsymmx + 1
        gushl(1:nshgmx) = -2

        pqn(1:3,1:nshgmx) = -1
        cup(1:3,1:nshgmx) = -1

        errorg = .false.

        nshlp(1:ndpmax)  =  0
        test(1:nconmx)   =  1
        nrcon(1:nconmx)  = -1
        nshcon(1:nconmx) =  0

        tcon(1,1:nshgmx,1:nconmx) = -1
        tcon(2,1:nshgmx,1:nconmx) =  0
        tcon(3,1:nshgmx,1:nconmx) =  0

    end subroutine wfnin


    !> \brief Defaults to be reset before every wfngrp
    !>
    !> This subroutine is called from \ref csfgen before every attempt to read next \c &wfngrp namelist.
    !> It erases all previous values to avait new input.
    !>
    !> \param nelecp  Number of electrons per electron set within the wave function group.
    !> \param defltc  Use of user-specified intermediate coupling.
    !> \param nerfg   Size of the array \c erfg.
    !> \param erfg    Array of error condition flags.
    !> \param gname   Arbitrary namelist text description (max 80 characters).
    !> \param qntar   Triplet of total quantum numbers for the wave function.
    !> \param errorg  Error flag (not used).
    !> \param ndpmax  Size of the array \c nelecp.
    !>
    subroutine wfnin0 (nelecp, defltc, nerfg, erfg, gname, qntar, errorg, ndpmax)

        integer :: defltc, ndpmax, nerfg
        logical :: errorg
        character(len=80) :: gname
        logical, dimension(*) :: erfg
        integer, dimension(*) :: nelecp
        integer, dimension(3) :: qntar
        intent (in) ndpmax, nerfg
        intent (out) defltc, erfg, gname, nelecp, qntar

        nelecp(1:ndpmax) = -1
        erfg(1:nerfg)    = .false.
        
        gname  = '        '
        defltc = 0

        qntar(1) = -1
        qntar(2) =  0
        qntar(3) =  0

    end subroutine wfnin0


    !> \brief Final text information output.
    !>
    subroutine wrnmlt (k, sname, nrerun, megul, symtyp, mgvn, s, sz, r, pin, ncsf, byproj, lcdi, lndi, lcdo, lndo, lcdt, lndt, &
                       nfto, ltri, idiag, npflg, thres, nelect, nsym, nob, refdet, nftw, iposit, nob0, nobl, nob0l, nx, nobe,  &
                       nobp, nobv)

        use precisn, only : wp

        integer :: byproj, idiag, iposit, k, lcdi, lcdo, lcdt, lndi, lndo, lndt, ltri, megul, mgvn, ncsf, nelect, nfto, nftw, &
                   nrerun, nsym, nx, symtyp
        real(kind=wp) :: pin, r, s, sz, thres
        character(len=80) :: sname
        integer, dimension(*) :: nob, nob0, refdet
        integer, dimension(nx) :: nob0l, nobl
        integer, dimension(nx) :: nobe, nobp, nobv
        integer, dimension(6) :: npflg
        intent (in) byproj, idiag, iposit, k, lcdi, lcdo, lcdt, lndi, lndo, lndt, ltri, megul, mgvn, ncsf, nelect, nfto, &
                    nob, nob0, nob0l, nobe, nobl, nobp, nobv, npflg, nrerun, nsym, nx, pin, r, refdet, s, sname, symtyp, &
                    sz, thres

        character(len=4) :: blank1 = '    '
        integer :: j

        write(k,'(" &INPUT")')
        write(k,'(" NAME=''",A80,"'',")') sname
        write(k,'(" NRERUN=",I3,",  MEGUL=",I3,",  SYMTYP=",I3,",")') nrerun, megul, symtyp
        write(k,'(" MGVN=",I3,",  S=",F6.1,",SZ=",F6.1,", R=",F6.1,",  PIN=",F6.1,",  NOCSF=",I6,",")') mgvn, s, sz, r, pin, ncsf
        write(k,'(" BYPROJ=",I2,",")') byproj
        write(k,'(" LCDI=",I15,",  LNDI=",I15,", LCDO=",I7,", LNDO=",I15,","," LCDT=",I7,",  LNDT=",I7,",")') &
            lcdi, lndi, lcdo, lndo, lcdt, lndt
        write(k,'(" NFTO=",I3,",  LTRI=",I5,", IDIAG=",I3,", NPFLG=",5(I2,",")," NPMSPD =",I2,",")') nfto, ltri, idiag, npflg
        write(k,'(" THRES=",1PD9.2,",")') thres
        write(k,'(" NELT=",I4,", NSYM=",I3,", NOB=",10(1x,I0,","))') nelect, nsym, (nob(j), j=1,nsym)
        write(k,'(" NDTRF=",12(1x,I0,",")/(8X,12(1x,I0,",")))') (refdet(j), j=1,nelect)
        write(k,'(" NOBL=",10(1x,I0,","))') nobl
        write(k,'(" NOB0L=",10(1x,I0,","))') nob0l

        if (iposit /= 0) then
            write(k,'(" IPOSIT=",I3)') iposit
            write(k,'(" NOB0=",10(1x,I0,","))') (nob0(j), j=1,nsym)
            write(k,'(" NOBE=",10(1x,I0,","))') (nobe(j), j=1,nsym)
            write(k,'(" NOBP=",10(1x,I0,","))') (nobp(j), j=1,nsym)
            write(k,'(" NOBV=",10(1x,I0,","))') (nobv(j), j=1,nsym)
        end if
        write(k,'(" /")')

    end subroutine wrnmlt


    !> \brief Display information about the wave function group
    !>
    subroutine wvwrit (nwfngp, gname, nelecg, defltc, irfcon, ndprod, symtyp, ntcon, navail, nrefog, nelecp, nshlp, qntar, &
                       nshlpt, nshgmx, mshl, gushl, pqn, cup, ncupp, npcupf, refdtg, nelp, ntconp, nshcon, nobt, reforb,   &
                       refgu, test, nrcon, tcon, refcon, nerfg, erfg, ndpp, nrefop, lsquare, errorg)

        use congen_data, only : nftw, rhead
        use congen_pagecontrol, only : addl, newpg, space

        integer :: defltc, irfcon, navail, ncupp, ndpp, ndprod, nelecg, nelp, nerfg, nobt, npcupf, nrefog, nrefop, &
                nshgmx, nshlpt, ntcon, ntconp, nwfngp, symtyp, lsquare
        logical :: errorg
        character(len=80) :: gname
        integer, dimension(3,*) :: cup, pqn
        logical, dimension(11) :: erfg
        integer, dimension(*) :: gushl, mshl, nelecp, nrcon, nshcon, nshlp, refcon, refdtg, refgu, test
        !     change refcon from integer to couble for consistency
        !      double precision  refcon(*) ! jmc the correponding actual arg is d.p. ???
        integer, dimension(3) :: qntar
        integer, dimension(5,*) :: reforb
        integer, dimension(3,nshgmx,*) :: tcon ! jmc changing the 2nd dimension from LG to NSHGMX for consistency with elsewhere.
                                            ! The dimensions could also be changed to (3,JX,JZ) if desired.
        intent (in) cup, defltc, erfg, gname, gushl, irfcon, mshl, navail, ncupp, ndpp, ndprod, nelecg, nelecp, nelp,         &
                    nerfg, nobt, npcupf, nrcon, nrefog, nrefop, nshcon, nshgmx, nshlp, nshlpt, ntcon, ntconp, nwfngp, pqn,    &
                    qntar, refcon, refdtg, refgu, reforb, symtyp, tcon, lsquare
        intent (inout) errorg, test

        character(len=32), dimension(11) :: ersntg
        character(len=4) :: hcup = 'CUP ', htcon = 'TCON', lpc = '  ( '
        integer :: i, i1, i2, ic, ii, iimax, imax, ip, ishp, it, jj, jnshl, nshcr, nit1 = 10, nit2 = 30

        data ersntg/'NELECG OUT OF RANGE             ', &
                    'NDPROD TOO LARGE                ', &
                    'NSHLP OUT OF RANGE              ', &
                    'NELECG NE SUM OF NELEP          ', &
                    'NREF OUT OF RANGE               ', &
                    'NELECG NE SUM OVER NELEC IN REFO', &
                    'ERROR IN REF ORB DATA           ', &
                    'ERROR IN SHELL DATA             ', &
                    'ERROR IN COUPLING DATA          ', &
                    'ERROR IN CONFIG CONSTRAINT DATA ', &
                    'NO SPACE FOR REFCON ARRAY       '/

        ! print wfngrp input data
        call newpg
        call addl (9)
        write(nftw,'(" WFN GROUP",I4,4X,A80)') nwfngp, gname
        write(nftw,'(" NELECG",I4,T15,"DEFLTC",I3,T27,"IRFCON",I6)') nelecg, defltc, irfcon
        write(nftw,'(" NDPROD",I4,T15,"NTCON ",I3,T27,"NAV  ",I10)') ndprod, ntcon, navail
        write(nftw,'(" LSQUARE",I3)') lsquare
        write(nftw,'(" NREFOG",I4,/)') nrefog
        write(nftw,'(9X,9I3)') (i, i=1,ndpp)
        write(nftw,'(" NELECP =",9I3)') (nelecp(i), i=1,ndpp)
        write(nftw,'(" NSHLP  =",9I3)') (nshlp(i), i=1,ndpp)

        if (qntar(1) /= -1) then
            call addl (1)
            write(nftw,'(5X,"TARGET MULTIPLICITY =",I5,5X,"TARGET SYMMETRY =",I5,5X,"TARGET INVERSION SYMMETRY =",I5)') &
                (qntar(i), i=1,3)
        end if

        call space (1)

        if (ndprod /= 0 .and. nshlpt <= nshgmx) then

            it = merge(1, 0, symtyp == 1)
            ishp = 0

            do ip = 1, ndprod
                jnshl = nshlp(ip)
                if (jnshl == 0) then
                    call addl (1)
                    write(nftw,'(I4," ORBITALS IN GROUP",I2)') jnshl, ip
                    call space (2)
                end if
                do i = 1, jnshl, nit1
                    imax = min(i + nit1 - 1, jnshl)
                    if (i == 1) then
                        call addl (4 + it)
                        write(nftw,'(I4," ORBITALS IN GROUP",I2)') jnshl, ip
                    else
                        call addl (3 + it)
                    end if
                    write(nftw,'(1X,A4,I8,9I12)') rhead(1), (ii, ii=i,imax)
                    write(nftw, '(1X,A4,I8,9I12)') rhead(2), (mshl(ishp+ii), ii=i,imax)
                    if (symtyp == 1) write(nftw,'(1X,A4,I8,9I12)') rhead(7), (gushl(ishp+ii), ii=i,imax)
                    write(nftw,'(1X,A4,10(A3,I2,",",I2,",",I2,")"))') rhead(3), (lpc,(pqn(jj,ishp+ii), jj=1,3), ii=i,imax)
                    call space (1)
                end do
                ishp = ishp + jnshl
            end do

            call space (1)

            if (ncupp * npcupf /= 0) then
                do i = 1, ncupp, nit1
                    imax = min(i + nit1 - 1, ncupp)
                    i1 = i + nshlpt
                    i2 = imax + nshlpt
                    if (i == 1) then
                        call addl (3)
                        write(nftw,'(" COUPLING DATA")')
                    else
                        call addl (2)
                    end if
                    write(nftw,'(1X,A4,I8,9I12)') rhead(1), (ii, ii=i1,i2)
                    write(nftw, '(1X,A4,10(A3,I2,",",I2,",",I2,")"))') hcup, (lpc,(cup(jj,ii), jj=1,3), ii=i,imax)
                    call space (1)
                end do
            end if

            it = merge(2, 1, symtyp == 1)

            do i = 1, nrefop, nit2
                imax = min(i + nit2 - 1, nrefop)
                call addl (6 + it)
                if (i == 1) then
                    write(nftw,'(" REFERENCE DETERMINANT INPUT DATA")')
                    it = it - 1
                end if
                write(nftw, '(1X,A4,I5,29I4)') rhead(1), (ip, ip=i,imax)
                do ii = 1, 5
                    write(nftw, '(1X,A4,I5,29I4)') rhead(ii+1), (reforb(ii,ip), ip=i,imax)
                end do
                if (symtyp == 1) write(nftw, '(1X,A4,I5,29I4)') rhead(7), (refgu(ip), ip=i,imax)
                call space (1)
            end do

            call space (1)

            it = (nelp + nit2 - 1) / nit2
            if (mod(nelp, nit2) == 0) it = it + 1
            if (nelp /= 0) then
                call addl (it)
                write(nftw,'(" REFDET =",30(I3,","),/,(9X,30(I3,",")))') (refdtg(ip), ip=1,nelp)
            end if

            call space (1)

            do ic = 1, ntconp
                nshcr = nshcon(ic)
                it = 1 + 2 * ((nshcr + nit1 - 1) / nit1) + 1 + (nobt + nit2 - 1) / nit2 ! jmc one too many here???
                if (mod(nobt, nit2) == 0) it = it + 1
                if (ic == 1) then
                    call addl (2 + it)
                    write(nftw,'(" EXITATION CONSTRAINTS / TCON(SYM,PQN,NE",/)')
                else
                    call addl(it)
                end if
                if (test(ic) /= 1) then
                    write(nftw,'(I3,10X,"GT",I3," REPLACEMENTS ALLOWED(INTERSECTION)")') ic, nrcon(ic)
                    test(ic) = 0
                else
                    write(nftw,'(I3,10X,"LE",I3," REPLACEMENTS ALLOWED(UNION)")') ic, nrcon(ic)
                end if
                if (nshcr /= 0) then
                    do ii = 1, nshcr, nit1
                        iimax = min(ii + nit1 - 1, nshcr)
                        write(nftw,'(1X,A4,I8,9I12)') rhead(1), (ip, ip=ii,iimax)
                        write(nftw,'(1X,A4,10(A3,I2,",",I2,",",I2,")"))') htcon, (lpc,(tcon(jj,ip,ic), jj=1,3), ip=ii,iimax)
                    end do
                    call space (1)
                    write(nftw,'(" REFCON =",30(I3,",")/(9X,30(I3,",")))') (refcon(ip), ip=1,nobt)
                end if
                call space (1)
            end do

        end if

        ! print wfngrp error messages

        do i = 1, nerfg
            if (.not. erfg(i)) cycle
            if (.not. errorg) then
                call space (2)
                call addl (2)
                write(nftw,'(" **** ERROR DATA FOR &WFNGRP FOLLOWS",/,12X,A32)') ersntg(i)
                errorg = .true.
                cycle
            end if
            call addl (1)
            write(nftw,'(12X,A32)') ersntg(i)
        end do

    end subroutine wvwrit

end module congen_driver
