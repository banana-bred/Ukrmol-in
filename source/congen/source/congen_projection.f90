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
!> \brief Projection on spin states.
!>
!> Routines in this module, most prominently the central subroutine \ref projec, read back all CSFs
!> generated in the previous part of CONGEN execution and recombine the determinants to satisfy spin
!> composition rules. While doing that, several new combinations of spin-orbitals (determinants) may
!> be added if there are any open shells where the electrons spins can be freely permuted. Also, the
!> routines apply a threshold for final selection of contributing determinants, so the output can be
!> even smaller than the input (though this mostly signalizes some error in setup).
!>
module congen_projection

    implicit none

    ! entry point of the module
    public projec

    ! support routines called from "projec"
    private cntrct
    private dophz
    private dophz0
    private iphase
    private mkorbs
    private pkwf
    private pmkorbs
    private popnwf
    private prjct
    private ptpwf
    private qsort
    private rdwf
    private rdwf_getsize
    private rfltn
    private snrm2
    private stmrg
    private wfgntr
    private wrnfto
    private wrwf

contains


    !> \brief Throw away determinants with negligible contribution.
    !>
    !> Scans the store of determinants, discard such whose contribution (multiplication factor) is below
    !> given tolerance, and bubbles out the vacated intervals from the storage arrays.
    !>
    subroutine cntrct (nelt, no, ndo, cdo, thres)

        use precisn, only : wp

        integer :: nelt, no
        real(kind=wp) :: thres
        real(kind=wp), dimension(*) :: cdo
        integer, dimension(*) :: ndo
        intent (in) nelt, thres
        intent (inout) cdo, ndo, no

        integer :: i, j, md, mdd, mov

        mov = 0     ! number of determinants discarded so far due to the threshold
        md  = 0     ! current position in determinant storage array
        mdd = 0     ! position of next determinant position vacated due to threshold

        ! loop over all determinants
        do i = 1, no
            if (abs(cdo(i)) <= thres) then
                mov = mov + 1                           ! increment number of discarded determinants
            else if (mov /= 0) then
                mdd = md - nelt * mov                   ! position of next determinant position vacated due to threshold
                cdo(i-mov) = cdo(i)                     ! move the next determinant factor to vacated space
                ndo(mdd+1:mdd+nelt) = ndo(md+1:md+nelt) ! move the next determinant data to vacated space
            end if
            md = md + nelt                              ! jump to the next determinant
        end do

        no = no - mov                                   ! update number of stored determinants

    end subroutine cntrct


    !> \brief Phase factor.
    !>
    !> Compute phase factor implied by placing continuum spin-orbital after all target spin-orbitals
    !>
    !> \note MAL 10/05/2011 Changes made are in order to bring dophz into line
    !> with the changes made in 'projec' in order to utilize dynamic memory
    !>
    subroutine dophz (nftw, nocsf, nelt, ndtrf, nconf, indo, ndo, lenndo, icdo, cdo, lencdo, iphz, leniphz, iphz0, &
                      leniphz0, nctarg, nctgt, notgt, mrkorb, mdegen, ntgsym, mcont, symtyp, npflg)

        use precisn, only : wp

        integer                      :: symtyp,i,n,ntci,nt,marked,nc,mb,md,m,na,mark1,iloc,inum,  &
                                        ntci1, ntci0,nftw,nctarg,iph,npflg,nocsf,nelt,     &
                                        leniphz,leniphz0,ntgsym,lenndo,lencdo
        integer, dimension(nelt)     :: ndtrf,nconf
        integer, dimension(lenndo)   :: ndo
        integer, dimension(nocsf)    :: indo,icdo
        real(kind=wp), allocatable   :: cdo(:)
        integer, dimension(leniphz)  :: iphz
        integer, dimension(ntgsym)   :: nctgt,mrkorb,mcont,notgt,mdegen
        integer, dimension(leniphz0) :: iphz0
        real(kind=wp), parameter     :: zero = 0.0_wp
        logical, parameter           :: zdebug = .false.

        ! Debug banner header
        write(nftw,'(/,5x,"Phase analysis for total wavefunction:",/)')
        write(nftw,'(  5x," ")')
        write(nftw,'(  5x,"  Number of CSFS           (nocsf) = ",I8)') nocsf
        write(nftw,'(  5x,"  Number of electrons       (nelt) = ",I8)') nelt
        write(nftw,'(  5x,"  Number of target states (ntgsym) = ",I8)') ntgsym
        write(nftw,'(  5x,"  Spatial group type      (symtyp) = ",I8)') symtyp
        write(nftw,'(  5x,"  Size of packed dets     (lenndo) = ",I8)') lenndo
        write(nftw,'(  5x,"  Size of cdo (#dets)     (lencdo) = ",I8)') lencdo
        write(nftw,'(  5x,"                          (nctarg) = ",I8,/)') nctarg
        write(nftw,'(  5x,"Structure of wavefunction:",/)')
        write(nftw,'(  5x,"Target   #CSFs    #Continuum  Spatial Sym ")')
        write(nftw,'(  5x,"State     targ     functions  continuum   ")')
        write(nftw,'(  5x,"------  -------   ----------  ----------- ")')
        do I = 1, ntgsym
            write(nftw,'(5x,I6,2x,I7,3x,I10,2x,I10)') i, nctgt(i), notgt(i), mcont(i)
        end do
        write(nftw,'(/,5x,"**** End of structure of wavefunction",/)') 

        n = 1
        ntci = 0

        ! Descend into loop over target states
        do nt = 1, ntgsym
            marked = mrkorb(nt)

            ! Descend into loop over number of continuum orbs  
            do nc = 1, nctgt(nt)

                ! First load reference determinant
                nconf(1:nelt) = ndtrf(1:nelt)

                ! Now make substitutions
                mb = indo(n)
                md = ndo(mb)
                
                Cyc1: do m = 1, md
                    na = ndo(mb + m)
                    do i = 1, nelt
                        if (na == ndtrf(i)) then
                            nconf(i) = ndo(mb + md + m)
                            cycle Cyc1
                        end if
                    end do
                    write(nftw,*) 'DOPHZ: help I should not have got here!!!                         na =', na
                    write(nftw,*) ' ndtrf ', ndtrf
                    write(nftw,*) ' nconf ', nconf
                end do Cyc1

                ! We have the present configuration, find the marked orbital
                m = merge(4, 2, symtyp <= 1 .or. mcont(nt) /= 0) ! number of tries (different spins)
                do mark1 = marked + 0, marked + m
                    if (mark1 == marked + m) then
                        ! Not found anywhere, stop!
                        write(nftw,*) ' Configuration is ', nconf
                        stop
                    else if (any(nconf(1:nelt) == mark1)) then
                        ! Found, let's continue!
                        exit
                    end if
                end do

                ! Phase depends on where the marked orbital is in the determinant
                ntci = ntci + 1
                inum = count(nconf(1:nelt) > mark1)
                if (mdegen(nt) >= 0) then
                    iph = merge(iphase(nconf,nelt), -iphase(nconf,nelt), cdo(icdo(n)) > zero)
                    iphz(ntci) = merge(iph, -iph, MOD(inum, 2) == 0)
                else
                    ! treat phase factor caused by coupling down rather than up for second
                    ! continua in degenerate symmetry/degenerate target as special case
                    iphz(ntci) = merge(iphz0(nc), -iphz0(nc), MOD(inum, 2) == 0)
                end if
                n = n + notgt(nt)

            end do ! End of loop over target states 
        end do

        ! Having completed the computation, print results.
        if (npflg > 0) then
            write(nftw,'(//," Phase factors for CI target states:")')
            ntci1 = 0
            do nt = 1, ntgsym
                ntci0 = ntci1 + 1
                ntci1 = ntci1 + nctgt(nt)
                write(nftw,'(/," Target symmetry",I3,/)') nt
                write(nftw,'(25I3)') (iphz(ntci), ntci=ntci0,ntci1)
            end do
            write(nftw,'(//)')
        end if

        ! Subroutine return point
        if (zdebug) then
            write(nftw, '(/,5x,"***** dophz() - completed ",/)')
        end if

    end subroutine dophz


    !> \brief Phase factor.
    !>
    !> Compute phase factor for target CSFs - given by reordering spin-orbitals in ascending order.
    !>
    !> \todo MAL 10/05/2011: Changes have been made to this subroutine to bring
    !>       it into line with the changes made to 'projec' and ensure to compliance
    !>       with F95 standards
    !>
    subroutine dophz0 (nftw, nocsf, nelt, ndtrf, nconf, indo, ndo, lenndo, icdo, cdo, lencdo, iphz, npflg)

        use precisn, only : wp

        integer                         :: n, i, mb, md, m, na, nftw, npflg, nocsf, nelt, lenndo, lencdo
        integer, dimension(nelt)        :: ndtrf,nconf
        integer, dimension(nocsf)       :: indo,icdo,iphz
        integer, dimension(lenndo)      :: ndo
        real(kind=wp), parameter        :: zero = 0.0_wp
        real(kind=wp),dimension(lencdo) :: cdo(lencdo)

        do n = 1, nocsf

            ! First load reference determinant
            nconf(1:nelt) = ndtrf(1:nelt)

            ! Then make substitutions
            mb = indo(n)
            md = ndo(mb)

            Cyc1: do m = 1,md
                na = ndo(mb + m)
                do i = 1, nelt
                    if (na == ndtrf(i)) then
                        nconf(i) = ndo(mb + md + m)
                        cycle Cyc1
                    end if
                end do
                write(nftw,*) 'DOPHZ0: help I should not have got here!!!                        na =', na
                write(nftw,*) ' ndtrf ', ndtrf
                write(nftw,*) ' nconf ', nconf
                stop
            end do Cyc1

            iphz(n) = merge(iphase(nconf,nelt), -iphase(nconf,nelt), cdo(icdo(n)) > zero)

        end do

        if (npflg > 0) then
            write(nftw,'(5x,"Phz factor, per target CSF (",I7,"), for future use:")') nocsf
            write(nftw,'((5x,15(I3,1x)))') (iphz(n), n=1,nocsf)
            write(nftw,'(/)')
        end if

    end subroutine dophz0


    !> \brief Sequence phase factor
    !>
    !> Compute phase factor (if any) due to out of sequence ordering of spin-orbitals in CSF stored in nconf.
    !>
    function iphase (nconf, nelt)

        integer :: nelt
        integer :: iphase
        integer, dimension(nelt) :: nconf
        intent (in) nconf, nelt

        integer :: iso, iswap, m, n, nst

        ! first check if there are any spin-orbitals out of sequence
        do n = 1, nelt
            if (n == nelt) then
                ! all spin-orbitals are in ascending order: no phase
                iphase = 1
                return
            else if (nconf(n) > nconf(n+1)) then
                ! there is a possible phase factor!
                exit
            end if
        end do

        ! can we eliminate some electrons from the phase computation?
        nst = 1
        do n = 1, nelt
            if (nconf(n) /= n) exit
            nst = nst + 1
        end do

        ! logic says you can't reach this statement
        iswap = 0
        do m = nelt, nst + 1, -1
            iso = nconf(m)
            do n = nst, m - 1
                if (iso < nconf(n)) iswap = iswap + 1
            end do
        end do

        ! phase given by whether number of swaps is odd or even
        iphase = merge(1, -1, mod(iswap, 2) == 0)

    end function iphase


    !> \brief Compute the orbital table.
    !>
    !> Computes the orbital table which is then used in the projection step. This is called from projec().
    !> \verbatim
    !>    Input data:
    !>       ISYMTYP  Switch for C-inf-v (=0 or 1) / Abelian point group (=2
    !>           NOB  Number of orbitals per symmetry
    !>          NSYM  Number of symmetries in the orbital set
    !>         NPFLG  Flag controlling printing of computed orbital table
    !>
    !>    Output data:
    !>             MN  Orbital number associated with each spin-orbital
    !>             MG  G/U designation for each spin-orbital (C-inf-v only)
    !>                 Actually this is always zero because C-inf-v does not
    !>                 distinguish between g/u. It exists because original
    !>                 version of Alchemy tried to use it for D-inf-h too;
    !>                 all CI evauation is doen in C-inf-v now because CONGE
    !>                 converts D-inf-h to C-inf-v data.
    !>             MM  Symmetry quantum number associated with each spin-orb
    !>             MS  Spin function ( alpha or beta ) associated with each
    !>                 spin orbital
    !>
    !>    Notes:
    !>
    !>     The orbital table establishes orbital and quantum number data for
    !>    each spin orbital in the set.
    !>
    !>    e.g. C-inf-v symmetry with NSYM=2, NOB=3,1, yields ten spin
    !>         orbitals which are designated as follows by this routine:
    !>
    !>       Spin orb.     MN  MG  MM  MS     Comments
    !>           1          1   0   0   0     1 sigma spin up
    !>           2          1   0   0   1     1 sigma spin down
    !>           3          2   0   0   0     2 sigma spin up
    !>           4          2   0   0   1     2 sigma spin down
    !>           5          3   0   0   0     3 sigma spin up
    !>           6          3   0   0   1     3 sigma spin down
    !>           7          4   0   1   0     1 pi(lambda=+1) spin up
    !>           8          4   0   1   1     1 pi(lambda=+1) spin down
    !>           9          4   0  -1   0     1 pi(lambda=-1) spin up
    !>          10          4   0  -1   1     1 pi(lambda=-1) spin down
    !> \endverbatim
    !> \note MAL 11/05/2011 : Changes made here are to bring the subroutine into
    !>       line with the changes that were made in 'projec' in order to utilize
    !>       dynamic memory and also to comply with the F95 standard
    !>
    subroutine mkorbs (nob, nsym, mn, mg, mm, ms, norb, nsrb_in, map, mpos, iposit, nobl, nob0l, symtyp)

        use precisn,     only : wp
        use congen_data, only : nftw

        integer :: nsym
        integer :: iposit
        integer :: symtyp
        integer :: norb
        integer :: nsrb_in
        integer :: nob(nsym)
        integer :: mn(nsrb_in)
        integer :: mg(nsrb_in)
        integer :: mm(nsrb_in)
        integer :: ms(nsrb_in)
        integer :: map(norb)
        integer :: mpos(nsrb_in)
        integer :: nobl(*)
        integer :: nob0l(nsym)

        integer  i, ik, ikp, ipos, is, ic, iso
        integer  j, k
        integer  m, ma, mb, m1, n, nep
        integer  ierr
        integer  len_noblj

        integer :: nsrb
        integer, allocatable :: noblj(:)
        integer, parameter :: iwrite = 6
        logical, parameter :: zdebug = .false.

        ! Debug banner header

        if (zdebug) then
            write(nftw,'(/,10x,"====> mkorbs() <====",/)')
            write(nftw,'(/,10x,"Input data: "     )')
            write(nftw,'(  10x,"  nsym    = ",I6  )') nsym
            write(nftw,'(  10x,"  symtyp  = ",I6  )') symtyp
            write(nftw,'(  10x,"  iposit  = ",I6  )') iposit
            write(nftw,'(  10x,"  norbs   = ",I6  )') norb
            write(nftw,'(  10x,"  nsrb_in = ",I6,/)') nsrb_in
            write(nftw,'(12x,"nob:  ",1x,20(I3,1x))') (nob(i), i=1,nsym)
            if (symtyp == 1) then
                write(nftw,'(12x,"nobl: ",1x,20(I3,1x))') (nobl(i), i=1,nsym)
            else
                write(nftw,'(12x,"nobl: ",1x,20(I3,1x))') (nobl(i), i=1,2*nsym)
            end if
            write(nftw,'(/,10x,"**** End of input data",/)')
        end if

        ! Copy the contents of input array nobl() to local storage noblj()
        !  - When we are working with D-inf-h we need to remember that
        !    nsym and nob refer to the C-inf-v representation, whereas
        !    nobl() holds the D-inf-h representation. Thus we need to
        !    double "nsym".

        if (symtyp == 1) then
            len_noblj = 2 * nsym
        else
            len_noblj = nsym
        end if

        allocate(noblj(len_noblj), stat = ierr)

        if (ierr /= 0) then
            write(nftw,'(/,10x,"**** Error in mkorbs() ",/)')
            write(nftw,'(/,10x,"Cannot allocate noblj - ierr = ",I6,/)') ierr
            stop
        end if

        noblj(1:len_noblj) = nobl(1:len_noblj)

        if (iposit /= 0) then
            do is = 1, nsym
                noblj(is) = nobl(is) / 2
            end do 
        end if

        !======================================================================
        !
        !     E L E C T R O N I C    O R B I T A L S
        !
        !======================================================================

        if (symtyp == 0) then
            ic = 1
            iso = 4
        else if (symtyp == 1) then
            ic = 2
            iso = 4
        else
            ic = 1
            iso = 2
        end if
        if (iposit /= 0) then
            do is = 1, nsym*ic
                noblj(is) = nobl(is) / 2
            end do
        else
            do is = 1, nsym*ic
                noblj(is) = nobl(is)
            end do
        end if

        ! First of all we loop over all non-degenerate electron orbitals and build the table of spin-orbitals for them.
        ! We set mpos() to be zero for "electron" orbitals.
        ! For linear molecules this is
        !    C-inf-v : Sigma type
        !    D-inf-h : Sigma_g and Sigma_u
        ! Actually the code also handle here the first IRR of Abelian point groups too.

        i  = 1
        ma = 0

        do j = 1, ic
            m1   = ma + 1
            ma   = ma + noblj(j)
            ipos = 0

            do n = m1, ma
                map(n) = n

                ! Spin orbital with spin-up

                mn(i)   = n
                mg(i)   = 0
                mm(i)   = 0
                ms(i)   = 0
                mpos(i) = ipos

                ! Spin orbital with spin-down

                i = i + 1

                mn(i)   = n
                mg(i)   = 0
                mm(i)   = 0
                ms(i)   = 1
                mpos(i) = ipos

                i = i + 1
            end do
        end do

        ! Process remaining orbitals
        !       C-inf-v:  Pi, Delta, ....
        !       D-inf-h:  Pi_u, Pi_g, Delta_g, ...
        !       Abelian:  irr = 2, 3, 4,

        k = ma + 1

        do m = ic + 1, nsym * ic
            ma = noblj(m)
            mb = (m - 1) / ic
            ipos = 0
            do n = 1, ma
                map(k) = k
                do j = 1, iso
                    mn(i) = k
                    mg(i) = 0
                    mm(i) = mb
                    ms(i) = 0
                    mpos(i) = ipos
                    i = i + 1
                end do
                k = k + 1
                ms(i - 1) = 1
                if (symtyp <= 1) then
                    mm(i - 1) = -mb
                    mm(i - 2) = -mb
                    ms(i - 3) = 1
                end if
            end do
        end do

        ! Compute the total number of electron type spin orbitals.

        nsrb = i - 1

        !======================================================================
        !
        !     P O S I T R O N I C    O R B I T A L S
        !
        !======================================================================

        if (iposit /= 0) then
            write(nftw,*) ' MAP : OLD, NEW'
            ik = k - 1
            ikp = ik + 1
            nep = ik + ik   ! total number of orbitals (electron + positron)
            do n = ikp, nep
                map(n) = map(n - ik)
            end do
            write(nftw,*) ' MAP : OLD, NEW'
            do n = 1, nep
                write(nftw,*) n,map(n)
            end do
            DO N = IKP, NEP
                MN(IK + IKP) = N
                IK = IK + 1
                MN(IK + IKP) = N
                IK = IK + 1
            end do
            do j=1, nep
                mg(j+nep)=mg(j)
                mm(j+nep)=mm(j)
                ms(j+nep)=ms(j)
                mpos(j+nep)=1
            end do
            nsrb = 2 * nsrb
        end if


        ! Now print-out a table of the spin-orbitals

        write(nftw,*) '                  I         N         G         M         S      MPOS      '

        do j = 1, nsrb
            write(iwrite,'(10x,7I10)') J, MN(J), MG(J), MM(J), MS(J), MPOS(J)
        end do

        if (nsrb /= nsrb_in) then
            write(iwrite,'(" HELP!!! MKORBS: NSRB, NSRBD = ",2I6)') nsrb, nsrb_in
            stop 999
        end if


        ! Return point
        !  - Release any allocated storage

        if (allocated(noblj)) then
            deallocate(noblj, stat = ierr)

            if (0 /= ierr) then
                write(iwrite,'(/,10x,"**** Error in mkorbs() ",/)')
                write(iwrite,'(/,10x,"Cannot deallocate noblj - ierr = ",I6,/)') ierr
                stop 999
            end if
        end if

        if (zdebug) then
            write(iwrite,'(/,10x,"***** Completed - mkorbs() ",/)')
        end if

    end subroutine mkorbs


    !> \brief Pack wave function.
    !>
    !> Reformats (packs) the CSF expression into the style used throughout the rest of Alchemy, that is as a set of
    !> replacements from the reference determinant. Adds this to the end of the array ntdo() from location "ndo".
    !>
    !> On entry to this routine we have the CSF defined for us as follows (from the projection step):
    !>  1. there are "nod" determinants 
    !>  2. each determinant is of length "ieltp". 
    !>  3. "cdo" contains the coefficient which multiplies each determinant. This is derived from the coupling process.
    !>  4. the determinants are stored in "mdo", as a list of spin orbitals - so it is of lenth nod*ieltp
    !>
    !> The above information is complemented by the analysis in the calling routine which classifies spin orbitals in this CSF 
    !> wrt the reference determinant:                      
    !>
    !>  5. "idopl"  is the number of spin orbs in the reference det but not present in this CSF;
    !>     "mdop()" is the list of those spin orbitals 
    !>  6. "idcpl"  is the number of spin orbs in this CSF but not present in the reference determinant;
    !>     "mdcp()" is the list of those spin orbitals.
    !>
    !> So, given all of the above information, the determinants in "mdo" are processed and each expressed in the format
    !>  - number of replacements from ref determinant
    !>  - list of replaced spin orbitals
    !>  - list of replacing spin orbitals
    !>
    !> The output is placed into array "ndto". The length available in "ndto" is passed into the routing in "len_ndto" and this is
    !> monitored to be sure we do not overflow it.
    !>
    !> During the process, it may be necessary to multiply "cdo" by -1 as we order spin orbitals. cdo() holds the coefficients 
    !> associated with each determinant. These were constructed earlier in the spin projection - note we only receive the relevant 
    !> "piece" of the cdo(0 array in the arglist, the bit for this CSF, not all of it for all CSFs, as is the case with "ndto()".
    !> "nftw" is the logical unit for the printer                      
    !>
    !> \note MAL 11/05/2011: Changes made here are to bring the subroutine into
    !>                       line with the changes that were made in 'projec' in order to utilize
    !>                       dynamic memory and to comply with the F95 standard.
    !>
    subroutine pkwf (nod, ieltp, cdo, mdo, idopl, mdop, idcpl, mdcp, nftw, ndo, ndto, len_ndto, ithis_csf)              

        use precisn, only : wp

        integer nod                 ! number of determinants
        integer ieltp               ! number of electrons in each determinant
        integer idopl, mdop(idopl)  ! S.O.s in refdet but not in this CSF.
        integer idcpl, mdcp(idcpl)  ! S.O.s in this CSF but not ref det.
        integer mdo(nod*ieltp)      ! the determinants on input 
        integer nftw                ! logical unit for the printer 
        integer len_ndto, ndto(len_ndto) ! the determinants on output 
        integer ndo                 ! On output points to the highest location used in array ndto
        integer ithis_csf           ! The present CSF index - helpful in error msgs.
    
        real(kind=wp)  :: cdo(nod)  ! coupling coeffs for each deterinant

        integer  i, j, k, n, nc, nd, md, mdopi
        integer  mdi(idcpl+ieltp)   ! Temporary workspace array
        integer  nt(idopl)          ! Temporary workspace array

        logical, parameter :: zdebug = .false.

        ! Debug banner header 
        if (zdebug) then
            write(nftw,'(/,10x,"====> PKWF() <====",/)')
            write(nftw,'(  10x,"Input data: ")')
            write(nftw,'(  10x,"  No. of determinants       (nod) = ",I5)') nod
            write(nftw,'(  10x,"  No of electrons per det (ieltp) = ",I5)') ieltp
            write(nftw,'(  10x,"  Input determinants: ")')
            md = 0
            do i = 1, nod
                write(nftw,'(/,10x,"  Determinant ",I5," Coeffcient = ",F13.6,/)') i, cdo(i)
                write(nftw,'(  10x,"  Spin orbs: ",20(I3,1x),/,(25x,20(I3,1x)))') (mdo(md+j), j=1,ieltp)
                md = md + ieltp
            end do
            write(nftw,'(/,10x,"  No. spin orbs in the reference det ")')
            write(nftw,'(  10x,"  but not present in this CSF (idopl) = ",I5,/)') idopl
            write(nftw,'(  10x,"  mdop: ",10(I3,1x),/,(16x,10(i3,1x)))') (mdop(i), i=1,idopl)
            write(nftw,'(/,10x,"  No. spin orbs in this CSF but ")')
            write(nftw,'(  10x,"  not present in ref det (idcpl) = ",I5,/)') idcpl
            write(nftw,'(  10x,"  mdcp: ",10(I3,1x),/,(16x,10(I3,1x)))') (mdcp(i), i=1,idcpl)
            write(nftw,'(/,10x,"Space available in ndto() = ",I10,/)') len_ndto              
        end if

        ! Special case is where idopl is 0. There is no work to be done.
        ! We have same spin orbitals as the reference determinant. 

        if (idopl == 0) then

            ndto(ndo) = 0
            ndo       = ndo + 1

        else

            ! Loop over all determinants in this CSF
            !
            !     "md" points at the location of the input determinant as we work though the list. Remember that each is 
            !          of length ieltp.

            md = 0

            do k = 1, nod

                ! Populate mdi() with the list of spin orbitals in this CSF but not in the reference determinant. 
                mdi(1:idcpl) = mdcp(1:idcpl)

                nd = idcpl

                ! Copy the current determinant onto the end of array mdi()
                do i = 1, ieltp
                    nd      = nd + 1
                    mdi(nd) = mdo(md + i)
                end do

                ! Extract the list of spin orbitals which have replaced those in the reference determinant but not present in 
                ! this CSF - as defined by array MDOP().  
                !
                ! We build "nt()" are a set pointers into mdi(). 

                nd = 0

                outer_loop: do i = 1, idopl
                    mdopi = mdop(i)
                    inner_loop: do j = 1, idopl
                        if (mdi(j) == mdopi) then
                            if (i /= j) then
                                cdo(k) = -cdo(k)
                                mdi(j) = mdi(i)
                            end if
                            mdi(i) = 0
                            cycle outer_loop
                        end if
                    end do inner_loop
                    nd = nd + 1
                    nt(nd) = i
                end do outer_loop

                ! Ok, we now know how long the "packed" determinant is.
                ! Let's check that we have enough space available in NDTO in which to store the data.
                !
                ! This data is: 
                !
                !     number of replacements
                !     list of replaced spin orbs
                !     list of replacing spin orbs

                if (ndo + 2 * nd > len_ndto) then
                    write(nftw,'(/,10x,"***** Error in: PKWF() ",/)')
                    write(nftw,'(  10x,"There is not enough space in NDTO to store the ")')
                    write(nftw,'(  10x,"present determinant (",I4," of ",I4," ). ")') k, nod
                    write(nftw,'(  10x,"Space needed = ",I8," Given (len_ndto) = ",I8)') ndo+2*nd, len_ndto
                    write(nftw,'(  10x,"This present CSF number = ",I10,/)') ithis_csf
                    stop 999
                end if

                ! Copy the determinant into place in NDTO
                ! It is useful to remember that "nd" is the number of replacements from the reference determinant

                ndto(ndo) = nd
                nc        = ndo + nd

                do i = 1, nd
                    n = nt(i)
                    ndto(ndo+i) = mdop(n)
                    ndto(nc+i)  = mdi(n)
                end do

                if (zdebug) then
                    write(nftw,'(/,10x,"Packed format for determinant ",I5,": ",/)') k
                    write(nftw,'(  10x,20(I3,1x))') (ndto(i), i =ndo,ndo+2*nd)
                end if

                ndo = ndo + nd + nd + 1

                ! Augment the pointer for the next determinant
                md = md + ieltp

            end do ! End of loop over determinants in this CSF

        end if

        ! Subroutine return point

        if (zdebug) then
            write(nftw,'(/,10x,"On output: ")')
            write(nftw,'(  10x,"   Highest location in ndto()  (ndo) = ",i10,/)') ndo
            write(nftw,'(/,10x,"**** PKWF() - completed",/)')
        end if

    end subroutine pkwf


    !> \brief ?
    !>
    subroutine pmkorbs (nob, nobe, nsym, mn, mg, mm, ms, nsrb, norb, nsrbd, map, mpos, iposit, symtyp)

        use congen_data, only : nftw

        integer :: symtyp, maxspin, imo, emo, ispin, iso, ipos, isym, j, jmo, maxmo, minmo, n, amo, nsrbd, &
                   nsym, iposit, norb, nsrb
        integer :: nob(nsym), nobe(nsym), mpos(nsrb), map(norb)
        integer :: mn(nsrb), mg(nsrb), mm(nsrb), ms(nsrb)

        !     Setting up the following arrays:
        !
        !       mn()   = orbital number
        !       mg()   = ??? distinguish gerade and ungerade
        !       mm()   = ??? m-quantum number for degenerate MOs
        !       ms()   = spin (for alpha=0, for beta=1)
        !       mpos() = flag for positron (for e-=0, for p+=1)
        !
        !     NOTE: only implemented for poly-atomic code (symtyp=2)  
        
        if (symtyp /= 2) then
            write(nftw,*) ' ERROR in PMKORBS: calculation with positrons'
            write(nftw,*) '                   only possible for SYMTYP=2'
            write(nftw,*) '                   (abelian groups).'
            write(nftw,*) ' here: SYMTYP=',SYMTYP         
            stop
        end if

        maxspin = 2

        !     imo = mo-number
        !     iso = so-number
        !     ipos = positron-flag 
        !          = 0 for 1..nobe(isym)
        !          = 1 for nobe(isym)..nob(isym)

        imo = 0
        iso = 0
        emo = 0

        do isym = 1, nsym

            ! electronic MOs
            ipos = 0
            maxmo = nobe(isym)
            amo = emo
            do jmo = 1, maxmo
                imo = imo + 1
                emo = emo + 1
                map(imo) = emo
                do ispin = 1, maxspin
                    iso = iso + 1
                    mn(iso) = imo
                    mg(iso) = 0
                    mm(iso) = isym - 1
                    ms(iso) = ispin - 1
                    mpos(iso) = ipos
                end do
            end do

            ! positronic MOs
            ipos = 1
            minmo = nobe(isym) + 1
            maxmo = nob(isym)
            !shift=ipos*nobe(isym)
            emo = amo
            do jmo = minmo, maxmo
                imo = imo + 1
                emo = emo + 1
                !map(imo) = imo - shift
                map(imo) = emo
                do ispin = 1, maxspin
                    iso = iso + 1
                    mn(iso) = imo
                    mg(iso) = 0
                    mm(iso) = isym - 1
                    ms(iso) = ispin - 1
                    mpos(iso) = ipos
                end do
            end do

        end do

        nsrb = iso

        ! output the labels
        write(nftw,*) ' MAP : OLD, NEW'
        do n = 1, norb
            write(nftw,*) n, map(n)
        end do

        write(6,*) '                  I         N         G         M         S      MPOS    '

        do j = 1, nsrb
            write(6,'(10X,7I10)') j, mn(j), mg(j), mm(j), ms(j), mpos(j)
        end do

        ! Control number of spin orbitals
        write(6,*) 'GIVEN NSRBD=', nsrbd
        write(6,*) 'CALCULATED NSRB=', nsrb

        if (nsrb /= nsrbd) then
            write(6,'(" HELP!!! MKORBS: NSRB, NSRBD = ",2I6)') nsrb, nsrbd
            stop
        end if

    end subroutine pmkorbs


    !> \brief Get open-shell part of determinants.
    !>
    !> Fills the array \c mop with subset of specification of every determinant for the current CSF. Only spin-orbitals
    !> that are in open shells are used, as the rest does not participate in spin composition done later in \ref prjct.
    !> The written spin-orbitals are also immediately sorted in non-descending order. The even or odd number of swaps needed
    !> for sorting is returned (as .false. and .true., respectively) via the logical array \c flip, so that the sign of
    !> determinant coefficients within this CSF can be adjusted to the used order.
    !>
    !> \verbatim
    !>     OUTPUT IDOP        NO OF SO IN DR BUT NOT IN DC
    !>            MDOP(NELT)        SO IN DR BUT NOT IN DC
    !>            IDCP        NO OF SO IN DC BUT NOT IN DR
    !>            MDCP(NELT)        SO IN DC BUT NOT IN DR
    !>            IELTP       NO OF SO IN OPEN SHELLS FOR A DTR
    !>            MOP(MOPMX)  SO
    !> \endverbatim
    !>
    !> \note In the original code there was a common block
    !>       \verbatim
    !>          /OWF/ IDOP,IDCP,IELTP
    !>       \endverbatim
    !>       which was used to pass three integer values back to the 
    !>       caling routine. these have now been placed into the
    !>       argument list; correspondingly the routine WFGNTR has
    !>       been modified. It is the only routine whihc calls this 
    !>       one. the purpose of the variables is as follows:
    !>       \verbatim
    !>            IDOP  holds the number of entries in MDOP 
    !>            IDCP  holds the number of entries in MDCP
    !>            IELTP is the number of electrons in open shell
    !>       \endverbatim
    !>
    subroutine popnwf (nsrb, nsrbs, nelt, ndtrf, mopmx, mdop, mdcp, mop, mdc, mdo, ndta, nod, nda, idop, idcp, ieltp, flip, nalm)

        use precisn, only : wp

        integer nsrb
        integer nsrbs
        integer nelt
        integer nod   ! Number of replacements from reference
        integer nalm  ! Output - return code
                      !   =0,    normal exit
                      !   =1,    different neltp
                      !   =2,    need more space for mop
                      !   =3,    neltp=0, but nod not =1

        integer :: ndtrf(nelt)  ! Input: the reference determinant
        integer :: mdc(nsrb)    ! Workspace: spin orbs in closed shell, common to all determinants 
        integer :: mdo(nsrb)    ! Workspace: Union of all spin-orbs in open shell  
        integer :: ndta(nsrb)   ! Workspace: to expand to full determinant (how many times each spin-orbital is used in CSF)
        integer :: mopmx        ! Maximum size of the mop() array
        integer :: mop(mopmx)   ! Spin-orbs 
        integer :: mdop(nelt)
        integer :: mdcp(nelt)
        logical :: flip(nod)    ! Sign-correction factor for each determinant
        integer :: nda(*)       ! Input - the array of packed determinants on which we work

        integer idop, idcp, ieltp

        integer i, k, m, md, me, n, na, nb, ndo, ndc, nod2, nd, no, no0
        integer ndop, ndcp, neltp

        integer, parameter :: nftw = 6 ! logical unit for printer

        logical, parameter :: zdebug = .false.

        ! Debug banner header
        if (zdebug) then
            write(nftw,'(/,25x,"====> POPNWF() <====",/)')
            write(nftw,'(/,25x,"Input data: ")')
            write(nftw,'(  25x,"  No. of spin orbitals          (nsrb) = ",I10)') nsrb
            write(nftw,'(  25x,"  No .of sigma-type spin orbs  (nsrbs) = ",I10)') nsrbs
            write(nftw,'(  25x,"  No. of electrons              (nelt) = ",I10)') nelt
            write(nftw,'(  25x,"  Units available in mop()     (mopmx) = ",I10)') mopmx
            write(nftw,'(  25x,"  No. of dets in this CSF        (nod) = ",I10,/)') nod
            write(nftw,'(  25x,"  Ref det = ",10(I5,1x),/,(37x,10(I5,1x)))') (ndtrf(i), i=1,nelt)
        end if

        ! Initialize the return data 
        idop  = 0 
        idcp  = 0
        ieltp = 0
        nalm  = 0

        ! Initialize local data 
        ndop  = 0 
        ndcp  = 0
        neltp = 0

        !======================================================================
        !
        !    S T E P :  1
        !
        !======================================================================

        ! We build NDTA() which repesents the spin-orbitals used in this CSF. 
        ! We start by initializing NDTA which is of length equal to the number of spin-orbitals in the system. 

        ndta(1:nsrb) = 0

        ! Now loop over the reference determinant and for every spin-orbital within it, we mark that spin orbital
        ! to be populated in EVERY determinant of this CSF. Thus we set its value in NDTA to be equal to the number
        ! of determinants in this CSF. 

        ndta(ndtrf(1:nelt)) = nod

        ! Loop over all determinants in this CSF and modify the count per spin-orb in ndta()

        md = 1  ! index in nda, where information about current determinant starts

        do i = 1, nod
            m = nda(md)  ! number of replacements (w.r.t. reference det) defining this determinant

            ! Loop over all replacement/replacing spin-orbitals
            !
            ! For each replaced spin-orb, we decrement its count in ndta() and for each replacing spin-orb, we increment 
            ! its count in ndta.

            do k = md + 1, md + m
                ndta( nda(k) )   = ndta( nda(k)   ) - 1     ! decrement use of "replaced" spin-orbital
                ndta( nda(k+m) ) = ndta( nda(k+m) ) + 1     ! increment use of "replacing" spin-orbital
            end do

            ! Update the pointer "md" to be at the start of the  next determinant. That is the value which defines
            ! the number of replacements in that determinant. 

            md = md + 2*m + 1  ! move to the next packed determinant
        end do

        if (zdebug) then
            write(nftw,'(/,25x,"Expanded determinant (NDTA) representation after")')
            write(nftw,'(  25x,"processing all (",I6,") dets within the ")') nod
            write(nftw,'(  25x,"present CSF.                            ",/)')
            write(nftw,'(  25x,"Spin Orb.   Count ")')
            write(nftw,'(  25x,"---------  -------")')
            write(nftw,'(  (25x,I9,2x,I7))') (i, ndta(i), i=1,nsrb)
        end if

        !==========================================================================
        !
        !    S T E P :  2
        !
        !==========================================================================

        ! We loop over all "orbitals" which are not lambda-degenerate
        !
        ! This means all sigma type orbitals when dealing with C_inf_v and ALL orbitals when dealing with D2h and sub-groups.
        !
        ! Given we have set the occupancy of an occupied spin-orbital to "nod" initially, then since each orbital is composed of 
        ! TWO spin orbitals, we will know that the orbital is full if its occupancy is 2*nod. It may or may have been processed
        ! in the list of determinants. Remember that NDTA() is summed over ALL determinants in the CSF. 
        !
        !     Following code: 
        !
        !         1. Examines each orbital 
        !         2. If an orbital is FULL, we store the constituent spin-orbs in MDC()
        !         3. Otherwise we store them in MDO() 
        !
        !     Note: this code works at the orbital level but produces output at the spin-orbital level

        nod2 = nod + nod

        ndo = 0
        ndc = 0

        do i = 1, nsrbs, 2
            if (ndta(i) + ndta(i+1) == nod2) then
                mdc(ndc+1) = i             ! Fully occupied
                mdc(ndc+2) = i + 1
                ndc        = ndc + 2
            else
                if (ndta(i) /= 0) then     ! Partially occupied
                    ndo      = ndo + 1
                    mdo(ndo) = i
                end if

                if (ndta(i+1) /= 0) then
                    ndo      = ndo + 1
                    mdo(ndo) = i + 1
                end if
            end if
        end do

        ! We do the same thing again but now for lambda-degenerate orbitals - if any exist.
        ! Of course the occupany is 4*nod.

        nod2 = nod2 + nod2

        do i = nsrbs + 1, nsrb, 4
            nd = ndta(I) + ndta(i+1) + ndta(i+2) + ndta(i+3)

            if (nd == nod2) then
                do k = i, i + 3
                    ndc      = ndc + 1
                    mdc(ndc) = k
                end do
            else
                do k = i, i + 3
                    if (ndta(k) /= 0) then
                        ndo      = ndo + 1
                        mdo(ndo) = k
                    end if
                end do
            end if
        end do

        if (zdebug) then
            write(nftw,'(/,25x,"After step 2 we have; ",/)')
            write(nftw,'(  25x,"  Number of closed orbtials (ndc) = ",I6)') ndc
            write(nftw,'(  25x,"  Number of open   orbitals (ndo) = ",I6,/)') ndo
            write(nftw,'(  25x,"Closed orbitals: ",20(I4,1x),/,(20x,20(I4,1x)))') (mdc(i), i=1,ndc)
            write(nftw,'(  25x,"Open   orbitals: ",20(I4,1x),/,(20x,20(I4,1x)))') (mdo(i), i=1,ndo)
        end if

        !==========================================================================
        !
        !    S T E P :  3
        !
        !==========================================================================

        ! ndta() is an array with one entry for each spin-orbitals

        ! First we zeroize the array + mark any spin-orbital in the reference determinant as being occuiped in ndta()
        ndta(1:nsrb) = 0
        ndta(ndtrf(1:nelt)) = 1

        if (ndc == 0) then
            mdop(1:nelt) = ndtrf(1:nelt)
            ndop = nelt
            ndcp = 0
        else
            do i = 1, ndc
                n = mdc(i)
                ndta(n) = ndta(n) + 1
            end do

            ndop = 0
            do i = 1, nelt
                n = ndtrf(i)
                if (ndta(n) == 2) cycle
                ndop = ndop + 1
                mdop(ndop) = n
            end do

            ndcp = 0
            do i = 1, ndc
                n = mdc(i)
                if (ndta(n) == 2) cycle
                ndcp = ndcp + 1
                mdcp(ndcp) = n
            end do
        end if

        neltp = nelt - ndc

        if (neltp /= 0) then

            ! Yet another use of ndta(), now: let it contain spin-orbitals from reference determinant
            ! and all occupied spin-orbitals in open shells.

            ndta(1:nsrb) = 0            ! clear the array
            ndta(ndtrf(1:nelt)) = 1     ! put 1 to every spin-orbital in reference determinant

            do i = 1, ndo               ! loop over all occupied spin-orbitals in open shells
                n = mdo(i)              ! index of the occupied spin-orbital
                ndta(n) = ndta(n) + 1   ! increment use of that spin-orbital
                mdc(i) = ndta(n)        ! backup the value of ndta(n)
            end do

            ! Loop over all determinants in the CSF, populating the array "mop" with subset of the specification of the determinants:
            ! only spin-orbitals in open shells are used here, as the rest (spin-orbitals in closed shells) is not interesting for
            ! the projection algorithm.

            no = 0
            md = 1                                          ! index in "mda" where the current determinant information starts

            do i = 1, nod                                   ! loop over all determinants
                m = nda(md)                                 ! number of replacements (wrt reference det) defining the current det
                me = md + m                                 ! index in "mda" where the to-be-replaced orbitals list ends

                do k = 1, m                                 ! loop over all replacements defining this determinant
                    na = nda(md + k)                        ! index of spin-orbital to be replaced ...
                    nb = nda(me + k)                        ! ... by this spin-orbital
                    ndta(na) = ndta(na) - 1                 ! decrease use of replaced spin-orbital
                    ndta(nb) = ndta(nb) + 1                 ! increase use of replacing spin-orbital
                end do

                no0 = no                                    ! remember size of "mop" before addition of spin-orbitals to it

                do k = 1, ndo                               ! loop over all occupied spin-orbitals in partially occupied orbitals
                    n = mdo(k)                              ! label of the occupied spin-orbital

                    if (ndta(n) == 2) then 

                        ! So, either already the reference determinant already uses a spin-orbital from open shell (and the current
                        ! determinant does not replace it), or the current determinat uses a spin-orbital from open shell as a replacement
                        ! of something else. In any case, this determinant ends with using the open-shell spin-orbital.

                        no = no + 1                         ! increment number of used open-shell spin-orbitals

                        if (no > mopmx) then                ! too much data - limits too tight
                            write(nftw,'(/,25x,"**** Error in; POPNWF() ",/)')
                            write(nftw,'(/,25x,"Exceeded size of mop() ",/)')
                            write(nftw,'(  25x,"   Determinant num (i) = ",I6)') i
                            write(nftw,'(  25x,"   spin orbital   (no) = ",I6)') no
                            write(nftw,'(  25x,"   mopmx               = ",I6,/)') mopmx
                            stop 999
                        end if

                        mop(no) = n                         ! store the open-shell spin-orbital label for further processing
                    end if

                    ndta(n) = mdc(k)                        ! restore "reference + open shells" spin-orbital use count ndta(n) from backup
                end do

                md = md + m + m + 1                         ! move on to the beginning of the next packed determinant

                ! sort the spin-orbitals in the open-shells-only subset of the determinant and store the corresponding permutation sign
                flip(i) = mod(qsort(no-no0, mop(no0+1:no)), 2) /= 0
            end do

            ! verify that each determinant contributed the same number of spin-orbitals from open shells
            if (mod(no, nod) /= 0) then
                write(nftw,'(/,25x,"**** Error in; POPNWF() ",/)')
                stop 999
            end if

        end if

        ! Finalize

        if (neltp == 0 .and. nod /= 1) then

            nalm = 3

        else

            ! We reach this point if all has gone successfully.
            ! We copy the work variables into the return variables and then set a return code of success (nalm=0)

            idop = ndop
            idcp = ndcp
            ieltp = neltp

            nalm = 0

        end if

        ! Subroutine return point

        if (zdebug) then
            write(nftw,'(/,25x,"Output data: ")')
            write(nftw,'(  25x,"                                (idop) = ",I10)') idop
            write(nftw,'(  25x,"                                (idcp) = ",I10)') idcp
            write(nftw,'(  25x,"  No. electrons in open shells (ieltp) = ",I10,/)') ieltp

            write(nftw,'(  27x,"Spin orbitals in DR but not DC: ",/)') 
            write(nftw,'(  27x,9(I4,1x))') (mdop(i), i=1,idop)
            write(nftw,'(/,27x,"Spin orbitals in DC but not DR: ",/)') 
            write(nftw,'(  27x,9(I4,1x))') (mdcp(i), i=1,idcp)

            write(nftw,'(/,25x,"Return code (nalm) = ",I10,/)') nalm
            write(nftw,'(/,25x,"**** Completed - POPNWF() ",/)')
        end if

    end subroutine popnwf


    !> \brief Apply Lowdin projection operator
    !>
    !> This routine applies the Lowdin projection operator. More details
    !> can be found in the literature at for example:
    !> \verbatim
    !>       Nelson F Beebe and Sten Lucil, J Phys B: At Mol Phys, Vol. 8, Issue 14, 1975, p2320
    !> \endverbatim
    !> This routine is called when a CSF is found to have two, or more electrons in open shells. Each pair of spin-orbitals in
    !> each determinant is examined and potentially used to create a  new determinant. Thus the output expression for the CSF 
    !> \verbatim
    !>       nodo, cdo(), ndo()
    !> \endverbatim
    !> may be much larger than the input 
    !> \verbatim
    !>       nodi, cdi(), ndo()                       
    !> \endverbatim
    !> Note reuse of ndo() here. cdi() is used an an extendable buffer too and must have more than "nodi" elements.
    !>
    !> The generated list of determinants is examined for any with very small coefficients (thres) and these are removed. Thus the
    !> the list may grow and shrink in this routine.                    
    !>
    !> Of course the nature of the projection process is controlled by the quantum numbers input.
    !>
    !> The routine terminates with an error message if any error conditions are found.
    !>
    subroutine prjct (nelt, mxss, nodi, ndo, cdi, nodo, cdo, maxcdo, mgvn, iss, isd, thres, r, ndtr, mm, ms, maxndo, symtyp, nsrb)

        use precisn,       only : wp
        use congen_bstree, only : det_tree

        real(kind=wp)              :: fcta,fctb,fctc,fctr,tmp
        real(kind=wp), parameter   :: zero = 0.0_wp      
        real(kind=wp), parameter   :: one  = 1.0_wp      
        real(kind=wp), parameter   :: four = 4.0_wp      

        integer  :: nelt    ! # of electrons in open shells (i.e. per det)
        integer  :: mxss    ! max S for projecttion operator
        integer  :: nodi    ! Number of determinants in CSF on input 
        integer  :: nodo    ! Number of determinants in CSF in output
        integer  :: maxcdo  ! Max size of cdo()
        integer  :: maxndo  ! Max size of ndo()
        integer  :: mgvn    ! Lambda or IRR for CSF
        integer  :: iss     ! Required S value
        integer  :: isd     ! Sz for determinants

        real(kind=wp)  :: r           ! +1.0 for Sigma(+), -1.0 for Sigma(-)
        real(kind=wp)  :: thres       ! Determinants with coefficients < thres are deleted
        real(kind=wp)  :: cdi(*)      ! Expansion coeffs of each det on input - grows as the routine adds new determinants to the CSF
        real(kind=wp)  :: cdo(maxcdo) ! On output holds the expansion coefficients for all determinants

        integer  :: symtyp      ! Designates C-inf-v, D-inf-h or Abelian
        integer  :: nsrb        ! Number of spin-orbitals in system

        integer, target :: ndo(maxndo) ! Input determinants - overwitten
        integer  :: ndtr(nsrb)  ! Workspace for building new determinants
        integer  :: mm(nsrb)    ! Lambda/IRR for each spin orbital
        integer  :: ms(nsrb)    ! Spin (Sz) of each spin prbital

        integer  :: i, ia, ib, id, idet, is, issp, istart, ma, mb, mga, mgb, nd, ninitial_dets 

        integer, pointer           :: ndo_ptr(:)
        logical, parameter         :: zdebug = .false.
        logical                    :: flip

        type(det_tree) :: bst     ! Binary search tree structure on top of the determinant storage with lexicographical compare

        ! Banner header
        if (zdebug) then
            write(6,'(/,20x,"====> prjct() - project wavefunction <====",/)')
            write(6,'(  20x,"Input data: ")')
            write(6,'(  20x,"  # electrons  open shell  (nelt) = ",I7)') nelt
            write(6,'(  20x,"  Maximum S for projection (mxss) = ",I7)') mxss
            write(6,'(  20x,"  Lambda value for Wavefn  (mgvn) = ",I7)') mgvn
            write(6,'(  20x,"  Required Spin value       (iss) = ",I7)') iss
            write(6,'(  20x,"  Required Sz   value       (isd) = ",I7)') isd
            write(6,'(  20x,"  Threshold               (thres) = ",D13.6)') thres
            write(6,'(  20x,"  Dimension of cdo()     (maxcdo) = ",I7)') maxcdo
            write(6,'(  20x,"  Dimension of ndo()     (maxndo) = ",I7)') maxndo
            write(6,'(  20x,"  Abelian/C-inf-v flag   (symtyp) = ",I7)') symtyp
            write(6,'(  20x,"  Number of spin orbs      (nsrb) = ",I7,/)') nsrb
            write(6,'(  20x,"No. of dets in this CSF (nodi) = ",I7)') nodi

            do idet = 1, nodi
                istart = (idet - 1) * nelt
                write(6,'(/,20x,"Det No = ",I7," Coeff (cdi) = ",D13.6,/)') idet, cdi(idet)
                write(6,'(  20x,"    Sp. orbs in open shells (ndi) = ",20(I4,1x))') (ndo(istart+i), i=1,nelt)
            end do
        end if

        ! Save the number of determinants in the CSF into a variable for use during diagnostic printing later (if used). 
        ! The value of "nodi"  may increase during execution of this routine.     
        ninitial_dets = nodi  

        ! Compute the pre-multiplication factors which depend only on overal values - number of electrons and overal Spin. 
        fcta = -nelt * (nelt - 4)
        fctb = iss * (iss + 2)
        issp = iss + 1

        ! We initialize the number of output determinants to be the same as the number input.
        nodo = nodi

        ! Copy CSF expansion coefficients from input array to output array before we launch into the loop over Sz components      
        cdo(1:nodi) = cdi(1:nodi)

        ! To speed up the searches, build a binary search tree on top of the determinant storage
        ndo_ptr => ndo(1:maxndo)
        call bst % init(ndo_ptr, nelt)
        do id = 1, nodi
            call bst % insert(id)
            if (zdebug) then
                write(6,'(/,20x,"Bstree after add det #",I0,": ")') id
                call bst % output(22)
            end if
        end do

        ! Loop over Sz components of spin 
        if (zdebug) then
            write(6,'(/,20x,"Entering loop over components of spin",/)')
            write(6,'(  20x,"  fcta = ",D13.6)') fcta
            write(6,'(  20x,"  fctb = ",D13.6)') fctb
            write(6,'(  20x,"  issp = ",I6,/)') issp
        end if

        ! As shown in equation (6) of the Beebe and Stencil 1975 paper, referenced above, the spin ptojector is product over 
        ! spin operators. This is instanciated here as the DO loop to line 180. 

        spin_loop: do is = isd + 1, mxss + 1, 2

            ! As shown in equation (6) we omit the case k = S
            if (is == issp) cycle

            ! Compute coefficient multiplication factors for this Sz 
            fctc = (is - 1) * (is + 1)
            fctr = (fcta - fctc) / (fctb - fctc)

            if (zdebug) then
                write(6,'(/,20x,"Working on Spin Iteration (is) = ",I5,/)') is
                write(6,'(  20x,"  Factor fctc = ",D13.6)') fctc
                write(6,'(  20x,"  Factor fctr = ",D13.6,/)') fctr
            end if

            ! Copy all existing determinant coefficients for this CSF into cdo() and in doing so multiply by the factor
            ! pertaining to this Sz ("nodo" is the associated length of "cdo").
            cdo(1:nodi) = fctr * cdi(1:nodi)
            nodo = nodi

            ! Descend into loop over all determinants currently in this CSF
            !  - This may include determinants that we have created in previous iterations of the loop to line 180 as well
            !    as those in the initial input list.
            !  - "nd" points at the location of the determinants in ndo() for each iteration

            nd = 0

            determinant_loop: do id = 1, nodi

                fctr = four * cdi(id) / (fctb - fctc)

                if (zdebug) then
                    write(6,'(23x,"Working on determinant (id) = ",I5," of ",I5)') id, nodi
                    write(6,'(23x,"Coefficient, cdi() = ",D13.6  )') cdi(id)
                    write(6,'(23x,"Factor (fctr)      = ",D13.6,/)') fctr
                    write(6,'(23x,"Current op-shl det: ",20(I3,1x),/,(20x,20(I3,1x)) )') (ndo(nd+i), i=1,nelt)
                end if

                ! Now descend into loop over all PAIRS of electrons
                !  - This is implemented as a double DO loop over electrons

                first_electron_loop: do ia = 2, nelt
                    second_electron_loop: do ib = 1, ia - 1

                        ! Copy the determinant from its location in NDO() to a temporary area in NDTR().
                        !  - here, the determinant is composed of "nelt" electrons and is expressed as a list of occupied spin-orbitals
                        ndtr(1:nelt) = ndo(nd+1:nd+nelt)

                        ! Consider the present pair of electrons and find their spin
                        !  - The debug here is extremely useful but can generate massive amounts of output - uncomment if really needed.

                        ma  = ndtr(ia)
                        mb  = ndtr(ib)

                        mga = ms(ma)
                        mgb = ms(mb)

                        if (zdebug) then
                            write(6,'(23x,"Evaluating electron pair: ")')
                            write(6,'(23x,"  #1: idx(ia) = ",I3,"sporb(ma) = ",I3," Sz (mga) = ",I3  )') ia, ma, mga
                            write(6,'(23x,"  #2: idx(ib) = ",I3,"sporb(mb) = ",I3," Sz (mgb) = ",I3,/)') ib, mb, mgb
                        end if

                        ! Now look at the spins of this pair.
                        !
                        !   Remember that 0 means spin-up, 1 means spin-down.                  
                        !    
                        !   Option #    Electron a       Electron b
                        !   --------    ----------       ----------
                        !      1            0                0
                        !      2            1                0                     
                        !      3            0                1
                        !      4            1                1
                        !    
                        !   When the spin-orbitals have the same value, that is
                        !   options 1 and 4, we can't change them so all we do 
                        !   is augment the (output) coefficient of this determinant 
                        !   by the required factor - no more work to be done, so
                        !   we can proceed to the next electron pair.
                        !
                        !   For options 2 and 3, we can create a new determinant
                        !   by switching the electrons around. We preserve the Sz
                        !   value by doing that. We rely in the following that the
                        !   spin orbitals are created (see ms() mn() ...) in a
                        !   particular order ( Orb N spin-up, Orb N spin down,
                        !   ...)
                        !

                        if (mga == mgb) then
                            cdo(id) = cdo(id) + fctr
                            cycle
                        end if

                        if (mga == 0) then
                            ma = ma + 1
                            mb = mb - 1
                        else
                            ma = ma - 1
                            mb = mb + 1
                        end if

                        ndtr(ia) = ma
                        ndtr(ib) = mb

                        ! sort the determinant, keeping an eye on the sign of the permutation
                        flip = (mod(qsort(nelt, ndtr(1:nelt)), 2) /= 0)

                        if (zdebug) then
                            write(6,'(23x,"New valid determinant produced by this pair",/)')
                            write(6,'(23x,"New open shell det: ",20(I3,1x),/,(20x,20(I3,1x)) )') (ndtr(i), i=1,nelt)
                        end if

                        ! We ee to screen to see that the spin-orbitals we have just added to the determinant (ma and mb) do not already
                        ! occur elsewhere in it. if so, the new determinant created is not valid and must be rejected. 

                        do i = 1, nelt
                            if (ndtr(i) == ma .and. i /= ia) cycle
                            if (ndtr(i) == mb .and. i /= ib) cycle
                        end do

                        ! Given the "new" determinant just constructed in "ndtr" with associated coefficient "fctr",
                        ! we examine the list of "nodo" determinants for this CSF, defined in ndo()/cdo(), and merge
                        ! the "new" determinant into the list. This may mean adding it onto the end of cdo()/ndo() - see
                        ! comments in the routine for explanation.
                        !
                        ! So, stmrg() needs to know the maximum dimensions of ndo()/cdo() to monitor the extension of these arrays.
                        !
                        ! "nodo" will be updated on return if we add ndtr() and fctr into the list.      

                        call stmrg (nelt, maxcdo, maxndo, ndo, cdo, nodo, ndtr, merge(-fctr, fctr, flip), bst)

                    end do second_electron_loop
                end do first_electron_loop ! End of loop(s) over pairs of electrons 

                ! End of loop over determinants in this CSF, update the pointer "nd" to start at next CSF in NDO(). 
                ! Remember that a determinant consists of "nelt" consequtive spin orbitals in ndi()     

                nd = nd + nelt

            end do determinant_loop ! End loop over "nodi" current dets in CSF

            ! Remove any determinants whihch have an expansion coefficient less than "thres".
            ! Number of elements "nodo" in cdo() and ndo() may actually decrease here as elemenst are removed.

            call cntrct (nelt, nodo, ndo, cdo, thres)

            ! If we find that we have no determinants left in the CSF, we have an error condition.

            if (nodo == 0) then
                write(6,'(/,10x,"**** Error in prjct() ",/)')
                write(6,'(  10x,"After removing all determinants with expansion ")')
                write(6,'(  10x,"coefficients less than (thres) = ",D13.6)') thres
                write(6,'(  10x,"there are no determinants left in this CSF",/)') 
                stop 999
            end if

            ! Ok, so now copy the full set of coefficients back to array cdi(). Note that originally cdi() was of length
            ! "nodi" but we have added to it. The variable "nodo" now holds the total number of determinants - so we need 
            ! to reset "nodi" to that value.
            !
            ! We do this because the next iteration of the loop to line over spins starts by copying from cdi().

            nodi = nodo

            cdi(1:nodi) = cdo(1:nodi)

        end do spin_loop

        ! For Sigma wavefunctions in C-inf-v/D-inf-h we needed to consider the reflection symmetry.
        ! They can be Sigma(+) or Sigma(-)

        if (symtyp <= 1 .and. mgvn == 0) then
            if (zdebug) write(6,'(/,20x,"Analysis of Reflection operator required",/)')

            nodi = nodo
            cdi(1:nodi) = cdo(1:nodi)

            call rfltn (nelt, nodi, ndo, cdi, r, maxcdo, maxndo, thres, nodo, cdo, ndtr, mm, bst)

            if (nodo <= 0) then
                write(6,'(/,10x,"**** Error in prjct() ",/)')
                write(6,'(  10x,"After reflection analysis there are no ")')
                write(6,'(  10x,"no determinants left in the CSF",/)')
                stop 999
            end if

        end if

        ! Debug printout of the coefficients after spin projection

        if (zdebug) then
            write(6,'(/,20x,"Coefficients after projection ")')
            write(6,'(  20x,"  (but not yet normalized)    ",/)')
            write(6,'(  20x,"  No. of determinants (nodi) = ",I6,/)') nodi
            write(6,'(20x,I4,2x,D13.6)') (i, cdi(i), i=1,nodi) 
        end if

        ! Compute inverse sum of squared coefficients and normalize output

        tmp = snrm2(nodo, cdo, 1)
        tmp = one / sqrt(tmp)

        do i = 1, nodo
            cdo(i) = cdo(i) * tmp
        end do

        ! Return point 

        continue

        if (zdebug) then
            write(6,'(20x,"Final number of determinants in CSF (nodo) = ",I5,/)') nodo

            if (nodo /= ninitial_dets) then
                ia = nodo - ninitial_dets
                write(6,'(20x,"The length of the CSF has changed: ")')
                write(6,'(20x,"    Initial number of dets = ",I6)') ninitial_dets
                write(6,'(20x,"    # of dets added   = ",I6," by projection",/)') ia
            else
                write(6,'(20x,"# of dets in the CSF has NOT changed ",/)')
            end if

            write(6,'(/,20x,"**** prjct() - completed",/)') 
        end if

    end subroutine prjct


    !> \brief Project the wave function
    !>
    !> The subroutine \c projec controls the projection of the wavefunctions
    !> and writes out the final wavefunctions plus header information for future use.
    !>
    !> \note MAL 06/05/11 PROJEC has been considerably modified to take advantage of dynamic memory allocation.
    !>
    subroutine projec (sname, megul, symtyp, mgvn, s, sz, r, pin, nocsf, byproj, idiag, npflg, thres, &
                       nelt, nsym, nob, ndtrf, nftw, iposit, nob0, nob1, nob01, iscat, ntgsym, notgt, &
                       nctgt, mcont, gucont, mrkorb, mdegen, mflag, nobe, nobp, nobv, maxtgsym)

        use precisn,      only : wp
        use global_utils, only : mprod

        integer, intent(in)        :: iscat, mflag
        integer, intent(inout)     :: ntgsym
        integer                    :: byproj, idiag, iposit, megul, mgvn, nelt, nftw, nocsf,  &
                                      nsym, symtyp, npflg(6), num_csfs_unproj, num_dets_unproj,  &
                                      len_pkd_dets_unproj
        real(kind=wp)              :: pin, r, s, sz, thres
        character(len=80)          :: sname
        integer, dimension(ntgsym) :: gucont, mcont, mdegen, mrkorb, nctgt, notgt
        integer, dimension(nsym)   :: nob, nob0, nob01, nobe, nobp, nobv
        integer, dimension(nelt)   :: ndtrf
        integer, dimension(*)      :: nob1

        integer                    :: i,nalm,nb,nctarg,nd,nl,norb,junk, isd,iss,n,k,msum,isum,m,ierr, &
                                      num_dets,nreps,maxndi,maxcdi, maxndo,maxcdo,lenndo,lencdo,      &
                                      leniphase,leniphase0,maxtgsym

        !-----------------------------------------------------------------------
        ! Following are used to store the "unprojected" wavefunction
        !-----------------------------------------------------------------------

        integer                    :: wfn_unproj_num_csfs,                &
                                      wfn_unproj_num_dets,                &
                                      wfn_unproj_len_pkd_dets
        integer, allocatable       :: wfn_unproj_dets_per_csf(:),         &
                                      wfn_unproj_packed_dets(:),          &
                                      wfn_unproj_indx_1st_det_per_csf(:), &
                                      wfn_unproj_indx_1st_coeff_per_csf(:)
        real(kind=wp), allocatable :: wfn_unproj_coefficients_per_det(:)

        !-----------------------------------------------------------------------
        ! Following are used to store the "projected" wavefunction
        !-----------------------------------------------------------------------

        integer                    :: wfn_proj_num_csfs,  &
                                      wfn_proj_num_dets,  &
                                      wfn_proj_len_pkd_dets
        integer, allocatable       :: wfn_proj_dets_per_csf(:),           &
                                      wfn_proj_packed_dets(:),            &
                                      wfn_proj_indx_1st_det_per_csf(:),   &
                                      wfn_proj_indx_1st_coeff_per_csf(:)
        real(kind=wp), allocatable :: wfn_proj_coefficients_per_det(:)

        !-----------------------------------------------------------------------
        ! Following are used to store the table of spin-orbitals
        !-----------------------------------------------------------------------

        integer                    :: nsrb,noarg
        integer, allocatable       :: itab_sporb_indx_in_sym(:),  &
                                      itab_sporb_gu_value(:),     &
                                      itab_sporb_sym(:),          &
                                      itab_sporb_isz(:),          &
                                      itab_sporb_mpos(:),         &
                                      map_orbitals(:)

        !------------------------------------------------------------------------
        ! Following are used during the phase analysis of the wavefunctio
        !------------------------------------------------------------------------

        integer, allocatable       :: nconf(:),iphase(:),iphase0(:)

        !-------------------------------------------------------------------------
        ! Following logical flags are used to control the computation
        !
        ! They are provided to make the code easier to read than trying to
        ! remember that "byproj .eq. 0" means "bypass the projection of the
        ! wavefunction".
        !-------------------------------------------------------------------------

        logical :: zbypass_wfn_projection,                &
                   zadjust_wfn_phase_for_scattering,      &
                   zpositrons,ztarget_state_calculation,  &
                   zscattering_calculation,zabelian

        integer :: num_csfs_proj, num_dets_proj, len_pkd_dets_proj ! MAL 12/05/11 : None of these were declared in CG code. Why?
                                                                   ! They are not assigned values either, so most likely a problem here

        if (symtyp >= 2) then
            write(nftw,'(" MOLECULE SYMMETRY CASE,  symtyp =",I2)') symtyp
            junk = mprod(1, 1, npflg(6), nftw)
        end if

        nalm = 0 ! JMC initialization
        iss = s + s
        isd = sz + sz
        if (iss < isd) then
            write(nftw, 40)
            write(nftw, 46)
            stop
        end if

        !-------------------------------------------------------------------------
        ! The following logicals are created to improve code readability
        !-------------------------------------------------------------------------

        zbypass_wfn_projection = byproj == 0
        zadjust_wfn_phase_for_scattering = iscat > 0
        zpositrons = iposit /= 0
        zabelian = symtyp == 2

        !-------------------------------------------------------------------------
        ! Compute the table of spin-orbitals
        !-------------------------------------------------------------------------

        nsrb = sum(nob(1:nsym))
        norb = nsrb
        ! JMC set but not used      NORBB=(NORB*(NORB+1))/2
        nsrb = 2 * nsrb

        !------------------------------------------------------------------------
        ! Memory allocation for table of spin-orbitals
        !------------------------------------------------------------------------

        allocate(itab_sporb_indx_in_sym(nsrb), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'itab_sporb_indx_in_sym', ierr
            stop
        end if

        allocate(itab_sporb_gu_value(nsrb),stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'itab_sporb_gu_value', ierr
            stop
        end if

        allocate(itab_sporb_sym(nsrb), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'itab_sporb_sym', ierr
            stop
        end if

        allocate(itab_sporb_isz(nsrb), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'itab_sporb_isz', ierr
            stop
        end if

        allocate(itab_sporb_mpos(nsrb), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'itab_sporb_mpos', ierr
            stop
        end if

        allocate(map_orbitals(norb), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'mpos_orbitals', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! Compute the table of spin orbitals
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        ! For positrons in Abelian point groups there is a separate
        ! processing routine
        !-------------------------------------------------------------------------

        if (zpositrons .and. zabelian) then
            call pmkorbs (nob, nobe, nsym,          &
                          itab_sporb_indx_in_sym,   &
                          itab_sporb_gu_value,      &
                          itab_sporb_sym,           &
                          itab_sporb_isz,           &
                          nsrb, norb, nsrb,         &
                          map_orbitals,             &
                          itab_sporb_mpos,          &
                          iposit, symtyp)
        else
            call mkorbs (nob, nsym,                &
                        itab_sporb_indx_in_sym,    &
                        itab_sporb_gu_value,       &
                        itab_sporb_sym,            &
                        itab_sporb_isz,            &
                        norb, nsrb, map_orbitals,  &
                        itab_sporb_mpos,           &
                        iposit, nob1, nob01, symtyp)
        end if

        !-------------------------------------------------------------------------
        ! Validate the reference determinant
        !-------------------------------------------------------------------------

        isum = 0
        if (.not. zabelian) then
            msum = 0
            do i = 1, nelt
                m = ndtrf(i)
                msum = msum + itab_sporb_sym(m)
                isum = isum + 1 - itab_sporb_isz(m) - itab_sporb_isz(m)
            end do
        else
            msum = 1
            do i = 1, nelt
                m = ndtrf(i)
                msum = mprod(msum, itab_sporb_sym(m) + 1, 0, nftw)   !mal 12/05/11 nftw was missing in cg code. why?
                isum = isum + 1 - itab_sporb_isz(m) - itab_sporb_isz(m)
            end do
            msum = msum - 1
        end if

        !-------------------------------------------------------------------------
        ! Cross check with input
        !-------------------------------------------------------------------------

        if (abs(msum) /= mgvn) then
            write(nftw,9900)
            write(nftw,1190)
            stop
        end if

        if (abs(isum) /= isd) then
            write(nftw,9900)
            write(nftw,1195) abs(isum), isd
            stop
        end if

        !-------------------------------------------------------------------------
        ! Read the unprojected wavefunctions from unit megul
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ! MAL 06/05/11 : Read the dimension information only, in order to see
        ! how many CSFs, how many determinants there are and also how long the
        ! determinant array is
        !-------------------------------------------------------------------------

        call rdwf_getsize (megul, wfn_unproj_num_csfs, wfn_unproj_num_dets, wfn_unproj_len_pkd_dets)

        !-------------------------------------------------------------------------
        ! Allocate the arrays dynamically for storage of the unprojected
        ! wavefunction
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ! (1). Number of determinants per CSF
        !-------------------------------------------------------------------------

        allocate(wfn_unproj_dets_per_csf(wfn_unproj_num_csfs), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_unproj_dets_per_csf', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (2). The coefficient for each and every determinant
        !-------------------------------------------------------------------------

        allocate(wfn_unproj_coefficients_per_det(wfn_unproj_num_dets), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_unproj_coefficients_per_det',ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (3). Each and every packed determinant
        !-------------------------------------------------------------------------

        allocate(wfn_unproj_packed_dets(wfn_unproj_len_pkd_dets), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_unproj_packed_dets', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (4). Index into the list of determinants to the location for the first
        !      determinant of each CSF
        !      Note how this has one extra entry.
        !-------------------------------------------------------------------------

        allocate(wfn_unproj_indx_1st_det_per_csf(wfn_unproj_num_csfs+1), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_unproj_indx_1st_det_per_csf',ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (5). Index into the list of coefficients to the location for the first
        !      coefficient of each CSF
        !      Note how this has one extra entry
        !-------------------------------------------------------------------------

        allocate(wfn_unproj_indx_1st_coeff_per_csf(wfn_unproj_num_csfs+1), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_unproj_indx_1st_coeff_per_csf', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! Now able to read the unprojected wavefunction from file
        ! Note that this does not include the indexing arrays
        ! These are computed next
        !-------------------------------------------------------------------------

        call rdwf (megul, num_csfs_unproj, wfn_unproj_dets_per_csf,     &
                   num_dets_unproj, wfn_unproj_coefficients_per_det,    &
                   len_pkd_dets_unproj, wfn_unproj_packed_dets)

        !-------------------------------------------------------------------------
        ! Check that the values read back match to those read earlier
        !-------------------------------------------------------------------------

        if (num_csfs_unproj /= wfn_unproj_num_csfs) then
            write(nftw,*) ' Error 1 '
            stop 999
        end if

        if (num_dets_unproj /= wfn_unproj_num_dets) then
            write(nftw,*) ' Error 2 '
            stop 999
        end if

        if (len_pkd_dets_unproj /= wfn_unproj_len_pkd_dets) then
            write(nftw,*) ' Error 3 '
            stop 999
        end if

        !-------------------------------------------------------------------------
        ! Given the input wavefunction, indexing vectors for it are built here
        !
        !    wnf_unproj_indx_1st_det_per_csf()
        !       points to the location of the first deteminant for each CSF
        !       within the array ndi() which holds all the packed determinants
        !       in the wavefunction.
        !
        !    wfn_unproj_indx_1st_coeff_per_csf()
        !       points to the location of the first coefficient for each CSF
        !       within the array cdi() which holds all the coefficients
        !       (one per determinant) in the wavefunction
        !
        !-------------------------------------------------------------------------

        wfn_unproj_indx_1st_coeff_per_csf(1) = 1

        do n = 2, num_csfs_unproj
            wfn_unproj_indx_1st_coeff_per_csf(n) =         &
              wfn_unproj_indx_1st_coeff_per_csf(n-1) +     &
              wfn_unproj_dets_per_csf(n-1)
        end do

        wfn_unproj_indx_1st_coeff_per_csf(num_csfs_unproj+1) =    &
            wfn_unproj_indx_1st_coeff_per_csf(num_csfs_unproj)  +   &
            wfn_unproj_dets_per_csf(num_csfs_unproj)

        !-------------------------------------------------------------------------
        ! Index determinants now
        !-------------------------------------------------------------------------

        wfn_unproj_indx_1st_det_per_csf(1) = 1
        k = 1
        do n = 1, num_csfs_unproj
            num_dets = wfn_unproj_dets_per_csf(n)
            do m = 1,num_dets
                nreps = wfn_unproj_packed_dets(k)
                k = k + (2*nreps + 1)
            end do
            if (n <= num_csfs_unproj) then
                wfn_unproj_indx_1st_det_per_csf(n+1) = k
            end if
        end do

        !-------------------------------------------------------------------------
        ! If the user has requested it, then print the wavefunction
        !-------------------------------------------------------------------------

        if (npflg(1) > 5) then
            write(nftw,1285) megul
            call ptpwf (nftw,num_csfs_unproj,nelt,ndtrf,   &
                        wfn_unproj_dets_per_csf,           &
                        wfn_unproj_indx_1st_det_per_csf,   &
                        wfn_unproj_indx_1st_coeff_per_csf, &
                        wfn_unproj_packed_dets,            &
                        wfn_unproj_coefficients_per_det)
        end if

        !-------------------------------------------------------------------------
        ! Allocate storage space for the projected wavefunction
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ! The number of CSFs in the projected wavefunction will be no greater
        ! than the number in the unprojected wavefunction. It is possible to
        ! even delete a few due to the THRESHOLD criterion on coefficients.
        ! However, the number of determinants may grow due to the projection
        ! for open shell determinants  (see prjct()).
        !
        ! A factor of 10 is allowed here as a guesstimate...this may need revision
        !
        !-------------------------------------------------------------------------
        ! If no projection is to be done, the data is just copied straight over
        ! which means the subsequent code does not need to be rewritten with
        ! different variable names
        !-------------------------------------------------------------------------

        if (zbypass_wfn_projection) then
            wfn_proj_num_csfs     = wfn_unproj_num_csfs
            wfn_proj_num_dets     = wfn_unproj_num_dets
            wfn_proj_len_pkd_dets = wfn_unproj_len_pkd_dets
        else
            wfn_proj_num_csfs     = wfn_unproj_num_csfs
            wfn_proj_num_dets     = 10 * wfn_unproj_num_dets
            wfn_proj_len_pkd_dets = 10 * wfn_unproj_len_pkd_dets
        end if

        !-------------------------------------------------------------------------
        ! Allocate space for the projected wavefunction
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ! (1). Number of determinants per CSF
        !-------------------------------------------------------------------------

        allocate(wfn_proj_dets_per_csf(wfn_proj_num_csfs), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_proj_dets_per_csf', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (2). The coefficient for each and every determinant
        !-------------------------------------------------------------------------

        allocate(wfn_proj_coefficients_per_det(wfn_proj_num_dets), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_proj_coefficients_per_det', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (3). Each and every packed determinant
        !-------------------------------------------------------------------------

        allocate(wfn_proj_packed_dets(wfn_proj_len_pkd_dets), stat=ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_proj_packed_dets', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (4). Index into the list of determinants to the location for the first
        !      determinant of each CSF
        !
        !      This (and 5 below) have to have a "fake" N+1 CSF, so the size
        !      of the last Nth CSF can be known
        !-------------------------------------------------------------------------

        allocate(wfn_proj_indx_1st_det_per_csf(wfn_proj_num_csfs+1), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_proj_indx_1st_det_per_csf', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! (5). Index into the list of coefficients to the location for the first
        !      coefficient of each CSF
        !-------------------------------------------------------------------------

        allocate(wfn_proj_indx_1st_coeff_per_csf(wfn_proj_num_csfs+1), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9950) 'wfn_proj_indx_1st_coeff_per_csf', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! Project the wavefunction
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ! There is a special case -- the wavefunctions read are already projected.
        ! In this case the input is just copied to the output
        !-------------------------------------------------------------------------

        if (zbypass_wfn_projection) then
            wfn_proj_dets_per_csf         = wfn_unproj_dets_per_csf
            wfn_proj_packed_dets          = wfn_unproj_packed_dets
            wfn_proj_coefficients_per_det = wfn_unproj_coefficients_per_det
            wfn_proj_indx_1st_det_per_csf = wfn_unproj_indx_1st_det_per_csf
            wfn_proj_indx_1st_coeff_per_csf = wfn_unproj_indx_1st_det_per_csf
        else

            maxndi = size(wfn_unproj_packed_dets)
            maxcdi = size(wfn_unproj_coefficients_per_det)
            maxndo = size(wfn_proj_packed_dets)
            maxcdo = size(wfn_proj_coefficients_per_det)

            if (maxndi <= 0 .or. maxcdi <= 0) stop 901
            if (maxndo <= 0 .or. maxcdo <= 0) stop 902
    
            call wfgntr (mgvn, iss, isd, thres, r, symtyp, nelt,    &
                         nsym, nob, nob1, nob01, nobe, norb, nsrb,  &
                         itab_sporb_indx_in_sym,                    &
                         itab_sporb_gu_value,                       &
                         itab_sporb_sym,                            &
                         itab_sporb_isz,                            &
                         iposit, map_orbitals, itab_sporb_mpos,     &
                         wfn_unproj_num_csfs,                       &
                         ndtrf,                                     &
                         wfn_unproj_dets_per_csf,                   &
                         wfn_unproj_packed_dets,                    &
                         wfn_unproj_coefficients_per_det,           &
                         wfn_unproj_indx_1st_det_per_csf,           &
                         wfn_unproj_indx_1st_coeff_per_csf,         &
                         maxndi, maxcdi,                            &
                         wfn_proj_dets_per_csf,                     &
                         wfn_proj_packed_dets,                      &
                         wfn_proj_coefficients_per_det,             &
                         wfn_proj_indx_1st_det_per_csf,             &
                         wfn_proj_indx_1st_coeff_per_csf,           &
                         maxndo, maxcdo, lenndo, lencdo,            &
                         npflg, byproj, nftw, nalm)

            if (nalm /= 0) then
                write(nftw,9900)
                write(nftw,46)
                stop
            end if

        end if ! end of bypass wavefunction projection switch

        !-------------------------------------------------------------------------
        ! Print the projected CSFs
        !-------------------------------------------------------------------------

        if (npflg(3) /= 0) then
            write(nftw,'("1 OUTPUT FUNCTIONS IN PACKED FORM")')
            call ptpwf (nftw, wfn_proj_num_csfs, nelt, ndtrf, &
                        wfn_proj_dets_per_csf,                &
                        wfn_proj_indx_1st_det_per_csf,        &
                        wfn_proj_indx_1st_coeff_per_csf,      &
                        wfn_proj_packed_dets,                 &
                        wfn_proj_coefficients_per_det)
        end if

        !-------------------------------------------------------------------------
        ! Clean up dynamic storage that is no longer needed
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ! Deallocate the arrays used to read the unprojected wavefunction
        !-------------------------------------------------------------------------

        deallocate(wfn_unproj_dets_per_csf, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'wfn_unproj_dets_per_csf', ierr
            stop
        end if

        deallocate(wfn_unproj_coefficients_per_det, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'wfn_unproj_coefficients_per_det', ierr
            stop
        end if

        deallocate(wfn_unproj_packed_dets, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'wfn_unproj_packed_dets', ierr
            stop
        end if

        deallocate(wfn_unproj_indx_1st_det_per_csf, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'wfn_unproj_indx_1st_det_per_csf', ierr
            stop
        end if

        deallocate(wfn_unproj_indx_1st_coeff_per_csf, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'wfn_unproj_1st_coeff_per_csf', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! Calculate phase correction for SCATCI
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        ! Figure out which calculation is being done
        !-------------------------------------------------------------------------

        ztarget_state_calculation = iscat == 1
        zscattering_calculation  = iscat > 1

        !-------------------------------------------------------------------------
        ! Allocate workspace arrays needed in this step
        !
        ! Depends on type of calculation
        !
        ! Note that the phase array is needed for one of the write to file options
        !-------------------------------------------------------------------------

        allocate(nconf(nelt), stat = ierr)

        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9950) 'nconf', ierr
            stop
        end if

        !-------------------------------------------------------------------------
        ! iphase()
        !-------------------------------------------------------------------------

        if (ztarget_state_calculation) then
            leniphase = nocsf
        else
            leniphase = sum(nctgt(1:ntgsym))
        end if

        allocate(iphase(leniphase), stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9950) 'iphase', ierr
            stop 999
        end if

        leniphase = size(iphase)

        !-------------------------------------------------------------------------
        ! iphase0()
        !-------------------------------------------------------------------------

        allocate(iphase0(3*nocsf), stat = ierr)

        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9950) 'iphase0', ierr
            stop 999
        end if

        !-------------------------------------------------------------------------
        ! Branch to appropriate phase handler -- watch for error condition
        ! (that is, iscat <= 0)
        !-------------------------------------------------------------------------

        if (ztarget_state_calculation) then

            write(nftw, 3000)

            nctarg = nocsf
            ntgsym = -1

            call dophz0 (nftw, nocsf, nelt, ndtrf, nconf,  &
                         wfn_proj_indx_1st_det_per_csf,    &
                         wfn_proj_packed_dets,             &
                         lenndo,                           &
                         wfn_proj_indx_1st_coeff_per_csf,  &
                         wfn_proj_coefficients_per_det,    &
                         lencdo,                           &
                         iphase,                           &
                         npflg(5))

        else if (zscattering_calculation) then

            write(nftw, 3100)
            write(nftw, 3110) ntgsym,notgt
            write(nftw, 3120) nctgt
            write(nftw, 3130) mcont

            if (symtyp == 1) write(nftw, 3140) gucont
            if (npflg(5) > 0) write(nftw, 3150) mrkorb
            if (symtyp <= 1 .and. mgvn > 0) write(nftw, 3160) mdegen

            !-------------------------------------------------------------------------
            ! Count the number of continuum functions, that is a sum over all target
            ! states
            !-------------------------------------------------------------------------

            if (iposit .NE. 0) then
              nctarg = sum(nctgt(1:maxtgsym))
            else
              nctarg = sum(nctgt(1:ntgsym))
            endif

            !-------------------------------------------------------------------------
            ! Execute the phase alignment routine
            !-------------------------------------------------------------------------

            call dophz (nftw, nocsf, nelt, ndtrf, nconf, &
                        wfn_proj_indx_1st_det_per_csf,   &
                        wfn_proj_packed_dets,            &
                        lenndo,                          &
                        wfn_proj_indx_1st_coeff_per_csf, &
                        wfn_proj_coefficients_per_det,   &
                        lencdo,                          &
                        iphase,                          &
                        leniphase,                       &
                        iphase0,                         &
                        leniphase0,                      &
                        nctarg, nctgt, notgt, mrkorb,    &
                        mdegen, ntgsym, mcont,           &
                        symtyp, npflg(5))

        else

            write(nftw,9900)
            write(nftw,*) 'neither a target state nor a scattering run'
            stop 999

        end if

        !------------------------------------------------------------------------
        ! Write the header and (projected) wavefunctions back to unit MEGUL
        !------------------------------------------------------------------------
        !
        !------------------------------------------------------------------------
        ! iscat > 0 => SCATCI, DENPROP format
        !              otherwise SPEEDY format
        !------------------------------------------------------------------------

        if (iscat > 0) then
            write(nftw, 3167) megul
            call wrnfto (sname, mgvn, s, sz, r, pin, norb, nsrb,        &
                         nocsf, nelt, idiag, nsym, symtyp,              &
                         nob, ndtrf, wfn_proj_dets_per_csf,             &
                         nocsf + 1,                                     &
                         wfn_proj_indx_1st_coeff_per_csf,               &
                         wfn_proj_indx_1st_det_per_csf,                 &
                         wfn_proj_packed_dets,                          &
                         lenndo,                                        &
                         wfn_proj_coefficients_per_det,                 &
                         lencdo,                                        &
                         megul, nob1, 2 * nsym,                         &
                         npflg, thres, iposit, nob0, nob01, nctarg,     &
                         ntgsym, notgt, nctgt, mcont, gucont, iphase,   &
                         nobe, nobp, nobv, maxtgsym)
            write(nftw, 3170) megul
        else
            write(nftw, 3168) megul
            call wrwf (megul,                          &
                       num_csfs_proj,                  & ! MAL 12/05/11 : Not declared in CG code 
                       wfn_proj_dets_per_csf,          &
                       num_dets_proj,                  & ! MAL 12/05/11 : Not declared in CG code
                       wfn_proj_coefficients_per_det,  &
                       len_pkd_dets_proj,              & ! MAL 12/05/11 : Not declared in CG code
                       wfn_proj_packed_dets)
        end if

        !-------------------------------------------------------------------------
        ! Deallocate workspace arrays used in dophz0()/dophz()
        !-------------------------------------------------------------------------

        if (allocated(nconf)) then
            deallocate(nconf, stat = ierr)
            if (ierr /= 0) then
                write(nftw, 9900)
                write(nftw, 9960) 'nconf', ierr
                stop
            end if
        end if

        if (allocated(iphase)) then
            deallocate(iphase, stat = ierr)
            if (ierr /= 0) then
                write(nftw, 9900)
                write(nftw, 9960) 'iphase', ierr
                stop
            end if
        end if

        if (allocated(iphase0)) then
            deallocate(iphase0, stat = ierr)
            if (ierr /= 0) then
                write(nftw, 9900)
                write(nftw, 9960) 'iphase0', ierr
                stop
            end if
        end if

        !-------------------------------------------------------------------------
        ! Subroutine common return point
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ! Deallocate storage used for the table of spin-orbitals
        !-------------------------------------------------------------------------

        deallocate(itab_sporb_indx_in_sym, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'itab_sporb_indx_in_sym', ierr
            stop
        end if

        deallocate(itab_sporb_gu_value, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'itab_sporb_gu_value', ierr
            stop
        end if

        deallocate(itab_sporb_sym, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'itab_sporb_sym', ierr
            stop
        end if

        deallocate(itab_sporb_isz, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'itab_sporb_isz', ierr
            stop
        end if

        deallocate(itab_sporb_mpos, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'itab_sporb_mpos', ierr
            stop
        end if

        deallocate(map_orbitals, stat = ierr)
        if (ierr /= 0) then
            write(nftw, 9900)
            write(nftw, 9960) 'mpos_orbitals', ierr
            stop
        end if

        write(nftw, 8000)

        !------------------------------------------------------------------------
        ! Format statements
        !------------------------------------------------------------------------

        40 format('  S.LT.SZ')
        1000 format(/,5x,'Projection and phase alignment of wavefunction ',/, &
                    5x,'============================================== ',//,&
                    5x,'Input data: ',/)
        1005 format(5x,'  Sname = ',a,/,     &
                    5x,'  Mgvn  = ',i10,/,   &
                    5x,'  S     = ',f10.4,/, &
                    5x,'  Sz    = ',f10.4,/, &
                    5x,'  R     = ',f10.4,/, &
                    5x,'  Pin   = ',f10.4,/, &
                    5x,'  Nocsf = ',i10,/,   &
                    5x,'  Idiag = ',i10,//)
        1007 format(5x,'  Number of electrons in system   (nelt) = ',i5,//, &
                    5x,'  Reference determinant: ',//, &
                    5x,'     (refdet) = ',10(i5,1x),/, &
                    (21x,10(i5,1x)))

        1020 format(5x,'  Point group (symmetry) of nuclear framework (symtyp)         = ',i3,/)
        1021 format(5x,'  This is the C-inf-v point group',/) 
        1022 format(5x,'  This is the D-inf-h point',/) 
        1023 format(5x,'  This is an Abelian point group ',/)

        1030 format(5x,'  Bypassing wavefunction projection (byproj .eq. 0)',/)
        1031 format(5x,'  Wavefunction will be projected (byproj .ne. 0)',/)

        1035 format(5x,'  Adjusting phase of wavefunction for scattering               (iscat .gt. 0)',/)
        1036 format(5x,'  Not adjusting phase of wavefunction for scattering               (iscat .le. 0)',/)

        1038 format(/,5x,'  Print flags (npflg) ',/,&
                    5x,'  ------------------- ',/,&
                    5x,'  1. Unprojec and projec wavefunctions  : ',i5,/,&
                    5x,'  2.                                    : ',i5,/,&
                    5x,'  3.                                    : ',i5,/,&
                    5x,'  4.                                    : ',i5,/,&
                    5x,'  5. Target or scattering phase compute : ',i5,/,&
                    5x,'  6. Abelian point grp multiplctn table : ',i5,/)

        1099 format(/,5x,'**** End of the input data',/)             

        1110 format(/,5x,'Total number of orbitals  (norb) = ',i8,/, &
                    5x,'Triangulation of norb    (norbb) = ',i8)
        1120 format(/,5x,'Total number of spin-orbs (nsrb) = ',i8,/)
        1130 format(5x,'Spin orbitals table of quantum numbers',//, &
                    5x,'  I      N      G      M      S     MPOS   ',/,&
                    5x,'-----  -----  -----  -----  -----  -----   ')
        1140 format((5x,6(i5,2x)))
        1145 format(/,5x,'**** End of table of spin-orbitals ',/)

        1160 format(/,5x,'User defined quantum numbers of ref determinant        :',//,&
                    5x,'  mgvn = ',i5,/,    &
                    5x,'  S    = ',f8.3,/,  &
                    5x,'  Sz   = ',f8.3,//, &
                    5x,'and locally computed vars for spin from S,Sz: ',//,&
                    5x,'  iss  = ',i5,/,    &
                    5x,'  isd  = ',i5)
        1170 format(/,5x,'For check, computed q-numbers of ref det:',//,&
                    5x,'  2*Sz + 1 = ',i5)
        1180 format(  5x,'  irreducible representation   = ',i5,//,&
                    5x,'   (Note: totally symmetric representation = 0) ',/)

        1185 format(5x,'  Lambda value = ',i5,/)
        1190 format(5x,'  Symmetry quantum number in refdet is not MGVN')
        1195 format(5x,'  Sz in refdet (',i0,') is not = SZ (',i0,')')

        1270 format(/,5x,'Starting to build indexes for the wavefunction',/)
        1280 format(/,5x,'Finished building indexes for the wavefunction',/)
        1285 format(/,5x,'Wavefunctions read from input file on unit ',i5,/)

        2000 format(/,5x,'The wavefunction will be projected',//, &
                    5x,'This means that the wavefunction on unit ',i5,' is an',/,&
                    5x,'unprojected wavefunction.',/)
        2010 format(5x,'Data read from the wavefunction on unit: ',i5,//,   &
                    5x,'   number of CSFs                = ',i10,/,  &
                    5x,'   number of determinants        = ',i10,/,  &
                    5x,'   length of packed determinants = ',i10,//, &
                    5x,'This data will be used to allocate dynamic storage',/,&
                    5x,'in which to hold the wavefunction.',/)
        2100 format(5x,'Projected CSFs in packed format:',/)
        3000 format(/,5x,'Computing phase for target state',/)
        3100 format(/,5x,'Performing phase correction for target',/,&
                    5x,'states in a scattering run',/)
        46 format('  DUE TO ALARM CONDITION THIS RUN WAS TERMINATED')
        3110 format(/' CI target data for SCATCI:', &
                    //' Number of target symmetries in expansion,   NTGSYM =',&
                    i5/' Number of continuum orbs for each state, NOTGT =',&
                    20I5,/,(' ',20I5))
        3120 format(' Number of CI components for each state,      NCTGT =',20I5,/,('  ',20I5))
        3130 format(' Continuum M projection  for each state,      MCONT =',20I5,/,('  ',20I5))
        3140 format(' Continuum G/U symmetry  for each state,     GUCONT =',20I5,/,('  ',20I5))
        3150 format(' Marked continuum orbital for each state,    MRKORB =',20I5,/,('  ',20I5))
        3160 format(' Degenerate coupling case flag               MDEGEN =',20I5,/,('  ',20I5))

        3167 format(5x,'Writing projected CSFs to unit ',i5,' in format',/,&
                    5x,'required by the SCATCI/DENPROP programs ',/)
        3168 format(5x,'Writing projected CSFs to unit ',i5,' in format',/,&
                    5x, 'required by the SPEEDY program ',/)
        3170 format(5x,'Data on CSFs has been written to file (MEGUL) = ',i5/)
        8000 format(5x,'***** Wavefn projection (projec()) - completed',/)
        9900 format(/,5x,'***** Error in: projec() ',//)
        9950 format(5x,'Cannot allocate space for array ',a,//,&
                    5x,'Return status from allocate() = ',i8,/)
        9960 format(5x,'Cannot de-allocate space for array ',a,//,&
                    5x,'Return status from deallocate() = ',i8,/)

    end subroutine projec


    !> \brief Print the CSFs.
    !>
    !> This provides a complete text dump of all CSFs in terms of their packed determinants.
    !>
    !> \param nftw      File unit for text output of the program.
    !> \param nocsf     Number of wave functions (CSFs).
    !> \param nelt      Number of electrons (and size of \c ndtrf).
    !> \param ndtrf     Reference determinant (spinorbitals per electron).
    !> \param nodi      Array with number of determinants per CSF.
    !> \param indi      Array with offsets in \c ndi per CSF.
    !> \param icdi      Array with offsets in \c cdi per CSF.
    !> \param ndi       Packed determinants.
    !> \param cdi       Deteminant coefficients.
    !>
    subroutine ptpwf (nftw, nocsf, nelt, ndtrf, nodi, indi, icdi, ndi, cdi)

        use precisn, only: wp

        integer :: nelt, nftw, nocsf
        real(kind=wp), dimension(*) :: cdi
        integer, dimension(*) :: icdi, indi, ndi, ndtrf, nodi
        intent (in) cdi, icdi, indi, ndi, ndtrf, nelt, nftw, nocsf, nodi

        integer :: i, k, ma, mb, mc, md, n

        write(nftw,'(" REFERENCE DETERMINANT"//(1X,20I5))') (ndtrf(i), i=1,nelt)
        write(nftw,'("  CSF",9X,"COEFFICIENT",2X,"NSO"/)')
        do n = 1, nocsf
            ma = nodi(n)
            mb = indi(n)
            mc = icdi(n) - 1
            md = ndi(mb)
            write(nftw,'(1x,i4,d20.10,i5,2x,20i5/(32x,20i5))') n, cdi(mc+1), md, (ndi(mb+i), i=1,2*md)
            mb = mb + md + md + 1
            do k = 2, ma
                md = ndi(mb)
                write(nftw,'(5x,d20.10,i5,2x,20i5/(32x,20i5))') cdi(mc+k), md, (ndi(mb+i), i=1,2*md)
                mb = mb + md + md + 1
            end do
        end do

    end subroutine ptpwf


    !> \brief   Sort integer array.
    !> \authors J Benda
    !> \date    2018
    !>
    !> Sort array in-place in non-descending order. Currently implemented as a stable insertion sort.
    !> This algorithm is advantageous for short and *almost sorted* arrays, which is the case for augmented
    !> open-shell-only subsets of determinants in CONGEN, typically resulting in asymptotic complexity of O(n).
    !>
    !> \param  n  Length of the array to sort.
    !> \param  a  Integer array to sort.
    !>
    !> \return Number of swaps done. This is needed elsewhere to update determinant signs.
    !>
    integer function qsort (n, a) result (swaps)

        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: a

        integer :: b, i, j, k

        swaps = 0

        do i = 2, n                         ! process all positions except for the first

            j = i - 1                       ! previous position
            if (a(j) <= a(i)) cycle         ! this is the correct order -> no change needed
            b = a(i)                        ! remember the current element a(i)

            do while (j > 1)                ! find the correct position for this element to the left of its current position
                if (b < a(j - 1)) then
                    j = j - 1
                else
                    exit
                end if
            end do

            do k = i, j + 1, -1             ! shift all elements left of i and larger than a(i) one position to the right, vacating j
                a(k) = a(k - 1)
            end do

            a(j) = b                        ! put the current element into the vacated space
            swaps = swaps + i - j           ! update the number of swaps needed to achieve the new order

        end do

    end function qsort


    !> \brief Read CSF.
    !>
    !> Routine \ref congen_distribution::wfn stores the CSF data in buffers of fixed size.
    !>
    !> When the buffers are full, they are emptied to disk
    !> and reused. This is why we can have several sets to read here.
    !>
    subroutine rdwf (nft, k1, nodi, k2, cdi, k3, ndi)

        use precisn, only : wp

        integer                     ::  k1, k2, k3, n1, n2, n3, i, nft, st
        integer, dimension(*)       ::  nodi, ndi
        real(kind=wp), dimension(*) ::  cdi
    
        rewind nft

        k1 = 0
        k2 = 0
        k3 = 0

        do
            read(nft, iostat = st) n1, (nodi(k1+i), i=1,n1)
            if (st /= 0) exit
            k1 = k1 + n1

            read(nft, iostat = st) n2, (cdi(k2+i), i=1,n2)
            if (st /= 0) exit
            k2 = k2 + n2

            read(nft, iostat = st) n3, (ndi(k3+i), i=1,n3)
            if (st /= 0) exit
            k3 = k3 + n3
        end do

        rewind nft

    end subroutine rdwf


    !> \brief Reads the size of the wavefunction
    !>
    !> Reads the information giving the size of the data used for the wavefunction, for example
    !> the number of determinants, but does not read the actual data,  such as the determinants.
    !>
    !> This is really just a stripped down version of the routine \ref rdwf which reads the full data.
    !>
    !> \param iunit     The logical unit on which the wavefucntion data is located.
    !> \param num_csfs  On return, number of CSFs stored in the file.
    !> \param num_dets  On return, number of determinants summed over all CSFs in the file.
    !> \param len_dets  On return, total summed length of all arrays of packed determinants in the file.
    !>
    !> \note MAL 10/05/2011: This subroutine is new to congen and is included to bring congen into
    !>                       line with the changes that were made in \ref projec in order to utilize
    !>                       dynamic memory
    !>
    subroutine rdwf_getsize (iunit, num_csfs, num_dets, len_dets)

        use precisn, only : wp

        integer :: iunit, num_csfs, num_dets, len_dets
        integer :: ncsfs, ndets, ldets, ios, ntemp

        rewind iunit

        ncsfs = 0
        ndets = 0
        ldets = 0

        ! The while loop reads records from the logical unit until the end of file, or some other error occurs.
        do

            ! First record is a CSF counter (followed by array with determinant counts per CSF)
            read(iunit, iostat = ios) ntemp
            if (0 /= ios) exit
            ncsfs = ncsfs + ntemp

            ! Second record is a counter for determinants (followed by coefficients per CSF per determinant)
            read(iunit, iostat = ios) ntemp
            if (0 /= ios) exit
            ndets = ndets + ntemp

            ! Third record is the counter of integer data comprising determinants themselves (followed by packed determinants)
            read(iunit, iostat = ios) ntemp
            if (0 /= ios) exit
            ldets = ldets + ntemp

        end do

        rewind iunit

        num_csfs = ncsfs 
        num_dets = ndets
        len_dets = ldets

    end subroutine rdwf_getsize


    !> \brief Add mirror-reflected spin-orbitals.
    !>
    !> Used only for C_infv and D_infh.
    !>
    subroutine rfltn (nelt, nodi, ndi, cdi, r, mxnd, ndmxp, thres, nodo, cdo, ndtr, mm, bst)

        use precisn,       only : wp
        use consts,        only : half => xhalf
        use congen_bstree, only : det_tree

        integer :: mxnd, ndmxp, nelt, nodi, nodo
        real(kind=wp) :: r, thres
        real(kind=wp), dimension(*) :: cdi, cdo
        integer, dimension(*) :: mm, ndi, ndtr
        type(det_tree) :: bst
        intent (in) cdi, mm, nodi, r
        intent (inout) nodo, bst

        real(kind=wp) :: cfd
        integer :: i, j, ma, nd

        cdo(1:nodi) = half * cdi(1:nodi)
        
        nd = 0
        nodo = nodi
        do i = 1, nodi
            cfd = half * r * cdi(i)
            do j = 1, nelt
                ma = ndi(nd+j)
                if (mm(ma) /= 0) ma = ma + sign(2, mm(ma))
                ndtr(j) = ma
            end do
            if (mod(qsort(nelt, ndtr(1:nelt)), 2) /= 0) cfd = -cfd
            call stmrg (nelt, mxnd, ndmxp, ndi, cdo, nodo, ndtr, cfd, bst)
            if (nodo < 0) return
            nd = nd + nelt
        end do

        call cntrct (nelt, nodo, ndi, cdo, thres)

    end subroutine rfltn


    !> \brief Real strided vector norm.
    !>
    !> This is actually an in-house version of the BLAS level 1 routine of the same name, included here to avoid
    !> dependency on BLAS library.
    !>
    !> \param n      Length of the vector
    !> \param array  Vector of real numbers
    !> \param istep  Stride
    !>
    !> \return Square of L2 norm of the vector.
    !>
    function snrm2 (n, array, istep)

        use precisn, only : wp
        use consts,  only : zero => xzero

        integer :: istep, n
        real(kind=wp), dimension(*) :: array
        real(kind=wp) :: snrm2
        intent (in) array, istep, n

        integer :: i, ii

        snrm2 = zero
        if (n <= 0) return

        ii = 1
        do i = 1, n
            snrm2 = snrm2 + array(ii)**2
            ii = ii + istep
        end do

    end function snrm2


    !> \brief Add determinant to a list of determinants
    !>
    !> Given an existing list of determinants in cdo()/ndo() and a new single determinant cdi/ndi(), the new determinant
    !> is merged into the list and the list extended if necessary.
    !>
    !> This operation is used during spin-projection of a CSF.
    !>
    !> \note MAL 16/05/11 : Modified to bring the subroutine into line with the changes made to projec
    !>
    !> \param nelt      Number of electrons in each det.
    !> \param maxcdo    Dimension of cdo.
    !> \param maxndo    Dimension of ndo.
    !> \param ndo       List of determinants each with "nelt" spin-orbitals.
    !> \param cdo       Coefficient for each det in ndo.
    !> \param nodo      On input is the number of determinants in cdo/ndo. Will be updated for output if the data in cdi/ndi
    !>                  new entry.
    !> \param ndi       A single det of "nelt" spin orbs which has to be merged into ndo.
    !> \param cdi       Single coefficient going with the single determinant defined in ndi.
    !> \param bst       Binary search tree used for fast localization of determinants.
    !>
    subroutine stmrg (nelt, maxcdo, maxndo, ndo, cdo, nodo, ndi, cdi, bst)

        use precisn,       only : wp
        use congen_bstree, only : det_tree

        real(kind=wp), parameter   :: one = 1.0_wp

        type(det_tree) :: bst

        integer :: nelt         ! Number of electrons in each det.
        integer :: maxcdo       ! Dimension of cdo
        integer :: maxndo       ! Dimension of ndo       
        integer :: ndo(maxndo)  ! List of determinants each with "nelt" spin-orbs
        integer :: nodo         ! On input is the number of dets in cdo()/ndo(). Will be updated
                                ! for output if the data in cdi/ndi() is merged into cdo()/ndo() as a new entry.
        integer :: ndi(nelt)    ! A single det of "nelt" spin orbs which has to be merged into ndo().

        real(kind=wp) :: cdo(maxcdo)  ! Coeff. for each det in ndo()
        real(kind=wp) :: cdi          ! Single coeff. going with the single determinant defined in ndi()

        integer  i, ibase, idet, n

        logical, parameter :: zdebug = .false.

        ! Debug banner header

        if (zdebug) then
            write(6,'(/,30x,"===> STMRG() <====",/)')
            write(6,'(  30x,"Input data: ")')
            write(6,'(  30x,"  nelt   = ",I6)') nelt
            write(6,'(  30x,"  maxcdo = ",I6)') maxcdo
            write(6,'(  30x,"  maxndo = ",I6)') maxndo
            write(6,'(  30x,"  cdi    = ",D13.6)') cdi
            write(6,'(/,30x,"# of dets in current (cdo,ndo) list (nodo) = ",I5,/)') nodo

            ibase = 0
            do idet = 1, nodo
                write(6,'( 30x,I5,2x,D13.6,20(I4,1x))') idet, cdo(idet), (ndo(ibase+i), i=1,nelt)
                ibase = ibase + nelt
            end do
        end if

        ! Try to find the input determinant in the list of determinants.
        !  - This should take around O(log2(nodo)) operations

        n = bst % locate_det(ndi)

        ! Now either update the coefficient (if det found), or insert new determinant.

        if (n > 0) then

            if (zdebug) then
                write(6,'(/,30x,"New determinant found in ""existing"" list:")')
                write(6,'(  30x,"  Found in exisiting list at n = ",I6)') n
                write(6,'(  30x,"  Coefficient cdo(n)           = ",D13.6)') cdo(n)
                write(6,'(  30x,"  Augment coeff,  cdi*sign     = ",D13.6,/)') cdi
            end if

            cdo(n) = cdo(n) + cdi

        else

            ! Ok, so we get here if we have to add the determinant onto the end of the exisitng list.

            if (zdebug) then
                write(6,'(30x,"Adding new determinant to end of list",/)')
            end if

            if (nodo + 1 > maxcdo .or. (nodo + 1) * nelt > maxndo) then
                write(6,'(/,10x,"**** Error in stmrg(): ",/)')
                write(6,'(  10x,"Insufficient space to add extra determinant onto")')
                write(6,'(  10x,"end of list.")')
                write(6,'(  10x,"  cdo() : required ",I10," available ",I10)') nodo + 1, maxcdo
                write(6,'(  10x,"  ndo() : required ",I10," available ",I10,/)') (nodo + 1) * nelt, maxndo
                stop 999
            end if

            cdo(nodo + 1) = cdi
            ndo(nodo * nelt + 1 : nodo * nelt + nelt) = ndi(1:nelt)
            nodo = nodo + 1

            call bst % insert(nodo)

            if (zdebug) then
                write(6,'(30x,"Bstree after stmrg: ")')
                call bst % output(32)
            end if

        end if

        ! Return point

        if (zdebug) then
            write(6,'(/,30x,"# of dets in current (cdo,ndo) list (nodo) = ",I5,/)') nodo

            ibase = 0
            do idet = 1,nodo
                write(6,'( 30x,I5,2x,D13.6,20(I4,1x))') idet, cdo(idet), (ndo(ibase+i), i=1,nelt)
                ibase = ibase + nelt
            end do

            write(6,'(/,30x,"***** STMRG() - completed ",/)')
        end if

    end subroutine stmrg


    !> \brief Form spin eigenstates
    !>
    !> Takes the wavefunction generated by CSFGEN and transforms it to be
    !> fully in accord with the spin quantum numbers of the system.
    !>
    !> \note MAL 10/05/2011: The changes made to wfgntr have made so as to make the
    !>       subroutine compatible with the changes made to its calling subroutine, \ref projec. 
    !>
    subroutine wfgntr (mgvn, iss, isd, thres, r, symtyp, nelt, nsyml, nob, nobl, nob0l, nobe, norb, nsrb,             &
                       mn, mg, mm, ms, iposit, map, mpos, nocsf, ndtrf, nodi, ndi, cdi, indil, icdil, maxndi, maxcdi, &
                       nodo, ndo, cdo, indo, icdo, maxndo, maxcdo, lenndo, lencdo, npflg, byproj, nftw, nalm)                            

        use precisn,      only : wp
        use global_utils, only : mprod

        ! Local temporary storage
        !  - len_cdit has to be >= the maximum number of determinants in any one CSF

        integer, parameter         :: len_cdit = 5000
        real(kind=wp)              :: cdit(len_cdit)

        ! Fixed constants

        real(kind=wp), parameter   :: verysmall = 1.0d-30

        ! Integer variables passed in the argument list 

        integer                    :: nftw  ! logical unit for the printer      
        integer                    :: symtyp, byproj 
        integer                    :: lencdo ! final usage in cdo()
        integer                    :: lenndo ! final usage in ndo()
        integer                    :: maxcdo ! maximum available in cdo()
        integer                    :: maxndo ! maximum available in ndo()
        integer                    :: mgvn
        integer                    :: iss 
        integer                    :: isd 
        integer                    :: norb 
        integer                    :: ndmx 
        integer                    :: ncmx
        integer                    :: ndmxp
        integer                    :: nsyml
        integer                    :: nelt
        integer                    :: iposit
        integer                    :: nocsf,nsrb
        integer                    :: maxndi,maxcdi
        integer, dimension(nsyml)  :: nob(nsyml),nobl(nsyml),nob0l(nsyml), nobe(nsyml) 
        integer, dimension(nsrb)   :: mn,mg,mm,ms,map,mpos 
        integer, dimension(nelt)   :: ndtrf
        integer                    :: npflg(6)
        real(kind=wp)              :: r,thres

        ! As elsewhere in the code, the wavefunction is represented as a 
        ! set of data packed into arrays. 
        !
        !     For the input wavefunction we have:
        !
        !            nodi()  is the number of dets in each CSF   
        !            ndi()   holds the dets in packed format, that is as    
        !                      number of replacements + replaced + replacing
        !            cdi()   is the coefficient for each determinant
        !
        !     We need two further indexing arrays, each of length,
        !     NOCSF+1 to handle the wavefunction data:
        !
        !            icdil(n) points at the first entry in the coefficients
        !                     array, cdi(), for CSF "n".
        !
        !            indil(n) points at the first entry in the determinants
        !                     array, ndi(), for CSF "n".
        !
        ! Remember these have to be one larger than the number of CSFs. We compute the storage size of CSF by looking at the location 
        ! of the "following one" - hence need to add one on at the end.

        integer, dimension(nocsf)           :: nodi 
        integer, dimension(maxndi)          :: ndi
        integer, dimension(nocsf+1)         :: indil, icdil
        real(kind=wp), dimension(maxcdi)    :: cdi

        ! Similar arrays exist for the output wavefunction and they are declared in the argument list too.

        integer, dimension(nocsf)           :: nodo 
        integer, dimension(maxndi)          :: ndo
        integer, dimension(nocsf+1)         :: indo, icdo
        real(kind=wp), dimension(maxcdi)    :: cdo

        ! Local integer variables
        integer                             :: n, i, num_dets_input, ipos_dets_input, ipos_coeffs_out, ipos_coeffs_in,   &
                                               ipos_dets_out, idet, isum, msum, nsrbs, needed, ierr, lenmop, nalm,       &
                                               ipos_this_det, nreps, me, mf, idop, idcp, ieltp, num_dets_final, mxss
        integer, dimension(nelt)            :: mdop, mdcp
        integer, dimension(nsrb)            :: mdc, mdo, ndta
        logical, allocatable                :: flip(:)

        ! Local real variables
        real(kind=wp)                       :: mysum

        integer, allocatable                :: mop(:) ! Holds all expanded open shell per CSF (see allocation)

        ! Following are used to analyze CSF dimensions on input
        integer                             :: icsf_with_max(1)
        integer                             :: max_num_dets_input
        integer                             :: max_num_dets_output

        ! Local fixed logical values
        logical, parameter                  :: zdebug = .false.

        ! Intrinsic Fortran functions used
        intrinsic                           :: sqrt, maxloc, minloc

        ! Debug banner header
        if (zdebug) then
            write(nftw,1000)
            write(nftw,1010) mgvn,iss,isd,thres,r,nsyml,nelt,nocsf,nsrb
            write(nftw,1020) norb, ndmx, ncmx, ndmxp
            write(nftw,1030) maxcdo, maxndo
            write(nftw,1035) 
            write(nftw,1036) (i,mn(i),mg(i),mm(i),ms(i),mpos(i),i=1,nsrb)
            write(nftw,1037) 
        end if

        ! We compute the number of spin-orbitals which are not 
        ! degenerate. For C-inf-v and D-inf-h this means sigma
        ! type. For Abelian point groups it is all orbitals.

        nsrbs = 0
        select case (symtyp)
            case (0) ; nsrbs = 2 * nob(1)
            case (1) ; nsrbs = 2 * (nob(1) + nob(2))
            case (2) ; nsrbs = 2 * sum(nob(1:nsyml))
            case default ; write(nftw, 9900) ; stop 
        end select 

        if (zdebug) then
            write(nftw,1040) nsrbs
        end if

        !=========================================================================
        !                                                 
        !     C O U N T   M A X I M A   F R O M   C S F s
        !
        !=========================================================================

        ! Looking at the number of determinants per CSF, we work out
        ! the CSF which has the maximum number of determinants.

        max_num_dets_input = maxval(nodi)
        icsf_with_max      = maxloc(nodi)

        ! mop() needs to hold the expanded determinants when processing each CSF one at a time (see popnwf().
        ! Each det is then a list of "ieltp" spin orbs, where "ieltp" <= "nelt".
        !
        ! So we can compute lenmop and allocate the array. We over allocate by a factor (10 ?) here as the 
        ! later processing requirements are not precisely known.

        lenmop = 7 * nelt * max_num_dets_input
        allocate(mop(lenmop), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900) 
            write(nftw,9925) lenmop
            stop 999
        end if

        ! There will be one permutation sign per determinant.

        allocate(flip(max_num_dets_input), stat = ierr)
        if (ierr /= 0) then
            write(nftw,9900)
            write(nftw,9926) max_num_dets_input
            stop 999
        end if

        !=========================================================================
        !                                                 
        !     L O O P   O V E R   C S F s
        !
        !=========================================================================

        if (zdebug) then
            write(nftw, 3000)
        end if

        nalm = 0 
        ipos_coeffs_out = 1
        ipos_dets_out   = 1

        do n = 1, nocsf 
            num_dets_input  = nodi(n) 
            ipos_dets_input = indil(n)
            icdo(n) = ipos_coeffs_out 
            indo(n) = ipos_dets_out

            if (zdebug) then
                write(nftw,3010) n, nocsf, num_dets_input, ipos_dets_input, ipos_coeffs_out, ipos_dets_out
            end if 

            !-------------------------------------------------------------------------
            ! Step 1: Validate spatial and spin quantum numbers 
            !         for all determinants in this CSF 
            !-------------------------------------------------------------------------

            ! Compute the overall change in spatial quantum number and in Sz component of spin for each determinant relative to the
            ! reference determinant. The overall change in spatial and in spin quantum numbers relative to the reference determinant
            ! should be zero. If not, then such a deterinant does not match the overall system quantum numbers and we have an error
            ! condition. Remember that we have already verified the quantum numbers of the reference determinant before calling this
            ! routine.

            ipos_this_det = ipos_dets_input 

            do idet = 1, num_dets_input
                nreps = ndi(ipos_this_det) 

                if (zdebug) then
                    write(nftw,3020) n, idet, num_dets_input, nreps 
                end if 

                if (nreps /=  0) then ! If 0 would be reference det 
                    if (zdebug) then
                        write(nftw,3030) (ndi(ipos_this_det+i), i=1,nreps)
                        write(nftw,3035) (ndi(ipos_this_det+nreps+i), i=1,nreps)
                    end if

                    isum = 0 ! Sum of Sz over spin-orbs 

                    select case (symtyp)

                        !... C-inf-v and D-inf-h
                        case (:1)
                            msum = 0 
                            do i = 1, nreps 
                                me   = ndi(ipos_this_det + i) 
                                mf   = ndi(ipos_this_det + nreps + i) 
                                isum = isum + ms(me) - ms(mf) 
                                msum = msum + mm(me) - mm(mf) 
                            end do

                        !... Abelian point groups 
                        case (2)
                            msum = 1 
                            do i = 1, nreps 
                                me   = ndi(ipos_this_det + i) 
                                mf   = ndi(ipos_this_det + nreps + i) 
                                isum = isum + ms(me) - ms(mf) 
                                msum = mprod(msum, mprod(mm(me)+1, mm(mf)+1, 0, nftw), 0, nftw)                                                       
                            end do
                            msum = msum - 1 

                        !... Erroneous "symtyp" value
                        case (3:) 
                            write(nftw,9900)
                            stop

                    end select 

                    if (msum /= 0) then 
                        write(nftw, 9900)
                        write(nftw, 9214) idet, n 
                        stop 999
                    end if

                    if (isum /= 0) then 
                        write(nftw,9900) 
                        write(nftw,9215) idet, n
                        stop 999
                    end if

                end if ! End of if test on zero replacements

                ! Align pointers for the next determinant. Remember, a determinant is stored as:
                !   - the number of replacements/replaced,
                !   - each replaced    spin-orb,
                !   - each replacement spin-orb,
                ! therefore this determinant was of length 2*md+1

                ipos_this_det = ipos_this_det + nreps + nreps + 1 

            end do 

            if (zdebug) then
                write(nftw, 3090)
            end if

            !-------------------------------------------------------------------------
            ! Step 2: Find the number of electrons which are in orbitals that are not fully occupied.
            !         These are known as "open shells". Also find the list of the spin-orbitals corresponding to these. 
            !-------------------------------------------------------------------------

            ! popnwf() operates on the complete set of determinants defining this CSF. It is worth remembering that in Alchemy
            ! a CSF is defined firstly as an assignment of orbital occupation numbers. One, or more, CSFs can then be formed 
            ! by distributing electrons to the spin-orbitals associated with this orbital assignment and 
            ! then generating the required coupling coefficients.
            !
            ! On successful return, "ieltp" will hold the number of electrons in open shells. This may even be zero.

            num_dets_input  = nodi(n) 
            ipos_dets_input = indil(n)

            call popnwf (nsrb, nsrbs, nelt, ndtrf, lenmop, mdop, mdcp, mop, mdc, mdo, ndta, &
                         num_dets_input, ndi(ipos_dets_input), idop, idcp, ieltp, flip, ierr)

            ! Test for error condition and print details 
            if (ierr == 0) then
                if (zdebug) then
                    write(nftw,4010) ieltp
                end if
            else
                write(nftw,9900)
                select case (ierr)
                    case (1)  ; write(nftw,233) n 
                    case (2)  ; write(nftw,235) n 
                    case (3:) ; write(nftw,237) n 
                end select
                stop 999
            end if

            !-------------------------------------------------------------------------
            ! Step 3: Normalize the expansion coefficients associated with the determinants in this CSF.
            !         These coefficients were obtained from products of the Clebsch-Gordan coupling coefficients.
            !-------------------------------------------------------------------------

            ! The computation depends on how many electrons are in open shells. 
            ! Normalized expansion coefficients are placed into the array "cdo()" during this process.

            ipos_coeffs_in = icdil(n) - 1

            ! update cdi with the sort permutation sign
            do i = 1, num_dets_input
                if (flip(i)) cdi(ipos_coeffs_in + i) = -cdi(ipos_coeffs_in + i)
            end do

            !-------------------------------------------------------------------------
            !        0 or 1 ELECTRONS in OPEN SHELLS
            !-------------------------------------------------------------------------

            if (ieltp <= 1) then  ! 0 or 1 electrons in open shells.

                if (zdebug) then
                    write (nftw,4510) 
                    write (nftw,4520) 
                    write (nftw,4522) ipos_coeffs_in, ipos_coeffs_out
                    write (nftw,4525) (i, cdi(ipos_coeffs_in+1+i-1), i=1,num_dets_input)
                end if

                mysum = snrm2(num_dets_input, cdi(ipos_coeffs_in+1), 1)

                if (zdebug) write(nftw,4530) mysum

                if (mysum < verysmall) then ! monitor for very small numbers
                    write(nftw, 9900)
                    write(nftw, 2222) mysum
                    stop 999
                end if

                mysum = 1.0_wp/SQRT(mysum)

                do i = 1, num_dets_input 
                    cdo(ipos_coeffs_out+i-1) = mysum * cdi(ipos_coeffs_in+1+i-1)
                end do

                num_dets_final = num_dets_input

                if (zdebug) then
                    write(nftw,4540) 
                    write(nftw,4525) (i, cdo(ipos_coeffs_out+i-1), i=1,num_dets_final)
                end if

            !-------------------------------------------------------------------------
            !         >= 2 ELECTRONS in OPEN SHELLS 
            !-------------------------------------------------------------------------

            else ! >= 2 electrons in open shells - project wfn

                if (zdebug) then
                    write(nftw,4610) 
                    write(nftw,4520) 
                    write(nftw,4522) ipos_coeffs_in, ipos_coeffs_out
                    write(nftw,4525) (i, cdi(ipos_coeffs_in+1+i-1), i=1,num_dets_input)
                end if

                ! We need to copy the expansion coefficients to temporary local storage 
                if (num_dets_input > len_cdit) then 
                    write(nftw,9900)
                    write(nftw,9935) n, num_dets_input, len_cdit
                    stop 999
                else
                    do i = 1, num_dets_input 
                        cdit(i) = cdi(ipos_coeffs_in + i)
                    end do
                end if

                if (zdebug) then
                    write(nftw,4625) 
                    write(nftw,4525) (i, cdit(i), i=1,num_dets_input)
                end if

                ! Call prjct() to perform the spin projection
                !
                !    ieltp = number of electrons 
                !    mxss  = maximum Spin for projection
                !    ma    = number of ddterminants
                !    mop   = input determinants (overwritten on output)
                !    cdit  = input coefficients with each determinant
                !    nod   = Return code (>0 means number of dets in output)
                !    cdo() = Coefficients of determinants after projection. These are stored into cdo() beginning at
                !            location ipos_coeffs_out, which is updated for every CSF.
                !   lencdo = space reaining in cdo() array - monitored 

                mxss = ieltp
                lencdo = maxcdo - ipos_coeffs_out + 1
                lenmop = size(mop)

                call prjct (ieltp, mxss, num_dets_input, mop, cdit, num_dets_final, cdo(ipos_coeffs_out), &
                            lencdo, mgvn, iss, isd, thres, r, ndta, mm, ms, lenmop, symtyp, nsrb)                   

                if (zdebug) then
                    write(nftw,4540) 
                    write(nftw,4525) (i, cdo(ipos_coeffs_out+i-1), i=1,num_dets_final)
                end if

            end if

            !-------------------------------------------------------------------------
            ! Step 4: Now we can move the determinants for this 
            !         CSF into place in the output array.
            !-------------------------------------------------------------------------

            if (zdebug) then
                write(nftw, 4060) 
                write(nftw, 4062) num_dets_final, ieltp, idop, idcp, ipos_dets_out, maxndo 
            end if

            call pkwf (num_dets_final, ieltp, cdo(ipos_coeffs_out), mop, idop, mdop, idcp, mdcp, &
                       nftw, ipos_dets_out, ndo, maxndo, n)


            ! Record the number of dets for this CSF in nodo()
            ! Update pointer to next location in cdo().

            nodo(n) = num_dets_final  
            ipos_coeffs_out = ipos_coeffs_out + num_dets_final 


            ! Screen for exhaustion of memory in ndo(), cdo()
            if (ipos_dets_out >= maxndo) then
                write(nftw,9900) 
                write(nftw,9981) n, maxndo
                stop 999
            end if

            if (ipos_coeffs_out >= maxcdo) then
                write(nftw,9900) 
                write(nftw,9982) n, maxcdo
                stop 999
            end if

            ! Finished with this CSF
            if (zdebug) then
                write (nftw,4690) n
            end if

        end do

        !=========================================================================
        !                                                 
        !     E N D   L O O P   O V E R   C S F s
        !
        !=========================================================================

        ! Finished with dynamic array mop()
        if (allocated(mop)) then
            deallocate(mop, stat = ierr)
            if (ierr /= 0) then
                write(nftw,9900)
                write(nftw,9927) lenmop 
                stop 999
            end if
        end if

        ! Finished with dynamic array flip()
        if (allocated(flip)) then
            deallocate(flip, stat = ierr)
            if (ierr /= 0) then
                write(nftw,9900)
                write(nftw,9928) max_num_dets_input
                stop 999
            end if
        end if

        ! Need starting values for the "N+1" th CSF (which does not exist). These also serve as an upper bound for the "N"-th CSF
        ! and are used in loops in the subsequent code which evaluates symbolic energy expressions.
        indo(nocsf+1) = ipos_dets_out 
        icdo(nocsf+1) = ipos_coeffs_out 

        ! High watermark in both arrays needs to be returned to caller.
        lenndo = ipos_dets_out   - 1
        lencdo = ipos_coeffs_out - 1

        ! Looking at the number of determinants per CSF, we work out the CSF which has the maximum number of determinants.
        max_num_dets_output = MAXVAL( nodo )
        icsf_with_max       = MAXLOC( nodo )

        write(nftw,2992) icsf_with_max(1), max_num_dets_output

        ! Subroutine return point
        write(nftw,7990) lenndo, maxndo, lencdo, maxcdo
        write(nftw,8000) 

        ! Format statements
        233  format(I6,' TH WF IN ERROR, (NO. OF OPEN SO NOT =)') 
        237  format(I6,' TH WF IN ERROR,(NELTP=0,BUT NOD GT. 1)') 
        235  format('  NEED MORE SPACE FOR MOP IN',I5,' TH WF',/, &
                    '  Increase parameter LNDT in input')         
        1000 format(/,10x,'====> WFGNTR() <====',/)
        1010 format(  10x,'Input data: ',/,      &
                      10x,'  mgvn  = ',i10,/,    &
                      10x,'  iss   = ',i10,/,    &
                      10x,'  isd   = ',i10,/,    &
                      10x,'  thres = ',f12.5,/,  &
                      10x,'  r     = ',f12.5,/,  &
                      10x,'  nsyml = ',i10,/,    &
                      10x,'  nelt  = ',i10,/,    &
                      10x,'  nocsf = ',i10,/,    &
                      10x,'  nsrb  = ',i10)
        1020 format(  10x,'  norb  = ',i10,/,    &
                      10x,'  ndmx  = ',i10,/,    &
                      10x,'  ncmx  = ',i10,/,    &
                      10x,'  ndmxp = ',i10,/)
        1030 format(  10x,'  maxcdo = ',i10,/,   &
                      10x,'  maxndo = ',i10,/)
        1035 format(5x,'Spin orbitals table of quantum numbers',//,     &
                    5x,'  I      N      G      M      S     MPOS   ',/, &
                    5x,'-----  -----  -----  -----  -----  -----   ')
        1036 format((5x,6(i5,2x)))
        1037 format(/,5x,'**** End of table of spin-orbitals ',/)
        1040 format(10x,'No. non-degenerate spin orbitals (nsrbs) = ',i6,/)
        1055 format(10x,'Allocating ',i8,' integers for array mop() ',/)
        1500 format(/,10x,'Allocated ',i10,' integer units to array nr() ',/)
        2222 format(/,' Sum IN WFGNTR =',E20.12,//) 
        2990 format(/,10x,'On input CSF ',i7,' has the largest ',/, &
                      10x,'number of determinants = ',i7,/)
        2992 format(/,10x,'On output CSF ',i7,' has the largest ',/, &
                      10x,'number of determinants = ',i7,/)
        3000 format(/,10x,'Entering loop over CSFs ',/)
        3010 format(/,10x,'>>>> Processing input CSF ',i10,' of ',i10,//,      &
                      10x,'Number of determinants (num_dets_input) = ',i7,/,   &
                      10x,'1st in pkd dets array (ipos_dets_input) = ',i7,//,  &
                      10x,'1st in Out coefs arry (ipos_coeffs_out) = ',i7,/,   &
                      10x,'1st in Output dets arry (ipos_dets_out) = ',i7)
        3020 format(/,15x,'CSF ',i7,' - Det. ',i6,' of ',i6,//, &
                      20x,'Number of replacements (nreps) = ',i5,/)
        3030 format(20x,'Replaced spin orbs    : ',20(i3,1x))
        3035 format(20x,'Replacments spin. orbs: ',20(i3,1x))
        3090 format(/,15x,'Quantum numbers for all determinants in this CSF',/, &
                      15x,'have been validated.',/)
        4010 format(15x,'Number of electrons in open shells (ieltp) = ',i5,/)
        4050 format(15x,'The expansion coefficients for each determinant ',/ &
                    15x,'have been normalized.',/)
        4060 format(/,15x,'This projected CSF will now be packed into',/  &
                      15x,'into the array holding all output packed CSF',/)
        4062 format(15x,'Data sent into PKWF() follows: (old name/new name)',/, &
                    15x,'    nod   (num_dets_final) = ',i10,/, &
                    15x,'    neltp (ieltp)          = ',i10,/, &
                    15x,'    idop  (idop)           = ',i10,/, &
                    15x,'    idcp  (idcp)           = ',i10,/, &
                    15x,'    no    (ipos_dets_out)  = ',i10,/, &
                    15x,'    ndmx  (maxndo)         = ',i10,/)
        4510 format(15x,'Normalizing CSF expansion coefficients using ',/, &
                    15x,'the SNRM2 function. Projection is not needed.'/)
        4520 format(/,15x,'CSF coefficients before normalization: ',/)
        4522 format(15x,'Storage locs for input and output coefficients: ',//, &
                    15x,'  For cdi(), ipos_coeffs_in  = ',i7,//, &
                    15x,'  For cdo(), ipos_coeffs_out = ',i7,/)
        4525 format(15x,i6,'.  ',2x,f13.7)
        4530 format(/,15x,'Sum - sqrs of coefficients (this CSF) = ',d13.6,/)
        4540 format(/,15x,'CSF coefficients after normalization: ',/)
        4610 format(15x,'Normalizing CSF expansion coefficients using ',/, &
                    15x,'projection because it has >= 2 electrons in  ',/, &
                    15x,'open shells.                                 '/)
        4625 format(/,15x,'Coefficients have been copied into CDIT(): ',/)
        4690 format(/,10x,'<<<< Completed processing of CSF number ',i6,/)
        7990 format(/,10x,'At end of wfgntr(), usage of memory in output',/,  &
                      10x,'arrays is as follows:                        ',//, &
                      10x,'  Array     Used     Max avail               ',/,  &
                      10x,'  -----  ---------   ---------               ',/,  &
                      10x,'  ndo()  ',i9,2x,i9,/,                             &
                      10x,'  cdo()  ',i9,2x,i9,/)
        8000 format(/,10x,'**** WFGNTR() - completed ',/)

        ! Error messages
        9214 format(5x,'Spatial symmetry for determinant ',i5,/, &
                    5x,'in  CSF ',i5,' does not match ref determinant',/)
        9215 format(5x,'Sz quantum number for determinant ',i5,/, &  
                    5x,'in  CSF ',i5,' does not match ref determinant',/)
        9900 format(/,5x,'***** Error in: wfngtr() ',//)
        9925 format(5x,'Cannot allocate array mop() of length (lenmop) ',i8,/)
        9926 format(5x,'Cannot allocate array flip() of length (max_num_dets_input) ',i8,/)
        9927 format(5x,'Cannot de-alloc array mop() of length (lenmop) ',i8,/)
        9928 format(5x,'Cannot de-alloc array flip() of length (max_num_dets_input) ',i8,/)
        9935 format(5x,'Insufficient space in cdit() for projection  ',/, &
                    5x,'CSF ',i10,' has ',i10,' determinants         ',/, &
                    5x,'cdit() holds the expansion coefficient       ',/, &
                    5x,'for each determinant and len_cdit must be    ',/, &
                    5x,'at least as big as the number of determinants',/, &
                    5x,'Currenttly it is fixed at ',i10,/)
        9955 format(5x,'WFNGTR() has received an error on return',/, &
                    5x,'fromPOPNWF(). Code is terminating now.',/) 
        9981 format(5x,'CSF ',i7,' has exhausted space available in ndo',/, &
                    5x,'Available (maxndo) = ',i10)    
        9982 format(5x,'CSF ',i7,' has exhausted space available in cdo',/, &
                    5x,'Available (maxcdo) = ',i10)    

    end subroutine wfgntr


    !> \brief WRNFTO - WRite wavefunction data to unit NFTO
    !>
    !> \note MAL 10/05/2011 : This subroutine has been changed to bring it into
    !>       line with the changes that have been made to 'projec' in order to
    !>       utilize dynamic memory. 'wrnfto' has also been modified in order to
    !>       comply to the F95 standard
    !>
    subroutine wrnfto (sname, mgvn, s, sz, r, pin, norb, nsrb, nocsf, nelt, idiag, nsym, symtyp,  &
                       nob, ndtrf, nodo, m, icdo, indo, ndo, lndi, cdo, lcdi, nfto, nobl, nx,     &
                       npflg, thres, iposit, nob0, nob0l, nctarg, ntgsym, notgt, nctgt, mcont,    &
                       gucont, iphz, nobe, nobp, nobv, maxtgsym)

        use precisn, only : wp

        integer                    :: symtyp, i, nfto, iposit, norb, nsrb, idiag, itg, maxtgsym, mgvn, &
                                      norbw, nsrbw, lcdi, lndi, nsym, nelt, m, nx, nctarg, ntgsym, nocsf
        integer, parameter         :: iwrite = 6
        integer, dimension(nsym)   :: nob, nob0, nobe, nobp, nobv
        integer, dimension(nelt)   :: ndtrf
        integer, dimension(nocsf)  :: nodo
        integer, dimension(m)      :: icdo, indo
        integer, dimension(lndi)   :: ndo
        integer, dimension(ntgsym) :: notgt, nctgt, mcont, gucont
        integer, dimension(nctarg) :: iphz
        integer, dimension(20)     :: nobw
        integer, dimension(nx)     :: nobl, nob0l
        integer, dimension(6)      :: npflg
        real(kind=wp),dimension(lcdi)   :: cdo
        real(kind=wp)              :: s, sz, r, pin, thres
        character(120)             :: name 
        character(80)              :: sname

        !-------------------------------------------------------------------------
        !     Debug banner header
        !-------------------------------------------------------------------------      

        write(iwrite,'(/,5x,"Writing final results to disk (wrnfto) ")')
        write(iwrite,'(  5x,"-------------------------------------- ")')
        write(iwrite,'(  5x,"  Logical unit         (nfto)   = ",I8)') nfto
        write(iwrite,'(  5x,"  Positrons present    (iposit) = ",I8)') iposit
        write(iwrite,'(  5x,"  Number of electrons  (nelt)   = ",I8)') nelt
        write(iwrite,'(  5x,"  Number of CSFs       (nocsf)  = ",I8)') nocsf
        write(iwrite,'(  5x,"  Number of symmetries (nsym)   = ",I8)') nsym
        write(iwrite,'(  5x,"  Point group flag     (symtyp) = ",I8)') symtyp
        write(iwrite,'(  5x,"  Orbitals per symm    (nob)    :  ",(20(1x,I0)))') (nob(i), i=1,nsym)
        write(iwrite,'(  5x,"  Total # of orbitals  (norb)   = ",I8)') norb
        write(iwrite,'(  5x,"  Total # of spin orbs (nsrb)   = ",I8)') nsrb

        write(iwrite,'(/,5x,"  Length of index  arrays - icdo,indo (m)    = ",I8)') m
        write(iwrite,'(  5x,"  Length of coeff. array  - cdo       (lcdi) = ",I8)') lcdi
        write(iwrite,'(  5x,"    (equals # dets in wfn)")')
        write(iwrite,'(  5x,"  Length of packed dets array - ndo   (lndi) = ",I8)') lndi

        write(iwrite,'(/,5x,"  Spin quantum number  (s)      = ",F8.2)') s
        write(iwrite,'(  5x,"  Z-projection of spin (sz)     = ",F8.2)') sz
        write(iwrite,'(  5x,"  Reflection quant. #  (r)      = ",F8.2)') r
        write(iwrite,'(  5x,"  Pin                  (pin)    = ",F8.2)') pin
        write(iwrite,'(  5x,"                       (idiag)  = ",I8,/)') idiag

        name = ' '
        name = sname

        norbw = norb
        nsrbw = nsrb

        nobw(1:nsym) = nob(1:nsym)

        !--------------------------------------------------------------------------
        !     Rewind the unit before writing on it
        !--------------------------------------------------------------------------

        rewind nfto

        !--------------------------------------------------------------------------
        !     Slightly different formats if positrons are involved
        !--------------------------------------------------------------------------

        if (iposit == 0) then
            write(nfto) name, mgvn, s, sz, r, pin, norbw, nsrbw,  nocsf, nelt, lcdi, idiag, &
                        nsym, symtyp, lndi, npflg, thres, nctarg, ntgsym
            if (ntgsym > 0) write(nfto) iphz, nctgt, notgt, mcont, gucont
            if (ntgsym <= 0) write(nfto) iphz
            write(nfto) (nobw(i), i=1,nsym), ndtrf, nodo, iposit, nob0, nobl, nob0l
        else
            write(nfto) name, mgvn, s, sz, r, pin, norbw, nsrbw, nocsf, nelt, lcdi, idiag, &
                        nsym, symtyp, lndi, npflg, thres, nctarg, maxtgsym
            if (maxtgsym > 0) then
                write(nfto) iphz
                write(nfto) (nctgt(itg),  itg=1,maxtgsym) 
                write(nfto) (notgt(itg),  itg=1,maxtgsym) 
                write(nfto) (mcont(itg),  itg=1,maxtgsym) 
                write(nfto) (gucont(itg), itg=1,maxtgsym)
            end if
            if (maxtgsym <= 0) write(nfto) iphz
            write(nfto) (nobw(i), i=1,nsym), ndtrf, nodo, iposit, nob0, nobl, nob0l
            write(nfto) (nobe(i), i=1,nsym)
            write(nfto) (nobp(i), i=1,nsym)
            write(nfto) (nobv(i), i=1,nsym)
        end if

        !-------------------------------------------------------------------------
        ! Now come the determinants per CSF
        !-------------------------------------------------------------------------

        write(nfto) icdo, indo
        write(nfto) ndo
        write(nfto) cdo

    end subroutine wrnfto


    !> \brief Write wave functions using SPEEDY format
    !>
    !> This subroutine writes all determinants in the current CSF to a binary file that has the following format:
    !> \verbatim
    !>    [WP] n1, [INTEGER] nodo(1:n1)
    !>    [WP] n2, [WP]      cdo(1:n2)
    !>    [WP] n3, [INTEGER] ndo(1:n3)
    !> \endverbatim
    !>
    !> On entry, the file will be rewindws to its beginning.
    !> On return, the output file is rewinded to its beginning, again.
    !>
    !> \param nft   An open binary file for output.
    !> \param n1    Length of \c nodo (number of determinants).
    !> \param nodo  Determinant sizes.
    !> \param n2    Length of \c cdo (number of determinants).
    !> \param cdo   Determinant factors within the CSF.
    !> \param n3    Length of \c ndo (number of integers defining the determinants).
    !> \param ndo   Packed determinants.
    !>
    subroutine wrwf (nft, n1, nodo, n2, cdo, n3, ndo)

        use precisn, only : wp

        integer :: n1, n2, n3, nft
        real(kind=wp), dimension(n2) :: cdo
        integer, dimension(n3) :: ndo
        integer, dimension(n1) :: nodo
        intent (in) cdo, n1, n2, n3, ndo, nft, nodo

        rewind nft

        write(nft) n1, nodo
        write(nft) n2, cdo
        write(nft) n3, ndo

        rewind nft

    end subroutine wrwf

end module congen_projection
