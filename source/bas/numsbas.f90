!> \brief   Evaluate given STO orbitals on radial grid
!> \authors J Benda
!> \date    2018
!>
!> NUMSBAS a variation on NUMCBAS that evaluates STO orbitals instead of the continuum wave functions.
!> This is useful for generating STO-like Gaussian basis sets e.g. when reproducing older STO-based
!> calculations.
!>
!> NUMSBAS uses the same input namelist as NUMCBAS, where many irrelevant entries are ignored. The only
!> entries that are used are: \c title, \c lunumb, \c lval, \c rlim and the new entries \c nval and \c expt
!> with the principal quantum number and the exponential coefficient of the STO function to sample.
!>
!> The STO function is defined as
!> \f[
!>      \xi_{nlm}(\zeta, \mathbf{r}) = S_n(\zeta) r^{n - 1} \mathrm{e}^{-\zeta r} Y_{lm}(\hat{\mathbf{r}})
!> \f]
!> with normalization constant
!> \f[
!>      S_n(\zeta) = \sqrt{\frac{(2\zeta)^{2n + 1}}{(2n)!}} \,.
!> \f]
!>
program numsbas

    use consts,       only: xzero
    use global_utils, only: print_ukrmol_header
    use numcbas_data, only: maxorb, rtol, tinyy, absacc, maxrx, nfta, hrx, nix, irx
    use precisn,      only: wp

    implicit none

    integer  :: i, ir, ira = 0, lunumb = 13, nval = 1, lval = 0, ibug(3) = 0
    real(wp) :: btol = 0.2_wp, charge = 0.0_wp, ecmax = 0.0_wp, rlim = 10.0_wp, tiny = 1.0E-11_wp, expt = 0
    real(wp), allocatable :: r(:), orb(:)
    character(len=120) :: title

    ! NUMCBAS-compatible namelist
    ! - many entries are not used
    ! - new entries "nval" and "expt" for the principal quantum number and exponent of the STO function
    namelist /input/ title, lunumb, nix, irx, hrx, nval, lval, ibug, btol, tiny, ecmax, rlim, charge, expt

    ! set up the default "integration grid"; here used simply as evaluation grid, because STO orbitals
    ! can be directly evaluated without solving any differential equations (copied from NUMCBAS)
    nix = 3
    hrx = (/ 1.E-02_wp, 2.E-02_wp, 2.605E-02_wp, (xzero, i = 4, maxrx) /)
    irx = (/ 30, 120, 500, (0, i = 4, maxrx) /)

    ! greet the user and read the input
    call print_ukrmol_header(nfta)
    write(nfta, '(//11x," Program NUMSBAS",//)')
    read(5, input)

    ! evaluate the radial grid and adjust the storage for evaluated orbital
    call evaluate_grid(nfta, ira, r, rlim)
    allocate(orb(ira))

    ! evaluate the STO function
    call evaluate_STO(expt, nval, lval, ira, r, orb)

    if (ibug(2) /= 0) then
        WRITE(NFTA, '(//11X,"FUNCTION NO.",I3/)') 1
        WRITE(NFTA, '(i5,f10.4,D16.8)') 0, xzero, xzero, (ir, r(ir), orb(ir), ir = 1, ira)
    end if

    ! write output to the file
    write(lunumb) title
    write(lunumb) 1, lval, ira
    write(lunumb) (  r(i), i = 1, ira + 1)
    write(lunumb) (orb(i), i = 1, ira + 1)
    write(lunumb) xzero

contains

    !> \brief   Create the evaluation grid
    !> \authors J Benda
    !> \date    2018
    !>
    !> Given the \c nix, \c hrx and \c irx variables, construct radial grid compatible with given \c rlim.
    !>
    !> \todo This is almost in its entirety directly copied from NUMCBAS. It would be better if a single version
    !>       of this subroutine existed and was shared by both programs.
    !>
    subroutine evaluate_grid (nfta, ira, r, rlim)

        use consts,       only: xzero
        use numcbas_data, only: maxrx, hrx, nix, irx

        integer,  intent(in)    :: nfta
        integer,  intent(inout) :: ira
        real(wp), intent(in)    :: rlim
        real(wp), intent(inout), allocatable :: r(:)

        integer  :: i, ir, ir1, ir2, ier, irxo, ira1, nix1
        real(wp) :: x, h, trange, rmat, rmat1

        rmat = rlim

        ! check for inappropriate or incomplete input data concerning NIX, IRX, HRX
        if (nix > maxrx) then
            write(nfta, 80) nix, maxrx
            stop
        end if
        do i = 1, nix
            if(irx(i) == 0 .or. hrx(i) == xzero) then
                write(nfta, 90) i
                stop
            end if
        end do

        write(nfta, 51) nix
        ier = 0
        irxo = 0
        trange = xzero
        do i = 1, nix
            trange = trange + (irx(i) - irxo) * hrx(i)
            if (i == nix .and. abs(trange - rmat) >= rtol) then
                irx(i) = int((rmat - trange) / hrx(i)) + irx(i)
                if(mod(irx(i), 2) /= 0) irx(i) = irx(i) + 1
                write(nfta, 52) irx(i)
            end if
            irxo = irx(i)
            write(nfta, 53) i, irx(i), hrx(i)
            if (mod(irx(i), 2) /= 0) ier = 1
        end do

        if (ier /= 0) then
            write(nfta, 70)
            stop
        end if
        ira = irx(nix)

        write(nfta, 68) lval, charge, rmat
        write(nfta, 69) ecmax, btol, tiny

        allocate(r(ira))

        x = xzero
        ir2 = 0
        do i = 1, nix
            ir1 = ir2 + 1
            ir2 = irx(i)
            h = hrx(i)
            do ir = ir1, ir2
                x = x + h
                r(ir) = x
            end do
        end do

        rmat = r(ira)
        if (nix <= 1) go to 7
        nix1 = nix - 1
        if (ira - irx(nix1) <= 100) go to 8
    7   nix1 = nix
        nix = nix + 1
        irx(nix) = ira
        irx(nix1) = ira - 100
        hrx(nix) = hrx(nix1)
    8   continue
        ira1 = irx(nix1)
        rmat1 = r(ira1)
        if (ibug(1) /= 0) then
            write(nfta, 57) nix
            write(nfta, 65) (i, irx(i), hrx(i), i = 1, nix)
            write(nfta, 58) ira, rmat, ira1, rmat1
        end if

    51  format(//11X,'INTEGRATION MESH INPUT DATA'//6X,'NIX   =',I3,5X,'No. of integration regions with different step-', &
               'sizes'//10X,'I',7X,'IRX',15X,'HRX')
    52  format(//5X,'Number of points (IRX) in last subrange has',' been changed to',5X,I5)
    53  format(1x,2(5X,I5),5X,D20.10,5X,I5)
    57  format(///11X,'NUMERICAL BASIS INTEGRATION MESH'//11X,'NIX   =',I5,5X,'NUMBER OF INTEGRATION REGIONS'//11X,'I',7X,'IRX', &
               15X,'HRX')
    58  format(//6X,'IRA   =',I5,10X,'TOTAL NUMBER OF INTEGRATION ','STEPS'/6X,'RMAT  =',F9.4,6X,'R-MATRIX BOUNDARY RADIUS'/6X, &
               'IRA1  =',I5,10X,'NUMBER OF OUTWARD INTEGRATION STEPS'/6X,'RMAT1 =',F9.4,6X,'MATCHING RADIUS')
    65  format(7x,i5,5x,I5,5X,D20.10)
    68  format(//11X,'NUMERICAL BASIS CALCULATION INPUT DATA'//6X,'LVAL   =',I8,5X,'Angular Momentum'/6x,'CHARGE =',F8.1,5X, &
               'Effective charge'/6x,'RMAT   =',F8.1,5X,'R-matrix boundary radius')
    69  format(//11X,'SEARCHING PROCEDURE PARAMETERS'//6X,'ECMAX  =',F8.2,5X,'Maximum energy for the ','eigensolutions'/6X, &
               'BTOL   =',D8.2,5X,'Iteration starting tolerance'/6X,'TINY   =',D8.2,5X,'Eigensolution convergence parameter')
    70  format(///11X,'***  INPUT DATA ERROR - JOB ABORTED  ***'//6X,'ERROR IN THE INTEGRATION MESH SPECIFICATION'/6X, &
               'THE NUMBER OF STEPS IN ONE OR MORE REGIONS IS NOT ','EXACTLY DIVISIBLE BY 2')
    80  format(///11X,'***  INPUT DATA ERROR - JOB ABORTED  ***'//6X,'ERROR IN THE INTEGRATION MESH SPECIFICATION'/6X,'NIX =', &
               I3,'; MUST BE <= ',I3)
    90  format(///11X,'***  INPUT DATA ERROR - JOB ABORTED  ***'//6X,'ERROR IN IRX AND/OR HRX STARTING AT'/6X,'COMPONENT ',I2, &
               ' - SOME VALUES ZERO')
    end subroutine evaluate_grid


    !> \brief   Evaluate a Slater-type orbital
    !> \authors J Benda
    !> \date    2018
    !>
    !> Given the STO exponent and quantum numbers (n,l), evaluate it on the radial grid reaching from
    !> the origin to the limiting value \c rlim. Only the radial part (i.e. excluding the spherical harmonic
    !> \f$ Y_{lm} \f$) is evaluated by this routine.
    !>
    subroutine evaluate_STO (expt, nval, lval, ira, r, orb)

        use precisn, only: wp

        integer,  intent(in) :: nval, lval
        real(wp), intent(in) :: expt
        integer,  intent(in) :: ira
        real(wp), intent(inout), allocatable :: r(:), orb(:)

        integer  :: i
        real(wp) :: norm

        ! compute normalization factor containing the factorial
        norm = 1
        do i = 1, 2*nval
            norm = norm * i
        end do
        norm = sqrt((2*expt)**(2*nval + 1) / norm)

        ! evaluate the radial part of the orbital at all distances
        do i = 1, ira + 1
            orb(i) = norm * r(i)**(nval - 1) * exp(-expt*r(i))
        end do

    end subroutine evaluate_STO

end program numsbas
