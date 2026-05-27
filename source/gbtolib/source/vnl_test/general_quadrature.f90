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
module general_quadrature
use precisn
use utils, only: xermsg, d1mach

 !> \class <bound_user_function>
 !> This is a class that defines an abstract function, whose specific implementation is deferred.
 !> The purpose is to use this abstract function in some bspline-related routines which require a user-defined function as a parameters.
 type, abstract :: bound_user_function
   !no data components
 contains
      !> \memberof bound_user_function
      procedure(user_function_interface), deferred :: eval
 end type bound_user_function
 abstract interface
      real(wp) function user_function_interface(data,x)
         import :: bound_user_function, wp
         class(bound_user_function) :: data
         real(kind=wp), intent(in) :: x
      end function user_function_interface
 end interface

!DQAGS is the routine that performs the numerical quadrature
public DQAGS

!auxiliary routines for DQAGS
private DQAGSE, DQELG, DQK21, DQPSRT

contains

!>***BEGIN PROLOGUE  DQAGS
!>***PURPOSE  The routine calculates an approximation result to a given
!>            Definite integral  I = Integral of F over (A,B),
!>            Hopefully satisfying following claim for accuracy
!>            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A1
!>***TYPE      real(kind=wp) (QAGS-S, DQAGS-D)
!>***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
!>             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
!>             QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>        Computation of a definite integral
!>        Standard fortran subroutine
!>        Double precision version
!>
!>
!>        PARAMETERS
!>         ON ENTRY
!>            F      - class(bound_user_function)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>            A      - Double precision
!>                     Lower limit of integration
!>
!>            B      - Double precision
!>                     Upper limit of integration
!>
!>            EPSABS - Double precision
!>                     Absolute accuracy requested
!>            EPSREL - Double precision
!>                     Relative accuracy requested
!>                     If  EPSABS.LE.0
!>                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
!>                     The routine will end with IER = 6.
!>
!>         ON RETURN
!>            RESULT - Double precision
!>                     Approximation to the integral
!>
!>            ABSERR - Double precision
!>                     Estimate of the modulus of the absolute error,
!>                     which should equal or exceed ABS(I-RESULT)
!>
!>            NEVAL  - Integer
!>                     Number of integrand evaluations
!>
!>            IER    - Integer
!>                     IER = 0 Normal and reliable termination of the
!>                             routine. It is assumed that the requested
!>                             accuracy has been achieved.
!>                     IER.GT.0 Abnormal termination of the routine
!>                             The estimates for integral and error are
!>                             less reliable. It is assumed that the
!>                             requested accuracy has not been achieved.
!>            ERROR MESSAGES
!>                     IER = 1 Maximum number of subdivisions allowed
!>                             has been achieved. One can allow more sub-
!>                             divisions by increasing the value of LIMIT
!>                             (and taking the according dimension
!>                             adjustments into account. However, if
!>                             this yields no improvement it is advised
!>                             to analyze the integrand in order to
!>                             determine the integration difficulties. If
!>                             the position of a local difficulty can be
!>                             determined (E.G. SINGULARITY,
!>                             DISCONTINUITY WITHIN THE INTERVAL) one
!>                             will probably gain from splitting up the
!>                             interval at this point and calling the
!>                             integrator on the subranges. If possible,
!>                             an appropriate special-purpose integrator
!>                             should be used, which is designed for
!>                             handling the type of difficulty involved.
!>                         = 2 The occurrence of roundoff error is detec-
!>                             ted, which prevents the requested
!>                             tolerance from being achieved.
!>                             The error may be under-estimated.
!>                         = 3 Extremely bad integrand behaviour
!>                             occurs at some points of the integration
!>                             interval.
!>                         = 4 The algorithm does not converge.
!>                             Roundoff error is detected in the
!>                             Extrapolation table. It is presumed that
!>                             the requested tolerance cannot be
!>                             achieved, and that the returned result is
!>                             the best which can be obtained.
!>                         = 5 The integral is probably divergent, or
!>                             slowly convergent. It must be noted that
!>                             divergence can occur with any other value
!>                             of IER.
!>                         = 6 The input is invalid, because
!>                             (EPSABS.LE.0 AND
!>                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)
!>                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
!>                             RESULT, ABSERR, NEVAL, LAST are set to
!>                             zero.  Except when LIMIT or LENW is
!>                             invalid, IWORK(1), WORK(LIMIT*2+1) and
!>                             WORK(LIMIT*3+1) are set to zero, WORK(1)
!>                             is set to A and WORK(LIMIT+1) TO B.
!>
!>         DIMENSIONING PARAMETERS
!>            LIMIT - Integer
!>                    DIMENSIONING PARAMETER FOR IWORK
!>                    LIMIT determines the maximum number of subintervals
!>                    in the partition of the given integration interval
!>                    (A,B), LIMIT.GE.1.
!>                    IF LIMIT.LT.1, the routine will end with IER = 6.
!>
!>            LENW  - Integer
!>                    DIMENSIONING PARAMETER FOR WORK
!>                    LENW must be at least LIMIT*4.
!>                    If LENW.LT.LIMIT*4, the routine will end
!>                    with IER = 6.
!>
!>            LAST  - Integer
!>                    On return, LAST equals the number of subintervals
!>                    produced in the subdivision process, determines the
!>                    number of significant elements actually in the WORK
!>                    Arrays.
!>
!>         WORK ARRAYS
!>            IWORK - Integer
!>                    Vector of dimension at least LIMIT, the first K
!>                    elements of which contain pointers
!>                    to the error estimates over the subintervals
!>                    such that WORK(LIMIT*3+IWORK(1)),... ,
!>                    WORK(LIMIT*3+IWORK(K)) form a decreasing
!>                    sequence, with K = LAST IF LAST.LE.(LIMIT/2+2),
!>                    and K = LIMIT+1-LAST otherwise
!>
!>            WORK  - Double precision
!>                    Vector of dimension at least LENW
!>                    on return
!>                    WORK(1), ..., WORK(LAST) contain the left
!>                     end-points of the subintervals in the
!>                     partition of (A,B),
!>                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
!>                     the right end-points,
!>                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
!>                     the integral approximations over the subintervals,
!>                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
!>                     contain the error estimates.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  DQAGSE, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>***END PROLOGUE  DQAGS
!>
!>
      SUBROUTINE DQAGS (F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      class(bound_user_function) :: F
      real(kind=wp) A,ABSERR,B,EPSABS,EPSREL,RESULT,WORK!,F
      INTEGER IER,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
!
      DIMENSION IWORK(*),WORK(*)
!
      !DECLARE F AS: real(kind=wp), EXTERNAL :: F
      !EXTERNAL F
!
!         CHECK VALIDITY OF LIMIT AND LENW.
!
!***FIRST EXECUTABLE STATEMENT  DQAGS
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
!
!         PREPARE CALL FOR DQAGSE.
!
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
!
      CALL DQAGSE(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
!
!         CALL ERROR HANDLER IF NECESSARY.
!
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAGS', 'ABNORMAL RETURN', IER, LVL)
      RETURN
      END SUBROUTINE DQAGS

!>***BEGIN PROLOGUE  DQAGSE
!>***PURPOSE  The routine calculates an approximation result to a given
!>            definite integral I = Integral of F over (A,B),
!>            hopefully satisfying following claim for accuracy
!>            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A1
!>***TYPE      real(kind=wp) (QAGSE-S, DQAGSE-D)
!>***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
!>             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
!>             QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>        Computation of a definite integral
!>        Standard fortran subroutine
!>        Double precision version
!>
!>        PARAMETERS
!>         ON ENTRY
!>            F      - class(bound_user_function)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>            A      - Double precision
!>                     Lower limit of integration
!>
!>            B      - Double precision
!>                     Upper limit of integration
!>
!>            EPSABS - Double precision
!>                     Absolute accuracy requested
!>            EPSREL - Double precision
!>                     Relative accuracy requested
!>                     If  EPSABS.LE.0
!>                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
!>                     the routine will end with IER = 6.
!>
!>            LIMIT  - Integer
!>                     Gives an upper bound on the number of subintervals
!>                     in the partition of (A,B)
!>
!>         ON RETURN
!>            RESULT - Double precision
!>                     Approximation to the integral
!>
!>            ABSERR - Double precision
!>                     Estimate of the modulus of the absolute error,
!>                     which should equal or exceed ABS(I-RESULT)
!>
!>            NEVAL  - Integer
!>                     Number of integrand evaluations
!>
!>            IER    - Integer
!>                     IER = 0 Normal and reliable termination of the
!>                             routine. It is assumed that the requested
!>                             accuracy has been achieved.
!>                     IER.GT.0 Abnormal termination of the routine
!>                             the estimates for integral and error are
!>                             less reliable. It is assumed that the
!>                             requested accuracy has not been achieved.
!>            ERROR MESSAGES
!>                         = 1 Maximum number of subdivisions allowed
!>                             has been achieved. One can allow more sub-
!>                             divisions by increasing the value of LIMIT
!>                             (and taking the according dimension
!>                             adjustments into account). However, if
!>                             this yields no improvement it is advised
!>                             to analyze the integrand in order to
!>                             determine the integration difficulties. If
!>                             the position of a local difficulty can be
!>                             determined (e.g. singularity,
!>                             discontinuity within the interval) one
!>                             will probably gain from splitting up the
!>                             interval at this point and calling the
!>                             integrator on the subranges. If possible,
!>                             an appropriate special-purpose integrator
!>                             should be used, which is designed for
!>                             handling the type of difficulty involved.
!>                         = 2 The occurrence of roundoff error is detec-
!>                             ted, which prevents the requested
!>                             tolerance from being achieved.
!>                             The error may be under-estimated.
!>                         = 3 Extremely bad integrand behaviour
!>                             occurs at some points of the integration
!>                             interval.
!>                         = 4 The algorithm does not converge.
!>                             Roundoff error is detected in the
!>                             extrapolation table.
!>                             It is presumed that the requested
!>                             tolerance cannot be achieved, and that the
!>                             returned result is the best which can be
!>                             obtained.
!>                         = 5 The integral is probably divergent, or
!>                             slowly convergent. It must be noted that
!>                             divergence can occur with any other value
!>                             of IER.
!>                         = 6 The input is invalid, because
!>                             EPSABS.LE.0 and
!>                             EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28).
!>                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
!>                             IORD(1) and ELIST(1) are set to zero.
!>                             ALIST(1) and BLIST(1) are set to A and B
!>                             respectively.
!>
!>            ALIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the left end points
!>                     of the subintervals in the partition of the
!>                     given integration range (A,B)
!>
!>            BLIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the right end points
!>                     of the subintervals in the partition of the given
!>                     integration range (A,B)
!>
!>            RLIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the integral
!>                     approximations on the subintervals
!>
!>            ELIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the moduli of the
!>                     absolute error estimates on the subintervals
!>
!>            IORD   - Integer
!>                     Vector of dimension at least LIMIT, the first K
!>                     elements of which are pointers to the
!>                     error estimates over the subintervals,
!>                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
!>                     form a decreasing sequence, with K = LAST
!>                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
!>                     otherwise
!>
!>            LAST   - Integer
!>                     Number of subintervals actually produced in the
!>                     subdivision process
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  D1MACH, DQELG, DQK21, DQPSRT
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>***END PROLOGUE  DQAGSE
!>
!>            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
!>            LIMEXP IN SUBROUTINE DQELG (RLIST2 SHOULD BE OF DIMENSION
!>            (LIMEXP+2) AT LEAST).
!>
!>            LIST OF MAJOR VARIABLES
!>            -----------------------
!>
!>           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
!>                       CONSIDERED UP TO NOW
!>           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
!>                       CONSIDERED UP TO NOW
!>           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
!>                       (ALIST(I),BLIST(I))
!>           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2 CONTAINING
!>                       THE PART OF THE EPSILON TABLE WHICH IS STILL
!>                       NEEDED FOR FURTHER COMPUTATIONS
!>           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
!>           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
!>                       ESTIMATE
!>           ERRMAX    - ELIST(MAXERR)
!>           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
!>                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
!>           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
!>           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
!>           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
!>                       ABS(RESULT))
!>           *****1    - VARIABLE FOR THE LEFT INTERVAL
!>           *****2    - VARIABLE FOR THE RIGHT INTERVAL
!>           LAST      - INDEX FOR SUBDIVISION
!>           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
!>           NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
!>                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
!>                       INTEGRAL HAS BEEN OBTAINED IT IS PUT IN
!>                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
!>                       BY ONE.
!>           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
!>                       TO NOW, MULTIPLIED BY 1.5
!>           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
!>                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
!>           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS
!>                       ATTEMPTING TO PERFORM EXTRAPOLATION I.E. BEFORE
!>                       SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO
!>                       DECREASE THE VALUE OF ERLARG.
!>           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
!>                       IS NO LONGER ALLOWED (TRUE VALUE)
!>
!>            MACHINE DEPENDENT CONSTANTS
!>            ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!>           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
      SUBROUTINE DQAGSE (F, A, B, EPSABS, EPSREL, LIMIT, RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
      class(bound_user_function) :: F
      real(kind=wp) A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,BLIST,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2, &
     &  DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,ERRMAX, &
     &  ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,OFLOW,RESABS,RESEPS,RESULT,&
     &  RES3LA,RLIST,RLIST2,SMALL,UFLOW
      INTEGER ID,IER,IERRO,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,KSGN,KTMIN,LAST,LIMIT,MAXERR,NEVAL,NRES,NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
!
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),RES3LA(3),RLIST(*),RLIST2(52)
!
!***FIRST EXECUTABLE STATEMENT  DQAGSE
      EPMACH = D1MACH(4)
!
!            TEST ON VALIDITY OF PARAMETERS
!            ------------------------------
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IF(EPSABS.LE.0.0D+00.AND.EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28)) IER = 6
      IF(IER.EQ.6) GO TO 999
!
!           FIRST APPROXIMATION TO THE INTEGRAL
!           -----------------------------------
!
      UFLOW = D1MACH(1)
      OFLOW = D1MACH(2)
      IERRO = 0
      CALL DQK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
!
!           TEST ON ACCURACY.
!
      DRES = ABS(RESULT)
      ERRBND = MAX(EPSABS,EPSREL*DRES)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      IF(ABSERR.LE.1.0D+02*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS).OR.ABSERR.EQ.0.0D+00) GO TO 140
!
!           INITIALIZATION
!           --------------
!
      RLIST2(1) = RESULT
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      ABSERR = OFLOW
      NRMAX = 1
      NRES = 0
      NUMRL2 = 2
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KSGN = -1
      IF(DRES.GE.(0.1D+01-0.5D+02*EPMACH)*DEFABS) KSGN = 1
!
!           MAIN DO-LOOP
!           ------------
!
      DO 90 LAST = 2,LIMIT
!
!           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR
!           ESTIMATE.
!
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        CALL DQK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
!
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!           AND ERROR AND TEST FOR ACCURACY.
!
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 15
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1D-04*ABS(AREA12).OR.ERRO12.LT.0.99D+00*ERRMAX) GO TO 10
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   10   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   15   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
!
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
!
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
!
!           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
!           EQUALS LIMIT.
!
        IF(LAST.EQ.LIMIT) IER = 1
!
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT A POINT OF THE INTEGRATION RANGE.
!
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)*(ABS(A2)+0.1D+04*UFLOW)) IER = 4
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
        IF(ERROR2.GT.ERROR1) GO TO 20
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 30
   20   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
!
!           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
!           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
!
   30   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
        IF(ERRSUM.LE.ERRBND) GO TO 115
! ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0) GO TO 100
        IF(LAST.EQ.2) GO TO 80
        IF(NOEXT) GO TO 90
        ERLARG = ERLARG-ERLAST
        IF(ABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 40
!
!           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
!           SMALLEST INTERVAL.
!
        IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
        EXTRAP = .TRUE.
        NRMAX = 2
   40   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 60
!
!           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
!           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
!           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
!
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 50 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
! ***JUMP OUT OF DO-LOOP
          IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
          NRMAX = NRMAX+1
   50   CONTINUE
!
!           PERFORM EXTRAPOLATION.
!
   60   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        CALL DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1D-02*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 70
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
! ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LE.ERTEST) GO TO 100
!
!           PREPARE BISECTION OF THE SMALLEST INTERVAL.
!
   70   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 100
        MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5D+00
        ERLARG = ERRSUM
        GO TO 90
   80   SMALL = ABS(B-A)*0.375D+00
        ERLARG = ERRSUM
        ERTEST = ERRBND
        RLIST2(2) = AREA
   90 CONTINUE
!
!           SET FINAL RESULT AND ERROR ESTIMATE.
!           ------------------------------------
!
  100 IF(ABSERR.EQ.OFLOW) GO TO 115
      IF(IER+IERRO.EQ.0) GO TO 110
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0D+00.AND.AREA.NE.0.0D+00) GO TO 105
      IF(ABSERR.GT.ERRSUM) GO TO 115
      IF(AREA.EQ.0.0D+00) GO TO 130
      GO TO 110
  105 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA)) GO TO 115
!
!           TEST ON DIVERGENCE.
!
  110 IF(KSGN.EQ.(-1).AND.MAX(ABS(RESULT),ABS(AREA)).LE.DEFABS*0.1D-01) GO TO 130
      IF(0.1D-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1D+03.OR.ERRSUM.GT.ABS(AREA)) IER = 6
      GO TO 130
!
!           COMPUTE GLOBAL INTEGRAL SUM.
!
  115 RESULT = 0.0D+00
      DO 120 K = 1,LAST
         RESULT = RESULT+RLIST(K)
  120 CONTINUE
      ABSERR = ERRSUM
  130 IF(IER.GT.2) IER = IER-1
  140 NEVAL = 42*LAST-21
  999 RETURN
      END SUBROUTINE DQAGSE

!>***BEGIN PROLOGUE  DQELG
!>***SUBSIDIARY
!>***PURPOSE  The routine determines the limit of a given sequence of
!>            approximations, by means of the Epsilon algorithm of
!>            P.Wynn. An estimate of the absolute error is also given.
!>            The condensed Epsilon table is computed. Only those
!>            elements needed for the computation of the next diagonal
!>            are preserved.
!>***LIBRARY   SLATEC
!>***TYPE      real(kind=wp) (QELG-S, DQELG-D)
!>***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Epsilon algorithm
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS
!>              N      - Integer
!>                       EPSTAB(N) contains the new element in the
!>                       first column of the epsilon table.
!>
!>              EPSTAB - Double precision
!>                       Vector of dimension 52 containing the elements
!>                       of the two lower diagonals of the triangular
!>                       epsilon table. The elements are numbered
!>                       starting at the right-hand corner of the
!>                       triangle.
!>
!>              RESULT - Double precision
!>                       Resulting approximation to the integral
!>
!>              ABSERR - Double precision
!>                       Estimate of the absolute error computed from
!>                       RESULT and the 3 previous results
!>
!>              RES3LA - Double precision
!>                       Vector of dimension 3 containing the last 3
!>                       results
!>
!>              NRES   - Integer
!>                       Number of calls to the routine
!>                       (should be zero at first call)
!>
!>***SEE ALSO  DQAGIE, DQAGOE, DQAGPE, DQAGSE
!>***ROUTINES CALLED  D1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>***END PROLOGUE  DQELG
!>
!>           LIST OF MAJOR VARIABLES
!>           -----------------------
!>
!>           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
!>           E1       ELEMENT IN THE EPSILON TABLE IS BASED
!>           E2
!>           E3                 E0
!>                        E3    E1    NEW
!>                              E2
!>           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
!>                    DIAGONAL
!>           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
!>           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
!>                    OF ERROR
!>
!>           MACHINE DEPENDENT CONSTANTS
!>           ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
!>           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
!>           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
!>           DIAGONAL OF THE EPSILON TABLE IS DELETED.
!>
      SUBROUTINE DQELG (N, EPSTAB, RESULT, ABSERR, RES3LA, NRES)
      real(kind=wp) ABSERR,DELTA1,DELTA2,DELTA3,EPMACH,EPSINF,EPSTAB,ERROR,ERR1,ERR2,ERR3,E0,E1,E1ABS,E2,E3,OFLOW,RES,RESULT,RES3LA,SS,TOL1,TOL2,TOL3
      INTEGER I,IB,IB2,IE,INDX,K1,K2,K3,LIMEXP,N,NEWELM,NRES,NUM
      DIMENSION EPSTAB(52),RES3LA(3)
!***FIRST EXECUTABLE STATEMENT  DQELG
      EPMACH = D1MACH(4)
      OFLOW = D1MACH(2)
      NRES = NRES+1
      ABSERR = OFLOW
      RESULT = EPSTAB(N)
      IF(N.LT.3) GO TO 100
      LIMEXP = 50
      EPSTAB(N+2) = EPSTAB(N)
      NEWELM = (N-1)/2
      EPSTAB(N) = OFLOW
      NUM = N
      K1 = N
      DO 40 I = 1,NEWELM
        K2 = K1-1
        K3 = K1-2
        RES = EPSTAB(K1+2)
        E0 = EPSTAB(K3)
        E1 = EPSTAB(K2)
        E2 = RES
        E1ABS = ABS(E1)
        DELTA2 = E2-E1
        ERR2 = ABS(DELTA2)
        TOL2 = MAX(ABS(E2),E1ABS)*EPMACH
        DELTA3 = E1-E0
        ERR3 = ABS(DELTA3)
        TOL3 = MAX(E1ABS,ABS(E0))*EPMACH
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
!
!           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
!           ACCURACY, CONVERGENCE IS ASSUMED.
!           RESULT = E2
!           ABSERR = ABS(E1-E0)+ABS(E2-E1)
!
        RESULT = RES
        ABSERR = ERR2+ERR3
! ***JUMP OUT OF DO-LOOP
        GO TO 100
   10   E3 = EPSTAB(K1)
        EPSTAB(K1) = E1
        DELTA1 = E1-E3
        ERR1 = ABS(DELTA1)
        TOL1 = MAX(E1ABS,ABS(E3))*EPMACH
!
!           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
!           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
!
        IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
        SS = 0.1D+01/DELTA1+0.1D+01/DELTA2-0.1D+01/DELTA3
        EPSINF = ABS(SS*E1)
!
!           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
!           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
!           OF N.
!
        IF(EPSINF.GT.0.1D-03) GO TO 30
   20   N = I+I-1
! ***JUMP OUT OF DO-LOOP
        GO TO 50
!
!           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
!           THE VALUE OF RESULT.
!
   30   RES = E1+0.1D+01/SS
        EPSTAB(K1) = RES
        K1 = K1-2
        ERROR = ERR2+ABS(RES-E2)+ERR3
        IF(ERROR.GT.ABSERR) GO TO 40
        ABSERR = ERROR
        RESULT = RES
   40 CONTINUE
!
!           SHIFT THE TABLE.
!
   50 IF(N.EQ.LIMEXP) N = 2*(LIMEXP/2)-1
      IB = 1
      IF((NUM/2)*2.EQ.NUM) IB = 2
      IE = NEWELM+1
      DO 60 I=1,IE
        IB2 = IB+2
        EPSTAB(IB) = EPSTAB(IB2)
        IB = IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      INDX = NUM-N+1
      DO 70 I = 1,N
        EPSTAB(I)= EPSTAB(INDX)
        INDX = INDX+1
   70 CONTINUE
   80 IF(NRES.GE.4) GO TO 90
      RES3LA(NRES) = RESULT
      ABSERR = OFLOW
      GO TO 100
!
!           COMPUTE ERROR ESTIMATE
!
   90 ABSERR = ABS(RESULT-RES3LA(3))+ABS(RESULT-RES3LA(2))+ABS(RESULT-RES3LA(1))
      RES3LA(1) = RES3LA(2)
      RES3LA(2) = RES3LA(3)
      RES3LA(3) = RESULT
  100 ABSERR = MAX(ABSERR,0.5D+01*EPMACH*ABS(RESULT))
      RETURN
      END SUBROUTINE DQELG

!>***BEGIN PROLOGUE  DQK21
!>***PURPOSE  To compute I = Integral of F over (A,B), with error
!>                           estimate
!>                       J = Integral of ABS(F) over (A,B)
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A2
!>***TYPE      real(kind=wp) (QK21-S, DQK21-D)
!>***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Integration rules
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS
!>            ON ENTRY
!>            F      - class(bound_user_function)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>              A      - Double precision
!>                       Lower limit of integration
!>
!>              B      - Double precision
!>                       Upper limit of integration
!>
!>            ON RETURN
!>              RESULT - Double precision
!>                       Approximation to the integral I
!>                       RESULT is computed by applying the 21-POINT
!>                       KRONROD RULE (RESK) obtained by optimal addition
!>                       of abscissae to the 10-POINT GAUSS RULE (RESG).
!>
!>              ABSERR - Double precision
!>                       Estimate of the modulus of the absolute error,
!>                       which should not exceed ABS(I-RESULT)
!>
!>              RESABS - Double precision
!>                       Approximation to the integral J
!>
!>              RESASC - Double precision
!>                       Approximation to the integral of ABS(F-I/(B-A))
!>                       over (A,B)
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  D1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>***END PROLOGUE  DQK21
!>
!>
!>           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!>           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!>           CORRESPONDING WEIGHTS ARE GIVEN.
!>
!>           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
!>                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
!>                    GAUSS RULE
!>                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!>                    ADDED TO THE 10-POINT GAUSS RULE
!>
!>           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
!>
!>           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
!>
!>
!> GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
!> AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
!> BELL LABS, NOV. 1981.
!>
!>
!>           LIST OF MAJOR VARIABLES
!>           -----------------------
!>
!>           CENTR  - MID POINT OF THE INTERVAL
!>           HLGTH  - HALF-LENGTH OF THE INTERVAL
!>           ABSC   - ABSCISSA
!>           FVAL*  - FUNCTION VALUE
!>           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
!>           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
!>           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!>                    I.E. TO I/(B-A)
!>
!>
!>           MACHINE DEPENDENT CONSTANTS
!>           ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!>
      SUBROUTINE DQK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
      class(bound_user_function) :: F
      real(kind=wp) A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
!
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
!
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.066671344308688137593568809893332D0 /
      DATA WG  (  2) / 0.149451349150580593145776339657697D0 /
      DATA WG  (  3) / 0.219086362515982043995534934228163D0 /
      DATA WG  (  4) / 0.269266719309996355091226921569469D0 /
      DATA WG  (  5) / 0.295524224714752870173892994651338D0 /
!
      DATA XGK (  1) / 0.995657163025808080735527280689003D0 /
      DATA XGK (  2) / 0.973906528517171720077964012084452D0 /
      DATA XGK (  3) / 0.930157491355708226001207180059508D0 /
      DATA XGK (  4) / 0.865063366688984510732096688423493D0 /
      DATA XGK (  5) / 0.780817726586416897063717578345042D0 /
      DATA XGK (  6) / 0.679409568299024406234327365114874D0 /
      DATA XGK (  7) / 0.562757134668604683339000099272694D0 /
      DATA XGK (  8) / 0.433395394129247190799265943165784D0 /
      DATA XGK (  9) / 0.294392862701460198131126603103866D0 /
      DATA XGK ( 10) / 0.148874338981631210884826001129720D0 /
      DATA XGK ( 11) / 0.000000000000000000000000000000000D0 /
!
      DATA WGK (  1) / 0.011694638867371874278064396062192D0 /
      DATA WGK (  2) / 0.032558162307964727478818972459390D0 /
      DATA WGK (  3) / 0.054755896574351996031381300244580D0 /
      DATA WGK (  4) / 0.075039674810919952767043140916190D0 /
      DATA WGK (  5) / 0.093125454583697605535065465083366D0 /
      DATA WGK (  6) / 0.109387158802297641899210590325805D0 /
      DATA WGK (  7) / 0.123491976262065851077958109831074D0 /
      DATA WGK (  8) / 0.134709217311473325928054001771707D0 /
      DATA WGK (  9) / 0.142775938577060080797094273138717D0 /
      DATA WGK ( 10) / 0.147739104901338491374841515972068D0 /
      DATA WGK ( 11) / 0.149445554002916905664936468389821D0 /
!
!***FIRST EXECUTABLE STATEMENT  DQK21
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
!
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
      RESG = 0.0D+00
      FC = F%eval(CENTR)
      RESK = WGK(11)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,5
        JTW = 2*J
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F%eval(CENTR-ABSC)
        FVAL2 = F%eval(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,5
        JTWM1 = 2*J-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F%eval(CENTR-ABSC)
        FVAL2 = F%eval(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(11)*ABS(FC-RESKH)
      DO 20 J=1,10
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00) ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END SUBROUTINE DQK21

!>***BEGIN PROLOGUE  DQPSRT
!>***SUBSIDIARY
!>***PURPOSE  This routine maintains the descending ordering in the
!>            list of the local error estimated resulting from the
!>            interval subdivision process. At each call two error
!>            estimates are inserted using the sequential search
!>            method, top-down for the largest error estimate and
!>            bottom-up for the smallest error estimate.
!>***LIBRARY   SLATEC
!>***TYPE      real(kind=wp) (QPSRT-S, DQPSRT-D)
!>***KEYWORDS  SEQUENTIAL SORTING
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Ordering routine
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS (MEANING AT OUTPUT)
!>              LIMIT  - Integer
!>                       Maximum number of error estimates the list
!>                       can contain
!>
!>              LAST   - Integer
!>                       Number of error estimates currently in the list
!>
!>              MAXERR - Integer
!>                       MAXERR points to the NRMAX-th largest error
!>                       estimate currently in the list
!>
!>              ERMAX  - Double precision
!>                       NRMAX-th largest error estimate
!>                       ERMAX = ELIST(MAXERR)
!>
!>              ELIST  - Double precision
!>                       Vector of dimension LAST containing
!>                       the error estimates
!>
!>              IORD   - Integer
!>                       Vector of dimension LAST, the first K elements
!>                       of which contain pointers to the error
!>                       estimates, such that
!>                       ELIST(IORD(1)),...,  ELIST(IORD(K))
!>                       form a decreasing sequence, with
!>                       K = LAST if LAST.LE.(LIMIT/2+2), and
!>                       K = LIMIT+1-LAST otherwise
!>
!>              NRMAX  - Integer
!>                       MAXERR = IORD(NRMAX)
!>
!>***SEE ALSO  DQAGE, DQAGIE, DQAGPE, DQAWSE
!>***ROUTINES CALLED  (NONE)
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>***END PROLOGUE  DQPSRT
!>
      SUBROUTINE DQPSRT (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)
      real(kind=wp) ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,NRMAX
      DIMENSION ELIST(*),IORD(*)
!
!           CHECK WHETHER THE LIST CONTAINS MORE THAN
!           TWO ERROR ESTIMATES.
!
!***FIRST EXECUTABLE STATEMENT  DQPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
!
!           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
!           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
!           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
!           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
!
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
!
!           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
!           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
!           SUBDIVISIONS STILL ALLOWED.
!
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
!
!           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
!           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
!
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
!
!           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
!
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
!
!           SET MAXERR AND ERMAX.
!
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END SUBROUTINE DQPSRT

end module general_quadrature
