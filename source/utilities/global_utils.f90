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
module global_utils
!
! Module containing global utility functions, and a parameter used in 
! one of the subroutines. 
! Joanne Carr, February 2010
!
  USE precisn, ONLY : wp ! for specifying the kind of reals
  implicit none
  private ! for now, for security

! Number of columns of vectors per line written out in subroutine PRINTV
  INTEGER, PARAMETER :: NLINE=4

! Number of columns written per line in MATTPT and WRECMT
  INTEGER, PARAMETER :: NCOL=6

  public :: initvr8 ! initializes a w.p real vector with zeroes
  public :: s15erfc ! calculates the complement of the error function erf(x)
  public :: mattpt ! prints a 2-dim symmetric matrix stored in lower triangle form
  public :: wrecmt ! prints an n by m matrix stored in a linear array 
  public :: printv ! prints out vectors in array VEC NMO in number with NBF elements per vector
  public :: cwbopn ! opens a sequential unit
  public :: utils_date_time ! returns the current date and time
  public :: indfunc ! returns a sequential index for element I,J in the lower half packed triangle.
  public :: utils_get_bytesizes ! returns the size in bytes of integers and reals in the compiled code
  public :: search ! search a file for a given header in Molecule-Sweden format
  public :: intape ! read a vector A of length IX from unit IU
  public :: mprod ! obtain and optionally print the group multiplication table for D2h
  public :: getin ! returns the character equivalent of the ND digit integer N
  public :: copyd ! copies the headers from one integrals file to another
  public :: print_ukrmol_header ! prints header at start-up of programs

! Laura Moore, May 2011
  public :: init_2d_wp_array_to_zero ! initializes a w.p. real 2D array with zeroes

  contains

!*==initvr8.spg  processed by SPAG 6.56Rc at 13:47 on 27 Jan 2010
      SUBROUTINE INITVR8(VEC,N)
!***********************************************************************
!
!     Initializes the components of a vector, VEC, of length N to
!     zero. This is for w.p. vectors only.
!
!***********************************************************************
      USE consts, ONLY : XZERO
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER, INTENT(IN) :: N ! dimension of the vector
      REAL(KIND=wp), DIMENSION(N), INTENT(OUT) :: VEC ! the vector 
!
! Local variables
!
      INTEGER :: I
!
!*** End of declarations rewritten by SPAG
!
      DO I=1, N
         VEC(I)=XZERO
      END DO
 
      RETURN
      END SUBROUTINE INITVR8

! Laura Moore, May 2011            
      SUBROUTINE INIT_2D_WP_ARRAY_TO_ZERO(ARRAY,M,N)
!***********************************************************************
!
!     Initializes the components of a 2D array, Array, of dim M x N to
!     zero. This is for w.p. arrays only.
!
!***********************************************************************
      USE consts, ONLY : XZERO
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER, INTENT(IN) :: M, N ! dimensions of the array
      REAL(KIND=wp), DIMENSION(M,N), INTENT(OUT) :: ARRAY ! the vector 
!
! Local variables
!
      INTEGER :: I, J
!
!*** End of declarations rewritten by SPAG
!
      DO I=1, N
        DO J=1, M
          ARRAY(J,I)=XZERO
        END DO
      END DO
 
      RETURN
      END SUBROUTINE INIT_2D_WP_ARRAY_TO_ZERO

!*==mattpt.spg  processed by SPAG 6.56Rc at 13:47 on 27 Jan 2010
      SUBROUTINE MATTPT(N,A,IPRINT,ZTEST)
!***********************************************************************
!
!     PRINTS 2-DIMENSIONAL SYMMETIC MATRIX OF DIMENSION N, THE LOWER
!     TRIANGLE OF WHICH IS STORED IN A
!     optionally modified so that zero lines are not printed
!     JMC (compare swedmos.f and gaustail.f)
!
!***********************************************************************
!
      USE consts, ONLY : XZERO
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER, INTENT(IN) :: IPRINT, N
      REAL(KIND=wp), DIMENSION(*), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: ZTEST
!
! Local variables
!
      INTEGER :: J, K, K1, K2, KL, KL1, KL2, L
!
!*** End of declarations rewritten by SPAG
!
      K1=1
      K2=NCOL
!
      DO ! jmc to replace a GO TO arrangement
         K=0
         DO L=1, N
            KL=K+L
            KL1=K+K1
            KL2=MIN(K+K2,KL)
            IF(ZTEST)THEN
               DO J=KL1, KL2
                  IF(A(J).NE.xzero)GO TO 7
               END DO
            END IF
            K=KL
         END DO
         IF(ZTEST)GO TO 2
!
! --- Write column labels
 7       WRITE(IPRINT,1002)(j,j=k1,k2)
!
! --- Write data
         K=0
         DO L=1, N
            KL=K+L
            KL1=K+K1
            KL2=MIN(K+K2,KL)
            IF(KL1.LE.KL2)WRITE(IPRINT,1001)l, (A(J),J=KL1,KL2)
            K=KL
         END DO
!
! --- Prepare new range
! jmc         IF(N.LE.K2)RETURN
 2       IF(N.LE.K2)EXIT ! jmc to replace a GO TO
         K1=K2+1
         K2=K2+NCOL
      END DO ! jmc to replace a GO TO
      RETURN
!
!     FORMAT STATEMENTS -
!
 1001 FORMAT(i5,2x,6E12.4)
 1002 FORMAT(/2x,6I12)
!
      END SUBROUTINE MATTPT

!*==wrecmt.spg  processed by SPAG 6.56Rc at 13:47 on 27 Jan 2010
      SUBROUTINE WRECMT(A,N,M,NN,MM,IWRITE)
!
!***********************************************************************
!
!     WRECMT PRINTS OUT A N*M MATRIX STORED IN A NN * MM ARRAY
!     modified so that zero lines are not printed
!
!***********************************************************************
!
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: IWRITE, M, MM, N, NN
      REAL(KIND=wp), DIMENSION(NN,MM) :: A
      INTENT (IN) A, IWRITE, M, MM, N, NN
!
! Local variables
!
      INTEGER :: I, J, K, NF, NI, NTIM
!
!*** End of declarations rewritten by SPAG
!
 1000 FORMAT(i5,2x,6D12.4)
 1010 FORMAT(/2x,6I12)
!
      NTIM=M/NCOL
      NF=0
! jmc      IF(NTIM.EQ.0)GO TO 30
! jmc the above statement is superfluous, as the contents of the do 
! jmc loop below will not be executed if upper_bound < lower_bound 
! jmc (with a +ve increment)
!
      DO I=1, NTIM
         NI=NF+1
         NF=NF+NCOL
         WRITE(IWRITE,1010)(k,k=ni,nf)
         DO J=1, N
            WRITE(IWRITE,1000)j, (A(J,K),K=NI,NF)
         END DO
      END DO
!
 30   NI=NF+1
      IF(NI.GT.M)RETURN
      WRITE(IWRITE,1010)(k,k=ni,m)
      DO J=1, N
         WRITE(IWRITE,1000)j, (A(J,K),K=NI,M)
      END DO
!
      RETURN
      END SUBROUTINE WRECMT

!*==printv.spg  processed by SPAG 6.56Rc at 14:46 on 23 Feb 2010
      SUBROUTINE PRINTV(VEC,NBF,NMO,ISYM,IWRITE)
!********************************************************************
!
!     PRINT - Prints out vectors in array VEC NMO in number with
!             NBF elements per vector
!
!********************************************************************
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: ISYM, NBF, NMO, IWRITE
      REAL(KIND=wp), DIMENSION(*) :: VEC
      INTENT (IN) ISYM, NBF, NMO, VEC, IWRITE
!
! Local variables
!
      INTEGER :: I, J, K, NDONE, NLEFT, NOW, NPASS, NPRINT
!
!*** End of declarations rewritten by SPAG
!
      NDONE=0
      NPASS=(NMO+NLINE-1)/NLINE
      NLEFT=NMO
      NOW=0
!
      WRITE(IWRITE,1)ISYM
!
      DO I=1, NPASS
         NPRINT=MIN(NLEFT,NLINE)
         WRITE(IWRITE,99)(NOW+K,K=1,NPRINT)
         DO J=1, NBF
            WRITE(IWRITE,100)J, (VEC(J+NDONE+(K-1)*NBF),K=1,NPRINT)
         END DO
         NOW=NOW+NLINE
         NDONE=NDONE+NPRINT*NBF
         NLEFT=NLEFT-NLINE
      END DO
!
      RETURN
!
 1    FORMAT(//' VECTORS FOR SYMMETRY',I4)
 99   FORMAT(/4X,I9,9I25)
 100  FORMAT(1X,I3,10(2X,E25.15))
!
      END SUBROUTINE PRINTV

!*==cwbopn.spg  processed by SPAG 6.56Rc at 14:46 on 23 Feb 2010
!
      SUBROUTINE CWBOPN(LUNIT)
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: LUNIT
      INTENT (IN) LUNIT
!
! Local variables
!
      LOGICAL :: OP
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     OPENS UP A SEQUENTIAL UNIT
!
!***********************************************************************
!
!
      INQUIRE(UNIT=LUNIT,OPENED=OP)
!
      IF(.NOT.OP)OPEN(UNIT=LUNIT,STATUS='UNKNOWN',FORM='UNFORMATTED',&
     &                ACCESS='SEQUENTIAL')
!
      REWIND LUNIT
!
      RETURN
      END SUBROUTINE CWBOPN

      SUBROUTINE UTILS_DATE_TIME(CURRENT_DATE, CURRENT_TIME)
      IMPLICIT NONE

      CHARACTER(LEN=8), OPTIONAL :: CURRENT_DATE
      CHARACTER(LEN=10), OPTIONAL :: CURRENT_TIME
 
      IF (PRESENT(CURRENT_DATE)) THEN
         IF (PRESENT(CURRENT_TIME)) THEN
            CALL DATE_AND_TIME(CURRENT_DATE, CURRENT_TIME)
         ELSE
            CALL DATE_AND_TIME(CURRENT_DATE)
         END IF
      ELSE
         IF (PRESENT(CURRENT_TIME)) THEN
            CALL DATE_AND_TIME(TIME=CURRENT_TIME)
         END IF
      END IF

      RETURN
      END SUBROUTINE UTILS_DATE_TIME

      SUBROUTINE utils_get_bytesizes(linteg, lreal)
      USE precisn, ONLY : rbyte ! for specifying the byte size of reals
      IMPLICIT NONE
 
! Dummy arguments
      INTEGER, INTENT(OUT) :: linteg, lreal
      
      linteg = bit_size(1) / 8 ! converting bit size of default integers to byte size
      lreal = rbyte

      RETURN
      END SUBROUTINE utils_get_bytesizes

!*==search.spg  processed by SPAG 6.56Rc at 11:57 on 13 May 2010
      SUBROUTINE SEARCH(IUNIT,A,ifail,iwrite)
!***********************************************************************
!
!     Utility to search a dataset IUNIT for a header A where the dataset
!     is assumed to have MOLECULE-SWEDEN convention headers. The header
!     convention is
!
!     '********', '        ', '        ', 'ABCDEFGH'
!
!     with ABCDEFGH being a character string such as ONEELINT etc.
!
!***********************************************************************
      USE params, ONLY : C8STARS
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      CHARACTER(LEN=8) :: A
      INTEGER :: IFAIL, IUNIT, IWRITE
      INTENT (IN) A, IUNIT, IWRITE
      INTENT (OUT) IFAIL
!
! Local variables
!
      INTEGER :: I
      CHARACTER(LEN=8), DIMENSION(4) :: B ! dimensioned according to the Sweden header format 
!
!*** End of declarations rewritten by SPAG
!
      REWIND IUNIT
 
      DO I = 1, HUGE(1) - 1 ! ensure the do loop is technically not infinite
         READ(IUNIT,END=990)B
         IF(B(1).EQ.C8STARS .AND. B(4).EQ.A)EXIT
      END DO
 
!..... We reach this point if the search has been successful
      ifail=0
      RETURN
 
!---- Process error condition namely, header not found by end of file.
 990  CONTINUE

      WRITE(IWRITE,9900)A, IUNIT
 9900 FORMAT(/' **** Error in SEARCH: ',//, &
             ' Attempt to find header (A) = ',A,' on unit = ',I3,/, &
             ' has failed.'/)
      ifail=1
      RETURN
 
      END SUBROUTINE SEARCH

!*==intape.spg  processed by SPAG 6.56Rc at 13:54 on  8 Oct 2010
      SUBROUTINE INTAPE(IU,A,IX)
!***********************************************************************
!
!     Utility routine to read a vector A of length IX from unit IU
!
!***********************************************************************
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: IU, IX
      REAL(KIND=wp), DIMENSION(IX) :: A
      INTENT (IN) IU, IX
      INTENT (INOUT) A
!
!*** End of declarations rewritten by SPAG
!
      READ(IU)A
!
      RETURN
!
      END SUBROUTINE INTAPE

!*==getin.spg  processed by SPAG 6.56Rc at 13:56 on 20 Aug 2010
      SUBROUTINE GETIN(N,ND,BCDIN,ICASE)
      IMPLICIT NONE
!
!     GETIN returns the character equivalent of the ND digit integer N
!      in the character*1 array BCDIN
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: ICASE, N, ND
      CHARACTER(LEN=1), DIMENSION(*) :: BCDIN
      INTENT (IN) ICASE, N, ND
      INTENT (OUT) BCDIN
!
! Local variables
!
      INTEGER :: I, IL, IS, NR
      CHARACTER(LEN=1), DIMENSION(10), PARAMETER :: LCHAR= &
     &         (/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
!
!*** End of declarations rewritten by SPAG
!
      IF(ICASE.EQ.1)THEN
         DO I=1, ND
            BCDIN(I)='0'
         END DO
      ELSE
         DO I=1, ND
            BCDIN(I)=' '
         END DO
      END IF
!
      IF(ND.LE.0)RETURN ! JMC rudimentary error handling
      NR=ABS(N)
!
      DO IS=ND, 1, -1
         IL=NR-(NR/10)*10
         BCDIN(IS)=LCHAR(IL+1)
         NR=NR/10
         IF(NR.EQ.0)EXIT
      END DO
      RETURN
      END SUBROUTINE GETIN

      INTEGER FUNCTION INDFUNC(I,J)
!
!       INDFUNC - Returns a sequential index for element I,J in the
!       lower half packed triangle.
!       Note: this is a little function to replace a statement function
!       (statement fns are obsolescent in Fortran95)
!
      IMPLICIT NONE

! Dummy arguments
      INTEGER, INTENT(IN) :: I, J

      INDFUNC=((I*(I-1))/2)+J ! Integer division is OK as i*(i-1) is always even

      END FUNCTION INDFUNC

!*==mprod.spg  processed by SPAG 6.56Rc at 13:56 on 20 Aug 2010
      FUNCTION MPROD(M1,M2,NPMULT,NFTW)
      USE PARAMS, ONLY : IPD2H, OS, MAXSYM
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: M1, M2, NPMULT, NFTW
      INTEGER :: MPROD
      INTENT (IN) M1, M2, NPMULT, NFTW
!
! Local variables
!
      CHARACTER(LEN=4) :: BLANK='    ', STAR12=' *  ', STAR2L='**  ', &
                          STAR4='****'
      INTEGER :: I, M, N, NCLASS
      CHARACTER(LEN=4), DIMENSION(MAXSYM+1) :: LIST
!
!*** End of declarations rewritten by SPAG
!
!     THIS FUNCTION RETURNS THE DIRECT PRODUCT OF M1 AND M2
!         PRESENT IMPLIMENTATION IS FOR D2H SYMMETRY
!     NPMULT .GT. 0 PERMITS PRINTING OF MULTIPLICATION TABLE
!
      MPROD=IPD2H(M1,M2)
      IF(NPMULT.LE.0)RETURN
!
      NCLASS=MAXSYM
      WRITE(NFTW, 15)
 15   FORMAT(//10X,'GROUP MULTIPLICATION TABLE FOR D2H SYMMETRY'/)
      WRITE(NFTW, 20) (I,I=1,NCLASS)
 20   FORMAT(16X,20I4)
      LIST(1)=BLANK
      DO I=1, NCLASS
         LIST(I+1)=OS(I)
      END DO
      WRITE(NFTW, 25) LIST
 25   FORMAT(13X,21A4)
      WRITE(NFTW, 30) (STAR4,I=1,NCLASS), STAR2L
 30   FORMAT(14X,' **',21A4)
      WRITE(NFTW, 35) (BLANK,I=1,NCLASS), STAR12
 35   FORMAT(14X,' * ',21A4)
      DO M=1, NCLASS
         LIST(1)=OS(M)
         DO N=1, NCLASS
            LIST(N+1)=OS(IPD2H(N,M))
         END DO
         WRITE(NFTW, 40) LIST, STAR12
 40      FORMAT(10X,A4,' * ',21A4)
         WRITE(NFTW, 35) (BLANK,I=1,NCLASS), STAR12
      END DO
      WRITE(NFTW, 30) (STAR4,I=1,NCLASS), STAR2L
      RETURN
      END FUNCTION MPROD

!*==s15erfc.spg  processed by SPAG 6.56Rc at 13:47 on 27 Jan 2010
!
      FUNCTION S15ERFC(X,IFAIL)
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
!     .. Scalar Arguments ..
      INTEGER :: IFAIL
      REAL(KIND=wp) :: X
      REAL(KIND=wp) :: S15ERFC
      INTENT (IN) X
      INTENT (OUT) IFAIL
!
! Local variables
!
!     .. Local Scalars ..
      REAL(KIND=wp) :: T, XHI, XLO, Y
!     .. Intrinsic Functions ..
      INTRINSIC ABS, EXP
!
!*** End of declarations rewritten by SPAG
!
!     COMPLEMENT OF ERROR FUNCTION ERFC(X)
!
!     .. Data statements ..
!     PRECISION DEPENDENT CONSTANTS
      DATA XLO/-6.25_wp/
!
!     RANGE DEPENDENT CONSTANTS
      DATA XHI/2.66E+1_wp/
!     XHI = LARGEST X SUCH THAT EXP(-X*X) .GT. MINREAL (ROUNDED DOWN)
!     .. Executable Statements ..
!
!     NO FAILURE EXITS
      IFAIL=0
!     TEST EXTREME EXITS
      IF(X.GE.XHI) THEN
         S15ERFC=0.0_wp
         RETURN
      END IF
      IF(X.LE.XLO) THEN
         S15ERFC=2.0_wp
         RETURN
      END IF
!
!     EXPANSION ARGUMENT
      T=1.0_wp-7.5_wp/(ABS(X)+3.75_wp)
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 16E
      Y = (((((((((((((((+3.328130055126039E-10_wp &
     &    *T-5.718639670776992E-10_wp)*T-4.066088879757269E-9_wp) &
     &    *T+7.532536116142436E-9_wp)*T+3.026547320064576E-8_wp) &
     &    *T-7.043998994397452E-8_wp)*T-1.822565715362025E-7_wp) &
     &    *T+6.575825478226343E-7_wp)*T+7.478317101785790E-7_wp) &
     &    *T-6.182369348098529E-6_wp)*T+3.584014089915968E-6_wp) &
     &    *T+4.789838226695987E-5_wp)*T-1.524627476123466E-4_wp) &
     &    *T-2.553523453642242E-5_wp)*T+1.802962431316418E-3_wp) &
     &    *T-8.220621168415435E-3_wp)*T+2.414322397093253E-2_wp
      Y = (((((Y*T-5.480232669380236E-2_wp)*T+1.026043120322792E-1_wp) &
     &    *T-1.635718955239687E-1_wp)*T+2.260080669166197E-1_wp) &
     &    *T-2.734219314954260E-1_wp)*T + 1.455897212750385E-1_wp
!
      S15ERFC=EXP(-X*X)*Y
      IF(X.LT.0.0_wp)S15ERFC=2.0_wp-S15ERFC
      RETURN
!
      END FUNCTION S15ERFC

!*==copyd.spg  processed by SPAG 6.56Rc at 10:41 on  8 Oct 2010
      SUBROUTINE COPYD(ITAPE,IWRITE,NFTTAIL,LSWORD)
!
!***********************************************************************
!
!     COPYD - Completely re-written version of the routine used in
!             several places throughout SWEDEN to copy integral
!             headers.This version can handle either Alchemy or
!             Cray version of MOLECULE on input but copies data to
!             Almlof's Cray format only! In essence this means that
!             character*4 data is converted to character*8.
!
!     Input data:
!          ITAPE  Logical unit for the MOLECULE output
!         IWRITE  Logical unit for the printer
!             IC  Workspace array for no. of primitives per contracted
!                 Gaussian function
!          ITRAN  Workspace array for integers
!          CTRAN  Workspace array for real*8 variables
!        NFTTAIL  Logical unit for the dataset which will hold the
!                 computed tail integrals
!
!     Notes:
!
!     This routine was developed as part of the tail integrals
!     package and then upgraded to general use. In particular it
!     should be noted that the SWEDEN package assumes that all data
!     records are at least 32 bytes long. Sometimes fillers must
!     added in this respect.
!
!     Code copies the header from the dataset of Gaussian integrals that
!     are computed by the code MOLECULE - this means that all data
!     prior to the records of integrals is copied. Note that some data
!     is returned to the caller here so that this is more than just a
!     utility routine - additionally, the records are copied onto
!     unit NFTTAIL which holds the file of integral tails, or in other
!     modules it would other things.
!
!***********************************************************************
      USE params, ONLY : CFILLERX, MAXSYM
      USE swmol3_data, ONLY : NHTT ! JMC replace this when the angular momentum parameter is moved to params
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: ITAPE, IWRITE, NFTTAIL
      LOGICAL :: LSWORD
      INTENT (IN) ITAPE, IWRITE, NFTTAIL, LSWORD
!
! Local variables
!
      CHARACTER(LEN=8) :: CL='        ', CNAKO='        ', CNAME
      REAL(KIND=wp) :: CH, CX, XX
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: CNUCNAM
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CTRAN
      INTEGER :: I, IABAS, IBBAS, II, J, K, KB, L, NMAX, NNUC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IC, ITRAN
      CHARACTER(LEN=8), DIMENSION(4) :: LABEL ! dimensioned according to the Sweden header format
      REAL(KIND=wp), DIMENSION(3) :: X
!
!*** End of declarations rewritten by SPAG
!
!---- Banner header is placed on the printer file
!
      IF(LSWORD)WRITE(IWRITE,1000)
!
!=======================================================================
!
!     1st set of data is on the nuclear configuration
!
!=======================================================================
!
      READ(ITAPE)(LABEL(I),I=1,4)
      READ(ITAPE)NNUC
!
      WRITE(NFTTAIL)(LABEL(I),I=1,4)
      WRITE(NFTTAIL)NNUC, CFILLERX, CFILLERX, CFILLERX, CFILLERX
!
      IF(LSWORD)WRITE(IWRITE,2000)(LABEL(I),I=1,4)
      IF(LSWORD)WRITE(IWRITE,2010)NNUC
!
!...... Loop over nuclei - obtain positions and charges
!
      IF(LSWORD)WRITE(IWRITE,2015)
!
      DO I=1, LEN(CNAME)
         CNAME(I:I)=' '
      END DO
!
      DO I=1, NNUC
         IF(LSWORD)THEN
            READ(ITAPE)CNAME(1:4), II, (X(J),J=1,3), CH
         ELSE
            READ(ITAPE)CNAME, II, (X(J),J=1,3), CH
         END IF
         WRITE(NFTTAIL)CNAME, II, (X(J),J=1,3), CH
         IF(LSWORD)WRITE(IWRITE,2020)CNAME, II, (X(J),J=1,3), CH
      END DO
!
!=======================================================================
!
!     2nd set of data is on Primitive and Contracted Gaussian functions
!
!=======================================================================
!
!---- IBBAS is the number of primitive functions; IABAS the contracted
!
      READ(ITAPE)IBBAS, IABAS
!
!     (this [JMC no longer] overallocates itran and ctran)
      ALLOCATE(ic(iabas),itran(MAXSYM*NHTT),ctran(MAXSYM*NHTT))
      WRITE(NFTTAIL)IBBAS, IABAS, CFILLERX, CFILLERX, CFILLERX, CFILLERX
      IF(LSWORD)WRITE(IWRITE,2500)IBBAS, IABAS
!
!---- Array IC has one entry for each contracted function. It gives the
!     number of primitives in that contraction.
!
      READ(ITAPE)(IC(I),I=1,IABAS)
      WRITE(NFTTAIL)(IC(I),I=1,IABAS)
!
!      WRITE(IWRITE,2550)
!      WRITE(IWRITE,2560) (I,IC(I),I=1,IABAS)
!
!---- For each contracted function we loop over the number of primitives
!     within it and obtain data on each.
!
!      WRITE(IWRITE,2600)
!
      DO I=1, IABAS
!      WRITE(IWRITE,2610) I,IC(I)
         DO J=1, IC(I)
            IF(LSWORD)THEN
               READ(ITAPE)CL(1:4), KB, CNAKO(1:4), XX, CX
            ELSE
               READ(ITAPE)CL, KB, CNAKO, XX, CX ! JMC bugfix, was CNAKO(1:4) previously
            END IF
            WRITE(NFTTAIL)CL, KB, CNAKO, XX, CX
!      WRITE(IWRITE,2620) J,CL,KB,CNAKO,XX,CX
         END DO
      END DO
!
!=======================================================================
!
!     3rd set of data is on symmetry information
!
!=======================================================================
!
      READ(ITAPE)(LABEL(I),I=1,4)
      READ(ITAPE)IABAS
!
      WRITE(NFTTAIL)(LABEL(I),I=1,4)
      WRITE(NFTTAIL)IABAS, CFILLERX, CFILLERX, CFILLERX, CFILLERX
!
!      WRITE(IWRITE,3000) (LABEL(I),I=1,4)
!      WRITE(IWRITE,3010) IABAS
!
      DO I=1, IABAS
         READ(ITAPE)J, (ITRAN(K),CTRAN(K),K=1,J)
         WRITE(NFTTAIL)J, (ITRAN(K),CTRAN(K),K=1,J), CFILLERX, &
     &                 CFILLERX, CFILLERX, CFILLERX
!      WRITE(IWRITE,3020) I,J,(ITRAN(K),CTRAN(K),K=1,J)
      END DO
      DEALLOCATE(CTRAN)
!
!=======================================================================
!
!     4th set of data is on MULLIKEN Population information
!
!     This contains information on the nuclear charges and on the
!     angular behaviour of the symmetry adaped basis functions.
!
!=======================================================================
!
      READ(ITAPE)(LABEL(I),I=1,4)
      READ(ITAPE)NMAX, IABAS
      ALLOCATE(CNUCNAM(NMAX),CTRAN(NMAX))
!
      WRITE(NFTTAIL)(LABEL(I),I=1,4)
      WRITE(NFTTAIL)NMAX, IABAS, CFILLERX, CFILLERX, CFILLERX, CFILLERX
!
!      WRITE(IWRITE,4000) (LABEL(I),I=1,4)
!      WRITE(IWRITE,4010) NMAX,IABAS
!
      READ(ITAPE)NMAX, (CNUCNAM(I),CTRAN(I),I=1,NMAX)
      WRITE(NFTTAIL)NMAX, (CNUCNAM(I),CTRAN(I),I=1,NMAX), CFILLERX,& 
     &              CFILLERX, CFILLERX, CFILLERX
!
!      WRITE(IWRITE,4020) NMAX
!      WRITE(IWRITE,4030) (I,CNUCNAM(I),CTRAN(I),I=1,NMAX)
!
!      WRITE(IWRITE,4040) IABAS
!
      DO I=1, IABAS
         READ(ITAPE)J, K, L
         WRITE(NFTTAIL)J, K, L, CFILLERX, CFILLERX, CFILLERX, CFILLERX
!      WRITE(IWRITE,4050) I,J,K,L
      END DO
!
      DEALLOCATE(ic,itran,ctran,cnucnam)
!
!---- Write a banner footer to the output
!
      IF(LSWORD)WRITE(IWRITE,8000)
!
      RETURN
!
!---- Format Statements
!
 203  FORMAT(1X,2I4,3X,A4,I2,3X,A3,2F15.8)
 204  FORMAT(/1X,2I4,(4(I3,F7.3)))
!
 500  FORMAT(5X,'Header Card (1:65) : ',A,/)
 510  FORMAT(5X,'No. of symmetries in basis set = ',I3,//,5X,&
     &       'No. of basis functions per symmetry: ',/,(5X,I2,1X,I3))
 520  FORMAT(/,5X,'Nuclear Potential Energy = ',F15.7,' (Hartrees) ',&
     &       //,5X,'Cartesian/Spherical flag = ',I8,/)
!
 1000 FORMAT(///,10X,'Atomic Integrals File: Header Records ',/)
 1010 FORMAT(//,10X,'These are copied from unit = ',I3,' to ',I3,/)
!
 2000 FORMAT(5X,'Section 1 Header = ',4A,/)
 2010 FORMAT(5X,'Number of nuclear centers = ',I3,/)
 2015 FORMAT(5X,' Symbol ',1X,'No.',' (X,Y,Z) Co-ordinates ',T50,&
     &       ' Charge ',/)
 2020 FORMAT(5X,A,1X,I2,1X,3(F10.6,1X),F6.3)
 2500 FORMAT(//,5X,'Number of primitive  Gaussian functions = ',I4,/,&
     &       5X,'Number "" contracted  " " "    " " " "  = ',I4,//)
 2550 FORMAT(5X,'Primitive functions per contracted function:',/)
 2560 FORMAT((5X,8(I3,'.',1X,I3,1X)))
 2600 FORMAT(/,5X,'A N A L Y S I S  of  each  C O N T R A C T I O N',/)
 2610 FORMAT(/,5X,'Contracted function no. = ',I3,/,5X,&
     &       'No. of primitives       = ',I3,/)
 2620 FORMAT(1X,'Primitive = ',I2,1X,'CL = ',A,1X,'KB = ',I3,1X,&
     &       'CNAKO = ',A,1X,'XX = ',D12.5,1X,'CX = ',D12.5,1X)
!
 3000 FORMAT(/,5X,'Section 2 Header = ',4A,/)
 3010 FORMAT(5X,'Number of symmetry adapted functions = ',I3,/)
 3020 FORMAT(/,5X,'Symmetry function no. = ',I3,/,5X,&
     &       'No. of components     = ',I3,/,5X,&
     &       ' (ITRAN,CTRAN) pairs: ',//,(5X,4(I4,1X,D12.5,1X)))
!
 4000 FORMAT(/,5X,'Section 3 Header = ',4A,/)
 4010 FORMAT(5X,'Number of atomic nuclei (NMAX) = ',I3,//,5X,&
     &       'Number ""   ""   basis functions (IABAS) = ',I5,//)
 4020 FORMAT(5X,'NMAX = ',I4,&
     &       ' Nuclei/Charge  pair definitions follow: ',/)
 4030 FORMAT((5X,3(I3,'.',1X,A,1X,D12.5,1X)))
 4040 FORMAT(/,5X,'No. of items (IABAS) = ',I5,/)
 4050 FORMAT(5X,'Item = ',I3,' J = ',I3,' K = ',I3,' L = ',I3)
!
 8000 FORMAT(//,10X,'***** Formatting of header records completed ',/)
!
      END SUBROUTINE COPYD

      SUBROUTINE print_ukrmol_header (unit)
      USE version_control, ONLY: print_git_revision
      INTEGER, INTENT(IN) :: unit
      WRITE(unit, '("+-------------------------------------------------------------------------+")')
      WRITE(unit, '("|              _              |                                           |")')
      WRITE(unit, '("| || || ||//  | \          |  | University College London (C) 1994 - 2020 |")')
      WRITE(unit, '("| || || ||    | / |/\/\ /\ |  | Open University           (C) 2007 - 2020 |")')
      WRITE(unit, '("| \\_// ||\\  | \ | | | \/ |  |                                           |")')
      WRITE(unit, '("|                             |                                           |")')
      CALL print_git_revision(unit)
      END SUBROUTINE print_ukrmol_header

end module global_utils
