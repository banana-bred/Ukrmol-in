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
!
! Data and utility module for swmol3
! JMC June 2010 
!
      module swmol3_data

      use precisn, only : wp ! for specifying the kind of reals
      use params, only : lbuf, maxsym
      implicit none
      public

      save
 
! Parameters
 
      INTEGER, PARAMETER :: LUT=6 ! unit for printing

      INTEGER, PARAMETER :: NHT=7 ! maximum value of (l + 1), where l is angular momentum

      INTEGER, PARAMETER :: MAXPRIM=25 ! 1. maximum number of basis function sub-blocks per symmetry (see input keyword jco)?
                                       ! 2. when squared, the limit on the product of any two numbers of primitives per 
                                       ! sub-block (the na's)

      INTEGER, PARAMETER :: MAXC=800 ! max of: ((na*nf, summed for a centre)*nont, summed from 1 to natype).
                                     ! See the docu for the meaning of the variables...

      INTEGER, PARAMETER :: MAXNRC=15 ! maximum number of contractions (basis functions per primitive?)

      INTEGER, PARAMETER :: IG=40000 ! where has this value come from?  Only used in ONEL, ONELH and EXCOEF

      INTEGER, PARAMETER :: NCMAX=100 ! maximum number of distinct symmetrically non-redundant nuclei

      INTEGER, PARAMETER :: NWMAX=1000 ! in readin, the value of norb is tested against this parameter; the program 
                                       ! will stop if norb > nwmax.  norb is then used to allocate some arrays.  No 
                                       ! idea where the value of 1000 has come from; maybe just "large".  nwmax is 
                                       ! also used below to dimension arrays - maybe make them dynamic too????

      INTEGER, PARAMETER :: NHTT=nht*(nht+1)/2

      INTEGER, PARAMETER :: NH4=4*NHT-3

      INTEGER, PARAMETER :: NHQ=nht*nht*nht*nht

      INTEGER, PARAMETER :: KWD=(NHT+1)*(NHT+2)*(NHT+3)/6

      INTEGER, PARAMETER :: MXP2=MAXPRIM*MAXPRIM

      INTEGER, PARAMETER :: MAXKM=ncmax*nht*maxprim

      INTEGER, PARAMETER :: MAXPOPC=((MAXSYM-1)*MAXSYM + MAXSYM-1)*MAXSYM + MAXSYM ! used to dimension the PC population 
                            ! count factor array.  See the value of k set above the call to LMUL in DRIVE for where this
                            ! expression has come from.

      REAL(KIND=wp), PARAMETER :: bohr_radius=0.529177_wp ! The Bohr radius in units of Angstrom

      ! Note: now taking LBUF from the params module in global_data.f90, 
      ! for consistency between different programs

! Data
 
      ! The Cartesian -> tesseral harmonics transformation matrices are stored 
      ! in this linear array. The values are hardwired [via DATA statements in the 
      ! subroutine READIN in swmol3.f, - JMC this is no longer true] for l values up to and including 6. 
      ! JMC the dimension, 1596, follows from nht=7
      REAL(KIND=wp), DIMENSION(1596) :: CHARM

      ! An array that goes with charm, where the ith component of iadq (starting from i=1) 
      ! is (index-1) for orbital angular momentum l=i-1, with index specifying the component 
      ! of charm (starting from 1) where the section of the array for this l starts.
      INTEGER, DIMENSION(nht+1) :: IADQ=(/0, 1, 10, 46, 146, 371, 812, 1596/)

      ! this should go into global_data.f90 when the angular momentum parameter goes too.
      CHARACTER(LEN=1), DIMENSION(nht), PARAMETER :: ISPD = (/'s', 'p', 'd', 'f', 'g', 'h', 'i'/)

      ! This array contains a 1 to indicate that the corresponding harmonic is a contaminant;
      ! 0 otherwise.  Components 1:n(n+1)/2 in the nth column match the components of the second 
      ! column of array KKTYP in the params module for angular momentum value (n-1).
      ! JCNTAM is initialized to zero here, with its 1's set in the subroutine set_jcntam below.
      INTEGER, DIMENSION(nhtt,nht) :: JCNTAM=0

      ! this is analogous to KKTYP in global_data.f90
      CHARACTER(LEN=4), DIMENSION(84), PARAMETER :: KWO=(/'S   ', 'X   ', 'Y   ', 'Z   ', 'XX  ', &
     &     'XY  ', 'XZ  ', 'YY  ', 'YZ  ', 'ZZ  ', 'F300', 'F210', &
     &     'F201', 'F120', 'F111', 'F102', 'F030', 'F021', 'F012', &
     &     'F003', 'G400', 'G310', 'G301', 'G220', 'G211', 'G202', &
     &     'G130', 'G121', 'G112', 'G103', 'G040', 'G031', 'G022', &
     &     'G013', 'G004', 'H500', 'H410', 'H401', 'H320', 'H311', &
     &     'H302', 'H230', 'H221', 'H212', 'H203', 'H140', 'H131', &
     &     'H122', 'H113', 'H104', 'H050', 'H041', 'H032', 'H023', &
     &     'H014', 'H005', 'I600', 'I510', 'I501', 'I420', 'I411', &
     &     'I402', 'I330', 'I321', 'I312', 'I303', 'I240', 'I231', &
     &     'I222', 'I213', 'I204', 'I150', 'I141', 'I132', 'I123', &
     &     'I114', 'I105', 'I060', 'I051', 'I042', 'I033', 'I024', &
     &     'I015', 'I006'/)

      ! Initialized in RFINIT; used in EXPA
      REAL(KIND=wp), DIMENSION(nh4) :: RFA

      ! From CMPIN and its ENTRY points: local variables with the SAVE attribute and variables 
      ! from the ICOMP common block (which was not used anywhere else in the code)
      INTEGER, DIMENSION(mxp2*mxp2) :: IBE
      INTEGER :: NN, IBEM

      ! Contents of common blocks

      ! These were all defined in the main program
      REAL(KIND=wp), DIMENSION(MXP2) :: ALP, BET, CAB, CCD, CSS1, CSS2, DEL, GAM, S1, S2
      INTEGER, DIMENSION(maxsym,maxsym) :: AND, EOR, OR ! JMC replaced hardwired dimension 8 with maxsym
      INTEGER, DIMENSION(maxsym,nwmax) :: DSTRT ! JMC replaced hardwired dimension 8 with maxsym
      LOGICAL :: HARMT
      INTEGER, DIMENSION(3) :: ISYTYP
      INTEGER :: ITAG, KMAX, LU2, MAXLOP, MAXLOT, NMAX, NRSS, NTAP, NUCZ
      INTEGER, DIMENSION(nht,nhtt) :: ITYPE
      INTEGER, DIMENSION(maxkm) :: JRS, JSTRT, MST, MUL, NHKT, NRCO, NSTRT, NUCO
      INTEGER, DIMENSION(nht) :: KHKT
      INTEGER, DIMENSION(maxsym) :: MSTOLD, MULT, NPAR, NPARSU ! JMC replaced hardwired dimension 8 with maxsym
      INTEGER, DIMENSION(ncmax) :: MULNUC
      INTEGER, DIMENSION(nwmax) :: NEWIND
      REAL(KIND=wp), DIMENSION(MAXPOPC) :: PC ! JMC replaced hardwired dimension 512
      REAL(KIND=wp), DIMENSION(MXP2,maxsym,3) :: PV, QV ! JMC replaced hardwired dimension 8 with maxsym
!      COMMON /INDX  / PC, DSTRT, NTAP, LU2, NRSS, NUCZ, ITAG, MAXLOP,
!     &                MAXLOT, KMAX, NMAX, KHKT, MULT, ISYTYP, ITYPE,
!     &                AND, OR, EOR, NPARSU, NPAR, MULNUC, NHKT, MUL,
!     &                NUCO, NRCO, JSTRT, NSTRT, MST, JRS
!      COMMON /ONE1  / ALP, BET, S1, CSS1, PV, CAB, GAM, DEL, S2, CSS2,
!     &                QV, CCD
!      COMMON /PRT   / HARMT
!      COMMON /REP   / NEWIND, MSTOLD

      ! These were not defined in the main program:

      REAL(KIND=wp) :: TLA, TLC
      REAL(KIND=wp), DIMENSION(maxkm) :: ALPHA
      REAL(KIND=wp), DIMENSION(3,maxkm) :: CENT
      REAL(KIND=wp), DIMENSION(ncmax) :: CHARGE
      REAL(KIND=wp), DIMENSION(maxc) :: CONT
      REAL(KIND=wp), DIMENSION(ncmax,3) :: CORD
      REAL(KIND=wp), DIMENSION(maxsym) :: FMULT ! JMC replaced hardwired dimension 8 with maxsym
      INTEGER, DIMENSION(nwmax) :: LABORB
!      COMMON /DAT   / ALPHA, CONT, CENT, CORD, CHARGE, FMULT, TLA, TLC
!      COMMON /SYMIND/ LABORB

      REAL(KIND=wp), DIMENSION(3,maxsym,nh4) :: XAND ! JMC replaced hardwired dimension 8 with maxsym
!      COMMON /XA    / XAND

!      REAL(KIND=wp) :: DISTPQ, DISTRS, RAX, RAY, RAZ, RPRX, RPRY, RPRZ, TOLR, TP52
      REAL(KIND=wp) :: TP52
      REAL(KIND=wp), DIMENSION(mxp2,maxsym) :: SABB, SCDD ! JMC replaced hardwired dimension 8 with maxsym
!      COMMON /CON   / RAX, RAY, RAZ, RPRX, RPRY, RPRZ, DISTPQ, DISTRS,
!     &                TOLR, TP52, SABB, SCDD

      REAL(KIND=wp), DIMENSION(4*LBUF) :: BUF
      INTEGER, DIMENSION(2,4*LBUF) :: IBUF
      LOGICAL :: IFD3
      INTEGER :: KPQR, MULA, MULB, MULC, MULD, NSTA, NSTB, NSTC, NSTD
      INTEGER, DIMENSION(nhtt) :: ITYA, ITYB, ITYC, ITYD
      INTEGER, DIMENSION(maxsym) :: JA, JB, JC, JD ! JMC replaced hardwired dimension 8 with maxsym
      INTEGER, DIMENSION(4) :: N2TAPE, NBLOCK, NUMSM, NUTTY
      INTEGER, DIMENSION(maxsym*maxsym) :: LEOR8 ! JMC replaced hardwired dimension 64 with maxsym*maxsym
!      COMMON /ITY   / ITYA, ITYB, ITYC, ITYD, KPQRS, NSTA, NSTB, NSTC,
!     &                NSTD, MUT, NUMSM, NBLOCK, JA, JB, JC, JD, BUF,
!     &                IBUF, NUTTY, N2TAPE, LEOR8, IFD3, MULA, MULB,
!     &                MULC, MULD

      INTEGER :: M2, M3, MIM, NMAS, MX, MRAF 
!      COMMON /MMMM  / M2, M3, MUM, MIM, NMAS, MX, MRAF

      INTEGER :: NNN
!      COMMON /TAB   / NNN

      INTEGER, DIMENSION(nhq) :: INDXY, INDZ
!      COMMON /DOR   / INDXY, INDZ

      LOGICAL :: IFD1X, IFD2X
      INTEGER :: IVA, JSTA, JSTQ, JSTR, JSTS, KQ, KR, KS, NAQRS, NCQRS, NRCA, &
                 NRCQ, NRCR, NRCS, NUCA, NUCAQ, NUCQ, NUCR, NUCRS, NUCS
!      COMMON /CONI  / NUCR, NUCS, NUCA, NUCQ, NRCR, NRCS, NRCA, NRCQ,
!     &                JSTR, JSTS, JSTA, JSTQ, IFD1X, IFD2X, NUCAQ,
!     &                NUCRS, NAQRS, NCQRS, KQ, KR, KS, IVA

      ! Blank common
      REAL(KIND=wp), DIMENSION(nh4) :: FACT, FACTM, RFACT, RFACTM
      LOGICAL :: IFD1, IFD2
      LOGICAL, DIMENSION(3) :: IFPL
      INTEGER :: KBCD, KCD, KHKTA, KHKTB, KHKTC, KHKTD, MAA, NHBCD, &
                 NHCD, NHKTA, NHKTB, NHKTC, NHKTD, NN1, NN2, NNB, NNC
!      COMMON  FACT, RFACT, FACTM, RFACTM, MAA, IFD1, IFD2, KCD, KBCD,
!     &        NHCD, NHBCD, NNC, IFPL, NN1, NN2, NHKTA, NHKTB, NHKTC,
!     &        NHKTD, KHKTA, KHKTB, KHKTC, KHKTD, NNB

! Subprograms

      CONTAINS

      INTEGER FUNCTION INEQV4(I,J,K,L)
      IMPLICIT NONE
!       IEOR returns the exclusive OR, bit by bit.
!       Note: this is a little function to replace a statement function
!       (statement fns are obsolescent in Fortran95)

! Dummy arguments
      INTEGER, INTENT(IN) :: I, J, K, L

      INEQV4=ieor(ieor(I,J),ieor(K,L))

      END FUNCTION INEQV4

      INTEGER FUNCTION IPK3(I,J,K)
      IMPLICIT NONE
!       Note: this is a little function to replace a statement function
!       (statement fns are obsolescent in Fortran95)
! MAXLOP is available via host association from the data part of this module

! Dummy arguments
      INTEGER, INTENT(IN) :: I, J, K

      IPK3=(I*MAXLOP+J)*MAXLOP+K

      END FUNCTION IPK3

      SUBROUTINE SET_JCNTAM

      JCNTAM(6,3)     = 1
      JCNTAM(6,4)     = 1 
      JCNTAM(9:10,4)  = 1 
      JCNTAM(9:11,5)  = 1
      JCNTAM(13:15,5) = 1
      JCNTAM(6,6)     = 1
      JCNTAM(12:13,6) = 1
      JCNTAM(15:28,6) = 1
      JCNTAM(10,7)    = 1
      JCNTAM(13,7)    = 1
      JCNTAM(15:16,7) = 1
      JCNTAM(18:28,7) = 1

      END SUBROUTINE SET_JCNTAM

      SUBROUTINE SET_HARMONIC_TRANS
! This subroutine solely defines the Cartesian -> tesseral harmonics 
! transformation matrices.  It's done this way because DATA statements 
! can't be used with objects accesssed via use association (i.e stored in a 
! module) so we can't keep the original initialization in swmol3's subroutine
! READIN.  Also, the array is so large it can't be initialized as part of the 
! declaration...
      charm(1:146) = (/1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 1.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, -1.0_wp, &
     &     0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 2.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp/)
      charm(147:371) = (/3.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, 0.0_wp, -24.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, -24.0_wp, 0.0_wp, 8.0_wp, 0.0_wp, &
     &     -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 6.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -3.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, -6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, -6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -3.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 2.0_wp/)
      charm(372:539) = (/1.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, -12.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, -12.0_wp, 0.0_wp, 8.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 2.0_wp, 0.0_wp, -12.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, -12.0_wp, 0.0_wp, 8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
     &     0.0_wp, 2.0_wp, 0.0_wp, 8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, &
     &     -24.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -10.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -6.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp/)
      charm(540:707) = (/0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -2.0_wp, 0.0_wp, 24.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, -8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 15.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 30.0_wp, 0.0_wp, -40.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 15.0_wp, 0.0_wp, -40.0_wp, 0.0_wp, 8.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     -10.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, &
     &     0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -1.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -1.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, &
     &     2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 2.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp/)
      charm(708:812) = (/0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     2.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 2.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     2.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, &
     &     0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -6.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 2.0_wp/)
      charm(813:980) = (/1.0_wp, 0.0_wp, 0.0_wp, -15.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 15.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, &
     &     -16.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     -16.0_wp, 0.0_wp, 16.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -10.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -1.0_wp, 0.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, 10.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     5.0_wp, 0.0_wp, -60.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -1.0_wp, 0.0_wp, 10.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -10.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, -16.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 16.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
     &     16.0_wp, 0.0_wp, -16.0_wp, 0.0_wp, 0.0_wp/)
      charm(981:1148) = (/0.0_wp, 6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, &
     &     0.0_wp, 8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 9.0_wp, 0.0_wp, &
     &     -24.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 40.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, -40.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -5.0_wp, 0.0_wp, 0.0_wp, &
     &     -15.0_wp, 0.0_wp, 90.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -15.0_wp, 0.0_wp, &
     &     180.0_wp, 0.0_wp, -120.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -5.0_wp, 0.0_wp, 90.0_wp, 0.0_wp, -120.0_wp, 0.0_wp, 16.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, -9.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -6.0_wp, &
     &     0.0_wp, 24.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     3.0_wp, 0.0_wp, -8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)
      charm(1149:1316) = (/1.0_wp, 0.0_wp, 0.0_wp, -5.0_wp, 0.0_wp, 1.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -5.0_wp, 0.0_wp, -6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 20.0_wp, 0.0_wp, -40.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 0.0_wp, -40.0_wp, 0.0_wp, &
     &     16.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, -5.0_wp, 0.0_wp, -6.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 5.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, &
     &     6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 10.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 20.0_wp, 0.0_wp, -40.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 0.0_wp, -40.0_wp, &
     &     0.0_wp, 16.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)
      charm(1317:1484) = (/0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, -6.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     -3.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, &
     &     2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 2.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
     &     0.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 3.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
     &     -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 2.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, &
     &     -1.0_wp, 0.0_wp, 0.0_wp/)
      charm(1485:1596) = (/0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -6.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     4.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 9.0_wp, 0.0_wp, -21.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 9.0_wp, 0.0_wp, -42.0_wp, 0.0_wp, -16.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, -21.0_wp, 0.0_wp, -16.0_wp, &
     &     0.0_wp, 8.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
     &     1.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     3.0_wp, 0.0_wp, 6.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
     &     0.0_wp, 1.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 1.0_wp/)

      END SUBROUTINE SET_HARMONIC_TRANS

      end module swmol3_data
