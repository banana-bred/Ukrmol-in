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
module cdenprop_aux

    use cdenprop_defs, only: idp

    implicit none

    public change_refdet, mkorbs, makemg, modrda, tmtma, tmtma_sparse, indexx, dryrun, setndj

    private zero

contains

    !> \brief   Transform packed determinants to a different reference determinant
    !> \authors A Harvey, J Benda
    !> \date    2011 - 2019
    !>
    !> Transform the configurations so that they correspond to the given reference determinant. This is advantageous
    !> for later fast comparison of determinants.
    !>
    subroutine change_refdet (csf_head, csf_body, ndtrf, no_spin_orbitals)

        use cdenprop_defs,  only: CSFheader, CSFbody
        use mpi_gbl,        only: local_rank, master, shared_communicator, mpi_xermsg, mpi_mod_bcast, shared_enabled
        use mpi_memory_gbl, only: mpi_memory_allocate_integer, mpi_memory_deallocate_integer, mpi_memory_synchronize

        type(CSFheader),      intent(inout) :: csf_head
        type(CSFbody),        intent(inout) :: csf_body
        integer, allocatable, intent(in)    :: ndtrf(:)
        integer,              intent(in)    :: no_spin_orbitals

        integer, allocatable :: imdi(:)
        integer              :: imdisz, isize_mdi, lndofj

        ! Only the node master will allocate and compute the new determinant offsets
        if (local_rank == master .or. .not. shared_enabled) then
            imdisz = 5 * (csf_head%lndof + max(csf_head%lndof / 2, 4))
            allocate (imdi(imdisz + 1000))

            call modrda(csf_head % nelt, no_spin_orbitals, ndtrf, csf_head % ndtrf, imdi, sum(csf_head % nodo), &
                        isize_mdi, csf_body % ndo, csf_body % cdo)

            if (isize_mdi > imdisz) then
                call mpi_xermsg('cdenprop_aux', 'change_refdet', 'IMDI too small.', isize_mdi - imdisz, 1)
            end if

            call setndj(csf_head % nocsf, imdi, csf_body % indo, csf_head % nodo, lndofj)
            csf_body % indo(csf_head%nocsf + 1) = lndofj
        end if

        ! Communicate the new determinant store size to other processes on the same node
        call mpi_mod_bcast(lndofj, master, shared_communicator)

        ! Save transformed determinants and associated data back in csf_head/body_j
        csf_head % ndtrf = ndtrf
        csf_head % lndof = lndofj
        if (size(csf_body % ndo) < lndofj) then
            call mpi_memory_deallocate_integer(csf_body % ndo, size(csf_body % ndo), csf_body % ndo_window, shared_communicator)
            csf_body % ndo_window = mpi_memory_allocate_integer(csf_body % ndo, lndofj, shared_communicator)
        end if
        if (local_rank == master .or. .not. shared_enabled) then
            csf_body % ndo(1:lndofj) = imdi(1:lndofj)
        end if

        ! Let other processes know that the shared memory was changed
        call mpi_memory_synchronize(csf_body % ndo_window,  shared_communicator)
        call mpi_memory_synchronize(csf_body % cdo_window,  shared_communicator)
        call mpi_memory_synchronize(csf_body % indo_window, shared_communicator)

    end subroutine change_refdet


!! OLD DENPROP ROUTINES

      SUBROUTINE MKORBS(ISYMTYP,NOB,NSYM,MN,MG,MM,MS,maxk,NPFLG)
! C***********************************************************************
! C
! C     MKORBS - Computes the orbital table which is then used in the
! C              evaluation of off diagonal density matrix elements
! C
! C     Input data:
! C        ISYMTYP  Switch for C-inf-v (=0 or 1) / Abelian point group (=2
! C            NOB  Number of orbitals per symmetry
! C           NSYM  Number of symmetries in the orbital set
! C          NPFLG  Flag controlling printing of computed orbital table
! C
! C     Output data:
! C              MN  Orbital number associated with each spin-orbital
! C              MG  G/U designation for each spin-orbital (C-inf-v only)
! C                  Actually this is always zero because C-inf-v does not
! C                  distinguish between g/u. It exists because original
! C                  version of Alchemy tried to use it for D-inf-h too;
! C                  all CI evauation is doen in C-inf-v now because CONGE
! C                  converts D-inf-h to C-inf-v data.
! C              MM  Symmetry quantum number associated with each spin-orb
! C              MS  Spin function ( alpha or beta ) associated with each
! C                  spin orbital
! C
! C     Notes:
! C
! C       The orbital table establishes orbital and quantum number data fo
! C     each spin orbital in the set.
! C
! C     e.g. C-inf-v symmetry with NSYM=2, NOB=3,1, yields ten spin
! C          orbitals which are designated as follows by this routine:
! C
! C        Spin orb.     MN  MG  MM  MS     Comments
! C            1          1   0   0   0     1 sigma spin up
! C            2          1   0   0   1     1 sigma spin down
! C            3          2   0   0   0     2 sigma spin up
! C            4          2   0   0   1     2 sigma spin down
! C            5          3   0   0   0     3 sigma spin up
! C            6          3   0   0   1     3 sigma spin down
! C            7          4   0   1   0     1 pi(lambda=+1) spin up
! C            8          4   0   1   1     1 pi(lambda=+1) spin down
! C            9          4   0  -1   0     1 pi(lambda=-1) spin up
! C           10          4   0  -1   1     1 pi(lambda=-1) spin down
! C
! C***********************************************************************
      !USE precisn, ONLY : wp ! for specifying the kind of reals
      IMPLICIT NONE
! C
! C*** Start of declarations rewritten by SPAG
! C
! C Dummy arguments
! C
      INTEGER :: ISYMTYP, MAXK, NPFLG, NSYM
      INTEGER, DIMENSION(*) :: MG, MM, MN, MS, NOB
      INTENT (IN) ISYMTYP, NOB, NPFLG, NSYM
      INTENT (OUT) MAXK
      INTENT (INOUT) MG, MM, MN, MS
! C
! C Local variables
! C
      INTEGER :: I, J, K, M, MA, MB, N, NSRB
! C
! C*** End of declarations rewritten by SPAG
! C
! C---- 'Afore we compute NSRB we set it to -1 so that we can error trap
! C     this routine. Should something happen such that the loops are
! C     not executed then NSRB is never reset to a proper value and so
! C     the test fails at the end of the routine.
! C
      NSRB=-1
! C
      IF(ISYMTYP.LT.2)THEN
! C
! C=======================================================================
! C
! C     Linear molecules i.e. C-inf-v point group
! C
! C=======================================================================
! C
! C---- Sigma type orbital are done first.
! C
         I=1
         MA=NOB(1)
         DO N=1, MA
            MN(I)=N
            MG(I)=0
            MM(I)=0
            MS(I)=0
            I=I+1
            MN(I)=N
            MG(I)=0
            MM(I)=0
            MS(I)=1
            I=I+1
         END DO
! C
! C---- Now do Pi, Delta, Phi orbitals etc.....
! C
         K=MA+1
         DO M=2, NSYM
            MA=NOB(M)
            MB=M-1
            DO N=1, MA
               DO J=1, 4
                  MN(I)=K
                  MG(I)=0
                  MM(I)=MB
                  MS(I)=0
                  I=I+1
               END DO
               K=K+1
               MM(I-1)=-MB
               MM(I-2)=-MB
               MS(I-1)=1
               MS(I-3)=1
            END DO
         END DO
! C
      ELSE
! C
! C=======================================================================
! C
! C     Abelian point group molecules (D2h and subgroups)
! C
! C=======================================================================
! C
! C---- Loop over all symmetries
! C          Loop over all orbitals within the symmetry
! C               Assign numbers for both spin orbitals associated
! C               with each orbital.
! C
          I=1
         K=1
         DO M=1, NSYM
            MA=NOB(M)
            MB=M-1
            DO N=1, MA
               MN(I)=K
               MG(I)=0
               MM(I)=MB
               MS(I)=0
               I=I+1
               MN(I)=K
               MG(I)=0
               MM(I)=MB
               MS(I)=1
               I=I+1
               K=K+1
            END DO
         END DO
! C
      END IF
! C
      NSRB=I-1
      maxk=k-1
! C
! C---- Error check for a stupid value of NSRB.
! C
! !       IF(NSRB.LE.0)THEN
! !          WRITE(6,9900)
! !          WRITE(6,9940)NSRB
! !          CALL TMTCLOS()
! !       END IF
! C
      IF(NPFLG.GE.1)THEN
         WRITE(6,505)
         WRITE(6,500)(I,MN(I),MG(I),MM(I),MS(I),I=1,NSRB)
         WRITE(6,510)NSRB
      END IF
! ! C
      RETURN
! C
! C---- Format Statements
! C
 500  FORMAT(8X,'I',4X,'N',4X,'G',4X,'M',4X,'S',/,(4X,5I5))
 505  FORMAT(//,5X,' THE ORBITAL TABLE FOR THIS CASE IS :',/)
 510  FORMAT(5X,' NO. OF SPIN ORBITALS FOR THIS RUN IS',I8,/)
! C
 9900 FORMAT(//,10X,'**** Error in MKORBS ',//)
 9940 FORMAT(10X,'NSRB is invalid = ',I10,/)
! C
      END SUBROUTINE MKORBS

      SUBROUTINE MODRDA(NELT,NSRB,NDTRI,NDTRJ,MDI,NODA,MDA,NDI,CA)
!       USE precisn, ONLY : wp ! for specifying the kind of reals
      use cdenprop_defs
      IMPLICIT NONE
! C
! C*** Start of declarations rewritten by SPAG
! C
! C Dummy arguments
! C
      INTEGER :: MDA, NELT, NODA, NSRB
      REAL(KIND=idp), DIMENSION(noda) :: CA
      INTEGER, DIMENSION(*) :: MDI, NDI
      INTEGER, DIMENSION(nelt) :: NDTRI, NDTRJ
      INTENT (IN) NDI, NDTRI, NDTRJ, NELT, NODA, NSRB
      INTENT (OUT) MDI
      INTENT (INOUT) CA, MDA
! C
! C Local variables
! C
      INTEGER :: I, J, M, MDB, MDC, N, ND, NDA, NDB, NDC, NP, NPS, NQ,&
     &           NR, NS, NT
      INTEGER, DIMENSION(nsrb) :: MA, MB, NA, NB, NTT
! C
! C*** End of declarations rewritten by SPAG
! C
! C***********************************************************************
! C
! C     MODRDA - MODIFY REFERENCE DETERMINANTS FOR THE DETERMINANTAL
! C              EXCITATION CODES.
! C
! C     IN CALCULATING THE OFF-DIAGONAL ELEMENTS OF THE DENSITY MATRIX
! C     IT IS SOMETIMES NECESSARY TO COMPARE DETERMINANTS WHICH HAVE BEEN
! C     CALCULATED WITH RESPECT TO DIFFERENT REFERENCE DETERMINANTS.SINCE
! C     THE DENSITY MATRIX GENERATION ALGORITHM EMPLOYED IN THIS CODE
! C     ASSUMES THAT THE DETERMINANTAL EXCITATION CODES FOR THE TWO GIVEN
! C     WAVEFUNCTIONS HAVE BEEN GENERATED FROM THE SAME REFERENCE
! C     DETERMINANT IT IS NECESSARY TO MODIFY THE DETERMINANTAL EXCITATION
! C     CODES OF ONE OF THE WAVEFUNCTIONS TO BE WITH RESPECT TO THE
! C     REFERENCE DETERMINANT OF THE OTHER. THIS SUBROUTINE CARRIES OUT
! C     THIS TASK THROUGH THE ENTRY POINT MODRDB.
! C
! C     INPUT DATA :
! C           NELT  NUMBER OF ELECTRONS
! C           NSRB  NUMBER OF SPIN ORBITALS
! C          NDTRI  REFERENCE DETERMINANT I
! C          NDTRJ  REFERENCE DETERMINANT J
! C           NODA  NUMBER OF DETERMINANTS IN WAVEFUNCTION J
! C            NDI  THE EXCITATION CODES FOR THE DETERMINANTS OF
! C                 WAVEFUNCTION J RELATIVE TO REFERENCE DETERMINANT I
! C             CA  COEFFICIENTS OF THE EXCITATION CODES IN WAVEFUNCTION J
! C
! C    OUTPUT DATA :
! C             NP  NUMBER OF UN-MATCHED ORBITALS SPIN ORBITALS
! C            NPS  PHASE FACTOR DETERMINED BY NUMBER OF UNMATCHED
! C                 ORBITALS
! C             NA  SPIN ORBITALS IN REFERENCE DETERMINANT I THAT ARE NOT
! C                 FOUND IN REFERENCE DETERMINANT J
! C             NB  SPIN ORBITALS IN REFERENCE DETERMINANT J THAT ARE NOT
! C                 FOUND IN REFERENCE DETERMINANT I
! C            MDI  THE NEW EXCITATION CODES FOR THE DETERMINANTS OF
! C                 WAVEFUNCTION J RELATIVE TO REFERENCE DETERMINANT I
! C            MDA  SIZE OF THE ARRAY MDI
! C             CA  THE COEFFICIENTS FOR EACH OF THE NEW EXCITATION CODES
! C                 THESE ARE THE SAME AS THE INPUT VALUES BUT POSSIBLY
! C                 MULTIPLIED BY -1.
! C
! C    TEMPORARY STORAGE AREAS:
! C            NTT  STORES THE REFERENCE DETERMINANT J DURING COMPARISON
! C             MA
! C             MB
! C
! C***********************************************************************
! C
! C-----REF DET J IS COPIED TO TEMPORARY AREA NTT;VARIABLES ARE INITIALIZE
! C
      DO I=1, NELT
         NTT(I)=NDTRJ(I)
      END DO
      NP=0
      NPS=1
! C
! C-----REF DET I IS SCANNED AND COMPARED TO REF DET J. AN ATTEMPT IS MADE
! C     TO MATCH EACH SPIN ORBITAL IN REF DET I WITH ITSELF IN REF DET J.
! C     NO MATCH IS FOUND THEN INFORMATION IS  STORED IN ARRAYS
! C     NA (THE SO NUMBER) AND NB (ITS POSITION IN REF DET I). IF A MATCH
! C     FOUND THEN IF THE SO OCCURS AT THE SAME POSITION IN EACH REF DET N
! C     ACTION IS NECESSARY. ON THE OTHER HAND IF THE SO OCCURS AT DIFFERE
! C     POSITIONS THEN REF DET J IS RE-ORDERED TO HAVE THE SO IN THE SAME
! C     POSITION AS IN REF DET I. A CHANGE OF PHASE IS REQUIRED. THIS
! C     APPROACH ALSO FACILITATES THE GENERATION OF AN EXCITATION CODE FOR
! C     REF DET J W.R.T. REF DET I. IT IS ONLY NECESSARY TO FIND THE SO'S
! C     IN REF DET I WHICH ARE UN-MATCHED AND THEN LOOK AT THE SO'S IN THE
! C     CORRESPONDING POSITIONS IN REF DET J.
! C
      DO I=1, NELT
         M=NDTRI(I)
         DO J=1, NELT
            IF(M.EQ.NTT(J))THEN
               IF(I.NE.J)THEN
                  N=NTT(J)
                  NTT(J)=NTT(I)
                  NTT(I)=N
                  NPS=-NPS
               END IF
               GO TO 30
            END IF
         END DO
         NP=NP+1
         NA(NP)=M
         NB(NP)=I
 30   END DO
! C
! C-----CHANGE NB TO BE SO'S RATHER THAN POINTERS TO STORAGE LOCATIONS
! C
      DO I=1, NP
         M=NB(I)
         NB(I)=NTT(M)
      END DO
! C
! C-----ARRAY NTT IS RE-INITIALIZED FOR SUBSEQUENT USE
! C
      DO I=1, Nsrb
         NTT(I)=0
      END DO
! C
! C-----RE-ORGANIZATION OF THE DETERMINANTS OF WAVEFUNCTION J BEGINS HERE
! C
      MDA=1
      NDA=1
! C
! C-----IF THE REFERENCE DETERMINANTS FOR BOTH WAVEFUNCTIONS ARE EQUAL,THE
! C     JUST COPY DETERMINANTS FROM NDI TO MDI.NO MODIFICATION IS REQUIRED
! C
      IF(NP.EQ.0)THEN
         DO ND=1, NODA
            NDB=NDI(NDA)
            NDA=NDA+2*NDB+1
         END DO
         MDA=NDA-1
         DO I=1, MDA
            MDI(I)=NDI(I)
         END DO
         RETURN
      END IF
! C
! C-----MODIFICATION IS REQUIRED.LOOP OVER DETERMINANTS AND MODIFY EACH ON
! C
      DO ND=1, NODA
         NDB=NDI(NDA)
         NDC=NDA+NDB
         NS=NPS
! C
! C--------THERE ARE NO REPLACEMENTS FROM REFERENCE DETERMINANT J FOR THIS
! C        PARTICULAR DETERMINANT. EXCITATION CODE IS JUST THAT FOR REFDET
! C        W.R.T. REF DET I.
! C
         IF(NDB.EQ.0)THEN
            MDI(MDA)=NP
            MDB=NP
            MDC=MDA+NP
            DO I=1, NP
               MDI(MDA+I)=NA(I)
               MDI(MDC+I)=NB(I)
            END DO
            GO TO 230
         END IF
! C
! C--------THERE ARE REPLACEMENTS FROM REFERENCE DETERMINANT J
! C        PERFORM INITIALIZATION FOR EACH DETERMINANT.
! C
         DO I=1, NP
            MA(I)=NA(I)
            MB(I)=NB(I)
         END DO
! C
         NR=NP
         DO I=1, NDB
            MA(NR+I)=NDI(NDA+I)
            MB(NR+I)=NDI(NDC+I)
         END DO
         NR=NP+NDB
! C
         DO I=1, NR
            NTT(MA(I))=I
         END DO
! C
! C--------MODIFY THE EXCITATION CODES. SPIN ORBITALS TO BE REPLACED ARE
! C        IN ARRAY MA WHILE THOSE THAT REPLACE ARE IN MB
! C
         MDB=0
         DO I=1, NP
            M=NTT(MB(I))
            IF(M.NE.0)THEN
               NTT(MA(I))=M
               NTT(MB(I))=0
               N=MA(M)
               MA(M)=MA(I)
               MA(I)=N
            ELSE
               MDB=MDB+1
            END IF
         END DO
! C
         NQ=NP+1
         DO I=NQ, NR
            M=NTT(MB(I))
            IF(M.NE.0)THEN
               IF(M.NE.I)THEN
                  NTT(MA(I))=M
                  NS=-NS
                  N=MA(M)
                  MA(M)=MA(I)
                  MA(I)=N
               END IF
               NTT(MB(I))=0
            ELSE
               MDB=MDB+1
            END IF
         END DO
! C
! C--------STORE THE NEW EXCITATION CODES INTO ARRAY MDI
! C
         MDI(MDA)=MDB
         MDC=MDA+MDB
         IF(MDB.NE.0)THEN
            NT=0
            DO I=1, NR
               M=NTT(MA(I))
               IF(M.NE.0)THEN
                  NT=NT+1
                  MDI(MDA+NT)=MA(I)
                  MDI(MDC+NT)=MB(I)
                  NTT(MA(I))=0
               END IF
            END DO
         END IF
! C
! C--------MODIFY THE PHASE FACTOR OF THE DETERMINANT
! C
 230     IF(NS.LT.0)CA(ND)=-CA(ND)
! C
! C--------AUGUMENT STORAGE POSITIONS IN THE DETERMINANT ARRAYS
! C        FOR THE NEXT PASS
! C
         MDA=MDC+MDB+1
         NDA=NDC+NDB+1
      END DO
! C
      MDA=MDA-1
! c
      RETURN
! C
      END SUBROUTINE MODRDA

     SUBROUTINE SETNDJ(NOCSF,JDETS,JDETPN,NODETJ,LNDOFN)
!       USE precisn, ONLY : wp ! for specifying the kind of reals
      IMPLICIT NONE
! C
! C*** Start of declarations rewritten by SPAG
! C
! C Dummy arguments
! C
      INTEGER :: LNDOFN, NOCSF
      INTEGER, DIMENSION(*) :: JDETPN, JDETS, NODETJ
      INTENT (IN) JDETS, NOCSF, NODETJ
      INTENT (OUT) JDETPN, LNDOFN
! C
! C Local variables
! C
      INTEGER :: I, IST, J, NDETS, NRP
! C
! C*** End of declarations rewritten by SPAG
! C
! C***********************************************************************
! C
! C     SETNDJ - SET NEW DETERMINANTS FOR WAVEFUNCTION J
! C
! C              THIS ROUTINE TAKES THE NEW ARRAY OF DETERMINANT CODES FOR
! C              WAVEFUNCTION J W.R.T. REFERENCE DETERMINANT I AND SETS UP
! C              A POINTER ARRAY GIVING THE START OF THE FIRST DETERMINANT
! C              FOR EACH CSF. THIS IS REQUIRED FOR FORMULA GENERATION
! C              SUBROUTINES.
! C
! C     INPUT DATA :
! C          NOCSF  NUMBER OF CSFS IN WAVEFUNCTION J
! C          JDETS  NEW DETERMINANT CODES FOR WAVEFUNCTION J W.R.T. THE
! C                 REFERENCE DETERMINANT FOR WAVEFUNCTION I.
! C         NODETJ  THE NUMBER OF DETERMINANTS PER CSF IN WAVEFUNCTION J
! C
! C     OUTPUT DATA :
! C         JDETPN   THE POSITION OF THE FIRST DETERMINANT FOR EACH CSF
! C                  IN THE ARRAY JDETS.
! C         LNDOFN   THE SIZE OF THE ARRAY HOLDING THE NEW DETERMINANTS
! C                  FOR WAVEFUNCTION J.
! C
! C***********************************************************************
! C
      IST=1
      DO J=1, NOCSF
         NDETS=NODETJ(J)
         JDETPN(J)=IST
         DO I=1, NDETS
            NRP=JDETS(IST)
            IST=IST+(2*NRP)+1
         END DO
      END DO
! C
      LNDOFN=IST-1
! C
      RETURN
! C
      END SUBROUTINE SETNDJ

      SUBROUTINE ZERO(NA,NB,NZ,NBL,MAXSO)
!       USE precisn, ONLY : wp ! for specifying the kind of reals
      IMPLICIT NONE
! C
! C*** Start of declarations rewritten by SPAG
! C
! C Dummy arguments
! C
      INTEGER :: MAXSO, NZ
      INTEGER, DIMENSION(*) :: NA, NB, NBL
      INTENT (IN) MAXSO, NA, NB, NBL
      INTENT (INOUT) NZ
! C
! C Local variables
! C
      INTEGER :: I, K, M, MA, MAXM, MB
      INTEGER, DIMENSION(maxso) :: NBB
! C
! C*** End of declarations rewritten by SPAG
! C
! C***********************************************************************
! C
! C     ZERO -   COMPARES THE FIRST PAIR OF DETERMINANTS FOR EACH CSF PAIR
! C              AND RETURNS IN NZ THE NUMBER OF ORBITAL DIFFERENCES
! C              BETWEEN THEM. THIS IS SIMILAR TO SUBROUTINE TO PZERO IN
! C              PROGRAM SPEEDY.
! C
! C***********************************************************************
! c
! C
      maxm=1 ! Was not initialised before use AH
      DO I=1, MAXSO
         NBB(I)=0
      END DO
! C
! C-----FOR THE FIRST DETERMINANT OF THE FIRST WAVEFUNCTION, SCAN ALL
! C     REPLACED ORBITALS AND DECREMENT POPULATION BY ONE. THEN
! C     SCAN ALL OF THE REPLACEMENTS AND AUGUMENT POPULATION BY ONE.
! C
      MA=NA(1)
      MB=NB(1)
! c
      DO I=1, MA
         M=NBL(NA(I+1))
         NBB(M)=NBB(M)-1
         M=NBL(NA(MA+I+1))
         NBB(M)=NBB(M)+1
      END DO
! C
! C-----SCAN FIRST DETERMINANT OF THE SECOND WAVEFUNCTION AND APPLY INVERS
! C     PROCESS TO THAT FOR FIRST WAVEFUNCTION.
! C
      DO I=1, MB
         M=NBL(NB(I+1))
         NBB(M)=NBB(M)+1
         M=NBL(NB(MB+I+1))
         NBB(M)=NBB(M)-1
      END DO
! C
! C-----COUNT ALL OF THE AFFECTED ORBITALS AND FIND OVERALL CHANGE.
! C
      NZ=0
      DO I=1, MA+MA
         M=NBL(NA(I+1))
!          maxm=max(maxm,m)
         K=NBB(M)
         NZ=NZ+ABS(K)
         NBB(M)=0
      END DO
! C
      DO I=1, MB+MB
         M=NBL(NB(I+1))
         K=NBB(M)
         NZ=NZ+ABS(K)
         NBB(M)=0
      END DO
! C
      RETURN
! C
      END SUBROUTINE ZERO

      SUBROUTINE DRYRUN(INDI,JNDI,ICSF,JCSF,EVALUE,MN,IDETS,JDETS,MAXSO)
!       USE precisn, ONLY : wp ! for specifying the kind of reals
      IMPLICIT NONE
! C
! C*** Start of declarations rewritten by SPAG
! C
! C Dummy arguments
! C
      LOGICAL :: EVALUE
      INTEGER :: ICSF, JCSF, MAXSO
      INTEGER, DIMENSION(*) :: IDETS, INDI, JDETS, JNDI, MN
      INTENT (IN) ICSF, INDI, JCSF, JNDI
      INTENT (OUT) EVALUE
! C
! C Local variables
! C
      INTEGER :: MD, ND, NZ
! C
! C*** End of declarations rewritten by SPAG
! C
! C***********************************************************************
! C
! C     DRYRUN - PERFORMS A 'DRYRUN' EVALUATION
! C
! C              EACH CSF CONTAINS THE SAME ORBITAL CONFIGURATION, ONLY
! C              THE ALLOCATION TO SPIN ORBITALS IS DIFFERENT. THUS EACH
! C              CSF PAIR REFERS TO A FIXED COMPARISON OF ORBITALS. THIS
! C              SUBROUTINE LOOKS AT THAT ORBITAL COMPARISON AND ALLOCATES
! C              THE VALUE .TRUE. TO VARIABLE EVALUE IF THE NUMBER OF
! C              ORBITAL DIFFERENCES IS TWO, OTHERWISE EVALUE=.FALSE.
! C              FROM THE BASIC PROPERTIES OF THE FIRST ORDER SPIN
! C              REDUCED DENSITY MATRIX, IT IS CLEAR WHY THIS CHOICE IS
! C              MADE. THIS SUBROUTINE CLOSELY PARALLELS SUBROUTINE DRYRUN
! C              IN PROGRAM SPEEDY.
! C
! C    INPUT DATA :
! C          INDI  STARTING POSITIONS FOR DETERMINANTS OF EACH CSF IN
! C                ARRAY IDETS FOR WAVEFUNCTION 1.
! C          JNDI  STARTING POSITIONS FOR DETERMINANTS OF EACH CSF IN
! C                ARRAY JDETS FOR WAVEFUNCTION 2.
! C          ICSF  NUMBER OF CSF FOR WAVEFUNCTION 1.
! C          JCSF  NUMBER OF CSF FOR WAVEFUNCTION 2.
! C            MN  THE ORBITAL TABLE RELATING SPIN ORBITALS TO ORBITALS
! C         IDETS  THE DETERMINANTS FOR WAVEFUNCTION 1
! C         JDETS  THE DETERMINANTS FOR WAVEFUNCTION 2
! C         MAXSO  PARAMETER DETERMINING MAXIMUM NUMBER OF SPIN ORBITALS
! C           NBB  WORK SPACE ARRAY
! C
! C    OUTPUT DATA :
! C         EVALUE  LOGICAL FLAG WHICH IS SET DEPENDING ON WHETHER OR NOT
! C                 THIS PAIR OF CSFS ACTUALLY CONTRIBUTES.
! C
! C***********************************************************************
! C
      MD=INDI(ICSF)
      ND=JNDI(JCSF)
      CALL ZERO(IDETS(MD),JDETS(ND),NZ,MN,MAXSO)
      IF(NZ.GT.2)THEN
         EVALUE=.FALSE.

      ELSE
         EVALUE=.TRUE.

      END IF
! C
      RETURN
! C
      END SUBROUTINE DRYRUN

      FUNCTION IHJSR(NSIZE,NSEQ,IST,ITARG)
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: IST, ITARG, NSIZE
      INTEGER :: IHJSR
      INTEGER, DIMENSION(NSIZE+1) :: NSEQ
      INTENT (IN) IST, ITARG, NSIZE
      INTENT (INOUT) NSEQ
!
! Local variables
!
      INTEGER :: ICURSE
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     IHJSR - INTEGER SEARCH FUNCTION
!
!             THIS FUNCTION IMPLEMENTS THE ALGORITHM GIVEN IN THE BOOK
!             'LEARNING TO PROGRAM' BY HOWARD JOHNSTON ON PAGE 389 FOR
!             SEARCHING AN UNORDERED LIST FOR THE OCCURRENCE OF A GIVEN
!             VALUE. THE METHOD OF SETTING A SENTINEL IS USED.
!
!     INPUT DATA :
!           NSIZE  NUMBER OF ELEMENTS IN THE SEQUENCE
!            NSEQ  ARRAY HOLDING THE SEQUENCE WITH ONE EXTRA POSITION AT
!                  THE END FOR HOLDING THE SENTINEL.
!             IST  INITIAL POSITION IN THE SEQUENCE TO START SEARCHING
!                  FROM. THIS IS USUALLY 1.
!           ITARG  INTEGER WHICH IS BEING SOUGHT. IE. THE TARGET.
!
!     OUTPUT DATA :
!           IHJSR  IS THE POSITION IN THE SEQUENCE OF THE REQUIRED
!                  TARGET. IF THE TARGET IS NOT IN THE SEQUENCE, THEN
!                  IHJSR IS GIVEN A VALUE ONE LARGER THAN THE ACTUAL
!                  SIZE OF THE SEQUENCE.
!
!***********************************************************************
!
      IHJSR=0
!
      NSEQ(NSIZE+1)=ITARG
      DO ICURSE=IST, NSIZE+1
         IF(NSEQ(ICURSE).NE.ITARG)CYCLE
         IHJSR=ICURSE
         GO TO 200

      END DO
!
 200  CONTINUE
!
      RETURN
!
      END FUNCTION IHJSR



      FUNCTION ISRCHE(N,NAR,INC,ITAR)
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: INC, ITAR, N
      INTEGER :: ISRCHE
      INTEGER, DIMENSION(n) :: NAR
      INTENT (IN) NAR
!
! Local variables
!
      INTEGER :: I
!       INTEGER :: IHJSR
      INTEGER, DIMENSION(n+1) :: IWORK
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     ISRCHE - PREPARE CALL TO INTEGER SEARCH FUNCTION
!
!              SETS UP A CALL TO THE FUNCTION IHJSR WHICH LOCATES THE
!              TARGET INTEGER IN THE SEQUENCE NAR.
!
!     INPUT DATA :
!              N  DIMENSION OF THE SEQUENCE TO BE SEARCHED
!            NAR  SEQUENCE OF INTEGERS TO BE SEARCHED. ITS ACTUAL
!                 DIMENSION IS ONE LARGER THAN THE SEQUENCE SIZE.
!            INC  STARTING POSITION IN SEQUENCE
!           ITAR  THE INTEGER BEING SOUGHT IN THE SEQUENCE IE. THE
!                 TARGET
!
!     OUTPUT DATA :
!          ISRCHE  THE POSITION OF THE TARGET IN THE SEQUENCE IS
!                  RETURNED VIA THE FUNCTION NAME. IF THE TARGET IS NOT
!                  FOUND THEN THE VALUE OF ISRCHE IS ONE LARGER THAN
!                  THE SIZE OF THE SEQUENCE.
!
!     LOCAL DATA :
!          IWORK  THIS ARRAY IS A WORKSPACE INTO WHICH IS COPIED THE
!                 ARRAY NAR WHICH IS BEING SEARCHED. IT IS IMPORTANT
!                 THAT THE ARRAY NAR IS NOT OVERWRITTEN AS THERE MAY
!                 BE ELEMENTS BEYOND THE LIMITS WHICH WE ARE SEARCHING
!
!***********************************************************************
!
      DO I=1, N
         IWORK(I)=NAR(I)

      END DO
!
      ISRCHE=IHJSR(N,IWORK,INC,ITAR)
!
      RETURN
!
      END FUNCTION ISRCHE



SUBROUTINE TMTMA(NODA,CA,NDA,NODB,CB,NDB,nsrb,NELT,NDTR,MN,MS,&
     &                 MDTR,&!NFTD,NELMT,NBLK,NEL,LENGTH,II,JJ,NPQ,CPQ,&
     &                 ICSF,JCSF,density_matrix_coefficients, density_matrix_orbital_pairs)
      USE cdenprop_defs, ONLY : idp ! for specifying the kind of reals
      IMPLICIT NONE
! C
! C*** Start of declarations rewritten by SPAG
! C
! C Dummy arguments
! C
      INTEGER :: ICSF, JCSF, LENGTH, NBLK, NEL, NELMT, NELT, NFTD, NODA, &
     &           NODB, NSRB
      REAL(KIND=idp), DIMENSION(*) :: CA, CB
      REAL(KIND=idp), DIMENSION(:,:) :: density_matrix_coefficients
      INTEGER, DIMENSION(:,:,:) :: density_matrix_orbital_pairs
      INTEGER, DIMENSION(*) ::  MDTR, MN, MS, NDA, NDB, NDTR
!       INTEGER, DIMENSION(2,*) :: NPQ
      INTENT (IN) CA, CB, MDTR, MN, MS, NDA, NDTR, NELT, NODA, NODB, &
     &            NSRB
! C
! C Local variables
! C
      REAL(KIND=idp) :: CDA, CFD
      REAL(KIND=idp), PARAMETER :: THRESH=1.0E-08_idp
      INTEGER :: I, IC, ID, IORB1, IORB2, J, JA, JB, M, MA, MAA, MB, &
     &           MBB, MDA, MDB, N
!       INTEGER :: ISRCHE
      INTEGER, DIMENSION(2) :: NDD,II, JJ,CPQ,NPQ
      INTEGER, DIMENSION(nsrb) :: NDTA
      INTEGER, DIMENSION(nsrb+1) :: NDTB
      INTEGER, DIMENSION(4) :: NDTC
      INTEGER, DIMENSION(nelt) :: NDTD
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     TMTMA -  COMPUTES OVERLAPS BETWEEN CSF PAIRS.
!
!              FOR OFF-DIAGONAL ELEMENTS OF THE DENSITY MATRIX, PAIRS OF
!              CSFS MUST BE CONSIDERED. WITHIN THIS OF COURSE IS IMPLIED
!              A SUMMATION OVER THE PAIRS OF DETERMINANTS WHICH
!              CONSTITUTE THE CSFS. THIS PROCEDURE IS ALSO REQUIRED IN
!              THE CONSTRUCTION OF THE HAMILTONIAN MATRIX AND THERE IS A
!              SUBROUTINE ENRGMA IN PROGRAM SPEEDY WHICH PERFORMS THIS
!              TASK. TMTMA IS A DIRECT ANALOGUE OF THAT ROUTINE BUT IS
!              NECESSARILY SIMPLER BECAUSE OF THE FACT THAT DETERMINANT
!              PAIRS MUST DIFFER BY NO MORE THAN TWO SPIN ORBITALS IN
!              ORDER TO CONTRIBUTE TO THE FIRST ORDER DENSITY MATRIX.
!
!     INPUT DATA :
!           NODA  NUMBER OF DETERMINANTS IN THE FIRST CSF
!             CA  COEFFICIENTS FOR EACH DETERMINANT IN THE FIRST CSF
!            NDA  DETERMINANTS IN THE FIRST CSF
!           NODB  NUMBER OF DETERMINANTS IN THE SECOND CSF
!             CB  COEFFICIENTS FOR EACH DETERMINANT IN THE SECOND CSF
!            NDB  DETERMINANTS IN THE SECOND CSF
!           NELT  NUMBER OF ELECTRON IN THE SYSTEM
!             MN  TABLE GIVING THE ORBITAL NUMBER FOR EACH SPIN ORBITAL
!             MS  TABLE GIVING SZ FOR EACH SPIN ORBITAL
!           MDTR  POINTER ARRAY LINKING THE REFERENCE DETERMINANT TO
!                 ALL SPIN ORBITALS
!
!     DATA FOR THE OUTPUT ROUTINE (SEE COMMENTS THEREIN):
!
!       NFTD,NELMT,NBLK,NEL,LENGTH,II,JJ,NPQ,CPQ,ICSF,JCSF,IWRITE
!
!     LOCAL DATA :
!           NDTC  WORKSPACE
!            NDD  WORKSPACE
!
!***********************************************************************
!
!-----COPY REFERENCE DETERMINANT INTO NDTD AND ZEROIZE THE COEFFICIENT
!     FOR THIS PAIR OF CSFS
!
      DO I=1, NELT
         NDTD(I)=NDTR(I)
      END DO
!
!-----OUTER LOOP OVER DETERMINANTS
!
      MDA=1
      DO IC=1, NODA
         CDA=CA(IC)
         MA=NDA(MDA)
         MAA=MDA+MA
         MDB=1
!
!--------INNER LOOP OVER DETERMINANTS
!
         DO ID=1, NODB
            CFD=CDA*CB(ID)
            MB=NDB(MDB)
            MBB=MDB+MB
!
!--------BASED UPON THE NUMBER OF REPLACEMENTS FROM THE REFERENCE
!        DETERMINANT, JUMP TO APPROPRIATE PIECE OF CODE.
!          A. BOTH DETERMINANTS ARE NOT THE REF. DET. GOTO 200
!          B. ONE IS THE REF DET. OTHER HAS > ONE REPLS. GOTO 300
!          C. BOTH ARE THE REF. DET. GOTO 270
!          D. ONE IS THE REF. DET. OTHER HAS ONE REPL. GOTO 260
!
            IF(MA*MB.NE.0)GO TO 200
            IF(MAX(MA,MB).GT.1)GO TO 302
            IF(MA.EQ.MB)GO TO 270

            IF(MA.EQ.0)THEN
               IF(MB.LE.1)THEN
                  NDTC(3)=NDB(MDB+1)
                  NDTC(4)=NDB(MBB+1)
                  GO TO 260

               END IF

            ELSE
               IF(MA.LE.1)THEN
                  NDTC(3)=NDA(MAA+1)
                  NDTC(4)=NDA(MDA+1)
                  GO TO 260

               END IF

            END IF
!
!---------PROCESS THE DIFFERENT SITUATIONS
!
!---------NEITHER IS THE REFERENCE DETERMINANT. THE FOLLOWING LINES OF C
!         IMPLEMENT AN ALGORITHM DESCRIBED BY NESBET IN HIS BOOK ON
!         VARIATIONAL METHODS IN ELECTRON-ATOM SCATTERING THEORY. THE
!         ALGORITHM EXPRESSES ON DETERMINANT AS AN EXCITATION OF THE OTH
!         RATHER THAN BOTH BEING EXCITATIONS RELATIVE TO THE REFERENCE
!         DETERMINANT. FOR A FULLER DESCRIPTION OF THE ALGORITHM SEE PAG
!         OF THE BOOK.
!
!
!--------- STEP ONE: COPY ALL REPLACED SPIN-ORBS INTO NDTB, AND
!                    ALL REPLACEMENT SPIN-ORBS INTO NDTA
!density_matrix_orbital_pairs
 200        DO I=1, MA
               NDTA(I)=NDA(MAA+I)
               NDTB(I)=NDA(MDA+I)
            END DO
!
            JA=MA
            DO I=1, MB
               J=ISRCHE(MA,NDTB,1,NDB(MDB+I))
               IF(J.LE.MA)THEN
                  NDTB(J)=NDB(MBB+I)

               ELSE
                  JA=JA+1
                  NDTA(JA)=NDB(MDB+I)
                  NDTB(JA)=NDB(MBB+I)

               END IF

            END DO
!
            IF((JA-MA).GT.1)GO TO 302
            JB=0
!
            DO I=1, JA
               J=ISRCHE(JA,NDTB,1,NDTA(I))
               IF(J.LE.JA)THEN
                  IF(I.NE.J)THEN
                     NDTB(J)=NDTB(I)
                     NDTB(I)=NDTA(I)
                     CFD=-CFD

                  END IF

               ELSE
                  JB=JB+1
                  IF(JB.GT.1)GO TO 302
                  NDD(JB)=I

               END IF

            END DO
!
            IF(JB.EQ.0)GO TO 270

            IF(JB.EQ.1)THEN
               NDTC(3)=NDTA(NDD(1))
               NDTC(4)=NDTB(NDD(1))
               GO TO 260

            END IF
!
!........ ONE PAIR IS DIFFERENT
!
!.......... CHECK FOR IDENTICAL SPIN VALUES
            IF(MS(NDTC(3)).NE.MS(NDTC(4)))GO TO 302
!
 260        IF(MN(NDTC(3)).LT.MN(NDTC(4)))THEN
               IORB1=MN(NDTC(3))
               IORB2=MN(NDTC(4))

            ELSE
               IORB1=MN(NDTC(4))
               IORB2=MN(NDTC(3))

            END IF

            if (abs(CFD).gt.THRESH) then
               density_matrix_coefficients(ICSF,JCSF)=density_matrix_coefficients(ICSF,JCSF)+CFD
               density_matrix_orbital_pairs(ICSF,JCSF,1) = IORB1
               density_matrix_orbital_pairs(ICSF,JCSF,2) = IORB2

            end if

            GO TO 302
!
!......... DETERMINANTS ARE IDENTICAL : THIS CASE DOES HAPPEN AFTER ALL!
!
 270        CONTINUE
!
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=NDA(MAA+I)

            END DO
!
            DO M=1, NELT
               N=NDTD(M)
               IORB1=MN(N)
               IORB2=IORB1

!                write(951,'(i4,i4,d20.5)') iorb1,iorb2,cfd !Alex test

              if (abs(CFD).gt.THRESH) then
                 density_matrix_coefficients(ICSF,JCSF)=density_matrix_coefficients(ICSF,JCSF)+CFD
                 density_matrix_orbital_pairs(ICSF,JCSF,1) = IORB1
                 density_matrix_orbital_pairs(ICSF,JCSF,2) = IORB2

              end if

            END DO
!
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=N

            END DO
!
!......... ASCEND IN THE COUPLING LOOPS
!
 302        MDB=MBB+MB+1

         END DO
         MDA=MAA+MA+1

      END DO
!
      RETURN
!
      END SUBROUTINE TMTMA




      SUBROUTINE TMTMA_SPARSE(NODA,CA,NDA,NODB,CB,NDB,nsrb,NELT,NDTR,MN,MS,&
     &                 MDTR,&!NFTD,NELMT,NBLK,NEL,LENGTH,II,JJ,NPQ,CPQ,&
     &                 ICSF,JCSF,density_matrix_coefficients, density_matrix_orbital_pairs_1,&
                       density_matrix_orbital_pairs_2,orbs_diag)
      USE cdenprop_defs, ONLY : idp ! for specifying the kind of reals
      use class_COOSparseMatrix_real
      use class_COOSparseMatrix_integer
      use omp_lib
      IMPLICIT NONE
! C
! C*** Start of declarations rewritten by SPAG
! C
! C Dummy arguments
! C
      INTEGER :: ICSF, JCSF, LENGTH, NBLK, NEL, NELMT, NELT, NFTD, NODA, &
     &           NODB, NSRB
      REAL(KIND=idp), DIMENSION(*) :: CA, CB

type(COOMatrix_real) :: density_matrix_coefficients
type(COOMatrix_integer) :: density_matrix_orbital_pairs_1 , density_matrix_orbital_pairs_2
integer, dimension(NELT,*):: orbs_diag


!       REAL(KIND=idp), DIMENSION(:,:) :: density_matrix_coefficients
!       INTEGER, DIMENSION(:,:,:) :: density_matrix_orbital_pairs
      INTEGER, DIMENSION(*) ::  MDTR, MN, MS, NDA, NDB, NDTR
!       INTEGER, DIMENSION(2,*) :: NPQ
!ZM added NDB to the intent(in) list
      INTENT (IN) CA, CB, MDTR, MN, MS, NDA, NDTR, NELT, NODA, NODB, &
     &            NSRB, NDB
! C
! C Local variables
! C

      REAL(KIND=idp) :: coefficient_temp
      INTEGER :: iorb1_temp,iorb2_temp

      REAL(KIND=idp) :: CDA, CFD
      REAL(KIND=idp), PARAMETER :: THRESH=1.0E-08_idp
      INTEGER :: I, IC, ID, IORB1, IORB2, J, JA, JB, M, MA, MAA, MB, &
     &           MBB, MDA, MDB, N
!       INTEGER :: ISRCHE
      INTEGER, DIMENSION(2) :: NDD,II, JJ,CPQ,NPQ
      INTEGER, DIMENSION(nsrb) :: NDTA
      INTEGER, DIMENSION(nsrb+1) :: NDTB
      INTEGER, DIMENSION(4) :: NDTC
      INTEGER, DIMENSION(nelt) :: NDTD
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     TMTMA -  COMPUTES OVERLAPS BETWEEN CSF PAIRS.
!
!              FOR OFF-DIAGONAL ELEMENTS OF THE DENSITY MATRIX, PAIRS OF
!              CSFS MUST BE CONSIDERED. WITHIN THIS OF COURSE IS IMPLIED
!              A SUMMATION OVER THE PAIRS OF DETERMINANTS WHICH
!              CONSTITUTE THE CSFS. THIS PROCEDURE IS ALSO REQUIRED IN
!              THE CONSTRUCTION OF THE HAMILTONIAN MATRIX AND THERE IS A
!              SUBROUTINE ENRGMA IN PROGRAM SPEEDY WHICH PERFORMS THIS
!              TASK. TMTMA IS A DIRECT ANALOGUE OF THAT ROUTINE BUT IS
!              NECESSARILY SIMPLER BECAUSE OF THE FACT THAT DETERMINANT
!              PAIRS MUST DIFFER BY NO MORE THAN TWO SPIN ORBITALS IN
!              ORDER TO CONTRIBUTE TO THE FIRST ORDER DENSITY MATRIX.
!
!     INPUT DATA :
!           NODA  NUMBER OF DETERMINANTS IN THE FIRST CSF
!             CA  COEFFICIENTS FOR EACH DETERMINANT IN THE FIRST CSF
!            NDA  DETERMINANTS IN THE FIRST CSF
!           NODB  NUMBER OF DETERMINANTS IN THE SECOND CSF
!             CB  COEFFICIENTS FOR EACH DETERMINANT IN THE SECOND CSF
!            NDB  DETERMINANTS IN THE SECOND CSF
!           NELT  NUMBER OF ELECTRON IN THE SYSTEM
!             MN  TABLE GIVING THE ORBITAL NUMBER FOR EACH SPIN ORBITAL
!             MS  TABLE GIVING SZ FOR EACH SPIN ORBITAL
!           MDTR  POINTER ARRAY LINKING THE REFERENCE DETERMINANT TO
!                 ALL SPIN ORBITALS
!
!     DATA FOR THE OUTPUT ROUTINE (SEE COMMENTS THEREIN):
!
!       NFTD,NELMT,NBLK,NEL,LENGTH,II,JJ,NPQ,CPQ,ICSF,JCSF,IWRITE
!
!     LOCAL DATA :
!           NDTC  WORKSPACE
!            NDD  WORKSPACE
!
!***********************************************************************
!
!-----COPY REFERENCE DETERMINANT INTO NDTD AND ZEROIZE THE COEFFICIENT
!     FOR THIS PAIR OF CSFS
!
      coefficient_temp =0._idp
      iorb1_temp = 0
      iorb2_temp = 0

      DO I=1, NELT
         NDTD(I)=NDTR(I)
      END DO
!
!-----OUTER LOOP OVER DETERMINANTS
!
      MDA=1
      DO IC=1, NODA
         CDA=CA(IC)
         MA=NDA(MDA)
         MAA=MDA+MA
         MDB=1
!
!--------INNER LOOP OVER DETERMINANTS
!
         DO ID=1, NODB
            CFD=CDA*CB(ID)
            MB=NDB(MDB)
            MBB=MDB+MB
!
!--------BASED UPON THE NUMBER OF REPLACEMENTS FROM THE REFERENCE
!        DETERMINANT, JUMP TO APPROPRIATE PIECE OF CODE.
!          A. BOTH DETERMINANTS ARE NOT THE REF. DET. GOTO 200
!          B. ONE IS THE REF DET. OTHER HAS > ONE REPLS. GOTO 300
!          C. BOTH ARE THE REF. DET. GOTO 270
!          D. ONE IS THE REF. DET. OTHER HAS ONE REPL. GOTO 260
!
            IF(MA*MB.NE.0)GO TO 200
            IF(MAX(MA,MB).GT.1)GO TO 302
            IF(MA.EQ.MB)GO TO 270
            IF(MA.EQ.0)THEN
               IF(MB.LE.1)THEN
                  NDTC(3)=NDB(MDB+1)
                  NDTC(4)=NDB(MBB+1)
                  GO TO 260

               END IF

            ELSE
               IF(MA.LE.1)THEN
                  NDTC(3)=NDA(MAA+1)
                  NDTC(4)=NDA(MDA+1)
                  GO TO 260

               END IF

            END IF
!
!---------PROCESS THE DIFFERENT SITUATIONS
!
!---------NEITHER IS THE REFERENCE DETERMINANT. THE FOLLOWING LINES OF C
!         IMPLEMENT AN ALGORITHM DESCRIBED BY NESBET IN HIS BOOK ON
!         VARIATIONAL METHODS IN ELECTRON-ATOM SCATTERING THEORY. THE
!         ALGORITHM EXPRESSES ON DETERMINANT AS AN EXCITATION OF THE OTH
!         RATHER THAN BOTH BEING EXCITATIONS RELATIVE TO THE REFERENCE
!         DETERMINANT. FOR A FULLER DESCRIPTION OF THE ALGORITHM SEE PAG
!         OF THE BOOK.
!
!
!--------- STEP ONE: COPY ALL REPLACED SPIN-ORBS INTO NDTB, AND
!                    ALL REPLACEMENT SPIN-ORBS INTO NDTA
!density_matrix_orbital_pairs
 200        DO I=1, MA
               NDTA(I)=NDA(MAA+I)
               NDTB(I)=NDA(MDA+I)

            END DO
!
            JA=MA
            DO I=1, MB
               J=ISRCHE(MA,NDTB,1,NDB(MDB+I))
               IF(J.LE.MA)THEN
                  NDTB(J)=NDB(MBB+I)

               ELSE
                  JA=JA+1
                  NDTA(JA)=NDB(MDB+I)
                  NDTB(JA)=NDB(MBB+I)

               END IF

            END DO
!
            IF((JA-MA).GT.1)GO TO 302
            JB=0
!
            DO I=1, JA
               J=ISRCHE(JA,NDTB,1,NDTA(I))
               IF(J.LE.JA)THEN
                  IF(I.NE.J)THEN
                     NDTB(J)=NDTB(I)
                     NDTB(I)=NDTA(I)
                     CFD=-CFD

                  END IF

               ELSE
                  JB=JB+1
                  IF(JB.GT.1)GO TO 302
                  NDD(JB)=I

               END IF

            END DO
!
            IF(JB.EQ.0)GO TO 270
            IF(JB.EQ.1)THEN
               NDTC(3)=NDTA(NDD(1))
               NDTC(4)=NDTB(NDD(1))
               GO TO 260

            END IF
!
!........ ONE PAIR IS DIFFERENT
!
!.......... CHECK FOR IDENTICAL SPIN VALUES
            IF(MS(NDTC(3)).NE.MS(NDTC(4)))GO TO 302
!
 260        IF(MN(NDTC(3)).LT.MN(NDTC(4)))THEN
               IORB1=MN(NDTC(3))
               IORB2=MN(NDTC(4))

            ELSE
               IORB1=MN(NDTC(4))
               IORB2=MN(NDTC(3))

            END IF

            if (abs(CFD).gt.THRESH) then
               coefficient_temp = coefficient_temp + CFD
               iorb1_temp = IORB1
               iorb2_temp = IORB2


            end if
            GO TO 302
!
!......... DETERMINANTS ARE IDENTICAL : THIS CASE DOES HAPPEN AFTER ALL!
!
 270        CONTINUE
!
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=NDA(MAA+I)
            END DO
!
            DO M=1, NELT
               N=NDTD(M)
               IORB1=MN(N)
               IORB2=IORB1

!                write(951,'(i4,i4,d20.5)') iorb1,iorb2,cfd

               ! needed when considering the case where icsf jcsf are equal.
               if (abs(CFD).gt.THRESH) then
                 if (icsf.eq.jcsf) then
                    orbs_diag(M,icsf) = IORB2
                 else
                    coefficient_temp=coefficient_temp+cfd
                 endif
              end if


            END DO
!
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=N

            END DO
!
!......... ASCEND IN THE COUPLING LOOPS
!
 302        MDB=MBB+MB+1

         END DO
         MDA=MAA+MA+1

      END DO
!
! SPARSE MATRIX INSERT
      if (abs(coefficient_temp).gt.THRESH) then
         !$OMP CRITICAL
         call density_matrix_coefficients%insert_single(ICSF,JCSF,coefficient_temp)
         call density_matrix_orbital_pairs_1%insert_single(ICSF,JCSF,iorb1_temp)
         call density_matrix_orbital_pairs_2%insert_single(ICSF,JCSF,iorb2_temp)
         !$OMP END CRITICAL
      end if

      RETURN
!
      END SUBROUTINE TMTMA_SPARSE




      SUBROUTINE MAKEMG(MG,NSRB,NELT,NDTRF)
!       USE precisn, ONLY : wp ! for specifying the kind of reals
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: NELT, NSRB
      INTEGER, DIMENSION(nsrb) :: MG
      INTEGER, DIMENSION(nelt) :: NDTRF
      INTENT (IN) NDTRF, NELT, NSRB
      INTENT (OUT) MG
!
! Local variables
!
      INTEGER :: I
!
!*** End of declarations rewritten by SPAG
!
!     ******************************************************************
!
!     ASSOCIATES WITH EACH SPIN ORBITAL IN THE REFERENCE DETERMINANT
!     A SEQUENCE NUMBER GIVING THE CONTOGUOUS POSITION OF THAT SPIN
!     ORBITAL IN THE REFERENCE DETERMINANT. NOTE THAT THIS ROUTINE IS
!     COPIED FROM SPEEDY.
!
!     INPUT DATA:
!           NSRB THE NUMBER OF SPIN ORBITALS IN THE WAVEFUNCTION SET
!           NELT THE NUMBER OF ELECTRONS IN THE WAVEFUNCTION
!          NDTRF THE REFERENCE DETERMINANT FOR THE WAVEFUNCTION
!
!     OUTPUT DATA:
!              MG THE POINTER ARRAY ITSELF.
!
!     ******************************************************************
!
      DO I=1, NSRB
         MG(I)=0

      END DO
!
      DO I=1, NELT
         MG(NDTRF(I))=I

      END DO
!
      RETURN
!
      END SUBROUTINE MAKEMG


      SUBROUTINE INDEXX(N,ARRIN,INDX)
!     ******************************************************************
!
!     INDEXX - Takes an array and produces a set of indices such
!              that ARRIN(INDX(J)) is in ascending order J=1,2,..,N
!
!     Input data:
!              N number of elements in the array to be ordered
!          ARRIN R*8 array which is to placed in ascending order.
!
!     Output data:
!            INDX a set of indices for ascending indices
!
!     Notes:
!
!     This routine is taken from the book Numerical Receipes by
!     Press, Flannery, Teukolsky and Vetterling Chapter 8 p. 233.
!     ISBN 0-521-30811-9 pub. Cambridge University Press (1986)
!     QA297.N866
!
!     This routines has been adapted by Charles J Gillan for use
!     in the R-matrix codes.
!
!     ******************************************************************
!       USE precisn, ONLY : wp ! for specifying the kind of reals
      IMPLICIT NONE
!
!     *** Start of declarations rewritten by SPAG
!
!     PARAMETER definitions
!
      REAL(KIND=idp), PARAMETER :: VSMALL=1.0E-20_idp
!
!     Dummy arguments
!
      INTEGER :: N
      REAL(KIND=idp), DIMENSION(n) :: ARRIN
      INTEGER, DIMENSION(n) :: INDX
      INTENT (IN) ARRIN, N
      INTENT (INOUT) INDX
!
!     Local variables
!
      INTEGER :: I, INDXT, IR, J, L
      REAL(KIND=idp) :: Q
!
!     *** End of declarations rewritten by SPAG
!
!
      if (N.eq.1) then
         INDX(1)=1
         return
      end if
!
!
!     Initialize the index array with consecutive integers
!
      DO J=1, N
         INDX(J)=J
      END DO
!
      L=N/2+1
      IR=N
!
!     From here on the algorithm is HEAPSORT wit indirect addressing
!     through INDX in all references to ARRIN
!
 10   CONTINUE
!
      IF(L.GT.1)THEN
         L=L-1
         INDXT=INDX(L)
         Q=ARRIN(INDXT)
      ELSE
         INDXT=INDX(IR)
         Q=ARRIN(INDXT)
         INDX(IR)=INDX(1)
         IR=IR-1
         IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            GO TO 800
         END IF
      END IF
!
      I=L
      J=L+L
!
 20   CONTINUE
!
      IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1))+VSMALL)J=J+1
         END IF
         IF(Q.LT.ARRIN(INDX(J))+VSMALL)THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         END IF
         GO TO 20
      END IF
      INDX(I)=INDXT
!
!     Loop back to process another element
!
      GO TO 10
!
 800  CONTINUE
!
      RETURN
!
      END SUBROUTINE INDEXX


!! END OF OLD DENPROP ROUTINES

end module cdenprop_aux
