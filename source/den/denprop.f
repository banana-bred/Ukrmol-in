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
!*==denprop.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      PROGRAM DENPROP
      USE GLOBAL_UTILS, ONLY: PRINT_UKRMOL_HEADER
      USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
      USE UKRMOL_INTERFACE_GBL, ONLY: START_MPI, FINALIZE_MPI
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      CHARACTER(LEN=8) :: CURDAT
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     DENPROP is a program to evaluate the following quantum
C             mechanical entities from a wavefunction or pair
C             of wavefunctions:
C
C        a) first order spin-reduced density matrices,
C        b) first order spin-reduced transition density matrices
C        c) wavefunction properties (moments)
C        d) transition moments
C
C     Computed properties and transition moments may be written out to
C     a target properties file. This file is read by the external region
C     program the data obtained therein used to construct the asymptotic
C     scattering potential.
C
C     Additionally the density matrix or matrices generated may be store
C     in a library dataset. They can be read by the pseudo natural
C     orbital code.
C
C     Author:
C
C        Charles J Gillan, Queen's University of Belfast, Feb 23rd 1989
C
C     Modifications:
C
C     1. Jonathan Tennyson   Nov 1991
C
C            Updated for better integration with the
C           (ALCHEMY I based) molecular R-matrix package
C
C     2. Jonathan Tennyson   May 1992
C
C            Subroutine POLAR added to calculate dipole polarizabitilies
C
C     3. Charles J Gillan    August 1993
C
C            Upgraded to handle Abelian point groups as well as C-inf-v
C            symmetry (linear) molecules; additionally more robust error
C            checking introduced in various subroutines.
c
C     4. Lesley Morgan       November 1997
C            Interfaced to SWEDEN based codes.  Some bugs fixed and o/p
C             rationalized
C
C     5. Zdenek Masin        July 2014
C            Interfaced to UKRmol+ codes.
C
C     5. Daniel Darby-Lewis        May 2019
C            Interfaced to new polarisiablity calculation in cdenprop.
C
C***********************************************************************
C
c
C---- Write the program title on the run output
C
      CALL PRINT_UKRMOL_HEADER(OUTPUT_UNIT)
      WRITE(*,12)
C
C----- Invoke subroutines to date stamp the run
C
      CALL DATE_and_TIME(CURDAT)
      WRITE(*,13)CURDAT
C
C---- Invoke the driver to perform the computation
C
      CALL START_MPI
      CALL DENXDR
C
C---- Notify the user of successful completion
C
      WRITE(*,8000)
C
      CALL FINALIZE_MPI
      STOP
C
C---- Format Statements
C
 12   FORMAT(//,20X,' Density Matrix package',/,20X,' PROGRAM DENPROP',
     &       /)
 13   FORMAT(20X,' DATE =',A)
C
 8000 FORMAT(//,10X,'D E N P R O P   has completed execution ',//)
C
      END PROGRAM DENPROP
!*==addptab.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE ADDPTAB(ISTATE,JSTATE,IPCDE,PROPS,NPROPS,MAXPROPS,IST,
     &                   JST,NUMBER,NEWCODES,XNEWPROP,IWRITE)
C***********************************************************************
C
C     ADDPTAB - Adds data to the Properties TABle
C
C     Input data:
C            IST  I state designation from unique target st. table
C            JST  J state designation from unique target st. table
C         NUMBER  Number of new peoperties to be added to Prp. table
C       NEWCODES  Operator codes for the new properties to be added
C       XNEWPROP  Expectation values for the new properties
C         IWRITE  Logical unit for the printer
C
C     Input/Output data:
C         ISTATE  Designation of state I - the bra vector
C         JSTATE  Designation of state J - the ket vector
C                 where both above are wrt the unique target state table
C          IPCDE  Slater property integral codes (8 per property)
C          PROPS  Expectation values for the wavefunction pair.
C         NPROPS  Number of properties in the table on entry and then
C                 on exit.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IST, IWRITE, JST, MAXPROPS, NPROPS, NUMBER
      INTEGER, DIMENSION(8,maxprops) :: IPCDE
      INTEGER, DIMENSION(maxprops) :: ISTATE, JSTATE
      INTEGER, DIMENSION(8,*) :: NEWCODES
      REAL(KIND=wp), DIMENSION(maxprops) :: PROPS
      REAL(KIND=wp), DIMENSION(*) :: XNEWPROP
      INTENT (IN) IST, IWRITE, JST, MAXPROPS, NEWCODES, NUMBER, XNEWPROP
      INTENT (OUT) IPCDE, ISTATE, JSTATE, PROPS
      INTENT (INOUT) NPROPS
C
C Local variables
C
      INTEGER :: I, ITEMP, J
C
C*** End of declarations rewritten by SPAG
C
c        WRITE(IWRITE,1000)
c        WRITE(IWRITE,1010) NPROPS,MAXPROPS,NUMBER
c        WRITE(IWRITE,1015)
c        DO 5 I=1,NUMBER
c        WRITE(IWRITE,1020) (NEWCODES(J,I),J=1,8),XNEWPROP(1,I)
c   5    CONTINUE
C
C---- Make sure that there is enough room in the table to add the new
C     entries
C
      ITEMP=NPROPS+NUMBER
C
      IF(ITEMP.GT.MAXPROPS)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9910)MAXPROPS, ITEMP
         STOP
      END IF
C
C---- Add the new data into the table
C
      DO I=1, NUMBER
         ISTATE(NPROPS+I)=IST
         JSTATE(NPROPS+I)=JST
      END DO
C
C---- There are NUMBER new Slater values to be added to the array
C
      DO I=1, NUMBER
         DO J=1, 8
            IPCDE(J,NPROPS+I)=NEWCODES(J,I)
         END DO
      END DO
C
C---- Similarly PROPS is augumented.
C
      DO I=1, NUMBER
         PROPS(NPROPS+I)=XNEWPROP(I)
      END DO
 
C
C---- Re-adjust the number of entries in the table now
C
      NPROPS=ITEMP
C
      RETURN
C
C---- Format statements
C
 1000 FORMAT(//,15X,'====> ADDPTAB - ADD PROPERTIES TO TABLE <====',/)
 1010 FORMAT(/,15X,'No. of properties already in table = ',I5,/,15X,
     &       'Maximum permitted properties       = ',I5,/,15X,
     &       'Number of properties to be added   = ',I5,/)
 1015 FORMAT(/,15X,'Property Op. Codes  Properties ',/,15X,
     &       '------------------  ---------- ',/)
 1020 FORMAT(17X,8I2,2X,F10.4)
C
 9900 FORMAT(//,10X,'**** Error in ADDPTAB:',//)
 9910 FORMAT(10X,'Not enough rows in Properties Table ',/,10X,
     &       'Given = ',I10,' Need = ',I10,/)
C
      END SUBROUTINE ADDPTAB
!*==addtst.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE ADDTST(NTGT,ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,IWRITE,
     &                  newst,MAXTGT,ISPINN,ISZN,MGVNN,GUTOTN,IREFLN,
     &                  ENERGYN)
C***********************************************************************
C
C     ADDTST - ADD details of a Target STate to the unique state table
C              and return the row number. In the event that the target
C              already exists then just return the row number.
C
C     Input data:
C          ISPIN 2*S+1 for each target state
C            ISZ 2*Sz for each target state
C           MGVN Lamda value (C-inf-v) or Irred. Rep of each state
C          GUTOT For D-inf-h only, the gerade or ungerade value
C          IREFL For C-inf-v only the +/- sigma character
C         ENERGY Eigen-energy in Hartrees for the each state
C         IWRITE Logical unit for the printer
C         MAXTGT Maximum nuber of rows in the target state table
C
C         Values to be added for the new state:
C
C         ISPINN
C           ISZN
C          MGVNN
C         GUTOTN
C         IREFLN
C        ENERGYN
C
C     Output data:
C            IROW The row which contains this new target state. Note
C                 that if the state existed already then IROW < NTGT.
C
C     Input/Output data:
C                  NTGT On input and output this is the number of rows
C                       used in the table. If a new target state is
C                       added then
C
C                          NTGT (output) = NTGT(input) + 1
C
C                       else it it already existed then
C
C                          NTGT(output) = NTGT(input)
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL(KIND=wp), PARAMETER :: VSMALL=1.0E-08_wp
C
C Dummy arguments
C
      REAL(KIND=wp) :: ENERGYN
      INTEGER :: GUTOTN, IREFLN, ISPINN, ISZN, IWRITE, MAXTGT, MGVNN, 
     &           NEWST, NTGT
      REAL(KIND=wp), DIMENSION(maxtgt) :: ENERGY
      INTEGER, DIMENSION(maxtgt) :: GUTOT, IREFL, ISPIN, ISZ, MGVN
      INTENT (IN) ENERGYN, GUTOTN, IREFLN, ISPINN, ISZN, IWRITE, MAXTGT, 
     &            MGVNN
      INTENT (OUT) NEWST
      INTENT (INOUT) ENERGY, GUTOT, IREFL, ISPIN, ISZ, MGVN, NTGT
C
C Local variables
C
      INTEGER :: I, K, NROW
      INTEGER, DIMENSION(ntgt) :: KROW
C
C*** End of declarations rewritten by SPAG
C
C---- Scan the table in existence already to see if the energy ENERGYN
C     already exists.
C
      nROW=0
C
 
!      write(172,*)(energy(i),i=1,maxtgt)
      DO I=1, NTGT
         IF(ABS(ENERGY(I)-ENERGYN).LT.VSMALL)THEN
            nROW=nROW+1
            krow(nrow)=i
         END IF
      END DO
C
C..... If the energy was not found then we proceed to add a new row
C
C      Otherwise proceed to check more carefully for degeneracy
C
      DO i=1, nrow
         k=krow(i)
C
C---- To match an existing row we must attempt to match all columns
c
         IF(ISPINN.EQ.ISPIN(k) .AND. ISZN.EQ.ISZ(k) .AND. 
     &      MGVNN.EQ.MGVN(k) .AND. GUTOTN.EQ.GUTOT(k) .AND. 
     &      IREFLN.EQ.IREFL(k))THEN
            newst=k
            GO TO 800
         END IF
      END DO
C
C---- Add a new row to the table, if space permits !
c
      NTGT=NTGT+1
      newst=NTGT
C
      IF(NTGT.GT.MAXTGT)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9910)MAXTGT, NTGT
         STOP
      END IF
C
      ISPIN(NTGT)=ISPINN
      ISZ(NTGT)=ISZN
      MGVN(NTGT)=MGVNN
      GUTOT(NTGT)=GUTOTN
      IREFL(NTGT)=IREFLN
      ENERGY(NTGT)=ENERGYN
C
!     Hemal Varambhia polarisability write statement
      WRITE(172,*)(mgvn(i),i=1,maxtgt)
      WRITE(172,*)(energy(i),i=1,maxtgt)
 800  RETURN
C
C---- Format Statements
C
 9900 FORMAT(//,10X,'**** Error in ADDTST :',/)
 9910 FORMAT(10X,'Need more rows in table',/,10X,'Current maximum   = ',
     &       I10,/,10X,'Trying to add row = ',I10,/)
C
      END SUBROUTINE ADDTST
!*==catche.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE CATCHE(NTGT,NFTSOR,NTGTF,NTGTS,NTGTL,IWRITE)
C***********************************************************************
C
C     CATCHE - CATCHes Errors in the namelist &INPUT and terminates
C              the run should any be found.
C
C     Input data: (This constitutes &INPUT excluding the NAME variable)
C           NTGT  Number of target states in the problem.
C         NFTSOR  Logical unit numbers for determinant information for
C                 all of the wavefunctions used.
C          NPFLG  Set of print flags for verious stages of the calculati
C         IWRITE  Logical unit for the printer
C           NFTD  Logical unit for the unsorted density matrix formulae
C           NFTA   " " "   ""   "   "  sorted density matrix formulae
C           NFDA   " " "   ""   "   "  direct access file used by the so
C           LCOF
C           LDAR  Number of inetegr words per record on unit NFDA
C           NFTG  Logical unit for the CI vectors defining each
C                 wavefunction.
C         NFTINT  Logical unit for the transfromed property interals tha
C                 were generated by TRANS.
C          NFTMT
C          NMSET  Set number at which the generated moments are to be ou
C          NFTDL  Logical unit for the density matrix library
C         NUCCEN  Number of the nucleus defined as the scattering center
C            ISW  Type of output moments requested
C           IPOL  Switch to carry out polarizibility calculation
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, NTGT
      INTEGER, DIMENSION(ntgt) :: NFTSOR, NTGTF, NTGTL, NTGTS
      INTENT (IN) IWRITE, NFTSOR, NTGT, NTGTF
      INTENT (INOUT) NTGTL, NTGTS
C
C Local variables
C
      INTEGER :: I, J
C
C*** End of declarations rewritten by SPAG
C
C---- Entries in NFTSOR must be unique - That means that all CSF
C     expansions are treated once only. CI vector pick up in module
C     DENMAK can't handle anything else and in any case it would be
C     a sever computational waste - double computation !
C
      DO I=1, NTGT
         DO J=I+1, NTGT
            IF(NFTSOR(J).EQ.NFTSOR(I))THEN
               WRITE(IWRITE,990)
               WRITE(IWRITE,988)I, J, NFTSOR(I)
               CALL TMTCLOS()
            END IF
         END DO
      END DO
C
C---- Make sure that all entries in NFTSOR, NTGTF and NTGTS are
C     positive
C
      DO I=1, NTGT
         IF(ntgts(i).NE.1)THEN
            WRITE(iwrite,800)i
            ntgts(i)=1
         END IF
         IF(NFTSOR(I).LE.0 .OR. NTGTF(I).LE.0 .OR. NTGTS(I).LE.0)THEN
            WRITE(IWRITE,990)
            WRITE(IWRITE,992)I, NFTSOR(I), NTGTF(I), NTGTS(I)
            CALL TMTCLOS()
         END IF
      END DO
C
C---- The user may not have set NTGTL at all. In this case the entries
C     will be at their default of zero and must be redefined to be
C     equal to NTGTS.
C
      IF(NTGTL(1).EQ.0)THEN
         DO I=1, NTGT
            NTGTL(I)=NTGTS(I)
         END DO
      END IF
C
C---- Now check that NTGTL(I) > or = to NTGTS(I) for all NTGT entries
C
      DO I=1, NTGT
         IF(NTGTL(I).LT.NTGTS(I))THEN
            WRITE(IWRITE,990)
            WRITE(IWRITE,994)I, NTGTS(I), NTGTL(I)
            CALL TMTCLOS()
         END IF
      END DO
C
      WRITE(6,8000)
C
      RETURN
C
C---- Error handler - shut down code due to error in data
C
      CONTINUE
C
      CALL TMTCLOS()
C
C---- Format Statements
C
 988  FORMAT(5X,'Double entry in NFTSOR array. ',//,5X,'Elements ',I5,
     &       ' and ',I5,' are equal to ',I5,//,5X,
     &       'This cannot be handled by code - Combine all  ',/,5X,
     &       'CI vectors associated with this CSF expansion ',/,5X,
     &       'into one input set.',/)
 990  FORMAT(//,5X,'**** Error in &INPUT data detected by routine',
     &       ' CATCHE',//)
 992  FORMAT(5X,'Entry number = ',I3,' in one of NFTSOR,NTGTF,NTGTS',/,
     &       5X,'less than 1. Vaues = ',3(I7,1X),/)
 994  FORMAT(5X,'Sequence number error for vector pickup. ',/,5X,
     &       'Entry number ',I3,' for NTGTS,NTGTL in error ',//,5X,
     &       'NTGTS (first vector required) = ',I6,/,5X,
     &       'NTGTL (last  vector required) = ',I6,/)
 8000 FORMAT(//,5X,'**** &INPUT data has been verified ',/)
 800  FORMAT(//5x,' ***** NTGTS.gt.1 not implemented',5x,' NTGTS(',i2,
     &       ') set to 1'//)
C
      END SUBROUTINE CATCHE
!*==chkwfn.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE CHKWFN(NSYM,NSYM0,NOB,NOB0,SYMTYP,SYMT2,S1,S2,IREFL1,
     &                  IREFL2,IWRITE,NALM)
C***********************************************************************
C
C     CHKWFN - Checks that the data defing the two wavefunctions is
C              compatible. If they are compatible then NALM=0 on
C              exit; otherwise NALM=1
C
C     Assuming that the input wavefunction data does not have a silly
C     mistake in that the orbital sets are different, this routine
C     enforces the following selection rules for non-relativistic
C     property operators:
C
C         (a) The two wavefunctions must have the same spin quant. nos
C         (b) The two wavefunctions, if both type sigma, must be both
C             + or both -.
C
C     Linkage:
C
C         TMTCLOS
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL(KIND=wp), PARAMETER :: EPS=0.0001_wp
C
C Dummy arguments
C
      INTEGER :: IREFL1, IREFL2, IWRITE, NALM, NSYM, NSYM0, SYMT2, 
     &           SYMTYP
      REAL(KIND=wp) :: S1, S2
      INTEGER, DIMENSION(*) :: NOB, NOB0
      INTENT (IN) IREFL1, IREFL2, IWRITE, NOB, NOB0, NSYM, NSYM0, S1, 
     &            S2, SYMT2, SYMTYP
      INTENT (OUT) NALM
C
C Local variables
C
      INTEGER :: I, J
C
C*** End of declarations rewritten by SPAG
C
C---- Default return code for NALM
C
      NALM=0
C
C---- Check that the two orbital sets are identical
C     If they are not then we have a critical condition and the code
C     input contained an error !
C
      DO I=1, MAX(NSYM,NSYM0)
         IF(NOB(I).NE.NOB0(I))THEN
            WRITE(IWRITE,9997)SYMTYP, NSYM, (NOB(J),J=1,20), SYMT2, 
     &                        NSYM0, (NOB0(J),J=1,20)
            CALL TMTCLOS()
         END IF
      END DO
C
C---- Make sure that the S**2 eigenvalues match for the two states
C
      IF(ABS(S1-S2).GT.EPS)THEN
         WRITE(IWRITE,9998)S1, S2
         NALM=1
      END IF
C
C---- For sigma states the +/- reflection symmetry must match too.
C
C     N.B.  The property operators are always of + type in the
C           non-relativistic code.
C
      IF(IREFL1*IREFL2.LT.0)THEN
         WRITE(IWRITE,9999)IREFL1, IREFL2
         NALM=1
      END IF
C
      RETURN
C
C---- Format Statements
C
 9997 FORMAT(/,5X,'Orbitals from configuration data inconsistant:'/,5X,
     &       'SYMTYP ',I2,' NSYM =',I3,', NOB =',20I3,/,5X,'SYMTYP ',I2,
     &       ' NSYM =',I3,', NOB =',20I3)
 9998 FORMAT(/,5X,'Attempt to compute moments between states of ',
     &       ' different spin',/,5X,'S1 =',F3.1,', S2 =',F3.1)
 9999 FORMAT(/,5X,'Attempt to compute moments between Sigma states',
     &       ' with different'/5x,'reflection symmetries',/,5X,
     &       'IREFL1 =',I3,', IREFL2 =',I3)
C
      END SUBROUTINE CHKWFN
!*==cird.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE CIRD(VECI,EIGI,VECJ,EIGJ,IDENT,IWRITE,NOVECI,NOVECJ,
     &                NCSFI,NCSFJ,ISYMTYP,NNUC,GEONUC,CHARG,cname,mgvni,
     &                mgvnj,NFTVI,NVSETI,NFTVJ,NVSETJ)
C***********************************************************************
C
C     CIRD  -  CONTROLS READING OF THE CI DUMPFILES
C
C              PICKS UP THE REQUIRED CI VECTORS AND STORES THEM INTO
C              THE BIG VECTOR. THE VECTORS FOR THE I STATE ARE STORED
C              AT THE BOTTOM OF CORE AND THE VECTORS FOR THE J STATE ARE
C              STORED AT THE NEXT PLACE IN CORE. THE HIGHER AREAS OF
C              CORE ARE USED FOR THE INPUT BUFFERS IN THIS ROUTINE.
C
C     INPUT DATA :
C              NR THE BIG VECTOR
C           IDENT FLAG INDICATING IF ONE OR TWO DIFFERENT SYMMETRIES OF
C                 WAVEFUNCTION ARE BEING USED. ONE WAVEFUNCTION SYMMETRY
C                 MEANS THAT THE VECTORS FROM ONLY ONE CI FILE ARE
C                 REQUIRED. EG. 1ST AND 2ND EIGENVECTORS OF THE
C                 HAMILTONIAN.
C          IWRITE THE LOGICAL UNIT FOR THE PRINTER
C          NOVECI NUMBER OF VECTORS FROM WAVEFUNCTION I
C          NOVECJ NUMBER OF VECTORS FROM WAVEFUNCTION J
C          NCSFI  NUMBER OF CSFS IN WAVEFUNCTION I
C          NCSFJ  NUMBER OF CSFS IN WAVEFUNCTION J
C
C     Linkage:
C
C         READER, WRTVEC
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IDENT, ISYMTYP, IWRITE, MGVNI, MGVNJ, NCSFI, NCSFJ, 
     &           NFTVI, NFTVJ, NNUC, NOVECI, NOVECJ, NVSETI, NVSETJ
      REAL(KIND=wp), DIMENSION(nnuc) :: CHARG
      CHARACTER(LEN=8), DIMENSION(nnuc) :: CNAME
      REAL(KIND=wp), DIMENSION(ncsfi) :: EIGI
      REAL(KIND=wp), DIMENSION(ncsfj) :: EIGJ
      REAL(KIND=wp), DIMENSION(3,nnuc) :: GEONUC
      REAL(KIND=wp), DIMENSION(noveci*ncsfi) :: VECI
      REAL(KIND=wp), DIMENSION(novecj*ncsfj) :: VECJ
      INTENT (IN) IDENT, ISYMTYP, NNUC
      INTENT (INOUT) EIGI, EIGJ
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp) :: E0, S, SZ
      INTEGER :: I, IFAIL, ITIME, NELT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPHZ
      CHARACTER(LEN=120) :: NAME
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(20) :: NHE
C
C*** End of declarations rewritten by SPAG
C
      ALLOCATE(iphz(max(ncsfi,ncsfj)))
C
C---- Call the subroutine to read the CI vector(s) for the first
C     wavefunction.
c
      IF(isymtyp.LT.2)THEN
         CALL READCID(nftvi,nvseti,NAME,NHE,NHD,DTNUC,ncsfi,noveci,
     &                ncsfi,EIGI,VECI,iphz,iwrite)
         e0=dtnuc(1)
         DO i=1, noveci
            eigi(i)=eigi(i)+e0
         END DO
 
         itime=1
      ELSE
         CALL READCIP(nftvi,nvseti,ncsfi,noveci,ncsfi,mgvni,s,sz,nelt,
     &                EIGI,VECI,iphz,cname,geonuc,charg,0,iwrite,1,
     &                ifail)
      END IF
C
C---- If there is a second wavefunction, i.e. not expanded in the same
C     set of CSFs, then read the CI vector(s) for it.
C
      IF(IDENT.EQ.1)THEN
         IF(isymtyp.LT.2)THEN
            CALL READCID(nftvj,nvsetj,NAME,NHE,NHD,DTNUC,ncsfj,novecj,
     &                   ncsfj,EIGj,VECj,iphz,iwrite)
            e0=dtnuc(1)
            DO i=1, novecj
               eigj(i)=eigj(i)+e0
            END DO
         ELSE
            CALL READCIP(nftvj,nvsetj,ncsfj,novecj,ncsfj,mgvnj,s,sz,
     &                   nelt,EIGj,VECj,iphz,cname,geonuc,charg,0,
     &                   iwrite,1,ifail)
         END IF
      END IF
c
      WRITE(iwrite,100)
 100  FORMAT(/5x,'Target vectors read successfully')
c
      DEALLOCATE(iphz)
C
      RETURN
C
      END SUBROUTINE CIRD
!*==cmakstn.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
!***********************************************************************
      FUNCTION CMAKSTN(ISTATE,MGVN,ISPIN,IGU,IREFL,ISYMTYP,NSYM,ksym)
C***********************************************************************
C
C     CMAKSTN - MAKe STate Name label
C
C     Takes the information on the state in binary machine form and
C     generates a character label - This can be used to make the
C     moments file much more accessible to the human eye
C
C     Input data:
C         ISTATE Sequential number assigned to this molecular state
C           MGVN Lamda value for this state or D2h Irr number.
C          ISPIN 2*S+1 for the state S=Total spin quantum number
C            IGU Gerade/Ungerade quantum number for D-inf-h molecules
C          IREFL +/- reflection q. number for Sigma states (Linear)
C        ISYMTYP Defines the molecular point group symmety type
C                =0 C-infinity-V
C                 1 D-infinity-H
C                 2 D2h (Abelian point group or subgroup)
C           NSYM Number of orbital symmetries (IRRs present in set)
C         IWRITE Logical unit for the printer
C
C     Author: Charles J Gillan, May 1994
C
C**********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER(LEN=1), PARAMETER :: CBLANK=' ', CPLUS='+', CMINUS='-', 
     &                           CLPAREN='(', CRPAREN=')'
      CHARACTER(LEN=2), PARAMETER :: CGERAD='-g', CUNGERAD='-u'
      CHARACTER(LEN=10), PARAMETER :: CSTATE='State No. '
C
C Dummy arguments
C
      INTEGER :: IGU, IREFL, ISPIN, ISTATE, ISYMTYP, KSYM, MGVN, NSYM
      CHARACTER(LEN=36) :: CMAKSTN
      INTENT (IN) IGU, IREFL, ISPIN, ISTATE, ISYMTYP, KSYM, MGVN, NSYM
C
C Local variables
C
      CHARACTER(LEN=4), DIMENSION(1) :: CC1
      CHARACTER(LEN=4), DIMENSION(2) :: CC2, CCI, CCS
      CHARACTER(LEN=4), DIMENSION(4) :: CC2H, CC2V, CD2
      CHARACTER(LEN=4), DIMENSION(8) :: CD2H
      CHARACTER(LEN=5), DIMENSION(5) :: CLAMDAS
      CHARACTER(LEN=10), DIMENSION(10) :: CTUPLETS
      INTEGER :: I, IMARK
      INTEGER :: LEN
C
C*** End of declarations rewritten by SPAG
C
C---- Defnitions of the arrays fixed for the duration of the run
C
      DATA(CTUPLETS(I),I=1,10)/'Singlet', 'Doublet   ', 'Triplet   ', 
     &     'Quartet   ', 'Quintuplet', 'Sextuplet ', 'Septuplet ', 
     &     'Octuplet  ', 'Novtuplet ', 'Dectuplet '/
C
      DATA(CLAMDAS(I),I=1,5)/'Sigma', 'Pi   ', 'Delta', 'Phi  ', 
     &     'Gamma'/
C
      DATA cd2h/' AG ', ' B3U', ' B2U', ' B1G', ' B1U', ' B2G', ' B3G', 
     &     ' AU '/
      DATA cc2v/' A1 ', ' B1 ', ' B2 ', ' A2 '/
      DATA ccs/' A'' ', ' A" '/
! 21/01/2008: other point groups added:
      DATA cc2h/' Ag ', ' Au ', ' Bu ', ' Bg '/
      DATA cc2/' A ', ' B '/
      DATA cd2/' A ', ' B3 ', ' B2 ', ' B1 '/
      DATA cci/' Ag ', ' Au '/
      DATA cc1/'A'/
 
C
C---- Initialize the string to all blanks
C
      DO I=1, LEN(CMAKSTN)
         CMAKSTN(I:I)=CBLANK
      END DO
c
C---- Step 1. Build the state number into the string
C
      CMAKSTN(1:10)=CSTATE
      WRITE(CMAKSTN(11:13),1110)ISTATE
C
C.... Step 2. The spin type
C
      IF(ISPIN.LE.10)CMAKSTN(15:24)=CTUPLETS(ISPIN)
C
C
C---- Branch to the section appropriate for the molecular point group
C
      IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
C
C=======================================================================
C
C     Linear Molecules : C-infinity-v and D-infinity-h point groups
C
C=======================================================================
C
C.... Step 3. The lamda value
C
         IF(MGVN.LE.4)CMAKSTN(26:30)=CLAMDAS(MGVN+1)
C
         IMARK=31
C
C..... Step 4. Gerade or Ungerade - if any !!!
C
         IF(IGU.EQ.-1)THEN
            CMAKSTN(31:32)=CUNGERAD
            IMARK=33
         ELSE IF(IGU.EQ.1)THEN
            CMAKSTN(31:32)=CGERAD
            IMARK=33
         END IF
C
C..... Step 5. Reflection symmetry if any !!!!
C
C      This will only be for sigma states
C
         IF(MGVN.EQ.0)THEN
            CMAKSTN(IMARK:IMARK)=CLPAREN
            IMARK=IMARK+1
C
            IF(IREFL.EQ.-1)THEN
               CMAKSTN(IMARK:IMARK)=CMINUS
            ELSE IF(IREFL.EQ.1)THEN
               CMAKSTN(IMARK:IMARK)=CPLUS
            END IF
            IMARK=IMARK+1
C
            CMAKSTN(IMARK:IMARK)=CRPAREN
         END IF
C
      ELSE
C
C=======================================================================
C
C     Non-Linear Molecules : D2h The Abelian point group and subgroups
C
C=======================================================================
C
C.... Step 3. The symmetry
C
! 21st January 2008 Hemal Varambhia: account for distinguishability
         IF(nsym.EQ.4)THEN
            IF(ksym.EQ.1)THEN
               CMAKSTN(26:30)=Cc2v(MGVN+1)
            ELSE IF(ksym.EQ.2)THEN
               CMAKSTN(26:30)=Cc2h(MGVN+1)
            ELSE IF(ksym.EQ.3)THEN
               CMAKSTN(26:30)=Cd2(MGVN+1)
            END IF
         ELSE IF(nsym.EQ.8)THEN
            CMAKSTN(26:30)=Cd2h(MGVN+1)
         ELSE IF(nsym.EQ.2)THEN
            IF(ksym.EQ.1)THEN
               cmakstn(26:30)=Ccs(mgvn+1)
            ELSE IF(ksym.EQ.2)THEN
               cmakstn(26:30)=Cc2(mgvn+1)
            ELSE IF(ksym.EQ.3)THEN
               cmakstn(26:30)=Cci(mgvn+1)
            END IF
         ELSE IF(nsym.EQ.1)THEN
            cmakstn(26:30)=Cc1(mgvn+1)
         END IF
c
      END IF
C
      RETURN
C
C---- Format Statements
C
 1110 FORMAT(I3)
C
      END FUNCTION CMAKSTN
!*==daopen.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DAOPEN(NFDA,LDAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LDAR, NFDA
      INTENT (IN) LDAR, NFDA
C
C Local variables
C
      INTEGER :: LBYTES
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     DAOPEN - Opens the direct access dataset on unit NFDA.
C
C     This unit is used as part of the density matrix formulae sorting
C     process.
C
C     Input data:
C           NFDA  Fortran logical unit number for the dataset
C           LDAR  Number of integer units for each record on the disk
C
C***********************************************************************
C
C---- Convert the number of integer units, LDAR, to bytes.
C     ZM: The original formula LBYTES=4*LDAR has been made platform independent
C         and it uses now the Fortran intrinsic bit_size function to determine
C         the number of bytes in an integer.
C
      LBYTES=(bit_size(LBYTES)/8)*LDAR
C
C---- Issue the open statement now
C
      OPEN(UNIT=NFDA,ACCESS='DIRECT',RECL=LBYTES,STATUS='UNKNOWN')
C
C---- Inform the user about the details of a successful opening
C     and then return to the caller
C
c      WRITE(6,5050) NFDA,LBYTES
C
      RETURN
C
C---- Format statements
C
 5050 FORMAT(/,5X,'The direct access dataset is now open on unit = ',I3,
     &       /,5X,'Its record length = ',I7,' (bytes) ',/)
C
      END SUBROUTINE DAOPEN
!*==denbld.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DENBLD(NFTA,IDENT,LDIAGBUF,LDM,DM,NOCSFI,NOCSFJ,CI,CJ,
     &                  LOFFDBUF,IWRITE,NCD,NPFLG)
C***********************************************************************
C
C     DENBLD - DENsity matrix BuiLDer routine
C
C     Reads the sorted density expression formulae and multiplies them
C     by the corresponding CI eigenvector components finally placing
C     the density elements in their correct place.
C
C     Output data:
C              DM The density matrix
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO, ITWO, ITHREE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IDENT, IWRITE, LDIAGBUF, LDM, LOFFDBUF, NCD, NFTA, 
     &           NOCSFI, NOCSFJ, NPFLG
      REAL(KIND=wp), DIMENSION(NOCSFI) :: CI
      REAL(KIND=wp), DIMENSION(NOCSFJ) :: CJ
      REAL(KIND=wp), DIMENSION(LDM) :: DM
      INTENT (IN) CI, CJ, IDENT, IWRITE, LDIAGBUF, LDM, LOFFDBUF, NCD, 
     &            NFTA, NOCSFI, NOCSFJ, NPFLG
      INTENT (INOUT) DM
C
C Local variables
C
      REAL(KIND=wp) :: CC, CI_I, CI_J
      REAL(KIND=wp), DIMENSION(LOFFDBUF) :: CPQ
      INTEGER :: I, IA, IB, IDBUFNO, IPASS, J, JA, JB, K, L, LBLK, LBOX, 
     &           LI, LJ, MEL, N, NBOX, NCODE, ND, NEL, NI, NL, NRECS, 
     &           H, R
      INTEGER, DIMENSION(LOFFDBUF) :: II, JJ, NPQ
      INTEGER, DIMENSION(LDIAGBUF) :: MP
      LOGICAL :: ZNPGT0
C
C*** End of declarations rewritten by SPAG
C
C---- Decide, once and for all, if the print flag is greater than zero.
C
      ZNPGT0=NPFLG.GT.0
c
C---- Banner header
C
      IF(ZNPGT0)THEN
         WRITE(IWRITE,1000)
         WRITE(IWRITE,1010)NFTA, IDENT, LDM, NOCSFI, NOCSFJ
         WRITE(IWRITE,1020)LDIAGBUF, LOFFDBUF, NCD
      END IF
C
 
!
C---- Initialize the density matrix
C
      DO I=1, LDM
         DM(I)=XZERO
      END DO
C
C---- Position the file of sorted density expression at the start of
C     the formulae. This means that we must skip the header records
C     which define the wavefunction symmetry information.
C
      REWIND NFTA
C
      READ(NFTA)
      IF(IDENT.EQ.1)READ(NFTA)
C
C---- Decide whether or not we have to consider diagonal elements
      IF(IDENT.EQ.1)GO TO 200
C
C=======================================================================
C
C     Process the diagonal elements, if they exist !
C
C=======================================================================
C
C..... Read the header record for the diagonal elements and error check
C      upon it !
C
      READ(NFTA)NCODE, LBOX
C
      IF(NCODE.NE.2 .OR. LBOX.NE.LDIAGBUF)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9950)NCODE, ITWO, LBOX, LDIAGBUF
         CALL TMTCLOS()
      END IF
C
C.... Read the first record of formulae - this will include CSF no 1.
C
C     Remember that IDBUFNO counts the buffer numbers
C
      IDBUFNO=1
C
      READ(NFTA)NEL, (MP(K),K=1,LDIAGBUF)
      MEL=ABS(NEL)
C
C.... Data for CSF 1 begins at element 1 as follows:
C
C          Element 1 - No of entries for CSF 1 formulae
C
C          Elements 2 through NL are the formulae
C
      NL=2*MP(1)+1
C
c        WRITE(IWRITE,1590) IONE,IDBUFNO
c        WRITE(IWRITE,1591) NEL,MEL,IONE,NL
c        WRITE(IWRITE,1595) (MP(N),N=2,NL)
C
C..... Due to the form of the algebra we need only enter the formulae
C      for CSF 1 into the density matrix if the kroenecker delta is 1;
C      that means we are dealing with a state with itself. For two
C      different states, i.e. same CSFs but different CI roots, we do
C      not have to add the CSF 1 formulae.
C
C      In both cases, despite the fact that diagonal formulae are
C      stored as D(ii)-D(11), we may then proceed to add these straight
C      to the density matrix.
C
      IF(NCD.EQ.1)THEN
         DO N=2, NL, 2
            DM(MP(N))=MP(N+1)
         END DO
      END IF
C
C---- If this is diagonal element processing and there is only one
C     CSF then we have finished.
C
      IF(NOCSFI.EQ.1)GO TO 500
C
C..... Process CSFs other than the first one
C
      ND=NL+1
      DO I=2, NOCSFI
         CC=CI(I)*CJ(I)
C
C....... Do we need to read another buffer or not ?
C
         IF(ND.GT.MEL)THEN
            IDBUFNO=IDBUFNO+1
            READ(NFTA)NEL, (MP(K),K=1,LDIAGBUF)
            ND=1
            MEL=ABS(NEL)
         END IF
         NI=ND+1
         NL=ND+2*MP(ND)
C
c           WRITE(IWRITE,1590) I,IDBUFNO
c           WRITE(IWRITE,1591) NEL,MEL,NI,NL
c           WRITE(IWRITE,1595) (MP(N),N=NI,NL)
C
         DO N=NI, NL, 2
            DM(MP(N))=DM(MP(N))+CC*MP(N+1)
         END DO
         ND=NL+1
      END DO
C
C=======================================================================
C
C     Off diagonal element processing
C
C=======================================================================
C
 200  IF(ident.EQ.2)GO TO 500
C
c        WRITE(IWRITE,2010)
C
C---- Now begin a loop over all boxes and passes of formulae
C
      READ(NFTA)NCODE, IPASS, NBOX, LBOX, LBLK
C
c        WRITE(IWRITE,2020) NCODE,IPASS,NBOX,LBOX,LBLK
C
      IF(ncode.EQ.0)GO TO 500
      IF(NCODE.NE.3 .OR. LBOX.NE.LOFFDBUF)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9955)NCODE, ITHREE, LBOX, LOFFDBUF
         CALL TMTCLOS()
      END IF
C
C...... Descend into a loop over boxes of formulae now. Each box
C       consists of a header record followed by records of formulae.
C
      DO L=1, NBOX
         READ(NFTA)IA, IB, JA, JB, NRECS
         LI=IA-1
         LJ=JA-1
C
c           WRITE(IWRITE,2030) L,IA,IB,JA,JB,NRECS,LI,LJ
C
C....... Descend into loop over records of formulae; read each
C        member of the box and proceed to build the density matrix
C        with it.
C
         DO J=1, NRECS
            READ(NFTA)NEL, (II(I),I=1,LBOX), (JJ(I),I=1,LBOX), 
     &                (NPQ(I),I=1,LBOX), (CPQ(I),I=1,LBOX)
            MEL=ABS(NEL)
C
c              WRITE(IWRITE,2110) MEL
c              WRITE(IWRITE,2120) (I,II(I),JJ(I),NPQ(I),CPQ(I),I=1,MEL)
C
C.......... Error check that the values of NPQ are valid. They must be
C           inside the range of the density matrix itself
C
C.......... ZM: Removed Hemal's nonsensical processing of the density matrix in case of equal wavefunctions. This was causing the permanent moments for target models other than HF to be calculated
C               erroneously. After removal of this bug the permanent moments are correct. This has been checked for uracil target CAS model and compared against the Molpro values.
C               I also implemented checking of the NPQ against the DM dimension as mentioned in the comment above.
            DO I=1, MEL
               IF (NPQ(I) > LDM) THEN
                  WRITE(IWRITE,5000) NPQ(I), LDM
                  CALL TMTCLOS()
               ENDIF
               DM(NPQ(I))=DM(NPQ(I))+CI(II(I)+LI)*CJ(JJ(I)+LJ)
     &                    *CPQ(I)
            END DO
c
         END DO
      END DO
C
 
C---- Normal completion : off-diagonal elements constructed
C
      IF(ZNPGT0)WRITE(IWRITE,520)
C
      RETURN
C
C---- Completion with no off-diagonal elements present
C
 500  IF(ZNPGT0)WRITE(IWRITE,550)
C
      RETURN
C
C---- Format Statements
C
 1983 FORMAT('ii(i)=',(i3),' jj(i)=',(i3))
 1000 FORMAT(///,5X,35('-'),//,8X,'Density Matrix Construction',//,5X,
     &       35('-'),//)
 1010 FORMAT(5X,'Sorted formulae are read from unit  = ',I7,/,5X,
     &       'The symmetry switch parameter is    = ',I7,/,5X,
     &       'Size of the density matrix in core  = ',I7,' words',/,5X,
     &       'No. of CSFs in diagonal wavefunct.  = ',I7,/,5X,
     &       'No. of CSFs in second wavefunction  = ',I7,/)
 1020 FORMAT(5X,'Diagonal buffer has maximum size    = ',I7,/,5X,
     &       'Off-diagonal buffers have max. size = ',I7,/,5X,
     &       'Value of the Kroenecker delta       = ',I3,/)
C
 1590 FORMAT(/,5X,'Csf No.',I10,' from buffer number ',I10)
 1591 FORMAT(/,5X,'No. of elements in buffer = ',I10,/,5X,
     &       'No. in absolute form      = ',I10,/,5X,
     &       'CSF  starts/finishes      = ',2I10,//,5X,
     &       'CSF  formulae follow : ',//)
 1595 FORMAT(5X,10(I5,1X,I3),/,(5X,10(I5,1X,I3)))
C
 2010 FORMAT(/,5X,'>>>> Beginning off-diagonal element processing ',/)
 2020 FORMAT(/,5X,'Header record for a new pass has been read ',//,5X,
     &       'Code Word    = ',I5,' Pass No. = ',I5,/,5X,
     &       'No. of boxes = ',I5,' Box size = ',I5,/,5X,
     &       'Blocking fac = ',I5,//)
 2030 FORMAT(//,15X,'Header for box No. ',I5,//,15X,'IA = ',I10,
     &       ' IB = ',I10,/,15X,'JA = ',I10,' JB = ',I10,/,15X,
     &       'No. of records  = ',I10,/,15X,'I bias = ',I10,
     &       ' J bias = ',I10,/)
 2110 FORMAT(/,5X,'Elements constituting the box = ',I10,/,5X,'  No.',
     &       1X,'   II',1X,'   JJ',1X,'  NPQ',1X,'  CPQ ',/,5X,'-----',
     &       1X,'   --',1X,'   --',1X,'  ---',1X,'  --- ',/)
 2120 FORMAT(5X,I5,1X,I5,1X,I5,1X,I5,1X,F15.10)
C
 9900 FORMAT(//,10X,'**** Error in DENBLD: ',//)
 9950 FORMAT(10X,'Header record for diagonal formulae is wrong.',//,10X,
     &       'Expression code read = ',I3,' expected = ',I3,/,10X,
     &       'Buffer size read in  = ',I9,' expected = ',I9,
     &       ' integer words.',/)
 9955 FORMAT(10X,'Header record for off-diagonal formulae is wrong.',//,
     &       10X,'Expression code read = ',I3,' expected = ',I3,/,10X,
     &       'Buffer size read in  = ',I9,' expected = ',I9,
     &       ' integer words.',/)
C
 520  FORMAT(/5X,'Off-Diagonal Elements have been added to the first',
     &       ' order spin-reduced'/5x,'density matrix.',/)
 550  FORMAT(/5X,'There are no Off-Diagonal Elements for this first',
     &       ' order spin-reduced'/5x,'density matrix',/)
 5000 FORMAT(/,5X,'NPQ(I) - density matrix element is wrong; it is ',
     &       'larger than the matrix dimension. NPQ(I) =',I9,' LDM = '
     &       ,I9,/)
C
      END SUBROUTINE DENBLD
!*==denget.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DENGET(LUNIT,INSET,INKEY,FORM)
C***********************************************************************
C
C     DENGET locates set number INSET on unit LUNIT with KEY = INKEY
C
C     If NSET = 0 file is positioned at end-of-information
C             = 1 the file is opened
C             = n file is positioned at the beginning of set number n
C
C     On return INSET = sequence number of current set
C
C***********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER(LEN=11) :: FORM
      INTEGER :: INKEY, INSET, LUNIT
      INTENT (IN) FORM, INKEY, LUNIT
      INTENT (INOUT) INSET
C
C Local variables
C
      INTEGER :: I, KEY, MSET, NREC, NSET
      LOGICAL :: OP
C
C*** End of declarations rewritten by SPAG
C
C---- Enquire about the status of the dataset on unit LUNIT,
C
      INQUIRE(UNIT=LUNIT,OPENED=OP)
C
C---- If the data set is not OPEN then we must open it and rewind.
C
      IF(.NOT.OP)OPEN(UNIT=LUNIT,ERR=99,FORM=FORM,STATUS='UNKNOWN')
      REWIND(UNIT=LUNIT,ERR=100)
C
C---- If we were simply asked to open the dataset, INSET=1, then we
C     may exit now.
C
      IF(INSET.EQ.1)GO TO 800
C
C---- Locate set number INSET
C
 5    CONTINUE
C
      IF(FORM.EQ.'FORMATTED')THEN
         READ(LUNIT,*,END=9)KEY, NSET, NREC
         IF(NSET.EQ.INSET .AND. KEY.EQ.INKEY)THEN
            BACKSPACE LUNIT
            GO TO 800
         ELSE
            DO I=1, NREC
               READ(LUNIT,*,END=199)
            END DO
         END IF
      ELSE
         READ(LUNIT,END=9)KEY, NSET, NREC
         IF(NSET.EQ.INSET .AND. KEY.EQ.INKEY)THEN
            BACKSPACE LUNIT
            GO TO 800
         ELSE
            DO I=1, NREC
               READ(LUNIT,END=199)
            END DO
         END IF
      END IF
C
      IF(NSET+1.EQ.INSET)GO TO 800
C
      MSET=NSET
C
      GO TO 5
C
C---- At END of file on read we branch here.
C
 9    CONTINUE
C
      IF(INSET.EQ.0)THEN
         BACKSPACE LUNIT
         INSET=MSET+1
         GO TO 800
      ELSE
         GO TO 99
      END IF
C
 800  CONTINUE
      RETURN
C
C---- Error handler - set number NSET not found on unit LUNIT
C
 99   CONTINUE
C
      WRITE(6,9900)
      WRITE(6,9920)LUNIT, INSET
C
      STOP
C
C---- Error on rewind of the file
C
 100  CONTINUE
C
      WRITE(6,9900)
      WRITE(6,9930)LUNIT
C
      STOP
C
C---- End of file while reading contents of a dataset
C
 199  CONTINUE
C
      WRITE(6,9900)
      WRITE(6,9940)LUNIT
C
C---- Format Statements
C
 9900 FORMAT(//,5X,'**** Error in DENGET ',//)
 9920 FORMAT(5X,'File, unit =',I3,', does not contain dataset NSET =',
     &       I5,' or data header error')
 9930 FORMAT(5X,'Rewind failed on unit = ',I3,/)
 9940 FORMAT(5X,'Error on read in middle of a set of data ',i5,/)
C
      END SUBROUTINE DENGET
!*==denmak.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DENMAK(IWRITE,NFTA,NPFLG,NFTDEN,IDENSET,NNUC,GEONUC,
     &                  CHARG,cname,NVSETI,NVSETJ,NVECI,NVECJ,NOVECI,
     &                  NOVECJ,NFTVI,NFTVJ,idiag)
C***********************************************************************
C
C     DENMAK - COMPUTES THE TRANSITION DENSITY MATRIX FOR PAIRS OF
C              STATES AND THEN PROCEDES TO EVALUATE THE TRANSITION
C              MOMENTS AND PSEUDO NATURAL ORBITALS FOR THESE PAIRS OF
C              STATES.
C
C     INPUT DATA :
C          IDENT  FLAG TO INDICATE THAT WAVEFUNCTIONS ARE OF THE SAME
C                 SAME/DIFFERENT SYMMETRY (0/1)
C         IWRITE  LOGICAL UNIT FOR THE PRINTER
C           NFTA  LOGICAL UNIT CONTAINING THE SORTED DENSITY MATRIX
C                 EXPRESSIONS
C          NCSFI  NUMBER OF CSFS IN THE FIRST WAVEFUNCTION
C          NCSFJ  NUMBER OF CSFS IN THE SECOND WAVEFUNCTION
C          NPFLG  PRINT FLAGS
C          NFTMT  LOGICAL UNIT FOR OUTPUT OF MOMENT EXPRESSIONS
C         NFTDEN  LOGICAL UNIT FOR OUTPUT OF DENSITY MATRIX. THIS
C                 IS LATER USED BY THE PSN PACKAGE
C
C     OUTPUT DATA :
C                  THERE IS NO OUTPUT DATA AS SUCH. THE CALCULATED
C                  DENSITY MATRICES ARE WRITTEN OUT TO THE PRINTER.
C
C     Linkage:
C             DENBLD
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE global_utils, ONLY : MPROD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: MAXSYM=20
C
C Dummy arguments
C
      INTEGER :: IDENSET, IDIAG, IWRITE, NFTA, NFTDEN, NFTVI, NFTVJ, 
     &           NNUC, NOVECI, NOVECJ, NPFLG, NVECI, NVECJ, NVSETI, 
     &           NVSETJ
      REAL(KIND=wp), DIMENSION(nnuc) :: CHARG
      CHARACTER(LEN=8), DIMENSION(nnuc) :: CNAME
      REAL(KIND=wp), DIMENSION(3,nnuc) :: GEONUC
      INTENT (IN) IDIAG, NVECI, NVECJ
      INTENT (INOUT) IDENSET
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: DEN, VECSTI, VECSTJ
      REAL(KIND=wp), ALLOCATABLE, TARGET, DIMENSION(:) :: ENSTI, ENSTJ
      INTEGER :: I, ICODE, IDENT, IDM, IEBIAS, IERR, IGUI, IGUJ, IREFLI, 
     &           IREFLJ, ISPINI, ISPINJ, ISYMTYP, ISZI, ISZJ, IVCSTR, J, 
     &           JEBIAS, JHIGH, JST, JVCSTR, LBOX, LCOF, LDENMAT, LDII, 
     &           MAXDM, MDEL, MGVNI, MGVNJ, NBLOCKS, NBOX, NCD, NCODE, 
     &           NDCOMP, NDMREC, NELT, NENDI, NENDJ, NOCSFI, NOCSFJ, 
     &           NORB, NSRB, NSYM
      INTEGER :: INT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MBAS, NBKPTL, NBKPTR
      CHARACTER(LEN=120) :: NAMESORT, NAMESORTI, NAMESORTJ
      INTEGER, DIMENSION(MAXSYM) :: NOB
      REAL(KIND=wp), POINTER, DIMENSION(:) :: PENJ
      REAL(KIND=wp) :: PIN, R, S, SZ
C
C*** End of declarations rewritten by SPAG
C
c        WRITE(IWRITE,1000)
c        WRITE(IWRITE,1010) NFTA,NFTDEN,NPFLG
c        WRITE(IWRITE,1020) NNUC
c        DO 1 I=1,NNUC
c           WRITE(IWRITE,1030) I,CHARG(I),
c     *                        GEONUC(1,I),GEONUC(2,I),GEONUC(3,I)
c   1    CONTINUE
C
C---- Prepare sorted moment expressions dataset for reading
C
      REWIND NFTA
C
C---- Obtain preliminary information from the file of sorted moment
C     expressions.
C
C...... (a) There is always one wavefunction and so at least one
C           header record. Note that we compute the g/u and +/-
C           symmetry numbers for C-inf-v/D-inf-h wavefunctions
C
      READ(NFTA)NAMESORTi, MGVNI, S, SZ, R, PIN, NORB, NSRB, NOCSFI, 
     &          NELT, NSYM, IDENT, ISYMTYP, (NOB(I),I=1,NSYM)
C
      ISPINI=2*S+1
      ISZI=INT(SZ)
C
 
      IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
         CALL REFLGU(IGUI,IREFLI,R,PIN)
      ELSE
         IGUI=0
         IREFLI=0
      END IF
C
C.......  Was local array NOB given enough core space ?
C
      IF(NSYM.GT.MAXSYM)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9930)NSYM, MAXSYM
         CALL TMTCLOS()
      END IF
C
C...... (b) If there is a second wavefunction then there are no
C           diagonal formulae and vice versa. Note that again we
C           must compute the g/u and +/- symmetries where needed.
C
C           Additionally we default the value of LDII to unity as it
C           will not be read in.
C
C
      IF(IDENT.EQ.1)THEN
         READ(NFTA)NAMESORTj, MGVNJ, S, SZ, R, PIN, NORB, NSRB, NOCSFJ, 
     &             NELT, NSYM, ISYMTYP
C
         ISPINJ=2*S+1
         ISZJ=INT(SZ)
C
         IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
            CALL REFLGU(IGUJ,IREFLJ,R,PIN)
         ELSE
            IGUj=0
            IREFLj=0
         END IF
C
         LDII=1
      ELSE
         ISPINJ=ISPINI
         ISZJ=ISZI
         MGVNJ=MGVNI
         NOCSFJ=NOCSFI
         IGUJ=IGUI
         IREFLJ=IREFLI
         READ(NFTA)NCODE, LDII
         IF(NCODE.NE.2)THEN
            WRITE(6,9900)
            CALL TMTCLOS()
         END IF
      END IF
C
C...... (c) There may or may not be off-diagonal formulae. If there are
C           then find out how large the buffers here are too.
C
C           This is easy if the wavefunctions involve different sets
C           of CSF expansions. If it is the same set then we must
C           loop over the diagonal expressions first
C
      LBOX=0
C
      IF(IDENT.EQ.1)READ(NFTA)NCODE, NBOX, NBOX, LBOX, LCOF
C
      IF(IDENT.EQ.0 .AND. NOCSFI.GT.1)THEN
 10      READ(NFTA)NCODE
         IF(NCODE.LT.0)THEN
            READ(NFTA)NCODE, NBOX, NBOX, LBOX, LCOF
            IF(NCODE.NE.3)THEN
               WRITE(IWRITE,9900)
               WRITE(IWRITE,9960)NCODE
               CALL TMTCLOS()
            END IF
         ELSE
            GO TO 10
         END IF
      END IF
C
C---- Debug printout of the wavefunctions read in
C
c        WRITE(IWRITE,2000) NAMESORT(1:50)
c        WRITE(IWRITE,2010) IDENT,ISYMTYP,NSYM,MGVNI,MGVNJ
c        WRITE(IWRITE,2020) NELT,S,SZ,NORB,NSRB
c        WRITE(IWRITE,2030) NOCSFI,NOCSFJ
c        WRITE(IWRITE,2040) (NOB(I),I=1,NSYM)
c        WRITE(IWRITE,2050) LDII,NBOX,LBOX,LCOF
C
C---- Establish the symmetry difference between the wavefunctions.
C
C     Remember that MGVN must be augumented by one before using the
C     routine MPROD
C
      IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
         MDEL=ABS(MGVNI-MGVNJ)
      ELSE IF(ISYMTYP.EQ.2)THEN
         MDEL=MPROD(MGVNI+1,MGVNJ+1,0,IWRITE)-1
      ELSE
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9910)ISYMTYP
         CALL TMTCLOS()
      END IF
C
      ALLOCATE(mbas(nsym+1),nbkptl(nsym),nbkptr(nsym))
      CALL SETMBA(MDEL,NSYM,NOB,MBAS,NBKPTL,NBKPTR,ISYMTYP,NBLOCKS,
     &            LDENMAT,0)
      DEALLOCATE(mbas,nbkptl,nbkptr)
C
C---- Data defining where, and how, we should pick up the CI
C     eigen vectors for density matrix construction has already
C     been passed in the argument list. Write a summary of this
C     to the printer
C
c      WRITE(IWRITE,4030) CHARI,NFTVI,CHARI,NVSETI
C
c      IF(IDENT.EQ.1)THEN
c        WRITE(IWRITE,4030) CHARJ,NFTVJ,CHARJ,NVSETJ
c        WRITE(IWRITE,4020)
c      ELSE
c        WRITE(IWRITE,4010)
c        NOVECJ=NOVECI
c      ENDIF
C
C---- We now may compute the number of different density matrices
C     that are to be built from the same set of sorted formulae.
C
C     If the CI vectors are from the same set then we build only
C     a lower half triangle of state pairs for that is the only
C     unique set.
C
      NDCOMP=0
      IF(IDENT.NE.1)THEN
         NDCOMP=NOVECI*(NOVECI+1)/2
      ELSE
         NDCOMP=NOVECI*NOVECJ
      END IF
C
c      WRITE(IWRITE,4012) NDCOMP
C
C..... Compute the number of the last vector needed; this is the
C      clearly given by the following formula
C
C         no of first vector read + number of vectors to read - 1
C
C      This we must do for I and J wavefunctions
C
      NENDI=NVECI+NOVECI-1
C
      IF(IDENT.NE.1)THEN
         NENDJ=NENDI
      ELSE
         NENDJ=NVECJ+NOVECJ-1
      END IF
C
      ALLOCATE(vecsti(noveci*nocsfi),ensti(nendi),stat=ierr)
      IF(ierr.NE.0)THEN
         WRITE(iwrite,910)noveci*nocsfi+nendi
         STOP
      END IF
C
      IF(IDENT.EQ.1)THEN
         ALLOCATE(vecstj(novecj*nocsfj),enstj(nendj),stat=ierr)
         IF(ierr.NE.0)THEN
            WRITE(iwrite,910)novecj*nocsfj+nendj
            STOP
         END IF
      ELSE
         ALLOCATE(vecstj(1),enstj(1))
      END IF
C
C---- Invoke the routine which controls the reading of the eigenvector
C     information
C
      CALL CIRD(vecsti,ensti,vecstj,enstj,IDENT,IWRITE,NOVECI,NOVECJ,
     &          NOCSFI,NOCSFJ,ISYMTYP,NNUC,GEONUC,CHARG,cname,mgvni,
     &          mgvnj,NFTVI,NVSETI,NFTVJ,NVSETJ)
C
C.... In the case that thereare no off diagonal elements it is
C     necessary to fudge the setting of the LBOX so that the
C     dynamic array allocation works on the subroutine call.
C
      IF(LBOX.EQ.0)LBOX=1
C
      IF(ident.EQ.1)THEN
         maxdm=noveci*novecj
      ELSE
         maxdm=noveci*(noveci+1)/2
      END IF
      maxdm=maxdm*ldenmat
      ALLOCATE(den(maxdm))
C
C---- Loop over all unique pairs of vectors and compute the density
C     matrix for each wavefunction set.
C
C     N.B. Variable NCD plays the important role of the Kroenecker delta
C          for wavefunctions of the same symmetry
C
C     We set JHIGH once and for all outside the loop as we know that
C     only the IDENT=0 case must be modified inside each loop.
C
      IDM=1
      IVCSTR=1
      JVCSTR=1
C
 
      IF(IDENT.NE.1)THEN
         DO I=NVECI, NENDI
            DO J=nvecj, I
 
               IF(mdel.EQ.0 .AND. i.EQ.j)THEN
                  ncd=1
               ELSE
                  ncd=0
               END IF
               CALL DENBLD(NFTA,IDENT,LDII,LDENMAT,den(IDM),NOCSFI,
     &                     NOCSFJ,vecsti(ivcstr),vecsti(jvcstr),LBOX,
     &                     IWRITE,NCD,NPFLG)
               JVCSTR=JVCSTR+NOCSFJ
               IF(idiag.EQ.0 .OR. i.EQ.j)IDM=IDM+LDENMAT
            END DO
            JVCSTR=1
            IVCSTR=IVCSTR+NOCSFI
         END DO
      ELSE
         DO I=NVECI, NENDI
            DO J=NVECJ, nendj
 
               IF(mdel.EQ.0 .AND. i.EQ.j)THEN
                  ncd=1
               ELSE
                  ncd=0
               END IF
               CALL DENBLD(NFTA,IDENT,LDII,LDENMAT,den(IDM),NOCSFI,
     &                     NOCSFJ,vecsti(ivcstr),vecstj(jvcstr),LBOX,
     &                     IWRITE,NCD,NPFLG)
               JVCSTR=JVCSTR+NOCSFJ
               IDM=IDM+LDENMAT
            END DO
            JVCSTR=1
            IVCSTR=IVCSTR+NOCSFI
         END DO
      END IF
      ndmrec=idm/ldenmat
C
 
C---- If requested by the user print all of the density matrices
C     that were calculated. In any case we write these to the
C     density matrix library on unit NFTDL starting at IDENSET.
C
      IDM=1
      IEBIAS=1
      DO I=NVECI, NENDI
         IF(IDENT.NE.1)THEN
            JHIGH=I
            namesort=namesorti(6:)
            penj=>ensti
            IF(idiag.EQ.0)THEN
               jst=nvecj
            ELSE
               jst=i
            END IF
         ELSE
            JHIGH=NOVECJ
            namesort=namesorti(6:40)//' / '//namesortj(6:40)
            penj=>enstj
            jst=nvecj
         END IF
         JEBIAS=1
         DO J=jst, JHIGH
            IF(npflg.NE.0)THEN
               WRITE(IWRITE,5025)I, J
               CALL DENPRT(ISYMTYP,MDEL,NSYM,NOB,den(idm),IWRITE)
            END IF
C
            IF(IDENT.NE.1)THEN
               IF(I.EQ.J)THEN
                  ICODE=1
               ELSE
                  ICODE=2
               END IF
            ELSE
               ICODE=3
            END IF
C
            CALL DENWRT(NFTDEN,IDENSET,den(idm),LDENMAT,ndmrec,ICODE,
     &                  NAMESORT,ISPINI,ISZI,MGVNI,IGUI,IREFLI,ISYMTYP,
     &                  ISPINJ,ISZJ,MGVNJ,IGUJ,IREFLJ,ISYMTYP,NSYM,NOB,
     &                  NOCSFI,NOCSFJ,I,J,ENSTi(IEBIAS),penj(JEBIAS),
     &                  NNUC,GEONUC)
            JEBIAS=JEBIAS+1
            IDM=IDM+LDENMAT
         END DO
         IEBIAS=IEBIAS+1
      END DO
      IDENSET=IDENSET+1
c
 
      DEALLOCATE(den)
      DEALLOCATE(vecsti,ensti,vecstj,enstj)
C
      RETURN
C
C---- Format statements
C
 910  FORMAT(//' Insufficient memory, need extra',i10)
 1000 FORMAT(///,5X,'====> DENMAK - BUILD THE DENSITY MATRICES <====',/)
 1010 FORMAT(/,5X,'Logical unit for sorted density expression = ',I3,/,
     &       5X,' " " "   ""   "  density matrix library    = ',I3,/,5X,
     &       'Set number to be used on density library   = ',I9,/,5X,
     &       'Integer words of core available            = ',I9,/,5X,
     &       'Print flag                                 = ',I2,/)
 1020 FORMAT(/,5X,'No. of atomic nuclei = ',I3,//,5X,
     &       'No.  Charge   X,Y,Z Co-ordinates ',/,5X,
     &       '---  ------   ------------------ ')
 1030 FORMAT(5X,I2,3X,F6.3,2X,3(F7.3,1X))
C
 2000 FORMAT(//,5X,'Wavefunction details read from NFTA ',//,5X,
     &       'Wavefunction I name card = ',A,/)
 2010 FORMAT(5X,'IDENT flag has value = ',I3,/,5X,
     &       'Abelian/Cinf-v flag  = ',I3,/,5X,
     &       'No. of orbital symms = ',I3,/,5X,
     &       'IRRed. reps          = ',I3,1X,I3,/)
 2020 FORMAT(5X,'No. of electrons     = ',I3,/,5X,
     &       'Total Spin Quant No. = ',F5.3,/,5X,
     &       'Z proj of Spin       = ',F5.3,/,5X,
     &       'Total no. of orbs    = ',I3,/,5X,
     &       'Total no. of sporbs  = ',I3,/)
 2030 FORMAT(5X,'CSFs per state       = ',I3,1X,I3,/)
 2040 FORMAT(5X,'NOB values = ',/,(5X,20I3,/))
 2050 FORMAT(5X,'Diagonal buffer size   = ',I10,' (integers) '/,5X,
     &       'No. of off-diag boxes  = ',I10,/,5X,
     &       'Size of boxes          = ',I10,' (integers) ',/,5X,
     &       'CI coeffs per box      = ',I10,/)
C
 4010 FORMAT(5X,'The wavefunctions have the same symmetry.')
 4012 FORMAT(5X,'No. of different density matrices to be built=',I3,/)
 4020 FORMAT(5X,'The wavefunctions have different symmetries.')
 4030 FORMAT(5X,'CI wavefunctions read from unit NFTV',A1,' =',I3,
     &       ', with set number NVSET',A1,' =',I3)
 5025 FORMAT(/' Density Matrix between State ',I3,' of manifold I',
     &       ' and State ',I3,' of manifold J',/)
C
 9900 FORMAT(/,10X,'**** Error in DENMAK: ',//)
 9910 FORMAT(10X,'ISymtyp must be 0,1 or 2 but is now = ',I10,/)
 9930 FORMAT(10X,'Local array NOB does not have enough space',//,10X,
     &       'Need = ',I5,' Given = ',I5,' words '/)
 9960 FORMAT(10X,'End of Sorted Diagonal Expressions Reached but ',/,
     &       10X,'next set does not have NCODE=3 rather NCODE=',I10,/)
C
      END SUBROUTINE DENMAK
!*==denprt.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DENPRT(ISYMTYP,MDEL,NSYM,NOB,DENMAT,IWRITE)
C***********************************************************************
C
C     DENPRT - Print out the density matrix in a suitably formatted
C              manner.
C
C     Input data:
C        ISYMTYP  Abelian/C-inf-v flag
C           MDEL  Delta lambda for C-inf-v wavefunctions and the
C                 direct product for Abelian wavefunctions
C           NSYM  Number of symmetries in the orbita set
C            NOB  Number of orbitals per symmetry.
C         DENMAT  The packed density matrix.
C         IWRITE  Logical unit for the printer
C
C     Linkage:
C
C         MPROD, OUTPAK, OUTPUT
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE global_utils, ONLY : MPROD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER(LEN=4), PARAMETER :: CPAIR='Pair'
      CHARACTER(LEN=6), PARAMETER :: CNUMB='Number'
C
C Dummy arguments
C
      INTEGER :: ISYMTYP, IWRITE, MDEL, NSYM
      REAL(KIND=wp), DIMENSION(*) :: DENMAT
      INTEGER, DIMENSION(*) :: NOB
      INTENT (IN) ISYMTYP, MDEL, NOB, NSYM
C
C Local variables
C
      INTEGER :: KDD, L, LC, LD, NA, NB
C
C*** End of declarations rewritten by SPAG
C
C---- Case 1. The density matrix is a set of blocks of lower half
C             triangles. This arises when we have
C
C                 C-inf-v and delta lamda (MDEL) is zero
C
C                 Abelian and MDEL is one
C
C     Case 2. The density matrix is a set of blocks of rectangles.
C
C..... Now branch to the appropriate section of the printing code
C
C
      kdd=1
      IF(MDEL.EQ.0)THEN
         DO L=1, NSYM
            NB=NOB(L)
            IF(NB.GT.0)THEN
               WRITE(IWRITE,520)CNUMB, L
               CALL OUTPAK(DENMAT(KDD),NB,1,IWRITE)
               KDD=KDD+(NB*(NB+1))/2
            END IF
         END DO
      ELSE
         DO L=1, NSYM
            IF(isymtyp.LE.1)THEN
               LC=ABS(L-1+MDEL)+1
            ELSE
               LC=MPROD(l,MDEL+1,0,IWRITE)
            END IF
            IF(LC.LT.L)CYCLE
            IF(LC.GT.NSYM)CYCLE
            NA=NOB(L)
            NB=NOB(LC)
            LD=NA*NB
            IF(LD.NE.0)THEN
               WRITE(IWRITE,520)CPAIR, L, LC
               CALL OUTPUT(DENMAT(KDD),1,Na,1,Nb,Na,Nb,1,IWRITE)
            END IF
            KDD=KDD+LD
         END DO
      END IF
C
      RETURN
C
C---- Format Statements
C
 520  FORMAT(/4X,'Density matrix sub-block for symmetry ',A,1X,I3,1X,I3)
C
      END SUBROUTINE DENPRT
!*==denrd.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DENRD(NFTDEN,IDENGET,DENMAT,LDENMAT,ICODE,CHEAD,ISPINI,
     &                 ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI,ISPINJ,ISZJ,
     &                 LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ,NSYM,NOB,NOCSFI,
     &                 NOCSFJ,NCI,NCJ,ECI,ECJ,NNUC,RGEOM,maxden)
C***********************************************************************
C
C     DENRD - Reads a Density Matrix from the library file on
C             unit NFTDEN
C
C     Input data:
C         NFTDEN Fortran logical unit number for the density matrix
C                library dataset
C        IDENSET Set number, on NFTDEN, at which this data will be
C                written
C         IWRITE Logical unit for the printer
C
C     Output data:
C         DENMAT The complete density matrix to be written out.
C        LDENMAT Total number of elements in the density matrix
C          ICODE Code defining the symmetry make up of the density
C                matrix:
C                = 1 means a state with itself i.e. ground state
C                    hence delta lambda=0 and matrix is triangular
C                = 2 means that states are of the same symmetry but
C                    are not the same state - matrix triangular.
C                = 3 means that states have different symmetries and
C                    so delta lambda is not zero.
C          CHEAD Character header describing the data
C         ISPINI 2*S+1 quantum number for first state (J for second)
C           ISZI Z-projection of S for first state    (J for second)
C         LAMDAI Z-projection of angular momentum (I=first,J=second)
C           IGUI G/U quantum number (if any)      (I=first,J=second)
C         IREFLI Reflection symmetry (if any)     (I=first,J=second)
C       ISYMTYPI Linear/Abelian flag for the molecular point group
C           NSYM Number of symmetries in the orbital set (C-inf-v)
C            NOB No. of orbitals per C-inf-v symmetry
C         NOCSFI No. of CSFs defining wavefunction of first state
C         NOCSFJ No. of CSFs  " " "    " " " " " " "  second state
C            NCI Which CI vector from the Hamiltonian was used for
C                the first state
C            NCJ Which CI vector from the Hamiltonian was used for
C                the second state.
C            ECI Absolute energy, in Hartrees, of the wavefunction I
C            ECJ   " "     " "    "    " " "   "   "    "  "  "    J
C           NNUC Number of nuclei in the system
C          RGEOM Nuclear configuration at which this density matrix was
C                generated: (X,Y,Z) co-ordinates for each nucleus.
C         IWRITE Logical unit for the printer
C
C     Linkage:
C
C         DENGET
C
C     Note:
C
C        The format of each member of the density matrix library is
C
C     Record 1:   Header defining the library member
C
C     Record 2:   Symmetry data for first wavefunction
C
C                 ISPIN,ISZ,LAMDA,IGU,IREFL,ISYMTYPI
C
C     Record 2a:  Symmetry data for the second wavefunction
C
C     Record 3:   Orbital set, CSF and eigenvector information
C
C                 NSYM,NOB,NOCSFI,NOCSFJ
C
C     Record 4:   Geometry information
C
C                 NNUC
C
C                 followed by NNUC records of the form
C
C                 RGEOM(1,I),RGEOM(2,I),RGEOM(3,I)
C
C     Record 5:   Density matrix elements
C
C                 NCI,NCJ,ECI,ECJ,LDENMAT,(DENMAT(I),I=1,LDENMAT)
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: IDENKEY=60
C
C Dummy arguments
C
      CHARACTER(LEN=80) :: CHEAD
      REAL(KIND=wp) :: ECI, ECJ
      INTEGER :: ICODE, IDENGET, IGUI, IGUJ, IREFLI, IREFLJ, ISPINI, 
     &           ISPINJ, ISYMTYPI, ISYMTYPJ, ISZI, ISZJ, LAMDAI, LAMDAJ, 
     &           LDENMAT, MAXDEN, NCI, NCJ, NFTDEN, NNUC, NOCSFI, 
     &           NOCSFJ, NSYM
      REAL(KIND=wp), DIMENSION(maxden) :: DENMAT
      INTEGER, DIMENSION(*) :: NOB
      REAL(KIND=wp), DIMENSION(3,*) :: RGEOM
      INTENT (IN) IDENGET, MAXDEN
      INTENT (OUT) CHEAD, DENMAT, ECI, ECJ, IGUJ, IREFLJ, ISPINJ, 
     &             ISYMTYPJ, ISZJ, LAMDAJ, NCI, NCJ, NOB, NOCSFI, 
     &             NOCSFJ, RGEOM
      INTENT (INOUT) ICODE, IGUI, IREFLI, ISPINI, ISYMTYPI, ISZI, 
     &               LAMDAI, LDENMAT, NNUC, NSYM
C
C Local variables
C
      INTEGER :: I, IRESET, J, NREC
      INTEGER, SAVE :: IDENSET
C
C*** End of declarations rewritten by SPAG
C
C---- Position the unit in order to read the set number IDENGET.
C
      ireset=0
c
 42   IF(idenget.LE.1 .OR. ireset.NE.0)THEN
         IF(idenget.LE.1)idenset=1
c
         CALL DENGET(NFTDEN,IDENSET,IDENKEY,'UNFORMATTED')
C
C---- Record 1:
C
C......... The header record
C
         READ(NFTDEN)I, IDENSET, NREC
         READ(NFTDEN)icode, CHEAD
C
C---- Record 2 and possibly 2a:
C
C......... Data defining the symmetries of the state(s) involved
C
         READ(NFTDEN)ISPINI, ISZI, LAMDAI, IGUI, IREFLI, ISYMTYPI
C
         IF(ICODE.EQ.3)THEN
            READ(NFTDEN)ISPINJ, ISZJ, LAMDAJ, IGUJ, IREFLJ, ISYMTYPJ
         ELSE
            ISPINJ=ISPINI
            ISZJ=ISZI
            LAMDAJ=LAMDAI
            IGUJ=IGUI
            IREFLJ=IREFLI
            ISYMTYPJ=ISYMTYPI
         END IF
C
C---- Record 3:
C
C.......... Data defining the orbital set and the CSF expansions
C           for each of the states
C
         READ(NFTDEN)NSYM, (NOB(I),I=1,NSYM), NOCSFI, NOCSFJ
C
C---- Record 4:
C
C.......... The nuclear geometry information
C
         READ(NFTDEN)NNUC
C
         DO I=1, NNUC
            READ(NFTDEN)(RGEOM(J,I),J=1,3)
         END DO
c
         ireset=0
      END IF
C
C---- Record 5:
C
C.......... The density matrix itself
C
      IF(idenget.EQ.0)RETURN
c
      ldenmat=0
      READ(NFTDEN,ERR=41)icode, NCI, NCJ, ECI, ECJ, LDENMAT, 
     &                   (DENMAT(I),I=1,LDENMAT)
C
C---- Print the data that has just been read.
C
c        WRITE(IWRITE,3000)
c        WRITE(IWRITE,3010) NFTDEN,IDENSET,NREC,ICODE
c        WRITE(IWRITE,3020) CINHEAD(1:40)
c        WRITE(IWRITE,3030) CI,ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI
c        WRITE(IWRITE,3030) CJ,ISPINJ,ISZJ,LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ
c        WRITE(IWRITE,3040) NSYM,(NOB(I),I=1,NSYM)
c        WRITE(IWRITE,3050) NOCSFI,NOCSFJ,NCI,NCJ,ECI,ECJ
c        WRITE(IWRITE,3060) NNUC
c        DO 1 J=1,NNUC
c           WRITE(IWRITE,3070) J,(RGEOM(I,J),I=1,3)
c   1    CONTINUE
C
      RETURN
c
 41   ireset=idenset
      idenset=ireset+1
      BACKSPACE nftden
      GO TO 42
C
C---- Format Statements
C
 3000 FORMAT(///,10X,'Data read from Density Matrix Library',//)
 3010 FORMAT(10X,'Logical unit number for the library   = ',I5,/,10X,
     &       'Density matrix set number to be read  = ',I5,/,10X,
     &       'Number of records in the set          = ',I5,/,10X,
     &       'Format code for type of matrix        = ',I5,//)
 3020 FORMAT(10X,'Set character header (1:40)           = ',A)
 3030 FORMAT(/,10X,'Symmetry details for wavefunction ',A,/,10X,
     &       '-----------------------------------',//,10X,
     &       'Total spin quantum number  = ',I5,/,10X,
     &       'Z-projection of spin       = ',I5,/,10X,
     &       'Lambda or Irreducible Rep. = ',I5,/,10X,
     &       'G/U flag for D-inf-h only  = ',I5,/,10X,
     &       'Sigma reflection (C-inf-v) = ',I5,/,10X,
     &       'C-inf-v or Abelian flag    = ',I5,/)
 3040 FORMAT(/,10X,'Orbital Set details: ',/,10X,
     &       '-------------------- ',//,10X,
     &       'No. of symmetries in the set = ',I5,/,10X,
     &       'Orbital per symmetry = ',20(I3,1X),/)
 3050 FORMAT(/,10X,'CSF details: ',/,10X,'-----------  ',//,10X,
     &       'No. of CSFs in wavefunction I = ',I5,/,10X,
     &       'No. of CSFs in wavefunction J = ',I5,/,10X,
     &       'Root number used for wfn I    = ',I5,/,10X,
     &       ' ""   """    ""   "  wfn J    = ',I5,/,10X,
     &       'Energy (Hartrees)    wfn I    = ',F20.12,/,10X,
     &       ' " "   (Hartrees)    wfn J    = ',F20.12,/)
 3060 FORMAT(/,10X,'Nuclear Configuration: ',/,10X,
     &       '---------------------- ',//,10X,'Number of nuclei = ',I5,
     &       //,10X,' No. ',1X,5X,'X',11X,'Y',11X,'Z',/)
 3070 FORMAT(10X,I5,1X,3(F10.8,1X))
C
      END SUBROUTINE DENRD
!*==denwrt.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
C***********************************************************************
C***********************************************************************
C
C     Density Matrix Library Manipulation Routines
C
C     Functions provided are:
C
C     1)  Read a single density matrix from a library  -  DENRD
C     2)  Write a single density matrix to a library   -  DENWRT
C
C     To locate a dataset on the library file, we use the utility
C     DENGET from the external region code.
C
C***********************************************************************
C***********************************************************************
      SUBROUTINE DENWRT(NFTDEN,IDENSET,DENMAT,LDENMAT,NDM,ICODE,CHEAD,
     &                  ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI,ISPINJ,
     &                  ISZJ,LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ,NSYM,NOB,
     &                  NOCSFI,NOCSFJ,NCI,NCJ,ECI,ECJ,NNUC,RGEOM)
C***********************************************************************
C
C     DENWRT - Writes Density Matrix to the library file on unit NFTDEN
C
C     Input data:
C         NFTDEN Fortran logical unit number for the density matrix
C                library dataset
C        IDENSET Set number, on NFTDEN, at which this data will be
C                written
C         DENMAT The complete density matrix to be written out.
C        LDENMAT Total number of elements in the density matrix
C          ICODE Code defining the symmetry make up of the density
C                matrix:
C                = 1 means a state with itself i.e. ground state
C                    hence delta lambda=0 and matrix is triangular
C                = 2 means that states are of the same symmetry but
C                    are not the same state - matrix triangular.
C                = 3 means that states have different symmetries and
C                    so delta lambda is not zero.
C          CHEAD Character header describing the data
C         ISPINI 2*S+1 quantum number for first state (J for second)
C           ISZI Z-projection of S for first state    (J for second)
C         LAMDAI Z-projection of angular momentum (I=first,J=second)
C           IGUI G/U quantum number (if any)      (I=first,J=second)
C         IREFLI Reflection symmetry (if any)     (I=first,J=second)
C       ISYMTYPI Linear/Abelian flag for the molecular point group
C           NSYM Number of symmetries in the orbital set (C-inf-v)
C            NOB No. of orbitals per C-inf-v symmetry
C         NOCSFI No. of CSFs defining wavefunction of first state
C         NOCSFJ No. of CSFs  " " "    " " " " " " "  second state
C            NCI Which CI vector from the Hamiltonian was used for
C                the first state
C            NCJ Which CI vector from the Hamiltonian was used for
C                the second state.
C            ECI Absolute energy, in Hartrees, of the wavefunction I
C            ECJ   " "     " "    "    " " "   "   "    "  "  "    J
C           NNUC Number of nuclei in the system
C          RGEOM Nuclear configuration at which this density matrix was
C                generated: (X,Y,Z) co-ordinates for each nucleus.
C         IWRITE Logical unit for the printer
C
C     Linkage:
C
C         DENGET
C
C     Note:
C
C        The format of each member of the density matrix library is
C
C     Record 1:   Header defining the library member
C
C     Record 2:   Symmetry data for first wavefunction
C
C                 ISPIN,ISZ,LAMDA,IGU,IREFL,ISYMTYPI
C
C     Record 2a:  Symmetry data for the second wavefunction
C
C     Record 3:   Orbital set, CSF and eigenvector information
C
C                 NSYM,NOB,NOCSFI,NOCSFJ
C
C     Record 4:   Geometry information
C
C                 NNUC
C
C                 followed by NNUC records of the form
C
C                 RGEOM(1,I),RGEOM(2,I),RGEOM(3,I)
C
C     Record 5:   Density matrix elements
C
C                 NCI,NCJ,LDENMAT,(DENMAT(I),I=1,LDENMAT)
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: IDENKEY=60
      CHARACTER(LEN=1), PARAMETER :: CBLANK=' '
C
C Dummy arguments
C
      CHARACTER(LEN=*) :: CHEAD
      REAL(KIND=wp) :: ECI, ECJ
      INTEGER :: ICODE, IDENSET, IGUI, IGUJ, IREFLI, IREFLJ, ISPINI, 
     &           ISPINJ, ISYMTYPI, ISYMTYPJ, ISZI, ISZJ, LAMDAI, LAMDAJ, 
     &           LDENMAT, NCI, NCJ, NDM, NFTDEN, NNUC, NOCSFI, NOCSFJ, 
     &           NSYM
      REAL(KIND=wp), DIMENSION(ldenmat) :: DENMAT
      INTEGER, DIMENSION(nsym) :: NOB
      REAL(KIND=wp), DIMENSION(3,*) :: RGEOM
      INTENT (IN) CHEAD, DENMAT, ECI, ECJ, ICODE, IGUI, IGUJ, IREFLI, 
     &            IREFLJ, ISPINI, ISPINJ, ISYMTYPI, ISYMTYPJ, ISZI, 
     &            ISZJ, LAMDAI, LAMDAJ, LDENMAT, NCI, NCJ, NDM, NNUC, 
     &            NOB, NOCSFI, NOCSFJ, NSYM, RGEOM
C
C Local variables
C
      CHARACTER(LEN=80) :: COUTHEAD
      INTEGER :: I, J, LENGTH, NREC
      INTEGER :: LEN
C
C*** End of declarations rewritten by SPAG
C
      IF(nci.EQ.1 .AND. ncj.EQ.1)THEN
C
c        WRITE(IWRITE,1000)
c        WRITE(IWRITE,1010) NFTDEN,IDENSET,LDENMAT,ICODE
c        WRITE(IWRITE,1020) CHEAD(1:LENGTH)
c        WRITE(IWRITE,1030) 'I',ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI
c        WRITE(IWRITE,1030) 'J',ISPINJ,ISZJ,LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ
c        WRITE(IWRITE,1040) NSYM,(NOB(I),I=1,NSYM)
c        WRITE(IWRITE,1050) NOCSFI,NOCSFJ,NCI,NCJ,ECI,ECJ
c        WRITE(IWRITE,1060) NNUC
c        DO 1 J=1,NNUC
c        WRITE(IWRITE,1070) J,(RGEOM(I,J),I=1,3)
c   1    CONTINUE
C
C---- Compute the number of records to be written on this set o
C
         NREC=ndm+4+NNUC
         IF(ICODE.EQ.3)NREC=NREC+1
C
C---- Copy the input header, CHEAD, to the fixed length (CH*80) header
C     which will be written out to disk. Truncate or pad with blanks
C     as necessary.
C
         LENGTH=LEN(CHEAD)
         IF(LENGTH.GT.80)LENGTH=80
C
         DO I=1, LENGTH
            COUTHEAD(I:I)=CHEAD(I:I)
         END DO
C
         DO I=LENGTH+1, 80
            COUTHEAD(I:I)=CBLANK
         END DO
C
C---- Record 1:
C
C......... The header record
C
         WRITE(NFTDEN)IDENKEY, IDENSET, NREC
         WRITE(NFTDEN)icode, COUTHEAD
C
C---- Record 2 and possibly 2a:
C
C......... Data defining the symmetries of the state(s) involved
C
         WRITE(NFTDEN)ISPINI, ISZI, LAMDAI, IGUI, IREFLI, ISYMTYPI
         IF(ICODE.EQ.3)WRITE(NFTDEN)ISPINJ, ISZJ, LAMDAJ, IGUJ, IREFLJ, 
     &                              ISYMTYPJ
C
C---- Record 3:
C
C.......... Data defining the orbital set and the CSF expansions
C           for each of the states
C
         WRITE(NFTDEN)NSYM, (NOB(I),I=1,NSYM), NOCSFI, NOCSFJ
C
C---- Record 4:
C
C.......... The nuclear geometry information
C
         WRITE(NFTDEN)NNUC
         DO I=1, NNUC
            WRITE(NFTDEN)(RGEOM(J,I),J=1,3)
         END DO
c
      END IF
C
C---- Record 5:
C
C.......... The density matrix itself
C
      WRITE(NFTDEN)icode, NCI, NCJ, ECI, ECJ, LDENMAT, 
     &             (DENMAT(I),I=1,LDENMAT)
C
 
      RETURN
C
C---- Format Statements
C
 1001 FORMAT('denmat(',i3,')= ',d20.12)
 1000 FORMAT(//,10X,'====> DENWRT - WRITE A DENSITY MATRIX <====',//)
 1010 FORMAT(10X,'Logical unit number for the library   = ',I5,/,10X,
     &       'Density matrix set number to be used  = ',I5,/,10X,
     &       'No. of elements in the density matrix = ',I5,/,10X,
     &       'Format code to be used for writing    = ',A,/,10X,
     &       'Density matrix type code              = ',I5)
 1020 FORMAT(10X,'Set character header (1:40)           = ',A)
 1030 FORMAT(/,10X,'Symmetry details for wavefunction ',A,/,10X,
     &       '-----------------------------------',//,10X,
     &       'Total spin quantum number  = ',I5,/,10X,
     &       'Z-projection of spin       = ',I5,/,10X,
     &       'Lambda or Irreducible Rep. = ',I5,/,10X,
     &       'G/U flag for D-inf-h only  = ',I5,/,10X,
     &       'Sigma reflection (C-inf-v) = ',I5,/,10X,
     &       'C-inf-v or Abelian flag    = ',I5,/)
 1040 FORMAT(/,10X,'Orbital Set details: ',/,10X,
     &       '-------------------- ',//,10X,
     &       'No. of symmetries in the set = ',I5,/,10X,
     &       'Orbital per symmetry = ',20(I3,1X),/)
 1050 FORMAT(/,10X,'CSF details: ',/,10X,'-----------  ',//,10X,
     &       'No. of CSFs in wavefunction I = ',I5,/,10X,
     &       'No. of CSFs in wavefunction J = ',I5,/,10X,
     &       'Root number used for wfn I    = ',I5,/,10X,
     &       ' ""   """    ""   "  wfn J    = ',I5,/,10X,
     &       'Energy (Hartrees)    wfn I    = ',F20.12,/,10X,
     &       ' " "   (Hartrees)    wfn J    = ',F20.12,/)
 1060 FORMAT(/,10X,'Nuclear Configuration: ',/,10X,
     &       '---------------------- ',//,10X,'Number of nuclei = ',I5,
     &       //,10X,' No. ',1X,5X,'X',11X,'Y',11X,'Z',/)
 1070 FORMAT(10X,I5,1X,3(F10.5,1X))
C
      END SUBROUTINE DENWRT
!*==denxdr.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DENXDR
C***********************************************************************
C
C     DENXDR - Density Matrix Driver
C
C     Controls data input and the loops over target states in the
C     calculation of target properties (i.e. moments) and transition
C     moments.
C
C     Input data:
C             NR  Workspace array
C           NCOR  Number of real*8 units in the array NR.
C
C     Linkage:
C
C        ERRORM, RDSPED, CHKWFN, MPROD, DMXIIS, DMXIJ,
C        SORTDX, DENMAK, CATCHE, DUMPNML, PRSDX
C        TMDVR1, PRTMOM, POLAR
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO, ITWO, ITHREE
      USE global_utils, ONLY : MPROD, CWBOPN
      USE ukrmol_interface_gbl, ONLY: READ_UKRMOLP_PROPERTY_INTS
      USE mpi_gbl, ONLY: mpi_mod_start, mpi_mod_finalize
      USE class_molecular_properties_data, ONLY :
     &          molecular_properties_data
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: MAXNOB=20, MAXNUC=100, MAXTGT=300, 
     &                      MAXDET=2000, MAXPFG=20
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(MAXNUC) :: CHARG
      CHARACTER(LEN=8), DIMENSION(maxnuc) :: CNAME
      REAL(KIND=wp), ALLOCATABLE, TARGET, DIMENSION(:) :: COEFSI, 
     &            COEFSJ, CPQ
      LOGICAL :: EVALUE, ZLAST, QMOLN
      REAL(KIND=wp) :: GEOCFN, GROUNDEN, PINI, PINJ, 
     &                R1, R2, S1, S1Z, S2, S2Z
      REAL(KIND=wp), DIMENSION(3,MAXNUC) :: GEONUC
                        ! get ground state sym for extracting right polarisability (groundsym)
      INTEGER :: GROUNDSPIN, GROUNDSYM, I, IDENSET, IDENT, IDIAG, 
     &           IDPFLAG, IGUI, IGUJ, IMAPST, IPLOTFG, IPOL, IPROP, 
     &           IQDFG, IREFLI, IREFLJ, ISTATE, ISW, IWRITE, JS1, 
     &           JSTATE, KEYTMT, KK, KSYM, LBOX, LCDOFKI, LCDOFKJ, LCOF, 
     &           LDAR, LENGTH, LNDOFKI, LNDOFKJ, LNGTH, MGVN1, MGVN2, 
     &           NALM, NEED, NELMT, NELT, NFDA, NFTA, NFTD, NFTDL, NFTG, 
     &           NFTI, NFTINT, NFTJ, NFTMOM, NFTMT, NMAPS, NMAPX, NMSET, 
     &           NMSET2, NMSETORI, NNUC, NOCSFKI, NOCSFKJ, NORB, NOVECI, 
     &           NOVECJ, NSRB, NSYM, NSYM0, NTGT, NUCCEN, NUMST, SYMT2, 
     &           SYMTYP
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: ICOEFPI, ICOEFPJ, 
     &       IDETPI, IDETPJ, IDETSI, IDETSJ, NODETI, NODETJ
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDMAP, NPQ
      INTEGER, ALLOCATABLE, DIMENSION(:) :: II, JJ
      CHARACTER(LEN=80) :: NAME
      CHARACTER(LEN=120) :: NAME1, NAME2
      INTEGER, DIMENSION(MAXDET) :: NDTRFI, NDTRFJ
      INTEGER, DIMENSION(MAXTGT) :: NFTSOR, NTGTF, NTGTL, NTGTS
      INTEGER, DIMENSION(MAXNOB) :: NOB, NOB0
      INTEGER, DIMENSION(MAXPFG) :: NPFLG
      LOGICAL :: ukrmolp_ints
      type(molecular_properties_data) :: mol_prop_data
C
C*** End of declarations rewritten by SPAG
C
C---- Namelist definition
C
C           NAME  A title for the run
C           NTGT  Number of different CSF expansions used. If there
C                 is only one eigenvector per expansion then this
C                 is also the number of target states.
C         NFTSOR  Logical units for the wavefunction CSF expansions
C                 (determinant information). These are produced by
C                 the module SORT.
C          NPFLG  Set of print flags for verious stages of the calculati
C         IWRITE  Logical unit for the printer
C          LNGTH  Number of integer words for buffer of diagonal formula
C         LENGTH  Number of words for each component of the records of
C                 off diagonal formulae
C           NFTD  Logical unit for the unsorted density matrix formulae
C           NFTA   " " "   ""   "   "  sorted density matrix formulae
C           NFDA   " " "   ""   "   "  direct access file used by the so
C           LCOF
C           LDAR  Number of integer words per record on unit NFDA
C          IDENT  Switch defining whether or not diagonal formulae have
C                 be computed for a CSF set.
C         NFTINT  Logical unit for the transfromed property interals tha
C                 were generated by TRANS.
C          NFTMT  Logical unit which contains the properties file that
C                 is read by the external region code
C          NMSET  Set number at which the generated moments are to be
C                 output
C          NFTDL  Logical unit for the density matrix library
C         NUCCEN  Number of the nucleus defined as the scattering center
C            ISW  Type of output moments requested
C           IPOL  Switch to carry out polarizibility calculation
C          IPROP  Switch on properties calculation
c          IDIAG  If non-zero, then only diagonal density matrices are
C                 calculated
C         NFTMOM  Logical unit for old style moment output
C         NMSET2  Set number at which to begin output on NFTMOM
C          QMOLN  Logical indicating whether we are running under Quantemol-N
C   UKRMOLP_INTS  Logical indicating whether the property integrals on input
C                 were calculated using UKRmol+.
C
      NAMELIST /INPUT / NAME, NTGT, NFTSOR, NFTG, NTGTF, NTGTS, NTGTL, 
     &   NPFLG, IWRITE, LNGTH, LENGTH, NFTD, NFTA, NFDA, LCOF, NFTMT, 
     &   NMSET, NFTDL, LDAR, NFTINT, iplotfg, iqdfg, ISW, IPOL, NFTMOM, 
     &   NMSET2, iprop, idiag, grounden, ZLAST, geocfn, ksym, QMOLN,
     &   ukrmolp_ints
C
C---- Variable initialization for each pass before reading &INPUT.
C
 2000 NTGT=0
      iprop=1
      idiag=0
      iplotfg=0
      iqdfg=0
      grounden=0.0_wp
      ksym=1
      qmoln=.false.
      groundsym=2
      groundspin=3
      ukrmolp_ints = .true.
C
      DO I=1, MAXTGT
         NFTSOR(I)=0
         NTGTF(I)=i
         NTGTS(I)=1
         NTGTL(I)=1
      END DO
C
C...... Logical unit numbers
C
      IWRITE=6
      NFTINT=17
      NFDA=20
      NFTD=21
      NFTG=25
      NFTA=30
      NFTMT=24
      NFTDL=60
      NFTMOM=61
C
C..... Record sizes for some logical units
C
      LDAR=620
      LNGTH=4095
      LENGTH=819
C
C..... Miscellaneous switches
C
      DO I=1, MAXPFG
         NPFLG(I)=0
      END DO
C
      DO I=1, MAXNOB
         NOB(I)=0
         NOB0(I)=0
      END DO
C
C...... Storage areas used for reference determinants of wavefunctions
C
      DO I=1, MAXDET
         NDTRFI(I)=0
         NDTRFJ(I)=0
      END DO
C
C...... No of CI coefficients per box
C
      LCOF=10000
C
      NMSET=1
      NMSET2=1
      IDENSET=1
      KEYTMT=50
      NUCCEN=-1
      ISW=0
C
C...... Nuclear data
C
      NNUC=-1
C
      DO I=1, MAXNUC
         GEONUC(1,I)=XZERO
         GEONUC(2,I)=XZERO
         GEONUC(3,I)=XZERO
         CHARG(I)=XZERO
      END DO
C
      IPOL=0
C
C..... End of a stack of &INPUT cards
C
      ZLAST=.FALSE.
C
C---- Read the namelist from the input stream
C
      READ(5,INPUT,END=5000,ERR=5500)
C
C---- Write title of this run to printer for user comprehension
C
      WRITE(IWRITE,5)NAME
C
      WRITE(IWRITE,7)NTGT*(NTGT+1)/2
C
C---- Allocate IDMAP for all pairs of target states (triangular)
C
      NMAPX=SUM(NTGTL(1:NTGT)-NTGTS(1:NTGT)+1)
      NEED=NMAPX*(NMAPX+1)/2
      ALLOCATE(IDMAP(3,NEED))
C
C---- DENPROP is no longer SWEDEN-compatible
C
      IF (.NOT.ukrmolp_ints) THEN
         WRITE(IWRITE,8)
         STOP 1
      END IF
C
C---- Calculate polarizabilities from previously constructed moments fil
C     If this branch is taken then no moments are evaluated from this
C     input at all. The code jumps to the polarizibility routine.
C
      IF(IPOL.LT.0)THEN
         WRITE(IWRITE,4)IPOL
         GO TO 4000
      END IF
c
C     IPROP=1 needs all density matrices
      IF(idiag.NE.0)iprop=0
C
C---- For a moments calculation write a summary of the input data that
C     was contained in namelist &INPUT
C
      CALL DUMPNML(NTGT,NFTSOR,NTGTF,NTGTS,NTGTL,NPFLG,IWRITE,NFTD,NFTA,
     &             NFDA,NFTG,NFTINT,NFTMT,NMSET,NUCCEN,ISW,NFTMOM,
     &             NMSET2,ukrmolp_ints)
C
C---- Read-in the UKRMol+ basis sets and property integrals (if needed).
C     All this data are stored in the ukrmolp_interface module. This
C     initialization must be done before the call to GETGEON which is the
C     first routine accessing the UKRMol+ data.
C
      IF (ukrmolp_ints) THEN
         CALL READ_UKRMOLP_PROPERTY_INTS(NFTINT,IWRITE)
      ENDIF
C
C---- Perform error checking on the namelist input data. The run will
C     terminate inside this routine if an error is found.
C
      CALL CATCHE(NTGT,NFTSOR,NTGTF,NTGTS,NTGTL,IWRITE)
C
C---- If the user has not defined NNUC in his namelist input then
C     we must interrogate the file of transformed property integrals
C     in order to obtain values of the variables:
C
C        NNUC, GEONUC and CHARG
C
      IF(NNUC.EQ.-1)CALL GETGEON(NFTINT,IWRITE,NNUC,GEONUC,CHARG,
     &                           npflg(1),ukrmolp_ints)
C
C=======================================================================
C=======================================================================
C
C     B E G I N  L O O P  O V E R  W A V E F U N C T I O N  P A I R S
C
C=======================================================================
C=======================================================================
C
      NMAPS=0
      IDPFLAG=1
      IMAPST=0
      numst=0
      ALLOCATE(ii(length),jj(length),npq(2,length),cpq(length))
      ii = -1; jj = -1; npq = -1; cpq = -1
C
      DO ISTATE=1, NTGT
         numst=numst+NTGTL(ISTATE)-NTGTS(ISTATE)+1
c
         IF(idiag.EQ.0)THEN
            js1=1
         ELSE
            js1=istate
         END IF
C
         DO JSTATE=js1, ISTATE
C
C---- Establish NFTI, NFTJ and IDENT for this pair
C
            NFTI=NFTSOR(ISTATE)
            NFTJ=NFTSOR(JSTATE)
C
            IF(NFTI.EQ.NFTJ)THEN
               IDENT=0
            ELSE
               IDENT=1
            END IF
C
C---- Read in the DETERMINANTAL form for the first and, optionally,
C     the second wavefunction.
C
C     Whenever a second wavefunction is read then some error checking
C     must be performed. It is necessary to make sure that the symmetry
C     properties are the same eg. both C-inf-v type, both with same spin
C     quantum numbers etc.... Additionally the two orbitals sets must
C     be the same !
C
            CALL CWBOPN(NFTD)
            CALL CWBOPN(NFTA)
c
            CALL RDSPED1(NFTI,NOCSFki,NORB,NSRB,NELT,LCDOFki,LNDOFki,
     &                   NSYM,MGVN1,S1,S1Z,R1,PINI,NAME1,SYMTYP)
c
            IF(SYMTYP.EQ.0 .OR. SYMTYP.EQ.1)THEN
               CALL REFLGU(IGUI,IREFLI,R1,PINi)
            ELSE
               IGUI=0
               IREFLI=0
            END IF
c
            ALLOCATE(nodeti(nocsfki),icoefpi(nocsfki),idetpi(nocsfki),
     &               idetsi(lndofki),coefsi(lcdofki))
c
            CALL RDSPED(NFTI,nftd,IWRITE,coefsi,nodeti,icoefpi,idetpi,
     &                  idetsi,NDTRFI,NOCSFki,NOB,LCDOFki,LNDOFki,
     &                  NPFLG(1))
c
            IF(IDENT.EQ.1)THEN
               CALL RDSPED1(NFTJ,NOCSFkj,NORB,NSRB,NELT,LCDOFkj,LNDOFkj,
     &                      NSYM0,MGVN2,S2,S2Z,R2,PINj,NAME2,SYMT2)
c
               IF(SYMT2.EQ.0 .OR. SYMT2.EQ.1)THEN
                  CALL REFLGU(IGUj,IREFLj,R2,PINj)
               ELSE
                  IGUj=0
                  IREFLj=0
               END IF
c
               ALLOCATE(nodetj(nocsfkj),icoefpj(nocsfkj),idetpj(nocsfkj)
     &                  ,idetsj(lndofkj),coefsj(lcdofkj))
c
               CALL RDSPED(NFTj,nftd,IWRITE,coefsj,nodetj,icoefpj,
     &                     idetpj,idetsj,NDTRFj,NOCSFkj,NOB0,LCDOFkj,
     &                     LNDOFkj,NPFLG(1))
c
               CALL CHKWFN(NSYM,NSYM0,NOB,NOB0,SYMTYP,SYMT2,S1,S2,
     &                     IREFLi,IREFLj,IWRITE,NALM)
               IF(NALM.NE.0)GO TO 3011
c        WRITE(IWRITE,112)NOCSFKJ
            ELSE
               NOCSFKJ=NOCSFKI
               MGVN2=MGVN1
               DO KK=1, MAXDET
                  NDTRFJ(KK)=NDTRFI(KK)
               END DO
C
            END IF
C
C---- Tell the user which type of point group this molecule belongs to.
C     Note that we print the D2h group table here as well if the
C     molecule is Abelian. Additionally we may catch an error in SYMTYP.
C
            IF(istate.EQ.1 .AND. jstate.EQ.1)THEN
               IF(SYMTYP.EQ.0 .OR. SYMTYP.EQ.1)THEN
                  WRITE(6,1977)
               ELSE IF(SYMTYP.EQ.2)THEN
                  WRITE(6,1978)NSYM
                  I=MPROD(MGVN1+1,MGVN2+1,1,6)
               ELSE
                  WRITE(6,9910)
                  WRITE(6,9912)SYMTYP
                  CALL TMTCLOS()
               END IF
            END IF
C
            IF(IDENT.NE.1)THEN
               WRITE(NFTD)ITWO, LNGTH
               CALL DMXIIS(idetpi,idetsi,NELT,NDTRFI,nocsfki,NFTD,LNGTH,
     &                     NORB,NSRB,NOB,IWRITE,LNDOFKI,NSYM,NPFLG(2),
     &                     SYMTYP)
            END IF
C
C---- Compute off-diagonal formulae, if required, and write to disk
C
C     The first record, on NFTD, for off-diagonal formulae is
C
C             ITHREE,LENGTH
C
C     where ITHREE is the code which identifies that these are diagonal
C     formulae and LENGTH gives the buffer size in integer words.
C
            IF(IDENT.EQ.1)THEN
               WRITE(NFTD)ITHREE, LENGTH
               CALL DMXIJ(nodeti,icoefpi,idetpi,idetsi,coefsi,nodetj,
     &                    icoefpj,idetpj,idetsj,coefsj,NDTRFI,NDTRFJ,
     &                    LCDOFKI,LCDOFKJ,LNDOFKI,LNDOFKJ,nocsfki,
     &                    nocsfkj,NELT,NOB,NSYM,NSRB,LENGTH,ii,jj,npq,
     &                    cpq,NFTD,IWRITE,IDENT,NPFLG(3),SYMTYP,EVALUE,
     &                    NELMT)
C
            ELSE IF(IDENT.EQ.0 .AND. NOCSFkI.GT.1)THEN
               WRITE(NFTD)ITHREE, LENGTH
               CALL DMXIJ(nodeti,icoefpi,idetpi,idetsi,coefsi,nodeti,
     &                    icoefpi,idetpi,idetsi,coefsi,NDTRFI,NDTRFi,
     &                    LCDOFKI,LCDOFKi,LNDOFKI,LNDOFKi,nocsfki,
     &                    nocsfki,NELT,NOB,NSYM,NSRB,LENGTH,ii,jj,npq,
     &                    cpq,NFTD,IWRITE,IDENT,NPFLG(3),SYMTYP,EVALUE,
     &                    NELMT)
C
            END IF
C
C
C----- END OF FORMULA GENERATION ---------------------------------------
C
C     IF REQUESTED, READ THE DENSITY MATRIX EXPRESSION FILE AND WRITE
C     SUMMARY OF IT TO UNIT 6
C
            IF(NPFLG(4).GE.1)CALL PRNEXP(NFTD,IDENT)
C
            IF(istate.EQ.1 .AND. jstate.EQ.1)THEN
               CALL DAOPEN(NFDA,LDAR)
               LBOX=(LDAR-2)/5
            END IF
C
C----- Sorting of the Density Matrix Expressions  ----------------------
c
            CALL SORTDX(NFTA,NFTD,NFDA,NPFLG(6),LCOF,IDENT,IWRITE,LBOX)
C
C..... Printing of sorted density matrix expressions
C
            IF(NPFLG(6).GE.1)CALL PRSDX(NFTA,IWRITE)
C
C---- Construct the density matrix from the sorted formulae ------------
C
C     First of all we must compute the number of vectors that are
C     required by the user for each wavefunction CSF expansion.
C
            NOVECI=NTGTL(ISTATE)-NTGTS(ISTATE)+1
            NOVECJ=NTGTL(JSTATE)-NTGTS(JSTATE)+1
C
C---- This ISTATE,JSTATE pair will generate a set of density matrix
C     pairs, this set must be flagged uniquely. This is done because
C     we can then easily produce an output file which is compatible with
C     the old code.
C
            IF(ISTATE.EQ.JSTATE)THEN
               NMAPS=(NOVECI*(NOVECI+1))/2
            ELSE
               NMAPS=NOVECI*NOVECJ
            END IF
C
            DO I=1, NMAPS
               IDMAP(1,IMAPST+I)=IDPFLAG
               IDMAP(2,IMAPST+I)=NOVECI
               IDMAP(3,IMAPST+I)=NOVECJ
            END DO
            IMAPST=IMAPST+NMAPS
            IDPFLAG=IDPFLAG+1
C
            CALL DENMAK(IWRITE,NFTA,NPFLG(7),NFTDL,IDENSET,NNUC,GEONUC,
     &                  CHARG,cname,NTGTF(ISTATE),NTGTF(JSTATE),
     &                  NTGTS(ISTATE),NTGTS(JSTATE),NOVECI,NOVECJ,NFTG,
     &                  NFTG,idiag)
c
 3011       DEALLOCATE(nodeti,icoefpi,idetpi,idetsi,coefsi)
            IF(ident.EQ.1)DEALLOCATE(nodetj,icoefpj,idetpj,idetsj,
     &                               coefsj)
C
C=======================================================================
C=======================================================================
C
C     E N D  O F  L O O P  O V E R  W A V E F U N C T I O N  P A I R S
C
C=======================================================================
C=======================================================================
C
         END DO
C
      END DO
C
      DEALLOCATE(ii,jj,npq,cpq)
      IF(iprop.EQ.0)RETURN
C
C---- Compute properties from the evaluated density matrices now
C     and at the same time build the target properties files.
C
      NMSETORI=NMSET2
C
      CALL TMDVR1(NFTDL,NFTINT,NFTMT,IWRITE,npflg(9),ISW,nelt,NUCCEN,
     &            NMSET,IDMAP,IMAPST,NFTMOM,NMSET2,numst,iplotfg,iqdfg,
     &            grounden,ksym,groundsym,groundspin,qmoln,ukrmolp_ints)
C
C---- If requsted by the user we read and print all sets of moments
C     generated in this run - these are the old style format. Here we
C     are reusing some old variables to save space.
C
      IF(NPFLG(10).GT.0)THEN
         DO I=NMSETORI, NMSET2-1
            ISTATE=I
            CALL PRTMOM(NFTMOM,IWRITE,ISTATE,'UNFORMATTED',KEYTMT)
         END DO
      END IF
c
      DEALLOCATE(idmap)
C
C---- Compute dipole polarizabilities if requested. This routine may
C     be invoked directly by after the reading of the namelist
C     or may be called after the run has constructed moments.
C
 4000 CONTINUE
c

      if ( ipol .ne. 0 ) then
         call mol_prop_data % clean
         call mol_prop_data % read_properties(nftmt,nmset)
         call mol_prop_data % calculate_polarisability(iwrite,qmoln)
      end if

c
C----- Alchemy style loop back to check for more input. If the user has
C      set the flag ZLAST to .true. then this means that we
C      automatically have the last namelist and can jump to the end of
C      the routine. There are two advantages to this practice:
C
C        (a) The system does not have to handle an end of file error
C            on read thereby saving a few hundred, maybe, CPU cycles
C
C        (b) Testing of long namelist input decks in step wise
C            refinement fashion is easily faciliated at input level.
C
      IF(ZLAST)THEN
         GO TO 5000
      ELSE
         GO TO 2000
      END IF
c
C----- End of file on input - this is reached when there are no more
C      &INPUTs left to read.
C
 5000 CONTINUE
C
      WRITE(IWRITE,15)
C
      RETURN
C
C---- Error condition handler - Error on read of &INPUT
C
 5500 CONTINUE
C
      WRITE(IWRITE,20)
C
      CALL TMTCLOS()
C
C---- Format Statements
C
 4    FORMAT(/,5x,'IPOL =',I3,' Calculate polarizabilities from',
     &       ' previously generated moments file')
 5    FORMAT(//,5X,A,//)
 6    FORMAT('In denxdr after read ksym= ',(i2))
 7    FORMAT(5X,'No. of different wavefunction pairs to be computed =',
     &       I6,/)
 8    FORMAT(/,5X,'+----------------------------------------------------
     &-------+',/,5X,'| ERROR:   Compatibility with SWEDEN integrals is 
     &no longer |',/,5X,'|          maintained in this version of DENPRO
     &P!!          |',/,5X,'+-------------------------------------------
     &----------------+',//)
 15   FORMAT(/5X,'END OF FILE ON READ. NO MORE INPUT')
 20   FORMAT(5X,'ERROR DURING READ OF &INPUT')
 112  FORMAT(/5x,'No. of CSF in second wavefunction=',I8)
 1977 FORMAT(/,5X,'Molecule is Linear and uses C-infinity-V symmetry',/)
 1978 FORMAT(/,5X,'Molecule is Abelian using D2h subgroup with ',I2,
     &       ' elements ',/)
C
 9910 FORMAT(//,5X,'**** Error in DENXDR ',//)
 9912 FORMAT(5X,'SYMTYP read for molecule not valid = ',I5,/)
C
      END SUBROUTINE DENXDR
!*==dmxii.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DMXII(NELT,NDTRF,NOCSF,NFD,MN,LNGTH,NDI,IWRITE,IDETP,
     &                 NORB,NPFLG,nsrb,lndof)
C***********************************************************************
C
C     DMXII COMPUTES THE DIAGONAL DENSITY MATRIX
C
C     EVALUATE DR
C
C     NC        NUMBER OF OCCUPIED ORBITALS IN REFERENCE DETERMINANT
C     JRON      OCCUPATION NUMBER OF ORBITAL
C     JROB      SEQUENCE NUMBER OF OCCUPIED ORBITAL IN REF. DET.
C
C     THE REFERENCE DETERMINANT IS ASSUMED TO BE ORDERED SO THAT
C     SPIN-ORBITALS CORRESPONDING TO A GIVEN ORBITAL ARE GROUPED
C     TOGETHER
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, LNDOF, LNGTH, NELT, NFD, NOCSF, NORB, NPFLG, 
     &           NSRB
      INTEGER, DIMENSION(nocsf) :: IDETP
      INTEGER, DIMENSION(nsrb) :: MN
      INTEGER, DIMENSION(lndof+1) :: NDI
      INTEGER, DIMENSION(nelt) :: NDTRF
      INTENT (IN) IDETP, LNDOF, NOCSF, NSRB
C
C Local variables
C
      INTEGER :: I, I1, L, M, N, NC, ND, NEL, NELM, NEMB, NN, NROB, NRP
      INTEGER, DIMENSION(norb) :: IROB, IRON, JROB, JRON
      INTEGER, DIMENSION(lngth) :: MP
C
C*** End of declarations rewritten by SPAG
C
C-----CALL SUBROUTINE TO SET UP JRON AND JROB FROM THE REFERENCE
C     DETERMINANT, SIZE IS RETURNED AS NC
C
      CALL INITJR(NELT,NDTRF,MN,NC,JROB,JRON,norb)
C
C-----INITIALIZE VARIABLES FOR THE OUTPUT ROUTINES OF THE
C     DIAGONAL DENSITY MATRIX CALCULATION
C
      NEL=0
      NELM=0
      NEMB=0
      MP=-1
C
C           ********** DENSITY MATRIX CALCULATION **********
C
C-----EVALUATE D11. ON CALLING SETD11 THE SIZE OF JRON AND JROB IS NC,
C     ON RETURN FROM SETD11 THEIR SIZE IS M.
C
      NRP=NDI(1)
      CALL SETD11(NRP,NDI,JROB,JRON,MN,NC,M,norb)
C
C-----STORE AWAY D11 FOR OUTPUT,BEFORE PROCEEDING TO NEXT STAGE OF
C     THE CALCULATION
C
      CALL WRTDII(NEL,M,LNGTH,MP,NFD,NELM,NEMB,JROB,JRON,IWRITE)
C
C-----EVALUATE -D11 AS DR-(DR+D11)
C
      CALL NEGD11(NORB,IRON,IROB,NDI,MN,NROB)
C
C-----EVALUATE DII-D11 FOR THE FIRST DETERMINANTS OF ALL CSFS
C     FROM I=2 TO NOCSF.
C
      DO L=2, NOCSF
C
C--------ESTABLISH THE STARTING POSITION FOR THE DETERMINANT AND
C        ALSO INITIALIZE JROB AND JRON FROM IRON AND IROB
C
         N=IDETP(L)
         NRP=NDI(N)
         I1=N+1
         NN=NROB
         IF(NN.NE.0)THEN
            DO I=1, NN
               JROB(I)=IROB(I)
               JRON(I)=IRON(I)
            END DO
         END IF
C
C--------LOOP OVER ALL REMOVALS OF ONE ELECTRON
C
         ND=-1
         CALL MODRON(JRON,JROB,NN,ND,NDI(I1),NRP,MN,norb)
C
C--------LOOP OVER ALL ADDITIONS OF ONE ELECTRON
C
         I1=I1+NRP
         ND=1
         CALL MODRON(JRON,JROB,NN,ND,NDI(I1),NRP,MN,norb)
C
C--------STORE CALCULATIONS FOR WRITING TO DISC
C
         CALL WRTDII(NEL,NN,LNGTH,MP,NFD,NELM,NEMB,JROB,JRON,IWRITE)
C
      END DO
C
C-----EMPTY CONTENTS OF STORAGE BUFFERS ONTO DISC IF THEY ARE NOT ALREAD
C     WRITTEN OUT
C
      CALL WRTDIX(NEL,LNGTH,MP,NFD,NELM,NEMB,IWRITE,NPFLG)
C
      RETURN
C
      END SUBROUTINE DMXII
!*==dmxiis.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DMXIIS(idetp,idets,NELT,NDTRF,NOCSF,NFD,LNGTH,NORB,
     &                  NSRB,NOB,IWRITE,LNDOF,NSYM,NPFLG,ISYMTYP)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ISYMTYP, IWRITE, LNDOF, LNGTH, NELT, NFD, NOCSF, NORB, 
     &           NPFLG, NSRB, NSYM
      INTEGER, DIMENSION(nocsf) :: IDETP
      INTEGER, DIMENSION(lndof+1) :: IDETS
      INTEGER, DIMENSION(nelt) :: NDTRF
      INTEGER, DIMENSION(*) :: NOB
C
C Local variables
C
      INTEGER :: I
      INTEGER, DIMENSION(nsrb) :: MN
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     DMXIIS - DIAGONAL ELEMENT PREPARATION
C
C     SETS UP THE SPIN ORBITAL TABLE AND ALLOCATES STORAGE IN DYNAMIC
C     ARRAY FOR USE BY ROUTINE DMXII, FOR THE CONSTRUCTION OF DIAGONAL
C     ELEMENTS OF THE DENSITY MATRIX.
C
C***********************************************************************
C
      IF(NPFLG.GT.0)WRITE(IWRITE,500)NELT, (NDTRF(I),I=1,NELT)
C
C     SET UP THE MN TABLE
C
      CALL MNTAB(ISYMTYP,NOB,NSYM,MN,NSRB,NPFLG)
C
C     SET UP THE DIAGONAL ELEMENTS OF THE DENSITY MATRIX
C
      CALL DMXII(NELT,NDTRF,NOCSF,NFD,MN,LNGTH,idets,IWRITE,IDETP,NORB,
     &           NPFLG,nsrb,lndof)
C
      RETURN
C
C---- Format Statements
C
 500  FORMAT(//5X,35('-'),//,8X,'DIAGONAL DENSITY ','CALCULATION'//,5X,
     &       35('-'),//5X,'Number of electrons is ',I3,//,5X,
     &       'Reference ',' Determinant : ',/,(5X,30I3))
C
      END SUBROUTINE DMXIIS
!*==dmxij.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DMXIJ(nodeti,icoefp,idetp,idets,coefsi,nodetj,jcoefp,
     &                 jdetp,jdets,coefsj,NDTRFI,NDTRFJ,LCDOFI,LCDOFJ,
     &                 LNDOFI,LNDOFJ,NOCSFI,NOCSFJ,NELT,NOB,NSYM,NSRB,
     &                 LENGTH,ii,jj,npq,cpq,NFTD,IWRITE,IDENT,NPFLG,
     &                 ISYMTYP,EVALUE,NELMT)
C***********************************************************************
C
C     DMXIJ - COMPUTE THE OFF DIAGONAL ELEMENTS OF THE DENSITY MATRIX
C
C             THIS SUBROUTINE DRIVES THE COMPUTATION OF FORMULAE FOR THE
C             OFF DIAGONAL ELEMENTS OF THE DENSITY MATRIX.
C
C     INPUT DATA :
C         NDTRFI  REFERENCE DETERMINANT OF WAVEFUNCTION I
C         NDTRFJ  REFERENCE DETERMINANT OF WAVEFUNCTION J
C         LCDOFI  NUMBER OF DETERMINANTS IN WAVEFUNCTION I
C         LCDOFJ  NUMBER OF DETERMINANTS IN WAVEFUNCTION J
C         LNDOFI  LENGTH OF THE DETERMINANT CODES FOR WAVEFN I
C         LNDOFJ  LENGTH OF THE DETERMINANT CODES FOR WAVEFN J
C         NOCSFI  NUMBER OF CSFS IN WAVEFUNCTION I
C         NOCSFJ  NUMBER OF CSFS IN WAVEFUNCTION J
C           NELT  NUMBER OF ELECTRONS IN THE WAVEFUNCTION
C            NOB  NUMBER OF ORBITALS IN EACH SYMMETRY
C           NSYM  NUMBER OF SYMMETRIES OVERALL
C           NORB  TOTAL NUMBER OF ORBITALS OVERALL
C           NSRB  TOTAL NUMBER OF SPIN ORBITALS OVERALL
C           NFTD  LOGICAL UNIT FOR OUTPUT OF DENSITY EXPRESSIONS
C         IWRITE  LOGICAL UNIT FOR THE PRINTER
C          NPFLG  PRINT FLAG
C        ISYMTYP  C-INF-V OR ABELIAN POINT GROUP FLAG
C
C    OUTPUT DATA :
C
C           THERE IS NO OUTPUT DATA AS SUCH. ALL FORMULAE ARE WRITTEN
C           ONTO DISC TO BE RE-READ BY OTHER SUBROUTINES.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL :: EVALUE
      INTEGER :: IDENT, ISYMTYP, IWRITE, LCDOFI, LCDOFJ, LENGTH, LNDOFI, 
     &           LNDOFJ, NELMT, NELT, NFTD, NOCSFI, NOCSFJ, NPFLG, NSRB, 
     &           NSYM
      REAL(KIND=wp), DIMENSION(lcdofi) :: COEFSI
      REAL(KIND=wp), DIMENSION(lcdofj) :: COEFSJ
      REAL(KIND=wp), DIMENSION(length) :: CPQ
      INTEGER, DIMENSION(nocsfi) :: ICOEFP, IDETP, NODETI
      INTEGER, DIMENSION(lndofi) :: IDETS
      INTEGER, DIMENSION(length) :: II, JJ
      INTEGER, DIMENSION(nocsfj) :: JCOEFP, JDETP, NODETJ
      INTEGER, DIMENSION(lndofj) :: JDETS
      INTEGER, DIMENSION(nelt) :: NDTRFI, NDTRFJ
      INTEGER, DIMENSION(*) :: NOB
      INTEGER, DIMENSION(2,length) :: NPQ
      INTENT (IN) LCDOFI, LNDOFI
      INTENT (INOUT) IDENT
C
C Local variables
C
      INTEGER :: ICA, ICB, ICSF, IMDISZ, INDA, INDB, JCSF, MAXSO, MDA, 
     &           NBLK, NEL, NODA, NODB
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMDI, JDETPN
      INTEGER, DIMENSION(nsrb) :: MG, MG2, MM, MN, MS
C
C*** End of declarations rewritten by SPAG
C
      IF(NPFLG.GT.0)WRITE(IWRITE,500)
C
C-----INITIALIZE PARAMETERS FOR THE OUTPUT ROUTINES WHICH WRITE FORMULAE
C     ONTO DISC. WRITE HEADER LINE GIVING CODE.
C
      NEL=0
      NELMT=0
      NBLK=0
c
C-----PRINT OUT THE TWO WAVEFUNCTIONS IN THEIR PRESENT FORM. ONLY ONE
C     EXISTS HOWEVER IF WE ARE LOOKING AT THE OFF DIAGONAL ELEMENTS
C     WITHIN A WAVEFUNCTION WITH > 1 CSF.
C
      IF(NPFLG.GT.0)THEN
         CALL PTPWF(IWRITE,NOCSFI,NELT,NDTRFI,NODETI,IDETP,ICOEFP,IDETS,
     &              COEFSi)
         IF(IDENT.EQ.1)CALL PTPWF(IWRITE,NOCSFJ,NELT,NDTRFJ,nodetj,
     &                            jdetp,jcoefp,jdets,coefsj)
      END IF
C
C-----COMPUTE THE SPIN ORBITAL TABLE
c
      CALL MKORBS(ISYMTYP,NOB,NSYM,MN,MG,MM,MS,maxso,NPFLG)
C
C-----IF THIS IS A RUN FOR OFF-DIAGONAL ELEMENTS WITHIN A SINGLE
C     WAVEFUNCTION THEN SET UP POINTERS AND SKIP TO THE FORMULA
C     GENERATION STAGE. NO DETERMINANT REARRANGEMENT IS NEEDED.
C
      IF(IDENT.EQ.1)THEN
C
C-----STORAGE MANAGEMENT - RE-ARRANGE J WAVEFN DETERMINANTS
C
c        allocate (jdetpn(nocsfj),imdi(lndofj+1))
c *** the above was not sufficient in some cases,
c     the line below is purely a guestimate
         imdisz=5*(lndofj+max(lndofj/2,4))+1000
         ALLOCATE(jdetpn(nocsfj),imdi(imdisz))
C
C-----CHANGE THE DETERMINANTS IN WAVEFUNCTION J TO BE DEFINED WITH
C     REPECT TO THE REFERENCE DETERMINANT OF WAVEFUNCTION I.
C
         CALL MODRDA(NELT,NSRB,NDTRFI,NDTRFJ,IMDI,LCDOFJ,MDA,jdets,
     &               coefsj)
         IF(mda.GT.imdisz)THEN
            WRITE(6,*)' Error in DMXIJ, IMDI too small '
            WRITE(6,*)' Need ', mda, ' given ', imdisz
            STOP
         END IF
         CALL SETNDJ(NOCSFJ,IMDI,jdetpn,nodetj,LNDOFj)
C
C-----WRITE OUT THE NEW DETERMINANTS TO THE PRINTER
C
         IF(NPFLG.GT.0)THEN
            WRITE(IWRITE,522)
            CALL PTPWF(IWRITE,NOCSFJ,NELT,NDTRFI,nodetj,jdetpn,jcoefp,
     &                 IMDI,coefsj)
         END IF
      END IF
C
C-----LOOP OVER ALL PAIRS OF CSFS FOR OFF DIAGONAL FORMULAE
C     SET LIMITS ON THE INNER LOOP OVER CSFS ACCORDING TO THE SITUATION
C     EITHER TWO DIFFERENT WAVEFUNCTIONS OF THE SAME WAVEFUNCTION. THE
C     INITIAL CALL TO MAKEMG SETS UP THE POINTER ARRAY IN MG.
C
      CALL MAKEMG(MG2,NSRB,NELT,NDTRFI)
c
      DO ICSF=1, NOCSFI
         IF(IDENT.NE.1)THEN
            DO JCSF=1, NOCSFi
                             !lower limit was icsf+1
C
C--------CONSIDER ALL PAIRS OF DETERMINANTS. PERFORM ORBITAL TEST
C        INITIALLY TO SEE IF THE PAIR CONTRIBUTE TO THE FORMULAE. IF
C        THE PAIR CONTRIBUTE THEN EVALUATE CORRESPONDING COEFFICIENT
C        AND STORE IN BUFFERS FOR WRITING TO DISC.
C
               IF(icsf.NE.jcsf)THEN
                             ! Hemal added this
                  CALL DRYRUN(IDETP,idetp,ICSF,JCSF,EVALUE,MN,IDETS,
     &                        idets,MAXSO)
! Hemal Varambhia DENPROP debug write statement 30th April
!          write(190,*)'subroutine tmtma:'
!          write(190,*)'EVALUE,icsf and jcsf=',EVALUE,icsf,jcsf
!
                  IF(EVALUE)THEN
                     NODA=NODETI(ICSF)
                     NODB=nodeti(JCSF)
                     ICA=ICOEFP(ICSF)
                     ICB=icoefp(JCSF)
                     INDA=IDETP(ICSF)
                     INDB=idetp(JCSF)
C
C-----CALL TMTMA TO EVALUATE THE COEFFICIENTS
C
                     CALL TMTMA(NODA,COEFSI(ica),IDETS(inda),NODB,
     &                          coefsi(ICB),idets(INDB),nsrb,NELT,
     &                          NDTRFi,MN,MS,Mg2,NFTD,NELMT,NBLK,NEL,
     &                          LENGTH,II,JJ,NPQ,CPQ,ICSF,JCSF)
c
                  END IF
               END IF
               ! Hemal added this
            END DO
         ELSE
            DO JCSF=1, NOCSFj
               CALL DRYRUN(IDETP,jdetpn,ICSF,JCSF,EVALUE,MN,IDETS,IMDI,
     &                     MAXSO)
               IF(EVALUE)THEN
                  NODA=NODETI(ICSF)
                  NODB=nodetj(JCSF)
                  ICA=ICOEFP(ICSF)
                  ICB=jcoefp(JCSF)
                  INDA=IDETP(ICSF)
                  INDB=jdetpn(JCSF)
C
C-----CALL TMTMA TO EVALUATE THE COEFFICIENTS
C
                  CALL TMTMA(NODA,COEFSI(ICA),IDETS(INDA),NODB,
     &                       coefsj(ICB),imdi(INDB),nsrb,NELT,NDTRFi,MN,
     &                       MS,Mg2,NFTD,NELMT,NBLK,NEL,LENGTH,II,JJ,
     &                       NPQ,CPQ,ICSF,JCSF)
c
               END IF
            END DO
         END IF
      END DO
C
C-----MAKE FINAL CALL TO EMPTY BUFFERS OF OFF-DIAGONAL FORMULAE.
C
      CALL WRTODX(NEL,LENGTH,II,JJ,NPQ,CPQ,NFTD,NELMT,NBLK,IWRITE,NPFLG)
c
c---- We must flag the case of diagonal density matrix of CI target
c     having no off diagonal elements.  This may happen if the full
c     symmetry of the molecule is not being used, eg homonuclear linear
c     molecule in D2h rather than D-inf-h
c
      IF(ident.EQ.0 .AND. nocsfi.GT.1 .AND. nel.EQ.0)ident=2
C
      IF(ident.EQ.1)DEALLOCATE(imdi,jdetpn)
c
      RETURN
C
 500  FORMAT(//,5X,25('-'),//,6X,'OFF DIAGONAL ELEMENTS',//,5X,25('-'),
     &       //)
 522  FORMAT(//,5X,'THE DETERMINANTS IN THE SECOND WAVEFUNCTION',
     &       ' HAVE BEEN REORGANIZED AND ARE NOW DEFINED W.R.T. THE',/,
     &       5X,'REFERENCE DETERMINANT OF THE FIRST ',
     &       'WAVEFUNCTION. A SUMMARY IS NOW PRINTED.',//)
C
      END SUBROUTINE DMXIJ
!*==dryrun.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DRYRUN(INDI,JNDI,ICSF,JCSF,EVALUE,MN,IDETS,JDETS,MAXSO)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      LOGICAL :: EVALUE
      INTEGER :: ICSF, JCSF, MAXSO
      INTEGER, DIMENSION(*) :: IDETS, INDI, JDETS, JNDI, MN
      INTENT (IN) ICSF, INDI, JCSF, JNDI
      INTENT (OUT) EVALUE
C
C Local variables
C
      INTEGER :: MD, ND, NZ
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     DRYRUN - PERFORMS A 'DRYRUN' EVALUATION
C
C              EACH CSF CONTAINS THE SAME ORBITAL CONFIGURATION, ONLY
C              THE ALLOCATION TO SPIN ORBITALS IS DIFFERENT. THUS EACH
C              CSF PAIR REFERS TO A FIXED COMPARISON OF ORBITALS. THIS
C              SUBROUTINE LOOKS AT THAT ORBITAL COMPARISON AND ALLOCATES
C              THE VALUE .TRUE. TO VARIABLE EVALUE IF THE NUMBER OF
C              ORBITAL DIFFERENCES IS TWO, OTHERWISE EVALUE=.FALSE.
C              FROM THE BASIC PROPERTIES OF THE FIRST ORDER SPIN
C              REDUCED DENSITY MATRIX, IT IS CLEAR WHY THIS CHOICE IS
C              MADE. THIS SUBROUTINE CLOSELY PARALLELS SUBROUTINE DRYRUN
C              IN PROGRAM SPEEDY.
C
C    INPUT DATA :
C          INDI  STARTING POSITIONS FOR DETERMINANTS OF EACH CSF IN
C                ARRAY IDETS FOR WAVEFUNCTION 1.
C          JNDI  STARTING POSITIONS FOR DETERMINANTS OF EACH CSF IN
C                ARRAY JDETS FOR WAVEFUNCTION 2.
C          ICSF  NUMBER OF CSF FOR WAVEFUNCTION 1.
C          JCSF  NUMBER OF CSF FOR WAVEFUNCTION 2.
C            MN  THE ORBITAL TABLE RELATING SPIN ORBITALS TO ORBITALS
C         IDETS  THE DETERMINANTS FOR WAVEFUNCTION 1
C         JDETS  THE DETERMINANTS FOR WAVEFUNCTION 2
C         MAXSO  PARAMETER DETERMINING MAXIMUM NUMBER OF SPIN ORBITALS
C           NBB  WORK SPACE ARRAY
C
C    OUTPUT DATA :
C         EVALUE  LOGICAL FLAG WHICH IS SET DEPENDING ON WHETHER OR NOT
C                 THIS PAIR OF CSFS ACTUALLY CONTRIBUTES.
C
C***********************************************************************
C
      MD=INDI(ICSF)
      ND=JNDI(JCSF)
      CALL ZERO(IDETS(MD),JDETS(ND),NZ,MN,MAXSO)
      IF(NZ.GT.2)THEN
         EVALUE=.FALSE.
      ELSE
         EVALUE=.TRUE.
      END IF
C
      RETURN
C
      END SUBROUTINE DRYRUN
!*==dumpnml.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE DUMPNML(NTGT,NFTSOR,NTGTF,NTGTS,NTGTL,NPFLG,IWRITE,
     &                   NFTD,NFTA,NFDA,NFTG,NFTINT,NFTMT,NMSET,NUCCEN,
     &                   ISW,NFTMOM,NMSET2,ukrmolp_ints)
C***********************************************************************
C
C     DUMPNML - Dumps out a list of the input parameters read via
C               the namelist &INPUT in the module TMTXDR. These
C               variables control the executon of the density matrix
C               and property evaluations. Note that the name of the
C               run, variable NAME, is printed bt TMTXDR and is not
C               included here.
C
C     Input data: (This constitutes &INPUT excluding the NAME variable)
C           NTGT  Number of target states in the problem.
C         NFTSOR  Logical unit numbers for determinant information for
C                 all of the wavefunctions used.
C          NPFLG  Set of print flags for verious stages of the calculati
C         IWRITE  Logical unit for the printer
C          LNGTH  Number of integer words for buffer of diagonal formula
C           NFTD  Logical unit for the unsorted density matrix formulae
C           NFTA   " " "   ""   "   "  sorted density matrix formulae
C           NFDA   " " "   ""   "   "  direct access file used by the so
C           LCOF
C           NFTG  Logical unit for the CI vectors defining the
C                 wavefunctions.
C         NFTINT  Logical unit for the transfromed property interals tha
C                 were generated by TRANS.
C          NFTMT  Logical unit number for the target properties file.
C          NMSET  Set number at which the generated moments are to be
C                 output.
C          NFTDL  Logical unit for the density matrix library
C         NUCCEN  Number of the nucleus defined as the scattering center
C            ISW  Type of output moments requested
C           IPOL  Switch to carry out polarizibility calculation
C         NFTMOM  Logical unit for the 'old style' moments output
C         NMSET2  Output set no. for NFTMOM at begining of output
C   ukrmolp_ints  Are we using the UKRmol+ property integrals?
C
C     Notes:
C
C       No error checking is performed by this routine. Instead the
C       routine CATCHE carries out that task.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ISW, IWRITE, NFDA, NFTA, NFTD, NFTG, NFTINT, NFTMOM, 
     &           NFTMT, NMSET, NMSET2, NTGT, NUCCEN
      INTEGER, DIMENSION(ntgt) :: NFTSOR, NTGTF, NTGTL, NTGTS
      INTEGER, DIMENSION(*) :: NPFLG
      LOGICAL :: ukrmolp_ints
      INTENT (IN) ISW, IWRITE, NFDA, NFTA, NFTD, NFTG, NFTINT, NFTMOM, 
     &            NFTMT, NFTSOR, NMSET, NMSET2, NPFLG, NTGT, NTGTF, 
     &            NTGTL, NTGTS, NUCCEN, ukrmolp_ints
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
C
C---- Source of the property integrals
C
      IF (ukrmolp_ints) THEN
         WRITE(IWRITE,'(/,10X,"Assuming UKRMol+ integrals on input.")')
      ELSE
         WRITE(IWRITE,'(/,10X,"Assuming SWEDEN integrals on input.")')
      ENDIF
C
C---- Wavefunction table defining location of determinants and CI
C     vectors for each wavefunction expansion.
C
      WRITE(IWRITE,9)NTGT
C
      DO I=1, NTGT
         WRITE(IWRITE,10)I, NFTSOR(I), NFTG, NTGTF(I), NTGTS(I), 
     &                   NTGTL(I)
      END DO
C
C---- Print flags defining level of printout
C
      WRITE(IWRITE,14)(NPFLG(I),I=1,11)
C
C---- Fortran units assigned to datasets
C
      WRITE(IWRITE,11)NFTD, NFTA, NFDA, NFTINT, NFTMT, NMSET, NFTMOM, 
     &                NMSET2
C
C---- Print information on the property operators for which moments
C     will be computed. In DENPROP we pick these up from the integrals
C     file so we just notify the user about this.
C
      WRITE(IWRITE,13)
C
C---- Construction of target properties file
C
      WRITE(IWRITE,20)
C
      WRITE(IWRITE,21)NFTMT, NMSET, NUCCEN, ISW
C
      RETURN
C
C---- Format Statements
C
 9    FORMAT(//,5X,'Wavefunction location table',/,5X,
     &       '---------------------------',//,5X,
     &       'Wavefunctions consist of two parts in the Alchemy',/,5X,
     &       'and R-matrix codes. The parts are: ',//,5X,
     &       '  (i)  Slater determinants and Clebsch-Gordan ',/,5X,
     &       '       coefficients i.e. CSF definitions ',/,5X,
     &       ' (ii)  CI expansion coefficient over CSFs     ',//,5X,
     &       'Data corresponding to (i) is in a single dataset ',/,5X,
     &       'for each different wavefunctions while data for (ii)',/,
     &       5X,'is stored, for all wavefunctions, on a single library',
     &       /,5X,'dataset',//,5X,'No. of target states = ',I3,//,5X,
     &       'No.  Determinants   CI vectors  Set. No.  Start/End',/,5X,
     &       '---  ------------   ----------  --------  ---------',/)
 10   FORMAT(5X,I3,2X,I7,5X,3X,I5,5X,2X,I6,2X,I3,1X,I3)
 11   FORMAT(/,5X,'Dataset information and logical unit association:',/,
     &       5X,'-------------------------------------------------',//,
     &       5X,'Un-sorted density matrix expressions (NFTD) = ',I2,/,
     &       5X,'Sorted density matrix expressions    (NFTA) = ',I2,/,
     &       5X,'Direct access file for sorting       (NFDA) = ',I2,/,
     &       5X,'Molecular property integrals       (NFTINT) = ',I2,/,
     &       5X,'Computed moments written to         (NFTMT) = ',I2,/,
     &       5X,'                      with set number NMSET = ',I2,/,
     &       5X,'Old style moments as in TMTJT/CJG  (NFTMOM) = ',I2,/,
     &       5X,'                     with set number NMSET2 = ',I2,/)
 13   FORMAT(/,5X,'Property operator information will be obtained ',/,
     &       5X,'from the file of transformed integrals. To control',/,
     &       5X,'the number evaluated modify your INTS input cards ',/)
 14   FORMAT(/,5X,'The print flags, NPFLG, are as follows: ',/,5X,
     &       '--------------------------------------- ',//,5X,
     &       '1. Wavefunction determinant data from SORT file = ',I2,/,
     &       5X,'2. Diagonal density matrix formulae evaluation  = ',I2,
     &       /,5X,'3. Off-diagonal density matrix form. evaluation = ',
     &       I2,/,5X,
     &       '4. Summary of unsorted density matrix formulae  = ',I2,/,
     &       5X,'5. Printing of transformed property integrals   = ',I2,
     &       /,5X,'6. Sorting procedure for density matrix form.   = ',
     &       I2,/,5X,
     &       '7.                                              = ',I2,/,
     &       5X,'8.                                              = ',I2,
     &       /,5X,'9.                                              = ',
     &       I2,/,5X,
     &       '10. Summary of properties for external region   = ',I2,/,
     &       5X,'11. Polarizibility calculation                  = ',I2,
     &       //)
 20   FORMAT(/,5X,'Target Data File Construction for External Region',/,
     &       5X,'-------------------------------------------------',/)
 21   FORMAT(/5X,'Logical unit for properties file     = ',I3,/,5X,
     &       'Set number to be used for properties = ',I3,/,5X,
     &       'Nuclear center (NUCCEN)  defined     = ',I3,/,5X,
     &       'Flag for add. of nuclear terms (ISW) = ',I3,/)
C
      END SUBROUTINE DUMPNML
!*==getgeon.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE GETGEON(NFT,IWRITE,NNUC,GEONUC,CHARG,ipflg,
     &                   ukrmolp_ints)
C***********************************************************************
C
C     GETGEON - GET GEONuc and other nuclear parameters from the dataset
C               of transformed one electron property integrals
C
C     Input data:
C            NFT Logical unit for the transformed integrals
C         IWRITE Logical unit for the printer
C   ukrmolp_ints Are we using UKRmol+ property integrals
C
C     Output data:
C            NNUC Number of nuclear centers in the problem
C          GEONUC Co-ordnates of each nuclear center
C           CHARG Charge on each nuclear center
C
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      USE GLOBAL_UTILS, ONLY : CWBOPN
      USE ukrmol_interface_gbl, ONLY: GET_GEOM, GET_NAME_SYM
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: MAXBF=255
C
C Dummy arguments
C
      INTEGER :: IPFLG, IWRITE, NFT, NNUC
      REAL(KIND=wp), DIMENSION(*) :: CHARG
      REAL(KIND=wp), DIMENSION(3,*) :: GEONUC
      LOGICAL :: ukrmolp_ints
      INTENT (IN) IPFLG, IWRITE, ukrmolp_ints
      INTENT (INOUT) CHARG, GEONUC, NNUC
C
C Local variables
C
      REAL(KIND=wp) :: EN
      INTEGER :: I, IFOUND, J, LTRB, NBFT, NCODT, NMVL, NLMQ
      INTEGER, DIMENSION(MAXBF) :: MBF, NBF, NPBF
      INTEGER, DIMENSION(8) :: NOB
      CHARACTER(LEN=132) :: NAME
      character(len=8), allocatable, dimension(:) :: CNAME
      real(kind=wp), allocatable, dimension(:) :: xnuc, ynuc, znuc
      real(kind=wp), allocatable, dimension(:) :: charge(:)
C
C*** End of declarations rewritten by SPAG
C
      IF (ukrmolp_ints) THEN
C
C------- Get the data from the objects in ukrmolp_interface module.
C
         write(IWRITE,'(/,10X,"Loading UKRmol+ data...")')
C
         CALL GET_GEOM(NNUC,CNAME,xnuc,ynuc,znuc,charge)
         CHARG(1:NNUC) = charge(1:NNUC)
         GEONUC(1,1:NNUC) = xnuc(1:NNUC)
         GEONUC(2,1:NNUC) = ynuc(1:NNUC)
         GEONUC(3,1:NNUC) = znuc(1:NNUC)
         deallocate(xnuc,ynuc,znuc,charge)
C
         CALL GET_NAME_SYM(NAME,NMVL,NOB,NLMQ)
C
      ELSE
C
C------- Read the header record
C
         write(IWRITE,'(/,10X,"Loading SWEDEN data...")')
C
         CALL CWBOPN(NFT)
         READ(NFT,ERR=40,END=40)NAME, NMVL, NBFT, NNUC, NCODT, LTRB, 
     &                 (NBF(I),I=1,NMVL), (NPBF(I),I=1,NBFT), 
     &                 (MBF(I),I=1,NBFT), (CHARG(I),I=1,NNUC), EN, 
     &                 ((GEONUC(j,I),I=1,NNUC),j=1,3)
         REWIND NFT
      ENDIF
C
C---- Error check the input for sillyness !
C
C...... Number of symmetries must be greater than 0 as must NNUC.
C
      IF(NMVL.LE.0 .OR. NNUC.LE.0)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9920)NMVL, NNUC
         CALL TMTCLOS()
      END IF
C
C---- Write a summary to the printer
C
      IF(ipflg.NE.0)THEN
         WRITE(IWRITE,2000)
         WRITE(IWRITE,2010)NAME(1:132), NMVL, NNUC
      END IF
C
      IF(ipflg.NE.0)WRITE(IWRITE,2020)
      ifound=0
      DO I=1, NNUC
         IF(charg(i).EQ.xzero)ifound=i
         IF(ipflg.NE.0)WRITE(IWRITE,2030)I, GEONUC(1,I), GEONUC(2,I), 
     &                                   GEONUC(3,I), CHARG(I)
      END DO
c
c---- If there is no scattering centre in the input data, invent one
c
      IF(ifound.EQ.0)THEN
         nnuc=nnuc+1
         charg(nnuc)=xzero
         DO i=1, 3
            geonuc(i,nnuc)=xzero
         END DO
      END IF
C
C
      RETURN
C
C---- Error condition handlers
C
 40   CONTINUE
C
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9910)NFT
C
      CALL TMTCLOS()
C
C---- Format Statements
C
 2000 FORMAT(//,10X,'Nuclear Framework Data from Transformed Integrals',
     &       /,10X,'-------------------------------------------------')
 2010 FORMAT(/,10X,'Name on dataset (1:132) = ',A,//,10X,
     &       'No. of symmetries in basis function set = ',I5,/,10X,
     &       'No. of nuclear centers = ',I3,/)
 2020 FORMAT(/,10X,'Data per center follows:',//,10X,' No. ',1X,
     &       '    X    ',1X,'    Y    ',1X,'    Z    ',1X,' Charge ',/,
     &       10X,' --- ',1X,'---------',1X,'---------',1X,'---------',
     &       1X,' ------ ')
 2030 FORMAT(10X,1X,I3,1X,1X,F9.5,1X,F9.5,1X,F9.5,1X,1X,F6.3)
C
 9900 FORMAT(//,10X,'**** Error in GETGEON : ',//)
 9910 FORMAT(10X,'Failed to read header of transformed integrals',/,10X,
     &       'logical unit = ',I3,/)
 9920 FORMAT(10X,'No. of symmetries or No. of nuclei is less than 1',//,
     &       10X,'Values read in = ',I10,1X,I10,' respectively ',/)
C
      END SUBROUTINE GETGEON
!*==getset.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE GETSET(LUNIT,INSET,INKEY,FORM)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER(LEN=11) :: FORM
      INTEGER :: INKEY, INSET, LUNIT
      INTENT (IN) FORM, INKEY, LUNIT
      INTENT (INOUT) INSET
C
C Local variables
C
      CHARACTER(LEN=10) :: CACCESS
      CHARACTER(LEN=11) :: CFORM
      INTEGER :: I, KEY, NREC, NSET
      LOGICAL :: OP, ZFORM
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     GETSET locates set number INSET on unit LUNIT with KEY = INKEY
C
C     If NSET = 0 file is positioned at end-of-information
C             = 1 the file is opened
C             = n file is positioned at the beginning of set number n
C
C     On return INSET = sequence number of current set
C
C     Notes:
C
C       Heavily modified from external region version by Charles J
C     Gillan, March 27th, in order to work properly in DENPROP.
C
C***********************************************************************
C
c        WRITE(6,1000)
c        WRITE(6,1010) LUNIT,INSET,INKEY,FORM
C
C---- Decide once and form all if the dataset is formatted or not.
C
C     If formatted the ZFORM assumes the value .TRUE.
C
      ZFORM=FORM.EQ.'FORMATTED'
C
C---- If necessary then OPEN the file. If already open we check its
C     attributes
C
      INQUIRE(UNIT=LUNIT,OPENED=OP,ACCESS=CACCESS,FORM=CFORM,ERR=10)
C
      IF(OP)THEN
         IF(CACCESS(1:3).EQ.'SEQ' .AND. CFORM(1:4).EQ.FORM(1:4))GO TO 22
         WRITE(6,18)
         WRITE(6,30)LUNIT, CACCESS, CFORM, FORM
         CALL TMTCLOS()
      ELSE
         OPEN(UNIT=LUNIT,ERR=99,FORM=FORM,STATUS='unknown')
         GO TO 22
      END IF
C
 10   CONTINUE
      WRITE(6,18)
      WRITE(6,35)LUNIT
      CALL TMTCLOS()
C
C---- If we are being asked to write the first set then we need only
C     rewind the file.
C
 22   REWIND lunit
      IF(INSET.EQ.1)RETURN
C
C---- Locate set number INSET
C
C     The set number is initialized to zero.
C
C     Successive loops over datasets return to line 5
C
      NSET=0
      NREC=0
C
 5    CONTINUE
C
C...... Read the header cards for this dataset using the appropriate
C       format.
C
      IF(ZFORM)THEN
         READ(LUNIT,*,END=9)KEY, NSET, NREC
      ELSE
         READ(LUNIT,END=9)KEY, NSET, NREC
      END IF
C
C...... If this is a hit then we need do no more work except return.
C
      IF(NSET.EQ.INSET .AND. KEY.EQ.INKEY)THEN
         BACKSPACE LUNIT
         RETURN
      END IF
C
C...... Loop over all elements of this dataset now
C
      IF(ZFORM)THEN
         DO I=1, NREC
            READ(LUNIT,*,END=99)
         END DO
      ELSE
         DO I=1, NREC
            READ(LUNIT,END=99)
         END DO
      END IF
C
C..... Loop back to read another dataset now
C
      GO TO 5
C
C---- The end of file has been reached and we have not found the set
C     number that we require. There are two options:
C
C       INSET = 0 which means that we must find the number of the
C                 next set and return that to the caller. Writing
C                 will commence there.
C
C       INSET > 0 this is valid only if INSET is the number of the
C                 next set that we would wish to write anyway.
C
 9    IF(INSET.EQ.0)THEN
         BACKSPACE LUNIT
         INSET=NSET+1
         RETURN
      ELSE IF(INSET.EQ.NSET+1)THEN
         RETURN
      ELSE
         WRITE(6,18)
         WRITE(6,29)
         WRITE(6,20)LUNIT, INSET
         WRITE(6,21)KEY, NSET, NREC, FORM
         CALL TMTCLOS()
      END IF
C
      RETURN
C
C---- Error handler - End of file during a scan over NREC records.
C
C     This is a serious error.
C
 99   CONTINUE
C
      WRITE(6,18)
      WRITE(6,19)
      WRITE(6,20)LUNIT, INSET
      WRITE(6,21)KEY, NSET, NREC, FORM
C
      CALL TMTCLOS()
C
C---- Format Statements
C
 18   FORMAT(//,10X,'**** Error in GETSET:',//)
 19   FORMAT(10X,' End of file on scan over NREC ',/)
 20   FORMAT(10X,' File on logical unit number   = ',I6,/,10X,
     &       ' does not contain dataset NSET = ',I6,//,10X,
     &       ' Data header error is possible: ',/)
 21   FORMAT(10X,' Last read of header gave ',//,10X,' KEY  = ',I6,
     &       ' NSET = ',I6,/,10X,' NREC = ',I6,/,10X,' FORM = ',A,/)
 29   FORMAT(10X,' Set number cannot be found ',/)
 30   FORMAT(5X,' Parameters for OPEN file on unit LUNIT',i5,' are not',
     &       ' correct. ',//,10X,'Access method: ',A,
     &       ' (obtained)  SEQUENTIAL (needed) ',/,10X,
     &       'Format status: ',A,' (obtained) ',A,' (needed) ',/)
 35   FORMAT(10X,'Inquire on unit = ',I6,' has failed ',/)
 1000 FORMAT(//,20X,'====> GETSET - LOCATE SET NUMBER <====',/)
 1010 FORMAT(/,20X,'Logical unit number      (LUNIT) = ',I6,/,20X,
     &       'Set number to be located (INSET) = ',I6,/,20X,
     &       'Dataset key number         (KEY) = ',I6,/,20X,
     &       'Format of the datset      (FORM) = ',A,/)
C
      END SUBROUTINE GETSET
!*==ihjsr.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      FUNCTION IHJSR(NSIZE,NSEQ,IST,ITARG)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IST, ITARG, NSIZE
      INTEGER :: IHJSR
      INTEGER, DIMENSION(NSIZE+1) :: NSEQ
      INTENT (IN) IST, ITARG, NSIZE
      INTENT (INOUT) NSEQ
C
C Local variables
C
      INTEGER :: ICURSE
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     IHJSR - INTEGER SEARCH FUNCTION
C
C             THIS FUNCTION IMPLEMENTS THE ALGORITHM GIVEN IN THE BOOK
C             'LEARNING TO PROGRAM' BY HOWARD JOHNSTON ON PAGE 389 FOR
C             SEARCHING AN UNORDERED LIST FOR THE OCCURRENCE OF A GIVEN
C             VALUE. THE METHOD OF SETTING A SENTINEL IS USED.
C
C     INPUT DATA :
C           NSIZE  NUMBER OF ELEMENTS IN THE SEQUENCE
C            NSEQ  ARRAY HOLDING THE SEQUENCE WITH ONE EXTRA POSITION AT
C                  THE END FOR HOLDING THE SENTINEL.
C             IST  INITIAL POSITION IN THE SEQUENCE TO START SEARCHING
C                  FROM. THIS IS USUALLY 1.
C           ITARG  INTEGER WHICH IS BEING SOUGHT. IE. THE TARGET.
C
C     OUTPUT DATA :
C           IHJSR  IS THE POSITION IN THE SEQUENCE OF THE REQUIRED
C                  TARGET. IF THE TARGET IS NOT IN THE SEQUENCE, THEN
C                  IHJSR IS GIVEN A VALUE ONE LARGER THAN THE ACTUAL
C                  SIZE OF THE SEQUENCE.
C
C***********************************************************************
C
      IHJSR=0
C
      NSEQ(NSIZE+1)=ITARG
      DO ICURSE=IST, NSIZE+1
         IF(NSEQ(ICURSE).NE.ITARG)CYCLE
         IHJSR=ICURSE
         GO TO 200
      END DO
C
 200  CONTINUE
C
      RETURN
C
      END FUNCTION IHJSR
!*==incent.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE INCENT(VLIST,vsymb,NOCMX,NON,NOADC,rmoi,ipflg)
C***********************************************************************
C
C     INCENT - Identify Nuclear CENTres
C
C     Given the list of nuclear charges this routine builds a table
C     which defines the moecule and then passes this to routine MOI
C     in order to compute the moments of inertia.
C
C     Input data:
C          NOCMX Maximum number of centers in the table allocated
C            NON Number of nuclear centers
C          NOADC Number of aditional centers
C         IWRITE Logical unit for the printer
C
C     Input/Output data:
C       VLIST(NOCMX,8) table containing one row of data per center
C
C     Notes:
C
C       This routine was obtained from the NYU properties code in the
C     form used in the Alchemy II package. It has been recoded in
C     standard Fortran 77 by Charles J Gillan at Queen's Belfast; there
C     is one exception to this however and that is in the use of the
C     EQUIVALENCE statement between a character and a real*8 variable.
C     This has been tested out howeveron IBM, Vax and Cray and it works.
C     The nuclear data table contained in VLIST would be best coded as
C     a structure but that is not available in Fortran 77.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO, HALF=>XHALF
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL(KIND=wp), PARAMETER :: TINY=1.0E-08_wp
      INTEGER, PARAMETER :: IWRITE=6
      CHARACTER(LEN=8), PARAMETER :: CBLANK='        '
C
C Dummy arguments
C
      INTEGER :: IPFLG, NOADC, NOCMX, NON
      REAL(KIND=wp), DIMENSION(3) :: RMOI
      REAL(KIND=wp), DIMENSION(NOCMX,8) :: VLIST
      CHARACTER(LEN=8), DIMENSION(nocmx) :: VSYMB
      INTENT (IN) NOADC, NON
      INTENT (INOUT) VLIST, VSYMB
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(103) :: AMASS, QUADN, SPIN
      CHARACTER(LEN=2), DIMENSION(103) :: ASYMB
      CHARACTER(LEN=8) :: CBLNK8, CSYMBOL
      INTEGER :: I, IATOM, J, NFIRST, NOC
      INTEGER :: INT
C
C*** End of declarations rewritten by SPAG
C
C---- Define the symbols for the elements in the periodic table
C
      DATA(ASYMB(I),I=1,103)/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 
     &     'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 
     &     'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 
     &     'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 
     &     'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 
     &     'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 
     &     'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 
     &     'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 
     &     'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 
     &     'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 
     &     'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'/
C
C---- Define the masses of the elements in the periodic table
C
      DATA(AMASS(I),I=1,92)/1.0078246_wp, 4.002601_wp, 7.01600_wp, 
     &     9.01218_wp, 11.009307_wp, 12.000000_wp, 14.0030738_wp, 
     &     15.9949141_wp, 18.9984022_wp, 19.992441_wp, 22.9898_wp, 
     &     23.98504_wp, 26.98153_wp, 27.976929_wp, 30.973764_wp, 
     &     31.9720727_wp, 34.9688531_wp, 39.962386_wp, 38.96371_wp, 
     &     39.96259_wp, 44.95592_wp, 48._wp, 50.9440_wp, 51.9405_wp, 
     &     54.9380_wp, 55.9349_wp, 58.9332_wp, 57.9353_wp, 62.9296_wp, 
     &     63.9291_wp, 68.9257_wp, 73.9219_wp, 74.9216_wp, 79.9165_wp, 
     &     78.91839_wp, 83.91151_wp, 84.9117_wp, 87.9056_wp, 88.9059_wp, 
     &     89.9043_wp, 92.9060_wp, 97.9055_wp, 98._wp, 101.9037_wp, 
     &102.9048_wp, 107.90389_wp, 106.90509_wp, 113.9036_wp, 114.9041_wp,
     &  120._wp, 120.9038_wp, 129.9067_wp, 126.90466_wp, 131.90416_wp, 
     &  132.9051_wp, 137.9050_wp, 138.9061_wp, 139.9053_wp, 140.9074_wp, 
     &     141.9075_wp, 145._wp, 151.9195_wp, 152.9209_wp, 157.9241_wp, 
     &  159.9250_wp, 163.9288_wp, 164.9303_wp, 165.9304_wp, 168.9344_wp, 
     &  173.9390_wp, 174.9409_wp, 179.9468_wp, 180.9480_wp, 183.9510_wp, 
     &     186.9560_wp, 192._wp, 192.9633_wp, 194.9648_wp, 196.9666_wp, 
     &  201.970625_wp, 204.9745_wp, 207.9766_wp, 208.9804_wp, 209._wp, 
     &  210._wp, 222._wp, 223._wp, 226._wp, 227._wp, 232._wp, 231._wp, 
     &     238._wp/
      DATA(AMASS(I),I=93,103)/237._wp, 244._wp, 243._wp, 247._wp, 
     &    247._wp, 251._wp, 252._wp, 257._wp, 258._wp, 259._wp, 260._wp/
C
C---- Define the nuclear spins and quadrupole moments for the elements
C
      DATA(SPIN(I),I=1,103)/103*0._wp/
      DATA(QUADN(I),I=1,103)/103*0._wp/
C
C---- Initialize CSYMBOL to all blanks
C
      CBLNK8=CBLANK
      CSYMBOL=CBLANK
C
C---- Header for this routine
C
c      WRITE(IWRITE,1010)
C
C---- Print details on the number of centers
C
c      WRITE(IWRITE,40) NON,NOADC
C
C---- Loop over the number of centers and identify each one. This
C     means that we pick up data from the periodic table that is
C     hardwired into the program. For each center a record is built
C     in the array VLIST. The components of VLIST are as follows:
C
C       VLIST(I,1) - X position of the nuclear center
C       VLIST(I,2) - Y position of the nuclear center
C       VLIST(I,3) - Z position of the nuclear center
C       VLIST(I,4) - Charge on this center
C       VLIST(I,5) - Symbol defining this atom
C       VLIST(I,6) - Atomic mass of this nucleus
C       VLIST(I,7) - Nuclear spin
C       VLIST(I,8) - Nuclear quadrupole
C
C     For dummy centers/scattering centers the data is defined as
C     zeros or blanks as appropriate.
C
      NOC=NON+NOADC
C
      DO I=1, NON
         IATOM=INT(VLIST(I,4))
C
C....... If the charge is not zero then we have a real nucleus.
C
         IF(IATOM.GT.0)THEN
            CSYMBOL(7:8)=ASYMB(IATOM)
            Vsymb(I)=cSYMBOL
            VLIST(I,6)=AMASS(IATOM)
            VLIST(I,7)=SPIN(IATOM)
            VLIST(I,8)=QUADN(IATOM)
         ELSE
            Vsymb(I)=cBLNK8
            VLIST(I,6)=XZERO
            VLIST(I,7)=XZERO
            VLIST(I,8)=XZERO
         END IF
      END DO
C
C---- Loop over the additional centers - if indeed any exist !
C
      IF(NOADC.LT.1)GO TO 30
C
      WRITE(6,50)
      WRITE(6,60)
C
      NFIRST=NON+1
      DO I=NFIRST, NOC
         VLIST(I,4)=XZERO
         Vsymb(I)=cBLNK8
         VLIST(I,6)=XZERO
         VLIST(I,7)=XZERO
         VLIST(I,8)=XZERO
      END DO
C
C---- Print a summary table defining each center in the molecule
C
 30   CONTINUE
C
      IF(ipflg.NE.0)THEN
         WRITE(IWRITE,1070)NOC
         WRITE(IWRITE,1075)
         DO I=1, NOC
            WRITE(IWRITE,1080)I, (VLIST(I,J),J=1,4), vsymb(i), 
     &                        (VLIST(I,J),J=6,8)
         END DO
      END IF
C
C---- Compute the moment of inertia tensor
C
      CALL MOI(rmoi,VLIST,NOC,NOCMX,IWRITE,ipflg)
C
      RETURN
C
C---- Format Statements
C
 1010 FORMAT(///,10X,'-----------------------------------------------',
     &       //,11X,'Moments of Inertia and Rotational Frequencies',//,
     &       10X,'-----------------------------------------------',//)
 40   FORMAT(/,10X,'The number of centers = ',I5,
     &       ' Additional Centers = ',I5,/)
 50   FORMAT(//,1X,'ADDITIONAL CENTRES',//)
 60   FORMAT(//,1X,'CENTRE',18X,'CO-ORDINATES'//)
 1070 FORMAT(/,10X,'Details of all ',I5,' centers follows: ',/)
 1075 FORMAT('  No.',2X,'   X    ',1X,'   Y    ',1X,'   Z    ',1X,
     &       ' Charge ',1X,' Symbol ',1X,'  Mass  ',1X,'  Spin  ',1X,
     &       '  Quad  ',/)
 1080 FORMAT(I5,1X,4(F8.3,1X),A8,1X,3(F8.3,1X))
C
      END SUBROUTINE INCENT
!*==index1.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE INDEX1(KT,NIP,NIP2,NOB,NSM,NRI,MBAS,IBKPTL,IBKPTR,
     &                  NBLOCKS,MDEL)
C***********************************************************************
C
C     INDEX1 - FORMS AN INDEX FROM THE ORBITAL LABELS
C
C              EACH OFF-DIAGONAL DENSITY MATRIX FORMULA ELEMENT
C              CONSISTS OF A PAIR OF ORBITALS. THIS ROUTINE CONVERTS
C              THIS PAIR INTO A SINGLE INDEX GIVING THE POSITION IN THE
C              DENSITY MATRIX ARRAY. NOTICE THAT THE INDEX IS COMPUTED
C              IN ORDER TO TIE IN WITH THE WAY THAT THE INTEGRALS ARE
C              PACKED. IE. ROW WISE - NOT COLUMN WISE !
C              THE CORRECT CHOICE OF INDEX IS CRUCIAL AND IS DISCUSSED
C              FURTHER IN THE DOCUMENTATION.
C
C     INPUT DATA :
C             KT  NUMBER OF ELEMENTS TO BE PROCESSED IN THE CURRENT CALL
C           NIP2  PAIR OF ORBITALS TO BE CONVERTED
C            NOB  NUMER OF ORBITALS PER C-INF-V SYMMETRY
C            NSM  LAMBDA VALUES OF ALL ORBITALS
C            NRI  SEQUENCE NUMBERS OF ORBITALS SYMMETRY BY SYMMETRY
C           MBAS  POINTER ARRAY TO START OF STORAGE LOCATIONS FOR EACH
C                 OVERLAP IN DENSITY MATRIX.
C
C     OUTPUT DATA :
C             NIP  SINGLE INDEX FOR EACH ORBITAL PAIR.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: KT, MDEL, NBLOCKS
      INTEGER, DIMENSION(nblocks) :: IBKPTL, IBKPTR
      INTEGER, DIMENSION(*) :: MBAS, NOB, NRI, NSM
      INTEGER, DIMENSION(kt) :: NIP
      INTEGER, DIMENSION(2,kt) :: NIP2
      INTENT (IN) IBKPTL, IBKPTR, KT, MBAS, MDEL, NBLOCKS, NIP2, NOB, 
     &            NRI, NSM
      INTENT (OUT) NIP
C
C Local variables
C
      INTEGER :: IA, IB, IDX, KK, MA, MB, NK, NP, NQ
C
C*** End of declarations rewritten by SPAG
C
C---- LOOP OVER THE ORBITAL PAIRS.
C
      DO NK=1, KT
         IA=NIP2(1,NK)
         IB=NIP2(2,NK)
         MA=NSM(IA)
         MB=NSM(IB)
         NP=NRI(IA)
         NQ=NRI(IB)
 
         IF(MDEL.LE.1)THEN
            IF(MA.EQ.MB)THEN
               IDX=(NP*(NP-1))/2+NQ+MBAS(MA+1)
            ELSE IF(MA.LT.MB)THEN
               IDX=(NQ-1)*NOB(MA+1)+NP+MBAS(MA+1)
            ELSE
               IDX=(NP-1)*NOB(MB+1)+NQ+MBAS(MB+1)
            END IF
            NIP(NK)=IDX
         ELSE
            DO KK=1, NBLOCKS
               IF((MA+1.EQ.IBKPTL(KK)) .AND. (MB+1.EQ.IBKPTR(KK)))THEN
                  IF(MA.EQ.MB)THEN
                     IDX=(NP*(NP-1))/2+NQ+MBAS(KK)
                  ELSE IF(MA.LT.MB)THEN
                     IDX=(NQ-1)*NOB(MA+1)+NP+MBAS(KK)
                  ELSE
                     IDX=(NP-1)*NOB(MB+1)+NQ+MBAS(KK)
                  END IF
                  NIP(NK)=IDX
               END IF
            END DO
         END IF
      END DO
C
      RETURN
C
 100  FORMAT(5x,i2,', nip=',i2,',ma=',i2,',mb=',i2,
     &       ',nip2(1,h)=',(i2),', nip2(2,h)=',(i2))
      END SUBROUTINE INDEX1
!*==indexx.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE INDEXX(N,ARRIN,INDX)
C***********************************************************************
C
C     INDEXX - Takes an array and produces a set of indices such
C              that ARRIN(INDX(J)) is in ascending order J=1,2,..,N
C
C     Input data:
C              N number of elements in the array to be ordered
C          ARRIN R*8 array which is to placed in ascending order.
C
C     Output data:
C            INDX a set of indices for ascending indices
C
C     Notes:
C
C     This routine is taken from the book Numerical Receipes by
C     Press, Flannery, Teukolsky and Vetterling Chapter 8 p. 233.
C     ISBN 0-521-30811-9 pub. Cambridge University Press (1986)
C     QA297.N866
C
C     This routines has been adapted by Charles J Gillan for use
C     in the R-matrix codes.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL(KIND=wp), PARAMETER :: VSMALL=1.0E-20_wp
C
C Dummy arguments
C
      INTEGER :: N
      REAL(KIND=wp), DIMENSION(n) :: ARRIN
      INTEGER, DIMENSION(n) :: INDX
      INTENT (IN) ARRIN, N
      INTENT (INOUT) INDX
C
C Local variables
C
      INTEGER :: I, INDXT, IR, J, L
      REAL(KIND=wp) :: Q
C
C*** End of declarations rewritten by SPAG
C
C
C
C
C
C---- Initialize the index array with consecutive integers
C
      DO J=1, N
         INDX(J)=J
      END DO
C
      L=N/2+1
      IR=N
C
C---- From here on the algorithm is HEAPSORT wit indirect addressing
C     through INDX in all references to ARRIN
C
 10   CONTINUE
C
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
C
      I=L
      J=L+L
C
 20   CONTINUE
C
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
C
C---- Loop back to process another element
C
      GO TO 10
C
 800  CONTINUE
C
      RETURN
C
      END SUBROUTINE INDEXX
!*==iniptab.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE INIPTAB(ISTATE,JSTATE,IPCDE,PROPS,MAXPROPS)
C***********************************************************************
C
C     INIPTAB - INItializes the Properties TABle
C
C     Input data:
C         ISTATE  Designation of state I - the bra vector
C         JSTATE  Designation of state J - the ket vector
C                 where both above are wrt the unique target state table
C          IPCDE  Slater property integral codes (8 per property)
C          PROPS  Expectation values for the wavefunction pair.
C
C     Example:
C
C         ISTATE   JSTATE         IPCDE         PROPS
C         ------   ------        -------        -----
C            1        1      2 0 0 0 0 0 0 0     1.0
C            2        1      2 0 0 0 0 0 0 0     0.0
C            2        2      2 0 0 0 0 0 0 0     1.0
C
C     where property 2 0 0 0 0 0 0 0 is the overlap, i.e.
C
C                 < State 1 | State 2 >  =  0.0
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: MAXPROPS
      INTEGER, DIMENSION(8*maxprops) :: IPCDE
      INTEGER, DIMENSION(maxprops) :: ISTATE, JSTATE
      REAL(KIND=wp), DIMENSION(maxprops) :: PROPS
      INTENT (IN) MAXPROPS
      INTENT (OUT) IPCDE, ISTATE, JSTATE, PROPS
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
C---- Loop over the rows of the table's columns and initialize them.
C
      DO I=1, MAXPROPS
         ISTATE(I)=0
         JSTATE(I)=0
      END DO
C
C---- To enhance vectorization we treat IPCDE as a single dimension
C     array in this case.
C
      DO I=1, 8*MAXPROPS
         IPCDE(I)=0
      END DO
C
C---- Similarly PROPS array.
C
      DO I=1, MAXPROPS
         PROPS(I)=xzero
      END DO
C
      RETURN
C
      END SUBROUTINE INIPTAB
!*==initjr.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE INITJR(NELT,NDTRF,MN,NC,JROB,JRON,norb)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NC, NELT, NORB
      INTEGER, DIMENSION(norb) :: JROB, JRON
      INTEGER, DIMENSION(*) :: MN
      INTEGER, DIMENSION(nelt) :: NDTRF
      INTENT (IN) MN, NDTRF, NELT, NORB
      INTENT (OUT) JROB
      INTENT (INOUT) JRON, NC
C
C Local variables
C
      INTEGER :: I, M, MC, N
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     SETS UP JRON AND JROB FROM THE REFERENCE DETERMINANT NDTRF OF
C     SPIN ORBITALS
C
C     INPUT :
C
C       NELT    NUMBER OF ELECTRONS IN SYSTEM
C      NDTRF    REFERENCE DETERMINANT OF SPIN ORBITALS     (ARRAY)
C         MN    SYMMETRY OF THE SPIN ORBITALS              (ARRAY)
C
C     OUTPUT :
C
C         NC   NUMBER OF ORBITALS IN THE REFERENCE DETERMINANT
C       JROB   LIST OF ORBITALS
C       JRON   LIST OF ORBITAL OCCUPATION NUMBERS
C
C***********************************************************************
C
C     INITIALIZE VARIABLES FOR EACH PASS THROUGH THE ROUTINE
C
      NC=0
      MC=0
C
C     LOOP OVER ELECTRONS TO CONSTRUCT JRON AND JROB FROM THE
C     REFERENCE DETERMINANT
C
      DO I=1, NELT
         N=NDTRF(I)
         M=MN(N)
         IF(M.NE.MC)THEN
            NC=NC+1
            JROB(NC)=M
            JRON(NC)=1
            MC=M
         ELSE
            JRON(NC)=JRON(NC)+1
         END IF
      END DO
C
      RETURN
C
      END SUBROUTINE INITJR
!*==initst.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE INITST(ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,MAXTGT)
C***********************************************************************
C
C     INITST - INItialize columns of the unique Target STate table
C
C     Input/Output data:
C          ISPIN 2*S+1 for each target state
C            ISZ 2*Sz for each target state
C           MGVN Lamda value (C-inf-v) or Irred. Rep of each state
C          GUTOT For D-inf-h only, the gerade or ungerade value
C          IREFL For C-inf-v only the +/- sigma character
C         ENERGY Eigen-energy in Hartrees for the each state
C         IWRITE Logical unit for the printer
C         MAXTGT Maximum nuber of rows in the target state table
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: MAXTGT
      REAL(KIND=wp), DIMENSION(maxtgt) :: ENERGY
      INTEGER, DIMENSION(maxtgt) :: GUTOT, IREFL, ISPIN, ISZ, MGVN
      INTENT (IN) MAXTGT
      INTENT (OUT) ENERGY, GUTOT, IREFL, ISPIN, ISZ, MGVN
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
C---- Zeroize the columns of the table
C
      DO I=1, MAXTGT
         ISPIN(I)=0
         ISZ(I)=0
         MGVN(I)=0
         GUTOT(I)=0
         IREFL(I)=0
         ENERGY(I)=XZERO
      END DO
C
      RETURN
C
      END SUBROUTINE INITST
!*==isrche.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
 
      FUNCTION ISRCHE(N,NAR,INC,ITAR)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: INC, ITAR, N
      INTEGER :: ISRCHE
      INTEGER, DIMENSION(n) :: NAR
      INTENT (IN) NAR
C
C Local variables
C
      INTEGER :: I
      INTEGER :: IHJSR
      INTEGER, DIMENSION(n+1) :: IWORK
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     ISRCHE - PREPARE CALL TO INTEGER SEARCH FUNCTION
C
C              SETS UP A CALL TO THE FUNCTION IHJSR WHICH LOCATES THE
C              TARGET INTEGER IN THE SEQUENCE NAR.
C
C     INPUT DATA :
C              N  DIMENSION OF THE SEQUENCE TO BE SEARCHED
C            NAR  SEQUENCE OF INTEGERS TO BE SEARCHED. ITS ACTUAL
C                 DIMENSION IS ONE LARGER THAN THE SEQUENCE SIZE.
C            INC  STARTING POSITION IN SEQUENCE
C           ITAR  THE INTEGER BEING SOUGHT IN THE SEQUENCE IE. THE
C                 TARGET
C
C     OUTPUT DATA :
C          ISRCHE  THE POSITION OF THE TARGET IN THE SEQUENCE IS
C                  RETURNED VIA THE FUNCTION NAME. IF THE TARGET IS NOT
C                  FOUND THEN THE VALUE OF ISRCHE IS ONE LARGER THAN
C                  THE SIZE OF THE SEQUENCE.
C
C     LOCAL DATA :
C          IWORK  THIS ARRAY IS A WORKSPACE INTO WHICH IS COPIED THE
C                 ARRAY NAR WHICH IS BEING SEARCHED. IT IS IMPORTANT
C                 THAT THE ARRAY NAR IS NOT OVERWRITTEN AS THERE MAY
C                 BE ELEMENTS BEYOND THE LIMITS WHICH WE ARE SEARCHING
C
C***********************************************************************
C
      DO I=1, N
         IWORK(I)=NAR(I)
      END DO
c
      ISRCHE=IHJSR(N,IWORK,INC,ITAR)
C
      RETURN
C
      END FUNCTION ISRCHE
!*==makemg.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE MAKEMG(MG,NSRB,NELT,NDTRF)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NELT, NSRB
      INTEGER, DIMENSION(nsrb) :: MG
      INTEGER, DIMENSION(nelt) :: NDTRF
      INTENT (IN) NDTRF, NELT, NSRB
      INTENT (OUT) MG
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     ASSOCIATES WITH EACH SPIN ORBITAL IN THE REFERENCE DETERMINANT
C     A SEQUENCE NUMBER GIVING THE CONTOGUOUS POSITION OF THAT SPIN
C     ORBITAL IN THE REFERENCE DETERMINANT. NOTE THAT THIS ROUTINE IS
C     COPIED FROM SPEEDY.
C
C     INPUT DATA:
C           NSRB THE NUMBER OF SPIN ORBITALS IN THE WAVEFUNCTION SET
C           NELT THE NUMBER OF ELECTRONS IN THE WAVEFUNCTION
C          NDTRF THE REFERENCE DETERMINANT FOR THE WAVEFUNCTION
C
C     OUTPUT DATA:
C              MG THE POINTER ARRAY ITSELF.
C
C***********************************************************************
C
      DO I=1, NSRB
         MG(I)=0
      END DO
C
      DO I=1, NELT
         MG(NDTRF(I))=I
      END DO
C
      RETURN
C
      END SUBROUTINE MAKEMG
!*==mkorbs.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE MKORBS(ISYMTYP,NOB,NSYM,MN,MG,MM,MS,maxk,NPFLG)
C***********************************************************************
C
C     MKORBS - Computes the orbital table which is then used in the
C              evaluation of off diagonal density matrix elements
C
C     Input data:
C        ISYMTYP  Switch for C-inf-v (=0 or 1) / Abelian point group (=2
C            NOB  Number of orbitals per symmetry
C           NSYM  Number of symmetries in the orbital set
C          NPFLG  Flag controlling printing of computed orbital table
C
C     Output data:
C              MN  Orbital number associated with each spin-orbital
C              MG  G/U designation for each spin-orbital (C-inf-v only)
C                  Actually this is always zero because C-inf-v does not
C                  distinguish between g/u. It exists because original
C                  version of Alchemy tried to use it for D-inf-h too;
C                  all CI evauation is doen in C-inf-v now because CONGE
C                  converts D-inf-h to C-inf-v data.
C              MM  Symmetry quantum number associated with each spin-orb
C              MS  Spin function ( alpha or beta ) associated with each
C                  spin orbital
C
C     Notes:
C
C       The orbital table establishes orbital and quantum number data fo
C     each spin orbital in the set.
C
C     e.g. C-inf-v symmetry with NSYM=2, NOB=3,1, yields ten spin
C          orbitals which are designated as follows by this routine:
C
C        Spin orb.     MN  MG  MM  MS     Comments
C            1          1   0   0   0     1 sigma spin up
C            2          1   0   0   1     1 sigma spin down
C            3          2   0   0   0     2 sigma spin up
C            4          2   0   0   1     2 sigma spin down
C            5          3   0   0   0     3 sigma spin up
C            6          3   0   0   1     3 sigma spin down
C            7          4   0   1   0     1 pi(lambda=+1) spin up
C            8          4   0   1   1     1 pi(lambda=+1) spin down
C            9          4   0  -1   0     1 pi(lambda=-1) spin up
C           10          4   0  -1   1     1 pi(lambda=-1) spin down
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ISYMTYP, MAXK, NPFLG, NSYM
      INTEGER, DIMENSION(*) :: MG, MM, MN, MS, NOB
      INTENT (IN) ISYMTYP, NOB, NPFLG, NSYM
      INTENT (OUT) MAXK
      INTENT (INOUT) MG, MM, MN, MS
C
C Local variables
C
      INTEGER :: I, J, K, M, MA, MB, N, NSRB
C
C*** End of declarations rewritten by SPAG
C
C---- 'Afore we compute NSRB we set it to -1 so that we can error trap
C     this routine. Should something happen such that the loops are
C     not executed then NSRB is never reset to a proper value and so
C     the test fails at the end of the routine.
C
      NSRB=-1
C
      IF(ISYMTYP.LT.2)THEN
C
C=======================================================================
C
C     Linear molecules i.e. C-inf-v point group
C
C=======================================================================
C
C---- Sigma type orbital are done first.
C
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
C
C---- Now do Pi, Delta, Phi orbitals etc.....
C
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
C
      ELSE
C
C=======================================================================
C
C     Abelian point group molecules (D2h and subgroups)
C
C=======================================================================
C
C---- Loop over all symmetries
C          Loop over all orbitals within the symmetry
C               Assign numbers for both spin orbitals associated
C               with each orbital.
C
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
C
      END IF
C
      NSRB=I-1
      maxk=k-1
C
C---- Error check for a stupid value of NSRB.
C
      IF(NSRB.LE.0)THEN
         WRITE(6,9900)
         WRITE(6,9940)NSRB
         CALL TMTCLOS()
      END IF
C
      IF(NPFLG.GE.1)THEN
         WRITE(6,505)
         WRITE(6,500)(I,MN(I),MG(I),MM(I),MS(I),I=1,NSRB)
         WRITE(6,510)NSRB
      END IF
C
      RETURN
C
C---- Format Statements
C
 500  FORMAT(8X,'I',4X,'N',4X,'G',4X,'M',4X,'S',/,(4X,5I5))
 505  FORMAT(//,5X,' THE ORBITAL TABLE FOR THIS CASE IS :',/)
 510  FORMAT(5X,' NO. OF SPIN ORBITALS FOR THIS RUN IS',I8,/)
C
 9900 FORMAT(//,10X,'**** Error in MKORBS ',//)
 9940 FORMAT(10X,'NSRB is invalid = ',I10,/)
C
      END SUBROUTINE MKORBS
!*==mntab.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE MNTAB(ISYMTYP,NOB,NSYM,MN,NSRB,NPFLG)
C***********************************************************************
C
C     CONSTRUCTS THE MN TABLE FROM THE VALUES OF NOB AND NSYM. THE
C     TOTAL NUMBER OF SPIN ORBITALS IS RETURNED VIA NSRB, ALSO.
C
C     THE MN TABLE GIVES THE ORBITAL SEQUENCE NUMBER ASSOCIATED WITH
C     EACH SPIN ORBITAL.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ISYMTYP, NPFLG, NSRB, NSYM
      INTEGER, DIMENSION(*) :: MN, NOB
      INTENT (IN) ISYMTYP, NOB, NPFLG, NSYM
      INTENT (INOUT) MN, NSRB
C
C Local variables
C
      INTEGER :: I, I1, I2, MP, MS, N
      LOGICAL :: ZFAIL
C
C*** End of declarations rewritten by SPAG
C
C---- Error checking of input data.
C
C     We must have that NSYM is positive and that the number of orbitals
C     per symmetry is positive too.
C
      IF(NSYM.LE.0)THEN
         WRITE(6,9900)
         WRITE(6,9910)NSYM
         CALL TMTCLOS()
      END IF
C
C..... Scan the orbital set now
C
      ZFAIL=.FALSE.
      DO I=1, NSYM
         IF(NOB(I).LT.0)ZFAIL=.TRUE.
      END DO
C
      IF(ZFAIL)THEN
         WRITE(6,9900)
         WRITE(6,9920)(NOB(I),I=1,NSYM)
         CALL TMTCLOS()
      END IF
C
C---- 'Afore we compute NSRB we set it to -1 so that we can error trap
C     this routine. Should something happen such that the loops are
C     not executed then NSRB is never reset to a proper value and so
C     the test fails at the end of the routine.
C
      NSRB=-1
C
C---- Depending upon the value of ISYMTYP we branch to an appropriate
C     in-line subroutine, or else we catch an error.
C
      IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)GO TO 1000
      IF(ISYMTYP.EQ.2)GO TO 2000
C
      WRITE(6,9900)
      WRITE(6,9930)ISYMTYP
      CALL TMTCLOS()
C
C---- Subroutine return point
C
 800  CONTINUE
C
C---- Back end error check for a stupid value of NSRB.
C
      IF(NSRB.LE.0)THEN
         WRITE(6,9900)
         WRITE(6,9940)NSRB
         CALL TMTCLOS()
      END IF
C
C---- Debug print out, if requested by the user
C
      IF(NPFLG.GT.0)THEN
         WRITE(6,500)
         WRITE(6,505)
         WRITE(6,506)(I,MN(I),I=1,NSRB)
         WRITE(6,510)NSRB
      END IF
C
      RETURN
C
C=======================================================================
C
C     In-line subroutine for C-inf-v (D-inf-h) Molecules
C
C=======================================================================
C
 1000 CONTINUE
C
C---- Count MS, the number of sigma orbitals, and MP, the number of
C     non-sigma orbitals. This is needed because sigma orbitals are
C     doubly degenerate while non sigma orbitals are quadruply
C     degenerate
C
      MS=NOB(1)
      MP=0
C
      DO I=2, NSYM
         MP=MP+NOB(I)
      END DO
C
      NSRB=2*MS+4*MP
C
C     CONSTRUCT THE SIGMA PART OF THE MN TABLE
C
      I2=0
      N=0
      IF(MS.GT.0)THEN
         I1=I2+1
         I2=I2+2*MS
         DO I=I1, I2, 2
            N=N+1
            MN(I)=N
            MN(I+1)=N
         END DO
      END IF
C
C     CONSTRUCT THE OTHER PART OF THE MN TABLE
C
      IF(MP.GT.0)THEN
         I1=I2+1
         I2=I2+4*MP
         DO I=I1, I2, 4
            N=N+1
            MN(I)=N
            MN(I+1)=N
            MN(I+2)=N
            MN(I+3)=N
         END DO
      END IF
C
C---- Return to the main section of this subroutine
C
      GO TO 800
C
C=======================================================================
C
C     In-line subroutine for Abelian point group Molecules
C
C     (D2h and subgroups)
C
C=======================================================================
C
 2000 CONTINUE
C
C---- Count the total number of orbitals
C
      MS=0
      DO I=1, NSYM
         MS=MS+NOB(I)
                   ! possibly problematic
      END DO
C
C---- Construct the orbital table
C
      NSRB=2*MS
C
      IF(MS.GT.0)THEN
         N=0
         DO I=1, NSRB, 2
            N=N+1
            MN(I)=N
            MN(I+1)=N
         END DO
      END IF
C
C---- Return to the main part of this subroutine
C
      GO TO 800
C
C---- Format Statements
C
 500  FORMAT(//,5X,'The Orbital Table : ',//,5X,
     &       'This equates a sequential orbital number with',/,5X,
     &       'every spin orbital. The orbital indices are  ',/,5X,
     &       'not reset to zero in a new symmetry.         ',//)
 505  FORMAT(5X,'Spin Orbital',1X,' Orbital  ',/,5X,' Sequential ',1X,
     &       'Sequential',/,5X,'    Index   ',1X,' Index    ',/5X,
     &       '------------',1X,'----------',/)
 506  FORMAT(5X,I6,6X,1X,I6)
 510  FORMAT(/,5X,'Total number of spin orbitals= ',I5,/)
 9900 FORMAT(//,10X,'**** Error in MNTAB ',//)
 9910 FORMAT(10X,'Value of NSYM input is invalid = ',I10,/)
 9920 FORMAT(10X,'NOB array has one or more errors:',//,
     &       (10X,10(I3,1X),/))
 9930 FORMAT(10X,'ISYMTYP not 0,1 or 2 but = ',I5,/)
 9940 FORMAT(10X,'Number of computed Spin-orbs .le. zero = ',I5,/)
C
      END SUBROUTINE MNTAB
!*==modrda.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE MODRDA(NELT,NSRB,NDTRI,NDTRJ,MDI,NODA,MDA,NDI,CA)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: MDA, NELT, NODA, NSRB
      REAL(KIND=wp), DIMENSION(noda) :: CA
      INTEGER, DIMENSION(*) :: MDI, NDI
      INTEGER, DIMENSION(nelt) :: NDTRI, NDTRJ
      INTENT (IN) NDI, NDTRI, NDTRJ, NELT, NODA, NSRB
      INTENT (OUT) MDI
      INTENT (INOUT) CA, MDA
C
C Local variables
C
      INTEGER :: I, J, M, MDB, MDC, N, ND, NDA, NDB, NDC, NP, NPS, NQ, 
     &           NR, NS, NT
      INTEGER, DIMENSION(nsrb) :: MA, MB, NA, NB, NTT
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     MODRDA - MODIFY REFERENCE DETERMINANTS FOR THE DETERMINANTAL
C              EXCITATION CODES.
C
C     IN CALCULATING THE OFF-DIAGONAL ELEMENTS OF THE DENSITY MATRIX
C     IT IS SOMETIMES NECESSARY TO COMPARE DETERMINANTS WHICH HAVE BEEN
C     CALCULATED WITH RESPECT TO DIFFERENT REFERENCE DETERMINANTS.SINCE
C     THE DENSITY MATRIX GENERATION ALGORITHM EMPLOYED IN THIS CODE
C     ASSUMES THAT THE DETERMINANTAL EXCITATION CODES FOR THE TWO GIVEN
C     WAVEFUNCTIONS HAVE BEEN GENERATED FROM THE SAME REFERENCE
C     DETERMINANT IT IS NECESSARY TO MODIFY THE DETERMINANTAL EXCITATION
C     CODES OF ONE OF THE WAVEFUNCTIONS TO BE WITH RESPECT TO THE
C     REFERENCE DETERMINANT OF THE OTHER. THIS SUBROUTINE CARRIES OUT
C     THIS TASK THROUGH THE ENTRY POINT MODRDB.
C
C     INPUT DATA :
C           NELT  NUMBER OF ELECTRONS
C           NSRB  NUMBER OF SPIN ORBITALS
C          NDTRI  REFERENCE DETERMINANT I
C          NDTRJ  REFERENCE DETERMINANT J
C           NODA  NUMBER OF DETERMINANTS IN WAVEFUNCTION J
C            NDI  THE EXCITATION CODES FOR THE DETERMINANTS OF
C                 WAVEFUNCTION J RELATIVE TO REFERENCE DETERMINANT I
C             CA  COEFFICIENTS OF THE EXCITATION CODES IN WAVEFUNCTION J
C
C    OUTPUT DATA :
C             NP  NUMBER OF UN-MATCHED ORBITALS SPIN ORBITALS
C            NPS  PHASE FACTOR DETERMINED BY NUMBER OF UNMATCHED
C                 ORBITALS
C             NA  SPIN ORBITALS IN REFERENCE DETERMINANT I THAT ARE NOT
C                 FOUND IN REFERENCE DETERMINANT J
C             NB  SPIN ORBITALS IN REFERENCE DETERMINANT J THAT ARE NOT
C                 FOUND IN REFERENCE DETERMINANT I
C            MDI  THE NEW EXCITATION CODES FOR THE DETERMINANTS OF
C                 WAVEFUNCTION J RELATIVE TO REFERENCE DETERMINANT I
C            MDA  SIZE OF THE ARRAY MDI
C             CA  THE COEFFICIENTS FOR EACH OF THE NEW EXCITATION CODES
C                 THESE ARE THE SAME AS THE INPUT VALUES BUT POSSIBLY
C                 MULTIPLIED BY -1.
C
C    TEMPORARY STORAGE AREAS:
C            NTT  STORES THE REFERENCE DETERMINANT J DURING COMPARISON
C             MA
C             MB
C
C***********************************************************************
C
C-----REF DET J IS COPIED TO TEMPORARY AREA NTT;VARIABLES ARE INITIALIZE
C
      DO I=1, NELT
         NTT(I)=NDTRJ(I)
      END DO
      NP=0
      NPS=1
C
C-----REF DET I IS SCANNED AND COMPARED TO REF DET J. AN ATTEMPT IS MADE
C     TO MATCH EACH SPIN ORBITAL IN REF DET I WITH ITSELF IN REF DET J.
C     NO MATCH IS FOUND THEN INFORMATION IS  STORED IN ARRAYS
C     NA (THE SO NUMBER) AND NB (ITS POSITION IN REF DET I). IF A MATCH
C     FOUND THEN IF THE SO OCCURS AT THE SAME POSITION IN EACH REF DET N
C     ACTION IS NECESSARY. ON THE OTHER HAND IF THE SO OCCURS AT DIFFERE
C     POSITIONS THEN REF DET J IS RE-ORDERED TO HAVE THE SO IN THE SAME
C     POSITION AS IN REF DET I. A CHANGE OF PHASE IS REQUIRED. THIS
C     APPROACH ALSO FACILITATES THE GENERATION OF AN EXCITATION CODE FOR
C     REF DET J W.R.T. REF DET I. IT IS ONLY NECESSARY TO FIND THE SO'S
C     IN REF DET I WHICH ARE UN-MATCHED AND THEN LOOK AT THE SO'S IN THE
C     CORRESPONDING POSITIONS IN REF DET J.
C
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
C
C-----CHANGE NB TO BE SO'S RATHER THAN POINTERS TO STORAGE LOCATIONS
C
      DO I=1, NP
         M=NB(I)
         NB(I)=NTT(M)
      END DO
C
C-----ARRAY NTT IS RE-INITIALIZED FOR SUBSEQUENT USE
C
      DO I=1, Nsrb
         NTT(I)=0
      END DO
C
C-----RE-ORGANIZATION OF THE DETERMINANTS OF WAVEFUNCTION J BEGINS HERE
C
      MDA=1
      NDA=1
C
C-----IF THE REFERENCE DETERMINANTS FOR BOTH WAVEFUNCTIONS ARE EQUAL,THE
C     JUST COPY DETERMINANTS FROM NDI TO MDI.NO MODIFICATION IS REQUIRED
C
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
C
C-----MODIFICATION IS REQUIRED.LOOP OVER DETERMINANTS AND MODIFY EACH ON
C
      DO ND=1, NODA
         NDB=NDI(NDA)
         NDC=NDA+NDB
         NS=NPS
C
C--------THERE ARE NO REPLACEMENTS FROM REFERENCE DETERMINANT J FOR THIS
C        PARTICULAR DETERMINANT. EXCITATION CODE IS JUST THAT FOR REFDET
C        W.R.T. REF DET I.
C
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
C
C--------THERE ARE REPLACEMENTS FROM REFERENCE DETERMINANT J
C        PERFORM INITIALIZATION FOR EACH DETERMINANT.
C
         DO I=1, NP
            MA(I)=NA(I)
            MB(I)=NB(I)
         END DO
C
         NR=NP
         DO I=1, NDB
            MA(NR+I)=NDI(NDA+I)
            MB(NR+I)=NDI(NDC+I)
         END DO
         NR=NP+NDB
C
         DO I=1, NR
            NTT(MA(I))=I
         END DO
C
C--------MODIFY THE EXCITATION CODES. SPIN ORBITALS TO BE REPLACED ARE
C        IN ARRAY MA WHILE THOSE THAT REPLACE ARE IN MB
C
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
C
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
C
C--------STORE THE NEW EXCITATION CODES INTO ARRAY MDI
C
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
C
C--------MODIFY THE PHASE FACTOR OF THE DETERMINANT
C
 230     IF(NS.LT.0)CA(ND)=-CA(ND)
C
C--------AUGUMENT STORAGE POSITIONS IN THE DETERMINANT ARRAYS
C        FOR THE NEXT PASS
C
         MDA=MDC+MDB+1
         NDA=NDC+NDB+1
      END DO
C
      MDA=MDA-1
c
      RETURN
C
      END SUBROUTINE MODRDA
!*==modron.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE MODRON(IOCC,IORB,NZ,ND,NDI,NRP,MN,norb)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ND, NORB, NRP, NZ
      INTEGER, DIMENSION(norb) :: IOCC, IORB
      INTEGER, DIMENSION(*) :: MN
      INTEGER, DIMENSION(nrp) :: NDI
      INTENT (IN) MN, ND, NDI, NORB, NRP
      INTENT (INOUT) IOCC, IORB, NZ
C
C Local variables
C
      INTEGER :: I, J, M, N
      LOGICAL :: USED
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     MODRON - MODIFIES OCCUPATION NUMBERS OF ORBITALS.
C
C              THE FUNDAMENTAL ALGORITHM BEHIND THE DIAGONAL DENSITY
C              MATRIX CONSTRUCTION RELIES ON THE RELATIVE OCCUPATION
C              NUMBERS OF ORBITALS IN DETERMINANT OVERLAPS. THIS ROUTINE
C              IS CALLED FORM SEVERAL PLACES WITHIN THE DIAGONAL ELEMENT
C              CONSTRUCTION ROUTINES.
C
C     INPUT DATA :
C             ND  THE NUMBER BY WHICH THE ORBITAL POPULATION IS TO BE
C                 CHANGED. WHEN CALLED FOR REPLACEMENTS ND=-1 IN GENERAL
C                 BUT WHEN CALLED FOR REPLACING ORBITALS ND=1.
C           IOCC  ORBITAL OCCUPATION TABLE
C           IORB  ORBITAL NUMBERS FOR ORBITAL OCCUPATION TABLE
C             NZ  CURRENT SIZE OF IN AND IR ARRAYS
C            NDI  ACTUAL DETERMINANTS ON WHICH TO WORK
C            NRP  SIZE OF THE DETERMINANTS
C             MN  ORBITAL TABLE
C
C    OUTPUT DATA :
C             IOCC,IORB AND NZ ARE SUITABLY MODIFIED
C
C***********************************************************************
C
C-----LOOP OVER THE SPECIFIED PART OF THE DETERMINANT AND MODIFY THE
C     ORBITAL OCCUPATION NUMBER.
C
      DO I=1, NRP
         N=NDI(I)
         M=MN(N)
         USED=.FALSE.
         IF(NZ.EQ.0)THEN
            NZ=NZ+1
            IORB(NZ)=M
            IOCC(NZ)=ND
            USED=.TRUE.
            CYCLE
         END IF
         DO J=1, NZ
            IF(IORB(J).NE.M)CYCLE
            IOCC(J)=IOCC(J)+ND
            USED=.TRUE.
         END DO
         IF(.NOT.USED)THEN
            NZ=NZ+1
            IORB(NZ)=M
            IOCC(NZ)=ND
         END IF
      END DO
C
      RETURN
C
      END SUBROUTINE MODRON
!*==moi.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE MOI(rmoi,VLIST,NOC,NOCMX,IWRITE,ipflg)
C***********************************************************************
C
C     MOI - MOment of Inertia tensor
C
C     Input data:
C          VLIST  Data on each nuclear center in compressed format
C            NOC  Number of centers in the problem
C          NOCMX  Maximum number of centers in the VLIST table
C         IWRITE  Logical unit for the printer
C
C     Linkage:
C
C        JACOBI
C
C     Notes:
C
C       This routine was obtained from the NYU properties code in the
C     form used in the Alchemy II package. It has been recoded in
C     standard Fortran 77 by Charles J Gillan at Queen's Belfast.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO, XONE
      USE blas_lapack_gbl, ONLY : blasint
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IPFLG, IWRITE, NOC, NOCMX
      REAL(KIND=wp), DIMENSION(3) :: RMOI
      REAL(KIND=wp), DIMENSION(NOCMX,8) :: VLIST
      INTENT (IN) IPFLG, IWRITE, NOC, NOCMX, VLIST
      INTENT (INOUT) RMOI
C
C Local variables
C
      REAL(KIND=wp) :: AMU, ANG, CM1, EV, TOTMSS, WT, WTR2, XM, YM, ZM
      REAL(KIND=wp), DIMENSION(3) :: COM
      INTEGER :: I, IJ, J, JJ, JK, K, KK
      INTEGER(blasint) :: INFO, THREE=3, TEN=10
      REAL(KIND=wp), DIMENSION(3,3) :: PATRAN
      REAL(KIND=wp), DIMENSION(10)  :: TENMOI
C
C*** End of declarations rewritten by SPAG
C
      DATA amu/1822.887_wp/, cm1/219474.6_wp/, eV/27.21165_wp/, 
     &     ang/0.2798_wp/
C
C---- Evaluate the center of mass
C
      TOTMSS=XZERO
      XM=XZERO
      YM=XZERO
      ZM=XZERO
C
      DO I=1, NOC
         TOTMSS=TOTMSS+VLIST(I,6)
         XM=XM+VLIST(I,6)*VLIST(I,1)
         YM=YM+VLIST(I,6)*VLIST(I,2)
         ZM=ZM+VLIST(I,6)*VLIST(I,3)
      END DO
C
      COM(1)=XM/TOTMSS
      COM(2)=YM/TOTMSS
      COM(3)=ZM/TOTMSS
C
      IF(ipflg.NE.0)WRITE(IWRITE,2050)TOTMSS, (COM(I),I=1,3)
C
C---- Now determine the moment of Inertia tensor and principal axis
C     of the system.
C
      DO IJ=1, 6
         TENMOI(IJ)=XZERO
      END DO
C
      DO I=1, NOC
         WT=VLIST(I,6)
         WTR2=WT*((VLIST(I,1)-COM(1))**2+(VLIST(I,2)-COM(2))**2
     &        +(VLIST(I,3)-COM(3))**2)
         DO J=1, 3
            JJ=J*(J+1)/2
            TENMOI(JJ)=TENMOI(JJ)+WTR2-WT*(VLIST(I,J)-COM(J))**2
            DO K=1, J-1
               JK=J*(J-1)/2+K
               TENMOI(JK)=TENMOI(JK)+WT*(VLIST(I,J)-COM(J))
     &                    *(VLIST(I,K)-COM(K))
            END DO
         END DO
      END DO
C
C.... Write out the Moment of Inertia Tensor for the user
C
      IF(ipflg.NE.0)THEN
         WRITE(IWRITE,5010)
         WRITE(IWRITE,7020)TENMOI(1)
         WRITE(IWRITE,7020)TENMOI(2), TENMOI(3)
         WRITE(IWRITE,7020)TENMOI(4), TENMOI(5), TENMOI(6)
      END IF
C
C.... Diagonalize the moment of inertia tensor to obtain the principal
C     axis.PATRAN now contains the moment of inertia tensor.
C
      KK=0
      DO I=1, 3
         DO J=1, I
            KK=KK+1
            PATRAN(I,J)=TENMOI(KK)
         END DO
      END DO
C
c     TENMOI is now a workspace. Value 10 for lwork was chosen
c     to satisfy minimal requirement of 3*N-1 elements with a small
c     additional margin.
      CALL DSYEV('V','L',THREE,PATRAN,THREE,RMOI,TENMOI,TEN,INFO)
C
C---- Print a summary of the Principal Moments of Inertia for the user
C     including the principal axes transformation matrix
C
      IF(ipflg.NE.0)THEN
         WRITE(IWRITE,7010)
         WRITE(IWRITE,7020)(amu*RMOI(I),I=1,3)
         WRITE(IWRITE,7011)
         WRITE(IWRITE,7020)(ang*RMOI(I),I=1,3)
C
         WRITE(IWRITE,7060)
         WRITE(IWRITE,7061)((PATRAN(I,J),I=1,3),J=1,3)
      END IF
C
C---- At this stage RMOI contains the principal moments of inertia;
C     take the reciprocals and convert to GigaHertz.
C
C     For a reference see the book
C
C         Molecular Quantum Mechanics by P W Atkins  2nd Edition p.289
C         Pub. Oxford University Press ISBN 0-19-855170-3
C
      DO I=1, 3
         IF(ABS(RMOI(I)).GT.1.0E-01_wp)THEN
            RMOI(I)=0.5_wp/(amu*RMOI(I))
         ELSE
            RMOI(I)=xzero
         END IF
      END DO
C
C---- Print a summary of the rotational frequencies
C
      IF(ipflg.NE.0)THEN
         WRITE(IWRITE,7050)
         WRITE(IWRITE,7020)(RMOI(I),I=1,3)
         WRITE(IWRITE,7051)
         WRITE(IWRITE,7020)(cm1*RMOI(I),I=1,3)
         WRITE(IWRITE,7052)
         WRITE(IWRITE,7020)(eV*RMOI(I),I=1,3)
      END IF
 
C
      RETURN
C
C---- Format Statements
C
 2050 FORMAT(/,10X,'Total Molecular Mass = ',F12.3,' Dalton (amu)',//,
     &       10X,'Co-ordinates of the Center of Mass in the present',
     &       ' axis system: ',//,10X,3(F12.5,1X),/)
 5010 FORMAT(/,10X,'Moment of Inertia Tensor (Lower Half Triangle):',/)
 7010 FORMAT(/,10X,'Principal Moments of Inertia (Atomic units): ',/)
 7011 FORMAT(/,10X,'Principal Moments of Inertia (amu angstrom^2): ',/)
 7020 FORMAT(10X,3(F20.10,1X),/)
 7050 FORMAT(/,10X,'Rotational Constants (Atomic Units): ',/)
 7051 FORMAT(/,10X,'Rotational Constants (cm**(-1)): ',/)
 7052 FORMAT(/,10X,'Rotational Constants (eV): ',/)
 7060 FORMAT(/,10X,'Principal Axes Transformation Matrix: ',/)
 7061 FORMAT(10X,F12.5,' x ',F12.5,' y ',F12.5,' z ')
C
      END SUBROUTINE MOI
!*==moidriv.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE MOIDRIV(NNUC,RGEOM,CHARGE,IWRITE)
C***********************************************************************
C
C     MOIDRIV - MOMents Of Inertia DRIVer routine
C
C     Collates the Alchemy style data into the format used by the
C     SWEDEN routine INCENT for computation of the moments of inertia
C     and rotational frequencies
C
C     Input data:
C           NNUC Number of nuclear centers in the molecule
C          RGEOM X,Y,Z co-ordinates for each nuclear center
C         CHARGE Charge on each nuclear center
C         IWRITE Logical unit for the printer
C             CR Real*8 workspace area for dynamic table building
C           NCOR Size of the core space in R*8 words
C
C     Linkage:
C
C          INCENT
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, NNUC
      REAL(KIND=wp), DIMENSION(nnuc) :: CHARGE
      REAL(KIND=wp), DIMENSION(3,nnuc) :: RGEOM
      INTENT (IN) CHARGE, IWRITE, RGEOM
C
C Local variables
C
      REAL(KIND=wp) :: CM1, Z
      REAL(KIND=wp), DIMENSION(8*nnuc) :: CR
      INTEGER :: DEND, I, LTABLE, NOADC, NOCMX
      REAL(KIND=wp), DIMENSION(3) :: RMOI
      CHARACTER(LEN=8), DIMENSION(nnuc) :: VSYMB
C
C*** End of declarations rewritten by SPAG
C
C---- Following saves space because the scope of the variables is
C     different
C
      EQUIVALENCE(LTABLE,NOADC)
C
C
      DATA cm1/219474.6_wp/
      DATA dend/288/
C
C---- If we have an atomic system then there is no Moment of Inertia
C     tensor to evaluate !
C
      IF(NNUC.EQ.1)THEN
         WRITE(IWRITE,9911)
         GO TO 800
      END IF
C
      LTABLE=NNUC*8
      DO I=1, LTABLE
         CR(I)=XZERO
      END DO
C
C---- Prepare the data for the SWEDEN routine INCENT.
C
C       It is useful to remember that:
C
C          NOADC is the no. of additional centers, always zero in Alchem
C          NOCMX is the maximum number of rows in the nuclear table
C
      NOADC=0
      NOCMX=NNUC
C
C..... The nuclear data table is stored in row order with NNUC entries p
C      column !
C
      DO I=1, NNUC
         CR(I)=RGEOM(1,I)
         CR(I+NNUC)=RGEOM(2,I)
         CR(I+NNUC*2)=RGEOM(3,I)
         CR(I+NNUC*3)=CHARGE(I)
      END DO
C
C---- Having prepared the data now invoke the routine from the SWEDEN
C     package which performs the computation
C
      CALL INCENT(CR(1),vsymb,NOCMX,NNUC,NOADC,rmoi,1)
!
      z=0._wp
! open borndat (fort.288) to write rotational constants for BORNCROS
      OPEN(UNIT=dend,FILE='borndat',STATUS='unknown')
      IF(RMOI(1).NE.z .AND. RMOI(2).NE.z .AND. RMOI(3).NE.z)THEN
         IF(RMOI(1).GT.RMOI(2) .AND. RMOI(1).GT.RMOI(3))THEN
            WRITE(dend,7021)cm1*RMOI(1)
         ELSE IF(RMOI(2).GT.RMOI(1) .AND. RMOI(2).GT.RMOI(3))THEN
            WRITE(dend,7021)cm1*RMOI(2)
         ELSE IF(RMOI(3).GT.RMOI(1) .AND. RMOI(3).GT.RMOI(2))THEN
            WRITE(dend,7021)cm1*RMOI(3)
         ELSE
            WRITE(dend,7021)cm1*RMOI(1)
         END IF
      ELSE IF(RMOI(1).EQ.z .AND. RMOI(2).NE.z .AND. RMOI(3).NE.z)THEN
         IF(RMOI(2).GT.RMOI(3))THEN
            WRITE(dend,7021)cm1*RMOI(2)
         ELSE IF(RMOI(3).GT.RMOI(2) .OR. RMOI(3).EQ.RMOI(2))THEN
            WRITE(dend,7021)cm1*RMOI(3)
         END IF
      ELSE IF(RMOI(2).EQ.z .AND. RMOI(1).NE.z .AND. RMOI(3).NE.z)THEN
         IF(RMOI(1).GT.RMOI(3))THEN
            WRITE(dend,7021)cm1*RMOI(1)
         ELSE IF(RMOI(3).GT.RMOI(1) .OR. RMOI(3).EQ.RMOI(1))THEN
            WRITE(dend,7021)cm1*RMOI(3)
         END IF
      ELSE IF(RMOI(3).EQ.z .AND. RMOI(1).NE.z .AND. RMOI(2).NE.z)THEN
         IF(RMOI(1).GT.RMOI(2))THEN
            WRITE(dend,7021)cm1*RMOI(1)
         ELSE IF(RMOI(2).GT.RMOI(1) .OR. RMOI(2).EQ.RMOI(1))THEN
            WRITE(dend,7021)cm1*RMOI(2)
         END IF
      ELSE IF(RMOI(3).NE.z .AND. RMOI(1).EQ.z .AND. RMOI(2).EQ.z)THEN
         WRITE(dend,7021)cm1*RMOI(3)
      ELSE IF(RMOI(3).EQ.z .AND. RMOI(1).NE.z .AND. RMOI(2).EQ.z)THEN
         WRITE(dend,7021)cm1*RMOI(1)
      ELSE IF(RMOI(3).EQ.z .AND. RMOI(1).EQ.z .AND. RMOI(2).NE.z)THEN
         WRITE(dend,7021)cm1*RMOI(2)
      END IF
!      close(dend)
 7021 FORMAT(F20.10)
!
C
 800  RETURN
C
C---- Format statements
C
 9911 FORMAT(/,10X,'No Moment of Inertia evaluation is done since ',/,
     &       10X,'there is only one nucleus: NNUC = 1',/)
C
      END SUBROUTINE MOIDRIV
!*==movew.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
!*********************************************************************
      SUBROUTINE MOVEW(NFT,NTH,NALM,NPFLG,NFT1)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE GLOBAL_UTILS, ONLY : CWBOPN
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: MEIG=20
C
C Dummy arguments
C
      INTEGER :: NALM, NFT, NFT1, NPFLG, NTH
      INTENT (IN) NPFLG
      INTENT (OUT) NALM
      INTENT (INOUT) NTH
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp), DIMENSION(MEIG) :: EIG
      INTEGER :: I, IDUM, J, K, M, NEIG, NKEEP, NSTAT, NT
      CHARACTER(LEN=120) :: NAME
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(20) :: NHE
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     MOVEW LOCATES THE NTH-DATASET CONTAINED ON THE ALCHEMY CI
C           DUMPFILE ASSOCIATED WITH UNIT NFT
C
C           NFT   =       UNIT ON WHICH CI VECTORS ARE HELD
C           NFT1  =       SUBSIDIARY UNIT USED BY PRTHD
C           NPFLG = 0     NO PRINTOUT DURING SEARCH
C                 = 1     PRINT HEADER LABELS ENCOUNTERED
C           NALM  = 0     NO ERRORS DETECTED
C                 = 1     ERRORS DETECTED
C
C***********************************************************************
C
      NALM=0
      NT=0
C
      CALL CWBOPN(NFT)
C
      IF(NTH.EQ.1)RETURN
      M=NTH-1
      IF(NTH.EQ.0)M=2000
      DO I=1, M
         READ(NFT,END=50)NT, NHD, NAME, NHE, DTNUC
         NKEEP=NHD(8)
         NSTAT=NHD(3)
         NEIG=MIN(MEIG,NSTAT)
         READ(NFT)(IDUM,J=1,NKEEP), (EIG(J),J=1,NEIG)
C
         IF(NPFLG.GT.0)CALL PRTHD(NT,NHD,NAME,NHE,DTNUC,NEIG,EIG,NFT1)
         DO K=1, NSTAT
            READ(NFT)
         END DO
C
      END DO
      RETURN
C
 50   IF(NTH.NE.0)THEN
         NALM=1
      ELSE
         NTH=NT+1
         BACKSPACE NFT
      END IF
      RETURN
      END SUBROUTINE MOVEW
!*==negd11.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE NEGD11(NORB,IRON,IROB,NDI,MN,M)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: M, NORB
      INTEGER, DIMENSION(NORB) :: IROB, IRON
      INTEGER, DIMENSION(*) :: MN, NDI
      INTENT (INOUT) IROB, IRON, M
C
C Local variables
C
      INTEGER :: I, IREPD, IREPS, L, ND, NRP, NZ
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     NEGD11 - COMPUTES THE NEGATIVE OF THE FIRST DETERMINANT OF THE
C              FIRST CSF.
C
C     INPUT DATA :
C           NORB  TOTAL NUMBER OF ORBITALS IN THE WAVEFUNCTION
C            NDI  DETERMINANTS CONSTITUTING THE WAVEFUNCTION
C             MN  THE ORBITAL TABLE
C
C    OUTPUT DATA :
C           IRON  NEGATIVE OF OCCUPATION NUMBER OF ORBITALS IN THE FIRST
C                 DETERMINANT OF THE FIRST CSF
C           IROB  ORBITAL NUMBERS FOR IRON ARRAY ELEMENTS
C              M  SIZE OF THE ARRAYS IRON AND IROB
C
C***********************************************************************
C
C-----ZEROIZE THE ARRAYS FIRST OF ALL
C
      DO I=1, NORB
         IRON(I)=0
         IROB(I)=0
      END DO
C
C-----LOOP OVER THE REPLACED ORBITALS IN THE FIRST DETERMINANTOF THE FIR
C     CSF AND SET OCCUPATION NUMBER TO +1. HOWEVER IF THE NUMBER OF
C     REPLACEMENTS IS ZERO THERE IS NO NEED TO DO THIS AND THE SUBROUTIN
C     EXITED.
C
      NRP=NDI(1)
      IF(NRP.EQ.0)THEN
         M=0
         RETURN
      END IF
      ND=1
      IREPD=2
      NZ=0
      CALL MODRON(IRON,IROB,NZ,ND,NDI(IREPD),NRP,MN,norb)
C
C-----LOOP OVER REPLACEMENT ORBITALS AND SET THEIR OCCUPATION NUMBER TO
C     -1.
C
      ND=-1
      IREPS=NRP+2
      CALL MODRON(IRON,IROB,NZ,ND,NDI(IREPS),NRP,MN,norb)
C
C-----PACK IRON AND IROB BY ELIMINATING ANY ORBITALS WITH ZERO OCCUPATIO
C     NUMBER
C
      M=0
      DO L=1, NZ
         IF(IRON(L).EQ.0)CYCLE
         M=M+1
         IROB(M)=IROB(L)
         IRON(M)=IRON(L)
      END DO
C
      RETURN
C
      END SUBROUTINE NEGD11
!*==ordtst.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE ORDTST(NTGT,ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,NPROPS,
     &                  ISTATE,JSTATE)
C***********************************************************************
C
C     TSTORD - Target STate ORDering by eigen-energy and property value
C              re-indexing
C
C     Input data:
C           NTGT Number of rows in the target state table. This is
C                obviously the number of unique target states.
C          INDEX Workspace area used for rank and indexing information
C         IWRITE Logical unit for the printer
C
C     Input/Output data:
C                 ISPIN  Columns of the unique target state table
C                   ISZ    "
C                  MGVN    "
C                 GUTOT    "
C                 IREFL    "
C                ENERGY    "
C                NPROPS  Number of property values
C                ISTATE  Pointers to the target state table for moments
C                JSTATE   " " "   ""  "   " "    " "   " "   "   " " "
C                INDEXV  Storage space for the index vector
C                IRANKV   " " "   " "   "   "  rank vector
C             IWRKSPACE  Integer workspace for ordering
C             WORKSPACE  R*8 workspace for ordering (same address)
C                IWRITE  Logical unit for the printer
C
C     Linkage:
C
C        INDEXX, IRANK, REORD8, REORD4
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NPROPS, NTGT
      REAL(KIND=wp), DIMENSION(ntgt) :: ENERGY
      INTEGER, DIMENSION(ntgt) :: GUTOT, IREFL, ISPIN, ISZ, MGVN
      INTEGER, DIMENSION(nprops) :: ISTATE, JSTATE
      INTENT (IN) NPROPS
      INTENT (INOUT) ISTATE, JSTATE
C
C Local variables
C
      INTEGER :: I
      INTEGER, DIMENSION(ntgt) :: INDEXV, IRANKV
      INTEGER, DIMENSION(nprops) :: IWRKSPACE
C
C*** End of declarations rewritten by SPAG
C
C---- Simple case NTGT=1 then no ordering required
C
      IF(NTGT.EQ.1)RETURN
c
C---- In the general case we first build an index vector
C     and from this a rank vector
C
      CALL INDEXX(NTGT,ENERGY,INDEXV)
      CALL RANK(NTGT,INDEXV,IRANKV)
C
C...... Now re-define all the state indices for the property values.
C
C       First of all the I-state vector
C
      DO I=1, NPROPS
         IWRKSPACE(I)=ISTATE(I)
      END DO
C
      DO I=1, NPROPS
         ISTATE(I)=IRANKV(IWRKSPACE(I))
      END DO
C
C..... Now the J-state vector
C
      DO I=1, NPROPS
         IWRKSPACE(I)=JSTATE(I)
      END DO
C
      DO I=1, NPROPS
         JSTATE(I)=IRANKV(iwrkspace(I))
      END DO
c
C
C...... Re-order all columns of the table using the index array.
C
      CALL REORD4(NTGT,ISPIN,INDEXV)
      CALL REORD4(NTGT,ISZ,INDEXV)
      CALL REORD4(NTGT,MGVN,INDEXV)
      CALL REORD4(NTGT,GUTOT,INDEXV)
      CALL REORD4(NTGT,IREFL,INDEXV)
      CALL REORD8(NTGT,ENERGY,INDEXV)
C
      RETURN
      END SUBROUTINE ORDTST
!*==outpak.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE OUTPAK(MATRIX,NROW,NCTL,IWRITE)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, NCTL, NROW
      REAL(KIND=wp), DIMENSION(*) :: MATRIX
      INTENT (IN) IWRITE, MATRIX, NCTL, NROW
C
C Local variables
C
      CHARACTER(LEN=4), DIMENSION(3) :: ASA
      INTEGER :: BEGIN, I, K, KCOL, KTOTAL, LAST, NCOL
      CHARACTER(LEN=4) :: BLANK, CTL
      CHARACTER(LEN=8) :: COLUMN
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     OUTPAK PRINTS A SYMMETRIC MATRIX OF REAL NUMBERS STORED IN ROW
C     PACKED LOWER TRIANGULAR FORM IN FORMATTED FORM WITH NUMBERED ROWS
C     AND COLUMNS.
C
C     INPUT DATA:
C         MATRIX PACKED MATRIX
C           NROW NUMBER OF ROWS TO BE OUTPUT
C           NCTL CARRAIGE CONTROL FLAG: 1 FOR SINGLE SPACE
C                                       2 FOR DOUBLE SPACE
C                                       3 FOR TRIPLE SPACE
C         IWRITE LOGICAL UNIT FOR THE PRINTER
C
C     THE MATRIX ELEMENTS ARE ARRANGED IN STORAGE AS FOLLOWS :
C
C         1
C         2   3
C         4   5   6
C         7   8   9  10
C        11  12  13  14  15
C        16  17  18  19  20  21
C        22  23  24  25  26  27  28 ETC.
C
C     OUTPAK IS SET UP TO HANDLE 8 COLUMNS PER PAGE WITH AN 8F15.8 FORMA
C     FOR THE COLUMNS. IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHA
C     FORMATS 1000 AND 2000, AND INITIALISE KCOL WITH THE NEW NUMBER OF
C     COLUMNS.
C
C     THIS ROUTINE WAS OBTAINED FROM THE QUANTUM THEORY PROJECT,
C     UNIV OF FLORIDA IN FORTRAN IV FORM AND UPGRADED TO FORTRAN 77 BY
C     CHARLES GILLAN AT QUEEN'S UNIV. BELFAST
C
C***********************************************************************
C
      DATA COLUMN/'Column'/, BLANK/'    '/, ASA/'    ', '0   ', '-   '/
      DATA KCOL/4/
C
      CTL=BLANK
      IF(NCTL.LE.3 .AND. NCTL.GT.0)CTL=ASA(NCTL)
C
C---- LAST IS THE LAST COLUMN IN THE ROW CURRENTLY BEING PRINTED
C
      LAST=MIN(NROW,KCOL)
C
C---- BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTE
C
      NCOL=1
      BEGIN=1
C
C---- START A NON STANDARD DO LOOP OVER THE COLUMNS
C
 1050 NCOL=1
      WRITE(IWRITE,1000)(COLUMN,I,I=BEGIN,LAST)
      DO K=BEGIN, NROW
         KTOTAL=(K*(K-1))/2+BEGIN-1
         WRITE(IWRITE,2000)CTL, K, (MATRIX(I+KTOTAL),I=1,NCOL)
         IF(K.LT.(BEGIN+KCOL-1))NCOL=NCOL+1
      END DO
      LAST=MIN(LAST+KCOL,NROW)
      BEGIN=BEGIN+NCOL
      IF(BEGIN.LE.NROW)GO TO 1050
C
      RETURN
C
 1000 FORMAT(/13X,7(3X,A8,I2,1X),(2X,A8,I2))
 2000 FORMAT(A4,' Row',I3,1X,8F14.7)
C
      END SUBROUTINE OUTPAK
!*==output.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE OUTPUT(MATRIX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,
     &                  NCTL,IWRITE)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: COLDIM, COLHI, COLLOW, IWRITE, NCTL, ROWDIM, ROWHI, 
     &           ROWLOW
      REAL(KIND=wp), DIMENSION(ROWDIM,COLDIM) :: MATRIX
      INTENT (IN) COLDIM, COLHI, COLLOW, IWRITE, MATRIX, NCTL, ROWDIM, 
     &            ROWHI, ROWLOW
C
C Local variables
C
      CHARACTER(LEN=8), DIMENSION(3) :: ASA
      INTEGER :: BEGIN, I, K, KCOL, LAST
      CHARACTER(LEN=8) :: BLANK, COLUMN, CTL
      REAL(KIND=wp) :: ZERO
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     OUTPUT PRINTS A MATRIX OF REAL NUMBERS IN FORMATTED FORM WITH
C     NUMBERED ROWS AND COLUMNS.
C
C     INPUT DATA:
C         MATRIX THE MATRIX TO BE OUTPUT
C         ROWLOW THE ROW AT WHICH THE OUTPUT IS TO BEGIN
C          ROWHI THE ROW AT WHICH OUTPUT IS TO TERMINATE
C         COLLOW THE COLUMN AT WHICH THE OUTPUT IS TO BEGIN
C          COLHI THE COLUMN AT WHICH OUTPUT IS TO TERMINATE
C         ROWDIM ROW DIMENSION OF THE MATRIX
C         COLDIM COLUMN DIMENSION OF THE MATRIX
C           NCTL CARRAIGE CONTROL FLAG: 1 FOR SINGLE SPACE
C                                       2 FOR DOUBLE SPACE
C                                       3 FOR TRIPLE SPACE
C         IWRITE LOGICAL UNIT FOR PRINTER
C
C     THE PROGRAM IS SET UP TO HANDLE 8 COLUMNS PER PAGE WITH A 8F14.7
C     FORMAT FOR THE COLUMNS. IF A DIFERENT NUMBER OF COLUMNS IS REQUIRE
C     CHANGE FORMATS 1000 AND 2000 AND INITIALIZE KCOL WITH THE NEW NUMB
C     OF COLUMNS.
C
C     THIS ROUTINE WAS OBTAINED FROM THE QUANTUM THEORY PROJECT,
C     UNIV OF FLORIDA IN FORTRAN IV FORM AND UPGRADED TO FORTRAN 77 BY
C     CHARLES GILLAN AT QUEEN'S UNIV. BELFAST
C
C***********************************************************************
C
      DATA COLUMN/' Column'/, ASA/'        ', '00000000', '--------'/
      DATA BLANK/'       '/
      DATA KCOL/4/, ZERO/0.0_wp/
C
C---- TEST TO SEE IF ALL OF THE ELEMENTS ARE ZERO. BRANCH TO END OF THE
C     ROUTINE IF THIS IS TRUE
C
c      DO 10 I=ROWLOW,ROWHI
c      DO 20 J=COLLOW,COLHI
c      IF (MATRIX(I,J).NE.ZERO) GOTO 15
c  20  CONTINUE
c  10  CONTINUE
c      WRITE(IWRITE,3000)
c      GOTO 3
C
C---- MATRIX IS TO BE PRINTED.
C
      CONTINUE
      CTL=BLANK
      IF(NCTL.LE.3 .AND. NCTL.GT.0)CTL=ASA(NCTL)
      IF(ROWHI.LT.ROWLOW .OR. COLHI.LT.COLLOW)GO TO 3
C
      LAST=MIN(COLHI,COLLOW+KCOL-1)
      DO BEGIN=COLLOW, COLHI, KCOL
         WRITE(IWRITE,1000)(COLUMN,I,I=BEGIN,LAST)
         DO K=ROWLOW, ROWHI
c      DO 4 I=BEGIN,LAST
c      IF (MATRIX(K,I).NE.ZERO) GOTO 5
c 4    CONTINUE
c      GOTO 1
            WRITE(IWRITE,2000)CTL, K, (MATRIX(K,I),I=BEGIN,LAST)
         END DO
         LAST=MIN(LAST+KCOL,COLHI)
      END DO
C
C---- END OF ROUTINE, RETURN POINT
C
 3    RETURN
C
 1000 FORMAT(/11X,7(3X,A8,I2,1X),(2X,A8,I2))
 2000 FORMAT(A4,'Row',I3,8F14.7)
 3000 FORMAT(3X,'ZERO MATRIX')
C
      END SUBROUTINE OUTPUT
!*==prnexp.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PRNEXP(NFD,IFLAGX)
C***********************************************************************
C
C     PRNEXP - PRINT DENSITY MATRIX EXPRESSIONS
C
C     PRINTS OUT THE FORMULAE THAT HAVE BEEN CALCULATED FOR THE
C     DIAGONAL AND OFF DIAGONAL ELEMENTS OF THE DENSITY MATRIX
C     THIS ROUTINE READS AND WRITES HEADERS, THEN CALLS SEPARATE
C     SUBROUTINES TO READ AND WRITE THE ACTUAL EXPRESSIONS
C
C     INPUT DATA:
C
C            NFD  THE LOGICAL UNIT CONTAINING THE FILE OF DENSITY
C                 MATRIX EXPRESSIONS
C             NR  WORKSPACE ARRAY
C         IFLAGX  FLAG WHICH DETERMINES THE SYMMETRY OF THE FUNCTIONS
C                 =1 WAVEFUNCTIONS HAVE DIFFERENT SYMMETRY
C                 =0 WAVEFUNCTIONS HAVE IDENTICAL SYMMETRY
C
C     OUTPUT DATA: NONE
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IFLAGX, NFD
      INTENT (IN) IFLAGX
C
C Local variables
C
      INTEGER :: I, IFLAG, INUM, LD, LDMII, NCODE
C
C*** End of declarations rewritten by SPAG
C
C---- Write out a header to the output stream
C
      WRITE(6,500)NFD
C
C---- Prepare the dataset of expressions for reading
C
      REWIND NFD
C
C---- Default some variables which may never be used before
C     they are printed
C
      LDMII=0
C
C---- Skip header records on unit NFD
C
      IFLAG=0
C
      IF(IFLAGX.EQ.1)IFLAG=1
      IF(IFLAGX.EQ.0)IFLAG=2
      INUM=6-IFLAG-IFLAG
      DO I=1, INUM
         READ(NFD)
      END DO
C
C---- READ AND PRINT DIAGONAL EXPRESSIONS IF THEY EXIST
C
      WRITE(6,510)
C
      IF(IFLAG.EQ.2)THEN
         READ(NFD)NCODE, LDMII
         WRITE(6,520)NCODE, LDMII
         CALL PTDII(NFD,LDMII)
      ELSE
         WRITE(6,530)IFLAG
      END IF
C
C---- Read and print off-diagonal expressions
C
      WRITE(6,540)
C
      READ(NFD,END=30)NCODE, LD
      WRITE(6,520)NCODE, LD
C
      CALL PTDIJ(NFD,LD)
C
C---- Start shut down processing for this routine
C
 20   CONTINUE
C
      REWIND NFD
C
      WRITE(6,550)
C
      RETURN
C
C---- Following code invoked for no off-diagonal elements
C
 30   CONTINUE
C
      WRITE(6,570)
      GO TO 20
C
C---- Format Statements
C
 500  FORMAT('1',//,5X,25('-'),//,6X,'Density Matrix Formulae',//,5X,
     &       25('-'),//,5X,
     &       'Expressions are read from logical unit number = ',I3,/)
 510  FORMAT(5X,'1. Diagonal Elements',//)
 520  FORMAT(5X,'The value of NCODE is',I5/5X,'The buffer size is',I5)
 530  FORMAT(5X,'There are no Diagonal Elements. IFLAG =',I3)
 540  FORMAT(//,5X,'2. Off-Diagonal Elements',//)
 550  FORMAT(//,5X,'All density matrix expressions have been printed')
 570  FORMAT(5X,' No Off-diagonal expressions exist',/)
C
      END SUBROUTINE PRNEXP
!*==prsdx.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PRSDX(NFTA,IWRITE)
C***********************************************************************
C
C     PRSDX - PRint the Sorted Density matrix eXpressions.
C
C     Input data:
C           NFTA Logical unit for the sorted expressions
C         IWRITE Logical unit for the printer
C             KI buffer for indices and, at another time, diagonal
C                formulae.
C             KJ buffer for indices in off-diagonal formulae.
C            KPQ buffer for density matrix position of off-diag elements
C            EPQ buffer for off-diagonal formulae coefficients
C           LBUF size the workspace arrays KI,KJ,KPQ and EPQ
C
C     Output data:
C                  None
C
C     Linkage:
C
C             PTDII
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      CHARACTER(LEN=1), PARAMETER :: CHARI='I', CHARJ='J'
      CHARACTER(LEN=10), PARAMETER :: CSEQU='SEQUENTIAL'
      CHARACTER(LEN=11), PARAMETER :: CMATED='UNFORMATTED'
C
C Dummy arguments
C
      INTEGER :: IWRITE, NFTA
      INTENT (IN) IWRITE
C
C Local variables
C
      CHARACTER(LEN=10) :: CACCESS
      CHARACTER(LEN=11) :: CFORM
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: EPQ
      INTEGER :: I, IA, IB, IDENT, IPASS, ISYMTYP, ITIME, J, JA, JB, L, 
     &           LBOX, LCOF, LDMII, MEL, MET, MGVNI, MGVNJ, MT, MTT, N, 
     &           NBOX, NCODE, NCSFI, NCSFJ, NELT, NORB, NSRB, NSYM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KI, KJ, KPQ
      CHARACTER(LEN=4), DIMENSION(30) :: NAME
      REAL(KIND=wp) :: PIN, R, S, SZ
      LOGICAL :: ZOP, ZTEST
C
C*** End of declarations rewritten by SPAG
C
C---- Banner header
C
      WRITE(IWRITE,500)NFTA
C
C---- Make sure that the file on unit NFTA has the proper attributes
C     for a file of unsorted density matrix expressions
C
      INQUIRE(UNIT=NFTA,ACCESS=CACCESS,FORM=CFORM,OPENED=ZOP)
C
C...... If the attributes are corret then ZTEST becomes .TRUE.
C
      ZTEST=CACCESS.EQ.CSEQU .AND. CFORM.EQ.CMATED
C
      IF(.NOT.ZTEST)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9905)NFTA, CACCESS, CMATED
         CALL TMTCLOS()
      END IF
C
      IF(.NOT.ZOP)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9906)NFTA
         CALL TMTCLOS()
      END IF
C
C---- Prepare the unit for reading
C
      REWIND NFTA
C
C---- Read first wavefunction header.
C
      READ(NFTA)(NAME(I),I=1,30), MGVNI, S, SZ, R, PIN, NORB, NSRB, 
     &          NCSFI, NELT, NSYM, IDENT, ISYMTYP
C
      WRITE(IWRITE,510)CHARI
      WRITE(IWRITE,515)MGVNI, S, SZ, NCSFI
      WRITE(IWRITE,517)NORB, NSRB, NELT, NSYM, IDENT, ISYMTYP
C
C---- If IDENT=0 then
C        read the diagonal formulae
C     else
C        read the header record for the second wavefunction
C     endif
C
      IF(IDENT.NE.1)THEN
         READ(NFTA)NCODE, LDMII
         ALLOCATE(ki(ldmii),kj(ldmii),kpq(ldmii),epq(ldmii))
         CALL PTDII(NFTA,LDMII)
         NCSFJ=NCSFI
      ELSE
         READ(NFTA)(NAME(I),I=1,30), MGVNJ, S, SZ, R, PIN, NORB, NSRB, 
     &             NCSFJ, NELT, NSYM, ISYMTYP
C
         WRITE(IWRITE,510)CHARJ
         WRITE(IWRITE,515)MGVNJ, S, SZ, NCSFJ
      END IF
C
C---- Read the off-diagonal elements
C
      WRITE(IWRITE,590)
C
      MTT=0
      MET=0
C
C...... This is a while loop implimented in Fortran-77
C
C       While NOT end of file
C             read a box of formulae and print them
C       end while
C
C       Successive iterations begin at line 100
C
      itime=0
C
 100  READ(NFTA,END=200)NCODE, IPASS, NBOX, LBOX, LCOF
      IF(NCODE.EQ.0)GO TO 180
      WRITE(IWRITE,520)NCODE, IPASS, NBOX, LBOX, LCOF
C
C..... Make sure that the size of the buffers is big enough
C
      IF(itime.EQ.0)THEN
         itime=itime+1
         IF(ident.NE.1)THEN
            IF(lbox.GT.ldmii)THEN
               DEALLOCATE(ki,kj,kpq,epq)
               ALLOCATE(ki(lbox),kj(lbox),kpq(lbox),epq(lbox))
            END IF
         ELSE
            ALLOCATE(ki(lbox),kj(lbox),kpq(lbox),epq(lbox))
         END IF
      END IF
C
      DO I=1, NCSFI, LCOF
         DO J=1, NCSFJ, LCOF
            READ(NFTA)IA, IB, JA, JB, MT
            MTT=MTT+MT
            WRITE(IWRITE,522)IA, IB, JA, JB, MT
            DO L=1, MT
               READ(NFTA)MEL, (KI(N),N=1,LBOX), (KJ(N),N=1,LBOX), 
     &                   (KPQ(N),N=1,LBOX), (EPQ(N),N=1,LBOX)
               MET=MET+MEL
               WRITE(IWRITE,524)MEL
               WRITE(IWRITE,526)(KI(N),KJ(N),KPQ(N),EPQ(N),N=1,MEL)
            END DO
         END DO
      END DO
C
C...... Loop back to read another record of formulae
C
      GO TO 100
C
C-----SPECIAL CASE : NO OFF DIAGONAL ELEMENTS
C
 180  WRITE(IWRITE,582)NCODE
      DEALLOCATE(ki,kj,kpq,epq)
      RETURN
C
C-----END OF FILE
C
 200  REWIND NFTA
      WRITE(IWRITE,585)MET, MTT
      DEALLOCATE(ki,kj,kpq,epq)
      RETURN
C
C---- Format Statements
C
 500  FORMAT(///,5X,'Printing of Sorted Density Matrix Expressions',/,
     &       5X,'---------------------------------------------',//,5X,
     &       'Expressions read from unit = ',I3,//)
 510  FORMAT(5X,'Wavefunction ',A,' details: ',/)
 515  FORMAT(7X,'Symmetry quantum number  = ',I10,/,7X,
     &       'Spin quantum number      = ',F10.4,/,7X,
     &       'Z projection of spin     = ',F10.4,/,7X,
     &       'Number of CSFs           = ',I10)
 517  FORMAT(7X,'Total number of orbitals = ',I10,/,7X,
     &       'Total spin orbitals      = ',I10,/,7X,
     &       'Number of electrons      = ',I10,/,7X,
     &       'No. of orbital symms     = ',I10,/,7X,
     &       'IDENT flag/switch        = ',I10,/,7X,
     &       'ISYMTYP - point group    = ',I10,//)
 520  FORMAT(5X,'Ncode =',I3,' Pass Number =',I3,' Number of Boxes =',
     &       I4,' Box Length =',I7,/,5X,'Csf blocking factor = ',I7,/)
 522  FORMAT(5X,'For this box :',/,5X,'I-Csfs bias =',I5,
     &       ' I-Csfs Largest Csf =',I5,/,5X,'J-Csfs bias =',I5,
     &       ' J-Csfs Largest Csf =',I5,/,5X,
     &       'Number of Records for the box =',I5,/)
 524  FORMAT(5X,' Number of Elements =',I5,/)
 526  FORMAT(3(1X,3I5,F10.6))
 582  FORMAT(5X,'NCODE =',I5)
 585  FORMAT(5X,'TOTAL NUMBER OF ELEMENTS =',I10,' TOTAL NUMBER',
     &       ' OF RECORDS =',I10,/)
 590  FORMAT(/5X,'SORTED OFF-DIAGONAL ELEMENTS',//)
C
 9900 FORMAT(/,10X,'**** Error in PRSDX ',/)
 9905 FORMAT(10X,'Logical unit for sorted expressions has wrong',/,10X,
     &       'attributes. ',//,10X,'Logical unit = ',I10,/,10X,
     &       'Access       = ',A,' should be SEQUENTIAL  ',/,10X,
     &       'Record form  = ',A,' should be UNFORMATTED ',/)
 9906 FORMAT(10X,'Logical unit = ',I3,' is not open ! ',/)
C
      END SUBROUTINE PRSDX
!*==prthd.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PRTHD(NSET,NHD,HEAD,NHE,DTNUC,NEIG,EIG,NFT)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER(LEN=120) :: HEAD
      INTEGER :: NEIG, NFT, NSET
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp), DIMENSION(neig) :: EIG
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(20) :: NHE
      INTENT (IN) DTNUC, EIG, HEAD, NEIG, NFT, NHD, NHE, NSET
C
C Local variables
C
      REAL(KIND=wp) :: E0
      INTEGER :: I, ND
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     PRTHD PRINTS THE HEADER RECORDS FROM AN ALCHEMY PERMANENT DATASET
C           FOR STORING CI EIGENVALUES AND EIGENVECTORS
C
C     NSET IS THE SET NUMBER
C     NHD,HEAD,NHE,E0,DTNUC,NEIG,EIG ARE READ BY THE ROUTINE WHICH CALLS
C     THIS ONE
C     NFT IS THE UNIT FOR OUTPUT
C
C***********************************************************************
C
      WRITE(NFT,100)NSET, HEAD
 100  FORMAT(/' SET',I4,4X,A)
      e0=dtnuc(1)
      WRITE(NFT,101)NHD, e0
 101  FORMAT(5X,'ICVC =',I5,4X,'NOCSF=',I5,4X,'NSTAT=',I5,4X,'NSYM =',
     &       I5,4X,'IDFLG=',2I5/5X,'IDEN =',I5,4X,'NKEEP=',I5,4X,
     &       'NNUC =',I5,4X,'LDEN =',I5,4X,'EN   =',F20.10)
      WRITE(NFT,103)(NHE(I),I=1,NHD(4))
 103  FORMAT(5X,'NOB  =',20I5)
      ND=21+NHD(9)
      WRITE(NFT,105)(DTNUC(I),I=22,ND)
 105  FORMAT(5X,'GEONUC =',10F10.4/(13X,10F10.4))
      ND=1+NHD(9)
      WRITE(NFT,104)(DTNUC(I),I=2,ND)
 104  FORMAT(5X,'CHARG=',10F10.4/(13X,10F10.4))
      WRITE(NFT,102)(EIG(I)+E0,I=1,NEIG)
 102  FORMAT(5X,'EIGEN-ENERGIES',/(3X,5F20.10))
C
      RETURN
      END SUBROUTINE PRTHD
!*==prtmom.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PRTMOM(NFTMT,IWRITE,NMSET,TMTFORM,KEYTMT)
C***********************************************************************
C
C     Takes the file of moment expressions which has been written by the
C     run and writes these out to the printer for checking
C
C     Linkage:
C
C        UNPAK9
C
C     Author: Charles J Gillan, March 1994
C
C***********************************************************************
C
      USE integer_packing, ONLY:UNPAK9
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, KEYTMT, NFTMT, NMSET
      CHARACTER(LEN=11) :: TMTFORM
      INTENT (IN) IWRITE
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(20) :: CHARG
      REAL(KIND=wp) :: EDIFF, XMOM
      REAL(KIND=wp), DIMENSION(3,20) :: GEONUC
      INTEGER :: GUTOT1, GUTOT2, I, IDENT, IOUTW1, IOUTW2, IREFL1, 
     &           IREFL2, ISET, ISPIN, IST, ISTI, J, JST, JSTJ, K, KEY, 
     &           MGVN1, MGVN2, NCSF1, NCSF2, NELT, NMOMS, NNUC, NOPREC, 
     &           NORB, NOVECI, NOVECJ, NPT, NREC, NSRB, NSTREC, NSYM, 
     &           NUCCEN, SYMTYP
      INTEGER, DIMENSION(8) :: IOPCDE
      CHARACTER(LEN=4), DIMENSION(30) :: NAME
      INTEGER, DIMENSION(20) :: NOB
C
C*** End of declarations rewritten by SPAG
C
C---- Banner header
C
      WRITE(IWRITE,1000)
C
C---- Position the moment file for reading the required set
C
      REWIND NFTMT
      CALL GETSET(NFTMT,NMSET,KEYTMT,TMTFORM)
C
C---- The set number is ISET
C
      READ(NFTMT,END=990,ERR=995)KEY, ISET, nrec, NOPREC, NOVECI, 
     &                           NOVECJ, IDENT, NUCCEN
C
      WRITE(IWRITE,1010)ISET, IDENT, NOVECI, NOVECJ, NUCCEN, NOPREC
C
C     NSTREC is the number of pairs of states with NOPREC moments per
C     pair. Thus the total number of moments per set is
C
C                 NSTREC*(NOPREC+1)
C
      IF(IDENT.NE.1)THEN
         NSTREC=(NOVECI*(NOVECI+1))/2
      ELSE
         NSTREC=NOVECI*NOVECJ
      END IF
C
      READ(NFTMT)NPT
      WRITE(IWRITE,1012)NPT
C
C.... 1. The header records
C
      READ(NFTMT,ERR=995)NORB, NSRB, NELT, NSYM, SYMTYP, ISPIN, 
     &                   (NOB(I),I=1,NSYM), NNUC, 
     &                   ((GEONUC(j,I),j=1,3),I=1,NNUC), 
     &                   (CHARG(I),I=1,NNUC)
 
      READ(NFTMT,ERR=995)(NAME(I),I=1,30), MGVN1, IREFL1, GUTOT1, NCSF1
C
      WRITE(IWRITE,1015)ISPIN, NORB, NSRB, NELT, NSYM, (NOB(I),I=1,NSYM)
      WRITE(IWRITE,1016)NNUC
      WRITE(IWRITE,1017)(I,(GEONUC(j,I),j=1,3),CHARG(I),I=1,NNUC)
C
      IF(IDENT.NE.1)THEN
         WRITE(IWRITE,1031)
      ELSE
         WRITE(IWRITE,1032)
      END IF
C
      WRITE(IWRITE,1040)(NAME(I),I=1,20), MGVN1, NCSF1
C
      IF(IDENT.EQ.0)THEN
         MGVN2=MGVN1
      ELSE
         READ(NFTMT,ERR=995)(NAME(I),I=1,30), MGVN2, IREFL2, GUTOT2, 
     &                      NCSF2
         WRITE(IWRITE,1033)
         WRITE(IWRITE,1040)(NAME(I),I=1,20), MGVN2, NCSF2
      END IF
C
C.... 2. The moment records themselves
C
      DO I=1, NSTREC
         READ(NFTMT,ERR=995)ISTI, JSTJ, EDIFF, nmoms
         WRITE(IWRITE,1055)ISTI, JSTJ, EDIFF
         DO J=1, nmoms
            READ(NFTMT,ERR=995)IOUTW2, IOUTW1, XMOM
c      WRITE(IWRITE,1060)  J,IOUTW2,IOUTW1,XMOM
            CALL UNPAK9(IOPCDE,MGVN1,MGVN2,IST,JST,IOUTW1,IOUTW2)
            WRITE(IWRITE,35)(IOPCDE(K),K=1,7), MGVN1, MGVN2, XMOM
         END DO
      END DO
C
C---- Subroutine return point
C
      WRITE(IWRITE,8000)
C
      RETURN
C
C---- Error condition handlers
C
C.... (a) End of file on read
C
 990  CONTINUE
C
      WRITE(IWRITE,9990)NFTMT
      CALL TMTCLOS()
C
C.... (b) Error during read of NFTMT
C
 995  CONTINUE
C
      WRITE(IWRITE,9995)NFTMT
      CALL TMTCLOS()
C
C---- Format Statements
C
 35   FORMAT(5X,7I2,I3,1X,I3,D20.13)
C
 1000 FORMAT(//,5X,'PRTMOM - Printing of Moment Expressions',/)
 1010 FORMAT(/,5X,'Moments file set number: ',I3,//,5X,
     &       'CSF mixing flag  (IDENT) = ',I5,/,5X,
     &       'No. of I states (NOVECI) = ',I5,/,5X,
     &       'No. of J states (NOVECJ) = ',I5,/,5X,
     &       'Nuclear Center  (NUCCEN) = ',I5,/,5X,
     &       'No. Properties  (NOPREC) = ',I5,/)
 1012 FORMAT(/,5X,'No. of property operators = ',I5,/)
 1015 FORMAT(/,5X,'Information from header record : ',//,5X,
     &       'Total Spin          (ISPIN) = ',I5,/,5X,
     &       'No. of orbitals      (NORB) = ',I5,/,5X,
     &       'No. of spin-orbitals (NSRB) = ',I5,/,5X,
     &       'Number of electrons  (NELT) = ',I5,/,5X,
     &       'No. of orbital syms  (NSYM) = ',I5,/,5X,
     &       'Orbitals per symmetry (NOB) : ',//,5X,20I4)
 1016 FORMAT(/,5X,'Number of nuclei (NNUC) = ',I5,/)
 1017 FORMAT(5X,' No.   Co-ordinates    Charge ',/,
     &       (5X,I3,5X,3F10.4,6X,F5.2))
 1031 FORMAT(/,5X,'Wavefunction details :',/)
 1032 FORMAT(/,5X,'Wavefunction I details :',/)
 1033 FORMAT(/,5X,'Wavefunction J details :',/)
 1040 FORMAT(/,5X,20A,/,5X,'Lambda value of state = ',I3,/,5X,
     &       'Number of CSFs        = ',I5)
 1055 FORMAT(/,5X,'I-state number =',I5,' J-state number =',I5,/,5X,
     &       'Energy difference (Hartree) =',D13.6,/)
 1060 FORMAT(5X,'Rec Number = ',I3,/,5X,'Word 1     = ',Z16,
     &       ' packed hex format '/,5X,'Word 2     = ',Z16,
     &       ' packed hex format '/,5X,'Moment     = ',D20.13)
 8000 FORMAT(/,5X,'**** Printing of Moment Expressions completed',/)
C
 9990 FORMAT(//,5X,'END OF FILE ON READ OF UNIT ',I3,//)
 9995 FORMAT(//,5X,'ERROR ON READ FOR UNIT ',I3)
C
      END SUBROUTINE PRTMOM
!*==prtptab.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PRTPTAB(ISTATE,JSTATE,IPCDE,PROPS,NPROPS,IWRITE)
C***********************************************************************
C
C     PRTPTAB - PRinTs the Properties TABle
C
C     Input data:
C         ISTATE  Designation of state I - the bra vector
C         JSTATE  Designation of state J - the ket vector
C                 where both above are wrt the unique target state table
C          IPCDE  Slater property integral codes (8 per property)
C          PROPS  Expectation values for the wavefunction pair.
C         NPROPS  No. of property elements in the table
C         IWRITE  Logical unit for the printer
C
C     Example:
C
C         ISTATE   JSTATE         IPCDE         PROPS
C         ------   ------        -------        -----
C            1        1      2 0 0 0 0 0 0 0     1.0
C            2        1      2 0 0 0 0 0 0 0     0.0
C            2        2      2 0 0 0 0 0 0 0     1.0
C
C     where property 2 0 0 0 0 0 0 0 is the overlap, i.e.
C
C                 < State 1 | State 2 >  =  0.0
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, NPROPS
      INTEGER, DIMENSION(8,nprops) :: IPCDE
      INTEGER, DIMENSION(nprops) :: ISTATE, JSTATE
      REAL(KIND=wp), DIMENSION(nprops) :: PROPS
      INTENT (IN) IPCDE, ISTATE, IWRITE, JSTATE, NPROPS, PROPS
C
C Local variables
C
      INTEGER :: I, J
C
C*** End of declarations rewritten by SPAG
C
      WRITE(IWRITE,1000)
      WRITE(IWRITE,1500)NPROPS
 
C
C---- Loop over the rows of the table's columns and print them.
C
      WRITE(IWRITE,2000)
C
      DO I=1, NPROPS
         WRITE(IWRITE,2010)I, ISTATE(I), JSTATE(I), (IPCDE(J,I),J=1,8), 
     &                     PROPS(I)
      END DO
C
 
      RETURN
C
C---- Format statements
C
 1000 FORMAT(//,15X,'========================',//,15X,
     &       'W A V E F U N C T I O N ',/,15X,
     &       '  P R O P E R T I E S   ',/,15X,
     &       '       T A B L E        ',//,15X,
     &       '========================',//)
 1500 FORMAT(//5x,'No. of Expectation Values (Rows) in table = ',I7)
 2000 FORMAT(/'     No.',1X,'I-State',1X,'J-State',1X,
     &       'Property Operator Code',1X,'Expectation Value',/,5X,'---',
     &       1X,'-------',1X,'-------',1X,'----------------------',1X,
     &       '-----------------',/)
 2010 FORMAT(1X,I5,1X,I5,1X,I5,1X,3X,8I2,3X,1X,D11.4)
C
      END SUBROUTINE PRTPTAB
!*==prttst.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PRTTST(NTGT,ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,ISYMTYP,
     &                  IWRITE)
C***********************************************************************
C
C     PRTTST - PRinT the Target STate table
C
C     Input data:
C           NTGT Number of rows in the target state table. This is
C                obviously the number of unique target states.
C          ISPIN 2*S+1 for each target state
C            ISZ 2*Sz for each target state
C           MGVN Lamda value (C-inf-v) or Irred. Rep of each state
C          GUTOT For D-inf-h only, the gerade or ungerade value
C          IREFL For C-inf-v only the +/- sigma character
C         ENERGY Eigen-energy in Hartrees for the state
C        ISYMTYP Switch for C-inf-v or Abelian point groups
C         IWRITE Logical unit for the printer
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ISYMTYP, IWRITE, NTGT
      REAL(KIND=wp), DIMENSION(ntgt) :: ENERGY
      INTEGER, DIMENSION(ntgt) :: GUTOT, IREFL, ISPIN, ISZ, MGVN
      INTENT (IN) ENERGY, GUTOT, IREFL, ISPIN, ISYMTYP, ISZ, IWRITE, 
     &            MGVN, NTGT
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
      WRITE(IWRITE,1000)NTGT
C
      IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
C
C=======================================================================
C
C     C-inf-v and D-inf-h type molecules
C
C=======================================================================
C
         WRITE(IWRITE,1005)
C
         DO I=1, NTGT
            WRITE(IWRITE,1010)I, ISPIN(I), ISZ(I), MGVN(I), GUTOT(I), 
     &                        IREFL(I), ENERGY(I)
         END DO
C
      ELSE
C
C=======================================================================
C
C     Abelian type molecules
C
C=======================================================================
C
         WRITE(IWRITE,2005)
C
         DO I=1, NTGT
            WRITE(IWRITE,2010)I, ISPIN(I), ISZ(I), MGVN(I), ENERGY(I)
         END DO
C
      END IF
c
      RETURN
C
C---- Format Statements
C
 1000 FORMAT(//,10X,'Unique Target State Table',/,10X,
     &       '-------------------------',//,10X,'There are ',I5,
     &       ' states in the table',/)
C
 1005 FORMAT(5X,'(2S+1)   Sz    Lamda   G/U     +/-   Energy (Hartrees)'
     &       ,/,5X,
     &       '------   --    -----   ---     ---   -----------------',/)
 1010 FORMAT(1X,I4,5(I6,1X),F17.10)
C
 2005 FORMAT(5X,'(2S+1)   Sz    Irrep   Energy (Hartrees)',/,5X,
     &       '------   --    -----   -----------------',/)
 2010 FORMAT(1X,I4,3(I6,1X),F17.10)
C
      END SUBROUTINE PRTTST
!*==ptdii.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PTDII(NFD,LDMII)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LDMII, NFD
      INTENT (IN) LDMII, NFD
C
C Local variables
C
      INTEGER :: II, J, K, N, NEL, NN
      INTEGER, DIMENSION(ldmii) :: MP
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     PTDII - PRINT DIAGONAL ELEMENTS OF DENSITY EXPRESSIONS
C
C     READS IN THE EXPRESSIONS FOR THE DIAGONAL ELEMENTS OF THE DENSITY
C     MATRIX AND PRINTS THESE OUT. END OF THE LIST IS SIGNALLED BY A
C     NEGATIVE VALUE OF NEL.
C
C     INPUT DATA:
C            NFD  LOGICAL UNIT CONTAINING THE EXPRESSIONS
C          LDMII  BUFFER SIZE. LDMII+1 IS THE MAXIMUM SIZE OF THE
C                 RECORD ON FILE NFD.
C             MP  BUFFER WHICH IS USED FOR THE INPUT AND OUTPUT
C
C     OUTPUT DATA: NONE
C
C***********************************************************************
c
C
C---- Header for the table following
C
      WRITE(6,500)
C
      II=0
 110  READ(NFD)NEL, MP
      N=ABS(NEL)
      K=1
 112  II=II+1
      NN=MP(K)*2
      IF(NN.NE.0)THEN
         WRITE(6,510)II, MP(K), (MP(K+J),J=1,NN)
      ELSE
         WRITE(6,510)II, MP(K)
      END IF
      K=K+NN+1
C
      IF(K.LE.N)GO TO 112
      IF(NEL.GT.0)GO TO 110
C
      RETURN
C
C---- Format Statements
C
 500  FORMAT(/,5X,'Diagonal density expressions are stored as the ',
     &       'difference',/,5X,
     &       'between D(II) and D(11). Of course D(11) is stored',/,5X,
     &       'in absolute form.',//,5X,'Csf',3X,'No. of Elements',3X,
     &       'Orbital Occupation No.'/)
 510  FORMAT(5X,I3,3X,I15,3X,I7,1X,I14,/,(29X,I7,1X,I14))
C
      END SUBROUTINE PTDII
!*==ptdij.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PTDIJ(NFD,LDMIJ)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LDMIJ, NFD
      INTENT (IN) LDMIJ, NFD
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(ldmij) :: CPQ
      INTEGER, DIMENSION(LDMIJ) :: II, JJ
      INTEGER, DIMENSION(2,LDMIJ) :: MPQ
      INTEGER :: N, NEL, NN
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     PTDIJ - PRINT THE OFF-DIAGONAL DENSITY EXPRESSIONS
C
C     READS IN THE EXPRESSIONS FOR THE OFF-DIAGONAL DENSITY MATRIX
C     ELEMENTS AND PRINTS THESE OUT. THE LAST RECORD IS SIGNALLED BY
C     A NEGATIVE VALUE OF NEL, UNLESS NEL= 0
C
C     INPUT DATA:
C            NFD  LOGICAL UNIT CONTAINING THE EXPRESSIONS
C          LDMIJ  SIZE OF THE BUFFERS II,JJ,MPQ AND CPQ
C             II  BUFFER
C             JJ    "
C            MPQ    "
C            CPQ    "
C
C     OUTPUT DATA: NONE
C
C***********************************************************************
c
C
      WRITE(6,500)
C
C---- Following is a while not end of file loop rading over the file
C     on unit NFD. Successive iterations begin at line 20.
C
 20   CONTINUE
C
C...... Read a buffer of expressions
C
      READ(NFD)NEL, II, JJ, MPQ, CPQ
C
      IF(NEL.NE.0)THEN
         NN=ABS(NEL)
         WRITE(6,510)(II(N),JJ(N),MPQ(1,N),MPQ(2,N),CPQ(N),N=1,NN)
      ELSE
         WRITE(6,520)
         STOP
      END IF
C
C..... Return to read another buffer unless this is the final one.
C      The last buffer is signified by NEL being negative.
C
      IF(NEL.GT.0)GO TO 20
C
      RETURN
C
C---- Format statements
C
 500  FORMAT(/,5X,'D(ij) EXPRESSIONS',//,5X,'I',5X,'J',5X,'P',5X,'Q',7X,
     &       'C(pq)')
 510  FORMAT(2(1X,4(I5,1X),F11.5))
 520  FORMAT(5X,'NUMBER OF ELEMENTS = 0')
C
      END SUBROUTINE PTDIJ
!*==ptpwf.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE PTPWF(NFTW,NOCSF,NELT,NDTRF,NODI,INDI,ICDI,NDI,CDI)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NELT, NFTW, NOCSF
      REAL(KIND=wp), DIMENSION(*) :: CDI
      INTEGER, DIMENSION(nocsf) :: ICDI, INDI, NODI
      INTEGER, DIMENSION(*) :: NDI
      INTEGER, DIMENSION(nelt) :: NDTRF
      INTENT (IN) CDI, ICDI, INDI, NDI, NDTRF, NELT, NFTW, NOCSF, NODI
C
C Local variables
C
      INTEGER :: I, K, MA, MB, MC, MD, N
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     PTPWF - PRINTS OUT DETERMINANTAL WAVEFUNCTIONS
C
C             THIS SUBROUTINE HAS BEEN TAKEN UNMODIFIED FORM THE CRAY
C             VERSION OF THE ALCHEMY PACKAGE. IT PRINTS OUT THE SLATER
C             DETERMINANTS, CONSTITUTING THE WAVEFUNCTION, IN A MORE
C             READABLE FORM.
C
C     INPUT DATA :
C           NFTW  LOGICAL UNIT FOR PRINTER
C          NOCSF  NUMBER OF CSFS IN THE WAVEFUNCTION
C           NELT  NUMBER OR ELECTRONS IN THE WAVEFUNCTION
C          NDTRF  REFERENCE DETERMINANT
C           NODI  NUMBER OF DETERMINANTS PER CSF
C           INDI  POSITION OF FIRST DETERMINANT OF EACH CSF IN ARRAY NDI
C           ICDI  POSITION OF FIRST COEFFICIENT OF EACH CSF IN ARRAY CDI
C            NDI  DETERMINANT CODES
C            CDI  COEFFICIENTS FOR EACH DETERMINANT
C
C***********************************************************************
C
      WRITE(NFTW,139)(NDTRF(I),I=1,NELT)
 139  FORMAT(' REFERENCE DETERMINANT'//(1X,15I5))
      WRITE(NFTW,137)
 137  FORMAT('  CSF',9X,'COEFFICIENT',2X,'NSO'/)
c
      DO N=1, NOCSF
         MA=NODI(N)
         MB=INDI(N)
         MC=ICDI(N)-1
         MD=NDI(MB)
         WRITE(NFTW,138)N, CDI(MC+1), MD, (NDI(MB+I),I=1,2*MD)
 138     FORMAT(1X,I4,D20.10,I5,2X,15I5/(32X,15I5))
         MB=MB+MD+MD+1
         DO K=2, MA
            MD=NDI(MB)
            WRITE(NFTW,140)CDI(MC+K), MD, (NDI(MB+I),I=1,2*MD)
 140        FORMAT(5X,D20.10,I5,2X,15I5/(32X,15I5))
            MB=MB+MD+MD+1
         END DO
      END DO
c
      RETURN
      END SUBROUTINE PTPWF
!*==rank.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE RANK(N,INDX,IRANK)
C***********************************************************************
C
C     RANK - Takes an array of indices as output by the routine
C            INDEXX and returns a table of ranks.
C
C     Input data:
C              N number of elements in the array to be ordered
C           INDX array of indices created by routine INDEXX.
C
C     Output data:
C           IRANK a set ranks corresponding to the indices
C
C     Notes:
C
C     This routine is taken from the book Numerical Receipes by
C     Press, Flannery, Teukolsky and Vetterling Chapter 8 p. 233.
C     ISBN 0-521-30811-9 pub. Cambridge University Press (1986)
C     QA297.N866
C
C     This routines has been adapted by Charles J Gillan for use
C     in the R-matrix codes.
C
C***********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: N
      INTEGER, DIMENSION(n) :: INDX, IRANK
      INTENT (IN) INDX, N
      INTENT (OUT) IRANK
C
C Local variables
C
      INTEGER :: J
C
C*** End of declarations rewritten by SPAG
C
C
C
C---- Convert the indices into ranks
C
      DO J=1, N
         IRANK(INDX(J))=J
      END DO
C
      RETURN
C
      END SUBROUTINE RANK
!*==rdsped.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE RDSPED(NFT,nftd,IWRITE,coefs,nodeti,icoefp,idetp,idets,
     &                  NDTRF,NOCSF,NOB,LCDOF,LNDOF,NPFLG)
C***********************************************************************
C
C     RDSPED  READS THE WAVEFUNCTION PRODUCED BY THE ALCHEMY FORMULA
C             GENERATOR SPEEDY. THE WAVEFUNCTION CONSISTS OF THE
C             REFERENCE DETERMINANT AND THE CSFS DEFINED IN TERMS OF
C             REPLACEMENTS FROM IT. IT IS UNBLOCKED.
C
C             THE DATA IS ACTUALLY READ FROM THE FILE PRODUCED BY SORT !
C             IT IS NOT PUT TOGETHER IN A SUITABLE FORM, IN THE OUTPUT
C             FROM SPEEDY, TO BE READ STRAIGHT IN HERE. SORT AMALGAMATES
C             THE NECESSARY PIECES FROM THE SPEEDY FILE READ BY IT.
C
C     INPUT DATA:
C            NFT LOGICAL UNIT CONTAINING THE SORTED FORMULAE
C         IWRITE LOGICAL UNIT FOR THE PRINTER
C
C     OUTPUT DATA:
C           LCDOF LENGTH OF THE COEFFICIENT LIST
C           LNDOF LENGTH OF THE DETERMINANT LIST
C            NSYM NUMBER OF SYMMETRIES IN THE WAVEFUNCTION
C           NPFLG PRINT FLAG
C             NOB NUMBER OF ORBITALS PER SYMMETRY
C            NELT NUMBER OF ELECTRONS IN THE WAVEFUNCTION
C           NOCSF NUMBER OD CSFS IN THE WAVEFUNCTION
C            NORB TOTAL NUMBER OF ORBITALS
C            NSRB TOTAL NUMBER OF SPIN ORBITAS IN THE WAVEFUNCTION
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, LCDOF, LNDOF, NFT, NFTD, NOCSF, NPFLG
      REAL(KIND=wp), DIMENSION(lcdof) :: COEFS
      INTEGER, DIMENSION(nocsf) :: ICOEFP, IDETP, NODETI
      INTEGER, DIMENSION(lndof) :: IDETS
      INTEGER, DIMENSION(500) :: NDTRF
      INTEGER, DIMENSION(20) :: NOB
      INTENT (IN) NFTD, NPFLG
      INTENT (OUT) NODETI
      INTENT (INOUT) LCDOF, LNDOF, NDTRF, NOB, NOCSF
C
C Local variables
C
      INTEGER :: GUTOT, I, IDIAG, IREFL, LNDI, LTRI, MGVN, NCTARG, NELT, 
     &           NORB, NPFLH, NSRB, NSYM, NTGSYM, SYMTYP
      CHARACTER(LEN=120) :: NAME
      REAL(KIND=wp) :: PIN, R, S, SZ, THRES
C
C*** End of declarations rewritten by SPAG
C
C     ------- READING SECTION AND TRANSFER TO NFTD ------------
C
      REWIND NFT
C
C     READ HEADER AND SECOND RECORD AND TRANSFER TO NFTD
C     (adjusted to take file from CONGEN as well as SORT)
C
      READ(NFT,END=50,ERR=50)NAME, MGVN, S, SZ, R, PIN, NORB, NSRB, 
     &                       NOCSF, NELT, LTRI, IDIAG, NSYM, SYMTYP, 
     &                       lndi, npflh, thres, nctarg, ntgsym
c     CONGEN-style header: skip record with phase information
      READ(nft)
      GO TO 60
c     This was a SORT-style header: read again with shorter record
 50   BACKSPACE nft
      READ(NFT)NAME, MGVN, S, SZ, R, PIN, NORB, NSRB, NOCSF, NELT, LTRI, 
     &         IDIAG, NSYM, SYMTYP
 60   CONTINUE
      READ(NFT)(NOB(I),I=1,NSYM), (NDTRF(I),I=1,NELT), 
     &         (nodeti(I),I=1,NOCSF)
c
      WRITE(NFTD)NAME, MGVN, S, SZ, R, PIN, NORB, NSRB, NOCSF, NELT, 
     &           NSYM, SYMTYP
      WRITE(NFTD)(NOB(I),I=1,NSYM)
C
C---- COMPUTE THE REFLECTION SYMMETRY AND G/U SYMMETRY FROM THE
C     INFORMATION READ IN.
C
      CALL REFLGU(GUTOT,IREFL,R,PIN)
C
C     READ THIRD RECORD
C
      READ(NFT)(ICOEFP(I),I=1,NOCSF), LCDOF, (IDETP(I),I=1,NOCSF), LNDOF
C
C     READ FOURTH RECORD
C
      LNDOF=LNDOF-1
      LCDOF=LCDOF-1
C
      READ(NFT)(IDETS(I),I=1,LNDOF)
      CALL RDSPX(NFT,COEFS,LCDOF)
C
C     ------- WAVEFUNCTION OUTPUT SECTION ------------
C
      IF(NPFLG.NE.0)THEN
         WRITE(iwrite,10)NAME, MGVN, S, SZ, R, PIN
         WRITE(iwrite,12)NORB, NSRB, NOCSF, NELT
         WRITE(iwrite,14)LTRI, IDIAG, NSYM, SYMTYP
         WRITE(iwrite,16)(NOB(I),I=1,NSYM)
         WRITE(iwrite,18)(NDTRF(I),I=1,NELT)
         CALL WFNDMP(ICOEFP,IDETP,COEFS,IDETS,LCDOF,LNDOF,NOCSF,IWRITE)
      END IF
C
      RETURN
C
C---- Format Statements
C
 10   FORMAT(///,5X,25('-'),//6X,'WAVEFUNCTION INPUT',//,5X,25('-'),
     &       //5X,30A4,//,5X,'MGVN =',I3,' S = ',F3.1,' SZ = ',F3.1,
     &       ' R = ',F3.1,' PIN = ',F4.1)
 12   FORMAT(5X,'NORB =',I4,' NSRB =',I4,' NOCSF =',I5,' NELT =',I4)
 14   FORMAT(5X,'LTRI =',I5,' IDIAG =',I5,' NSYM =',I3,' SYMTYP =',I2)
 16   FORMAT(5X,'NOB =',(20(I5,1X)))
 18   FORMAT(/5X,'Reference Determinant',//,(5X,20(I5,1X)))
C
      END SUBROUTINE RDSPED
!*==rdsped1.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE RDSPED1(NFT,NOCSF,NORB,NSRB,NELT,LCDOF,LNDOF,NSYM,MGVN,
     &                   S,SZ,R,PIN,NAME,SYMTYP)
C***********************************************************************
C
C     RDSPED  READS THE WAVEFUNCTION PRODUCED BY THE ALCHEMY FORMULA
C             GENERATOR SPEEDY. THE WAVEFUNCTION CONSISTS OF THE
C     INPUT DATA:
C            NFT LOGICAL UNIT CONTAINING THE SORTED FORMULAE
C         IWRITE LOGICAL UNIT FOR THE PRINTER
C
C     OUTPUT DATA:
C           LCDOF LENGTH OF THE COEFFICIENT LIST
C           LNDOF LENGTH OF THE DETERMINANT LIST
C            NSYM NUMBER OF SYMMETRIES IN THE WAVEFUNCTION
C            NELT NUMBER OF ELECTRONS IN THE WAVEFUNCTION
C           NOCSF NUMBER OD CSFS IN THE WAVEFUNCTION
C            NORB TOTAL NUMBER OF ORBITALS
C            NSRB TOTAL NUMBER OF SPIN ORBITAS IN THE WAVEFUNCTION
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE GLOBAL_UTILS, ONLY : CWBOPN
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LCDOF, LNDOF, MGVN, NELT, NFT, NOCSF, NORB, NSRB, NSYM, 
     &           SYMTYP
      CHARACTER(LEN=120) :: NAME
      REAL(KIND=wp) :: PIN, R, S, SZ
      INTENT (OUT) LCDOF, LNDOF, MGVN, NAME, NELT, NORB, NSRB, NSYM, 
     &             PIN, R, S, SYMTYP, SZ
      INTENT (INOUT) NOCSF
C
C Local variables
C
      INTEGER :: I, ICOEFP, IDETP, IDIAG, LNDI, LTRI, NCTARG, NPFLH, 
     &           NTGSYM
      REAL(KIND=wp) :: THRES
C
C*** End of declarations rewritten by SPAG
C
      CALL cwbopn(nft)
C
C     READ HEADER AND SECOND RECORD
C     (adjusted to take file from CONGEN as well as SORT)
C
      READ(NFT,END=50,ERR=50)NAME, MGVN, S, SZ, R, PIN, NORB, NSRB, 
     &                       NOCSF, NELT, LTRI, IDIAG, NSYM, SYMTYP, 
     &                       lndi, npflh, thres, nctarg, ntgsym
 
c     CONGEN-style header: skip record with phase information
      READ(nft)
      GO TO 60
c
c     This was a SORT-style header: read again with shorter record
 50   BACKSPACE nft
      READ(NFT)NAME, MGVN, S, SZ, R, PIN, NORB, NSRB, NOCSF, NELT, LTRI, 
     &         IDIAG, NSYM, SYMTYP
 60   CONTINUE
C
      READ(NFT)
      READ(NFT)(ICOEFP,I=1,NOCSF), LCDOF, (IDETP,I=1,NOCSF), LNDOF
C
      RETURN
      END SUBROUTINE RDSPED1
!*==rdspx.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE RDSPX(NFT,CR,LC)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LC, NFT
      REAL(KIND=wp), DIMENSION(LC) :: CR
      INTENT (IN) LC, NFT
      INTENT (OUT) CR
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     RDSPX - READS THE COEFFICIENTS ACCOMPANYING EACH SLATER
C             DETERMINANT INTO CORE
C
C***********************************************************************
      READ(NFT)CR
      RETURN
      END SUBROUTINE RDSPX
!*==reflgu.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE REFLGU(IGU,IREFL,R,PIN)
C***********************************************************************
C
C     REFLGU - Computes the G/U symmetry and reflection symmetries
C              for C-inf-v/D-inf-h point groups from the R and PIN
C              numbers.
C
C     Input data:
C              R R*8 packed reflection symmetry information
C            PIN R*8 packed G/U information
C
C     Output data:
C             IGU G or U flag if appropriate:
C                 =+1 means gerade symmetry of D-inf-h
C                  -1 means ungerade symmetry of D-inf-h
C                   0 the molecule, of the CSFs, are in C-inf-v
C           IREFL + or - reflection symmetry for sigma wavefunctions
C                 =+1 means + reflection
C                  -1 means - reflection
C                   0 symmetry is pi or higher i.e. no reflection
C
C     Notes:
C
C       Care has been taken to make this machine portable
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL, PARAMETER :: POINT01=0.01E+00
C
C Dummy arguments
C
      INTEGER :: IGU, IREFL
      REAL(KIND=wp) :: PIN, R
      INTENT (IN) PIN, R
      INTENT (OUT) IGU, IREFL
C
C Local variables
C
      INTEGER :: INT
      REAL :: PINTMP, RTEMP
      REAL :: REAL
C
C*** End of declarations rewritten by SPAG
C
      RTEMP=REAL(R)
      PINTMP=REAL(PIN)
C
      IREFL=INT(RTEMP+SIGN(POINT01,RTEMP))
      IGU=INT(PINTMP+SIGN(POINT01,PINTMP))
C
      RETURN
C
      END SUBROUTINE REFLGU
!*==reord8.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE REORD8(N,ARRIN,INDX)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: N
      REAL(KIND=wp), DIMENSION(n) :: ARRIN
      INTEGER, DIMENSION(n) :: INDX
      INTENT (IN) INDX, N
      INTENT (INOUT) ARRIN
C
C Local variables
C
      INTEGER :: J
      REAL(KIND=wp), DIMENSION(n) :: WORKSPACE
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     REORD8  - Given an array and a set of index vectors this
C               routine places the array in the order specified
C               by the index vector. The vector is treated as
C               as a real*8 data type.
C
C     Input data:
C              N number of elements in the array to be ordered
C           INDX array of indices created by routine INDEXX.
C      WORKSPACE workspace array used for a temporary copy of ARRIN
C
C     Input/Output data:
C                 ARRIN On input this is not ordered but is overwritten
C                       by the ordered array on output.
C
C***********************************************************************
C
C
C---- Make a copy of the input array ARRIN into the WORKSPACE
      DO J=1, N
         WORKSPACE(J)=ARRIN(J)
      END DO
C
C---- Create the vector ARRIN in order specified by array INDEX
      DO J=1, N
         ARRIN(J)=WORKSPACE(INDX(J))
      END DO
C
      RETURN
      END SUBROUTINE REORD8
!*==reord4.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE REORD4(N,IARRIN,INDX)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: N
      INTEGER, DIMENSION(n) :: IARRIN, INDX
      INTENT (IN) INDX, N
      INTENT (INOUT) IARRIN
C
C Local variables
C
      INTEGER, DIMENSION(n) :: IWRKSPACE
      INTEGER :: J
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     REORD4  - Given an array and a set of index vectors this
C               routine places the array in the order specified
C               by the index vector. The vector is treated as
C               as an integer data type.
C
C     Input data:
C              N number of elements in the array to be ordered
C           INDX array of indices created by routine INDEXX.
C      IWRKSPACE workspace array used for a temporary copy of ARRIN
C
C     Input/Output data:
C                IARRIN On input this is not ordered but is overwritten
C                       by the ordered array on output.
C
C***********************************************************************
C
C
C---- Make a copy of the input array ARRIN into the WORKSPACE
      DO J=1, N
         IWRKSPACE(J)=IARRIN(J)
      END DO
C
C---- Create the vector ARRIN in order specified by array INDEX
      DO J=1, N
         IARRIN(J)=IWRKSPACE(INDX(J))
      END DO
C
      END SUBROUTINE REORD4
!*==rwdii.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE RWDII(NFD,NFT,LDMII,MN,MM,MBAS)
C***********************************************************************
C
C     RWDII - Read and sort the diagonal density matrix expressions
C
C     By sorting it is meant that one changes the orbital label
C     into a position in the array which holds density matrix.
C
C     INPUT DATA:
C            NFD  LOGICAL UNIT CONTAINING THE DENSITY MATRIX
C                 EXPRESSIONS
C            NFT  LOGICAL UNIT TO HOLD THE OUTPUT EXPRESSIONS
C          LMDII  LENGTH OF BUFFER USED FOR INPUT
C             MN  THE ORBITAL TABLE SYMMETRY BY SYMMETRY
C             MM  M-VALUES FOR ALL OF THE ORBITALS
C           MBAS  NUMBER OF ORBITALS FOR EACH SYMMETRY (C-INF-V)
C
C     OUTPUT DATA:
C                 THERE IS NO OUTPUT DATA AS SUCH. SORTED DENSITY
C                 EXPRESSIONS ARE WRITTEN TO A DATASET.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LDMII, NFD, NFT
      INTEGER, DIMENSION(*) :: MBAS, MM, MN
      INTENT (IN) LDMII, MBAS, MM, MN, NFD, NFT
C
C Local variables
C
      INTEGER :: I, K, KP, M, MP, N, ND, NDD, NI, NL, NP
      INTEGER, DIMENSION(ldmii) :: LBUF
C
C*** End of declarations rewritten by SPAG
C
C-----INITIALIZE THE COUNTERS FOR NUMBER OF ELEMENTS, M, AND
C     NUMBER OF RECORDS WRITTEN TO NEW FILE, K.
C
      M=0
      K=0
C
C-----READ ONE RECORD FROM THE INPUT FILE, NFD, AND BEGIN
C     PROCESSING IT. CHANGE THE ORBITAL NUMBERS TO POSITIONS IN THE
C     DENSITY MATRIX ARRAY.
C
 10   READ(NFD)N, LBUF
      NDD=ABS(N)
      ND=1
C
 12   NI=ND+1
      NL=ND+2*LBUF(ND)
C
      IF(LBUF(ND).NE.0)THEN
         DO I=NI, NL, 2
            NP=LBUF(I)
            MP=MM(NP)+1
            KP=MN(NP)
            LBUF(I)=(KP*(KP+1))/2+MBAS(MP)
         END DO
      END IF
      ND=NL+1
C
C-----CONTINUE TO PROCESS THIS RECORD IF NOT ALREADY DONE.
C     OTHERWISE WRITE NEW RECORD AND AUGMENT COUNTERS.
C
      IF(ND.LE.NDD)GO TO 12
      WRITE(NFT)N, LBUF
      M=M+ABS(N)
      K=K+1
C
C-----IF THIS IS NOT THE LAST RECORD GO BACK AND READ ANOTHER FROM
C     THE FILE. THE LAST RECORD IS SIGNALLED BY N BEING ZERO OR
C     NEGATIVE.
C
      IF(N.GT.0)GO TO 10
C
C-----WRITE SUMMARY INFORMATION TO OUTPUT STREAM
C     FOR THE ENTIRE PROCESS
C
C      WRITE(6,500)M,K,LDMII
C
      RETURN
C
C---- Format Statements
C
 500  FORMAT(//,5X,'Statistics from the sorting of D(II)',/,5X,
     &       ' (On diagonal density expressions)  ',/,5X,
     &       '------------------------------------',//,5X,
     &       'Number of coefficients written = ',I7,/,5X,
     &       'Number of records written      = ',I7,/,5X,
     &       'Buffer size (integer words)    = ',I7,//)
C
      END SUBROUTINE RWDII
!*==rwdij.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE RWDIJ(NOCSF,LBLK,NFTI,LDM,NFTX,NBOX,LBOX,NFTO,NBOXJ,
     &                 NOB,NSM,NRI,MBAS,IBKPTL,IBKPTR,NBLOCKS,MDEL)
C***********************************************************************
C
C     RWDIJ - Sorts off-diagonal density matrix expressions. Orbital
C             indices are converted into a sequential index and
C             expressions sorted into boxes defined by the partitioning
C             of the CI vectors.
C
C     INPUT :
C       NOCSF NUMBER OF CSFS PER WAVEFUNCTION (AN ARRAY SIZE 2)
C       LBLK  BLOCK SIZE FOR THE CSF ARRAY USED FOR PARTITIONING
C             THE CSF ARRAY INTO BOXES.
C       NFTI  UNIT CONTAINING THE DENSITY MATRIXC ELEMENTS TO BE
C             SORTED
C       LDM   SIZE OF BUFFER FOR THE OFF DIAGONAL ELEMENTS
C       NI    BUFFERS FOR CSF NUMBERS
C       NJ    BUFFERS FOR CSF NUMBERS
C       NPQ   AN ARRAY WHICH WILL HOLD THE SEQUENTIAL INDEX FOR EACH
C             PAIR OF ORBITALS
C       CPQ   BUFFER FOR COEFFICIENT WITH EACH ELEMENT
C       NFTX  LOGICAL UNIT NUMBER OF THE DIRECT ACCESS FILE USED IN THE
C             SORTING STEP
C       NBOX  NUMBER OF BOXES IN THE SORTING STEP
C       LBOX  SIZE OF A BOX
C       NCHNA USED IN THE SORTING STEP TO BACK CHAIN DIRECT ACCESS
C             RECORDS FOR EACH BOX
C       NCHNB NUMBER OF ELEMENTS IN EACH BOX
C       NCHNC
C       MI    BOX FOR CSF INDEX USED IN SORT STEP
C       MJ     "   "   "    "     "   "  "    "
C       MPQ   BOX FOR SEQUENTIAL ORBITAL INDEX USED IN SORT STEP
C       DPQ    "   "  COEFFICIENTS USED IN SORT STEP
C       NFTO  SEQUENTIAL FILE TO WHICH SORTED EXPRESSIONS ARE COPIED
C       KI    OUTPUT BUFFER FOR TRANSFER TO SEQUENTIAL FILE
C       KJ      "      "     "     "     "     "         "
C       KPQ     "      "     "     "     "     "         "
C       EPQ     "      "     "     "     "     "         "
C       NBOXJ J-BOX SIZE
C       NPQ2  BUFFER FOR THE INPUT ORBITAL INDICES
C
C       THIS SUBROUTINE RETURNS NO VALUES AS SUCH BUT CREATES A FILE OF
C       SORTED DENSITY EXPRESSIONS WHICH IS PASSED ON TO THE ACTUAL
C       DENSITY MATRIX CONSTRUCTION ROUTINES.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : ITHREE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LBLK, LBOX, LDM, MDEL, NBLOCKS, NBOX, NBOXJ, NFTI, 
     &           NFTO, NFTX
      INTEGER, DIMENSION(*) :: IBKPTL, IBKPTR, MBAS, NOB, NRI, NSM
      INTEGER, DIMENSION(2) :: NOCSF
      INTENT (IN) LBLK, LDM, NBOX, NBOXJ, NFTI, NOCSF
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(ldm) :: CPQ
      REAL(KIND=wp), DIMENSION(lbox,nbox) :: DPQ
      REAL(KIND=wp), DIMENSION(lbox) :: EPQ
      INTEGER :: I, IA, IB, IBX, IPASS, J, JA, JB, JDSK, K, KDSK, M, 
     &           MEL, MT, MTT, N, NCSFI, NCSFJ, NEL, NN, NT, NTEL
      INTEGER, DIMENSION(lbox) :: KI, KJ, KPQ
      INTEGER, DIMENSION(lbox,nbox) :: MI, MJ, MPQ
      INTEGER, DIMENSION(nbox) :: NCHNA, NCHNB, NCHNC
      INTEGER, DIMENSION(ldm) :: NI, NJ, NPQ
      INTEGER, DIMENSION(2,ldm) :: NPQ2
C
C*** End of declarations rewritten by SPAG
C
C----- INITIALIZATION OF VARIABLES
C
      IPASS=0
      NCSFI=NOCSF(1)
      NCSFJ=NOCSF(2)
      MI=-1
      MJ=-1
      MPQ=-1
      CPQ=-1
      DPQ=-1
      EPQ=-1
C
C----- PASSES OVER THE FILE OF MOMENT EXPRESSIONS START HERE
C
 100  CONTINUE
C
      JDSK=0
      IPASS=IPASS+1
      NTEL=0
      MTT=0
      DO N=1, NBOX
         NCHNA(N)=0
         NCHNB(N)=0
         NCHNC(N)=0
      END DO
C
c      IF (IPASS .EQ. 1) WRITE(6,510)
C
C----- READ MOMENT EXPRESSIONS AND SORT THEM
C
 120  READ(NFTI)NEL, NI, NJ, NPQ2, CPQ
C
C      IF THERE ARE NO MOMENT EXPRESSIONS THEN SKIP SORT
C
      IF(NEL.EQ.0)THEN
         WRITE(NFTO)NEL, IPASS, NBOX, LBOX, LBLK
         WRITE(6,500)
         RETURN
      END IF
C
      NT=ABS(NEL)
C
C----- CONVERT PACKED ORBITAL LABELS INTO A SEQUENTIAL INDEX
C
 
      CALL INDEX1(NT,NPQ,NPQ2,NOB,NSM,NRI,MBAS,IBKPTr,IBKPTl,NBLOCKS,
     &            MDEL)
 
C
C----- LOOP OVER THE EXPRESSIONS AND DROP INTO APPROPRIATE BOX.
C      UPON FILLING UP A BOX WRITE IT OUT TO THE DIRECT ACCESS FILE
C
      DO N=1, NT
         IA=(NI(N)-1)/LBLK+1
         IBX=NBOXJ*(IA-1)+(NJ(N)-1)/LBLK+1
         M=NCHNB(IBX)+1
         NCHNB(IBX)=M
         MI(M,IBX)=NI(N)
         MJ(M,IBX)=NJ(N)
         MPQ(M,IBX)=NPQ(N)
         DPQ(M,IBX)=CPQ(N)
         IF(M.GE.LBOX)THEN
            JDSK=JDSK+1
            CALL WRTDA(NFTX,JDSK,NCHNA(IBX),LBOX,M,MI(1,IBX),MJ(1,IBX),
     &                 MPQ(1,IBX),DPQ(1,IBX))
            NCHNA(IBX)=JDSK
            NCHNB(IBX)=0
            NCHNC(IBX)=NCHNC(IBX)+1
         END IF
      END DO
C
C----- HAVING SORTED ONE RECORD OF EXPRESSIONS, IS THIS THE LAST ONE
C      OR ARE THERE MORE ? IF THERE ARE MORE, DO WE NEED TO START A
C      NEW PASS OVER THE FILE OF MOMENT EXPRESSIONS ? IF SPACE IS STILL
C      AVAILABLE ON THE DA FILE THEN GO BACK AND READ MORE RECORDS OF
C      MOMENT EXPRESSIONS.
C
      IF(NEL.GT.0)GO TO 120
C
C----- IF THERE IS NO MORE DA SPACE FOR THIS PASS OR IF WE HAVE REACHED
C      THE LAST RECORD OF MOMENT EXPRESSIONS, THEN WRITE OUT D(IJ) TO A
C      SEQUENTIAL FILE.
C
      WRITE(NFTO)ITHREE, IPASS, NBOX, LBOX, LBLK
c
      N=0
      DO I=1, NCSFI, LBLK
         IA=I-1
         IB=MIN(I+LBLK-1,NCSFI)
         DO J=1, NCSFJ, LBLK
            JA=J-1
            JB=MIN(J+LBLK-1,NCSFJ)
            N=N+1
            MT=NCHNC(N)
            IF(NCHNB(N).NE.0)MT=MT+1
            WRITE(NFTO)I, IB, J, JB, MT
            MTT=MTT+MT
c
            IF(NCHNB(N).NE.0)THEN
               NN=NCHNB(N)
               DO K=1, NN
                  MI(K,N)=MI(K,N)-IA
                  MJ(K,N)=MJ(K,N)-JA
               END DO
               CALL WRTDB(NFTO,LBOX,NN,MI(1,N),MJ(1,N),MPQ(1,N),DPQ(1,N)
     &                    )
               NTEL=NTEL+NN
            END IF
c
 
 150        IF(NCHNA(N).EQ.0)CYCLE
            KDSK=NCHNA(N)
            READ(NFTX,REC=KDSK)NCHNA(N), MEL, KI, KJ, KPQ, EPQ
            DO K=1, MEL
               KI(K)=KI(K)-IA
               KJ(K)=KJ(K)-JA
            END DO
            WRITE(NFTO)MEL, KI, KJ, KPQ, EPQ
            NTEL=NTEL+MEL
            GO TO 150
         END DO
      END DO
C
C----- END OF THIS PASS. PRINT STATISTICS.
C
C      WRITE(6,550) IPASS,JDSK
C
C----- IS ANOTHER PASS REQUIRED OR HAVE WE REACHED THE END OF THE
C      MOMENT EXPRESSIONS ?
C
      IF(NEL.GT.0)GO TO 100
c      WRITE(6,600) NTEL,MTT,LBOX,IPASS
C
 
      RETURN
C
C---- Format Statements
C
 500  FORMAT(//,5X,'THERE ARE NO OFF-DIAGONAL ELEMENTS. THE SORT ',
     &       'HAS BEEN SKIPPED',/)
 510  FORMAT(/,5X,'Sorting transition moment expressions:')
 550  FORMAT(/,5X,'PASS NUMBER',I3,3X,'DIRECT ACCESS RECORDS USED =',I7)
 600  FORMAT(/,5X,'THE SORTING PROCEDURE HAS BEEN COMPLETED. ',/,5X,
     &       'THE FOLLOWING STATISTICS SUMMARIZE THE PROCESS ;',/,5X,
     &       'COEFFICIENTS WRITTEN =',I10,/,5X,'RECORDS WRITTEN =',I10,
     &       /,5X,'LBOX =',I10,/,5X,
     &       'TOTAL NUMBER OF PASSES OVER FILE =',I10,//)
C
      END SUBROUTINE RWDIJ
!*==setd11.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE SETD11(NRP,NDI,JROB,JRON,MN,NC,M,norb)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: M, NC, NORB, NRP
      INTEGER, DIMENSION(norb) :: JROB, JRON
      INTEGER, DIMENSION(*) :: MN, NDI
      INTENT (IN) MN, NDI, NORB, NRP
      INTENT (INOUT) JROB, JRON, M, NC
C
C Local variables
C
      INTEGER :: I, J, L, N, ND
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     LOOKS THROUGH THE FIRST DETERMINANT OF THE FIRST CSF AND USING
C     THE REPLACEMENTS FROM THE REFERENCE DETERMINANT, MODIFIES JRON
C     AND JROB ACCORDINGLY, THUS SETTING UP D11
C
C     INPUT :
C
C       NRP  NUMBER OF REPLACEMENTS IN 1ST DETERMINANT OF 1ST CSF
C       NDI  THE PACKED DETERMINANT OF 1ST CSF
C      JROB  LIST OF ORBITALS IN REFERENCE DETERMINANT
C      JRON  OCCUPATION OF THE ORBITALS IN REFERENCE DETERMINANT
C        MN  THE ORBITAL TABLE
C        NC  NUMBER OF OCCUPIED ORBITALS IN THE REFERENCE DET.
C
C     OUTPUT :
C       JROB  LIST OF ORBITALS IN THE 1ST DETERMINANT OF THE 1ST CSF
C       JRON  LIST OF ORBITALS IN THE 1ST DETERMINANT OF THE 1ST CSF
C          M  SIZE OF ARRAYS JROB AND JRON
C
C***********************************************************************
C
C-----LOOP OVER REPLACEMENTS FROM REF DET. INCLUDED WITHIN LOOP ALSO ARE
C     THE ORBITALS WHICH REPLACE THESE. THIS SECTION IS SKIPPED IF NRP I
C     ZERO UNDER FORTRAN 77
C
      DO I=1, NRP
         ND=-1
         N=NDI(I+1)
         M=MN(N)
 130     DO J=1, NC
            IF(JROB(J).NE.M)CYCLE
            JRON(J)=JRON(J)+ND
            GO TO 160
         END DO
C
         NC=NC+1
         JROB(NC)=M
         JRON(NC)=ND
 160     IF(ND.LE.0)THEN
            N=NDI(I+NRP+1)
            M=MN(N)
            ND=1
            GO TO 130
         END IF
C
      END DO
C
C-----PACK D11 IE. REMOVE ANY ORBITALS WITH ZERO OCCUPATION NUMBER.
C     M COUNTS THE SIZE OF THE FINAL ARRAY.
C
      M=0
      DO L=1, NC
         IF(JRON(L).NE.0)THEN
            M=M+1
            JROB(M)=JROB(L)
            JRON(M)=JRON(L)
         END IF
      END DO
C
      RETURN
      END SUBROUTINE SETD11
!*==setmba.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE SETMBA(MDEL,NSYM,NOB,MBAS,NBKPTL,NBKPTR,ISYMTYP,
     &                  NBLOCKS,LDM,iprnt)
C***********************************************************************
C
C     SETMBA - SETS UP THE POINTER ARRAY MBAS
C
C              GENERATES AN ARRAY OF POINTERS TO THE FIRST STORAGE
C              LOCATION IN THE DENSITY MATRIX ARRAY, FOR EACH OVERLAP
C              ORBITAL PAIR.
C
C     INPUT DATA :
C           MDEL  ABSOLUTE VALUE OF THE LAMBDA DIFFERENCE BETWEEN THE
C                 TWO WAVEFUNCTIONS.
C           NSYM  NUMBER OF C-INF-V SYMMETRIES IN THE ORBITAL SET.
C            NOB  NUMBER OF ORBITALS IN EACH C-INF-V SYMMETRY.
C
C     OUTPUT DATA :
C            MBAS  ARRAY OF POINTERS.
C             LDM  number of elements in the density matrix
C
C     Linkage:
C
C        MPROD, TMTCLOS
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE global_utils, ONLY : MPROD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IPRNT, ISYMTYP, LDM, MDEL, NBLOCKS, NSYM
      INTEGER, DIMENSION(nsym+1) :: MBAS
      INTEGER, DIMENSION(nsym) :: NBKPTL, NBKPTR, NOB
      INTENT (IN) IPRNT, ISYMTYP, MDEL, NOB, NSYM
      INTENT (OUT) LDM
      INTENT (INOUT) MBAS, NBKPTL, NBKPTR, NBLOCKS
C
C Local variables
C
      CHARACTER(LEN=5), DIMENSION(4) :: CC2V
      CHARACTER(LEN=5), DIMENSION(8) :: CD2H
      CHARACTER(LEN=5), DIMENSION(20) :: CINFV
      INTEGER :: I, IBAS, ICOUNT, ISYM, JJ, JSYM, L, L2, LAST, M, N, 
     &           NBGN
C
C*** End of declarations rewritten by SPAG
C
C---- Data statements
C
C...... Names of Irreducible representations for C-inf-v
C
C       Alchemy I and II are dimensioned to deal with up to 20 of these
C
      DATA(CINFV(I),I=1,20)/'SIGMA', 'PI   ', 'DELTA', 'PHI  ', 'GAMMA', 
     &     'H    ', 'I    ', 'J    ', 'K    ', 'L    ', 'M    ', 
     &     'N    ', 'O    ', 'P    ', 'Q    ', 'R    ', 'S    ', 
     &     'T    ', 'U    ', 'V    '/
C
C...... Names of Irreducible representations for D2h
C
      DATA(CD2H(I),I=1,8)/' AG  ', ' B3U ', ' B2U ', ' B1G ', ' B1U ', 
     &     ' B2G ', ' B3G ', ' AU  '/
      DATA(Cc2v(I),I=1,4)/' A1  ', ' B1  ', ' B2  ', ' A2 '/
C
C---- Data initialization
C
      NBLOCKS=0
      MBAS(1)=0
      LDM=0
C
C---- Depending upon the value of ISYMTYP branch to an appropriate
C     in-line subroutine and evaluate the number of blocks in the
C     density matrix and their sizes.
C
C     We may also trap an error in ISYMTYP at this point.
C
      IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
C
C=======================================================================
C
C     for C-infinity-v molecules
C
C=======================================================================
C
C---- Set up of pointers depends upon the symmetry difference of the
C     two wavefunctions
C
         IF(MDEL.EQ.0)THEN
            DO L=1, NSYM
               MBAS(L+1)=MBAS(L)+(NOB(L)*(NOB(L)+1))/2
               NBKPTL(L)=L
               NBKPTR(L)=L
            END DO
            NBLOCKS=NSYM
            LDM=MBAS(NSYM+1)
         ELSE IF(MDEL.EQ.1)THEN
            DO L=1, NSYM
               L2=L+1
               IF(l2.LE.nsym)THEN
                  IBAS=NOB(L)*NOB(L2)
               ELSE
                  ibas=0
               END IF
               MBAS(L2)=MBAS(L)+IBAS
               NBKPTL(L)=L
               NBKPTR(L)=L2
            END DO
C         NBLOCKS=NSYM-MDEL
            NBLOCKS=NSYM-1
            LDM=MBAS(NSYM)
         END IF
C
         IF(MDEL.GE.2)THEN
            ICOUNT=1
            MBAS(1)=0
            NBGN=(MDEL+1)/2+1
            DO JJ=NBGN, NSYM
               N=JJ
               ICOUNT=ICOUNT+1
               M=ABS(MDEL+1-JJ)+1
               IF(N.EQ.M)THEN
                  MBAS(ICOUNT)=MBAS(ICOUNT-1)+NOB(M)*(NOB(M)+1)/2
               ELSE
                  MBAS(ICOUNT)=MBAS(ICOUNT-1)+NOB(M)*NOB(N)
               END IF
               NBKPTL(ICOUNT-1)=m
               NBKPTR(ICOUNT-1)=n
C
            END DO
            NBLOCKS=ICOUNT-1
            LDM=MBAS(ICOUNT)
C
         END IF
c
      ELSE
C
C=======================================================================
C
C     Abelian point group molecules (D2h and subgroups)
C
C=======================================================================
C
C---- If the product of the two wavefunctions transforms as the totally
C     symmetric representation, i.e. MDEL=1, then wavefunctions have
C     the same symmetry. This means that we have NSYM blocks of
C     lower half triangles in the density matrix.
C
         IF(MDEL.EQ.0)THEN
            DO L=1, NSYM
               MBAS(L+1)=MBAS(L)+(NOB(L)*(NOB(L)+1))/2
               NBKPTL(L)=L
               NBKPTR(L)=L
            END DO
         ELSE
C
            last=0
            DO isym=1, nsym
               mbas(isym)=last
               jsym=MPROD(isym,MDEL+1,0,6)
               ! watch out - NOB may be shorter than max sym if it ends with zeros
               if (jsym <= nsym) then
                  last=last+nob(isym)*nob(jsym)
               end if
               NBKPTL(isym)=isym
               NBKPTR(isym)=jsym
            end do
            mbas(nsym+1)=last
         END IF
         NBLOCKS=NSYM
         LDM=MBAS(NSYM+1)
      END IF
c
C
C..... Loop over all blocks of the density matrix and print their
C      specifications
C
      IF(iprnt.NE.0)THEN
         WRITE(6,8010)ISYMTYP, NBLOCKS
C
         IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
            DO I=1, NBLOCKS
               WRITE(6,8020)I, CINFV(NBKPTL(I)), CINFV(NBKPTR(I)), 
     &                      MBAS(I)
            END DO
         ELSE
            IF(nsym.EQ.4)THEN
               DO I=1, NBLOCKS
                  WRITE(6,8020)I, Cc2v(NBKPTL(I)), Cc2v(NBKPTR(I)), 
     &                         MBAS(I)
               END DO
            ELSE
               DO I=1, NBLOCKS
                  WRITE(6,8020)I, CD2H(NBKPTL(I)), CD2H(NBKPTR(I)), 
     &                         MBAS(I)
               END DO
            END IF
         END IF
      END IF
C
      RETURN
C
C---- Format Statements
C
 8010 FORMAT(/,5X,'Density Matrix Morphology:',/,5X,
     &       '--------------------------',//,5X,
     &       'Value of ISYMTYP defining point group  = ',I2,/,5X,
     &       'Number of blocks in the density matrix = ',I2,//,5X,'No.',
     &       1X,'Sym 1',1x,'Sym 2',1X,'Pointer'/,5X,'---',1X,'-----',1x,
     &       '-----',1X,'-------',/)
 8020 FORMAT(5X,I3,1X,A,1x,A,4x,I7)
C
      END SUBROUTINE SETMBA
!*==setndj.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE SETNDJ(NOCSF,JDETS,JDETPN,NODETJ,LNDOFN)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LNDOFN, NOCSF
      INTEGER, DIMENSION(*) :: JDETPN, JDETS, NODETJ
      INTENT (IN) JDETS, NOCSF, NODETJ
      INTENT (OUT) JDETPN, LNDOFN
C
C Local variables
C
      INTEGER :: I, IST, J, NDETS, NRP
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     SETNDJ - SET NEW DETERMINANTS FOR WAVEFUNCTION J
C
C              THIS ROUTINE TAKES THE NEW ARRAY OF DETERMINANT CODES FOR
C              WAVEFUNCTION J W.R.T. REFERENCE DETERMINANT I AND SETS UP
C              A POINTER ARRAY GIVING THE START OF THE FIRST DETERMINANT
C              FOR EACH CSF. THIS IS REQUIRED FOR FORMULA GENERATION
C              SUBROUTINES.
C
C     INPUT DATA :
C          NOCSF  NUMBER OF CSFS IN WAVEFUNCTION J
C          JDETS  NEW DETERMINANT CODES FOR WAVEFUNCTION J W.R.T. THE
C                 REFERENCE DETERMINANT FOR WAVEFUNCTION I.
C         NODETJ  THE NUMBER OF DETERMINANTS PER CSF IN WAVEFUNCTION J
C
C     OUTPUT DATA :
C         JDETPN   THE POSITION OF THE FIRST DETERMINANT FOR EACH CSF
C                  IN THE ARRAY JDETS.
C         LNDOFN   THE SIZE OF THE ARRAY HOLDING THE NEW DETERMINANTS
C                  FOR WAVEFUNCTION J.
C
C***********************************************************************
C
      IST=1
      DO J=1, NOCSF
         NDETS=NODETJ(J)
         JDETPN(J)=IST
         DO I=1, NDETS
            NRP=JDETS(IST)
            IST=IST+(2*NRP)+1
         END DO
      END DO
C
      LNDOFN=IST-1
C
      RETURN
C
      END SUBROUTINE SETNDJ
!*==sortdx.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE SORTDX(NFTA,NFTD,NFDA,NPFLG,LCOF,IDENT,IWRITE,LBOX)
C***********************************************************************
C
C     SORTDX - SORTS the density matrix expessions from being defined
C              as pairs of orbitals into positions in the density matrix
C              array. Obviously this parallels the program SORT used to
C              order Hamiltonian elements form SPEEDY.
C
C     Input data:
C           NFTA  Logical unit on which sorted expressions are written
C           NFTD   " " "   ""  from which unsorted expression are read
C           NFDA   " " "   ""  for the Fortran Direct Access file
C          NPFLG  Print flag
C           LCOF  Box sort size for Hamiltonian eigen vectors
C          IDENT  moment/transition moment switch
C         IWRITE  Logical unit for the printer
C           MDEL  Value of mdel under consideration
C           LBOX
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE global_utils, ONLY : MPROD, CWBOPN
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IDENT, IWRITE, LBOX, LCOF, NFDA, NFTA, NFTD, NPFLG
      INTENT (IN) IDENT
C
C Local variables
C
      INTEGER :: I, ISYMTYP, LDM, LDMII, LDMIJ, MDEL, MGVN1, MGVN2, 
     &           NBLOCKS, NBOX, NBOXI, NBOXJ, NCODE, NCSF1, NCSF2, NELT, 
     &           NOBT, NORB, NSRB, NSYM, NTOT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MBAS, NBKPTL, NBKPTR, NRI, 
     &       NSM
      INTEGER, DIMENSION(2) :: MQN, NOCSF
      CHARACTER(LEN=4), DIMENSION(30) :: NAME
      INTEGER, DIMENSION(20) :: NOB
      REAL(KIND=wp) :: PIN, R, S, SZ
C
C*** End of declarations rewritten by SPAG
C
C---- Write banner header to the output stream
C
      IF(NPFLG.GT.0)WRITE(IWRITE,5000)
C
C---- Prepare the unsorted density matrix expression file, NFTD, for
C     reading from the top and the new scratch file, NFTA, to be written
C     to - it holds sorted density matrix expressions
C
      REWIND NFTD
      CALL CWBOPN(nfta)
C
C---- Read the header records from the density matrix expression dataset
C     and copy them to the dataset of sorted density matrix formulae.
C
C     For diagnostic purposes some details are written to the printer.
C
C...... First wavefunction
C
      READ(NFTD)(NAME(I),I=1,30), MGVN1, S, SZ, R, PIN, NORB, NSRB, 
     &          NCSF1, NELT, NSYM, ISYMTYP
      READ(NFTD)(NOB(I),I=1,NSYM)
C
      WRITE(NFTA)(NAME(I),I=1,30), MGVN1, S, SZ, R, PIN, NORB, NSRB, 
     &           NCSF1, NELT, NSYM, IDENT, ISYMTYP, (NOB(I),I=1,NSYM)
C
      IF(NPFLG.GT.0)THEN
         WRITE(IWRITE,5800)(NAME(I),I=1,30)
         WRITE(IWRITE,5810)MGVN1, S, SZ, NELT, NCSF1, NSYM, 
     &                     (NOB(I),I=1,NSYM)
      END IF
C
C...... Second wavefunction, if applicable !
C
C       Note that if there is no second wavefunction, that is we are
C       dealing with a wavefunction with itself, then we must set the
C       appropriate values into the NOCSF and MQN arrays.
C
      IF(IDENT.EQ.1)THEN
         READ(NFTD)(NAME(I),I=1,30), MGVN2, S, SZ, R, PIN, NORB, NSRB, 
     &             NCSF2, NELT, NSYM, ISYMTYP
         READ(NFTD)(NOB(I),I=1,NSYM)
         WRITE(NFTA)(NAME(I),I=1,30), MGVN2, S, SZ, R, PIN, NORB, NSRB, 
     &              NCSF2, NELT, NSYM, ISYMTYP
         IF(NPFLG.GT.0)THEN
            WRITE(IWRITE,5800)(NAME(I),I=1,30)
            WRITE(IWRITE,5810)MGVN2, S, SZ, NELT, NCSF2, NSYM, 
     &                        (NOB(I),I=1,NSYM)
         END IF
         NOCSF(1)=NCSF1
         NOCSF(2)=NCSF2
         MQN(1)=MGVN1
         MQN(2)=MGVN2
      ELSE
         NOCSF(1)=NCSF1
         NOCSF(2)=NOCSF(1)
         MQN(1)=MGVN1
         MQN(2)=MQN(1)
      END IF
C
C---- Count the total number of orbitals in the problem
C
      NTOT=0
      DO I=1, NSYM
         NTOT=NTOT+NOB(I)
      END DO
C
C---- Evaluate pointers for dynamic storage allocation.
C
C     It is useful to remember that:
C
C     MBAS holds the biases for each symmetry block of the
C               density matrix (c.f. routine SETMBA)
C     NRI  is first column of the orbital table
C     NSM  the second column of the orbital table
C
C     Note that we over estimate the number of density matrix
C     subblocks (2*NSYM) but this wastes little storage and is
C     quicker than an accurate calculation.
C
      ALLOCATE(mbas(nsym+1),nbkptl(nsym),nbkptr(nsym),nri(ntot),
     &         nsm(ntot))
C
C---- Establish the symmetry difference between the wavefunctions.
C
C     For D2h, and subgroups, remember that we must augument MGVN
C     values by one in order to get a correct result from MPROD.
C
C     We may verify that ISYMTYP is sensible here too.
C
      IF(ISYMTYP.EQ.0 .OR. ISYMTYP.EQ.1)THEN
c        WRITE(IWRITE,1092)
         MDEL=ABS(MQN(1)-MQN(2))
      ELSE
c        WRITE(IWRITE,1093)
         mdel=mprod(mqn(1)+1,mqn(2)+1,0,IWRITE)-1
      END IF
C
C..... Tell the user about the wavefunction symmetries and their
C      product
C
      WRITE(IWRITE,1090)MQN(1), MQN(2)
      WRITE(IWRITE,1095)MDEL
C
C---- Establish the bias values for density matrix sub-blocks.
C
C     Error check that MDEL has indeed been set correctly
C
      IF(MDEL.LT.0)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9924)MDEL
         CALL TMTCLOS()
      END IF
C
c ** changed
      CALL SETMBA(MDEL,NSYM,NOB,MBAS,NBKPTL,NBKPTR,ISYMTYP,NBLOCKS,LDM,
     &            npflg)
c
C---- Compute the orbital table for the sorting procedure and then
C     write it to the printer. The unsorted density matrix
C     expressions use contiguous labels for the orbitals, i.e. the
C     label does not change through a symmetry.
C
C     e.g. if NOB=3,2 then orbitals are numbered 1,2,3,4,5,
C
C     After NOBID one has the table
C
C         Orig. label    NRI         NSM
C              1          1           0
C              2          2           0
C              3          3           0
C              4          1           1
C              5          2           1
C
C     which is used to sort expressions into symmetry blocks.
C
      CALL NOBID(NOBT,NSYM,NOB,NRI,NSM,IWRITE,NPFLG)
C
C---- Read back the diagonal elements, if these exist, and sort them.
C
 
      LDMII=0
C
      IF(IDENT.NE.1)THEN
         READ(NFTD)NCODE, LDMII
         WRITE(NFTA)NCODE, LDMII
         CALL RWDII(NFTD,NFTA,LDMII,NRI,NSM,MBAS)
      END IF
C
C---- Sort the off-diagonal elements. Obviously we do not sort
C     off-diagonal elements if none exist - that occurs in the
C     single configuration case which is quantified by the
C     following variable values:
C
C       (IDENT=0 and NOCSF(1)=1) or IDENT=2
C
      IF((IDENT.EQ.0 .AND. NOCSF(1).EQ.1) .OR. ident.EQ.2)GO TO 120
C
C---- Read the code defining that these are off-diagonal elements and
C     the buffer size for them.
C
C     We do not check the buffer size yet because it is part and parcel
C     of the whole storage allocation - that allocation is checked at
C     the end. However we check that NCODE is correct and that LDMIJ
C     is sensible - failure of the test means that there is probably
C     an error either on the data file or in the code's handling of it.
C
      READ(NFTD)NCODE, LDMIJ
C
      IF(NCODE.NE.3)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9920)NCODE
         CALL TMTCLOS()
      END IF
C
      IF(LDMIJ.LE.0)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9930)LDMIJ
         CALL TMTCLOS()
      END IF
C
C---- Off-diagonal sorting - storage management. Note that the box size
C     chosen depends upon the length of the direct access record.
C
C     Further explanation can be had in the routine WRTDA.
C
      NBOXI=(NOCSF(1)-1)/LCOF+1
      NBOXJ=(NOCSF(2)-1)/LCOF+1
      NBOX=NBOXI*NBOXJ
C
C-----PERFORM THE SORT ON OFF-DIAGONAL EXPRESSIONS
C
 
      CALL RWDIJ(NOCSF,LCOF,NFTD,LDMIJ,NFDA,NBOX,LBOX,NFTA,NBOXJ,NOB,
     &           NSM,NRI,MBAS,NBKPTL,NBKPTR,NBLOCKS,MDEL)
C
C---- Finish with all datasets used in the SORTing procedure
C
 120  CONTINUE
C
      DEALLOCATE(mbas,nbkptl,nbkptr,nri,nsm)
      REWIND NFTA
      REWIND NFTD
C
      RETURN
C
C---- Format Statements
C
 1090 FORMAT(//5X,'Wavefunction symmetry numbers = ',I3,1X,I3)
 1092 FORMAT(5X,'Molecule belongs to C-infinity-v point group',/)
 1093 FORMAT(5X,'Molecule belongs to an Abelian point group',/)
 1095 FORMAT(/,5X,'Wavefunction symmetry difference = ',I3)
 5000 FORMAT(//,5X,45('-'),//6X,'SORTING OF MOMENT ','EXPRESSIONS',//,
     &       5X,45('-'),//)
 5800 FORMAT(' NAME = ',30A)
 5810 FORMAT(5X,' Wavefunction Information :',//,5X,
     &       ' Symmetry Quantum Number      (MGVN)  =',I10,/,5X,
     &       ' Spin Quantum Number          (S)     = ',F10.1,/,5X,
     &       ' Z-projection of Spin         (SZ)    = ',F10.1,/,5X,
     &       ' Number of Electrons          (NELT)  =',I10,/,5X,
     &       ' Number of CSFs               (NOCSF) =',I10,/,5X,
     &       ' Number of orbital symmetries (NSYM)  =',I10,//,5X,
     &       ' No. of orbitals per symmetry (NOB)   =',(8(I3,1X),/))
C
 9900 FORMAT(/,10X,'**** ERROR in SORTDX ',/)
 9920 FORMAT(10X,'Value of NCODE, the key for off-diagonal unsorted ',/,
     &       10X,'density matrix expressions is not 3 but = ',I10,/)
 9924 FORMAT(10X,'ABEND - MDEL',i5,' NOT SET BY CODE - BUG',/)
 9930 FORMAT(10X,'Value of LDMIJ, the buffer size for off-diagonal ',/,
     &       10X,'unsorted density matrix expressions is silly = ',I10)
C
      END SUBROUTINE SORTDX
!*==nobid.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE NOBID(NOBT,NSYM,NOB,NRI,NSM,IWRITE,NPFLG)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, NOBT, NPFLG, NSYM
      INTEGER, DIMENSION(*) :: NOB, NRI, NSM
      INTENT (IN) IWRITE, NOB, NPFLG, NSYM
      INTENT (INOUT) NOBT, NRI, NSM
C
C Local variables
C
      INTEGER :: I, ISYM, L, N
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     NOBID - Identify the symmetry of the orbitals in the array NOB
C             and produce an orbital table for use in the sorting
C             step.
C
C     Input data:
C           NSYM  Number of symmetries in array NOB
C            NOB  The number of orbitals in each symmetry
C         IWRITE  Logical unit for the printer
C
C     Output data:
C            NOBT  Size of the arrays NRI and NSM, in other words
C                  the total number of orbitals in the set.
C             NRI  Table of orbital numbers beginning at 1 for each
C                  symmetry in the orbital set.
C             NSM  The symmetry value for each orbital in NRI.
C
C     Example:
C
C     If NSYM=2 and NOB=3,2 then orbitals are numbered 1,2,3,4,5,
C
C     After NOBID one has the table
C
C         Orig. label    NRI         NSM
C              1          1           0
C              2          2           0
C              3          3           0
C              4          1           1
C              5          2           1
C
C     Note: This routine is independent of point group symmetry because
C           it deals with orbitals not spin-orbitals !
C
C***********************************************************************
C
C
C---- Initialize the orbital counter to zero
C
      N=0
C
C---- Outer loop over the number of symmetries
C
      DO L=1, NSYM
         ISYM=L-1
         DO I=1, NOB(L)
            N=N+1
            NRI(N)=I
            NSM(N)=ISYM
         END DO
      END DO
      NOBT=N
C
C---- If requested, write out the orbital table to the printer.
C
      IF(NPFLG.GT.0)WRITE(IWRITE,500)(I,NRI(I),NSM(I),I=1,NOBT)
C
      RETURN
C
C---- Format Statements
C
 500  FORMAT(5X,'ORBITAL TABLE FOR THE SORTING PROCEDURE',//,5X,
     &       'N(OLD)',4X,'N(NEW)',9X,'M',//,(1X,I9,I10,2X,I10))
C
      END SUBROUTINE NOBID
!*==tmdvr1.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE TMDVR1(NFTDEN,NFTINT,NFTMT,IWRITE,npflg,ISW,nelt,
     &                  NUCCEN,NMSET,IDMAP,IDMAPL,LUMOM,NMSET2,maxtgt,
     &                  iplotfg,iqdfg,grounden,ksym,groundsym,
     &                  groundspin,qmoln,ukrmolp_ints)
C***********************************************************************
C
C     TMDVR1 - Transition Moments DriVer 1
C
C     Evaluates expectation values of one electron operators for
C     wavefunctions by multiplying density and transition density
C     matrices by transformed property integrals. This is the
C     supervisor program for the process. This module decides whether
C     Gaussian or Slater type functions have been used as the expansion
C     basis for the wavefunction and then passes control to an
C     appropriate lower level routine.
C
C     Input data:
C
C           NAME Character string defining the run
C         NFTDEN Logical unit for the density matrix library
C         NFTINT Logical unit for the transformed property integrals
C          NFTMT Logical unit for the target properties file
C             NR Integer workspace array
C           NCOR Size of the integer workspace array
C         IWRITE Logical unit for the printer
C            ISW Moment type switch
C         NUCCEN Sequence code designating which nuclear center is the
C                scattering center.
C          NMSET Set number at which the computed target properties are
C                to be written out.
C          IDMAP Defines unique sets of density matrix pairs which are
C                used to write an old style moments file - WRTMTH
C         IDMAPL Size of array IDMAP
C          LUMOM Logical unit for output of old style moments file as
C                in the original TMTJT series of codes.
C         NMSET2 set number for output on unit LUMOM for the first
C                set of data.
c         NOELEM Number of elements in density matrix
c                 (for lambda doubling)
c        ILOWINT Lower limit on integrals
c         IHINIT Upper limit on integrals
c        ILOWDEN Lower limit on density matrix
c         IHIDEN Upper limit on density matrix
c       SUBTRACT Matrix containing target density
c         KCOUNT Kth CSF from first wavefunction
c         JCOUNT Jth CSF from second wavefunction
c        MVALUES Contains info on mvalues to be considered
C          QMOLN Logical indicating whether we are running under Quantemol-N
C   ukrmolp_ints Are we using UKRmol+ property integrals
C
C     Linkage:
C
C           TMGP, TMSP, WRTMTH, WRTMT
C
C     Authors: Doug McLean, Bowen Liu and Megumu Yoshimine,
C              in the MOTECC-90 release of Alchemy II from
C              IBM San Jose, California, USA.
C
C     Change Log:
C
C       (i)  October 1993
C
C            Charles J Gillan, Queen's University of Belfast
C            modified to run with UK R-matrix scattering codes
C
C      (ii)  March 1994
C
C            Charles J Gillan, Queen's University of Belfast
C            added WRTMTH and WRTMT to generate old style output
C            which is used by polarizibility modules POLAR,GETPOL.
C
C     (iii)  July 2014
C            ZM interfaced with UKRmol+ integral code.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE global_utils, ONLY : MPROD, INITVR8
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL(KIND=wp), PARAMETER :: VSMALL=1.0E-10_wp
C---- NPMAX defines the maximum number of properties that we can
C     evaluate.
      INTEGER, PARAMETER :: NPMAX=17, MAXNUC=100, IPLOTFILE=44, 
     &                      IDENKEY=60, POTV=45
C
C Dummy arguments
C
      REAL(KIND=wp) :: GROUNDEN
      INTEGER :: GROUNDSPIN, GROUNDSYM, IDMAPL, IPLOTFG, IQDFG, ISW, 
     &           IWRITE, KSYM, LUMOM, MAXTGT, NELT, NFTDEN, NFTINT, 
     &           NFTMT, NMSET, NMSET2, NPFLG, NUCCEN
      INTEGER, DIMENSION(3,IDMAPL) :: IDMAP
      INTENT (IN) IDMAP, IDMAPL, IPLOTFG, NPFLG
      INTENT (OUT) GROUNDSYM
      INTENT (INOUT) GROUNDSPIN, NMSET2, NUCCEN
      LOGICAL :: QMOLN, ukrmolp_ints
      INTENT (IN) QMOLN, ukrmolp_ints
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(MAXNUC) :: CHARG
      CHARACTER(LEN=80) :: CHEAD
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: DEN, PROPS, XPROP
      REAL(KIND=wp), DIMENSION(2) :: DUM
      REAL(KIND=wp) :: ECI, ECJ
      REAL(KIND=wp), DIMENSION(maxtgt) :: ENERGY
      INTEGER, DIMENSION(maxtgt) :: GUTOT, IREFL, ISPIN, ISZ, MGVN
      INTEGER :: I, ICODE, IDENSET, IDENT, IGUI, IGUJ, IOLDMAP, IREFLI, 
     &           IREFLJ, ISPINI, ISPINJ, IST, ISTYLE, ISYMTYPI, 
     &           ISYMTYPJ, ISZI, ISZJ, IZERO, J, JST, LAMDAI, LAMDAJ, 
     &           LDENMAT, MAXDEN, MAXPTAB, MDEL, NCI, NCJ, NEWPROP, 
     &           NNUC, NOCSFI, NOCSFJ, NOPREC, NORB, NPROPS, NREC, NSRB, 
     &           NSTATI, NSTATJ, NSYM, NTGT
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISLATAB, KSLATCDE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISTATE, JSTATE
      INTEGER :: LEN
      CHARACTER(LEN=132) :: NAMED
      INTEGER, DIMENSION(30) :: NOB
      REAL(KIND=wp), DIMENSION(3,MAXNUC) :: RGEOM
      LOGICAL :: ZGAUSS
C
C*** End of declarations rewritten by SPAG
C
      DATA izero/0/
C
c        WRITE(IWRITE,1000)
c        WRITE(IWRITE,1010) NFTDEN,NFTINT,NFTMT,IWRITE,NCOR
c        WRITE(IWRITE,1020) IDMAPL
c        WRITE(IWRITE,1030)(IDMAP(1,I),IDMAP(2,I),IDMAP(3,I),I=1,IDMAPL)
C
C---- Write a title to the output
C
C      WRITE(6,8001) NAME
C
C---- If either of the output set numbers for the properties,
C     NMSET or NMSET2, is zero then we must locate the end of the
C     dataset and replace those numbers by the next set numbers
C     values. Also there is no point going on to compute moments
C     if the set numbers are invalid. These operations are carried
C     out by invoking the GETSET routine with appropriate
C     parameters.
C
C     Key numbers are:  9  for the formatted moments in WRTARG
C                      50  for the unformatted moments in WRTMTH
C
      DO i=1, len(named)
         named(i:i)=' '
      END DO
c
c     NFTMT holds finaL properties data
      CALL GETSET(NFTMT,NMSET,9,'FORMATTED  ')
      CALL GETSET(LUMOM,NMSET2,50,'UNFORMATTED')
C
      CALL GETGEON(NFTINT,IWRITE,NNUC,RGEOM,CHARG,izero,ukrmolp_ints)
C
C---- If NUCCEN is less than 1 then we scan the CHARG array to find the
C     scattering center; it has charge zero and is at the end of the
C     list! Should no charge center be found then NUCCEN is set to 1.
C
      IF(NUCCEN.LE.0)THEN
         DO I=NNUC, 1, -1
            IF(ABS(CHARG(I)).LT.VSMALL)NUCCEN=I
         END DO
      END IF
      IF(NUCCEN.LE.0)NUCCEN=1
C
C---- Read the first member on the density matrix library and thereby
C     determine the symmetry group of the density matrices. Of course it
C     is assumed that the library has been constructed from a single
C     run of the DENPROP code and that all members of the library are
C     constructed from the same orbital set, i.e. same point group too.
C
C     At this stage we also pick up the NOB array.
C
      IDENSET=0
      maxden=1
      CALL DENRD(NFTDEN,IDENSET,dum,LDENMAT,ICODE,CHEAD,ISPINI,ISZI,
     &           LAMDAI,IGUI,IREFLI,ISYMTYPI,ISPINJ,ISZJ,LAMDAJ,IGUJ,
     &           IREFLJ,ISYMTYPJ,NSYM,NOB,NOCSFI,NOCSFJ,NCI,NCJ,ECI,ECJ,
     &           NNUC,RGEOM,maxden)
C
C---- Compute absolute maximum value for the size of the density
C     matrix
C
      MAXDEN=0
      norb=0
      DO I=1, NSYM
         NORB=NORB+NOB(I)
         MAXDEN=MAXDEN+NOB(I)*NOB(I)
      END DO
c
      ALLOCATE(den(maxden))
C
      IF(ISYMTYPI.EQ.2)THEN
         NSRB=NORB*2
      ELSE
         NSRB=NOB(1)*2+(NORB-NOB(1))*4
      END IF
C
C---- From the information on the nuclear framework we may now compute
C     moment of inertia information as well as rotational constants.
C
      CALL MOIDRIV(NNUC,RGEOM,CHARG,IWRITE)
C
C---- If ISYMTYPI = 0 or 1 then we have C-inf-V type molecular point
C     group and this means that we have used Slater type functions in
C     basis set. On the other hand if ISYMTYPI=2 then the molecule has
C     and Abelian point group (D2h or subgroup) symmetry meaning that
C     the expansion basis is Gaussian functions.
C
C     We take this opportunity to trap possible errors too.
C
      IF(ISYMTYPI.EQ.0 .OR. ISYMTYPI.EQ.1)THEN
         ZGAUSS=.FALSE.
      ELSE IF(ISYMTYPI.EQ.2)THEN
         ZGAUSS=.TRUE.
      ELSE
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9960)ISYMTYPI
         STOP
      END IF
C
      IF(ISYMTYPI.NE.ISYMTYPJ)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9965)ISYMTYPI, ISYMTYPJ
         STOP
      END IF
C
      CALL INITST(ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,MAXTGT)
C
C---- Allocate storage for the properties table and initialize it too.
C
      MAXPTAB=NPMaX*MAXTGT*MAXTGT
      ALLOCATE(istate(maxptab),jstate(maxptab),islatab(8,maxptab),
     &         props(maxptab))
C
      CALL INIPTAB(ISTATE,JSTATE,ISLATAB,PROPS,MAXPTAB)
C
C---- Make dynamic storage allocation in the big vector to read the
C     density matrix
C
C     It is useful to remember that:
C     KSLATCDE is a buffer for the Slater property codes
C        XPROP is a buffer for final expectation values of the
C              wavefunction for each property operator in the Slater
C              case.
C
 
      ALLOCATE(kslatcde(8,npmax),xprop(npmax))
C
C---- Loop over all density matrices on the dataset and compute the
C     wavefunction properties associated with them. Note that we
C     dynamically compute the number of target states so we begin by
C     initializing the number to zero; the same is true of the number
C     entries in the properties table.
C
      NTGT=0
      NPROPS=0
      IOLDMAP=-100
C
 
      DO I=1, idmapl
         IDENSET=I
         CALL DENRD(NFTDEN,IDENSET,den,LDENMAT,ICODE,CHEAD,ISPINI,ISZI,
     &              LAMDAI,IGUI,IREFLI,ISYMTYPI,ISPINJ,ISZJ,LAMDAJ,IGUJ,
     &              IREFLJ,ISYMTYPJ,NSYM,NOB,NOCSFI,NOCSFJ,NCI,NCJ,ECI,
     &              ECJ,NNUC,RGEOM,maxden)
C
C....... Assign the state(s) to the unique state table and find
C        their positions. Whether we have one state, or two, is
C        controlled by the value of ICODE returned (cf. DENWRT)
C
         CALL ADDTST(NTGT,ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,IWRITE,IST,
     &               MAXTGT,ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ECI)
C
         IF(ICODE.EQ.2 .OR. ICODE.EQ.3)THEN
            CALL ADDTST(NTGT,ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,IWRITE,
     &                  JST,MAXTGT,ISPINJ,ISZJ,LAMDAJ,IGUJ,IREFLJ,ECJ)
         ELSE IF(ICODE.EQ.1)THEN
            JST=IST
         ELSE
            WRITE(IWRITE,9900)
            WRITE(IWRITE,9970)ICODE
            STOP
         END IF
C
C....... Zero-ize the auxiliary properties table for this pass
C
         CALL INITVR8(xPROP,Npmax)
C
C....... Compute properties now and add them to the properties table
C
 
         IF(ZGAUSS)THEN
            MDEL=MPROD(lamdaI+1,lamdaJ+1,0,IWRITE)-1
            CALL TMGP(NSYM,NOB,NFTINT,den,KSLATCDE,xprop,newprop,
     &                ldenmat,NPMaX,nuccen,mdel,ukrmolp_ints)
         ELSE
            MDEL=ABS(LAMDAJ-LAMDAI)
            CALL TMSP(IST,JST,nsym,nob,NFTINT,den,KSLATCDE,xprop,
     &                NEWPROP,npmax,MDEL,IWRITE)
         END IF
c
         CALL ADDPTAB(ISTATE,JSTATE,ISLATAB,PROPS,NPROPS,MAXPTAB,IST,
     &                JST,NEWPROP,KSLATCDE,xPROP,IWRITE)
 
C
C....... Place the properties upon the 'old style' output dataset !
C
C        If this wavefunction pair corresponds to the same symmetries
C        as the previous then the corresponding entry in IDMAP will be
C        the same. We need only write the moments ! However if it is a
C        new symmetry set then we must write a header first !
C
C        In order to write a header we must account for all the
C        records in the set and so on.
C
         IF(IDMAP(1,I).NE.IOLDMAP)THEN
            IOLDMAP=IDMAP(1,I)
            NREC=0
            DO J=1, IDMAPL
               IF(IDMAP(1,J).EQ.IOLDMAP)NREC=NREC+1
            END DO
            NREC=NREC*(NEWPROP+1)+3
            NAMED(1:LEN(CHEAD))=CHEAD
            IF(ICODE.EQ.1 .OR. ICODE.EQ.2)THEN
               IDENT=0
            ELSE IF(ICODE.EQ.3)THEN
               IDENT=1
               NREC=NREC+1
            END IF
            NOPREC=NEWPROP
            NSTATI=IDMAP(2,I)
            NSTATJ=IDMAP(3,I)
c
            CALL WRTMTH(LUMOM,NMSET2,nrec,NOPREC,NSTATI,NSTATJ,IDENT,
     &                  NUCCEN,NEWPROP,KSLATCDE,NORB,NSRB,NELT,NSYM,
     &                  ISYMTYPI,ISPINI,NOB,NNUC,RGEOM,CHARG,NAMED,
     &                  LAMDAI,IREFLI,IGUI,NOCSFI,NAMED,LAMDAJ,IREFLJ,
     &                  IGUJ,NOCSFJ)
            NMSET2=NMSET2+1
         END IF
C
         CALL WRTMT(LUMOM,ECI,ECJ,NCI,NCJ,LAMDAI,LAMDAJ,KSLATCDE,xPROP,
     &              NEWPROP)
C
C
      END DO
C
C
      IF(iplotfg.EQ.1)THEN
         istyle=-1
         IF(istyle.EQ.-1)THEN
            IF(isymtypi.LT.2)THEN
               istyle=0
            ELSE
               istyle=1
            END IF
         END IF
         IF(istyle.EQ.0)THEN
            CALL plotdpoten(iplotfile,NMSET,NNUC,NTGT,ISW,NUCCEN,RGEOM,
     &                      CHARG,ISPIN,MGVN,GUTOT,IREFL,ENERGY,ISTATE,
     &                      JSTATE,IWRITE,ISYMTYPI,NSYM,grounden,ksym)
         ELSE
            CALL plotppoten(iplotfile,NMSET,NNUC,NTGT,NUCCEN,RGEOM,
     &                      CHARG,ISPIN,MGVN,ENERGY,ISTATE,JSTATE,
     &                      IWRITE,ISYMTYPI,NSYM,grounden,ksym)
         END IF
      END IF
C---- Now that the unique state table is completed re-order it into
C     increasing energy order and write the target properties file
C
      CALL ORDTST(NTGT,ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,NPROPS,ISTATE,
     &            JSTATE)
C
!     Hemal Varambhia polarisability write statement
!      write(172,*) (mgvn(i),i=1,ntgt)
      groundsym=mgvn(1)
      groundspin=ispin(1)
      IF (QMOLN) WRITE(412,*)'ground state spin: ', groundspin
      IF(npflg.GT.0)THEN
         WRITE(IWRITE,4000)
         CALL PRTTST(NTGT,ISPIN,ISZ,MGVN,GUTOT,IREFL,ENERGY,ISYMTYPI,
     &               IWRITE)
         CALL PRTPTAB(ISTATE,JSTATE,ISLATAB,PROPS,NPROPS,IWRITE)
      END IF
C
 
      istyle=-1
      IF(istyle.EQ.-1)THEN
         IF(isymtypi.LT.2)THEN
            istyle=0
         ELSE
            istyle=1
         END IF
      END IF
      IF(istyle.EQ.0)THEN
         CALL WRDTARG(NFTMT,NMSET,NNUC,NTGT,NPROPS,ISW,NUCCEN,RGEOM,
     &                CHARG,ISPIN,MGVN,GUTOT,IREFL,ENERGY,ISTATE,JSTATE,
     &                ISLATAB,PROPS,IWRITE,ISYMTYPI,NSYM,iqdfg,grounden,
     &                ksym)
      ELSE
         CALL WRPTARG(NFTMT,NMSET,NNUC,NTGT,NPROPS,ISW,NUCCEN,RGEOM,
     &                CHARG,ISPIN,MGVN,ENERGY,ISTATE,JSTATE,ISLATAB,
     &                PROPS,IWRITE,ISYMTYPI,NSYM,iqdfg,grounden,ksym,
     &                ukrmolp_ints)
      END IF
c
      DEALLOCATE(den)
      DEALLOCATE(istate,jstate,islatab,props)
C
      RETURN
C
C---- Error condition handlers
C
C       (a) Error on open of the Density matrix library dataset
C
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9990)NFTDEN
C
      STOP
C
C---- Format statements
C
 6    FORMAT('xprop(',(i2),')= ',(i3))
 1000 FORMAT(10X,'====> TMDVR1 - CONTROL T.M. EVALUATION <====',//)
 1010 FORMAT(/,10X,'Logical unit for density matrices   = ',I3,/,10X,
     &       ' " " "   ""   "  property integrals = ',I3,/,10X,
     &       ' " " "   ""   "  target properties  = ',I3,/,10X,
     &       ' " " "   ""   "  printer            = ',I3,/,10X,
     &       'Size of integer workspace available = ',I10,/)
 1020 FORMAT(/,10X,'Length of mapping table for old style output = ',I7,
     &       //,10X,' IDMAP flag  NOVECI  NOVECJ ')
 1030 FORMAT(10X,I7,7X,I5,3X,I5)
 4000 FORMAT(///,' The following target state and properties tables are'
     &       ,' now in energy order')
 8001 FORMAT(///,5X,'CI PROPERTIES',/,5X,A,//,10X,
     &       'Moments printed in this output do NOT include the',/,10X,
     &       'nuclear terms. If the user has requested nuclear ',/,10X,
     &       'terms to beincluded (ISW=0) they are included in ',/,10X,
     &       'the output dump files written (in WRTARG) for the',/,10X,
     &       'external region codes. ',//)
C
 9900 FORMAT(/,5X,'**** Error in TMDVR1: ',/)
C
 9960 FORMAT(10X,'Value of ISYMTYPI read from density matrix library',/,
     &       10X,'is invalid. Must be 0,1 or 2 but is = ',I10,/)
 9965 FORMAT(10X,'ISYMTYPI and ISYMTYPJ read from density library ',/,
     &       10X,'are not identical but are = ',I10,1X,I10,/)
 9970 FORMAT(10X,'ICODE read from density library is not valid ',/,10X,
     &       'Should be 1,2 or 3 but is = ',I10,/)
 9990 FORMAT(10X,'Opening of the Density matrix library dataset',/,10X,
     &       'has failed on unit number = ',I3,/)
C
      END SUBROUTINE TMDVR1
!*==tmgp.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE TMGP(NSYM,NOB,NFT,DEN,ncod,PR,np,lden,NPMX,nuccen,mdel,
     &                ukrmolp_ints)
C***********************************************************************
C
C     TMGP - Supervises the computation of properties for wavefunctions
C            expanded in Gaussian type sets of functions
C
C     Input data:
C           NSYM Number of symmetries in the orbital set
C            NOB Number of orbitals per symmetry
C              D The density matrix for thwe wavefunction pair
C           NFT  Logical unit number for the property integrals
C           NCOR Size of R*8 core available for workspace
C             CR R*8 core workspace available
C           NPMX Maximum number of property operators that can be
C                handled
C   ukrmolp_ints Are we using UKRmol+ property integrals
C
C     Linkage:
C
C           TMGQ, GET_NAME_SYM, TMG_UKPLUS
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      USE global_utils, ONLY: mprod
      USE ukrmol_interface_gbl, ONLY: GET_NAME_SYM, TMG_UKPLUS
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL(KIND=wp), PARAMETER :: THRESH=1.E-10_wp
C
C Dummy arguments
C
      INTEGER :: LDEN, MDEL, NFT, NP, NPMX, NSYM, NUCCEN
      REAL(KIND=wp), DIMENSION(lden) :: DEN
      INTEGER, DIMENSION(8,npmx) :: NCOD
      INTEGER, DIMENSION(8) :: NOB
      REAL(KIND=wp), DIMENSION(npmx) :: PR
      LOGICAL :: ukrmolp_ints
      INTENT (IN) MDEL, NPMX, NUCCEN, ukrmolp_ints
      INTENT (OUT) NCOD
      INTENT (INOUT) NOB, NP, NSYM, PR
C
C Local variables
C
      INTEGER :: I, IPSTART, ISYM, JSYM, L, LAST, LBUF, LMQ, M, NALM, 
     &           NCODT, NINTS, NLMQ, NNSYM, NNUC, NOBT, Q, NOBSYM, J
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDV
      INTEGER, DIMENSION(9) :: ISTART
      CHARACTER(LEN=132) :: NAME
      INTEGER, DIMENSION(8) :: NNOB
      CHARACTER(LEN=8) :: PNAME
      REAL(KIND=wp) :: PRLMQ
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: XBUF
      LOGICAL :: sweden_ints
      INTEGER, ALLOCATABLE :: block_data(:,:)
      INTEGER :: block
C
C*** End of declarations rewritten by SPAG
C
      NALM=0
      NOBSYM=NSYM  ! actual length of incoming NOB
      np=0
C
      IF (.not.(ukrmolp_ints)) THEN
         sweden_ints = .true.
      ELSE
         sweden_ints = .false.
      ENDIF
C
      IF (ukrmolp_ints) THEN
C     
C------- Read header information from the UKRmol+ data file
C  
         CALL GET_NAME_SYM(NAME,NNSYM,NNOB,NLMQ)
C
      ENDIF
C
      IF (sweden_ints) THEN
C
C------- Read header information from the SWEDEN property integrals file
C
         REWIND NFT
         READ(NFT)NAME, nnsym, nobt, NNUC, NCODT, lbuf,
     &            (nnob(I),I=1,nnsym)
C
      ENDIF
C 
C---- Check that input is compatible with that from CONGEN.
C     This check must be the same regardless of where the integrals come from.
C
      IF(nsym.NE.nnsym)THEN
         IF(nsym.GT.nnsym)THEN
            IF(ANY(nob(nnsym+1:nsym).NE.0))GO TO 9
         ELSE
            IF(ANY(nnob(nsym+1:nnsym).NE.0))GO TO 9
            nob(nsym+1:)=0
         END IF
         nsym=nnsym
      END IF
      DO i=1, nsym
         IF(nnob(i).NE.nob(i))GO TO 9
      END DO
C
      ALLOCATE(block_data(3, NSYM*NSYM))
c
c.... Set up pointers to start of each block of transformed integrals.
c     Note that for each block (m) the values isym,jsym give the symmetries
c     of the orbitals involved in the corresponding block of the density matrix
c     and therefore determine exactly which property integrals are to be combined
c     with the density matrix.
c
      last=0
      block = 0
      DO m=1, nsym
         istart(m)=last
         DO isym=1, nsym
            jsym = mprod(isym, m, 0, 0)
            ! watch out - NOB may be shorter than congen's NSYM if it ends with zeros
            IF (MAX(isym,jsym) <= NOBSYM) THEN
               IF(isym.EQ.jsym)THEN
                  last=last+nob(isym)*(nob(isym)+1)/2
               ELSE
                  last=last+nob(isym)*nob(jsym)
               END IF
            END IF
            block = block + 1
            block_data(1:3,block) = (/m,isym,jsym/)
         END DO
      END DO
      istart(nsym+1)=last
c
      IF (sweden_ints) THEN
         READ(NFT)nlmq, nints
         IF(nlmq.GT.npmx)THEN
            WRITE(6,301)nlmq, npmx
            STOP 999
         END IF
C
C---- Storage allocation
C
         ALLOCATE(xbuf(nints))
         ALLOCATE(indv(nints))
C
      ENDIF !sweden_ints
C
C---- Read the Gaussian property integrals into core. If the process
C     goes wrong then the code terminates on receipt of the warning
C     signal in NALM. Note that NALM is not used in TMG_UKPLUS since all
C     errors are handled there immediately.
C
      DO lmq=1, nlmq
c
         IF (sweden_ints) THEN
            ipstart=istart(mdel+1)+1
            CALL TMGQ(NFT,lbuf,nints,indv,xbuf,DEN,lden,prlmq,l,m,q,
     &                pname,ipstart,nalm)
         ENDIF
C
         IF (ukrmolp_ints) THEN
            CALL TMG_UKPLUS(mdel+1,block_data,block,NOB,DEN,lden,
     &                      prlmq,lmq,l,m,q)
         ENDIF
C
         IF(nalm.NE.0)THEN
            WRITE(6,300)lmq, pname
            STOP 999
         END IF
c
c---- Set up Alchemy style property labels
         IF(lmq.EQ.1)THEN
            np=np+1
            ncod(1,np)=2
            ncod(2,np)=0
            ncod(3,np)=0
            ncod(4,np)=0
            ncod(5,np)=0
            ncod(6,np)=0
            ncod(7,np)=0
            ncod(8,np)=0
            pr(np)=prlmq+pr(np)
         ELSE IF(abs(prlmq).GT.thresh)THEN
            np=np+1
            ncod(1,np)=4
            ncod(2,np)=l
            ncod(3,np)=0
            ncod(4,np)=0
            ncod(5,np)=nuccen
            ncod(6,np)=l
            ncod(7,np)=m
            ncod(8,np)=q
            PR(np)=prlmq+PR(np)
            !ZM: the relation of the Molpro quadrupoles to the solid harmonic ones can be worked out from: http://www.molpro.net/pipermail/molpro-user/2005-July/001440.html
            !For example: (20) = QMZZ, (22) = (QMXX-QMYY)/sqrt(3)
         END IF
      END DO
c
      IF (sweden_ints) THEN
         DEALLOCATE(xbuf)
         DEALLOCATE(indv)
      ENDIF
C
      RETURN
c
 9    WRITE(6,900)
      WRITE(6,901)nsym, (nob(i),i=1,nsym)
      WRITE(6,902)nnsym, (nnob(i),i=1,nnsym)
      STOP
c
C---- Format Statements
C
 300  FORMAT(' ****** TMGP',5X,
     &       'ERROR READING PROPERTY INTEGRALS'/'         LMQ=',i3,
     &       '   Property ',a)
 301  FORMAT(' ****** TMGP',5X,'Too many properties',2I5)
 900  FORMAT(/' ****** Mismatch in symmetry data *****'/)
 901  FORMAT(' CONGEN   output ',9I7)
 902  FORMAT(' Property integrals output ',9I7)
C
      END SUBROUTINE TMGP
!*==tmgq.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE TMGQ(NFtprop,Lbuf,maxnin,indexv,xbuf,DEN,lden,pr,l,m,q,
     &                PNAME,ipstart,nalm)
C***********************************************************************
C
C     TMGQ -
C
C     Input data:
C           NFTPROP
C           LBUF
C             PR
C           NALM
C          PNAME
C
C     Linkage:
C
C         none
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
c.... rootp=1/sqrt(2*pi)
      REAL(KIND=wp), PARAMETER :: ROOTP=0.398942280401433_wp
C
C Dummy arguments
C
      INTEGER :: IPSTART, L, LBUF, LDEN, M, MAXNIN, NALM, NFTPROP, Q
      CHARACTER(LEN=8) :: PNAME
      REAL(KIND=wp) :: PR
      REAL(KIND=wp), DIMENSION(lden) :: DEN
      INTEGER, DIMENSION(maxnin) :: INDEXV
      REAL(KIND=wp), DIMENSION(maxnin) :: XBUF
      INTENT (IN) DEN, IPSTART, LBUF, LDEN, MAXNIN, NFTPROP
      INTENT (OUT) L, M, NALM, PNAME, PR, Q
      INTENT (INOUT) INDEXV, XBUF
C
C Local variables
C
      INTEGER :: I, IFINISH, ISTART, J, LMQ, N, NBUF, NINTS, NLEFT, NREC
      REAL(KIND=wp) :: SUM
C
C*** End of declarations rewritten by SPAG
C
      ISTART=1
      nalm=0
      pr=xzero

c
      READ(nftprop,ERR=900)PNAME, lmq, nrec, nints, l, m, q
c      WRITE(6,5001) pname,nints
c
      IF(nints.GT.0)THEN
C
C---- Loop over records
C
         nleft=nints
 
         DO i=1, nrec
            IF(nleft.LT.lbuf)THEN
               ifinish=istart+nleft-1
            ELSE
               IFINISH=ISTART+LBUF-1
               nleft=nleft-lbuf
            END IF
            READ(nftprop,ERR=900)(Xbuf(J),J=ISTART,IFINISH), 
     &                           (indexv(j),J=ISTART,IFINISH), nBUF
            ISTART=IFINISH+1
c
         END DO
         IF(nbuf.LE.0)GO TO 900
c
 
         IF(ipstart.GT.0)THEN
            SUM=XZERO
            DO N=1, lden
               DO i=1, nints
                  IF(n.EQ.indexv(i)-ipstart+1)THEN
                     SUM=SUM+DEN(N)*XBUF(i)
                  END IF
               END DO
            END DO
 
c
c          if(l.gt.0) sum = SUM*rootp
            pr=sum
c          write(6,120) pr
         END IF
C
      END IF
C
      RETURN
c
 900  nalm=1
      RETURN
C
C---- Format statements
C
 120  FORMAT(' Final property',d20.12)
! Hemal Varambhia de-bug
 121  FORMAT('DEN(',i3,')=',d20.12,' xbuf(',i3,')=',d20.12)
!
 
 5001 FORMAT(/' Number of property integrals of type ',a8,' read is',
     &       i10)
C
      END SUBROUTINE TMGQ
!*==tmrs19.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE TMRS19(NFTINT,XINT,NINT,NCOD,mxint,IWRITE,IERR,IEND)
C***********************************************************************
C
C     TMRS19 - Transition Moments Read Slater property integrals
C              from unit 19.
C
C     Reads a record of transformed Slater property integrals from
C     unit NFTINT
C
C     Input data:
C         NFTINT Logical unit for the property integrals
C         IWRITE Logical unit for the printer
C
C     Output data:
C            NINT  Number of integrals read
C            IERR  Signal for error on read of integrals
C                  =0  no error, or problem, occurred on read
C                   1  there has been an error during read
C            IEND  Signal for end of file on read of integrals
C                  =0  an attempt was made to read a record of
C                      integrals (see IERR)
C                   1  the end of file was encountered on trying
C                      to read the integrals record
C
C     Input/Output data:
C                  NCOD Storage location where integral codes are
C                       placed upon reading
C                  XINT Storage location where integrals are placed
C                       upon reading
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IEND, IERR, IWRITE, MXINT, NFTINT, NINT
      INTEGER, DIMENSION(8) :: NCOD
      REAL(KIND=wp), DIMENSION(mxint) :: XINT
      INTENT (IN) IWRITE, MXINT, NFTINT
      INTENT (OUT) IEND, IERR, NCOD, XINT
      INTENT (INOUT) NINT
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
C---- Initialize the return codes
C
      IERR=0
      IEND=0
C
C---- Read a record of integrals
C
      READ(NFTINT,ERR=80,END=85)(NCOD(I),I=1,8), NINT, 
     &                          (XINT(I),I=1,NINT)
C
c        WRITE(IWRITE,2010) (NCOD(I),I=1,8),NINT
c        WRITE(IWRITE,2020) (XINT(I),I=1,NINT)
C
      RETURN
C
C---- Error condition handlers
C
C...... (a) Error on read of data
C
 80   CONTINUE
C
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9920)
C
      IERR=1
      RETURN
C
C...... (b) End of file on read - Not strictly an error within the
C                                 context of routine TMGP.
C
 85   CONTINUE
C
      IEND=1
      RETURN
C
C---- Format Statements
C
 2010 FORMAT(/,20X,'Property Integral Code = ',8I3,/,20X,
     &       'No. of Integrals       = ',I10,/)
 2020 FORMAT(20X,'Integrals follow: ',/,(20X,3(F15.7,1X),/))
 9900 FORMAT(/,10X,'**** Error in TMRS19: ',/)
 9920 FORMAT(10X,'Failed to read Slater property integrals ',
     &       'from unit = ',I3,/)
C
      END SUBROUTINE TMRS19
!*==tmsp.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE TMSP(IST,JST,nsym,nob,NFTINT,DEN,NCOD,pr,NEWPROP,npmax,
     &                MDEL1,IWRITE)
C***********************************************************************
C
C     TMSP - Computes expectation values of property operators for
C            wavefunctions expanded in a basis of Slater type
C            functions.
C
C     Input data:
C
C            IST Designation of the I state
C            JST Designation ""  "  J  " "
C           NSYM No. of C-inf-V symmetries in the basis set
C            NOB No. of orbitals per symmetry
C         NFTINT Logical unit for the transformed property integrals
C            DEN The density matrix for the wavefunction pair
C           XINT Buffer for reading the property integrals
C           NCOD Buffer for reading the property integral codes
C             PR Buffer for storing computed properties
C          MDEL1 The delta lambda for the wavefunction pair making up
C                the density matrix. (NO LAMBDA DOUBLING)
C          MDEL2 The delta lambda for the wavefunction pair making up
C                the density matrix. (LAMBDA DOUBLING)
C        ISYMTYP The symmetry type of molecule
C         IWRITE Logical unit for the printer
c******************************************************************
C*****See the input to TMDVR1 to find definitions of remaining input
C*****Variables
c****************************************************************
C     Output data:
C         NEWPROP The number of computed properties
C
C                 NCOD and PR are used as storage areas to pass data
C                 forward.
C
C     Linkage:
C
C          TMRS19
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO, TWO=>XTWO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IST, IWRITE, JST, MDEL1, NEWPROP, NFTINT, NPMAX, NSYM
      REAL(KIND=wp), DIMENSION(*) :: DEN
      INTEGER, DIMENSION(8,npmax) :: NCOD
      INTEGER, DIMENSION(*) :: NOB
      REAL(KIND=wp), DIMENSION(npmax) :: PR
      INTENT (IN) DEN, IST, JST, MDEL1, NOB, NPMAX, NSYM
      INTENT (OUT) NEWPROP
      INTENT (INOUT) PR
C
C Local variables
C
      INTEGER :: I, IEND, IERR, ISYM, J, M, MAXDEN, N, NINT, NP
      REAL(KIND=wp) :: SUM, TOL
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: XINT
C
C*** End of declarations rewritten by SPAG
C
      DATA tol/1.E-10_wp/
C
c        WRITE(IWRITE,1000)
c        WRITE(IWRITE,1010) NFTINT,NSYM,(NOB(I),I=1,NSYM)
c        WRITE(IWRITE,1020) IST,JST
C
C---- Prepare the transformed property integrals for reading
C
      REWIND NFTINT
C
C---- Define the number of properties computed to be zero here.
C
      NEWPROP=0
C
C---- Compute absolute maximum number of integrals
C
      MAXDEN=0
      DO I=1, NSYM
         MAXDEN=MAXDEN+NOB(I)*NOB(I)
      END DO
c
      ALLOCATE(xint(maxden))
C
      READ(NFTINT)
      READ(NFTINT)
C
C---- Loop over all of the property operators on the dataset and
C     compute the expectation value for the wavefunction. The
C     following is a repeat until end of file loop beginning at line
C     100 terminating with a branch to line 200. Not all property
C     integral codes, on NFTINT, will be used to compute expectation
C     values thus we must count the ones that are using variable NP.
C
      NP=0
C
 100  NP=NP+1
      IF(np.GT.1)THEN
         IF(abs(pr(np-1)).LT.tol .AND. ncod(1,np-1).EQ.4)np=np-1
      END IF
C
C...... Read a record of property integrals in preparation for
C       computation of the expectation value. The end of file is
C       signalled by having IEND=0 on return.
C
      CALL TMRS19(NFTINT,XINT,NINT,NCOD(1,NP),maxden,IWRITE,IERR,IEND)
C
      IF(IERR.NE.0 .OR. IEND.NE.0)GO TO 200
c
      IF(np.GT.npmax)THEN
         WRITE(iwrite,301)np, npmax
         STOP 999
      END IF
C
      M=NCOD(7,NP)
c
C...... Should the operator have a code not equal to that of the
C       wavefunction pair difference/sum then we cannot compute the
C       expectation value for it.
c
      IF(M.NE.MDEL1)GO TO 100
C
C..... Compute the expectation value by multiplying the density matrix
C      by the vector of property integrals.
C
      SUM=XZERO
      n=0
      IF(mdel1.EQ.0)THEN
         DO isym=1, nsym
            DO i=1, nob(isym)
               DO j=1, i-1
                  n=n+1
                  SUM=SUM+two*DEN(N)*XINT(N)
               END DO
               n=n+1
               SUM=SUM+DEN(N)*XINT(N)
            END DO
         END DO
      ELSE
         DO n=1, nint
            SUM=SUM+DEN(N)*XINT(N)
         END DO
      END IF
C
      IF(IST.NE.JST)PR(NP)=XZERO
      PR(NP)=SUM+PR(NP)
C
      GO TO 100
C
C---- Branch here at the end of loop
C
 200  CONTINUE
C
      NEWPROP=NP-1
C
      REWIND NFTINT
      DEALLOCATE(xint)
C
      RETURN
C
C---- Format Statements
C
 301  FORMAT(' ****** TMSP',5X,'Too many properties',2I5)
 1000 FORMAT(//,20X,'====> TMSP - COMPUTE SLATER PROPERTIES <====',//)
 1010 FORMAT(/,20X,'Logical unit for the property integrals = ',I3,/,
     &       20X,'No. of symmetries in the orbital set    = ',I3,/,20X,
     &       'Orbital set = ',3(10I3,/,30X))
 1020 FORMAT(/,20X,'I state designation = ',I5,/,20X,
     &       'J state designation = ',I5,/)
C
      END SUBROUTINE TMSP
!*==tmtclos.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE TMTCLOS()
C***********************************************************************
C
C     TMTCLOS - Shuts down TMT program when an error is found
C
C     All TMT errors should invoke this routine as the shutdown process,
C     and not the system statement STOP alone.
C
C***********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C*** End of declarations rewritten by SPAG
C
C
      WRITE(6,1000)
C
C---- Write the message
C
      WRITE(6,1010)
      WRITE(6,1020)
      WRITE(6,1030)
      WRITE(6,1020)
      WRITE(6,1010)
C
      STOP
C
C---- Format Statements
C
 1000 FORMAT(///)
 1010 FORMAT(10X,'***********************************************')
 1020 FORMAT(10X,'*                                             *')
 1030 FORMAT(10X,'*  DENPROP IS SHUTTING DOWN DUE TO AN ERROR   *')
C
      END SUBROUTINE TMTCLOS
!*==tmtma.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE TMTMA(NODA,CA,NDA,NODB,CB,NDB,nsrb,NELT,NDTR,MN,MS,
     &                 MDTR,NFTD,NELMT,NBLK,NEL,LENGTH,II,JJ,NPQ,CPQ,
     &                 ICSF,JCSF)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ICSF, JCSF, LENGTH, NBLK, NEL, NELMT, NELT, NFTD, NODA, 
     &           NODB, NSRB
      REAL(KIND=wp), DIMENSION(*) :: CA, CB, CPQ
      INTEGER, DIMENSION(*) :: II, JJ, MDTR, MN, MS, NDA, NDB, NDTR
      INTEGER, DIMENSION(2,*) :: NPQ
      INTENT (IN) CA, CB, MDTR, MN, MS, NDA, NDTR, NELT, NODA, NODB, 
     &            NSRB
C
C Local variables
C
      REAL(KIND=wp) :: CDA, CFD
      INTEGER :: I, IC, ID, IORB1, IORB2, J, JA, JB, M, MA, MAA, MB, 
     &           MBB, MDA, MDB, N
      INTEGER :: ISRCHE
      INTEGER, DIMENSION(2) :: NDD
      INTEGER, DIMENSION(nsrb) :: NDTA
      INTEGER, DIMENSION(nsrb+1) :: NDTB
      INTEGER, DIMENSION(4) :: NDTC
      INTEGER, DIMENSION(nelt) :: NDTD
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     TMTMA -  COMPUTES OVERLAPS BETWEEN CSF PAIRS.
C
C              FOR OFF-DIAGONAL ELEMENTS OF THE DENSITY MATRIX, PAIRS OF
C              CSFS MUST BE CONSIDERED. WITHIN THIS OF COURSE IS IMPLIED
C              A SUMMATION OVER THE PAIRS OF DETERMINANTS WHICH
C              CONSTITUTE THE CSFS. THIS PROCEDURE IS ALSO REQUIRED IN
C              THE CONSTRUCTION OF THE HAMILTONIAN MATRIX AND THERE IS A
C              SUBROUTINE ENRGMA IN PROGRAM SPEEDY WHICH PERFORMS THIS
C              TASK. TMTMA IS A DIRECT ANALOGUE OF THAT ROUTINE BUT IS
C              NECESSARILY SIMPLER BECAUSE OF THE FACT THAT DETERMINANT
C              PAIRS MUST DIFFER BY NO MORE THAN TWO SPIN ORBITALS IN
C              ORDER TO CONTRIBUTE TO THE FIRST ORDER DENSITY MATRIX.
C
C     INPUT DATA :
C           NODA  NUMBER OF DETERMINANTS IN THE FIRST CSF
C             CA  COEFFICIENTS FOR EACH DETERMINANT IN THE FIRST CSF
C            NDA  DETERMINANTS IN THE FIRST CSF
C           NODB  NUMBER OF DETERMINANTS IN THE SECOND CSF
C             CB  COEFFICIENTS FOR EACH DETERMINANT IN THE SECOND CSF
C            NDB  DETERMINANTS IN THE SECOND CSF
C           NELT  NUMBER OF ELECTRON IN THE SYSTEM
C             MN  TABLE GIVING THE ORBITAL NUMBER FOR EACH SPIN ORBITAL
C             MS  TABLE GIVING SZ FOR EACH SPIN ORBITAL
C           MDTR  POINTER ARRAY LINKING THE REFERENCE DETERMINANT TO
C                 ALL SPIN ORBITALS
C
C     DATA FOR THE OUTPUT ROUTINE (SEE COMMENTS THEREIN):
C
C       NFTD,NELMT,NBLK,NEL,LENGTH,II,JJ,NPQ,CPQ,ICSF,JCSF,IWRITE
C
C     LOCAL DATA :
C           NDTC  WORKSPACE
C            NDD  WORKSPACE
C
C***********************************************************************
C
C-----COPY REFERENCE DETERMINANT INTO NDTD AND ZEROIZE THE COEFFICIENT
C     FOR THIS PAIR OF CSFS
C
      DO I=1, NELT
         NDTD(I)=NDTR(I)
      END DO
C
C-----OUTER LOOP OVER DETERMINANTS
C
      MDA=1
      DO IC=1, NODA
         CDA=CA(IC)
         MA=NDA(MDA)
         MAA=MDA+MA
         MDB=1
C
C--------INNER LOOP OVER DETERMINANTS
C
         DO ID=1, NODB
            CFD=CDA*CB(ID)
            MB=NDB(MDB)
            MBB=MDB+MB
C
C--------BASED UPON THE NUMBER OF REPLACEMENTS FROM THE REFERENCE
C        DETERMINANT, JUMP TO APPROPRIATE PIECE OF CODE.
C          A. BOTH DETERMINANTS ARE NOT THE REF. DET. GOTO 200
C          B. ONE IS THE REF DET. OTHER HAS > ONE REPLS. GOTO 300
C          C. BOTH ARE THE REF. DET. GOTO 270
C          D. ONE IS THE REF. DET. OTHER HAS ONE REPL. GOTO 260
C
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
C
C---------PROCESS THE DIFFERENT SITUATIONS
C
C---------NEITHER IS THE REFERENCE DETERMINANT. THE FOLLOWING LINES OF C
C         IMPLEMENT AN ALGORITHM DESCRIBED BY NESBET IN HIS BOOK ON
C         VARIATIONAL METHODS IN ELECTRON-ATOM SCATTERING THEORY. THE
C         ALGORITHM EXPRESSES ON DETERMINANT AS AN EXCITATION OF THE OTH
C         RATHER THAN BOTH BEING EXCITATIONS RELATIVE TO THE REFERENCE
C         DETERMINANT. FOR A FULLER DESCRIPTION OF THE ALGORITHM SEE PAG
C         OF THE BOOK.
C
C
C--------- STEP ONE: COPY ALL REPLACED SPIN-ORBS INTO NDTB, AND
C                    ALL REPLACEMENT SPIN-ORBS INTO NDTA
C
 200        DO I=1, MA
               NDTA(I)=NDA(MAA+I)
               NDTB(I)=NDA(MDA+I)
            END DO
C
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
C
            IF((JA-MA).GT.1)GO TO 302
            JB=0
C
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
C
            IF(JB.EQ.0)GO TO 270
            IF(JB.EQ.1)THEN
               NDTC(3)=NDTA(NDD(1))
               NDTC(4)=NDTB(NDD(1))
               GO TO 260
            END IF
C
C........ ONE PAIR IS DIFFERENT
C
C.......... CHECK FOR IDENTICAL SPIN VALUES
            IF(MS(NDTC(3)).NE.MS(NDTC(4)))GO TO 302
C
 260        IF(MN(NDTC(3)).LT.MN(NDTC(4)))THEN
               IORB1=MN(NDTC(3))
               IORB2=MN(NDTC(4))
            ELSE
               IORB1=MN(NDTC(4))
               IORB2=MN(NDTC(3))
            END IF
            CALL WRTOFD(NFTD,NELMT,NBLK,NEL,LENGTH,ICSF,JCSF,CFD,IORB1,
     &                  IORB2,II,JJ,NPQ,CPQ)
            GO TO 302
C
C......... DETERMINANTS ARE IDENTICAL : THIS CASE DOES HAPPEN AFTER ALL!
C
 270        CONTINUE
C
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=NDA(MAA+I)
            END DO
C
            DO M=1, NELT
               N=NDTD(M)
               IORB1=MN(N)
               IORB2=IORB1
 
               CALL WRTOFD(NFTD,NELMT,NBLK,NEL,LENGTH,ICSF,JCSF,CFD,
     &                     IORB1,IORB2,II,JJ,NPQ,CPQ)
            END DO
C
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=N
            END DO
C
C......... ASCEND IN THE COUPLING LOOPS
C
 302        MDB=MBB+MB+1
         END DO
         MDA=MAA+MA+1
      END DO
C
      RETURN
C
      END SUBROUTINE TMTMA
!*==jacobi.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
 
      SUBROUTINE JACOBI(A,B,C,N)
C***********************************************************************
C
C     JACOBI - Use the JACOBI method to diagonalize a real symmetric
C              matrix
C
C     Input data:
C              A  lower triangle of a real symmetric matrix
C              N  order of the matrix
C
C     Output data:
C               C array of eigenvalues of the matrix
C               B matrix of the eigenvectors
C
C     Notes:
C
C       This routine was obtained from the NYU properties code in the
C     form used in the Alchemy II package.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: N
      REAL(KIND=wp), DIMENSION(*) :: A, B, C
      INTENT (IN) N
      INTENT (OUT) C
      INTENT (INOUT) A, B
C
C Local variables
C
      REAL(KIND=wp) :: ALPHA, AMAX, BETA, CC, CSQ, EPS, ONE, ROOT, S, 
     &                SSQ, SUM, TEMPA, TERM, THREE, THRESH, THRSHG, 
     &                TWOSC, ZERO
      INTEGER :: I, ID, II, J, JD, JJ, JJJ, K, KA, KB, KC, KD, KI, KJ, 
     &           L, LL, LOOPC, LOOPMX, M, NA, NN, NNA
C
C*** End of declarations rewritten by SPAG
C
C---- Data initialization
C
      DATA EPS/1.0E-8_wp/, LOOPMX/5000/, 
     &     zero, one, three/0._wp, 1._wp, 3._wp/
C
      TERM=zero
C
      LOOPC=0
      NA=N
      NN=(NA*(NA+1))/2
      K=1
      DO I=1, NA
         DO J=1, NA
            B(K)=zero
            IF(I.EQ.J)B(K)=one
            K=K+1
         END DO
      END DO
      SUM=zero
      NNA=NA-1
      IF(NA.EQ.1)GO TO 170
      K=1
      AMAX=zero
      DO I=1, NA
         DO J=1, I
            IF(I.EQ.J)GO TO 30
            IF(ABS(A(K)).GT.AMAX)AMAX=ABS(A(K))
 30         TERM=A(K)*A(K)
            SUM=SUM+TERM+TERM
            K=K+1
         END DO
         SUM=SUM-TERM
      END DO
      SUM=SQRT(SUM)
      THRESH=SUM/SQRT(real(NA,kind=wp))
      THRSHG=THRESH*EPS
      IF(THRSHG.GE.AMAX)GO TO 170
      THRESH=AMAX/three
      IF(THRESH.LT.THRSHG)THRESH=THRSHG
 60   K=2
      M=0
      JD=1
      KJ=0
      DO J=2, NA
         KJ=KJ+NA
         JD=JD+J
         JJ=J-1
         JJJ=JJ-1
         ID=0
         KI=-NA
         DO I=1, JJ
            KI=KI+NA
            ID=ID+I
            IF(ABS(A(K)).LE.THRESH)CYCLE
            M=M+8
            ALPHA=(A(JD)-A(ID))/(A(K)+A(K))
            BETA=0.25_wp/(1.0_wp+ALPHA**2)
            ROOT=0.5_wp+ABS(ALPHA)*SQRT(BETA)
            IF(ALPHA.GE.0._wp)GO TO 70
            CSQ=BETA/ROOT
            SSQ=ROOT
            GO TO 80
 70         SSQ=BETA/ROOT
            CSQ=ROOT
 80         CC=SQRT(CSQ)
            S=-SQRT(SSQ)
            TWOSC=CC*(S+S)
            TEMPA=CSQ*A(ID)+TWOSC*A(K)+SSQ*A(JD)
            A(JD)=CSQ*A(JD)-TWOSC*A(K)+SSQ*A(ID)
            A(ID)=TEMPA
            A(K)=0._wp
            KA=JD-J
            KB=ID-I
            KC=KI
            KD=KJ
            II=I-1
            IF(II.EQ.0)GO TO 100
            DO L=1, II
               KC=KC+1
               KD=KD+1
               TEMPA=CC*B(KC)+S*B(KD)
               B(KD)=-S*B(KC)+CC*B(KD)
               B(KC)=TEMPA
               KB=KB+1
               KA=KA+1
               TEMPA=CC*A(KB)+S*A(KA)
               A(KA)=-S*A(KB)+CC*A(KA)
               A(KB)=TEMPA
            END DO
 100        KC=KC+1
            KD=KD+1
            TEMPA=CC*B(KC)+S*B(KD)
            B(KD)=-S*B(KC)+CC*B(KD)
            B(KC)=TEMPA
            KB=KB+1
            KA=KA+1
            IF(I.EQ.JJ)GO TO 120
            DO L=I, JJJ
               KC=KC+1
               KD=KD+1
               TEMPA=CC*B(KC)+S*B(KD)
               B(KD)=-S*B(KC)+CC*B(KD)
               B(KC)=TEMPA
               KB=KB+L
               KA=KA+1
               TEMPA=CC*A(KB)+S*A(KA)
               A(KA)=-S*A(KB)+CC*A(KA)
               A(KB)=TEMPA
            END DO
 120        KC=KC+1
            KD=KD+1
            TEMPA=CC*B(KC)+S*B(KD)
            B(KD)=-S*B(KC)+CC*B(KD)
            B(KC)=TEMPA
            KB=KB+JJ
            KA=KA+1
            IF(J.EQ.NA)CYCLE
            DO L=J, NNA
               KC=KC+1
               KD=KD+1
               TEMPA=CC*B(KC)+S*B(KD)
               B(KD)=-S*B(KC)+CC*B(KD)
               B(KC)=TEMPA
               KB=KB+L
               KA=KA+L
               TEMPA=CC*A(KB)+S*A(KA)
               A(KA)=-S*A(KB)+CC*A(KA)
               A(KB)=TEMPA
            END DO
            K=K+1
         END DO
         K=K+1
      END DO
      LOOPC=LOOPC+1
      IF(LOOPC.GT.LOOPMX)GO TO 200
      IF(M.GT.NN)GO TO 60
      IF(THRESH.EQ.THRSHG)GO TO 160
      THRESH=THRESH/3._wp
      IF(THRESH.GE.THRSHG)GO TO 60
      THRESH=THRSHG
      GO TO 60
 160  IF(M.NE.0)GO TO 60
 170  LL=0
      DO L=1, NA
         LL=LL+L
         C(L)=A(LL)
      END DO
C
C---- Subroutine return point
C
      RETURN
C
C---- Error condition handler - No convergence
C
 200  PRINT 210, LOOPMX
C
      STOP
C
C---- Format Statements
C
 210  FORMAT(/' NO CONVERGENCE IN',I4,' ITERATIONS')
C
      END SUBROUTINE JACOBI
!*==wfndmp.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WFNDMP(ICOEFP,IDETP,COEFS,IDETS,LC,LN,NOCSF,IWRITE)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, LC, LN, NOCSF
      REAL(KIND=wp), DIMENSION(LC) :: COEFS
      INTEGER, DIMENSION(NOCSF) :: ICOEFP, IDETP
      INTEGER, DIMENSION(LN) :: IDETS
      INTENT (IN) COEFS, ICOEFP, IDETP, IDETS, IWRITE, LC, LN, NOCSF
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     WFNDMP - WAVEFUNCTION DUMP
C
C     WRITES OUT THE DETERMINANT INFORMATION FROM THE WAVEFUNCTION
C     FILES TO THE OUTPUT STREAM.
C
C***********************************************************************
C
 10   FORMAT(/,5X,'Determinant Positions',//,(5X,20(I5,1X)))
 20   FORMAT(/,5X,'Determinants',//,(5X,20(I5,1X)))
 30   FORMAT(/,5X,'Coefficient Positions',//,(5X,20(I5,1X)))
 40   FORMAT(/,5X,'Coefficients',//,(5X,8(E14.7,1X)))
C
      WRITE(IWRITE,10)IDETP
      WRITE(IWRITE,20)IDETS
      WRITE(IWRITE,30)ICOEFP
      WRITE(IWRITE,40)COEFS
C
      RETURN
      END SUBROUTINE WFNDMP
!*==wrdtarg.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRDTARG(NFTMT,NMSET,NNUC,NTARG,NPROPS,ISW,NUCCEN,RGEOM,
     &                   CHARG,ISPIN,MGVN,GUTOT,IREFL,ENERGY,ISTATE,
     &                   JSTATE,ICODES,PROPS,IWRITE,ISYMTYP,NSYM,iqdfg,
     &                   grounden,ksym)
C***********************************************************************
C
C     WRDTARG - WRite a TARGet properties file for linear molecule
C
C     Linkage:
C
C        GETSET,CMAKSTN
C
C     Notes:
C
C     This is the partner of the routine RDTARG in the external region
C     codes. The properties file is always formatted and all records
C     have a fixed length of 80 bytes. They are written (and read) with
C     the format
C
C                   I1, 7I3, D20.12, 2X, A36
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: KEYH=9, KEYG=8, KEYT=5, KEYP=1, NULL=0
      CHARACTER(LEN=11), PARAMETER :: TFORM='FORMATTED  '
C
C Dummy arguments
C
      REAL(KIND=wp) :: GROUNDEN
      INTEGER :: IQDFG, ISW, ISYMTYP, IWRITE, KSYM, NFTMT, NMSET, NNUC, 
     &           NPROPS, NSYM, NTARG, NUCCEN
      REAL(KIND=wp), DIMENSION(*) :: CHARG, ENERGY, PROPS
      INTEGER, DIMENSION(*) :: GUTOT, IREFL, ISPIN, ISTATE, JSTATE, MGVN
      INTEGER, DIMENSION(8,*) :: ICODES
      REAL(KIND=wp), DIMENSION(3,*) :: RGEOM
      INTENT (IN) CHARG, ENERGY, GROUNDEN, ICODES, IQDFG, ISTATE, ISW, 
     &            IWRITE, JSTATE, NPROPS, NTARG, NUCCEN, RGEOM
      INTENT (INOUT) PROPS
C
C Local variables
C
      CHARACTER(LEN=36), EXTERNAL :: CMAKSTN
      INTEGER :: DEND, I, ICHARG, ICOUNT, ICOUNT1, ICOUNT2, INUC, ITYPE, 
     &           J, K, MDIFF1, MDIFF2, MNUC, NOADC, NRECS, 
     &           ORDPOT=55
      REAL(KIND=wp) :: DIP1, DIP2, DIP3, GS, R, RAB, TOL=1.E-10_wp, Z
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDEX1, INDEX2
      INTEGER :: INT, LEN
      INTEGER, DIMENSION(3) :: MASS
      REAL(KIND=wp), DIMENSION(3) :: RMOI
      CHARACTER(LEN=36) :: RNAME
      REAL(KIND=wp), DIMENSION(nnuc,8) :: VLIST
      CHARACTER(LEN=8), DIMENSION(nnuc) :: VSYMB
      LOGICAL :: ZISWIS0
C
C*** End of declarations rewritten by SPAG
C
      DATA dend/288/, gs/0._wp/
C
      WRITE(IWRITE,1000)
      WRITE(IWRITE,1010)NFTMT, NMSET, TFORM
      WRITE(IWRITE,1020)NNUC, NUCCEN, NTARG, NPROPS, ISW
c        WRITE(IWRITE,1030)
c        DO 101 I=1,NNUC
c        WRITE(IWRITE,1040) CHARG(I),(RGEOM(J,I),J=1,3)
c  101   CONTINUE
C
C---- Error check the input as much as possible
C
C..... Must have a +ve number of nuclei and target states too.
C
      IF(NNUC.LT.1 .OR. NTARG.LT.1)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9910)NNUC, NTARG
         STOP
      END IF
C
C..... ISW must be zero or unity.
C
      IF(ISW.NE.0 .AND. ISW.NE.1)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9920)ISW
         STOP
      END IF
C
C..... NUCCEN must be in the range 1 <= NUCCEN <= NNUC
C
      IF(NUCCEN.LT.1 .OR. NUCCEN.GT.NNUC)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9930)NUCCEN, NNUC
         STOP
      END IF
C
C---- Initialize the output character string to blanks
C
      DO I=1, LEN(RNAME)
         RNAME(I:I)=' '
      END DO
C
C---- First position NFTMT at the begining of set NTSET
C
      CALL GETSET(NFTMT,NMSET,KEYH,TFORM)
C
C---- Must scan the properties records an extract only those that are
C     going to be needed in an external region run. They are
C
C     r**l Y_{lm}   centered on the NUCCEN nucleus
C
C     Also we must make sure that delta lamda matches for the target
C     state pair and the operator.
C
C     Operators with M > 1 are ignored at present in this routine !!!!!
C
C     Once the density matrix section has been amended this will change
C
C     All others are ignored in this routine !
C
C     An index vector is built which is later used to extract the
C     elements from the appropriate arrays; a second index vector
C     corresponding to the sub set of these selected moments having
C     target states with themselves is also built. The second vector
C     is used to subtract nuclear moment contributions.
C
      ICOUNT=0
      ICOUNT2=0
      ZISWIS0=ISW.EQ.0
      ALLOCATE(index1(nprops))
      ALLOCATE(index2(nprops))
C
      DO I=1, NPROPS
         IF(ICODES(1,I).EQ.4 .AND. ICODES(2,I).EQ.ICODES(6,I) .AND. 
     &      ICODES(3,I).EQ.0 .AND. ICODES(4,I).EQ.0 .AND. ICODES(5,I)
     &      .EQ.NUCCEN .AND. ICODES(7,I).LT.2)THEN
            MDIFF1=ABS(MGVN(ISTATE(I))-MGVN(JSTATE(I)))
            MDIFF2=MGVN(ISTATE(I))+MGVN(JSTATE(I))
            IF(MDIFF1.EQ.ICODES(7,I) .OR. MDIFF2.EQ.ICODES(7,I))THEN
               ICOUNT=ICOUNT+1
               INDEX1(ICOUNT)=I
               IF(ZISWIS0 .AND. ISTATE(I).EQ.JSTATE(I))THEN
                  ICOUNT2=ICOUNT2+1
                  INDEX2(ICOUNT2)=I
               END IF
            END IF
         END IF
      END DO
C
C---- Must now loop over the external region moments and convert the
C     moments to be inclusive of the nuclear term - if requested
C     by the user (ISW=0). The following is highly specific to
C     diatomic molecules in which one has
C
C       Nucleus, Nucleus, Scattering center
C
C     specified all the way through from INTS.
C
      IF(ZISWIS0)THEN
c
C     Get nuclear data (this has been done before, but its easier to
C     repeat the calculation than to store the output)
c
         noadc=0
         DO I=1, NNUC
            vlist(I,1)=RGEOM(1,I)
            vlist(I,2)=RGEOM(2,I)
            vlist(I,3)=RGEOM(3,I)
            vlist(I,4)=CHARG(I)
         END DO
C
         CALL INCENT(VLIST,vsymb,nnuc,NNUC,NOADC,rmoi,0)
C
         DO inuc=1, nnuc
            IF(inuc.EQ.nuccen)CYCLE
c
            IF(rgeom(1,inuc).NE.xzero .OR. rgeom(2,inuc).NE.xzero)THEN
               WRITE(iwrite,9940)
               STOP
            END IF
            R=RGEOM(3,inuc)-RGEOM(3,NUCCEN)
            Z=charg(inuc)
            mass(inuc)=vlist(inuc,6)
C
            icount1=0
            DO I=1, ICOUNT2
               k=index2(i)
               PROPS(k)=PROPS(k)-Z*R**ICODES(6,k)
               IF(abs(props(k)).LT.tol)CYCLE
               icount1=icount1+1
            END DO
         END DO
      END IF
C
      WRITE(iwrite,2111)nftmt
C
C=======================================================================
C
C     SET HEADER with KEY = 9 :
C     FIELD
C       2   SET NUMBER
C       3   NUMBER OF RECORDS IN SET
C       4   NUMBER OF GEOMETRY RECORDS ( NNUC )
C       5   NUMBER OF RECORDS OF TARGET DATA
C       6   NUMBER OF RECORDS OF MOMENT DATA
C       7   MOMENT TYPE SWITCH, ISW
C       8
C       9   INTERNUCLEAR SEPARATION, RAB  ( DIATOMICS ONLY )
C      10   36-CHARACTER HEADER FIELD
C
C=======================================================================
C
      RAB=XZERO
      mnuc=0
      DO I=1, NNUC
         RAB=RAB+ABS(RGEOM(3,I))
         ICHARG=INT(ABS(CHARG(I)))
         IF(I.EQ.NUCCEN .AND. icharg.EQ.0)CYCLE
         mnuc=mnuc+1
      END DO
C
      NRECS=NTARG+mNUC+ICOUNT1
C
      WRITE(NFTMT,1100)KEYH, NMSET, NRECS, mNUC, NTARG, ICOUNT1, ISW, 
     &                 NULL, RAB, RNAME
      WRITE(iwrite,2100)KEYH, NMSET, NRECS, mNUC, NTARG, ICOUNT1, ISW, 
     &                  NULL, RAB, RNAME
C
C=======================================================================
C
C     DATA DEFINING MOLECULAR GEOMETRY with KEY = 8 :
C     FIELD
C       2   NUCLEAR SEQUENCE NUMBER (I)
C       3   SET TO 0 FOR A REAL TARGET NUCLEUS, TO 1 FOR THE SCATTERING
C       4   NUCLEAR CHARGE ( USE NAMELIST INPUT FOR NONINTEGRAL VALUES )
C       5   NUCLEAR MASS ( IN ATOMIC UNITS )
C       6
C       7
C       8
C       9   NUCLEAR POSITION, GEONUC(I)
C      10   26-CHARACTER HEADER FIELD
C
C=======================================================================
C
      DO I=1, NNUC
         ICHARG=INT(ABS(CHARG(I)))
         IF(I.NE.NUCCEN)THEN
            ITYPE=0
            WRITE(NFTMT,1100)KEYG, I, ITYPE, ICHARG, MASS(I), NULL, 
     &                       NULL, NULL, RGEOM(3,I), RNAME
            WRITE(iwrite,2100)KEYG, I, ITYPE, ICHARG, MASS(I), NULL, 
     &                        NULL, NULL, RGEOM(3,I), RNAME
         END IF
      END DO
C
C=======================================================================
C
C     TARGET DATA with KEY = 5  :
C
C     FIELD
C       2   STATE INDEX         (I)
C       3   MANIFOLD INDEX
C       4   INDEX WITHIN MANIFOLD
C       5   |M|            MTARG(I)
C       6   2*S+1          STARG(I)
C       7   +/- INDEX      ITARG(I)
C       8   G/U INDEX      GTARG(I)
C       9   E IN AU        ETARG(I)
C
C=======================================================================
C
      DO I=1, NTARG
         RNAME=CMAKSTN(I,MGVN(I),ISPIN(I),GUTOT(I),IREFL(I),ISYMTYP,
     &         NSYM,ksym)
         WRITE(NFTMT,1100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I), 
     &                    IREFL(I), GUTOT(I), ENERGY(I), RNAME
         WRITE(iwrite,2100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I), 
     &                     IREFL(I), GUTOT(I), ENERGY(I), RNAME
      END DO
      IF(iqdfg.EQ.1)THEN
         WRITE(ordpot,2222)((ENERGY(I)-grounden)/2,i=1,ntarg)
      END IF
 2222 FORMAT(7(f12.4))
C
C=======================================================================
C
C     Properties records
C
C=======================================================================
!
      DO I=1, ICOUNT
         K=INDEX1(I)
         IF(abs(props(k)).LT.tol)CYCLE
         WRITE(NFTMT,1100)KEYP, ISTATE(K), MGVN(ISTATE(K)), JSTATE(K), 
     &                    MGVN(JSTATE(K)), (ICODES(J,K),J=5,7), PROPS(K)
         WRITE(iwrite,2100)KEYP, ISTATE(K), MGVN(ISTATE(K)), JSTATE(K), 
     &                     MGVN(JSTATE(K)), (ICODES(J,K),J=5,7), 
     &                     PROPS(K)
 
! write g/s dipole moment into borndat for borncros
         IF(ISTATE(K).EQ.JSTATE(K) .AND. ISTATE(K).EQ.1)THEN
            IF(ICODES(6,K).EQ.1 .AND. ICODES(7,K).EQ.0)THEN
               IF(PROPS(k).GE.0._wp)THEN
                  dip1=PROPS(k)
               ELSE
                  dip1=-1._wp*PROPS(k)
               END IF
            ELSE IF(ICODES(6,K).EQ.1 .AND. ICODES(7,K).EQ.1)THEN
               IF(PROPS(k).GE.0._wp)THEN
                  dip2=PROPS(k)
               ELSE
                  dip2=-1._wp*PROPS(k)
               END IF
            ELSE IF(ICODES(6,K).EQ.1 .AND. ICODES(7,K).EQ.-1)THEN
               IF(PROPS(k).GE.0._wp)THEN
                  dip3=PROPS(k)
               ELSE
                  dip3=-1._wp*PROPS(k)
               END IF
            END IF
         END IF
!
      END DO
!
      WRITE(dend,654)dip1
      WRITE(dend,654)dip2
      WRITE(dend,654)dip3
      CLOSE(dend)
 654  FORMAT(F20.10)
 
c
      DEALLOCATE(index1)
      DEALLOCATE(index2)
C
      RETURN
C
C---- Format Statement
C
 12   FORMAT(1X,I1,7I3,D20.12,2X,A26)
C
 1000 FORMAT(////,10X,' ====> Create/Update target property file <====',
     &       //)
 1010 FORMAT(/,10X,'Logical unit for target properties file   = ',I3,/,
     &       10X,'Set number at which to write target data  = ',I3,/,
     &       10X,'Format style for the target properties    = ',A,/)
 1020 FORMAT(/,10X,'No. of nuclear centers in molecular frame = ',I3,/,
     &       10X,'Sequence number of the scattering center  = ',I3,/,
     &       10X,'Number of target states in this dataset   = ',I3,/,
     &       10X,'The number of properties computed         = ',I3,/,
     &       10X,'Nuclear contribution switch (ISW) is set  = ',I3,/)
 1030 FORMAT(/,10X,'Nuclear Statistics Table: ',/,10X,
     &       '------------------------- ',//,10X,'  Charge ',1X,5X,'X',
     &       10X,'Y',10X,'Z')
 1040 FORMAT(10X,F8.3,1X,F10.5,1X,F10.5,1X,F10.5)
C
 1100 FORMAT(I1,7I3,D20.12,2X,A36)
 1110 FORMAT(A)
 2111 FORMAT(/' Output to unit ',i0,' is as follows:'/)
 2100 FORMAT(1x,I1,7I3,D20.12,2X,A35)
C
 9900 FORMAT(/10X,'**** Error in WRTARG: ',/)
 9910 FORMAT(10X,'Either, or both NNUC and NTARG are less than zero:',/,
     &       10X,'NNUC = ',I5,' NTARG= ',I5,/)
 9920 FORMAT(10X,'ISW should be 0 or 1. It is input as = ',I10,/)
 9930 FORMAT(10X,'NUCCEN must lie in the range 0 < NUCCEN <= ',I3,//,
     &       10X,'NUCCEN has been input as = ',I5,/)
 9940 FORMAT(10X,'Routine has been hardwired for linear case')
C
      END SUBROUTINE WRDTARG
!*==wrptarg.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRPTARG(NFTMT,NMSET,NNUC,NTARG,NPROPS,ISW,NUCCEN,RGEOM,
     &                   CHARG,ISPIN,MGVN,ENERGY,ISTATE,JSTATE,ICODES,
     &                   PROPS,IWRITE,ISYMTYP,NSYM,iqdfg,grounden,ksym,
     &                   ukrmolp_ints)
C***********************************************************************
C
C     WRPTARG - WRite a TARGet properties file for non-linear molecule
C
C     Linkage:
C
C        GETSET,CMAKSTN
C
C     Notes:
C
C     This is the partner of the routine RDTARG in the external region
C     codes. The properties file is always formatted and all records
C     have a fixed length of 80 bytes. They are written (and read) with
C     the format
C
C                   I1, 7I3, D20.12, 2X, A36
C
c     ICODES(1,*) = 4
c     ICODES(2,*) = power of r
c     ICODES(3,*) = power of sin(theta)
c     ICODES(4,*) = power of cos(theta)
c     ICODES(5,*) = centre
c     ICODES(6,*) = l of associated Legendre P_lm (cos(theta))
c     ICODES(7,*) = m of associated Legendre P_lm (cos(theta))
c     ICODES(8,*) = q
c
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      USE ukrmol_interface_gbl, only: GET_ECP_CORE_CHARGES
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: KEYH=6, KEYG=8, KEYT=5, KEYP=1, NULL=0
      CHARACTER(LEN=11), PARAMETER :: TFORM='FORMATTED  '
C
C Dummy arguments
C
      REAL(KIND=wp) :: GROUNDEN
      INTEGER :: IQDFG, ISW, ISYMTYP, IWRITE, KSYM, NFTMT, NMSET, NNUC, 
     &           NPROPS, NSYM, NTARG, NUCCEN
      REAL(KIND=wp), DIMENSION(*) :: CHARG, ENERGY, PROPS
      INTEGER, DIMENSION(8,*) :: ICODES
      INTEGER, DIMENSION(*) :: ISPIN, ISTATE, JSTATE, MGVN
      REAL(KIND=wp), DIMENSION(3,*) :: RGEOM
      LOGICAL :: ukrmolp_ints
      INTENT (IN) CHARG, ENERGY, GROUNDEN, IQDFG, ISTATE, ISW, IWRITE, 
     &            JSTATE, NPROPS, NTARG, NUCCEN, RGEOM, ukrmolp_ints
      INTENT (INOUT) ICODES, PROPS
C
C Local variables
C
      CHARACTER(LEN=36), EXTERNAL :: CMAKSTN
      CHARACTER(LEN=8) :: CSYMBOL
      INTEGER :: DEND, I, ICHARG, ICOUNT, ICOUNT2, INUC, J, K, L, M, 
     &           NNUCT, NOADC, NRECS, ORDPOT=55, Q
      REAL(KIND=wp) :: DIP1, DIP2, DIP3, GS, RMASS, X, Y, Z, ZZ
      INTEGER, DIMENSION(nprops) :: INDEX1, INDEX2
      INTEGER :: INT, LEN
      REAL(KIND=wp), DIMENSION(3) :: RMOI
      CHARACTER(LEN=36) :: RNAME
      LOGICAL :: SCATI
      REAL(KIND=wp), DIMENSION(nnuc,8) :: VLIST
      CHARACTER(LEN=8), DIMENSION(nnuc) :: VSYMB
      INTEGER, ALLOCATABLE :: core_charges(:)
C
C*** End of declarations rewritten by SPAG
C
      DATA dend/288/, gs/0._wp/
C
      WRITE(IWRITE,1000)
      WRITE(IWRITE,1010)NFTMT, NMSET, TFORM
      WRITE(IWRITE,1020)NNUC, NUCCEN, NTARG, NPROPS, ISW
c
      scati=.FALSE.
      DO I=1, NNUC
         IF(i.EQ.nuccen)scati=.TRUE.
      END DO
C
C---- Initialize the output character string to blanks
C
      DO I=1, LEN(RNAME)
         RNAME(I:I)=' '
      END DO
C
C---- First position NFTMT at the begining of set NTSET
C
      CALL GETSET(NFTMT,NMSET,KEYH,TFORM)
C
C---- Must scan the properties records an extract only those that are
C     going to be needed in an external region run. They are
C
C     r**l S_{lmq}   centered on the NUCCEN nucleus
C
C     Also we must make sure that delta lamda matches for the target
C     state pair and the operator.
C
C     An index vector is built which is later used to extract the
C     elements from the appropriate arrays; a second index vector
C     corresponding to the sub set of these selected moments having
C     target states with themselves is also built. The second vector
C     is used to subtract nuclear moment contributions.
C
      ICOUNT=0
      ICOUNT2=0
C
      DO I=1, NPROPS
         IF(ICODES(1,I).EQ.4 .AND. ICODES(2,I).EQ.ICODES(6,I) .AND. 
     &      ICODES(3,I).EQ.0 .AND. ICODES(4,I).EQ.0 .AND. ICODES(5,I)
     &      .EQ.NUCCEN)THEN
            ICOUNT=ICOUNT+1
            INDEX1(ICOUNT)=I
            IF(ISW.EQ.0 .AND. ISTATE(I).EQ.JSTATE(I))THEN
               ICOUNT2=ICOUNT2+1
               INDEX2(ICOUNT2)=I
            END IF
         END IF
      END DO
C
C     Add the nuclear term - if requested (ISW=0).
c
      IF(ISW.EQ.0)THEN

         if (ukrmolp_ints) then
            CALL GET_ECP_CORE_CHARGES(core_charges)
         else
            allocate(core_charges(nnuc))
            core_charges = 0
         endif

         DO inuc=1, nnuc
            IF(inuc.EQ.nuccen)CYCLE
c
            x=RGEOM(1,inuc)-RGEOM(1,NUCCEN)
            y=RGEOM(2,inuc)-RGEOM(2,NUCCEN)
            z=RGEOM(3,inuc)-RGEOM(3,NUCCEN)
            ZZ=charg(inuc)-core_charges(inuc)
C
            DO I=1, ICOUNT2
               l=ICODES(6,INDEX2(I))
               m=ICODES(7,INDEX2(I))
               q=ICODES(8,INDEX2(I))
c
               IF(l.EQ.1)THEN
                  IF(m.EQ.0)THEN
                     PROPS(INDEX2(I))=PROPS(INDEX2(I))-ZZ*z
                  ELSE IF(q.EQ.1)THEN
                     PROPS(INDEX2(I))=PROPS(INDEX2(I))-ZZ*x
                  ELSE
                     PROPS(INDEX2(I))=PROPS(INDEX2(I))-ZZ*y
                  END IF
               ELSE IF(l.EQ.2)THEN
                  IF(m.EQ.0)THEN
                     PROPS(INDEX2(I))=PROPS(INDEX2(I))
     &                                -0.5_wp*zz*(2._wp*z*z-x*x-y*y)
                  ELSE IF(m.EQ.1)THEN
                     IF(q.EQ.1)THEN
                        PROPS(INDEX2(I))=PROPS(INDEX2(I))-sqrt(3._wp)
     &                                   *zz*x*z
                     ELSE
                        PROPS(INDEX2(I))=PROPS(INDEX2(I))-sqrt(3._wp)
     &                                   *zz*y*z
                     END IF
                  ELSE
                     IF(q.EQ.1)THEN
                        PROPS(INDEX2(I))=PROPS(INDEX2(I))
     &                                  -0.5_wp*sqrt(3._wp)*zz*(x*x-y*y)
                     ELSE
                        PROPS(INDEX2(I))=PROPS(INDEX2(I))-sqrt(3._wp)
     &                                   *zz*x*y
                     END IF
                  END IF
               ELSE
                  WRITE(iwrite,9900)
                  WRITE(iwrite,9910)l, m, q
                  STOP 999
               END IF
            END DO
         END DO
      END IF
c
C     Get nuclear data (this has been done before, but its easier to
C     repeat the calculation than to store the output)
c
      noadc=0
      DO I=1, NNUC
         vlist(I,1)=RGEOM(1,I)
         vlist(I,2)=RGEOM(2,I)
         vlist(I,3)=RGEOM(3,I)
         vlist(I,4)=CHARG(I)
      END DO
C
      CALL INCENT(VLIST,vsymb,nnuc,NNUC,NOADC,rmoi,0)
C
C=======================================================================
C
C     SET HEADER with KEY = 6 :
C     FIELD
C       2   SET NUMBER
C       3   NUMBER OF RECORDS IN SET
C       4   NUMBER OF GEOMETRY RECORDS ( NNUC )
C       5   NUMBER OF RECORDS OF TARGET DATA
C       6   NUMBER OF RECORDS OF MOMENT DATA
C       7   MOMENT TYPE SWITCH, ISW
C       8   Rotational constant AX (a.u.)
C       9   Rotational constant BY (a.u.)
C      10   Rotational constant CZ (a.u.)
C
C=======================================================================
C
C --- If the input is from a scattering run, then adjust nnuc
      IF(scati)THEN
         nnuct=nnuc-1
      ELSE
         nnuct=nnuc
      END IF
c
      WRITE(iwrite,2111)nftmt
c
      NRECS=NTARG+NNUC+ICOUNT
      WRITE(NFTMT,1101)KEYH, NMSET, NRECS, NNUCt, NTARG, ICOUNT, ISW, 
     &                 rmoi
      WRITE(iwrite,2101)KEYH, NMSET, NRECS, NNUCt, NTARG, ICOUNT, ISW, 
     &                  rmoi
C
C=======================================================================
C
C     DATA DEFINING MOLECULAR GEOMETRY with KEY = 8 :
C     FIELD
C       2   NUCLEAR SEQUENCE NUMBER (I)
C       3   Name of atom
C       4   NUCLEAR CHARGE
C       5   NUCLEAR MASS
C       6   Nuclear x-coordinate
C       7   Nuclear y-coordinate
C       8   Nuclear z-coordinate
C
C=======================================================================
C
      DO I=1, NNUC
         IF(I.NE.NUCCEN)THEN
            csymbol=vsymb(i)
            rmass=vlist(i,6)
            ICHARG=INT(ABS(CHARG(I)))
            WRITE(NFTMT,1102)KEYG, I, csymbol(6:), ICHARG, rMASs, 
     &                       (RGEOM(j,I),j=1,3)
            WRITE(iwrite,2102)KEYG, I, csymbol(6:), ICHARG, rMASs, 
     &                        (RGEOM(j,I),j=1,3)
         END IF
      END DO
C
C=======================================================================
C
C     TARGET DATA with KEY = 5  :
C
C     FIELD
C       2   STATE INDEX         (I)
C       3
C       4
C       5   Symmetry number (Alchemy style) MTARG(I)
C       6   2*S+1          STARG(I)
C       7
C       8
C       9   E IN AU        ETARG(I)
C
C=======================================================================
C
      DO I=1, NTARG
         RNAME=CMAKSTN(I,MGVN(I),ISPIN(I),null,null,ISYMTYP,NSYM,ksym)
         if (ntarg .lt. 100) then 
            WRITE(NFTMT,1100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I), 
     &                     null, null, ENERGY(I), RNAME
            WRITE(iwrite,2100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I),
     &                      null, null, ENERGY(I), RNAME
         else
            WRITE(NFTMT,3100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I), 
     &                     null, null, ENERGY(I), RNAME
            WRITE(iwrite,5100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I),
     &                      null, null, ENERGY(I), RNAME
         end if
c~          WRITE(NFTMT,1100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I), null, 
c~      &                    null, ENERGY(I), RNAME
c~          WRITE(iwrite,2100)KEYT, I, NULL, NULL, MGVN(I), ISPIN(I), null, 
c~      &                     null, ENERGY(I), RNAME
      END DO
      IF(iqdfg.EQ.1)THEN
         WRITE(ordpot,2222)((ENERGY(I)-grounden)/2,i=1,ntarg)
      END IF
 2222 FORMAT(7(f12.4))
C
C=======================================================================
C
C     Properties records
C
C=======================================================================
C
C ZM: zero dip1,dip2,dip3 first
      dip1 = 0.0_wp
      dip2 = 0.0_wp
      dip3 = 0.0_wp
      DO I=1, ICOUNT
         K=INDEX1(I)
         IF(icodes(8,k).EQ.-1)icodes(7,k)=-icodes(7,k)
         if (ntarg .lt. 100) then
            WRITE(NFTMT,1100)KEYP, ISTATE(K), MGVN(ISTATE(K)), 
     &                   JSTATE(K),MGVN(JSTATE(K)), (ICODES(J,K),J=5,7),
     &                   PROPS(K)
            WRITE(iwrite,2100)KEYP, ISTATE(K),MGVN(ISTATE(K)),JSTATE(K), 
     &                    MGVN(JSTATE(K)), (ICODES(J,K),J=5,7), 
     &                    PROPS(K)
         else
            WRITE(NFTMT,4100)KEYP, ISTATE(K), MGVN(ISTATE(K)), 
     &                   JSTATE(K),MGVN(JSTATE(K)), (ICODES(J,K),J=5,7),
     &                   PROPS(K)
            WRITE(iwrite,6100)KEYP, ISTATE(K),MGVN(ISTATE(K)),JSTATE(K), 
     &                    MGVN(JSTATE(K)), (ICODES(J,K),J=5,7), 
     &                    PROPS(K)

         end if
c~          WRITE(NFTMT,1100)KEYP, ISTATE(K), MGVN(ISTATE(K)), JSTATE(K), 
c~      &                    MGVN(JSTATE(K)), (ICODES(J,K),J=5,7), PROPS(K)
c~          WRITE(iwrite,2100)KEYP, ISTATE(K), MGVN(ISTATE(K)), JSTATE(K), 
c~      &                     MGVN(JSTATE(K)), (ICODES(J,K),J=5,7), 
c~      &                     PROPS(K)
 
!
! write g/s dipole moment into borndat for borncros
         IF(ISTATE(K).EQ.JSTATE(K) .AND. ISTATE(K).EQ.1)THEN
            IF(ICODES(6,K).EQ.1 .AND. ICODES(7,K).EQ.0)THEN
               IF(PROPS(k).GE.0._wp)THEN
                  dip1=PROPS(k)
               ELSE
                  dip1=-1._wp*PROPS(k)
               END IF
            ELSE IF(ICODES(6,K).EQ.1 .AND. ICODES(7,K).EQ.1)THEN
               IF(PROPS(k).GE.0._wp)THEN
                  dip2=PROPS(k)
               ELSE
                  dip2=-1._wp*PROPS(k)
               END IF
            ELSE IF(ICODES(6,K).EQ.1 .AND. ICODES(7,K).EQ.-1)THEN
               IF(PROPS(k).GE.0._wp)THEN
                  dip3=PROPS(k)
               ELSE
                  dip3=-1._wp*PROPS(k)
               END IF
            END IF
         END IF
!
      END DO
!
      WRITE(dend,654)dip1
      WRITE(dend,654)dip2
      WRITE(dend,654)dip3
      CLOSE(dend)
 654  FORMAT(F20.10)
c
      RETURN
C
C---- Format Statement
 6    FORMAT('In subroutine wrptarg ksym is',(i2))
 1000 FORMAT(////,10X,' ====> Create/Update target property file <====',
     &       //)
 1010 FORMAT(/,10X,'Logical unit for target properties file   = ',I3,/,
     &       10X,'Set number at which to write target data  = ',I3,/,
     &       10X,'Format style for the target properties    = ',A,/)
 1020 FORMAT(/,10X,'No. of nuclear centers in molecular frame = ',I3,/,
     &       10X,'Sequence number of the scattering center  = ',I3,/,
     &       10X,'Number of target states in this dataset   = ',I3,/,
     &       10X,'The number of properties computed         = ',I6,/,
     &       10X,'Nuclear contribution switch (ISW) is set  = ',I3,/)
C
 1100 FORMAT(I1,7I3,D20.12,2X,A36)
 1101 FORMAT(I1,I3,I6,I3,I4,I6,I3,1x,3D20.12)
 1102 FORMAT(I1,I3,a3,i3,f10.4,3F20.10)
 2111 FORMAT(/' Output to unit ',i0,' is as follows:'/)
 2100 FORMAT(I1,7I3,D20.12,2X,A36)
 2101 FORMAT(I1,I3,I6,I3,I4,I6,I3,1x,3D20.12)
 2102 FORMAT(I1,I3,a3,i3,f10.4,3F20.10)
 3100 format(I1,i8,6I3,D20.12,2X,A36) !A.Harvey 2012 for nstat > 99
 4100 format(I1,I8,I3,I8,4I3,D20.12,2X,A35) !A.Harvey 2012 for nstat > 99
 5100 format(1x,I1,i8,6I3,D20.12,2X,A36) !A.Harvey 2012 for nstat > 99
 6100 format(1x,I1,I8,I3,I8,4I3,D20.12,2X,A35) !A.Harvey 2012 for nstat > 99
C
 9900 FORMAT(/10X,'**** Error in WRTARG: ',/)
 9910 FORMAT(10X,' Error in property operator code ',3I5)
c
      END SUBROUTINE WRPTARG
!*==wrtda.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTDA(NFTX,JDSK,NCHA,LBOX,NEL,MI,MJ,MPQ,DPQ)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: JDSK, LBOX, NCHA, NEL, NFTX
      REAL(KIND=wp), DIMENSION(LBOX) :: DPQ
      INTEGER, DIMENSION(LBOX) :: MI, MJ, MPQ
      INTENT (IN) DPQ, JDSK, LBOX, MI, MJ, MPQ, NCHA, NEL, NFTX
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     WRITES OUT A DIRECT ACCESS RECORD ONTO THE FILE NFTX. THIS IS
C     A RECORD OF INDEXED MOMENT EXPRESSIONS.
C
C***********************************************************************
C
C
      WRITE(NFTX,REC=JDSK)NCHA, NEL, MI, MJ, MPQ, DPQ
C
      RETURN
C
      END SUBROUTINE WRTDA
!*==wrtdb.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTDB(NFTO,LBOX,NEL,MI,MJ,MPQ,DPQ)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LBOX, NEL, NFTO
      REAL(KIND=wp), DIMENSION(LBOX) :: DPQ
      INTEGER, DIMENSION(LBOX) :: MI, MJ, MPQ
      INTENT (IN) DPQ, LBOX, MI, MJ, MPQ, NEL, NFTO
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     WRTDB - WRITES OUT WHAT IS ALREADY STORED IN BOXES ONTO THE
C             SEQUENTIAL FILE DURING THE SORT PROCEDURE. THIS IS THE
C             DATA THAT WAS NOT WRITTEN ONTO THE DIRECT ACCESS FILE
C
C***********************************************************************
C
C
      WRITE(NFTO)NEL, MI, MJ, MPQ, DPQ
 
C
      RETURN
      END SUBROUTINE WRTDB
!*==wrtdii.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTDII(NEL,NN,LNGTH,MP,NFD,NELM,NEMB,JROB,JRON,IWRITE)
C***********************************************************************
C
C     WRTDII - WRITE DIAGONAL ELEMENTS.
C
C              WRITES THE FORMULAE FOR THE DIAGONAL DENSITY MATRIX
C              ELEMENTS TO THE FILE ON UNIT NFD. ACTUALLY THIS ROUTINE
C              STORES THE FORMULAE IN A BUFFER UNTIL IT IS FULL. THE
C              BUFFER IS EMPTIED EACH TIME IT BECOMES FULL. A FINAL
C              CALL IS REQUIRED TO SUBROUTINE WRTDIX TO EMPTY THE
C              BUFFER AFTER THE DIAGONAL CALCULATION IS COMPLETE.
C
C     INPUT DATA :
C            NEL  NUMBER OF FORMULAE ELEMENTS ALREADY IN BUFFER
C             NN  NUMBER OF FORMULAE ELEMENTS TO BE ADDED
C          LNGTH  SIZE OF THE BUFFER
C             MP  BUFFER AREA
C            NFD  LOGICAL UNIT TO WHICH BUFFER IS WRITTEN
C           NELM  TOTAL NUMBER OF FORMULAE ELEMENTS ALREADY WRITTEN
C           NEMB  TOTAL NUMBER OF TIMES THE BUFFER HAS BEEN EMPTIED SO
C                 FAR. AS EACH WRITE PRODUCES ONLY ONE RECORD THEN THIS
C                 COUNTS THE NUMBER OF RECORDS WRITTEN.
C           JROB  ORBITAL NUMBERS. THIS IS PART OF THE FORMULAE
C           JRON  ORBITAL OCCUPATION NUMBERS. THE OTHER PART OF THE
C                 FORMULAE.
C          IWRITE LOGICAL UNIT FOR THE PRINTER
C
C   OUTPUT DATA :
C                 ON OUTPUT NELM,NEMB,NEL AND MP ARE UPDATED
C
C     Notes:
C
C       The routine writes full records to the disk, the length of each
C     being
C
C           (1 + LNGTH) integer words
C
C     Thus LNGTH, an input parameter to the code, should be chosen to
C     optimize the disk I/O transfer rate; this means that full buffers
C     should be written out.
C
C     Some possible values for LNGTH and associated record size are:
C
C     Cray Unicos - set LNGTH = 4095 giving 4096 words (32768 bytes)
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, LNGTH, NEL, NELM, NEMB, NFD, NN
      INTEGER, DIMENSION(*) :: JROB, JRON
      INTEGER, DIMENSION(LNGTH) :: MP
      INTENT (IN) IWRITE, JROB, JRON, LNGTH, NFD, NN
      INTENT (INOUT) MP, NEL, NELM, NEMB
C
C Local variables
C
      INTEGER :: I, ICOUNTNZ, ISTART
C
C*** End of declarations rewritten by SPAG
C
C-----IF (2*NN+1), WHICH ARE THE NUMBER OF ELEMENTS TO BE ADDED, IS
C     GREATER THAN THE BUFFER SPACE AVAILABLE THEN CALL TMTCLOS() THE CO
C     AND WARN USER.
C
      IF((2*NN+1).GT.LNGTH)THEN
         WRITE(IWRITE,9010)
         WRITE(IWRITE,9020)LNGTH, (2*NN+1)
         CALL TMTCLOS()
      END IF
C
C-----ADD THE NEW ELEMENTS TO THE BUFFER MP,UNLESS THERE IS NOT
C     ENOUGH SPACE LEFT. IN THAT CASE EMPTY THE BUFFER ONTO DISK
C     FIRST AND THEN ADD NEW ELEMENTS. NEMB COUNTS THE NUMBER OF
C     RECORDS AND NELM THE NUMBER OF ELEMENTS.
C
      IF((NEL+2*NN+1).GT.LNGTH)THEN
         WRITE(NFD)NEL, MP
         NELM=NELM+NEL
         NEMB=NEMB+1
         NEL=0
      END IF
C
C...... Move to the next free storage location. It will contain the
C       number of elements that are stored. However we do not know
C       that until we have excluded the zero ones during storage.
C       Keep the position of this free location in ISTART.
C
      NEL=NEL+1
      ISTART=NEL
C
C...... Loop over all the elements and add the non-zero ones only.
C       i.e. non zero occupation numbers JRON.
C
      ICOUNTNZ=0
C
      DO I=1, NN
         IF(JRON(I).NE.0)THEN
            ICOUNTNZ=ICOUNTNZ+1
            NEL=NEL+1
            MP(NEL)=JROB(I)
            NEL=NEL+1
            MP(NEL)=JRON(I)
         END IF
      END DO
C
C...... Go back and fill in position MP(ISTART) now
C
      MP(ISTART)=ICOUNTNZ
C
      RETURN
C
C---- Format Statements
C
 9010 FORMAT(//,10X,'**** Error in WRTDII: ',//)
 9020 FORMAT(10X,'Insufficient space in buffer for output of diagonal',
     &       ' density matrix elements',//,10X,'Given space = ',I10,
     &       ' Need = ',I10,/)
C
      END SUBROUTINE WRTDII
!*==wrtdix.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTDIX(NEL,LNGTH,MP,NFD,NELM,NEMB,IWRITE,NPFLG)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, LNGTH, NEL, NELM, NEMB, NFD, NPFLG
      INTEGER, DIMENSION(LNGTH) :: MP
      INTENT (IN) IWRITE, LNGTH, MP, NFD, NPFLG
      INTENT (INOUT) NEL, NELM, NEMB
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     WRTDIX - WRITE THE FINAL ELEMENTS OF THE DENSITY MATRIX
C              EXPRESSIONS.
C
C     WRITE OUT THE LAST ELEMENTS OF THE DENSITY MATRIX EXPRESSIONS
C     TO THE FILE NFD AND WRITE MESSAGE TO THE OUTPUT STREAM.
C     INPUT PARAMETERS ARE THE SAME AS FOR INITIAL CALL.
C
C***********************************************************************
C
 
C
      NEMB=NEMB+1
      NELM=NELM+NEL
      NEL=-NEL
      WRITE(NFD)NEL, MP
C
      IF(NPFLG.GT.0)WRITE(IWRITE,20)NELM, NEMB
C
      RETURN
C
C---- Format Statements
C
 20   FORMAT(//,5X,'The final diagonal formulae elements have been ',
     &       'written out.',//,5X,
     &       'The following parameters summarize the process: ',/,5X,
     &       'Coefficients written = ',I10,/,5X,
     &       'Records required     = ',I10,/)
C
      END SUBROUTINE WRTDIX
!*==wrtmt.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTMT(LUMOM,ECI,ECJ,ISTATE,JSTATE,LAMDAI,LAMDAJ,IOPCDS,
     &                 XMOM,NMOMS)
C***********************************************************************
C
C     WRTMT - WRite Transition MomenTs in old style
C
C     Input data:
C
C     Linkage:
C
C         IPACK9
C
C     Author: Charles J Gillan, March 1994
C
C     Copyright 1994 (c) Charles J Gillan
C     All rights reserved
C
C***********************************************************************
      USE integer_packing, ONLY:IPACK9
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: ECI, ECJ
      INTEGER :: ISTATE, JSTATE, LAMDAI, LAMDAJ, LUMOM, NMOMS
      INTEGER, DIMENSION(8,*) :: IOPCDS
      REAL(KIND=wp), DIMENSION(*) :: XMOM
      INTENT (IN) ECI, ECJ, IOPCDS, LUMOM, NMOMS, XMOM
C
C Local variables
C
      REAL(KIND=wp) :: EDIFF
      INTEGER :: I, IOUTW1, IOUTW2, K
      INTEGER, DIMENSION(7) :: IOPLOC
C
C*** End of declarations rewritten by SPAG
C
c        WRITE(IWRITE,1000)
c        WRITE(IWRITE,1010) LUMOM,ISTATE,JSTATE,ECI,ECJ
c        WRITE(IWRITE,1020) LAMDAI,LAMDAJ
c        WRITE(IWRITE,1030) NMOMS
c        WRITE(IWRITE,1040) ((IOPCDS(J,I),J=1,7),XMOM(1,I),I=1,NMOMS)
C
C---- Write the state definition information at the head of the moment
C     records
C
      EDIFF=ECJ-ECI
C
      WRITE(LUMOM)ISTATE, JSTATE, EDIFF, nmoms
C
C---- Begin loop over all moment values and
C
C       (i)  pack the operator and state to state information
C      (ii)  write the data to disk
C
C     The style of packing cannot handle operator codes with negative
C     numbers. These are never needed in R-matrix anyway at present
C     so they are simply set to zeros !
C
      DO K=1, NMOMS
         DO I=1, 7
            IOPLOC(I)=IOPCDS(I,K)
         END DO
         IOUTW1=IPACK9(IOPLOC,ISTATE,LAMDAI,JSTATE,LAMDAJ,IOUTW2)
         WRITE(LUMOM)IOUTW2, IOUTW1, XMOM(K)
      END DO
C
      RETURN
C
C---- Format Statements
C
 1000 FORMAT(//,25X,'====> WRTMT - OLD STYLE TMT OUTPUT <====',//)
 1010 FORMAT(/,25X,'Data will be written to logical unit = ',I10,/,25X,
     &       'ISTATE number                        = ',I10,/,25X,
     &       'JSTATE number                        = ',I10,/,25X,
     &       'Energy in Hartrees of I state        = ',F10.5,/,25X,
     &       'Energy in Hartrees of J state        = ',F10.5)
 1020 FORMAT(/,25X,'Lambda value of state I              = ',I10,/,25X,
     &       'Lambda value of state J              = ',I10)
 1030 FORMAT(/,25X,'Number of moment values = ',I5,//,25X,
     &       ' Operator codes       Moment Value ',/)
 1040 FORMAT(25X,8I3,1X,F20.10)
C
      END SUBROUTINE WRTMT
!*==wrtmth.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTMTH(LUMOM,NMSET,nrec,NOPREC,NSTATI,NSTATJ,ISYM,
     &                  NUCCEN,NPT,IOPCDS,NORB,NSRB,NELT,NSYM,ISYMTYP,
     &                  ISPIN,NOB,NNUC,GEONUC,CHARG,NAMEI,MGVNI,IREFLI,
     &                  GUTOTI,NCSFI,NAMEJ,MGVNJ,IREFLJ,GUTOTJ,NCSFJ)
C***********************************************************************
C
C     WRTMTH - WRites Transition MomenT Header in the old style used
C              in the TMTCJG and TMTJT series of codes. This routine is
C              provided because the users require backwards
C              compatibility. The partner routine is RDTMTH in the
C              external region modules of the R-matrix suite.
C
C     Author: Charles J Gillan, Queen's University Belfast, March 1994
C
C     Copyright 1994 (c) Charles J Gillan
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: KEYT=50
      CHARACTER(LEN=11), PARAMETER :: TFORM='UNFORMATTED'
C
C Dummy arguments
C
      INTEGER :: GUTOTI, GUTOTJ, IREFLI, IREFLJ, ISPIN, ISYM, ISYMTYP, 
     &           LUMOM, MGVNI, MGVNJ, NCSFI, NCSFJ, NELT, NMSET, NNUC, 
     &           NOPREC, NORB, NPT, NREC, NSRB, NSTATI, NSTATJ, NSYM, 
     &           NUCCEN
      CHARACTER(LEN=120) :: NAMEI, NAMEJ
      REAL(KIND=wp), DIMENSION(*) :: CHARG
      REAL(KIND=wp), DIMENSION(3,*) :: GEONUC
      INTEGER, DIMENSION(8,*) :: IOPCDS
      INTEGER, DIMENSION(*) :: NOB
      INTENT (IN) CHARG, GEONUC, GUTOTI, GUTOTJ, IOPCDS, IREFLI, IREFLJ, 
     &            ISPIN, ISYM, ISYMTYP, LUMOM, MGVNI, MGVNJ, NAMEI, 
     &            NAMEJ, NCSFI, NCSFJ, NELT, NMSET, NNUC, NOB, NOPREC, 
     &            NORB, NPT, NREC, NSRB, NSTATI, NSTATJ, NSYM, NUCCEN
C
C Local variables
C
      INTEGER :: I, IFAIL, J, NSET
C
C*** End of declarations rewritten by SPAG
C
      DATA nset/1/
C
C---- Position the dataset at its end now ready for writing.
C
      IFAIL=0
C
c      CALL GETSET(LUMOM,NSET,KEYT,TFORM)
c      nset = 0
C
C---- Header records for this output set of properties
C
      WRITE(LUMOM)KEYT, NMSET, nrec, NOPREC, NSTATI, NSTATJ, ISYM, 
     &            NUCCEN
C
      WRITE(LUMOM)NPT, ((IOPCDS(I,J),I=1,8),J=1,NPT)
C
      WRITE(LUMOM)NORB, NSRB, NELT, NSYM, ISYMTYP, ISPIN, 
     &            (NOB(I),I=1,NSYM), NNUC, 
     &            ((GEONUC(j,I),j=1,3),I=1,NNUC), (CHARG(I),I=1,NNUC)
C
      WRITE(LUMOM)NAMEI, MGVNI, IREFLI, GUTOTI, NCSFI
C
C..... A second wave function header is only needed in some cases !
C
      IF(ISYM.EQ.1)WRITE(LUMOM)NAMEJ, MGVNJ, IREFLJ, GUTOTJ, NCSFJ
C
      RETURN
C
C---- Format Statements
C
      END SUBROUTINE WRTMTH
!*==wrtodx.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTODX(NEL,LENGTH,II,JJ,NPQ,CPQ,NFD,NELMT,NBLK,IWRITE,
     &                  NPFLG)
C***********************************************************************
C
C     WRTODX - WRITE THE FINAL OFF-DIAGONAL ELEMENTS OF THE DENSITY
C              MATRIX EXPRESSIONS
C
C     WRITES OUT THE REMAINING CONTENTS OF THE BUFFERS FOR OFF-DIAGONAL
C     ELEMENTS AFTER THEIR COMPUTTATION HAS BEEN COMPLETED. OUTPUT GOES
C     TO FILE NFD. A MESSAGE IS WRITTEN TO THE PRINTER STREAM AS WELL.
C     INPUT PARAMETERS ARE SAME AS FOR WRTOFD.
C
C     Notes:
C
C       The routine writes full records to the disk, the length of each
C     being
C
C       (1 + 4*LENGTH) integer words + LENGTH double precision words
C
C     Thus LENGTH, an input parameter to the code, should be chosen to
C     optimize the disk I/O transfer rate; this means that full buffers
C     should be written out.
C
C     Some possible values for LNGTH and associated record size are:
C
C     Cray Unicos - set LNGTH = 819 giving 4096 words (32768 bytes)
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, LENGTH, NBLK, NEL, NELMT, NFD, NPFLG
      REAL(KIND=wp), DIMENSION(LENGTH) :: CPQ
      INTEGER, DIMENSION(LENGTH) :: II, JJ
      INTEGER, DIMENSION(2,LENGTH) :: NPQ
      INTENT (IN) CPQ, II, IWRITE, JJ, LENGTH, NFD, NPFLG, NPQ
      INTENT (INOUT) NBLK, NEL, NELMT
C
C Local variables
C
      INTEGER :: KK
C
C*** End of declarations rewritten by SPAG
C
      NELMT=NELMT+NEL
      NBLK=NBLK+1
      NEL=-NEL
      WRITE(NFD)NEL, II, JJ, NPQ, CPQ
C
      IF(npflg.GT.0)WRITE(IWRITE,20)NELMT, NBLK
C
      IF(NELMT.GT.1 .AND. NPFLG.GT.0)WRITE(IWRITE,21)
     &                                     (CPQ(KK),KK=1,NELMT)
C
      RETURN
C
C---- Format Statements
C
 20   FORMAT(/,5x,'Number of off-diagonal elements =',I8,5x,'in block',
     &       i5)
 21   FORMAT(/,5x,'Coefficients=',6F10.5/(18x,6F10.5))
      END SUBROUTINE WRTODX
!*==wrtofd.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE WRTOFD(NFD,NELMT,NBLK,NEL,LENGTH,ICSF,JCSF,COEFF,IORB1,
     &                  IORB2,II,JJ,NPQ,CPQ)
C***********************************************************************
C
C     WRTOFD - WRITE THE OFF-DIAGONAL DENSITY MATRIX FORMULAE TO DISC
C
C              STORES THE OFF-DIAGONAL DENSITY MATRIX FORMULAE IN TO
C              BUFFERS AND WHEN THESE ARE FULL WRITES THEM OUT TO DISC.
C              A FINAL CALL IS REQUIRED TO WRTODX TO WRITE THE FINAL
C              BUFFER TO DISC. THIS ROUTINE IS SIMILAR TO WRTDII FOR THE
C              DIAGONAL ELEMENTS, HOWEVER IN THE PRESENT CASE ONLY ONE
C              ELEMENT IS ADDED PER CALL. IT ASSUMED THAT THE SIZE OF
C              BUFFERS, GIVEN BY VARIABLE LENGTH, HAS BEEN SET TO A NON
C              ZERO VALUE.
C
C              IF THE VALUE OF COEFF IS LESS THAN A CERTAIN THRESHOLD
C              THEN THE FORMULA ELEMENT IS NOT STORED. THIS ELIMINATES
C              ANY ZERO ELEMENTS THAT ARE PASSED. THESE CAN OCCUR DUE TO
C              THE WAY IN WHICH SUBROUTINE ZERO WORKS. SUCH ZERO ELEMENT
C              ARE NOT INCORRECT BUT MERELY WASTE TIME.
C
C     INPUT DATA :
C            NFD  LOGICAL UNIT ON WHICH FORMULAE ARE WRITTEN
C          NELMT  COUNTS THE NUMBER OF FORMULAE ELEMENTS WRITTEN OUT
C                 ALREADY.
C           NBLK  COUNT NUMBERS OF TIMES THAT THE BUFFERS HAVE BEEN
C                 EMPTIED ALREADY. AS EACH WRITE OPERATION CREATES ONE
C                 RECORD THEN THIS COUNTS THE NUMBER OF RECORDS WRITTEN
C            NEL  NUMBER OF ELEMENTS ALREADY IN BUFFER
C         LENGTH  SIZE OF OFF-DIAGONAL BUFFERS
C           ICSF  FIRST CSF NUMBER
C           JCSF  SECOND CSF NUMBER
C          COEFF  COEFFICIENT
C          IORB1  FIRST ORBITAL INDEX
C          IORB2  SECOND ORBITAL INDEX
C             II  BUFFER FOR CSF NUMBERS
C             JJ  BUFFER FOR CSF NUMBERS
C            NPQ  BUFFER FOR ORBITAL INDICES
C            CPQ  BUFFER FOR COEFFICIENTS
C         IWRITE  LOGICAL UNIT FOR THE PRINTER
C
C     OUTPUT DATA :
C                   ON OUTPUT NELMT,NBLK,NEL ARE UPDATED
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
C     Note that THRESH defines the deletion threshold for density
C     matrix formulae.
      REAL(KIND=wp), PARAMETER :: THRESH=1.0E-08_wp
C
C Dummy arguments
C
      REAL(KIND=wp) :: COEFF
      INTEGER :: ICSF, IORB1, IORB2, JCSF, LENGTH, NBLK, NEL, NELMT, NFD
      REAL(KIND=wp), DIMENSION(LENGTH) :: CPQ
      INTEGER, DIMENSION(LENGTH) :: II, JJ
      INTEGER, DIMENSION(2,LENGTH) :: NPQ
      INTENT (IN) COEFF, ICSF, JCSF, LENGTH, NFD
      INTENT (INOUT) CPQ, II, IORB1, IORB2, JJ, NBLK, NEL, NELMT, NPQ
C
C Local variables
C
      INTEGER :: I, IHITROW, NTEMP
C
C*** End of declarations rewritten by SPAG
C
C---- Check the buffer space available. If there is not enough space
C     left then write out the buffer and update the pointers.
C
      IF(NEL+1.GT.LENGTH)THEN
         WRITE(NFD)NEL, II, JJ, NPQ, CPQ
         NELMT=NELMT+NEL
         NBLK=NBLK+1
         NEL=0
      END IF
C
C---- Place the orbitals in contiguous order
C
      IF(IORB2.GT.IORB1)THEN
         NTEMP=IORB1
         IORB1=IORB2
         IORB2=NTEMP
      END IF
C
C---- Add this data to the buffer but only if it passes the threshold
C     test.
C
      IF(ABS(COEFF).LT.THRESH)RETURN
C
C..... Let's scan the existing buffer to see if this element already
C      exists maybe - Then we just add the coefficient ! We set IHITROW
C      initially to zero and then scan the existing buffers for a hit
C      on the CSF and Orbital pairs.
C
      IHITROW=0
C
      DO I=1, NEL
         IF(II(I).EQ.ICSF .AND. JJ(I).EQ.JCSF .AND. NPQ(1,I)
     &      .EQ.IORB1 .AND. NPQ(2,I).EQ.IORB2)IHITROW=I
      END DO
C
C..... If we have a hit then just add the coefficient otherwise we must
C      enter a new element in the buffer.
C
C      When adding the new element there are two possibilities
C
C        (i) find a zero location and add over write it
C
C       (ii) add the element to the end of the buffer
C
      IF(IHITROW.NE.0)THEN
         CPQ(IHITROW)=CPQ(IHITROW)+COEFF
      ELSE
         IHITROW=0
         DO I=1, NEL
            IF(ABS(CPQ(I)).LT.THRESH)IHITROW=I
         END DO
C
         IF(IHITROW.NE.0)THEN
            II(IHITROW)=ICSF
            JJ(IHITROW)=JCSF
            NPQ(1,IHITROW)=IORB1
            NPQ(2,IHITROW)=IORB2
            CPQ(IHITROW)=COEFF
         ELSE
            NEL=NEL+1
            II(NEL)=ICSF
            JJ(NEL)=JCSF
            NPQ(1,NEL)=IORB1
            NPQ(2,NEL)=IORB2
            CPQ(NEL)=COEFF
         END IF
      END IF
C
      RETURN
C
      END SUBROUTINE WRTOFD
!*==zero.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE ZERO(NA,NB,NZ,NBL,MAXSO)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: MAXSO, NZ
      INTEGER, DIMENSION(*) :: NA, NB, NBL
      INTENT (IN) MAXSO, NA, NB, NBL
      INTENT (INOUT) NZ
C
C Local variables
C
      INTEGER :: I, K, M, MA, MAXM, MB
      INTEGER, DIMENSION(maxso) :: NBB
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     ZERO -   COMPARES THE FIRST PAIR OF DETERMINANTS FOR EACH CSF PAIR
C              AND RETURNS IN NZ THE NUMBER OF ORBITAL DIFFERENCES
C              BETWEEN THEM. THIS IS SIMILAR TO SUBROUTINE TO PZERO IN
C              PROGRAM SPEEDY.
C
C***********************************************************************
c
C
      DO I=1, MAXSO
         NBB(I)=0
      END DO
C
C-----FOR THE FIRST DETERMINANT OF THE FIRST WAVEFUNCTION, SCAN ALL
C     REPLACED ORBITALS AND DECREMENT POPULATION BY ONE. THEN
C     SCAN ALL OF THE REPLACEMENTS AND AUGUMENT POPULATION BY ONE.
C
      MA=NA(1)
      MB=NB(1)
c
      DO I=1, MA
         M=NBL(NA(I+1))
         NBB(M)=NBB(M)-1
         M=NBL(NA(MA+I+1))
         NBB(M)=NBB(M)+1
      END DO
C
C-----SCAN FIRST DETERMINANT OF THE SECOND WAVEFUNCTION AND APPLY INVERS
C     PROCESS TO THAT FOR FIRST WAVEFUNCTION.
C
      DO I=1, MB
         M=NBL(NB(I+1))
         NBB(M)=NBB(M)+1
         M=NBL(NB(MB+I+1))
         NBB(M)=NBB(M)-1
      END DO
C
C-----COUNT ALL OF THE AFFECTED ORBITALS AND FIND OVERALL CHANGE.
C
C
      NZ=0
      DO I=1, MA+MA
         M=NBL(NA(I+1))
C ZM maxm is not needed anywhere in this routine
C         maxm=max(maxm,m)
         K=NBB(M)
         NZ=NZ+ABS(K)
         NBB(M)=0
      END DO
C
      DO I=1, MB+MB
         M=NBL(NB(I+1))
         K=NBB(M)
         NZ=NZ+ABS(K)
         NBB(M)=0
      END DO
C
      RETURN
C
      END SUBROUTINE ZERO
!*==readcip.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE READCIP(nftw,nciset,nocsf,noveci,nctgt,irrep,s,sz,nelt,
     &                   EIG,VEC,iphase,cname,geonuc,charge,nfti,NFT,
     &                   itime,ifail)
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IFAIL, IRREP, ITIME, NCISET, NCTGT, NELT, NFT, NFTI, 
     &           NFTW, NOCSF, NOVECI
      REAL(KIND=wp) :: S, SZ
      REAL(KIND=wp), DIMENSION(*) :: CHARGE
      CHARACTER(LEN=8), DIMENSION(*) :: CNAME
      REAL(KIND=wp), DIMENSION(nocsf) :: EIG
      REAL(KIND=wp), DIMENSION(3,*) :: GEONUC
      INTEGER, DIMENSION(nocsf) :: IPHASE
      REAL(KIND=wp), DIMENSION(nocsf*noveci) :: VEC
      INTENT (IN) IRREP, ITIME, NCISET, NCTGT, NFT, NOCSF, NOVECI
      INTENT (OUT) NELT, S, SZ
      INTENT (INOUT) CHARGE, CNAME, EIG, GEONUC
C
C Local variables
C
      REAL(KIND=wp) :: CHARGE1, DELTA, E0, EI, EJ, XNUC1, YNUC1, ZNUC1
      CHARACTER(LEN=8) :: CNAME1
      INTEGER :: I, II, J, MGVN, NALM, NNUC, NNUC1, NOCSFI, NOVECR, 
     &           NREC, NSET, NSTAT
      CHARACTER(LEN=120) :: NAME
C
C*** End of declarations rewritten by SPAG
C
C**********************************************************************
C
C     READCIP reads CI data from unit NFTW in format used for
C     polyatomic targets
C     Note that MGVN is the value used by CONGEN, not Molecule
C
C**********************************************************************
C
c
c.... Find file
      CALL SEARCH(nftw,'CIDATA  ',ifail)
      IF(ifail.EQ.1)GO TO 1800
C
 5    READ(NFTW)Nset, nrec, NAME, nnuc, nocsfi, nstat, mgvn, s, sz, 
     &          nelt, e0
      IF(nset.EQ.nciset)THEN
         IF(nocsfi.NE.nocsf)GO TO 1920
         IF(mgvn.NE.irrep)GO TO 1800
         novecr=noveci
      ELSE
         novecr=nstat
      END IF
      DO i=1, nnuc
         READ(nftw)cname(i), geonuc(1,i), geonuc(2,i), geonuc(3,i), 
     &             charge(i)
      END DO
C
      IF(itime.EQ.0 .AND. nset.EQ.nciset)THEN
C....   Read integral header and check that the geometries match
         CALL SEARCH(nfti,'POLYAINP',ifail)
         READ(nfti)nnuc1
         IF(nnuc.NE.nnuc1)GO TO 1900
         DO ii=1, nnuc1
            READ(nfti)cname1, i, xnuc1, ynuc1, znuc1, charge1
            IF(cname1.NE.cname(ii) .OR. xnuc1.NE.geonuc(1,ii) .OR. 
     &         ynuc1.NE.geonuc(2,ii) .OR. znuc1.NE.geonuc(3,ii) .OR. 
     &         charge1.NE.charge(ii))GO TO 1900
         END DO
c
      ELSE IF(nset.EQ.nciset)THEN
C....   READ CI COEFFICIENTS
         CALL CIVIO(NFTW,1,NOCSFi,novecr,EIG,VEC,NALM,iphase)
         IF(NALM.NE.0)GO TO 2900
         IF(NOCSF.NE.NCTGT)THEN
            WRITE(NFT,2920)NOCSF, NCTGT
            STOP
         END IF
      ELSE
         READ(nftw)
         DO i=1, nstat
            READ(nftw)
         END DO
      END IF
      READ(nftw,END=4)
 4    IF(nset.NE.nciset)GO TO 5
c
c....   Add nuclear potential to electronic eigenvalues
      DO i=1, noveci
         eig(i)=eig(i)+e0
      END DO
c
! Hemal Varambhia DENPROP bug fix
! Here one checks for repeated roots for a Hamiltonian block of a certain
! symmetry. Denprop does appears to omit repeated roots when writing to the
! target properties file. This ought not to happen but as yet, no solution
! has been
      delta=1.0E-08_wp
      OUTER:DO i=1, noveci
         ei=eig(i)
         INNER:DO j=1, i
            ej=eig(j)
!             write(173,*)'ei=',ei,' ej=',ej,' i=',i,' j=',j
!             write(173,*)'ei-ej<delta=',ei-ej.lt.delta
!             write(173,*)'i.ne.j',i.ne.j
            IF(ei-ej.LT.delta .AND. i.NE.j)THEN
 
               WRITE(6,*)' WARNING: this molecule has repeated roots'
               WRITE(6,*)
     &                  ' DENPROP may be not be able work with repeated'
               WRITE(6,*)' hence proceed with care'
               GO TO 200
            END IF
         END DO INNER
      END DO OUTER
! end of DENPROP bug fix
 
 200  RETURN
c
 1800 WRITE(nft,1804)irrep, nftw
 1804 FORMAT(/' CIDATA for MGVN=',i1,' NOT FOUND ON UNIT',I3)
      STOP
c
 1900 WRITE(NFT,1910)
 1910 FORMAT(/' Target and integral data are inconsistent ')
      WRITE(nft,1911)cname1, cname(i), nnuc, nnuc1, xnuc1, geonuc(1,i), 
     &               ynuc1, geonuc(2,i), xnuc1, geonuc(3,i), charge1, 
     &               charge(i)
 1911 FORMAT(1x,2A8,2I5,4(/2F10.6))
      STOP
c
 1920 WRITE(NFT,1921)
 1921 FORMAT(/' Target and formula data are inconsistent ')
      WRITE(nft,1922)nocsfi, nocsf
 1922 FORMAT(' NOCSF(target) =',i6,'    NOCSF(formula) =',i6)
      STOP
c
 2900 WRITE(NFT,2910)
 2910 FORMAT(/' UNABLE TO GET CI-TARGET VECTOR ')
      STOP
c
 2920 FORMAT(/' HAMILTONIAN TRANSFORMATION DATA INCONSISTENT ',2I10)
c
      END SUBROUTINE READCIP
!*==readcid.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE READCID(nftw,nciset,NAME,NHE,NHD,DTNUC,nocsf,NSTAT,
     &                   nctgt,EIG,VEC,iphase,NFT)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER(LEN=120) :: NAME
      INTEGER :: NCISET, NCTGT, NFT, NFTW, NOCSF, NSTAT
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp), DIMENSION(*) :: EIG, VEC
      INTEGER, DIMENSION(*) :: IPHASE
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(20) :: NHE
      INTENT (IN) NCTGT
      INTENT (INOUT) NHD, NOCSF
C
C Local variables
C
      INTEGER :: NALM, NT
C
C*** End of declarations rewritten by SPAG
C
C**********************************************************************
C
C     READCID reads CI data from unit NFTW in format used for diatomic
C     targets
C
C**********************************************************************
C
c
      CALL MOVEW(NFTW,nciset,NALM,0,NFT)
      IF(NALM.NE.0)GO TO 2900
C
      READ(NFTW,END=2900)NT, NHD, NAME, NHE, DTNUC
      NOCSF=NHD(2)
      IF(nstat.GT.nhd(3))THEN
         WRITE(nft,2980)nstat, nhd(3)
 2980    FORMAT(/' Target state number',i4,' requested',
     &          /' Only               ',i4,' on file: STOP')
         STOP
      END IF
c
C     READ CI COEFFICIENTS
C
      CALL CIVIO(NFTW,1,NOCSF,NSTAT,EIG,VEC,NALM,iphase)
      IF(NALM.NE.0)GO TO 2900
      IF(NOCSF.EQ.NCTGT)THEN
         CALL PRTHD(nciset,NHD,NAME,NHE,DTNUC,NSTAT,EIG,NFT)
      ELSE
         WRITE(NFT,2920)NOCSF, NCTGT
 2920    FORMAT(/' HAMILTONIAN TRANSFORMATION DATA INCONSISTENT ',2I10)
         STOP
      END IF
c
      RETURN
c
 2900 WRITE(NFT,2910)
 2910 FORMAT(/' UNABLE TO GET CI-TARGET VECTOR ')
      STOP
      END SUBROUTINE READCID
!*==search.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE SEARCH(IUNIT,A,ifail)
      USE GLOBAL_UTILS, ONLY : CWBOPN
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: IWRITE=6
      CHARACTER(LEN=8), PARAMETER :: A1='********'
C
C Dummy arguments
C
      CHARACTER(LEN=8) :: A
      INTEGER :: IFAIL, IUNIT
      INTENT (IN) A
      INTENT (OUT) IFAIL
C
C Local variables
C
      CHARACTER(LEN=32) :: B
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     Utility to search a dataset IUNIT for a header A where the dataset
C     is assumed to have MOLECULE-SWEDEN convention headers. The header
C     convention is
C
C     '********', '        ', '        ', 'ABCDEFG'
C
C     with ABCDEFG being a character string such as ONELINT etc.
C
C***********************************************************************
C
      CALL cwbopn(iunit)
      ifail=0
c      IREC=0
C
 1    READ(IUNIT,END=990)B
C
c      IREC=IREC+1
c      WRITE(IWRITE,1500) IREC,B
C
      IF(B(1:8).EQ.A1 .AND. B(25:32).EQ.A)RETURN
c
      GO TO 1
C
C---- Process error condition namely, header not found by end of file.
C
 990  WRITE(IWRITE,9900)A, IUNIT
      ifail=1
      STOP 999
C
C---- Format Statements
C
 1500 FORMAT(10X,'Record No. = ',I6,' Data = ',A)
 9900 FORMAT(/,10X,'**** Error in SEARCH: ',//,10X,
     &       'Attempt to find header (A) = ',A,' on unit = ',I3,/,10X,
     &       'has failed.',/)
C
      END SUBROUTINE SEARCH
!*==civio.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE CIVIO(NFT,NRW,NK,NS,EI,CV,NALM,iphz)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NALM, NFT, NK, NRW, NS
      REAL(KIND=wp), DIMENSION(*) :: CV, EI
      INTEGER, DIMENSION(*) :: IPHZ
C
C Local variables
C
      INTEGER :: I, IB, M
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     CIVIO CONTROLS THE I/O OF CI COEFFICIENTS AND STATE DATA
C
C***********************************************************************
C
C
C     READ ENERGY AND CSF SPECIFICATION FOR EACH CI STATE
C
      CALL CIVIOA(NFT,NRW,NK,NS,EI,NALM,iphz)
      IF(NALM.NE.0)GO TO 90
C
C     READ COEFFICIENTS FOR EACH CI STATE
C
      IB=1
      DO I=1, NS
         M=I
         CALL CIVIOB(NFT,NRW,M,NK,CV(IB),NALM)
         IF(NALM.NE.0)EXIT
         IB=IB+NK
      END DO
C
 90   RETURN
      END SUBROUTINE CIVIO
!*==civioa.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE CIVIOA(NFT,NRW,NOCSF,NSTAT,EI,NALM,iphz)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NALM, NFT, NOCSF, NRW, NSTAT
      REAL(KIND=wp), DIMENSION(NSTAT) :: EI
      INTEGER, DIMENSION(nocsf) :: IPHZ
      INTENT (IN) NFT, NOCSF, NRW, NSTAT
      INTENT (OUT) NALM
      INTENT (INOUT) EI, IPHZ
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     CIVIOA READS OR WRITES CI DUMPFILES
C     IKEEP is no longer used.  A dummy variable is written for format
C     consistency
C
C***********************************************************************
C
C
      NALM=0
      IF(NRW.EQ.0)THEN
         WRITE(NFT,ERR=200)iphz, EI
      ELSE
         READ(NFT,ERR=200)iphz, EI
      END IF
      RETURN
C
 200  NALM=1
      RETURN
      END SUBROUTINE CIVIOA
!*==civiob.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE CIVIOB(NFT,NRW,NTH,NOCSF,VC,NALM)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NALM, NFT, NOCSF, NRW, NTH
      REAL(KIND=wp), DIMENSION(NOCSF) :: VC
      INTENT (IN) NFT, NOCSF, NRW
      INTENT (OUT) NALM
      INTENT (INOUT) NTH, VC
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     CIVIOB READS OR WRITES CI COEFFICIENTS IN THE DUMPFILE FORMAT
C
C***********************************************************************
C
C
      IF(NRW.EQ.0)THEN
         WRITE(NFT,ERR=200)NTH, VC
      ELSE
         READ(NFT,ERR=200)NTH, VC
      END IF
      NALM=0
      RETURN
C
 200  NALM=1
      RETURN
      END SUBROUTINE CIVIOB
!*==plotdpoten.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE plotdpoten(iplotfile,NMSET,NNUC,NTARG,ISW,NUCCEN,RGEOM,
     &                      CHARG,ISPIN,MGVN,GUTOT,IREFL,ENERGY,ISTATE,
     &                      JSTATE,IWRITE,ISYMTYP,NSYM,grounden,ksym)
C***********************************************************************
C     PLOTPOTEN: PLOTPOTentialENergies.
C     This routine includes part of routine WRDTARG in order to write
C     potential energy in a data format suitable for plots of energy versus
C     internuclear distance. If calculations are performed for several
C     internuclear distances, one can swich on the plot data writer setting
C     flag iplotfg=1 in input data and then cat units ipotv (fort.45) from
C     each internuclear distance denprop run in one data file.
C     This means that the target energies are written in ipotv in the order
C     given by the CI data and not in ascending order as in the target
C     property file. The energies are given in hartree and are relative to
C     the ground target energy at equilibrium geometry, if grounden is
C     provided in input when iplotfg=1. Also a write statment depending
C     on flag iqdfg=1 in namelist has been included in subroutine WRDTARG
C     in order to write the ordered target energies on a file. These
C     data may be used in the outer region for resonance parameters from
C     complex quantum defect calculations (using subroutine respars in
C     mcqd, NV 03).
C
C     WRDTARG - WRite a TARGet properties file for linear molecule
C     WRPTARG - WRite a TARGet properties file for non-linear molecule
C
C     Linkage:
C
C        GETSET,CMAKSTN
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: KEYH=9, NULL=0
      CHARACTER(LEN=11), PARAMETER :: TFORM='FORMATTED  '
C
C Dummy arguments
C
      INTEGER :: IPLOTFILE, ISW, ISYMTYP, IWRITE, NMSET, NNUC, NSYM, 
     &           NTARG, NUCCEN, KSYM
      REAL(KIND=wp) :: GROUNDEN
      REAL(KIND=wp), DIMENSION(*) :: CHARG, ENERGY
      INTEGER, DIMENSION(*) :: GUTOT, IREFL, ISPIN, ISTATE, JSTATE, MGVN
      REAL(KIND=wp), DIMENSION(3,*) :: RGEOM
      INTENT (IN) CHARG, ENERGY, ISW, IWRITE, NNUC, NTARG, NUCCEN, 
     &            RGEOM, GROUNDEN
C
C Local variables
C
      CHARACTER(LEN=36), EXTERNAL :: CMAKSTN
      REAL(KIND=wp) :: RAB, TOL=1.E-10_wp
      INTEGER :: I, ICHARG, INUC, IPOTV, MNUC
      INTEGER :: INT, LEN
      INTEGER, DIMENSION(3) :: MASS
      REAL(KIND=wp), DIMENSION(3) :: RMOI
      CHARACTER(LEN=36) :: RNAME
      LOGICAL :: ZISWIS0
C
C*** End of declarations rewritten by SPAG
C
!      PARAMETER(XZERO=0.0D+00,ipotv=45)
C
C---- Error check the input as much as possible
C
C..... Must have a +ve number of nuclei and target states too.
C
      IF(NNUC.LT.1 .OR. NTARG.LT.1)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9910)NNUC, NTARG
         STOP
      END IF
C
C..... ISW must be zero or unity.
C
      IF(ISW.NE.0 .AND. ISW.NE.1)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9920)ISW
         STOP
      END IF
C
C..... NUCCEN must be in the range 1 <= NUCCEN <= NNUC
C
      IF(NUCCEN.LT.1 .OR. NUCCEN.GT.NNUC)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9930)NUCCEN, NNUC
         STOP
      END IF
C
C---- Initialize the output character string to blanks
C
      DO I=1, LEN(RNAME)
         RNAME(I:I)=' '
      END DO
C
C---- First position iplotfile at the begining of set NTSET
C
      CALL GETSET(iplotfile,NMSET,KEYH,TFORM)
C
      ZISWIS0=ISW.EQ.0
C
      IF(ZISWIS0)THEN
c
C     Get nuclear data (this has been done before, but its easier to
C     repeat the calculation than to store the output)
c
         DO inuc=1, nnuc
            IF(inuc.EQ.nuccen)CYCLE
c
            IF(rgeom(1,inuc).NE.xzero .OR. rgeom(2,inuc).NE.xzero)THEN
               WRITE(iwrite,9940)
               STOP
            END IF
C
         END DO
      END IF
C
      RAB=XZERO
      mnuc=0
      DO I=1, NNUC
         RAB=RAB+ABS(RGEOM(3,I))
         ICHARG=INT(ABS(CHARG(I)))
         IF(I.EQ.NUCCEN .AND. icharg.EQ.0)CYCLE
         mnuc=mnuc+1
      END DO
C
      DO I=1, NTARG
         RNAME=CMAKSTN(I,MGVN(I),ISPIN(I),GUTOT(I),IREFL(I),ISYMTYP,
     &         NSYM,ksym)
         WRITE(iplotfile,2100)energy(i), rname
      END DO
      WRITE(ipotv,2220)rab, ((energy(i)-grounden),i=1,ntarg)
C
      RETURN
C
C---- Format Statement
C
 2100 FORMAT(D20.12,2x,A35)
 2220 FORMAT(f8.4,7(f8.4))
C
 9900 FORMAT(/10X,'**** Error in WRTARG: ',/)
 9910 FORMAT(10X,'Either, or both NNUC and NTARG are less than zero:',/,
     &       10X,'NNUC = ',I5,' NTARG= ',I5,/)
 9920 FORMAT(10X,'ISW should be 0 or 1. It is input as = ',I10,/)
 9930 FORMAT(10X,'NUCCEN must lie in the range 0 < NUCCEN <= ',I3,//,
     &       10X,'NUCCEN has been input as = ',I5,/)
 9940 FORMAT(10X,'Routine has been hardwired for linear case')
C
      END SUBROUTINE PLOTDPOTEN
!*==plotppoten.spg  processed by SPAG 6.56Rc at 16:19 on 19 Nov 2010
      SUBROUTINE plotppoten(iplotfile,NMSET,NNUC,NTARG,NUCCEN,RGEOM,
     &                      CHARG,ISPIN,MGVN,ENERGY,ISTATE,JSTATE,
     &                      IWRITE,ISYMTYP,NSYM,grounden,ksym)
C***********************************************************************
C     PLOTPPOTEN: PLOT Polyatomics POTential ENergies.
C     This routine includes part of routine WRPTARG in order to write
C     potential energy in a data format suitable for plots of energy versus
C     geometry. A numbering of geometry configuration is set in namelist (geocfn).
C     If calculations are performed for several geometry, one can swich on the plot
C     data writer setting flag iplotfg=1 and geocfn= n in the input data. Here n is
C     the number assigned to current fixed geometry calculation, or, in other words,
C     to the  current coordinate configuration of the molecular target. It can be
C     replaced with one of the internal coordinates when one performs calculations
C     keeping  fixed the other two coordinates. in input data and then cat units
C     ipotv (fort.45) from each internuclear distance denprop run in one data file.
C     This means that the target enegies are written in ipotv in the order
C     given by the CI data and not in ascending order as in the target
C     property file. The energies are given in hartree and are relative to
C     the ground target energy at equilibrium geometry, if grounden is
C     provided in input when iplotfg=1. Also a write statment depending
C     on flag iqdfg=1 in namelist has been included in subroutine WRDTARG
C     in order to write the ordered target energies on a file. These
C     data may be used in the outer region for resonance parameters from
C     complex quantum defect calculations (see respars, NV).
C
C     Linkage:
C
C        GETSET,CMAKSTN
C
C***********************************************************************
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: KEYH=6, NULL=0, IPOTV=45
      CHARACTER(LEN=11), PARAMETER :: TFORM='FORMATTED  '
C
C Dummy arguments
C
      REAL(KIND=wp) :: GROUNDEN
      INTEGER :: IPLOTFILE, ISYMTYP, IWRITE, KSYM, NMSET, NNUC, NSYM, 
     &           NTARG, NUCCEN
      REAL(KIND=wp), DIMENSION(*) :: CHARG, ENERGY
      INTEGER, DIMENSION(*) :: ISPIN, ISTATE, JSTATE, MGVN
      REAL(KIND=wp), DIMENSION(3,*) :: RGEOM
      INTENT (IN) ENERGY, GROUNDEN, NNUC, NTARG, NUCCEN
C
C Local variables
C
      CHARACTER(LEN=36), EXTERNAL :: CMAKSTN
      CHARACTER(LEN=8) :: CSYMBOL
      INTEGER :: GEOCFN, I, NNUCT
      INTEGER :: LEN
      REAL(KIND=wp), DIMENSION(3) :: RMOI
      CHARACTER(LEN=36) :: RNAME
      LOGICAL :: SCATI
C
C*** End of declarations rewritten by SPAG
C
      scati=.FALSE.
      DO I=1, NNUC
         IF(i.EQ.nuccen)scati=.TRUE.
      END DO
C
C---- Initialize the output character string to blanks
C
      DO I=1, LEN(RNAME)
         RNAME(I:I)=' '
      END DO
C
C---- First position nplotfile at the begining of set NTSET
C
      CALL GETSET(iplotfile,NMSET,KEYH,TFORM)
C
C --- If the input is from a scattering run, then adjust nnuc
      IF(scati)THEN
         nnuct=nnuc-1
      ELSE
         nnuct=nnuc
      END IF
c
      WRITE(iplotfile,2221)
      DO i=1, ntarg
         rname=cmakstn(i,mgvn(i),ispin(i),null,null,isymtyp,nsym,ksym)
         WRITE(iplotfile,2223)energy(i), rname
      END DO
      IF(ntarg.GT.1)THEN
         WRITE(iplotfile,2224)
      END IF
 2221 FORMAT(/,'STATE ENERGIES FOR THE MOLECULAR TARGET ARE: ',/,/,3x,
     &       'Energy (Hartree) ',20x,'State',/)
 2223 FORMAT((f18.12),10x,A36)
 2224 FORMAT(/,'See file results/SelectedTargetStates.txt',
     &       ' for vertical excitation energies',/,
     &       ' of the selected states in energy order.',/,/)
!       write (ipotv,2220) (rgeom(1,j),rgeom(2,j),rgeom(3,j),j=1,nnuct),
!     * ((energy(i)-grounden),i=1,ntarg)
!I THINK THAT THIS ABOVE SHOULD BE ANYWAY SOMETHING LIKE THIS:
      WRITE(ipotv,2220)geocfn, ((energy(i)-grounden),i=1,ntarg)
!
C
      RETURN
C
C---- Format Statement
C
 2100 FORMAT(D20.12,2X,A36)
 2220 FORMAT(i5,10(f14.4))
 
C
      END SUBROUTINE PLOTPPOTEN
