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
!> This module contains some routines taken from the UKRmol-in code GAUSPROP. These routines read the basis set from the integral files.
module ukrmol_routines

contains
!
!---- Read and print the header information on the file of atomic
!     basis function integrals. It is also copied to the header of
!     the file of property integrals.
!
!     It is useful to remember that:
!
!        LTRAN is the maximum number of transformation functions from
!              contracted Gaussians to symmetry adapted basis functions
!              Thus an accurate value would be
!
!                Sum over all sym. adapted functions of the
!                             no. of contracted functions
!
!           IC pointer to storage for the first column of the symmetry
!              adapted basis function table.
!        ITRAN pointer to storage for the second column of the symmetry
!              adapted basis function table
!       ICTRAN pointer to storage for the third column of the symmetry
!              adapted basis function array
!
!        NSYMT is the number of Irreducible representations in this
!              molecular point group i.e. symmetries
!         NBFT is the number of symmetry adapted basis functions
!              per molecular symmetry.
!
!     For the purposes of the PRTINF routine these pointers are used
!     solely as pointers to workspaces. It is in the routine GTGETAB
!     that the symmetry adapted basis function table is actuall built.
!     iclab will hold the labels of the symmetry adapted functions
!     and must be saved.  LTRAN is arbitrary and could be changed.
!

      SUBROUTINE PRTINF(ITAPE,IWRITE,NSYMT,NBFT,CLMQTYPE,jsum)
!***********************************************************************
!
!     Prints the header from the dataset of Gaussian integrals that
!     are computed by the code MOLECULE - this means that all data
!     prior to the records of integrals is printed. Note that some data
!     is returned to the caller here so that this is more than just a
!     utility routine.
!
!     Input data:
!          ITAPE  Logical unit for the MOLECULE output
!         IWRITE  Logical unit for the printer
!             IC  Workspace array for no. of primitives per contracted
!                 Gaussian function
!          ITRAN  Workspace array for integers
!          CTRAN  Workspace array for real*8 variables
!       CLMQTYPE  Angular behaviour array which lists all l,m,q
!                 functions.
!
!     Output data:
!           NSYMT Number of symmetries (IRRs) in this molecular point
!                 group.
!            NBFT Number of basis functions per symmetry
!
!***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE GLOBAL_UTILS, ONLY : CWBOPN
      USE GAUSPROP_DATA, ONLY : MAXDIM
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: ITAPE, IWRITE, JSUM, NSYMT
      CHARACTER(LEN=4), DIMENSION(*) :: CLMQTYPE
      INTEGER, DIMENSION(*) :: NBFT
      INTENT (IN) CLMQTYPE, IWRITE
      INTENT (INOUT) JSUM, NBFT, NSYMT
!
! Local variables
!
      REAL(KIND=wp) :: CH, CX, POTNUC, XX
      CHARACTER(LEN=192) :: CHEADER
      CHARACTER(LEN=8) :: CL, CNAKO, CNAME
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: CNUCNAM
      REAL(KIND=wp), DIMENSION(maxdim) :: CTRAN
      INTEGER :: I, IABAS, IBBAS, IDOSPH, II, J, JEND, K, KB, L, NMAX, NNUC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IC
      INTEGER, DIMENSION(maxdim) :: ITRAN
      CHARACTER(LEN=8), DIMENSION(4) :: LABEL
      REAL(KIND=wp), DIMENSION(3) :: X
!
!*** End of declarations rewritten by SPAG
!
! Initialization with blanks, replacing a DATA statement
      CL='        '
      CNAKO='        '
!
!---- Banner header is placed on the printer file
!
      WRITE(IWRITE,1000)
!
!---- Rewind the file and read the identifier record at the start of
!     it.
!
!
      CALL cwbopn(ITAPE)
!
!---- Data on this record is as follows:
!
!      CHEADER - Character header from input to MOLECULE code
!        NSYMT - No. of Irreducible Representations (Symmetries) in the
!                set of basis functions
!         NBFT - No. of atomic basis functions per symmetry
!       POTNUC - Nuclear potential energy
!       IDOSPH - Flag stating whether or not we have spherical harmonic
!                basis
!
      READ(ITAPE)CHEADER, NSYMT, (NBFT(I),I=1,NSYMT), POTNUC, IDOSPH
!
      WRITE(IWRITE,500)CHEADER(1:65)
      WRITE(IWRITE,510)NSYMT, (I,NBFT(I),I=1,NSYMT)
      WRITE(IWRITE,520)POTNUC, IDOSPH
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
      WRITE(IWRITE,2000)(LABEL(I),I=1,4)
      WRITE(IWRITE,2010)NNUC
!
!...... Loop over nuclei - obtain positions and charges
!
      WRITE(IWRITE,2015)
!
      DO I=1, LEN(CNAME)
         CNAME(I:I)=' '
      END DO
!
      DO I=1, NNUC
         READ(ITAPE)CNAME(1:4), II, (X(J),J=1,3), CH
         WRITE(IWRITE,2020)CNAME, II, (X(J),J=1,3), CH
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
      WRITE(IWRITE,2500)IBBAS, IABAS
!
!---- Array IC has one entry for each contracted function. It gives the
!     number of primitives in that contraction.
!
      ALLOCATE(ic(iabas))
      READ(ITAPE)(IC(I),I=1,IABAS)
!      write(6,*)  'ic(iabas) is',ic(iabas),iabas
!
!     WRITE(IWRITE,2550)
!     WRITE(IWRITE,2560) (I,IC(I),I=1,IABAS)
!
!---- For each contracted function we loop over the number of primitives
!     within it and obtain data on each.
!
!     WRITE(IWRITE,2600)
!
      DO I=1, IABAS
         JEND=IC(I)
!     WRITE(IWRITE,2610) I,JEND
         DO J=1, JEND
            READ(ITAPE)CL(1:4), KB, CNAKO(1:4), XX, CX
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
!      WRITE(IWRITE,3000) (LABEL(I),I=1,4)
!      WRITE(IWRITE,3010) IABAS
!
      jsum=0
      DO I=1, IABAS
         READ(ITAPE)J, (ITRAN(K),CTRAN(K),K=1,J)
!         WRITE(6,3020) I,J,(ITRAN(K),CTRAN(K),K=1,J)
         jsum=jsum+j
      END DO
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
!
      ALLOCATE(cnucnam(nmax))
!
!      WRITE(IWRITE,4000) (LABEL(I),I=1,4)
!      WRITE(IWRITE,4010) NMAX,IABAS
!
      READ(ITAPE)NMAX, (CNUCNAM(I),CTRAN(I),I=1,NMAX)
!
!      WRITE(IWRITE,4020) NMAX
!      WRITE(IWRITE,4030) (I,CNUCNAM(I),CTRAN(I),I=1,NMAX)
!
!      WRITE(IWRITE,4040) IABAS
!
      DO I=1, IABAS
         READ(ITAPE)J, K, L
         WRITE(IWRITE,4050) I,J,K,L,CLMQTYPE(L)
      END DO
!
      DEALLOCATE(ic,cnucnam)
!
      RETURN
!
!---- Format Statements
!
 500  FORMAT(5X,'Header Card (1:65) : ',A,/)
 510  FORMAT(5X,'No. of symmetries in basis set = ',I3,//,5X,'No. of basis functions per symmetry: ',/,(5X,I2,1X,I3))
 520  FORMAT(/,5X,'Nuclear Potential Energy = ',F15.7,' (Hartrees) ',//,5X,'Cartesian/Spherical flag = ',I8,/)
!
 1000 FORMAT(///,10X,'Atomic Integrals File: Header Records ',/)
!
 2000 FORMAT(5X,'Section 1 Header = ',4A,/)
 2010 FORMAT(5X,'Number of nuclear centers = ',I3,/)
 2015 FORMAT(5X,' Symbol ',1X,'No.',' (X,Y,Z) Co-ordinates ',T50,' Charge ',/)
 2020 FORMAT(5X,A,1X,I2,1X,3(F10.6,1X),F6.3)
 2500 FORMAT(//,5X,'Number of primitive  Gaussian functions = ',I4,/,5X,'Number "" contracted  " " "    " " " "  = ',I4,//)
 2550 FORMAT(5X,'Primitive functions per contracted function:',/)
 2560 FORMAT((5X,8(I3,'.',1X,I3,1X)))
 2600 FORMAT(/,5X,'A N A L Y S I S  of  each  C O N T R A C T I O N',/)
 2610 FORMAT(/,5X,'Contracted function no. = ',I3,/,5X,'No. of primitives       = ',I3,/)
 2620 FORMAT(1X,'Primitive = ',I2,1X,'CL = ',A,1X,'KB = ',I3,1X,'CNAKO = ',A,1X,'XX = ',D12.5,1X,'CX = ',D12.5,1X)
!
 3000 FORMAT(/,5X,'Section 2 Header = ',4A,/)
 3010 FORMAT(5X,'Number of symmetry adapted functions = ',I3,/)
 3020 FORMAT(/,5X,'Symmetry function no. = ',I3,/,5X,'No. of components     = ',I3,/,5X,' (ITRAN,CTRAN) pairs: ',//,(5X,4(I4,1X,D12.5,1X)))
!
 4000 FORMAT(/,5X,'Section 3 Header = ',4A,/)
 4010 FORMAT(5X,'Number of atomic nuclei (NMAX) = ',I3,//,5X,'Number ""   ""   basis functions (IABAS) = ',I5,//)
 4020 FORMAT(5X,'NMAX = ',I4,' Nuclei/Charge  pair definitions follow: ',/)
 4030 FORMAT((5X,3(I3,'.',1X,A,1X,D12.5,1X)))
 4040 FORMAT(/,5X,'No. of items (IABAS) = ',I5,/)
 4050 FORMAT(5X,'Item = ',I3,' J = ',I3,' K = ',I3,' L = ',I3,1X,'Angular function = ',A)
!
      END SUBROUTINE PRTINF

      SUBROUTINE GTGETAB(IWRITE,ITAPE,NCONTRAS,NPRCONT,ICONTRNO,CSYMBOL,ISYMBNUM,&
     &                   ngnuc,CANGULAR,cscatter,nnuc,cnucname,XPONENT,&
     &                   XCONTCOF,XCONTCOFNN,NSYMABF,ISYMABF,ITRAN,CTRAN,LTRAN,&
     &                   LSABFTAB,ZSCATCFN,ZSCATSFN,MAXSFN)
!***********************************************************************
!
!     GTGETAB - GeT primitive and symmetry adapted Gaussian TABles
!
!     This routine reads data from the atomic integrals file into a
!     tabular format that is used later by the code. The first table
!     relates primitive Gaussian functions to contracted ones while the
!     second relates contracted Gaussian functions to symmetry adapted
!     basis functions.
!     ZM: I've adjusted this routine to return also the unnormalized contraction
!     coefficients XCONTCOFNN.
!
!     Input data:
!         IWRITE Logical unit for the printer
!          ITAPE Logical unit holding the data on the functions
!                i.e. Output by the MOLECULE code
!       NCONTRAS Number of contracted Gaussian functions
!        NPRCONT Table of the no. of primitives per contracted function
!          LTRAN Size of the arrays ISYMABF,ITRAN and CTRAN
!
!     Output data:
!        ICONTRNO Sequence number of the contracted Gaussian function
!                 to which this primitive belongs
!         CSYMBOL Symbol which identifies the symmetry independent
!                 atom to which this basis function belongs
!        ISYMBNUM For each atom MOLECULE generates CSYMBOL and the number ISYMBNUM which 
!                 together can be used to identify all atoms in the molecule. The ISYMBNUM is an
!                 array which contains this number for each primitive GTO in the basis and can be
!                 used (together with CSYMBOL) to identify the nucleus on which the GTO is centered.
!        CANGULAR Designation of the angular term for this primitive
!                 Gaussian
!         XPONENT The exponent for this primitive Gaussian function
!        XCONTCOF The contraction coefficient for this primitive Gaussian (normalized)
!      XCONTCOFNN The contraction coefficient for this primitive Gaussian (unnormalized)
!         NSYMABF Number of symmetry adapted basis functions
!         ISYMABF Pointer array defining symmetry adapted basis functio
!           ITRAN Pointer to contracted Gaussian function which belongs
!                 to the symmetry adapted basis function in the
!                 corresponfing entry in ISYMABF
!           CTRAN Coefficient partnering the corresponding entries in
!                 ISYMABF and ITRAN.
!        LSABFTAB Length of the symmetry adapted basis function table
!        ZSCATCFN List of values set .true. for scattering contracted
!                 functions and .false. for target contractions.
!        ZSCATSFN List of values set .true. for scattering symmetry
!                 adapted basis functions and .false. for target ones.
!          MAXSFN Maximum space in the array ZSCATSFN
!
!     Notes:
!
!        The contracted Gaussian expansions are normalized after
!     being read in. This is because the output from MOLECULE may be
!     un-normalized.
!
!     Linkage:
!
!        GTINDEX, GTNORM
!
!***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE gaus_shared, ONLY : GTINDEX, GTNORM
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      CHARACTER(LEN=4) :: CSCATTER
      INTEGER :: ITAPE, IWRITE, LSABFTAB, LTRAN, MAXSFN, NCONTRAS, NNUC, NSYMABF
      CHARACTER(LEN=4), DIMENSION(*) :: CANGULAR, CSYMBOL
      CHARACTER(LEN=4), DIMENSION(nnuc) :: CNUCNAME
      REAL(KIND=wp), DIMENSION(*) :: CTRAN, XCONTCOF, XCONTCOFNN, XPONENT
      INTEGER, DIMENSION(*) :: ICONTRNO, ISYMABF, ITRAN, NGNUC, ISYMBNUM
      INTEGER, DIMENSION(ncontras) :: NPRCONT
      LOGICAL, DIMENSION(ncontras) :: ZSCATCFN
      LOGICAL, DIMENSION(maxsfn) :: ZSCATSFN
      INTENT (IN) CNUCNAME, CSCATTER, ITAPE, IWRITE, LTRAN, MAXSFN, NCONTRAS, NNUC
      INTENT (OUT) CTRAN, ICONTRNO, ISYMABF, LSABFTAB, NGNUC, ZSCATSFN, ISYMBNUM
      INTENT (INOUT) CSYMBOL, ITRAN, NSYMABF, ZSCATCFN
!
! Local variables
!
      INTEGER :: I, ILAST, INUC, IXPOWER, IYPOWER, IZPOWER, J, JCODE, JNC, JNLAST, JNUC, K, KBIAS, L, MULT
!
!*** End of declarations rewritten by SPAG
!
!---- Print banner header
!
!        WRITE(IWRITE,1000)
!        WRITE(IWRITE,1010) ITAPE,NCONTRAS
!
!---- We assume that ITAPE is positioned so that we may begin reading.
!
      K=0
      jnlast=0
      ilast=0
!
!      WRITE(IWRITE,'(/,10X,"Contracted and primitive GTO table")')
!      WRITE(IWRITE,'(  10X,"==================================",/)')
!
      DO I=1, NCONTRAS
         DO J=1, NPRCONT(I)
            K=K+1
            ICONTRNO(K)=I
            READ(ITAPE,ERR=900)CSYMBOL(K), ISYMBNUM(K), CANGULAR(K), XPONENT(K), XCONTCOF(K)
!
!            WRITE(IWRITE,'(10X,i4,1X,i4,1X,a,1X,i4,1X,a,1X,2e25.15)') K, I, CSYMBOL(K), ISYMBNUM(K), CANGULAR(K), XPONENT(K), XCONTCOF(K)
            XCONTCOFNN(K)=XCONTCOF(K)
!
            IF(j.EQ.1)THEN
               IF(CSYMBOL(K).EQ.CSCATTER)THEN
                  ZSCATCFN(I)=.TRUE.
               ELSE
                  ZSCATCFN(I)=.FALSE.
               END IF
               IF(i.GT.1)THEN
                  IF(csymbol(k).NE.csymbol(k-1))THEN
                     ilast=i-1
                     jnlast=jnuc
                  END IF
               END IF
               DO inuc=1, nnuc
                  IF(csymbol(k).EQ.cnucname(inuc))jnuc=inuc
               END DO
!
!---- Must check for 'equivalent atoms' with non-unique labels.
               IF(jnuc.NE.jnlast .AND. jnuc.NE.jnlast+1)THEN
                  mult=jnuc-jnlast
                  jnc=mult-1-mod(i-ilast-1,mult)
                  jnuc=jnuc-jnc
               END IF
            END IF
            ngnuc(k)=jnuc
         END DO
!
!.... Build the JCODE value by looking at the angular terms and
!     then invoke the normalization routine for this contraction.
!
!        At this stage K is a high water marker and will point to the
!        last primitive read. Since all primitives within a contraction
!        are of the same angular behaviour then we may use this to work
!        out the JCODE value. However a new pointer must be set to
!        give access to the exponents and coefficients.
!
         CALL GTINDEX(CANGULAR(K),IXPOWER,IYPOWER,IZPOWER)
         JCODE=IXPOWER+IYPOWER+IZPOWER+1
!
         L=K-NPRCONT(I)+1
         CALL GTNORM(XPONENT(L),XCONTCOF(L),NPRCONT(I),JCODE)
      END DO
!
!---- Having read the contracted Gaussian information we must now
!     skip the header that follows it and read the symmetry adapted
!     basis function information. It defines linear combinations of
!     contracted Gaussians which make symmetry adapted functions
!
      READ(ITAPE,ERR=910)
      READ(ITAPE,ERR=920)NSYMABF
!
      IF(NSYMABF.GT.MAXSFN)THEN
         WRITE(IWRITE,9900)
         WRITE(IWRITE,9970)MAXSFN, NSYMABF
         STOP 992
      END IF
!
      KBIAS=0
!
      DO I=1, NSYMABF
         READ(ITAPE,ERR=930)J, (ITRAN(KBIAS+K),CTRAN(KBIAS+K),K=1,J)
         DO K=1, J
            ISYMABF(KBIAS+K)=I
         END DO
!
!.... We need only test one of the contracted functions to see
!     if it is a scattering function or not.
!
         ZSCATSFN(I)=.FALSE.
         IF(ZSCATCFN(ITRAN(KBIAS+1)))ZSCATSFN(I)=.TRUE.
!
         KBIAS=KBIAS+J
!
         IF(KBIAS.GT.LTRAN)THEN
            WRITE(IWRITE,9900)
            WRITE(IWRITE,9950)KBIAS, LTRAN
            STOP 999
         END IF
      END DO
!
!---- Print out the symmetry adapted basis function table
!
      WRITE(IWRITE,5000)
!
      DO K=1,KBIAS
         WRITE(IWRITE,5010) K,ISYMABF(K),ITRAN(K),CTRAN(K),ZSCATSFN(ISYMABF(K))
      ENDDO
!
!---- Store the length of the table for later use in the transform
!     code.
!
      LSABFTAB=KBIAS
!
      RETURN
!
!---- Error condition handlers
!
!..... Error during of a record of primitive basis function data
!
 900  CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9910)K
      STOP 999
!
!..... Error reading header of SYMTRANS data
!
 910  CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9920)
      STOP 999
!
!..... Error reading number of symmetry adapted functions
!
 920  CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9930)
      STOP 999
!
!..... Error reading a record for one symmetry adapted function
!
 930  CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9940)I
      STOP 999
!
!---- Format Statements
!
 1000 FORMAT(//,15X,'====> GTGETAB - READ PRIMITIVE FN. TABLE <====',/)
 1010 FORMAT(/,15X,'Logical unit for atomic ints     (ITAPE) = ',I5,/,15X,'Total number of contracted fns. (NCONTRAS) = ',I5)
!
 5000 FORMAT(/,10X,'Symmetry Adapted Function Definition Table',/,10X,&
     &       '==========================================',//,10X,&
     &       'This relates contracted Gaussian functions',/,10X,&
     &       'to symmetry adapted basis functions for the',/,10X,&
     &       'molecular point group. ',///,10X,&
     &       'Row  Basis Function  Cont. Gaussian    Coefficient',&
     &       '    Scatter ',/,10X,&
     &       '---  --------------  --------------    -----------',&
     &       '    ------- ',/)
 5010 FORMAT(5X,I7,3X,I10,8X,I10,7X,F10.6,1X,L7)
 9900 FORMAT(//,10X,'**** Error in GTGETAB: ',//)
 9910 FORMAT(10X,'Reading of record number = ',I5,' has failed ',/)
 9930 FORMAT(10X,'Reading NSYMABF has failed ',/)
 9920 FORMAT(10X,'Reading header for SYMTRANS has failed',/)
 9940 FORMAT(10X,'Reading of SYMTRANS record no. = ',I5,' failed',/)
 9950 FORMAT(10X,'Not enough space in the symmetry adapted function',' table.',//,10X,'Need at present = ',I10,/,10X,'Given value     = ',I10,/)
 9970 FORMAT(10X,'Not enough space in ZSCATSFN array ',//,10X,'MAXSFN = ',I10,' NSYMABF = ',I10,/)
!
      END SUBROUTINE GTGETAB
!*==gpsplit.spg  processed by SPAG 6.56Rc at 14:46 on 23 Feb 2010
      SUBROUTINE GPSPLIT(NBFT,NRI,NSM,ZSCATSFN,NBFTARG,NSYM)
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: NSYM
      INTEGER, DIMENSION(*) :: NBFT, NBFTARG, NRI, NSM
      LOGICAL, DIMENSION(*) :: ZSCATSFN
      INTENT (IN) NBFT, NSYM, ZSCATSFN
      INTENT (OUT) NBFTARG, NRI, NSM
!
! Local variables
!
      INTEGER :: I, ISCAT, ISYM, ITARG, L, N
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     GPSPLIT - GT SPLITs the basis functions into target and
!               continuum pieces
!
!***********************************************************************
!
!
!
      N=0
      DO L=1, NSYM
         ISYM=L-1
         iscat=0
         itarg=0
         DO I=1, NBFT(L)
            N=N+1
            NSM(N)=ISYM
            IF(ZSCATSFN(N))THEN
               iscat=iscat+1
               NRI(N)=iscat
            ELSE
               itarg=itarg+1
               NRI(N)=itarg
            END IF
         END DO
         nbftarg(l)=itarg
      END DO
!
      RETURN
!
!---- Format Statements
!
 1000 FORMAT(/10X,'Splitting of Basis Functions',/,10X,&
     &       '----------------------------',//,10X,&
     &       'Number of symmetry adapted B.Fns = ',I5,/,10X,&
     &       'Number of symmetries in set      = ',I5,/)
 2000 FORMAT(10X,'Sym No.  Total  Target ',/)
 2010 FORMAT(10X,I4,5X,I4,3X,I4)
!
      END SUBROUTINE GPSPLIT

      SUBROUTINE NORMAL(NSYM,NBFT,XNORM,IWRITE,NORDINTS,LTRI,ITAPE,NFTTAIL,ZTAIL)
!***********************************************************************
!
!     NORMAL - Computes normalization factors for the symmetry adapted
!              basis functions in each integral.
!
!***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE CONSTS, ONLY : XZERO, XONE
      USE params, ONLY : COVERLAP, MAXSYM
      USE global_utils, ONLY : search
!      USE sword_data, ONLY : VSMALL
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      REAL(KIND=wp), PARAMETER :: VSMALL=1.0D-07 ! used as a tolerance on overlap matrix elements
      INTEGER :: ITAPE, IWRITE, LTRI, NFTTAIL, NORDINTS, NSYM
      LOGICAL :: ZTAIL
      INTEGER, DIMENSION(*) :: NBFT
      REAL(KIND=wp), DIMENSION(MAXSYM,*) :: XNORM
      INTENT (IN) IWRITE, LTRI, NBFT, NORDINTS, NSYM, ZTAIL
      INTENT (INOUT) XNORM
!
! Local variables
!
      INTEGER :: I, IFAIL, IFUNCT, INT, J, L, NUMB, KK
      INTEGER, DIMENSION(ltri) :: INDEXV
      REAL(KIND=wp), DIMENSION(nordints) :: OLAP
      REAL(KIND=wp), DIMENSION(ltri) :: XBUF
!
!*** End of declarations rewritten by SPAG
!
!---- Debug banner header
!
!        WRITE(IWRITE,1000)
!        WRITE(IWRITE,1010) NSYM,(L,NBFT(L),L=1,NSYM)
!        WRITE(IWRITE,1020) NORDINTS,LTRI,ITAPE
!
!---- Initialization of ordered integral buffer
!
      DO I=1, NORDINTS
         OLAP(I)=XZERO
      END DO
!
!---- Find the overlap integrals on the dataset written by MOLECULE.
!     Then read and order these.
!
      CALL SEARCH(ITAPE,COVERLAP,ifail,iwrite)
      IF(ifail.NE.0)STOP 999
! JMC replacing GO TO structure with the following do loop
      DO KK = 1, HUGE(1) ! make sure this loop is technically not infinite
         READ(ITAPE)(XBUF(I),I=1,LTRI), (INDEXV(I),I=1,LTRI), NUMB
!
         IF(NUMB.LT.0)EXIT ! i.e. go to the IF(ZTAIL)THEN line
         DO I=1, NUMB
            OLAP(INDEXV(I))=XBUF(I)
         END DO
      END DO
!
!---- Find the overlap tail integrals and do the subtraction sum if
!     necessary.
!
      IF(ZTAIL)THEN
         REWIND nfttail
         CALL SEARCH(NFTTAIL,COVERLAP,ifail,iwrite)
         IF(ifail.NE.0)STOP 999
! JMC replacing GO TO structure with the following do loop
         DO KK = 1, HUGE(1) ! make sure this loop is technically not infinite
            READ(NFTTAIL)(XBUF(I),I=1,LTRI), (INDEXV(I),I=1,LTRI), NUMB
            IF(NUMB.LT.0)EXIT ! i.e. go to the INT=0 line below
            DO I=1, NUMB
               OLAP(INDEXV(I))=OLAP(INDEXV(I))-XBUF(I)
            END DO
         END DO
      END IF
!
!---- Begin loop over the ordered lower half triangles of integrals
!
!     N.B. this is an assumption about the input !
!
!     We copy out the overlap matrix elements at this stage.
!
      INT=0
      DO L=1, NSYM
         IFUNCT=0
         DO I=1, NBFT(L)
            DO J=1, I
               INT=INT+1
               IF(I.EQ.J)THEN
                  IFUNCT=IFUNCT+1
                  XNORM(L,IFUNCT)=OLAP(INT)
                  IF(OLAP(INT).LT.VSMALL)THEN
                     WRITE(IWRITE,9900)
                     WRITE(IWRITE,9910)L, I, OLAP(INT)
                     STOP 999
                  END IF
               END IF
            END DO
         END DO
      END DO
!
!---- Compute the normalization factors
!
      DO L=1, NSYM
         DO I=1, NBFT(L)
            XNORM(L,I)=XONE/SQRT(XNORM(L,I))
         END DO
      END DO
!
        WRITE(IWRITE,2000)
        DO L=1,NSYM
           WRITE(IWRITE,2010) L
           WRITE(IWRITE,2020) (I,XNORM(L,I),I=1,NBFT(L))
        END DO
!
      RETURN
!
!---- Format Statements
!
 1000 FORMAT(/,10X,'====> NORMAL - COMPUTE FACTORS <====',//)
 1010 FORMAT(10X,'No. of symmetries is orbital set (NSYM) = ',I3,/,10X,'Basis functions per symmetry : ',//,(10X,I3,1X,I3))
 1020 FORMAT(/,10X,'Length of ordered integral vector = ',I6,/,10X,&
     &       'Number of ints per record (LTRI)  = ',I6,/,10X,&
     &       'Unit for overlap ints     (ITAPE) = ',I6,/)
 2000 FORMAT(//,10X,'Normalization factors follow : ',/)
 2010 FORMAT(/,10X,'Symmetry number = ',I3,/)
 2020 FORMAT(10X,I3,1X,F20.10)
 9900 FORMAT(/,10X,'**** Error in NORMAL : ',//)
 9910 FORMAT(/,10X,'One of the sym. adapted basis function ',/,10X,&
     &       'overlap integrals is too small. This   ',/,10X,&
     &       'is probably an error ',//,10X,'Symmetry = ',I3,1X,&
     &       'Function no. = ',I3,/,10X,'Overlap integral = ',F20.10,/)
 9935 FORMAT(/,10X,'Error in an integral index : ',//,10X,&
     &       'Within record sequence = ',I5,/,10X,&
     &       'INDEXV value is        = ',I10,/,10X,&
     &       'Integral itself        = ',F25.13,/)
!
      END SUBROUTINE NORMAL

   !todo supply all allocate statements with error checking.
   !todo replace all calls to 'stop' with calls to 'xermsg'.
   !> This routine reads-in the MOLECULE basis set of symmetry adapted contracted spherical GTOs. The type and amount of output data has been chosen to allow (in ukrmol_basis_data) the conversion of
   !> molecular orbital coefficients from the MOLECULE format to any other format (in our case we convert to the basis of contracted spherical GTOs). However, this routine performs more than just reading of
   !> the MOLECULE basis set. It also normalizes the coefficients for the symmetry adapted functions calculating the self-overlaps of the symmetry-adapted functions. These self-overlaps take into account
   !> the possible presence of tails which are subtracted in that case. It is these normalized symmetry-adapted coefficients which together with the GTO basis set data allow for the conversion of molecular
   !> orbital coefficients between MOLECULE and any other basis. This routine has been written assembling together some routines taken from SWORD and GAUSPROP.
   !> \param[in] iunit Unit number for the unformatted file containing the basis set and raw integrals produced by swmol3 (typically fort.2).
   !> \param[in] nfttail Unit number for the tail integrals produced by gaustail (typically fort.96).
   !> \param[in] ztail Logical variable indicating whether we need to subtract tails or not.
   !> \param[out] CANGULAR Character array specifying the angular type of each primitive GTO in the basis.
   !> \param[out] CHARGNUC Real array containing the charge on each nucleus.
   !> \param[out] CNUCNAME Character array containing the name of the atom and other data concerning the symmetry.
   !> \param[out] CONTCOFNN Unnormalized contraction coefficients for each primitive GTO in the basis.
   !> \param[out] CTRANNN Normalized symmetry-adapted coefficients for the transformation from the contracted GTO basis to the basis of symmetry-adapted functions. See IC for details.
   !> \param[out] EXPONT Exponents of each primitive GTO in the basis.
   !> \param[out] IC Integer array of length LSABFTAB. MOLECULE specifies the transformation from the contracted GTO basis to the basis of symmetry-adapted functions in terms of a table of length LSABFTAB.
   !>                Each row of this table gives the index of the symmetry adapted function (IC), the index of the contracted cartesian GTO contributing to this function (ITRAN) and the coefficient with
   !>                which the GTO contributes to the symmetry adapted function (CTRANNN). This table is printed by the routine GTGETAB.
   !> \param[out] ICONTNO For each primitive GTO this array gives the index of the contracted GTO to which the primitive belongs.
   !> \param[out] IPRCNPT For each contracted GTO this array gives the number of primitive GTOs.
   !> \param[out] ITRAN Integer array of length LSABFTAB. See IC for details.
   !> \param[out] LSABFTAB Integer value giving the dimensions of the arrays IC, ITRAN, CTRANNN.
   !> \param[out] MAX_L Maximum GTO L in the basis.
   !> \param[out] NCONTRAS Number of contracted GTOs in the basis.
   !> \param[out] NGREAD Number of primitive GTOs in the basis.
   !> \param[out] NNUC Total number of nuclei.
   !> \param[out] NSYMABF Total number of symmetry-adapted functions in the basis.
   !> \param[out] NUCIND For each primitive GTO this array gives the index of the nucleus on which the GTO is sitting.
   !> \param[out] XNUC X-coordinate for each nucleus.
   !> \param[out] YNUC Y-coordinate for each nucleus.
   !> \param[out] ZNUC Z-coordinate for each nucleus.
   subroutine READ_MOLECULE_BASIS(iunit,nfttail,ztail,CANGULAR,CHARGNUC,CNUCNAME,CONTCOFNN,CTRANNN,EXPONT,IC,ICONTNO,IPRCNPT,ITRAN,LSABFTAB,MAX_L,NCONTRAS,NGREAD,NNUC,NSYMABF,NUCIND,XNUC,YNUC,ZNUC)
      use precisn
      use const, only: stdout
      use params, only : maxsym, kktyp
      use consts, only : xzero
      use params, only : maxsym, kktyp, LRECL=>LBUF
      use gaus_shared, only : gtindex, gtphiint, maxlp1, ntable
      use global_utils, only : initvr8
      use gausprop_data, only : maxig, maxsfn, iwrite
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: nfttail
      logical, intent(in) :: ztail
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:), intent(out) :: CANGULAR, CNUCNAME
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:), intent(out) :: CHARGNUC, CONTCOFNN, CTRANNN, EXPONT, XNUC, YNUC, ZNUC
      INTEGER, ALLOCATABLE, DIMENSION(:), intent(out) :: IC, ICONTNO, IPRCNPT, ITRAN, NUCIND
      INTEGER, INTENT(OUT) :: LSABFTAB, MAX_L, NCONTRAS, NGREAD, NNUC, NSYMABF
!
! Local variables
!
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: CSYMBOL
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CONTCOF, CTRAN
      CHARACTER(LEN=4) :: CSCATTER
      INTEGER :: I, J, LTRAN, NPRIMS, NSYMT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INUCNO, IXPOWER, IYPOWER, IZPOWER, NRI, NSMI
      INTEGER, DIMENSION(MAXSYM) :: NBFT, NBFTARG
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUCNUMB
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: XNORM
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: ZSCATCFN
      LOGICAL, DIMENSION(maxsfn) :: ZSCATSFN
!ZM
      real(kind=wp), allocatable :: sa_norms(:)
      integer :: nbtot, n1, k
      integer, allocatable :: isymbnum(:)
!
!   
!   ---- Read and print the header information on the file of atomic
!        basis function integrals.
!   
!        It is useful to remember that:
!   
!           LTRAN is the maximum number of transformation functions from
!                 contracted Gaussians to symmetry adapted basis functions
!                 Thus an accurate value would be
!   
!                   Sum over all sym. adapted functions of the
!                                no. of contracted functions
!   
!           NSYMT is the number of Irreducible representations in this
!                 molecular point group i.e. symmetries
!            NBFT is the number of symmetry adapted basis functions
!                 per molecular symmetry.
!   
         CALL PRTINF(IUNIT,STDOUT,NSYMT,NBFT,KKTYP(:,2),LTRAN)
!   
!        NSYMABF counts the total number of symmetry adapted basis
!        functions here.
!   
         NSYMABF=0
         DO I=1, NSYMT
            NSYMABF=NSYMABF+NBFT(I)
         END DO
!   
!   ---- Read the nuclear location table from the file of atomic integrals
!   
!        This identifies all of the nuclei in the system. The designation
!        of the nuclei is related to the MOLECULE input format i.e. in terms
!        of symmetry independent atoms. Remember that each symmetry
!        independent atom must be given a unique label such as C1 and C2 to
!        distinguish between two Carbon atoms which cannot be related to
!        each other by a symmetry operation of the molecular point group
!        (D2h or lower). MOLECULE then generates all actual nuclei from the
!        groups operations and assigns a sequential index to each case for
!        each symmetry independent nucleus.  It is is this information that
!        is stored in the columns
!   
!           CNUCNAME and CNUCNUMB
!   
!           XNUC,YNUC and ZNUC are the X,Y,Z co-ordinates while
!   
!           CHARGNUC is the charge on the nucleus
!   
!        Later in the code each primitive Gaussian function is identified
!        with one row in this table.
!   
         REWIND IUNIT
!   
         READ(IUNIT)
         READ(IUNIT)
         READ(IUNIT)NNUC
!   
         ALLOCATE(CNUCNAME(NNUC),NUCNUMB(NNUC), XNUC(NNUC),YNUC(NNUC),ZNUC(NNUC),CHARGNUC(NNUC))
!   
         DO I=1, NNUC
            READ(IUNIT)CNUCNAME(I), NUCNUMB(I), XNUC(I), YNUC(I), ZNUC(I), CHARGNUC(I)
         END DO
!   
!   ---- Read the Gaussian basis function information from the integrals
!        file now.
!   
!            NPRIMS - Number of primitive Gaussian functions
!          NCONTRAS - Number of contracted Gaussian functions
!   
!        Starting at position IPRCNPT of the big vector is a table of
!        the number of primitives per contracted function.
!   
         READ(IUNIT)NPRIMS, NCONTRAS
         ALLOCATE(iprcnpt(ncontras),zscatcfn(ncontras))
         READ(IUNIT)(IPRCNPT(I),I=1,NCONTRAS)
!   
!   .... Find the total number of primitive Gaussians
         NGREAD=0
         DO I=1, NCONTRAS
            NGREAD=NGREAD+IPRCNPT(I)
         END DO
!   
!   ...... Allocate storage space from the big vector for the primitive
!          basis function table.
!   
         ALLOCATE(icontno(ngread),inucno(ngread),expont(ngread),&
        &         contcof(ngread),contcofnn(ngread),nri(nsymabf),nsmi(nsymabf),ic(ltran),&
        &         itran(ltran),ctran(ltran),ctrannn(ltran),csymbol(ngread),isymbnum(ngread),&
        &         cangular(ngread))
!   
!   ...... Invoke the routine which reads the table of primitive Gaussians
!          and their identifications into core
!   
         CALL GTGETAB(stdout,IUNIT,NCONTRAS,IPRCNPT,ICONTNO,CSYMBOL,ISYMBNUM,InucNO,&
        &             CANGULAR,cscatter,nnuc,cnucname,EXPONT,CONTCOF,CONTCOFNN,&
        &             NSYMABF,IC,ITRAN,CTRAN,LTRAN,LSABFTAB,ZSCATCFN,&
        &             ZSCATSFN,MAXSFN)
!   
         CALL GPSPLIT(NBFT,NRI,NSMI,ZSCATSFN,NBFTARG,NSYMT)
!   
!-----   Each primitive Gaussian has associated with it an angular term in
!        the form x^i * y^j * z^k and this is represented in MOLECULE by
!        a character string. The string for each basis function is stored
!        in the array beginning at IANGULR. This is now converted into
!        an integer representation.
!   
!        It is useful to remember that:
!   
!           IXPOWER is the starting location for the storage of i powers
!           IYPOWER is  "   " " "    " " "    "   "  storage of j powers
!           IZPOWER is  "   " " "    " " "    "   "  storage of k powers
!   
         ALLOCATE(ixpower(ngread),iypower(ngread),izpower(ngread))
!   
         max_l = 0
         DO I=1, NGREAD
            CALL GTINDEX(CANGULAR(I),IXPOWER(I),IYPOWER(I),IZPOWER(I))
            if (IXPOWER(I)+IYPOWER(I)+IZPOWER(I) > max_l) max_l = IXPOWER(I)+IYPOWER(I)+IZPOWER(I)
         END DO
!   
!
!---- Determine absolute (i.e. 1..NNUC) indices of the nuclei where each of the primitive GTOs are sitting.
!
      allocate(NUCIND(NGREAD))
!
      do i=1,NGREAD
         do k=1,NNUC
            if  (CNUCNAME(k) .eq. CSYMBOL(i) .and. NUCNUMB(k) .eq. ISYMBNUM(i)) then
                NUCIND(i) = k
            endif
         enddo !k
      enddo !i
!
!---- Compute : the total number of basis function, NBTOT;
!               the total storage for all lower half triangles, N1;
!               pointers into each symmetries lower triangle, NZBF. (JMC set but not used so removed)
!
      NBTOT=0
      N1=0
!
      DO I=1, NSYMT
         N1=N1+NBFT(I)*(NBFT(I)+1)/2
         NBTOT=NBTOT+NBFT(I) 
      END DO

!
!     There are LRECL/2 integrals and LRECL/2 indices per record ! JMC ??? don't believe that...
!
      ALLOCATE(xnorm(maxsym,nbtot),sa_norms(NBTOT))
!
!---- Invoke the subroutine which computes symmetry adapted basis
!     function normalization factors. These are obtained from the
!     overlap integrals.
!
      CALL NORMAL(NSYMT,NBFT,XNORM,stdout,N1,lrecl,IUNIT,NFTTAIL,ZTAIL)
!     
!---- Normalize the symmetry adapted functions using the XNORM factors calculated by NORMAL. The orbital coefficients in MOLECULE are expressed in terms of 
!     coefficients for the NORMALIZED symmetry adapted fns.: the normalization of all molecular integrals is performed in SWORD using the XNORM values. This normalization is equivalent to normalizing the 
!     symmetry adapted functions themselves. This is what we need to do if we want to convert the molecular orbital coefficients from the MOLECULE basis to any other basis.
!
      K = 0
      DO I=1, NSYMT
         DO J=1, NBFT(I)
            K = K + 1
            sa_norms(k) = XNORM(I,J)
         ENDDO
      END DO
!
!     Normalize the symmetry-adapted expansion coefficients
!
      do k=1,LSABFTAB
         CTRANNN(K) = CTRAN(K)*sa_norms(IC(k))
      enddo
! 
   end subroutine READ_MOLECULE_BASIS

end module ukrmol_routines
