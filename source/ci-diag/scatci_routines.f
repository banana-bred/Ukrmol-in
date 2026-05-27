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
      MODULE SCATCI_ROUTINES

      public

      CONTAINS

!*==check.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE CHECK(NH,THRES,CH)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NH
      REAL(KIND=wp) :: THRES
      REAL(KIND=wp), DIMENSION(*) :: CH
      INTENT (IN) CH, THRES
      INTENT (INOUT) NH
C
C Local variables
C
      INTEGER :: MH, N
C
C*** End of declarations rewritten by SPAG
C
      MH=0
      DO N=1, NH
         IF(ABS(CH(N)).GE.THRES)MH=MH+1
      END DO
      NH=MH
      RETURN
      END SUBROUTINE CHECK
!*==chn2e.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE chn2e(nftin,iwrite,xint2e,nint2e)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE params, ONLY : LRECIN=>LRECLX, ctrans2
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, NFTIN, NINT2E
      REAL(KIND=wp) :: XINT2E(:)
      INTENT (IN) IWRITE, NINT2E
      INTENT (OUT) XINT2E
C
C Local variables
C
      INTEGER :: I, ICOUNT, IFAIL, INBOX, LAST
C
C*** End of declarations rewritten by SPAG
C
c***********************************************************************
c
c     Read the two electron part of the file of unformatted transformed
c     integrals from MOLECULE-SWEDEN and copy these to a new format
c     which is Alchemy I compatible.
c
c     Input Data:
c          nftin  logical unit holding MOLECULE-SWEDEN integrals
c         iwrite  logical unit for the printer
c     Output Data:
c         xint2e  array for the integrals
c
c     Notes: Alchemy requires to read integrals in boxes of preset
c            size because that is the way that the formulae are sorted.
c            Sweden on the other hand writes all integrals in one long
c            vector split into equal length records. This routine
c            is a simple interface between the two codes which
c            rearrangesthe Sweden integrals into records which are the
c            same size as an Alchemy box. Thus routine procbs can simpl
c            use a vectorized read of a single vector. It is clearly
c            possible to enhance the speed of this coding at the expense
c            of complicating the do loops.
c
c***********************************************************************
C
c---- Rewind transformed integrals file and start of the two electron
c     integrals
c
      REWIND nftin
c
      CALL search(nftin,ctrans2,ifail)
c
c---- Loop until end of file reading records of length lrecin
c     words each time.
c
c.... Initialize for the loop
      icount=1
      inbox=0
c
c.... Begin a Fortran implementation of a repeat until loop
c     which returns to 300 until end file forces jump to 400
c
      last=min(lrecin,nint2e)
 300  READ(nftin,END=400,ERR=990)(xint2e(inbox+i),i=1,last)
      icount=icount+1
      inbox=inbox+lrecin
      IF(inbox+lrecin.GT.nint2e)last=nint2e-inbox
      GO TO 300
c
c---- End of file condition reached, jump out of the loop
c
 400  WRITE(iwrite,4000)icount-1
 4000 FORMAT(/,10x,'End of file reached. Swed.Record count = ',i10,//)
c
      REWIND nftin
c
      RETURN
c
c---- Abnormal termination on two electron integral reading
c
 990  WRITE(iwrite,9900)icount
 9900 FORMAT(/,10x,'Error on read for 2-e ints. Record = ',i10,//)
      STOP
c
      END SUBROUTINE CHN2E
!*==cirmat.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE CIRMAT(NTGsym,NFTG,NTGTF,CTGT,etgt,NTGTS,NCTGT,NFT,
     &                  numtgt,iphz,nctarg,symtyp)
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : ZERO=>XZERO, ONE=>XONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NCTARG, NFT, NFTG, NTGSYM, SYMTYP
      REAL(KIND=wp), DIMENSION(*) :: CTGT, ETGT
      INTEGER, DIMENSION(*) :: IPHZ, NCTGT, NTGTF, NTGTS, NUMTGT
      INTENT (IN) NCTARG, NTGSYM, NTGTF, NUMTGT, SYMTYP
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(41) :: DTNUC1
      REAL(KIND=wp) :: S, SZ
      INTEGER :: MGVN, NELT
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: EIG, VC
      INTEGER :: IFL, IFNUM, II, IST, ITGT, JTGT, LST, NHDIM, NOCSF1, 
     &           NPHZ, NSTAT, ERR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KPHZ
      CHARACTER(LEN=120) :: NAME
      INTEGER, DIMENSION(10) :: NHD1
      INTEGER, DIMENSION(20) :: NHE1
C
C*** End of declarations rewritten by SPAG
C
C     PROCESS DATA FOR HAMILTONIAN CONTRACTION IN R-MATRIX CI
C
      IFNUM=0
      ist=0
      lst=0
      ifl=0
      nphz=1
C
C     OBTAIN CI COEFFICIENTS FROM DUMPFILE FOR H-TRANSFORMATION
C
      DO Itgt=1, ntgsym
         nstat=numtgt(itgt)
         DO ii=1, numtgt(itgt)
            nstat=max(nstat,NTGTS(IFL+ii))
         END DO
         DO jtgt=1, numtgt(itgt)
            ifl=ifl+1
C
C     DETERMINE FILE NUMBER AND POSITION AT START OF DATASET
C
            IF(NTGTF(IFL).GT.0)THEN
               WRITE(NFT,2182)IFL, NTGTF(IFL), NTGTS(IFL)
c
C     Do we need to find another set of target vectors ?
               IF(IFNUM.NE.NTGTF(IFL))THEN
                  IFNUM=NTGTF(IFL)
C
                  nhdim=nctgt(itgt)
                  IF (ALLOCATED(kphz)) DEALLOCATE(kphz)
                  IF (ALLOCATED(eig)) DEALLOCATE(eig)
                  IF (ALLOCATED(vc)) DEALLOCATE(vc)
                  ALLOCATE(kphz(nhdim),eig(nhdim),vc(nhdim*nstat),
     &                     stat=ERR)
                  IF (ERR.ne.0) THEN
                     print *,'CIRMAT: memory allocation error'
                     stop
                  ENDIF
C
C     READ CI COEFFICIENTS
C
                  IF(symtyp.LT.2)THEN
                     CALL READCID(nftg,ifnum,NAME,NHE1,NHD1,DTNUC1,
     &                            NOCSF1,NSTAT,nctgt(itgt),EIG,VC,kphz,
     &                            NFT)
                  ELSE
                     CALL READCIP(nftg,ifnum,NOCSF1,NSTAT,nctgt(itgt),
     &                            mgvn,s,sz,nelt,EIG,VC,kphz,NFT)
                  END IF
               END IF
C            ELSE
C
C     THIS IS AN SCF WAVEFUNCTION
C               WRITE(NFT,2183)IFL
C               IF (ALLOCATED(kphz)) DEALLOCATE(kphz)
C               IF (ALLOCATED(eig)) DEALLOCATE(eig)
C               IF (ALLOCATED(vc)) DEALLOCATE(vc)
C               ALLOCATE(kphz(1),eig(1),vc(1))
C               vc(1)=ONE
C               NOCSF1=1
C               eig(1)=zero
C
            END IF
C
C     SELECT AND PACK CI COEFFICIENTS FOR H-TRANSFORMATION
C
            IF(nctarg.LE.0)THEN
               CALL CPAK(NTGTS(IFL),NOCSF1,VC,eig,CTGT,etgt,ist,lst)
            ELSE
C        .... also introducing phase correction from CONGEN
C
               CALL CPAKPZ(NTGTS(IFL),NOCSF1,VC,eig,CTGT,etgt,ist,lst,
     &                     iphz(nphz),kphz)
            END IF
C
         END DO
         nphz=nphz+nocsf1
      END DO
C     successful completion
      IF(nphz-1.EQ.nctarg .OR. nctarg.EQ.0)RETURN
C
      WRITE(nft,2180)nphz-1, nctarg
 2180 FORMAT(/' Target data from CI file not consistent with data',
     &       ' from CONGEN',/' Total number of target CSFs here =',i6,
     &       ' from CONGEN =',i6)
C
      STOP
 2182 FORMAT(/' TRANSFORMATION VECTOR ',I3,' FROM SET',I3,' STATE',I3)
 2183 FORMAT(/' TRANSFORMATION VECTOR ',I3,' IS FOR AN SCF STATE')
C
      END SUBROUTINE CIRMAT
!*==civio.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE CIVIO(NFT,NRW,NK,NS,EI,CV,NALM,iphz,dg)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NALM, NFT, NK, NRW, NS
      REAL(KIND=wp), DIMENSION(nk,ns) :: CV
      REAL(KIND=wp), DIMENSION(nk) :: DG
      REAL(KIND=wp), DIMENSION(ns) :: EI
      INTEGER, DIMENSION(nk) :: IPHZ
      INTENT (IN) DG, NFT, NK, NRW, NS
      INTENT (OUT) NALM
      INTENT (INOUT) CV, EI, IPHZ
C
C Local variables
C
      INTEGER :: I, J, M
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
C     Read/Write ENERGY AND CSF SPECIFICATION FOR EACH CI STATE
C     annd then Read/Write COEFFICIENTS FOR EACH CI STATE
      NALM=0
      IF(NRW.EQ.0)THEN
         WRITE(NFT,ERR=200)iphz, EI, dg
         DO i=1, ns
            WRITE(NFT,ERR=200)i, (CV(j,i),j=1,nk)
         END DO
      ELSE
         READ(NFT,ERR=200)iphz, EI
         DO i=1, ns
            READ(NFT,ERR=200)M, (CV(j,i),j=1,nk)
         END DO
      END IF
      RETURN
C
 200  nalm=1
      RETURN
      END SUBROUTINE CIVIO
!*==comprs.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE COMPRS(nci,nold,nnew,XMJK,CJK,xmnew,cnew,coef,nmax)
c
c     Comprs takes latest integral contribution to CI target and
c     compresseswith previous ones: more than one coefficent is require
c     if there is more than one target state of the same symmetry.
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NCI, NMAX, NNEW, NOLD
      REAL(KIND=wp), DIMENSION(nmax,nci) :: CJK
      REAL(KIND=wp), DIMENSION(nnew) :: CNEW
      REAL(KIND=wp), DIMENSION(nci) :: COEF
      INTEGER(longint), DIMENSION(2,nmax) :: XMJK
      INTEGER(longint), DIMENSION(2,nnew) :: XMNEW
      INTENT (IN) CNEW, COEF, NCI, NMAX, NNEW, XMNEW
      INTENT (INOUT) CJK, NOLD
C
C Local variables
C
      INTEGER :: I, IJK, N, NJK
      INTEGER(longint), DIMENSION(2) :: XNA
C
C*** End of declarations rewritten by SPAG
C
      njk=nold
      ijk=0
c     loop over new integrals
      DO i=1, nnew
         xna(1)=xmnew(1,i)
         xna(2)=xmnew(2,i)
c     is this one already stored?
         ijk=ijkpqrs(xna,xmjk,nold)
         IF(IJK.EQ.0)THEN
C     if not present: add to the end of the list
            NJK=NJK+1
            IF(NJK.LE.NMAX)THEN
               XMJK(1,NJK)=XNA(1)
               XMJK(2,NJK)=XNA(2)
               DO n=1, nci
                  CJK(NJK,n)=cnew(i)*coef(n)
               END DO
            ELSE
               WRITE(6,999)nmax
 999           FORMAT(/
     &             ' **** Attempt to store too many integrals in COMPRS'
     &             ,/,5x,'Present limit of',i7,' integrals exceeded')
               STOP
            END IF
         ELSE
C     if present: add coefficients
            DO n=1, nci
               CJK(IJK,n)=CJK(IJK,n)+cnew(i)*coef(n)
            END DO
            IJK=0
         END IF
      END DO
c     the number of integrals now stored:
      nold=njk
      RETURN
      END SUBROUTINE COMPRS
!*==cpak.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE CPAK(K,NOSCF,CIV,eig,C,e,i,l)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: I, K, L, NOSCF
      REAL(KIND=wp), DIMENSION(*) :: C, E, EIG
      REAL(KIND=wp), DIMENSION(NOSCF,*) :: CIV
      INTENT (IN) CIV, EIG, K, NOSCF
      INTENT (OUT) C, E
      INTENT (INOUT) I, L
C
C Local variables
C
      INTEGER :: J
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     CPAK ACCUMULATES THE CI EIGENVECTOR EXPANSION COEFFICIENTS
C          FOR THE TARGET EIGENSTATES IN A PACKED FORM SUITABLE
C          for CI contraction in subroutine ENRGMX
C
C***********************************************************************
C
C
      l=l+1
      e(l)=eig(k)
C
      DO J=1, NOSCF
         I=I+1
         C(I)=CIV(J,K)
      END DO
C
      RETURN
C
      END SUBROUTINE CPAK
!*==cpakpz.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE CPAKPZ(K,NOSCF,CIV,eig,C,e,i,l,iphz,iphz0)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: I, K, L, NOSCF
      REAL(KIND=wp), DIMENSION(*) :: C, E, EIG
      REAL(KIND=wp), DIMENSION(NOSCF,*) :: CIV
      INTEGER, DIMENSION(noscf) :: IPHZ, IPHZ0
      INTENT (IN) CIV, EIG, IPHZ, IPHZ0, K, NOSCF
      INTENT (OUT) C, E
      INTENT (INOUT) I, L
C
C Local variables
C
      INTEGER :: J
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     CPAKPZ ACCUMULATES THE CI EIGENVECTOR EXPANSION COEFFICIENTS
C          FOR THE TARGET EIGENSTATES IN A PACKED FORM SUITABLE
C          for CI contraction in subroutine ENRGMX
C     Also include the phase correction contained in array iphz
C
C***********************************************************************
C
C
      l=l+1
      e(l)=eig(k)
      DO J=1, NOSCF
         I=I+1
         C(I)=CIV(J,K)*REAL(iphz(j)*iphz0(j),KIND=wp)
      END DO
C
      RETURN
C
      END SUBROUTINE CPAKPZ
!*==diff.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE DIFF(NCH0,NCH,XNPQ,CPQ,CH)
C
C     INPUT
C      NCH  = NO. OF NON-ZERO INTEGRALS OF THE REFERENCE SET
C      NPQ  = INDICES OF INTEGRALS
C      CPQ  = VALUES OF INTEGRALS
C      THRES= THRESHOLD FOR ZERO
C      CH   =ORDERED SET OF CURRENT ELEMENTS
C
C      OUTPUT
C      NCHP =NO. OF NON/ZERO INTEGRALS OF THE CURRENT SET
C      NPQP =INDICES OF INTEGRALS
C      CPQP =VALUES OF INTEGRALS
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NCH, NCH0
      REAL(KIND=wp), DIMENSION(*) :: CH, CPQ
      INTEGER(longint), DIMENSION(2,*) :: XNPQ
      INTENT (IN) CPQ, NCH, NCH0
      INTENT (INOUT) CH
C
C Local variables
C
      INTEGER :: I, IC, JA, JB
      INTEGER, DIMENSION(8) :: ITEMP
C
C*** End of declarations rewritten by SPAG
C
      DO I=NCH0+1, NCH
         CALL unpack8ints(XNPQ(1,I),ITEMP)
         JA=ITEMP(1)
         JB=ITEMP(4)
         IC=(JA*(JA-1))/2+JB
         CH(IC)=CH(IC)-CPQ(I)
      END DO
C
      RETURN
      END SUBROUTINE DIFF
!*==enrgmb.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE ENRGMB(CJA,CJB,CKA,CKB,MN,MM,MS,MPOS,NI,CFD,XMJK,CJK,
     &                  SYMTYP, NOEX)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE consts, ONLY : ONE=>XONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: LPOSIT, MFLG, NFLG, NH
      COMMON /ENA   / NH, NFLG, MFLG
      COMMON /GD4   / LPOSIT
C
C Dummy arguments
C
      REAL(KIND=wp) :: CFD
      INTEGER :: SYMTYP, NOEX
      REAL(KIND=wp), DIMENSION(*) :: CJA, CJB, CJK, CKA, CKB
      INTEGER, DIMENSION(*) :: MM, MN, MPOS, MS
      INTEGER, DIMENSION(4) :: NI
      INTEGER(longint), DIMENSION(2,*) :: XMJK
      INTENT (IN) CFD, MM, MN, MPOS, MS, NI, SYMTYP
      INTENT (INOUT) CJA, CJB, CKA, CKB
C
C Local variables
C
      INTEGER, SAVE :: I, KA, KB, M1, M2, M3, M4, MLA, MLB, MP1, MP2, 
     &                 MP3, MP4, N, NP, NQ, NR, NS, NS1, NS2, NS3, NS4, 
     &                 NSFA, NSFB
      INTEGER, DIMENSION(4), SAVE :: MJ, MPS, NJ, NSJ
      REAL(KIND=wp), SAVE :: SIGN
C
C*** End of declarations rewritten by SPAG
C
      EQUIVALENCE(NJ(1),NP)
      EQUIVALENCE(NJ(2),NQ)
      EQUIVALENCE(NJ(3),NR)
      EQUIVALENCE(NJ(4),NS)
      EQUIVALENCE(MJ(1),M1)
      EQUIVALENCE(MJ(2),M2)
      EQUIVALENCE(MJ(3),M3)
      EQUIVALENCE(MJ(4),M4)
      EQUIVALENCE(NSJ(1),NS1)
      EQUIVALENCE(NSJ(2),NS2)
      EQUIVALENCE(NSJ(3),NS3)
      EQUIVALENCE(NSJ(4),NS4)
      EQUIVALENCE(MPS(1),MP1)
      EQUIVALENCE(MPS(2),MP2)
      EQUIVALENCE(MPS(3),MP3)
      EQUIVALENCE(MPS(4),MP4)
C
C     SET UP INPUT PARAMETERS
C
      DO I=1, 4
         N=NI(I)
         MPS(I)=MPOS(N)
         NJ(I)=MN(N)
         MJ(I)=MM(N)
         NSJ(I)=MS(N)
      END DO
C
C     TEST SPIN AND M VALUE
C
      IF(NS1+NS2.EQ.1 .OR. NS3+NS4.EQ.1 .OR. MP1.NE.MP2)THEN
         NSFA=0
      ELSE
         NSFA=1
         MLA=0
         IF(SYMTYP.LT.2)THEN
            KA=MAX(NI(1),NI(2))
            KB=MAX(NI(3),NI(4))
            IF(KA.LT.KB)THEN
               IF(M3*M4.LT.0)MLA=1
            ELSE
               IF(M1*M2.LT.0)MLA=1
            END IF
         END IF
      END IF
C NO EXCHANGE FOR POSITRONS
      IF(NS1+NS4.EQ.1 .OR. NS2+NS3.EQ.1 .OR. MP1.NE.MP4)THEN
         NSFB=0
      ELSE
         NSFB=1
         MLB=0
         IF(SYMTYP.LT.2)THEN
            KA=MAX(NI(1),NI(4))
            KB=MAX(NI(3),NI(2))
            IF(KA.LT.KB)THEN
               IF(M2*M3.LT.0)MLB=1
            ELSE
               IF(M1*M4.LT.0)MLB=1
            END IF
         END IF
      END IF
C NO COULOMB OR EXCHANGE -- RETURN
      IF(NSFA.EQ.0 .AND. NSFB.EQ.0)RETURN
C CHANGE SIGN FOR POSITRONS
      SIGN=ONE
      IF(LPOSIT.EQ.1)SIGN=-ONE
      if (NOEX .eq. 1) sign=ONE
C
C     TEST FOR TYPE
C
      IF(MFLG.EQ.1)GO TO 400
      IF(NP.EQ.NQ .AND. NR.EQ.NS)GO TO 100
      IF(NP.EQ.NQ .OR. NR.EQ.NS)GO TO 400
      IF(NP.EQ.NR .AND. NQ.EQ.NS)GO TO 300
      IF(NP.EQ.NR)GO TO 400
      IF(NP.EQ.NS .AND. NQ.EQ.NR)GO TO 200
      GO TO 400
C
C     CASE 1   (PP/RR) AND (PR/PR)
C
 100  IF(NP.LT.NR)THEN
         N=(NR*(NR-1))/2+NP
      ELSE
         N=(NP*(NP-1))/2+NR
      END IF
      IF(NSFA.NE.0)THEN
         IF(MLA.EQ.0)THEN
            CJA(N)=CJA(N)+SIGN*CFD
         ELSE
            CJB(N)=CJB(N)+SIGN*CFD
         END IF
      END IF
C
      IF(NSFB.EQ.0)RETURN
      IF(MLB.EQ.0)THEN
         CKA(N)=CKA(N)-SIGN*CFD
      ELSE
         CKB(N)=CKB(N)-SIGN*CFD
      END IF
      RETURN
C
C     CASE 2   (PQ/QP) AND (PP/QQ)
C
 200  IF(NP.LT.NQ)THEN
         N=(NQ*(NQ-1))/2+NP
      ELSE
         N=(NP*(NP-1))/2+NQ
      END IF
      IF(NSFA.NE.0)THEN
         IF(MLA.NE.0)THEN
            CKB(N)=CKB(N)+SIGN*CFD
         ELSE
            CKA(N)=CKA(N)+SIGN*CFD
         END IF
      END IF
C
      IF(NSFB.EQ.0)RETURN
      IF(MLB.NE.0)THEN
         CJB(N)=CJB(N)-SIGN*CFD
      ELSE
         CJA(N)=CJA(N)-SIGN*CFD
      END IF
      RETURN
C
C     CASE 3   (PQ/PQ)
C
 300  IF(NP.LT.NQ)THEN
         N=(NQ*(NQ-1))/2+NP
      ELSE
         N=(NP*(NP-1))/2+NQ
      END IF
      IF(NSFA.NE.0)THEN
         IF(MLA.NE.0)THEN
            CKB(N)=CKB(N)+SIGN*CFD
         ELSE
            CKA(N)=CKA(N)+SIGN*CFD
         END IF
      END IF
C
      IF(NSFB.EQ.0)RETURN
      IF(MLB.NE.0)THEN
         CKB(N)=CKB(N)-SIGN*CFD
      ELSE
         CKA(N)=CKA(N)-SIGN*CFD
      END IF
      RETURN
C
C     THE REST OF THE CASES
C
 400  IF(NSFA.NE.0)CALL PQRS(NP,NQ,NR,NS,MLA,0,SIGN*CFD,XMJK,CJK)
      IF(NSFB.NE.0)CALL PQRS(NP,NS,NR,NQ,MLB,1,SIGN*CFD,XMJK,CJK)
      RETURN
      END SUBROUTINE ENRGMB
!*==enrgmc.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE ENRGMC(CJK,NORBB)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : ZERO=>XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NORBB
      REAL(KIND=wp), DIMENSION(norbb,4) :: CJK
      INTENT (IN) NORBB
      INTENT (OUT) CJK
C
C Local variables
C
      INTEGER :: I, J
C
C*** End of declarations rewritten by SPAG
C
      DO j=1, 4
         DO I=1, norbb
            CJK(I,j)=ZERO
         END DO
      END DO
C
      RETURN
      END SUBROUTINE ENRGMC
!*==enrgms.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE ENRGMS(NORB,MN,MM,MS,NELT,NDTRF,MDTR,NOCSF,NODT,INDT,
     &                  NDT,ICDT,CDT,THRES,NALM,IDIAG,nelm,NFTW,IPOSIT,
     &                  MAP,MPOS,mcon,mcorb,NSRB,SYMTYP,ncont,ncont2,
     &                  nfte,lembf,xint1e,nint1e,xint2e,thrhm,IODR,NRI,
     &                  NSM,NOB,MBAS,NBAS,NBASH,IPAIR,NCORB,kpt,mocsf,
     &                  iexpc,eig,notgt,nummx,LUSME,enh_factor,
     &                  enh_factor_type,ukrmolp_ints,NOEX)
c
c    ENRGMS controls the construction and evaluation of matrix elements
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE consts, ONLY : ZERO=>XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: MALM, MFLG, NFLG, NH, NJK, NJKMX
      COMMON /ENA   / NH, NFLG, MFLG
      COMMON /STPQ  / NJK, MALM, NJKMX
C
C Dummy arguments
C
      INTEGER :: IDIAG, IEXPC, IODR, IPOSIT, LEMBF, LUSME, MOCSF, NALM, 
     &           NCONT, NCONT2, NCORB, NELM, NELT, NFTE, NFTW, NINT1E, 
     &           NOCSF, NORB, NOTGT, NSRB, NUMMX, SYMTYP,
     &           enh_factor_type, NOEX
      REAL(KIND=wp) :: THRES, THRHM, enh_factor
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CDT, EIG
      REAL(KIND=wp), ALLOCATABLE :: XINT1E(:), XINT2E(:)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDT, NDT, ICDT, MAP, MCON,
     &                                      MCORB, MDTR, MM, MN, MPOS,
     &                                      MS, NODT, NRI
      INTEGER, DIMENSION(:) :: MBAS, NBAS, NBASH, NDTRF, NOB, NSM
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KPT
      LOGICAL :: ukrmolp_ints
      INTENT (IN) EIG, ICDT, IEXPC, INDT, MOCSF, NCONT, NCONT2, NFTW, 
     &            NOCSF, NSRB, NUMMX, THRHM, enh_factor,enh_factor_type,
     &            ukrmolp_ints
      INTENT (OUT) NALM
      INTENT (INOUT) NELM
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CCJK, CHC, CHR, ELEM
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CH, CHT, EM
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: CJK
      REAL(KIND=wp) :: EONE
      CHARACTER(LEN=32) :: HEADERP
      INTEGER :: I, ICH, IEMBF, IERR, IJ, IKK, ISTEP, ISTEP1, ISTEP2, J, 
     &           JDIAG, KCC, M, MA0, MA00, MA000, MC, MD, MMM, MST, 
     &           MXCIN, MXCINT, N, NA0, NA00, NBEG, NC, 
     &           NCI, NCL, ND, NEED, NELEM, 
     &           NELM0, NELMX, NHC, NHR, NHT, NINT, NINT0, NINT00, 
     &           NORBB, NORBL, NTLIC, NTLIR, NTLIT, NZ
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBB
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NEMBF
      INTEGER, DIMENSION(5) :: NRG
      INTEGER(longint), ALLOCATABLE, DIMENSION(:,:) :: XMH
      INTEGER(longint), ALLOCATABLE, DIMENSION(:,:) :: XMHC, XMHC1, 
     &                XMHR, XMJK
C
C*** End of declarations rewritten by SPAG
C
      EQUIVALENCE(NRG(1),NHR)
      DATA headerp/'************************NUMOFP  '/
C
      ALLOCATE(CH(norb*(norb+1)/2), CHT(norb), CJK(norb*(norb+1)/2,4))
      ALLOCATE(EM(lembf), NBB(nsrb), NEMBF(2,LEMBF), XMH(2,NORB))
C
      ierr=0
      NORBB=(NORB*(NORB+1))/2
      NORBL=NORB
      NALM=0
C-----------------------------------------------------------------------
      mxcin=4*NORBB+NORBL
      mxcint=2*mxcin
c
c *** xmhc needs to be bigger than mxcint in EXPODL
c
      ALLOCATE(chc(mxcin),xmhc(2,mxcint),xmhc1(2,mxcint),stat=ierr)
      need=mxcin+2*mxcint
      IF(ierr.NE.0)THEN
         WRITE(6,99)need
 99      FORMAT(/' UNABLE TO ALLOCATE MEMORY IN ENRGMX, need extra',i10,
     &          ' (real) words')
         STOP
      END IF
C
C    This is a best guess, although errors are flagged in lower level
C    routines, they are not always trapped in this rotuine
      NJKMX=mxcint
      nelmx=1
      ALLOCATE(elem(nelmx),xmjk(2,njkmx),ccjk(njkmx),xmhr(2,mxcin),
     &         chr(mxcin))
C-----------------------------------------------------------------------
      NFLG=0
      MFLG=0
C     ZERO ARRAY FOR PZERO
      NBB=0
C-----------------------------------------------------------------------
C
c     First construct the (1,1) diagonal element
c
C-----------------------------------------------------------------------
      NJK=0
      MALM=0
      nelm=0
      nelem=1
      ma0=0
      m=1
      iembf=0
      IF(iexpc.EQ.0)THEN
         istep=1
         istep1=1
         istep2=1
      ELSE
         istep=2
         istep1=2
         istep2=2
      END IF
      CALL ENRGMC(CJK,NORBB)
      CALL ENRGMZ(CH,xMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,IPOSIT,
     &            MPOS,NODT(1),CDT(1),NDT(1),NODT(1),CDT(1),NDT(1),CJK,
     &            xMJK,CCJK,SYMTYP,mcon,idiag,NOEX)
      CALL CHECK(NH,THRES,CHT)
      CALL CHECK(NJK,THRES,CCJK)
      IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
         MALM=0
         WRITE(NFTW,9010)
 9010    FORMAT('    THE 1ST CSF IS IN ERROR')
         WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
 9020    FORMAT('    NH =',I3,' NJK =',I3,' MALM =',I3,' NJKMX =',I7)
         NALM=1
         RETURN
      END IF
C
      CALL STDIAG(CH,CHR,xMHR,CHR,xMHR,NORB,NORBL,NHR,NHC,THRES,1,mpos)
      NTLIR=NHR
C
      DO N=1, 4
         CALL STORE(NTLIR,xMHR,CHR,NORB,THRES,CJK(1,n),N)
         NRG(N+1)=NTLIR
      END DO
C
      NJK=NTLIR-NHR
      nint0=ntlir
      ntlit=ntlir
      nci=1
      ncl=1
      IF(iexpc.NE.0)THEN
c
c        expand prototype diagonal CSF
c
         nelem=notgt
         CALL expdg(kpt(m),nelem,ntlit,xmhr,xmhr,ntlit)
      END IF
      nint00=ntlit
c
c     evaluate the integrals
c
      kcc=1
      ikk=1
      iembf=0
      DO i=1, nelem
         ma0=i-notgt
         IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, nci
         CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlit,
     &               xmhr(1,ikk),Chr,mxcint,xint1e,nint1e,xint2e,elem,
     &               symtyp,map,ncorb,mcorb,nci,LUSME,enh_factor,
     &               enh_factor_type,ukrmolp_ints,NOEX)
         KCC=KCC+1
         IF(LUSME.NE.0)WRITE(LUSME,*)i, i, 0
         IF(i.EQ.1 .AND. idiag.EQ.0)THEN
            eone=elem(1)
            elem(1)=elem(1)-eone
         END IF
         ij=1
         IF(idiag.EQ.0)THEN
            elem(ij)=elem(ij)+eone
         ELSE IF(idiag.EQ.2)THEN
            elem(ij)=elem(ij)+eig(1)
         END IF
         ma0=ma0+notgt
         iembf=iembf+1
         na0=ma0
         nembf(1,iembf)=na0
         nembf(2,iembf)=ma0
         em(iembf)=elem(ij)
         ij=ij+1
         IF(iembf.EQ.lembf)THEN
            CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
            nelm=nelm+iembf
            iembf=0
         END IF
         ikk=ikk+ntlit
      END DO
c
      NINT=NTLIT*nelem
      nelm0=1
c
C-----------------------------------------------------------------------
c
c     Other diagonal elements in target CI case (skip otherwise)
c
C-----------------------------------------------------------------------
      IF(iexpc.EQ.0)THEN
         mst=2
         ma0=1
      ELSE
         mst=ncont+1
         ma00=notgt
         NFLG=0
         MFLG=0
      END IF
c
C-----------------------------------------------------------------------
c
c     Now make all the other diagonal elements
c
C-----------------------------------------------------------------------
      NFLG=0
      MFLG=0
      jdiag=idiag
      DO M=mst, NOCSF
         MC=ICDT(M)
         MD=INDT(M)
         IF(m.GT.ncont)jdiag=0
         NJK=0
         CALL ENRGMC(CJK,NORBB)
         CALL ENRGMZ(CH,xMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &               IPOSIT,MPOS,NODT(M),CDT(MC),NDT(MD),NODT(M),CDT(MC)
     &               ,NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,jdiag,NOEX)
         CALL CHECK(NH,THRES,CHT)
         CALL CHECK(NJK,THRES,CCJK)
         IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
            MALM=0
            WRITE(NFTW,9030)M
 9030       FORMAT(I6,' TH CSF IN ERROR')
            WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
            NALM=1
            RETURN
         END IF
C
         CALL STDIAG(CH,CHC,xMHC,CHR,xMHR,NORB,NORBL,NHC,NHR,THRES,
     &               IDIAG,mpos)
         NTLIC=NHC
C
         DO N=1, 4
            IF(IDIAG.EQ.0 .AND. NRG(N+1).GT.NRG(N))
     &         CALL DIFF(NRG(N),NRG(N+1),xMHR,CHR,CJK(1,N))
            CALL STORE(NTLIC,xMHC,CHC,NORB,THRES,CJK(1,N),N)
         END DO
C
         NJK=NTLIC-NHC
         nelm0=nelm0+1
         ma0=ma0+1
c
c     evaluate the integrals
c
         IF(LUSME.NE.0)WRITE(LUSME,*)ntliC, 1
         CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlic,xmhc,
     &               CHC,mxcint,xint1e,nint1e,xint2e,elem,symtyp,map,
     &               ncorb,mcorb,1,LUSME,enh_factor,enh_factor_type,
     &               ukrmolp_ints,NOEX)
         KCC=KCC+1
         IF(idiag.EQ.0)elem(1)=elem(1)+eone
         IF(LUSME.NE.0)WRITE(LUSME,*)MA0, MA0, 0
         IF(abs(elem(1)).GT.thrhm)THEN
            iembf=iembf+1
            nembf(1,iembf)=ma0
            nembf(2,iembf)=ma0
            em(iembf)=elem(1)
            IF(iembf.EQ.lembf)THEN
               CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
               nelm=nelm+iembf
               iembf=0
            END IF
         END IF
         NINT=NINT+NTLIC
         NINT0=NINT0+NTLIC
         NINT00=NINT00+NTLIC
      END DO
c
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c     The diagonal elements have now all been constructed
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C     Start the true OFF-DIAGONAL TERMS
C
c     The first portion of this is adapted from the scheme given by
c     L.A. Morgan & J. Tennyson, J Phys B, 26, 2429 (1993)
c     to expand prototype CSFs: numbering refers to appendix of
C
      NFLG=1
      MFLG=1
      CALL ENRGMC(CJK,NORBB)
      NJK=0
      MALM=0
C
c     loops 2000-2999 are for expanding prototype off-diagonal matrix
c     elements and/or for CI target contractions
c     skip if no neither to be done
      IF(iexpc.EQ.0)THEN
         nbeg=1
         ma0=0
         GO TO 4000
      ELSE
         WRITE(nftw,*)'This calculation is not allowed in the current co
     &de. Please contact Jonathan Tennyson and his group'
         STOP
      END IF
      mmm=1
      ma000=0
C-----------------------------------------------------------------------
c
c     EXPCSF / CI target contraction option: case (3)
c     compute remaining diagonal elements with continuum functions
c
C-----------------------------------------------------------------------
c
c     first: off-diagonal continuum case for present target state
c     (only one prototype)
      nci=1
      ncl=1
      mst=mmm-istep
      ntlit=0
      m=mst
c     then loop over CI componants
c     first do off-diagonal term within same target state
c     starting with those for the same target CSF
c
      m=m+istep
      MC=ICDT(M)
      MD=INDT(M)
      n=m+1
      NC=ICDT(N)
      ND=INDT(N)
      CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
      IF(NZ.GT.4)THEN
         WRITE(nftw,9050)m, n, nz
 9050    FORMAT(/' Surely some mistake here:'/' Matrix element (',i5,
     &          ',',i5,') has',i3,
     &          ' differences'/' but come from same CI target CSF: STOP'
     &          )
         STOP
      END IF
      NJK=0
      CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &            IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),NODT(M),CDT(MC),
     &            NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,idiag,NOEX)
C
      NHT=0
      IF(NFLG.NE.1)THEN
         NHT=NORBL
         CALL CHECK(NHT,THRES,CH)
      END IF
      NHC=0
      IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,THRES)
      IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
         NALM=1
         WRITE(NFTW,9040)N, M
 9040    FORMAT(I8,I7,5X,'FNS ARE NOT ORTHOGONAL')
         RETURN
      END IF
      NFLG=1
      ICH=NHC
      IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,THRES)
      IF(nhc.NE.0)THEN
c      IJK=NHC-ICH
         nint0=nint0+nhc
         ntlit=nhc
      END IF
      ma00=ma000
      nelm0=nelm0+1
      nint00=nint00+ntlit
c
c     (3) off diagonal case within main block
      nelem=notgt*(notgt-1)/2
      CALL EXPODL(kpt(m),kpt(n),notgt,nelem,ntlit,xmhc1,xmhc1,ntlit)
      ikk=1
c     evaluate the integrals
      DO i=1, notgt
         ma00=ma00+1
         na00=ma00
         DO j=i+1, notgt
            na00=na00+1
            IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
            CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlit,
     &               xmhc1(1,ikk),ch,mxcint,xint1e,nint1e,xint2e,
     &               elem,symtyp,map,ncorb,mcorb,nci,LUSME,enh_factor,
     &               enh_factor_type,ukrmolp_ints,NOEX)
            kcc=kcc+1
            ij=0
            ma0=ma00
            na0=na00
            ij=ij+1
            IF(LUSME.NE.0)WRITE(lusme,*)na00, ma00, 0
            IF(abs(elem(ij)).GT.thrhm)THEN
               iembf=iembf+1
               nembf(1,iembf)=na0
               nembf(2,iembf)=ma0
               em(iembf)=elem(ij)
               IF(iembf.EQ.lembf)THEN
                  CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                  NELM=NELM+iembf
                  iembf=0
               END IF
            END IF
            ikk=ikk+ntlit
         END DO
      END DO
      NINT=nint+NTLIT*nelem
C-----------------------------------------------------------------------
c
c      Case (2) Off diagonal element with L**2 CSF
c
C-----------------------------------------------------------------------
      na0=ncont2
      DO N=ncont+1, nocsf
         NC=ICDT(N)
         ND=INDT(N)
         na0=na0+1
c     then loop over CI componants
         m=mmm-istep
         ntlit=0
         m=m+istep
         MC=ICDT(M)
         MD=INDT(M)
         CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
         IF(NZ.GT.4)GO TO 2110
         NJK=0
         CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &               IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),NODT(M),CDT(MC)
     &               ,NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,0,NOEX)
C
         NHT=0
         IF(NFLG.NE.1)THEN
            NHT=NORBL
            CALL CHECK(NHT,THRES,CH)
         END IF
         NHC=0
         IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,THRES)
         IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
            NALM=1
            WRITE(NFTW,9040)N, M
            RETURN
         END IF
         NFLG=1
         ICH=NHC
         IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,THRES)
         IF(NHC.NE.0)THEN
            nint0=nint0+nhc
            ntlit=nhc
         END IF
 2110    CONTINUE
         nint00=nint00+ntlit
         nelm0=nelm0+1
         nelem=notgt
         CALL expcor(kpt(m),nelem,ntlit,xmhc1,xmhc1,ntlit)
         ikk=1
         ma00=ma000
c     evaluate the integrals
         DO i=1, nelem
            ma00=ma00+1
            ma0=ma00
            IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCL
            CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlit,
     &              xmhc1(1,ikk),ch,mxcint,xint1e,nint1e,xint2e,
     &              elem,symtyp,map,NCORB,mcorb,ncl,LUSME,enh_factor,
     &              enh_factor_type,ukrmolp_ints,NOEX)
            kcc=kcc+1
            IF(LUSME.NE.0)WRITE(lusme,*)na0, ma0, 0
            IF(abs(elem(1)).GT.thrhm)THEN
               iembf=iembf+1
               nembf(1,iembf)=na0
               nembf(2,iembf)=ma0
               em(iembf)=elem(1)
               IF(iembf.EQ.lembf)THEN
                  CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                  NELM=NELM+iembf
                  iembf=0
               END IF
            END IF
            ikk=ikk+ntlit
         END DO
         NINT=nint+NTLIT*nelem
      END DO
      ma000=ma000+notgt
      mmm=mmm+2
c
      ma0=ncont2
      nbeg=ncont+1
C-----------------------------------------------------------------------
C
C     Off-diagonal L**2 matrix elements: no expansion possible
C
C-----------------------------------------------------------------------
 4000 CONTINUE
      DO M=nbeg, NOCSF-1
         MC=ICDT(M)
         MD=INDT(M)
         ma0=ma0+1
         na0=ma0
         jdiag=idiag
         DO N=M+1, NOCSF
            NC=ICDT(N)
            ND=INDT(N)
            na0=na0+1
            CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
            IF(NZ.GT.4)CYCLE
            NJK=0
            IF(n.GT.ncont)jdiag=0
            CALL ENRGMZ(CH,xMHC,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &                  IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),NODT(M),
     &                  CDT(MC),NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,jdiag,
     &                  NOEX)
            IF(MALM.NE.0)THEN
               MALM=0
               WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
               NALM=1
               RETURN
            END IF
C
            NHT=0
            IF(NFLG.NE.1)THEN
               NHT=NORBL
               CALL CHECK(NHT,THRES,CH)
            END IF
            NHC=0
            IF(NH.GT.0)CALL MVDIAG(CHR,xMHC,CHR,xMHC,NH,NHC,THRES)
            IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
               NALM=1
               WRITE(NFTW,9040)N, M
               RETURN
            END IF
            NFLG=1
            ICH=NHC
            IF(NJK.GT.0)CALL MVDIAG(CHR,xMHC,CCJK,xMJK,NJK,NHC,THRES)
            IF(NHC.EQ.0)CYCLE
            nelm0=nelm0+1
c
            IF(LUSME.NE.0)WRITE(LUSME,*)NHC, 1
            CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,NHC,xMHC,
     &                  CHR,mxcint,xint1e,nint1e,xint2e,elem,symtyp,map,
     &                  NCORB,mcorb,1,LUSME,enh_factor,enh_factor_type,
     &                  ukrmolp_ints,NOEX)
            kcc=kcc+1
            IF(LUSME.NE.0)WRITE(lusme,*)na0, ma0, 0
            IF(abs(elem(1)).GT.thrhm)THEN
               iembf=iembf+1
               nembf(1,iembf)=na0
               nembf(2,iembf)=ma0
               em(iembf)=elem(1)
               IF(iembf.EQ.lembf)THEN
                  CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                  NELM=NELM+iembf
                  iembf=0
               END IF
            END IF
            NINT=NINT+NHC
            NINT0=NINT0+NHC
            NINT00=NINT00+NHC
         END DO
      END DO
      IF(LUSME.NE.0)THEN
         WRITE(lusme,*)headerp
         WRITE(lusme,*)KCC-1
      END IF
c     write out last buffer of the Hamiltonian matrix
      IF(iembf.NE.0)THEN
         CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
         NELM=NELM+iembf
         iembf=0
      END IF
c     write dummy record at end of file
      CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
c
      N=(NOCSF*(NOCSF+1))/2
      m=(mOCSF*(mOCSF+1))/2
      WRITE(NFTW,6000)m, NELM
 6000 FORMAT(//' TOTAL H Matrix     ELEMENTS      =',
     &       I10/' NON-ZERO ELEMENTS evaluated      =',I10)
      IF(iexpc.NE.0)WRITE(nftw,6010)n, nelm0
 6010 FORMAT(' Total    prototype ELEMENTS      =',
     &       I10/' Non-zero prototype ELEMENTS      =',I10)
      IF(iexpc.NE.0)WRITE(nftw,6020)nint0, nint00
 6020 FORMAT(/' Number (prototype) integrals     =',
     &       I10/' Compressed (prototype) integrals =',I10)
      WRITE(nftw,6030)nint
 6030 FORMAT(' Number integrals evaluated       =',I10)
c
      DEALLOCATE(xmjk,ccjk,xmhr,chr,chc,xmhc,xmhc1,elem)
      RETURN
      END SUBROUTINE ENRGMS
!*==enrgmt.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE ENRGMT(NORB,MN,MM,MS,NELT,NDTRF,MDTR,NOCSF,NODT,INDT,
     &                  NDT,ICDT,CDT,THRES,NALM,NFTW,IPOSIT,MAP,MPOS,
     &                  mcon,mcorb,NSRB,SYMTYP,idiag,nelm,nfte,lembf,
     &                  xint1e,nint1e,xint2e,thrhm,IODR,NRI,NSM,NOB,
     &                  MBAS,NBAS,NBASH,IPAIR,istart,jump,ukrmolp_ints,
     &                  NOEX)
c
c    ENRGMT controls the construction and evaluation of matrix elements
c    for target wavefunctions that need to be phase corrected
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: MALM, MFLG, NFLG, NH, NJK, NJKMX
      COMMON /ENA   / NH, NFLG, MFLG
      COMMON /STPQ  / NJK, MALM, NJKMX
C
C Dummy arguments
C
      INTEGER :: IDIAG, IODR, IPOSIT, ISTART, JUMP, LEMBF, NALM, NELM, 
     &           NELT, NFTE, NFTW, NINT1E, NOCSF, NORB, NSRB, SYMTYP,
     &           NOEX
      REAL(KIND=wp) :: THRES, THRHM
      REAL(KIND=wp), DIMENSION(*) :: CDT
      REAL(KIND=wp), ALLOCATABLE :: XINT1E(:), XINT2E(:)
      INTEGER, DIMENSION(*) :: ICDT, INDT, MAP, MBAS, MCON, MCORB, MDTR, 
     &                         MM, MN, MPOS, MS, NBAS, NBASH, NDT, 
     &                         NDTRF, NOB, NODT, NRI, NSM
      INTEGER, DIMENSION(0:*) :: IPAIR
      LOGICAL :: ukrmolp_ints
      INTENT (IN) ICDT, INDT, ISTART, JUMP, NFTW, NOCSF, NSRB, THRHM
      INTENT (IN) ukrmolp_ints
      INTENT (OUT) NALM
      INTENT (INOUT) NELM
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CCJK, CHC, CHR
      REAL(KIND=wp), DIMENSION(norb) :: CHT
      REAL(KIND=wp), DIMENSION(norb*(norb+1)/2,4) :: CJK
      REAL(KIND=wp), DIMENSION(1) :: ELEM
      REAL(KIND=wp), DIMENSION(lembf) :: EM
      REAL(KIND=wp) :: EONE
      INTEGER :: I, IEMBF, M, MA0, MC, MD, MXCIN, N, NA0, NC, 
     &           NCORB, ND, NHC, NHR, NHT, NORBB, NORBL, NTLIC, NTLIR, 
     &           NZ
      INTEGER, DIMENSION(nsrb) :: NBB
      INTEGER, DIMENSION(2,lembf) :: NEMBF
      INTEGER, DIMENSION(5) :: NRG
      INTEGER(longint), ALLOCATABLE, DIMENSION(:,:) :: XMH, XMHC, XMHR, 
     &                XMJK
C
C*** End of declarations rewritten by SPAG
C
      DATA NA0/0/
      EQUIVALENCE(NRG(1),NHR)
C
      ncorb=0
      NORBB=(NORB*(NORB+1))/2
      NORBL=NORB
      njkmx=norbb
      NALM=0
C-----------------------------------------------------------------------
      mxcin=4*NORBB+NORBL
      ALLOCATE(xmhr(2,mxcin),xmh(2,norb),xmhc(2,mxcin),xmjk(2,norbb),
     &         chr(mxcin),chc(mxcin),ccjk(norbb))
C-----------------------------------------------------------------------
      NFLG=0
      MFLG=0
C     ZERO ARRAY FOR PZERO
      DO I=1, NSRB
         NBB(I)=0
      END DO
C-----------------------------------------------------------------------
C
c     First construct the (1,1) diagonal element
c
C-----------------------------------------------------------------------
      nelm=0
      m=istart
      MC=ICDT(M)
      MD=INDT(M)
      NJK=0
      MALM=0
      CALL ENRGMC(CJK,NORBB)
      CALL ENRGMZ(CHc,xMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &            IPOSIT,MPOS,NODT(M),CDT(MC),NDT(MD),NODT(M),CDT(MC),
     &            NDT(MD),CJK,xmJK,CCJK,SYMTYP,mcon,idiag,NOEX)
      CALL CHECK(NH,THRES,CHT)
      CALL CHECK(NJK,THRES,CCJK)
      IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
         MALM=0
         WRITE(NFTW,9010)
 9010    FORMAT('    THE 1ST CSF IS IN ERROR')
         WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
 9020    FORMAT('    NH =',I3,' NJK =',I3,' MALM =',I3,' NJKMX =',I7)
         NALM=1
         RETURN
      END IF
C
      CALL STDIAG(CHc,CHR,xMHR,CHR,xMHR,NORB,NORBL,NHR,NHC,THRES,1,mpos)
      NTLIR=NHR
C
      DO N=1, 4
         CALL STORE(NTLIR,xMHR,CHR,NORB,THRES,CJK(1,N),N)
         NRG(N+1)=NTLIR
      END DO
C
      NJK=NTLIR-NHR
      CALL TPINDE(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlir,xmhr,chr,
     &            ntlir,xint1e,nint1e,xint2e,elem,symtyp,map,ncorb,
     &            mcorb,1,ukrmolp_ints)
      IF(idiag.EQ.0)eone=elem(1)
      nembf(1,1)=1
      nembf(2,1)=1
      em(1)=elem(1)
      iembf=1
C
      IF(NOCSF.LE.1)GO TO 5000
C-----------------------------------------------------------------------
c
c     Now make the other diagonal elements
c
C-----------------------------------------------------------------------
      DO MA0=2, NOCSF
         m=m+jump
         MC=ICDT(M)
         MD=INDT(M)
         NJK=0
         CALL ENRGMC(CJK,NORBB)
         CALL ENRGMZ(CHc,XMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &               IPOSIT,MPOS,NODT(M),CDT(MC),NDT(MD),NODT(M),CDT(MC)
     &               ,NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,idiag,NOEX)
         CALL CHECK(NH,THRES,CHT)
         CALL CHECK(NJK,THRES,CCJK)
         IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
            MALM=0
            WRITE(NFTW,9030)M
 9030       FORMAT(I6,' TH CSF IN ERROR')
            WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
            NALM=1
            RETURN
         END IF
C
         CALL STDIAG(CHc,CHC,xMHC,CHR,xMHR,NORB,NORBL,NHC,NHR,THRES,
     &               IDIAG,mpos)
         NTLIC=NHC
C
         DO N=1, 4
            IF(IDIAG.EQ.0 .AND. NRG(N+1).GT.NRG(N))
     &         CALL DIFF(NRG(N),NRG(N+1),xMHR,CHR,CJK(1,N))
            CALL STORE(NTLIC,xMHC,CHC,NORB,THRES,CJK(1,n),N)
         END DO
C
         NJK=NTLIC-NHC
c
c     evaluate the integrals
c
         CALL TPINDE(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlic,xmhc,
     &               CHC,ntlic,xint1e,nint1e,xint2e,elem,symtyp,map,
     &               ncorb,mcorb,1,ukrmolp_ints)
         IF(idiag.EQ.0)elem(1)=elem(1)+eone
         IF(abs(elem(1)).GT.thrhm)THEN
            iembf=iembf+1
            nembf(1,iembf)=ma0
            nembf(2,iembf)=ma0
            em(iembf)=elem(1)
            IF(iembf.EQ.lembf)THEN
               CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
               nelm=nelm+iembf
               iembf=0
            END IF
         END IF
      END DO
C-----------------------------------------------------------------------
C
C     Off-diagonal matrix elements
C
C-----------------------------------------------------------------------
      NFLG=1
      MFLG=1
      CALL ENRGMC(CJK,NORBB)
      NJK=0
      MALM=0
C     ZERO ARRAY FOR PZERO
      DO I=1, NSRB
         NBB(I)=0
      END DO
      m=istart-jump
      DO MA0=1, NOCSF-1
         m=m+jump
         MC=ICDT(M)
         MD=INDT(M)
         n=m
         DO NA0=MA0+1, NOCSF
            n=n+jump
            NC=ICDT(N)
            ND=INDT(N)
            CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
            IF(NZ.GT.4)CYCLE
            NJK=0
            CALL ENRGMZ(CHc,xMHC,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,MS,
     &                  NORBL,IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),
     &                  NODT(M),CDT(MC),NDT(MD),CJK,xMJK,CCJK,SYMTYP,
     &                  mcon,idiag,NOEX)
C
            NHT=0
            IF(NFLG.NE.1)THEN
               NHT=NORBL
               CALL CHECK(NHT,THRES,CHc)
            END IF
            NHC=0
            IF(NH.GT.0)CALL MVDIAG(CHR,xMHC,CHR,xMHC,NH,NHC,THRES)
            IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
               NALM=1
               WRITE(NFTW,9040)N, M
 9040          FORMAT(I8,I7,5X,'FNS ARE NOT ORTHOGONAL')
               RETURN
            END IF
            NFLG=1
            IF(NJK.GT.0)CALL MVDIAG(CHR,xMHC,CCJK,xMJK,NJK,NHC,THRES)
            IF(NHC.EQ.0)CYCLE
c
            CALL TPINDE(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,NHC,xMHC,
     &                  CHR,nhc,xint1e,nint1e,xint2e,elem,symtyp,map,
     &                  NCORB,mcorb,1,ukrmolp_ints)
            IF(abs(elem(1)).GT.thrhm)THEN
               iembf=iembf+1
               nembf(1,iembf)=na0
               nembf(2,iembf)=ma0
               em(iembf)=elem(1)
               IF(iembf.EQ.lembf)THEN
                  CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                  nelm=nelm+iembf
                  iembf=0
               END IF
            END IF
         END DO
      END DO
c
c     write out last buffer of the Hamiltonian matrix
 5000 IF(iembf.NE.0)THEN
         CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
         nelm=nelm+iembf
         iembf=0
      END IF
c     write dummy record at end of file
      CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
c      WRITE(NFTW,6000) (NOCSF*(NOCSF+1))/2,NELM
c 6000 FORMAT(//' TOTAL H Matrix     ELEMENTS      =',I10
c     *      /' NON-ZERO ELEMENTS evaluated      =',I10)
C
      DEALLOCATE(xmhr,xmh,xmhc,xmjk,chr,chc,ccjk)
c
      RETURN
      END SUBROUTINE ENRGMT
!*==enrgmx.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE ENRGMX(NORB,MN,MM,MS,NELT,NDTRF,MDTR,NOCSF,NODT,INDT,
     &                  NDT,ICDT,CDT,THRES,NALM,IDIAG,nelm,NFTW,IPOSIT,
     &                  MAP,MPOS,mcon,mcorb,NSRB,SYMTYP,ncont,ncont2,
     &                  nfte,lembf,xint1e,nint1e,xint2e,thrhm,IODR,NRI,
     &                  NSM,NOB,MBAS,NBAS,NBASH,IPAIR,NCORB,kpt,mocsf,
     &                  iexpc,icitg,ctgt,eig,numtgt,ntgsym,nctgt,notgt,
     &                  nummx,gucont,mcont,lusme,enh_factor,
     &                  enh_factor_type,ukrmolp_ints,NOEX)
c
c    ENRGMX controls the construction and evaluation of matrix elements
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: MALM, MFLG, NFLG, NH, NJK, NJKMX
      COMMON /ENA   / NH, NFLG, MFLG
      COMMON /STPQ  / NJK, MALM, NJKMX
C
C Dummy arguments
C
      INTEGER :: ICITG, IDIAG, IEXPC, IODR, IPOSIT, LEMBF, LUSME, MOCSF, 
     &           NALM, NCONT, NCONT2, NCORB, NELM, NELT, NFTE, NFTW, 
     &           NINT1E, NOCSF, NORB, NSRB, NTGSYM, NUMMX, SYMTYP,
     &           enh_factor_type,NOEX
      REAL(KIND=wp) :: THRES, THRHM, enh_factor
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CDT, CTGT, EIG
      REAL(KIND=wp), ALLOCATABLE :: XINT1E(:), XINT2E(:)
      INTEGER, DIMENSION(ntgsym) :: GUCONT, MCONT, NCTGT, NOTGT, NUMTGT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ICDT, INDT, MAP, MCON,
     &                                      MCORB, MDTR, MM, MN, MPOS,
     &                                      MS, NDT, NODT, NRI, NSM
      INTEGER, DIMENSION(*) :: NOB, NDTRF, MBAS, NBAS, NBASH
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KPT
      LOGICAL :: ukrmolp_ints
      INTENT (IN) EIG, GUCONT, ICDT, ICITG, IEXPC, INDT, MCONT, MOCSF, 
     &            NCONT, NCONT2, NFTW, NOCSF, NOTGT, NSRB, NTGSYM, 
     &            NUMMX, THRHM, enh_factor,enh_factor_type, ukrmolp_ints
      INTENT (OUT) NALM
      INTENT (INOUT) NELM
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CCJK, CHC, CHR, ELEM,
     &                                            CH, CHT, CDEF, EM,
     &                                            COEF
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: CJK
      REAL(KIND=wp) :: EONE
      CHARACTER(LEN=32) :: HEADERP
      INTEGER :: I, ICH, ICTGT, ICTGT1, ICTGT2, ICTGTP, IEMBF, IERR, II, 
     &           IJ, IKK, IOTGT, IOTGT1, IOTGT2, ISTEP, ISTEP1, ISTEP2, 
     &           ITARG, J, JCTGT, JCTGTP, JDIAG, JJ, KCC, M, MA0, MA00, 
     &           MA000, MC, MD, MMM, MST, MXCBIG, MXCIN, MXCINT, MXCORB, 
     &           N, NA0, NA00, NA000, NBEG, NC, NCHAN, NCHANA, NCHANB, 
     &           NCI, NCL, ND, NEED, NELEM, NELM0, NELMX, NHC, NHR, NHT, 
     &           NINT, NINT0, NINT00, NNN, NORBB, NORBL, NOTOP, NOTOP1, 
     &           NOTOP2, NTGT1, NTGT2, NTGTS, NTLIC, NTLIR, NTLIT, NZ
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBB
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NEMBF
      INTEGER, DIMENSION(5) :: NRG
      INTEGER(longint), ALLOCATABLE, DIMENSION(:,:) :: XMH
      INTEGER(longint), ALLOCATABLE, DIMENSION(:,:) :: XMHC, XMHC1, 
     &                XMHR, XMJK
C
C*** End of declarations rewritten by SPAG
C
      EQUIVALENCE(NRG(1),NHR)
      DATA NA0/0/, headerp/'************************NUMOFP  '/
C
      ALLOCATE(CH(norb*(norb+1)/2), CHT(norb), CJK(norb*(norb+1)/2,4))
      ALLOCATE(COEF(nummx*nummx), EM(lembf), NBB(nsrb), NEMBF(2,lembf))
      ALLOCATE(XMH(2,norb))
C
      ierr=0
      NORBB=(NORB*(NORB+1))/2
      NORBL=NORB
      NALM=0
C-----------------------------------------------------------------------
      mxcin=4*NORBB+NORBL
      mxcint=2*mxcin
      mxcbig=mxcint*nummx**2
c
c *** xmhc needs to be bigger than mxcint in EXPODL
c     the following seems to be a gross over estimate
      mxcorb=0
      DO i=1, ntgsym
         mxcorb=max(mxcorb,norbb*notgt(i)*(notgt(i)-1)/2)
      END DO
      mxcorb=100000000
      WRITE(6,'(5X,"MXCORB SET TO: ",i15)') mxcorb
c
      IF(icitg.EQ.0)THEN
         ALLOCATE(chc(mxcin),xmhc(2,mxcint),xmhc1(2,mxcint),stat=ierr)
         need=mxcin+2*mxcint
      ELSE
         ALLOCATE(chc(mxcbig),xmhc(2,mxcorb),xmhc1(2,mxcint),stat=ierr)
         need=mxcbig+mxcorb+mxcint
      END IF
      IF(ierr.NE.0)THEN
         WRITE(6,99)need
 99      FORMAT(/' UNABLE TO ALLOCATE MEMORY IN ENRGMX, need extra',i10,
     &          ' (real) words')
         STOP
      END IF
C
C    This is a best guess, although errors are flagged in lower level
C    routines, they are not always trapped in this routine
      NJKMX=mxcint
 
c*****************************************************************************
c The definition of nelmx has been modified in the following statment and in
c loop 1 below. This modification takes into account the off diagonal case
c with continuum and different target states (see loop 3000 in ENRGMX). Read
c note below for details.(N. Vinci, 20/11/01)
c
c The matrix elements are stored in the array elem that is allocated in
c ENRGMX (and also in ENRGMS) as elem(nelmx). Thus, its dimension, nelmx,
c must be consistent with the value of variable nci according to the dimension
c statment for elem(nci) in subroutine PINDEX. Nelmx must be big enough to
c provide as many as nci elements to be stored in elem-dynamical array
c dimension inconsistence causing execution problems for same input values
c of numtgt (number of target states per symmetry).
c****************************************************************************
 
      nelmx=numtgt(1)**2
 
      DO i=2, ntgsym
         nelmx=max(nelmx,numtgt(i)**2)
      END DO
      ALLOCATE(elem(nelmx),xmjk(2,njkmx),ccjk(njkmx),xmhr(2,mxcin),
     &         chr(mxcin))
C-----------------------------------------------------------------------
      NFLG=0
      MFLG=0
C     ZERO ARRAY FOR PZERO
      DO I=1, NSRB
         NBB(I)=0
      END DO
C-----------------------------------------------------------------------
C
c     First construct the (1,1) diagonal element
c
C-----------------------------------------------------------------------
      NJK=0
      MALM=0
      nelm=0
      nelem=1
      ma0=0
      notop=0
      m=1
      iembf=0
      IF(iexpc.EQ.0)THEN
         istep=1
         istep1=1
         istep2=1
         nchan=1
         nchana=1
         nchanb=1
      ELSE
         istep=2
         istep1=2
         istep2=2
      END IF
      CALL ENRGMC(CJK,NORBB)
      CALL ENRGMZ(CH,xMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,IPOSIT,
     &            MPOS,NODT(1),CDT(1),NDT(1),NODT(1),CDT(1),NDT(1),CJK,
     &            xMJK,CCJK,SYMTYP,mcon,idiag,NOEX)
      CALL CHECK(NH,THRES,CHT)
      CALL CHECK(NJK,THRES,CCJK)
      IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
         MALM=0
         WRITE(NFTW,9010)
 9010    FORMAT('    THE 1ST CSF IS IN ERROR')
         WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
 9020    FORMAT('    NH =',I3,' NJK =',I3,' MALM =',I3,' NJKMX =',I7)
         NALM=1
         RETURN
      END IF
C
      CALL STDIAG(CH,CHR,xMHR,CHR,xMHR,NORB,NORBL,NHR,NHC,THRES,1,mpos)
      NTLIR=NHR
C
      DO N=1, 4
         CALL STORE(NTLIR,xMHR,CHR,NORB,THRES,CJK(1,n),N)
         NRG(N+1)=NTLIR
      END DO
C
      NJK=NTLIR-NHR
      nint0=ntlir
c
c     CI target contraction option:
c     compute other diagonal elements, add and compress
c
      KCC=1
      IF(icitg.NE.0)THEN
         nci=numtgt(1)*(numtgt(1)+1)/2
         ntlit=0
         IF(idiag.EQ.0)THEN
c        IDIAG=0 option treats first matrix element as a special case
            IF(LUSME.NE.0)WRITE(LUSME,*)ntlir, 1
            CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlir,
     &                  xmhr,CHR,mxcint,xint1e,nint1e,xint2e,elem,
     &                  symtyp,map,ncorb,mcorb,1,lusme,enh_factor,
     &                  enh_factor_type,ukrmolp_ints,NOEX)
            eone=elem(1)
            KCC=KCC+1
            IF(LUSME.NE.0)WRITE(LUSME,*)1, 1, 0
         ELSE
c        first: account for first term in CI target expansion
            CALL mkcfd(coef,ctgt(1),numtgt(1),nctgt(1))
            CALL Comprs(nci,ntlit,ntlir,xMHC,CHC,xMHR,CHR,coef,mxcint)
         END IF
c     then loop over the other CSFs for first target state
         DO ictgt=2, nctgt(1)
            m=m+istep
            MC=ICDT(M)
            MD=INDT(M)
            njk=0
            NFLG=0
            MFLG=0
C
            CALL ENRGMC(CJK,NORBB)
            CALL ENRGMZ(CH,xMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &                  IPOSIT,MPOS,NODT(M),CDT(MC),NDT(MD),NODT(M),
     &                  CDT(MC),NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,idiag,
     &                  NOEX)
            CALL CHECK(NH,THRES,CHT)
            CALL CHECK(NJK,THRES,CCJK)
            IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
               MALM=0
               WRITE(NFTW,9030)m
               WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
               NALM=1
               RETURN
            END IF
C
            CALL STDIAG(CH,CH,xMHC1,CHR,xMHR,NORB,NORBL,NHC,NHR,THRES,
     &                  idiag,mpos)
            NTLIC=NHC
C
            DO N=1, 4
               IF(IDIAG.EQ.0 .AND. NRG(N+1).GT.NRG(N))
     &            CALL DIFF(NRG(N),NRG(N+1),xMHR,CHR,CJK(1,N))
               CALL STORE(NTLIC,xMHC1,CH,NORB,THRES,CJK(1,n),N)
            END DO
C
            NJK=NTLIC-NHC
            nint0=nint0+ntlic
            CALL mkcfd(coef,ctgt(ictgt),numtgt(1),nctgt(1))
            CALL Comprs(nci,ntlit,ntlic,xMHC,CHC,xMHC1,CH,coef,mxcint)
C
C        OFF-DIAGONAL TERMS
C
            NFLG=1
            MFLG=1
            CALL ENRGMC(CJK,NORBB)
            NJK=0
            MALM=0
            n=1-istep
            DO jctgt=1, ictgt-1
               n=n+istep
               NC=ICDT(N)
               ND=INDT(N)
               CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
               IF(NZ.GT.4)CYCLE
               NJK=0
               CALL ENRGMZ(CH,xMHC1,CH,NDTRF,MDTR,NELT,NORB,MN,MM,MS,
     &                     NORBL,IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),
     &                     NODT(M),CDT(MC),NDT(MD),CJK,xMJK,CCJK,SYMTYP,
     &                     mcon,idiag,NOEX)
C
               NHT=0
               IF(NFLG.NE.1)THEN
                  NHT=NORBL
                  CALL CHECK(NHT,THRES,CH)
               END IF
               NHC=0
               IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CH,xMHC1,NH,NHC,THRES)
               IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                  NALM=1
                  WRITE(NFTW,9040)N, M
                  RETURN
               END IF
               NFLG=1
               ICH=NHC
               IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,THRES)
               IF(NHC.EQ.0)CYCLE
c         IJK=NHC-ICH
               nint0=nint0+nhc
               CALL mkcfo(coef,ctgt(ictgt),ctgt(jctgt),numtgt(1),
     &                    nctgt(1))
               CALL Comprs(nci,ntlit,nhc,xMHC,CHC,xMHC1,CH,coef,mxcint)
            END DO
         END DO
c
         IF(iexpc.NE.0)THEN
c
c        expand prototype diagonal CSF
c
            nelem=notgt(1)
            CALL expdg(kpt(m),nelem,ntlit,xmhc,xmhc,ntlit)
            nint00=ntlit
c
c     evaluate the integrals
c
            ikk=1
            iembf=0
            DO i=1, nelem
               ma0=i-notgt(1)
               IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
               IF (ikk+ntlit > mxcorb) then
                  write(*,6040) ikk+ntlit
                  stop
               endif
               CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlit,
     &               xmhc(1,ikk),chc,mxcint,xint1e,nint1e,xint2e,
     &               elem,symtyp,map,ncorb,mcorb,nci,lusme,enh_factor,
     &               enh_factor_type,ukrmolp_ints,NOEX)
               KCC=KCC+1
               IF(icitg.EQ.0 .AND. i.EQ.1 .AND. idiag.EQ.0)THEN
                  eone=elem(1)
                  elem(1)=elem(1)-eone
               END IF
               ij=1
               DO ii=1, numtgt(1)
                  IF(idiag.EQ.0)THEN
                     elem(ij)=elem(ij)+eone
                  ELSE IF(idiag.EQ.2)THEN
                     elem(ij)=elem(ij)+eig(ii)
                  END IF
                  ma0=ma0+notgt(1)
                  DO jj=ii, numtgt(1)
                     iembf=iembf+1
                     na0=ma0+(jj-ii)*notgt(1)
                     nembf(1,iembf)=na0
                     nembf(2,iembf)=ma0
                     IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                     em(iembf)=elem(ij)
                     ij=ij+1
                     IF(iembf.EQ.lembf)THEN
                        CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                        nelm=nelm+iembf
                        iembf=0
                     END IF
                  END DO
               END DO
               ikk=ikk+ntlit
            END DO
         END IF
      ELSE
         ntlit=ntlir
         nci=1
         ncl=1
         IF(iexpc.NE.0)THEN
c
c        expand prototype diagonal CSF
c
            nelem=notgt(1)
            CALL expdg(kpt(m),nelem,ntlit,xmhr,xmhr,ntlit)
         END IF
         nint00=ntlit
c
c     evaluate the integrals
c
         ikk=1
         iembf=0
         DO i=1, nelem
            ma0=i-notgt(1)
            IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
            CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlit,
     &              xmhr(1,ikk),Chr,mxcint,xint1e,nint1e,xint2e,
     &              elem,symtyp,map,ncorb,mcorb,nci,lusme,enh_factor,
     &              enh_factor_type,ukrmolp_ints,NOEX)
            kcc=kcc+1
            IF(icitg.EQ.0 .AND. i.EQ.1 .AND. idiag.EQ.0)THEN
               eone=elem(1)
               elem(1)=elem(1)-eone
            END IF
            ij=1
            DO ii=1, numtgt(1)
               IF(idiag.EQ.0)THEN
                  elem(ij)=elem(ij)+eone
               ELSE IF(idiag.EQ.2)THEN
                  elem(ij)=elem(ij)+eig(ii)
               END IF
               ma0=ma0+notgt(1)
               DO jj=ii, numtgt(1)
                  iembf=iembf+1
                  na0=ma0+(jj-ii)*notgt(1)
                  nembf(1,iembf)=na0
                  nembf(2,iembf)=ma0
                  IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                  em(iembf)=elem(ij)
                  ij=ij+1
                  IF(iembf.EQ.lembf)THEN
                     CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                     nelm=nelm+iembf
                     iembf=0
                  END IF
               END DO
            END DO
            ikk=ikk+ntlit
         END DO
      END IF
c
      NINT=NTLIT*nelem
      nelm0=1
C
      IF(NOCSF.LE.numtgt(1))GO TO 5000
C-----------------------------------------------------------------------
c
c     Other diagonal elements in target CI case (skip otherwise)
c
C-----------------------------------------------------------------------
      IF(icitg.EQ.0 .AND. iexpc.EQ.0)THEN
         mst=2
         ma0=1
         GO TO 1500
      END IF
      mst=ncont+1
      ictgtp=nctgt(1)*numtgt(1)
      itarg=numtgt(1)
      ma00=notgt(1)*numtgt(1)
c
c     EXPCSF / CI target contraction option:
c     compute remaining diagonal elements with continuum functions
c
      DO ntgts=2, ntgsym
         IF(icitg.NE.0)THEN
            nci=numtgt(ntgts)*(numtgt(ntgts)+1)/2
         ELSE
            nci=1
         END IF
         IF(iexpc.EQ.0)THEN
            notop=notgt(ntgts)-1
            istep=notgt(ntgts)
            m=m+1-istep
         END IF
         ntlit=0
         DO iotgt=0, notop
            IF(iotgt.GT.0)m=m-notgt(ntgts)*nctgt(ntgts)+1-istep
            nnn=m
c     then loop over CI componants (skip out if icitg=0)
            DO ictgt=1, nctgt(ntgts)
               m=m+istep
               MC=ICDT(M)
               MD=INDT(M)
               njk=0
               NFLG=0
               MFLG=0
C
               CALL ENRGMC(CJK,NORBB)
               CALL ENRGMZ(CH,xMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,
     &                     NORBL,IPOSIT,MPOS,NODT(M),CDT(MC),NDT(MD),
     &                     NODT(M),CDT(MC),NDT(MD),CJK,xMJK,CCJK,SYMTYP,
     &                     mcon,idiag,NOEX)
               CALL CHECK(NH,THRES,CHT)
               CALL CHECK(NJK,THRES,CCJK)
               IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
                  MALM=0
                  WRITE(NFTW,9030)m
                  WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
                  NALM=1
                  RETURN
               END IF
C
               CALL STDIAG(CH,CH,xMHC1,CHR,xMHR,NORB,NORBL,NHC,NHR,
     &                     THRES,idiag,mpos)
               NTLIC=NHC
C
               DO N=1, 4
                  IF(IDIAG.EQ.0 .AND. NRG(N+1).GT.NRG(N))
     &               CALL DIFF(NRG(N),NRG(N+1),xMHR,CHR,CJK(1,N))
                  CALL STORE(NTLIC,xMHC1,CH,NORB,THRES,CJK(1,n),N)
               END DO
C
               NJK=NTLIC-NHC
               nint0=nint0+ntlic
               IF(icitg.EQ.0)THEN
                  ntlit=ntlic
                  EXIT
               END IF
               CALL mkcfd(coef,ctgt(ictgtp+ictgt),numtgt(ntgts),
     &                    nctgt(ntgts))
               CALL Comprs(nci,ntlit,ntlic,xMHC,CHC,xMHC1,CH,coef,
     &                     mxcint)
C
C     OFF-DIAGONAL TERMS
C
               NFLG=1
               MFLG=1
               CALL ENRGMC(CJK,NORBB)
               NJK=0
               MALM=0
               n=nnn
               DO jctgt=1, ictgt-1
                  n=n+istep
                  NC=ICDT(N)
                  ND=INDT(N)
                  CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
                  IF(NZ.GT.4)CYCLE
                  NJK=0
                  CALL ENRGMZ(CH,xMHC1,CH,NDTRF,MDTR,NELT,NORB,MN,MM,MS,
     &                        NORBL,IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),
     &                        NODT(M),CDT(MC),NDT(MD),CJK,xMJK,CCJK,
     &                        SYMTYP,mcon,idiag,NOEX)
C
                  NHT=0
                  IF(NFLG.NE.1)THEN
                     NHT=NORBL
                     CALL CHECK(NHT,THRES,CH)
                  END IF
                  NHC=0
                  IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CH,xMHC1,NH,NHC,THRES)
                  IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                     NALM=1
                     WRITE(NFTW,9040)N, M
                     RETURN
                  END IF
                  NFLG=1
                  ICH=NHC
                  IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,
     &                                    THRES)
                  IF(NHC.EQ.0)CYCLE
c      IJK=NHC-ICH
                  nint0=nint0+nhc
                  CALL mkcfo(coef,ctgt(ictgtp+ictgt),ctgt(ictgtp+jctgt),
     &                       numtgt(ntgts),nctgt(ntgts))
                  CALL Comprs(nci,ntlit,nhc,xMHC,CHC,xMHC1,CH,coef,
     &                        mxcint)
               END DO
            END DO
c
            nint00=nint00+ntlit
            IF(iexpc.NE.0)THEN
c
c        expand prototype diagonal CSF
c
               nelem=notgt(ntgts)
               CALL expdg(kpt(m),nelem,ntlit,xmhc,xmhc,ntlit)
            END IF
c
c     evaluate the integrals
c
            nelm0=nelm0+1
            ikk=1
            DO i=1, nelem
               ma0=i-notgt(ntgts)+ma00+iotgt
               IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
               IF (ikk+ntlit > mxcorb) then
                  write(*,6040) ikk+ntlit
                  stop
               endif
               CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlit,
     &               xmhc(1,ikk),CHC,mxcint,xint1e,nint1e,xint2e,
     &               elem,symtyp,map,ncorb,mcorb,nci,lusme,enh_factor,
     &               enh_factor_type,ukrmolp_ints,NOEX)
               kcc=kcc+1
               ij=1
               DO ii=1, numtgt(ntgts)
                  IF(idiag.EQ.0)THEN
                     elem(ij)=elem(ij)+eone
                  ELSE IF(idiag.EQ.2)THEN
                     elem(ij)=elem(ij)+eig(itarg+ii)
                  END IF
                  ma0=ma0+notgt(ntgts)
                  DO jj=ii, numtgt(ntgts)
                     iembf=iembf+1
                     na0=ma0+(jj-ii)*notgt(ntgts)
                     nembf(1,iembf)=na0
                     nembf(2,iembf)=ma0
                     IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                     em(iembf)=elem(ij)
                     ij=ij+1
                     IF(iembf.EQ.lembf)THEN
                        CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                        nelm=nelm+iembf
                        iembf=0
                     END IF
                  END DO
               END DO
               ikk=ikk+ntlit
            END DO
c
            NINT=nint+NTLIT*nelem
         END DO
         itarg=itarg+numtgt(ntgts)
         ma00=ma00+notgt(ntgts)*numtgt(ntgts)
         ictgtp=ictgtp+nctgt(ntgts)*numtgt(ntgts)
      END DO
c
C-----------------------------------------------------------------------
c
c     Now make all the other diagonal elements
c
C-----------------------------------------------------------------------
      NFLG=0
      MFLG=0
 1500 jdiag=idiag
      DO M=mst, NOCSF
         MC=ICDT(M)
         MD=INDT(M)
         IF(m.GT.ncont)jdiag=0
         NJK=0
         CALL ENRGMC(CJK,NORBB)
c
c*************************************************************************
c The variable mcon (dimension NSRB) in SUBROUTINE ENRGMZ is set equal to
c mcon in the call statment below. This was set to mcorb (dimension NORB)
c in a previous version of scatci. (N. Vinci, 10/11/01).
c*************************************************************************
c
         CALL ENRGMZ(CH,xMH,CHT,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &               IPOSIT,MPOS,NODT(M),CDT(MC),NDT(MD),NODT(M),CDT(MC)
     &               ,NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,jdiag,NOEX)
         CALL CHECK(NH,THRES,CHT)
         CALL CHECK(NJK,THRES,CCJK)
         IF(NH.NE.0 .OR. NJK.NE.0 .OR. MALM.NE.0)THEN
            MALM=0
            WRITE(NFTW,9030)M
 9030       FORMAT(I6,' TH CSF IN ERROR')
            WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
            NALM=1
            RETURN
         END IF
C
         CALL STDIAG(CH,CHC,xMHC,CHR,xMHR,NORB,NORBL,NHC,NHR,THRES,
     &               IDIAG,mpos)
         NTLIC=NHC
C
         DO N=1, 4
            IF(IDIAG.EQ.0 .AND. NRG(N+1).GT.NRG(N))
     &         CALL DIFF(NRG(N),NRG(N+1),xMHR,CHR,CJK(1,N))
            CALL STORE(NTLIC,xMHC,CHC,NORB,THRES,CJK(1,N),N)
         END DO
C
         NJK=NTLIC-NHC
         nelm0=nelm0+1
         ma0=ma0+1
c
c     evaluate the integrals
c
         IF(LUSME.NE.0)WRITE(LUSME,*)ntlic, 1
         IF (ntlic > mxcorb) then
            write(*,6040) ntlic
            stop
         endif
         CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,ntlic,xmhc,
     &               CHC,mxcint,xint1e,nint1e,xint2e,elem,symtyp,map,
     &               ncorb,mcorb,1,lusme,enh_factor,enh_factor_type,
     &               ukrmolp_ints,NOEX)
         kcc=kcc+1
         IF(idiag.EQ.0)elem(1)=elem(1)+eone
         IF(LUSME.NE.0)WRITE(LUSME,*)ma0, ma0, 0
         IF(abs(elem(1)).GT.thrhm)THEN
            iembf=iembf+1
            nembf(1,iembf)=ma0
            nembf(2,iembf)=ma0
            em(iembf)=elem(1)
            IF(iembf.EQ.lembf)THEN
               CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
               nelm=nelm+iembf
               iembf=0
            END IF
         END IF
         NINT=NINT+NTLIC
         NINT0=NINT0+NTLIC
         NINT00=NINT00+NTLIC
      END DO
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c     The diagonal elements have now all been constructed
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C     Start the true OFF-DIAGONAL TERMS
C
c     The first portion of this is adapted from the scheme given by
c     L.A. Morgan & J. Tennyson, J Phys B, 26, 2429 (1993)
c     to expand prototype CSFs: numbering refers to appendix of
C
      NFLG=1
      MFLG=1
      CALL ENRGMC(CJK,NORBB)
      NJK=0
      MALM=0
C
c     loops 2000-2999 are for expanding prototype off-diagonal matrix
c     elements and/or for CI target contractions
c     skip if no neither to be done
      IF(iexpc.EQ.0 .AND. icitg.EQ.0)THEN
         nbeg=1
         ma0=0
         GO TO 4000
      END IF
      ictgtp=0
      mmm=1
      ma000=0
C-----------------------------------------------------------------------
c
c     EXPCSF / CI target contraction option: case (3)
c     compute remaining diagonal elements with continuum functions
c
C-----------------------------------------------------------------------
      DO ntgts=1, ntgsym
c
c     first: off-diagonal continuum case for present target state
c     (only one prototype)
         IF(icitg.NE.0)THEN
            nci=numtgt(ntgts)*(numtgt(ntgts)+1)/2
            ncl=numtgt(ntgts)
         ELSE
            nci=1
            ncl=1
         END IF
         IF(iexpc.EQ.0)THEN
            notop=notgt(ntgts)-1
            istep=notgt(ntgts)
         END IF
         DO iotgt1=0, notop
            mst=mmm+iotgt1-istep
            DO iotgt2=iotgt1+1, max(notop,1)
               ntlit=0
               m=mst
c     then loop over CI componants
c     first do off-diagonal term within same target state
c     starting with those for the same target CSF
               DO ictgt=1, nctgt(ntgts)
c
                  m=m+istep
                  MC=ICDT(M)
                  MD=INDT(M)
                  n=m+(iotgt2-iotgt1)
                  NC=ICDT(N)
                  ND=INDT(N)
                  CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
                  IF(NZ.GT.4)THEN
                     WRITE(nftw,9050)m, n, nz
 9050                FORMAT(/' Surely some mistake here:'/
     &                      ' Matrix element (',i5,',',i5,') has',i3,
     &                      ' differences'/
     &                      ' but come from same CI target CSF: STOP')
                     STOP
                  END IF
                  NJK=0
                  CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,MS,
     &                        NORBL,IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),
     &                        NODT(M),CDT(MC),NDT(MD),CJK,xMJK,CCJK,
     &                        SYMTYP,mcon,idiag,NOEX)
C
                  NHT=0
                  IF(NFLG.NE.1)THEN
                     NHT=NORBL
                     CALL CHECK(NHT,THRES,CH)
                  END IF
                  NHC=0
                  IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,THRES)
                  IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                     NALM=1
                     WRITE(NFTW,9040)N, M
 9040                FORMAT(I8,I7,5X,'FNS ARE NOT ORTHOGONAL')
                     RETURN
                  END IF
                  NFLG=1
                  ICH=NHC
                  IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,
     &                                    THRES)
                  IF(nhc.NE.0)THEN
c      IJK=NHC-ICH
                     nint0=nint0+nhc
                     IF(icitg.NE.0)THEN
                        CALL mkcfd(coef,ctgt(ictgtp+ictgt),numtgt(ntgts)
     &                             ,nctgt(ntgts))
                        CALL Comprs(nci,ntlit,nhc,xMHC,CHC,xMHC1,CH,
     &                              coef,mxcint)
                     ELSE
                        ntlit=nhc
                     END IF
                  END IF
C     now the terms for same target but different CSF
c     (only present for a CI target)
                  n=mst+(iotgt2-iotgt1)
                  DO jctgt=1, nctgt(ntgts)
                     n=n+istep
                     IF(jctgt.EQ.ictgt) CYCLE
                     NC=ICDT(N)
                     ND=INDT(N)
                     CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
                     IF(NZ.GT.4)CYCLE
                     NJK=0
                     CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,
     &                           MS,NORBL,IPOSIT,MPOS,NODT(N),CDT(NC),
     &                           NDT(ND),NODT(M),CDT(MC),NDT(MD),CJK,
     &                           xMJK,CCJK,SYMTYP,mcon,idiag,NOEX)
C
                     NHT=0
                     IF(NFLG.NE.1)THEN
                        NHT=NORBL
                        CALL CHECK(NHT,THRES,CH)
                     END IF
                     NHC=0
                     IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,
     &                  THRES)
                     IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                        NALM=1
                        WRITE(NFTW,9040)N, M
                        RETURN
                     END IF
                     NFLG=1
                     ICH=NHC
                     IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,
     &                  THRES)
                     IF(NHC.EQ.0)CYCLE
c      IJK=NHC-ICH
                     nint0=nint0+nhc
                     IF(icitg.NE.0)THEN
                        CALL mkcft(coef,ctgt(ictgtp+ictgt),
     &                             ctgt(ictgtp+jctgt),numtgt(ntgts),
     &                             nctgt(ntgts))
                        CALL Comprs(nci,ntlit,nhc,xMHC,CHC,xMHC1,CH,
     &                              coef,mxcint)
                     END IF
                  END DO
               END DO
            END DO
            ma00=ma000
            nelm0=nelm0+1
            nint00=nint00+ntlit
            IF(iexpc.NE.0)THEN
c        (3) off diagonal case within main block
               nchan=notgt(ntgts)
               nelem=nchan*(nchan-1)/2
               CALL EXPODL(kpt(m),kpt(n),nchan,nelem,ntlit,xmhc,xmhc,
     &                     ntlit)
            END IF
            ikk=1
c     evaluate the integrals
            DO i=1, nchan
               ma00=ma00+1
               na00=ma00
               DO j=i+1, nchan
                  na00=na00+1
                  IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
                  IF (ikk+ntlit > mxcorb) then
                     write(*,6040) ikk+ntlit
                     stop
                  endif
                  CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,
     &                        ntlit,xmhc(1,ikk),chc,mxcint,xint1e,
     &                        nint1e,xint2e,elem,symtyp,map,ncorb,mcorb,
     &                        nci,lusme,enh_factor,enh_factor_type,
     &                        ukrmolp_ints,NOEX)
                  kcc=kcc+1
                  ij=0
                  ma0=ma00-notgt(ntgts)
                  DO ii=1, numtgt(ntgts)
                     ma0=ma0+notgt(ntgts)
                     na0=na00+(ii-2)*notgt(ntgts)
                     DO jj=ii, numtgt(ntgts)
                        na0=na0+notgt(ntgts)
                        IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                        ij=ij+1
                        IF(abs(elem(ij)).GT.thrhm)THEN
                           iembf=iembf+1
                           nembf(1,iembf)=na0
                           nembf(2,iembf)=ma0
                           em(iembf)=elem(ij)
                           IF(iembf.EQ.lembf)THEN
                              CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                              NELM=NELM+iembf
                              iembf=0
                           END IF
                           IF(ii.EQ.jj)CYCLE
c        allow for upper triangle of diagonal block
                           iembf=iembf+1
                           nembf(1,iembf)=na0-j+i
                           nembf(2,iembf)=ma0+j-i
                           IF(LUSME.NE.0)WRITE(lusme,*)na0-j+i, ma0+j-i, 
     &                        1
                           em(iembf)=elem(ij)
                           IF(iembf.EQ.lembf)THEN
                              CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                              NELM=NELM+iembf
                              iembf=0
                           END IF
                        END IF
                     END DO
                  END DO
                  ikk=ikk+ntlit
               END DO
            END DO
            NINT=nint+NTLIT*nelem
C-----------------------------------------------------------------------
c
c      Case (2) Off diagonal element with L**2 CSF
c
C-----------------------------------------------------------------------
            IF(icitg.NE.0)nci=numtgt(ntgts)
            na0=ncont2
            DO N=ncont+1, nocsf
               NC=ICDT(N)
               ND=INDT(N)
               na0=na0+1
c     then loop over CI componants
               m=mmm+iotgt1-istep
               ntlit=0
               DO ictgt=1, nctgt(ntgts)
c
                  m=m+istep
                  MC=ICDT(M)
                  MD=INDT(M)
                  CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
                  IF(NZ.GT.4)CYCLE
                  NJK=0
                  CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,MS,
     &                        NORBL,IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),
     &                        NODT(M),CDT(MC),NDT(MD),CJK,xMJK,CCJK,
     &                        SYMTYP,mcon,0,NOEX)
C
                  NHT=0
                  IF(NFLG.NE.1)THEN
                     NHT=NORBL
                     CALL CHECK(NHT,THRES,CH)
                  END IF
                  NHC=0
                  IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,THRES)
                  IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                     NALM=1
                     WRITE(NFTW,9040)N, M
                     RETURN
                  END IF
                  NFLG=1
                  ICH=NHC
                  IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,
     &                                    THRES)
                  IF(NHC.EQ.0)CYCLE
                  nint0=nint0+nhc
                  IF(icitg.NE.0)THEN
                     CALL mkcfl(coef,ctgt(ictgtp+ictgt),ncl,nctgt(ntgts)
     &                          )
                     CALL Comprs(ncl,ntlit,nhc,xMHC,CHC,xMHC1,CH,coef,
     &                           mxcint)
                  ELSE
                     ntlit=nhc
                  END IF
               END DO
               nint00=nint00+ntlit
               nelm0=nelm0+1
               IF(iexpc.NE.0)THEN
                  nelem=notgt(ntgts)
                  CALL expcor(kpt(m),nelem,ntlit,xmhc,xmhc,ntlit)
               END IF
               ikk=1
               ma00=ma000
c     evaluate the integrals
               DO i=1, nelem
                  ma00=ma00+1
                  ma0=ma00-notgt(ntgts)
                  IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCL
                  IF (ikk+ntlit > mxcorb) then
                     write(*,6040) ikk+ntlit
                     stop
                  endif
                  CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,
     &                        ntlit,xmhc(1,ikk),chc,mxcint,xint1e,
     &                        nint1e,xint2e,elem,symtyp,map,NCORB,mcorb,
     &                        ncl,lusme,enh_factor,enh_factor_type,
     &                        ukrmolp_ints,NOEX)
                  kcc=kcc+1
                  DO ii=1, ncl
                     ma0=ma0+notgt(ntgts)
                     IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                     IF(abs(elem(ii)).GT.thrhm)THEN
                        iembf=iembf+1
                        nembf(1,iembf)=na0
                        nembf(2,iembf)=ma0
                        em(iembf)=elem(ii)
                        IF(iembf.EQ.lembf)THEN
                           CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                           NELM=NELM+iembf
                           iembf=0
                        END IF
                     END IF
                  END DO
                  ikk=ikk+ntlit
               END DO
               NINT=nint+NTLIT*nelem
            END DO
         END DO
         ictgtp=ictgtp+nctgt(ntgts)*numtgt(ntgts)
         ma000=ma000+notgt(ntgts)*numtgt(ntgts)
         IF(iexpc.EQ.0)THEN
            mmm=mmm+notgt(ntgts)*nctgt(ntgts)
         ELSE
            mmm=mmm+2*nctgt(ntgts)
         END IF
      END DO
c
c
C-----------------------------------------------------------------------
c     Off diagonal case with continuum in both bra and ket
c     and different target states: have to worry whether or not
c     the continuum is the same
c
C-----------------------------------------------------------------------
c
      ictgtp=0
      mmm=1
      ma000=0
      notop1=0
      notop2=0
      DO ntgt1=1, ntgsym-1
         jctgtp=ictgtp
         na000=ma000
         nnn=mmm
         DO ntgt2=ntgt1+1, ntgsym
            jctgtp=jctgtp+nctgt(ntgt2-1)*numtgt(ntgt2-1)
            na000=na000+notgt(ntgt2-1)*numtgt(ntgt2-1)
            IF(iexpc.EQ.0)THEN
               nnn=nnn+notgt(ntgt2-1)*nctgt(ntgt2-1)
               notop1=notgt(ntgt1)-1
               istep1=notgt(ntgt1)
               notop2=notgt(ntgt2)-1
               istep2=notgt(ntgt2)
               istep=min(istep1,istep2)
               notop=min(notop1,notop2)
            ELSE
               nnn=nnn+2*nctgt(ntgt2-1)
               nchana=notgt(ntgt1)
               nchanb=notgt(ntgt2)
               nchan=min(nchana,nchanb)
            END IF
            IF(icitg.NE.0)nci=numtgt(ntgt1)*numtgt(ntgt2)
c     if continua different jump to loop 3500
            IF(mcont(ntgt1).NE.mcont(ntgt2))GO TO 3500
            IF(gucont(ntgt1).NE.gucont(ntgt2))GO TO 3500
c
C-----------------------------------------------------------------------
c
c     EXPCSF / CI target contraction option: case (4)
c     continuum in bra and ket the same but target states different
c     (be careful because the expansion lengths may not be the same)
c
C-----------------------------------------------------------------------
c
C-----------------------------------------------------------------------
c           (4a) continuum orbital in bra and ket the same: expand
c           along diagonal of off-diagonal block
C-----------------------------------------------------------------------
            DO iotgt=0, notop
               ntlit=0
               m=mmm+iotgt-istep1
               DO ictgt1=1, nctgt(ntgt1)
c
                  m=m+istep1
                  MC=ICDT(M)
                  MD=INDT(M)
                  n=nnn+iotgt-istep2
                  DO ictgt2=1, nctgt(ntgt2)
                     n=n+istep2
                     NC=ICDT(N)
                     ND=INDT(N)
                     CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
                     IF(NZ.GT.4)CYCLE
                     NJK=0
                     CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,
     &                           MS,NORBL,IPOSIT,MPOS,NODT(N),CDT(NC),
     &                           NDT(ND),NODT(M),CDT(MC),NDT(MD),CJK,
     &                           xMJK,CCJK,SYMTYP,mcon,idiag,NOEX)
C
                     NHT=0
                     IF(NFLG.NE.1)THEN
                        NHT=NORBL
                        CALL CHECK(NHT,THRES,CH)
                     END IF
                     NHC=0
                     IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,
     &                  THRES)
                     IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                        NALM=1
                        WRITE(NFTW,9040)N, M
                        RETURN
                     END IF
                     NFLG=1
                     ICH=NHC
                     IF(NJK.GT.0)CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,
     &                  THRES)
                     IF(nhc.EQ.0)CYCLE
c      IJK=NHC-ICH
                     nint0=nint0+nhc
                     IF(icitg.NE.0)THEN
                        CALL mkcfg(coef,ctgt(ictgtp+1),ictgt1,
     &                             ctgt(jctgtp+1),ictgt2,numtgt(ntgt1),
     &                             numtgt(ntgt2),nctgt(ntgt1),
     &                             nctgt(ntgt2))
                        CALL Comprs(nci,ntlit,nhc,xMHC,CHC,xMHC1,CH,
     &                              coef,mxcint)
                     ELSE
                        ntlit=nhc
                     END IF
                  END DO
               END DO
               nelm0=nelm0+1
               nint00=nint00+ntlit
               IF(iexpc.NE.0)THEN
                  nelem=nchan
                  CALL expdg(kpt(m),nelem,ntlit,xmhc,xmhc,ntlit)
               END IF
               ikk=1
               ma00=ma000
               na00=na000
c     evaluate the integrals
               DO i=1, nchan
                  ma00=ma00+1
                  na00=na00+1
                  IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
                  IF (ikk+ntlit > mxcorb) then
                     write(*,6040) ikk+ntlit
                     stop
                  endif
                  CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,
     &                        ntlit,xmhc(1,ikk),CHC,mxcint,xint1e,
     &                        nint1e,xint2e,elem,symtyp,map,NCORB,mcorb,
     &                        nci,lusme,enh_factor,enh_factor_type,
     &                        ukrmolp_ints,NOEX)
                  kcc=kcc+1
                  ij=0
                  ma0=ma00-notgt(ntgt1)
                  DO ii=1, numtgt(ntgt1)
                     ma0=ma0+notgt(ntgt1)
                     na0=na00-notgt(ntgt2)
                     DO jj=1, numtgt(ntgt2)
                        na0=na0+notgt(ntgt2)
                        IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                        ij=ij+1
                        IF(abs(elem(ij)).GT.thrhm)THEN
                           iembf=iembf+1
                           nembf(1,iembf)=na0
                           nembf(2,iembf)=ma0
                           IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                           em(iembf)=elem(ij)
                           IF(iembf.EQ.lembf)THEN
                              CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                              NELM=NELM+iembf
                              iembf=0
                           END IF
                        END IF
                     END DO
                  END DO
                  ikk=ikk+ntlit
               END DO
               NINT=nint+NTLIT*nelem
            END DO
C-----------------------------------------------------------------------
c        (4b) off diagonal continuum lower triangle case: expand
c           along upper & lower triangle of off-diagonal block
C-----------------------------------------------------------------------
            DO iotgt1=0, notop1
               m=mmm+iotgt1-istep1
               DO iotgt2=1, max(1,notop2)
                  ntlit=0
                  DO ictgt1=1, nctgt(ntgt1)
c
                     m=m+istep1
                     MC=ICDT(M)
                     MD=INDT(M)
                     n=nnn+iotgt2-istep2
                     DO ictgt2=1, nctgt(ntgt2)
                        n=n+istep2
                        NC=ICDT(N)
                        ND=INDT(N)
                        CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
                        IF(NZ.GT.4)CYCLE
                        NJK=0
                        CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,
     &                              MM,MS,NORBL,IPOSIT,MPOS,NODT(N),
     &                              CDT(NC),NDT(ND),NODT(M),CDT(MC),
     &                              NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,
     &                              idiag,NOEX)
C
                        NHT=0
                        IF(NFLG.NE.1)THEN
                           NHT=NORBL
                           CALL CHECK(NHT,THRES,CH)
                        END IF
                        NHC=0
                        IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,
     &                     THRES)
                        IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                           NALM=1
                           WRITE(NFTW,9040)N, M
                           RETURN
                        END IF
                        NFLG=1
                        ICH=NHC
                        IF(NJK.GT.0)
     &                     CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,THRES)
                        IF(nhc.EQ.0)CYCLE
c      IJK=NHC-ICH
                        nint0=nint0+nhc
                        IF(icitg.NE.0)THEN
                           CALL mkcfg(coef,ctgt(ictgtp+1),ictgt1,
     &                                ctgt(jctgtp+1),ictgt2,
     &                                numtgt(ntgt1),numtgt(ntgt2),
     &                                nctgt(ntgt1),nctgt(ntgt2))
                           CALL Comprs(nci,ntlit,nhc,xMHC,CHC,xMHC1,CH,
     &                                 coef,mxcint)
                        ELSE
                           ntlit=nhc
                        END IF
                     END DO
                  END DO
                  nelm0=nelm0+1
                  nint00=nint00+ntlit
                  IF(iexpc.NE.0)THEN
                     nelem=nchana*nchanb-nchan
                     CALL expndg(kpt(m),kpt(n),nchana,nchanb,nelem,
     &                           ntlit,xmhc,xmhc,ntlit,ipair)
                  END IF
                  ikk=1
                  ma00=ma000
c     evaluate the integrals
                  DO i=1, nchana
                     ma00=ma00+1
                     na00=na000
                     DO j=1, nchanb
                        na00=na00+1
                        IF(i.EQ.j)CYCLE
                        IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
                        IF (ikk+ntlit > mxcorb) then
                           write(*,6040) ikk+ntlit
                           stop
                        endif
                        CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,
     &                          IPAIR,ntlit,xmhc(1,ikk),CHC,mxcint,
     &                          xint1e,nint1e,xint2e,elem,symtyp,
     &                          map,NCORB,mcorb,nci,lusme,enh_factor,
     &                          enh_factor_type,ukrmolp_ints,NOEX)
                        kcc=kcc+1
                        ij=0
                        ma0=ma00-notgt(ntgt1)
                        DO ii=1, numtgt(ntgt1)
                           ma0=ma0+notgt(ntgt1)
                           na0=na00-notgt(ntgt2)
                           DO jj=1, numtgt(ntgt2)
                              na0=na0+notgt(ntgt2)
                              IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                              ij=ij+1
                              IF(abs(elem(ij)).GT.thrhm)THEN
                                 iembf=iembf+1
                                 nembf(1,iembf)=na0
                                 nembf(2,iembf)=ma0
                                 em(iembf)=elem(ij)
                                 IF(iembf.EQ.lembf)THEN
                                    CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,
     &                                 EM)
                                    NELM=NELM+iembf
                                    iembf=0
                                 END IF
                              END IF
                           END DO
                        END DO
                        ikk=ikk+ntlit
                     END DO
                  END DO
                  NINT=nint+NTLIT*nelem
               END DO
            END DO
c
            CYCLE
c
C-----------------------------------------------------------------------
c
c     EXPCSF / CI target contraction option: case (5)
c     continuum in bra and ket and target states both different
c
C-----------------------------------------------------------------------
c
 3500       DO iotgt1=0, notop1
               DO iotgt2=0, notop2
                  ntlit=0
                  m=mmm+iotgt1-istep1
                  DO ictgt1=1, nctgt(ntgt1)
c
                     m=m+istep1
                     MC=ICDT(M)
                     MD=INDT(M)
                     n=nnn+iotgt2-istep2
                     DO ictgt2=1, nctgt(ntgt2)
                        n=n+istep2
                        NC=ICDT(N)
                        ND=INDT(N)
                        CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
                        IF(NZ.GT.4)CYCLE
                        NJK=0
                        CALL ENRGMZ(CH,xMHR,CHR,NDTRF,MDTR,NELT,NORB,MN,
     &                              MM,MS,NORBL,IPOSIT,MPOS,NODT(N),
     &                              CDT(NC),NDT(ND),NODT(M),CDT(MC),
     &                              NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,
     &                              idiag,NOEX)
C
                        NHT=0
                        IF(NFLG.NE.1)THEN
                           NHT=NORBL
                           CALL CHECK(NHT,THRES,CH)
                        END IF
                        NHC=0
                        IF(NH.GT.0)CALL MVDIAG(CH,xMHC1,CHR,xMHR,NH,NHC,
     &                     THRES)
                        IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
                           NALM=1
                           WRITE(NFTW,9040)N, M
                           RETURN
                        END IF
                        NFLG=1
                        ICH=NHC
                        IF(NJK.GT.0)
     &                     CALL MVDIAG(CH,xMHC1,CCJK,xMJK,NJK,NHC,THRES)
                        IF(nhc.EQ.0)CYCLE
c      IJK=NHC-ICH
                        nint0=nint0+nhc
                        IF(icitg.NE.0)THEN
                           CALL mkcfg(coef,ctgt(ictgtp+1),ictgt1,
     &                                ctgt(jctgtp+1),ictgt2,
     &                                numtgt(ntgt1),numtgt(ntgt2),
     &                                nctgt(ntgt1),nctgt(ntgt2))
                           CALL Comprs(nci,ntlit,nhc,xMHC,CHC,xMHC1,CH,
     &                                 coef,mxcint)
                        ELSE
                           ntlit=nhc
                        END IF
                     END DO
                  END DO
                  nelm0=nelm0+1
                  nint00=nint00+ntlit
                  IF(iexpc.NE.0)THEN
                     nelem=nchana*nchanb
                     CALL expoff(kpt(m),kpt(n),nchana,nchanb,nelem,
     &                           ntlit,xmhc,xmhc,ntlit)
                  END IF
                  ikk=1
                  ma00=ma000
c     evaluate the integrals
                  DO i=1, nchana
                     ma00=ma00+1
                     na00=na000
                     DO j=1, nchanb
                        na00=na00+1
                        IF(LUSME.NE.0)WRITE(LUSME,*)ntlit, NCI
                        IF (ikk+ntlit > mxcorb) then
                           write(*,6040) ikk+ntlit
                           stop
                        endif
                        CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,
     &                          IPAIR,ntlit,xmhc(1,ikk),CHC,mxcint,
     &                          xint1e,nint1e,xint2e,elem,symtyp,
     &                          map,NCORB,mcorb,nci,lusme,enh_factor,
     &                          enh_factor_type,ukrmolp_ints,NOEX)
                        kcc=kcc+1
                        ij=0
                        ma0=ma00-notgt(ntgt1)
                        DO ii=1, numtgt(ntgt1)
                           ma0=ma0+notgt(ntgt1)
                           na0=na00-notgt(ntgt2)
                           DO jj=1, numtgt(ntgt2)
                              na0=na0+notgt(ntgt2)
                              IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
                              ij=ij+1
                              IF(abs(elem(ij)).GT.thrhm)THEN
                                 iembf=iembf+1
                                 nembf(1,iembf)=na0
                                 nembf(2,iembf)=ma0
                                 em(iembf)=elem(ij)
                                 IF(iembf.EQ.lembf)THEN
                                    CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,
     &                                 EM)
                                    NELM=NELM+iembf
                                    iembf=0
                                 END IF
                              END IF
                           END DO
                        END DO
                        ikk=ikk+ntlit
                     END DO
                  END DO
                  NINT=nint+NTLIT*nelem
c
               END DO
            END DO
         END DO
         ictgtp=ictgtp+nctgt(ntgt1)*numtgt(ntgt1)
         ma000=ma000+notgt(ntgt1)*numtgt(ntgt1)
         IF(iexpc.EQ.0)THEN
            mmm=mmm+notgt(ntgt1)*nctgt(ntgt1)
         ELSE
            mmm=mmm+2*nctgt(ntgt1)
         END IF
      END DO
      nbeg=ncont+1
      ma0=ncont2
C-----------------------------------------------------------------------
C
C     Off-diagonal L**2 matrix elements: no expansion possible
C
C-----------------------------------------------------------------------
 4000 CONTINUE
      DO M=nbeg, NOCSF-1
         MC=ICDT(M)
         MD=INDT(M)
         ma0=ma0+1
         na0=ma0
         jdiag=idiag
         DO N=M+1, NOCSF
            NC=ICDT(N)
            ND=INDT(N)
            na0=na0+1
            CALL PZERO(MN,NDT(ND),NDT(MD),NZ,NBB)
            IF(NZ.GT.4)CYCLE
            NJK=0
            IF(n.GT.ncont)jdiag=0
            CALL ENRGMZ(CH,xMHC,CHR,NDTRF,MDTR,NELT,NORB,MN,MM,MS,NORBL,
     &                  IPOSIT,MPOS,NODT(N),CDT(NC),NDT(ND),NODT(M),
     &                  CDT(MC),NDT(MD),CJK,xMJK,CCJK,SYMTYP,mcon,jdiag,
     &                  NOEX)
            IF(MALM.NE.0)THEN
               MALM=0
               WRITE(NFTW,9020)NH, NJK, MALM, NJKMX
               NALM=1
               RETURN
            END IF
C
            NHT=0
            IF(NFLG.NE.1)THEN
               NHT=NORBL
               CALL CHECK(NHT,THRES,CH)
            END IF
            NHC=0
            IF(NH.GT.0)CALL MVDIAG(CHR,xMHC,CHR,xMHC,NH,NHC,THRES)
            IF(NFLG.NE.1 .AND. NHT+NHC.NE.0)THEN
               NALM=1
               WRITE(NFTW,9040)N, M
               RETURN
            END IF
            NFLG=1
            ICH=NHC
            IF(NJK.GT.0)CALL MVDIAG(CHR,xMHC,CCJK,xMJK,NJK,NHC,THRES)
            IF(NHC.EQ.0)CYCLE
c      IJK=NHC-ICH
            nelm0=nelm0+1
c
            IF(LUSME.NE.0)WRITE(LUSME,*)NHC, 1
            IF (nhc > mxcorb) then
               write(*,6040) nhc
               stop
            endif
            CALL PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,NHC,xMHC,
     &                  CHR,mxcint,xint1e,nint1e,xint2e,elem,symtyp,map,
     &                  NCORB,mcorb,1,lusme,enh_factor,enh_factor_type,
     &                  ukrmolp_ints,NOEX)
            kcc=kcc+1
            IF(LUSME.NE.0)WRITE(LUSME,*)na0, ma0, 0
            IF(abs(elem(1)).GT.thrhm)THEN
               iembf=iembf+1
               nembf(1,iembf)=na0
               nembf(2,iembf)=ma0
               em(iembf)=elem(1)
               IF(iembf.EQ.lembf)THEN
                  CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
                  NELM=NELM+iembf
                  iembf=0
               END IF
            END IF
            NINT=NINT+NHC
            NINT0=NINT0+NHC
            NINT00=NINT00+NHC
         END DO
      END DO
      IF(LUSME.NE.0)THEN
         WRITE(lusme,*)headerp
         WRITE(lusme,*)KCC-1
      END IF
c     write out last buffer of the Hamiltonian matrix
 5000 IF(iembf.NE.0)THEN
         CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
         NELM=NELM+iembf
         iembf=0
      END IF
c     write dummy record at end of file
      CALL WRTEM(NFTE,iembf,LEMBF,NEMBF,EM)
c
      N=(NOCSF*(NOCSF+1))/2
      m=(mOCSF*(mOCSF+1))/2
      WRITE(NFTW,6000)m, NELM
 6000 FORMAT(//' TOTAL H Matrix     ELEMENTS      =',
     &       I10/' NON-ZERO ELEMENTS evaluated      =',I10)
      IF(iexpc.NE.0)WRITE(nftw,6010)n, nelm0
 6010 FORMAT(' Total    prototype ELEMENTS      =',
     &       I10/' Non-zero prototype ELEMENTS      =',I10)
      IF(iexpc.NE.0)WRITE(nftw,6020)nint0, nint00
 6020 FORMAT(/' Number (prototype) integrals     =',
     &       I10/' Compressed (prototype) integrals =',I10)
      WRITE(nftw,6030)nint
 6030 FORMAT(' Number integrals evaluated       =',I10)
 6040 FORMAT(' ENRGMX error: MXCORB too small, need: ',I15)
c
      DEALLOCATE(xmjk,ccjk,xmhr,chr,chc,xmhc,xmhc1,elem)
c
c**************************************************************************
c     Deallocate statment for the array elem has been introduced in
c     subroutines ENRGMX (see stament just above) and ENRGMS where this
c     array is allocated. Also the if statment "if(icitg.ne.0) deallocate
c     (xmhc1)" in ENRGMX was replaced by a deallocate statment only.
c     xmhc1(2,mxcint) is allocated for both cases,(icitg.ne.0) and
c     (icitg.eq.0), in this subroutine. (N.Vinci, 15/11/01).
c**************************************************************************
c
      RETURN
      END SUBROUTINE ENRGMX
!*==enrgmz.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE ENRGMZ(CH,XMHPQ,CHPQ,NDTR,MDTR,NELT,NORB,MN,MM,MS,
     &                  NORBL,IPOSIT,MPOS,NODA,CA,NDA,NODB,CB,NDB,CJKAB,
     &                  XMJK,CJK,SYMTYP,mcon,jdiag,NOEX)
C
C     NODA          NO OF DTRS IN WF(A)
C     NODB          NO OF DTRS IN WF(B)
C     CA            COEFFICINTS FOR WF(A)
C     CB            COEFFICINTS FOR WF(B)
C     NDA           PACKED DTRS FOR WF(A)
C     NDB           PACKED DTRS FOR WF(B)
C     MDTR(NSRB)    REF. DTR IN EXPONDED FORM
C     NFLG          =0, DIAGONAL, =1,OFF DIAGONAL
C
C     TEMP  NDTA,NDTB,NDTD(NELT)
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE consts, ONLY : ZERO=>XZERO
      USE integer_packing, ONLY : pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: LPOSIT, MFLG, NFLG, NHPQ
      COMMON /ENA   / NHPQ, NFLG, MFLG
      COMMON /GD4   / LPOSIT
C
C Dummy arguments
C
      INTEGER :: IPOSIT, JDIAG, NELT, NODA, NODB, NORB, NORBL, SYMTYP,
     &           NOEX
      REAL(KIND=wp), DIMENSION(*) :: CA, CB, CH, CHPQ, CJK
      REAL(KIND=wp), DIMENSION(norb*(norb+1)/2,4) :: CJKAB
      INTEGER, DIMENSION(*) :: MCON, MDTR, MM, MN, MPOS, MS, NDA, NDB, 
     &                         NDTR
      INTEGER(longint), DIMENSION(2,*) :: XMHPQ, XMJK
      INTENT (IN) CA, CB, JDIAG, MCON, MDTR, NDA, NDB, NDTR, NELT, NODA, 
     &            NODB, NORBL
      INTENT (INOUT) CH, CHPQ
C
C Local variables
C
      REAL(KIND=wp), SAVE :: CDA, CFD
      INTEGER, SAVE :: I, IA1, IA2, IC, ID, J, JA, JB, LNADD, M, MA, 
     &                 MAA, MB, MBB, MDA, MDB, N, NB, NCHK, NDBMPI, 
     &                 NDTAI
      INTEGER, DIMENSION(2), SAVE :: NDD
      INTEGER, DIMENSION(4), SAVE :: NDIAG, NDTC
      INTEGER, DIMENSION(nelt) :: NDTA, NDTB, NDTD
      INTEGER(longint), DIMENSION(2), SAVE :: XIAB
C
C*** End of declarations rewritten by SPAG
C
      LPOSIT=0
      LNADD=0
      DO I=1, NELT
         NDTD(I)=NDTR(I)
      END DO
C
      NHPQ=0
      DO I=1, NORBL
         CH(I)=ZERO
      END DO
C
      MDA=1
      DO IC=1, NODA
         CDA=CA(IC)
         MA=NDA(MDA)
         MAA=MDA+MA
         MDB=1
C
         DO ID=1, NODB
            LPOSIT=0
            CFD=CDA*CB(ID)
            MB=NDB(MDB)
            MBB=MDB+MB
C
CDIR% NOVECTOR
            IF(MA*MB.NE.0)GO TO 200
            IF(MAX(MA,MB).GT.2)GO TO 301
            IF(MA.EQ.MB)GO TO 270
            IF(MA.EQ.0)THEN
               IF(MB.LE.1)THEN
                  NDTC(3)=NDB(MDB+1)
                  NDTC(4)=NDB(MBB+1)
                  GO TO 260
               ELSE
                  NDTC(1)=NDB(MDB+1)
                  NDTC(2)=NDB(MBB+1)
                  NDTC(3)=NDB(MDB+2)
                  NDTC(4)=NDB(MBB+2)
               END IF
            ELSE
               IF(MA.LE.1)THEN
                  NDTC(3)=NDA(MAA+1)
                  NDTC(4)=NDA(MDA+1)
                  GO TO 260
               ELSE
                  NDTC(1)=NDA(MAA+1)
                  NDTC(2)=NDA(MDA+1)
                  NDTC(3)=NDA(MAA+2)
                  NDTC(4)=NDA(MDA+2)
               END IF
            END IF
            GO TO 250
C
 200        DO I=1, MA
               NDTA(I)=NDA(MAA+I)
               NDTB(I)=NDA(MDA+I)
            END DO
C
            JA=MA
            DO I=1, MB
               NDBMPI=NDB(MDB+I)
               DO J=1, MA
                  IF(NDTB(J).EQ.NDBMPI)THEN
                     NDTB(J)=NDB(MBB+I)
                     GO TO 208
                  END IF
               END DO
               JA=JA+1
               NDTA(JA)=NDBMPI
               NDTB(JA)=NDB(MBB+I)
 208        END DO
C
            IF(JA-MA.GT.2)GO TO 301
            JB=0
C
            DO I=1, JA
               NDTAI=NDTA(I)
               DO J=1, JA
                  IF(NDTB(J).EQ.NDTAI)THEN
                     IF(I.NE.J)THEN
                        NDTB(J)=NDTB(I)
                        NDTB(I)=NDTAI
                        CFD=-CFD
                     END IF
                     GO TO 220
                  END IF
               END DO
               JB=JB+1
               IF(JB.GT.2)GO TO 301
               NDD(JB)=I
 220        END DO
C
            IF(JB.EQ.0)GO TO 270
            IF(JB.EQ.1)THEN
               NDTC(3)=NDTA(NDD(1))
               NDTC(4)=NDTB(NDD(1))
               GO TO 260
            ELSE
               NDTC(1)=NDTA(NDD(1))
               NDTC(2)=NDTB(NDD(1))
               NDTC(3)=NDTA(NDD(2))
               NDTC(4)=NDTB(NDD(2))
            END IF
C
C     TWO PAIRS ARE DIFFERENT
C
 250        IF(jdiag.GT.1)THEN
               IF(min(mcon(ndtc(1)),mcon(ndtc(2)),mcon(ndtc(3)),
     &            mcon(ndtc(4))).GT.0)GO TO 301
            END IF
            IF(IPOSIT.NE.0)CALL POSADD(LPOSIT,NDTC,IPOSIT,1,4,MPOS)
            CALL ENRGMB(CJKAB(1,1),CJKAB(1,2),CJKAB(1,3),CJKAB(1,4),MN,
     &                  MM,MS,MPOS,NDTC,CFD,XMJK,CJK,SYMTYP,NOEX)
            GO TO 301
C
C     ONE PAIR IS DIFFERENT
C
 260        DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=NDA(MAA+I)
            END DO
C
            IF(jdiag.GT.1 .AND. min(mcon(ndtc(3)),mcon(ndtc(4))).GT.0)
     &         THEN
               nchk=1
            ELSE
               nchk=0
            END IF
            DO I=1, NELT
               IF(NDTD(I).EQ.NDTC(3))CYCLE
               IF(nchk.EQ.1 .AND. mcon(ndtd(i)).GT.0)CYCLE
               NDTC(1)=NDTD(I)
               NDTC(2)=NDTD(I)
               IF(IPOSIT.NE.0)CALL POSADD(LPOSIT,NDTC,IPOSIT,1,4,MPOS)
               CALL ENRGMB(CJKAB(1,1),CJKAB(1,2),CJKAB(1,3),CJKAB(1,4),
     &                     MN,MM,MS,MPOS,NDTC,CFD,XMJK,CJK,SYMTYP,NOEX)
            END DO
C
            IF(IPOSIT.NE.0)CALL POSADD(LNADD,NDTC,NORB,3,4,MPOS)
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=N
            END DO
            IF(nchk.EQ.1)GO TO 301
            IF(MM(NDTC(3)).NE.MM(NDTC(4)) .OR. MS(NDTC(3))
     &         .NE.MS(NDTC(4)))GO TO 301
            IF(MN(NDTC(3)).LT.MN(NDTC(4)))THEN
               IA1=MN(NDTC(4))
               IA2=MN(NDTC(3))
            ELSE
               IA1=MN(NDTC(3))
               IA2=MN(NDTC(4))
            END IF
C
            IF(IA1.EQ.IA2 .AND. NFLG.EQ.0)THEN
C BELOW IS AN EXAMPLE OF USE OF ARRAY DIMENSIONS TO LABEL A POSITRONIC
C INTEGRAL - CONVERTED TO SIGN LABELLING IN STDIAG.
               CH(IA1+LNADD)=CH(IA1+LNADD)+CFD
               GO TO 301
            END IF
CDIR% VECTOR
            CALL pack8ints(0,ia1,0,ia2,lnadd,0,0,0,xiab)
c     IAB=IPACK2(IA1,IA2)
            i=ijkpqrs(xiab,xmhpq,nhpq)
c
            IF(I.GT.0)THEN
               CHPQ(I)=CHPQ(I)+CFD
            ELSE
               NHPQ=NHPQ+1
               XMHPQ(1,NHPQ)=XIAB(1)
               XMHPQ(2,NHPQ)=XIAB(2)
               CHPQ(NHPQ)=CFD
            END IF
            GO TO 301
C
C     THE SAME DETERMINANTS
C
 270        DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=NDA(MAA+I)
            END DO
            DO M=2, NELT
               NDTC(1)=NDTD(M)
               NDTC(2)=NDTC(1)
               DO N=1, M-1
                  NDTC(3)=NDTD(N)
                  NDTC(4)=NDTC(3)
                  IF(jdiag.GT.1)THEN
                     IF(min(mcon(ndtc(1)),mcon(ndtc(3))).GT.0)CYCLE
                  END IF
                  IF(IPOSIT.NE.0)
     &               CALL POSADD(LPOSIT,NDTC,IPOSIT,1,4,MPOS)
                  CALL ENRGMB(CJKAB(1,1),CJKAB(1,2),CJKAB(1,3),
     &                        CJKAB(1,4),MN,MM,MS,MPOS,NDTC,CFD,XMJK,
     &                        CJK,SYMTYP,NOEX)
               END DO
            END DO
C
            IF(NFLG.NE.0)NFLG=2
            DO M=1, NELT
               N=NDTD(M)
               IF(jdiag.GT.1 .AND. mcon(n).GT.0)CYCLE
               NDIAG(1)=N
               NDIAG(2)=N
               IF(IPOSIT.NE.0)CALL POSADD(LPOSIT,NDIAG,NORB,1,2,MPOS)
               NB=MN(N)
               CH(NB+LNADD)=CH(NB+LNADD)+CFD
            END DO
C
            DO I=1, MA
               N=NDA(MDA+I)
               M=MDTR(N)
               NDTD(M)=N
            END DO
C
 301        MDB=MBB+MB+1
         END DO
C
         MDA=MAA+MA+1
      END DO
      RETURN
      END SUBROUTINE ENRGMZ
!*==exdrvf.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE EXDRVF
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE global_utils, ONLY : MPROD, CWBOPN
      USE scatci_data, ONLY : NTGTMX
      USE mpi_gbl, ONLY: mpi_mod_start, mpi_mod_finalize
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      CHARACTER(LEN=8) :: BLANK
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CDO, EIG, VECS, X1E, 
     &            X2E
      REAL(KIND=wp) :: CRITC, CRITE, CRITES, CRITR, ORTHO, PIN, R, S, 
     &                SCALEM, SZ, THRES, THRHM, enh_factor
      INTEGER, DIMENSION(ntgtmx) :: GUCONT, MCONT, NCTGT, NOTEMP, NOTGT, 
     &                              NTGTF, NTGTS, NUMTGT
      INTEGER :: I, ICIDG, ICITG, IDIAG, IDIAGT, IEXPC, IGH, IODR, 
     &           IPOSIT, ISD, ISS, ISTART, ISYM, IU, JUMP, JUNK, KEIGEN, 
     &           KVECS, LCDOF, LEMBF, LNDOF, LUSME, MAXITER, MEGUL, 
     &           MGVN, MOCSF, MXCSF, MXRW, NALM, NCITOT, NCONT, NCONT2, 
     &           NCORB, NCTARG, NELM, NELT, NFTE, NFTG, NFTI, NFTW, 
     &           NINT1E, NINT2E, NIPAIR, NMPRS, NN0, NOBEP, NOBMX, NOBT, 
     &           NOCSF, NORB, NORBP, NPQRS, NSRB, NSYM, NSYM0, 
     &           NTGCON, NTGSYM, NTGT, NUMMX, SCFUSE, SYMTYP, j,k,m,
     &           size_kphz, enh_factor_type, NOEX
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KICDO, KINDO, KMCON, KMCORB, 
     &       KMG, KMM, KMN, KMS, KNDO, KNORB, KNSRB, KPAIR, KPHZ, KPT, 
     &       NODO
      INTEGER, DIMENSION(5,255) :: LPQRS
      INTEGER, DIMENSION(21) :: MBAS, NBASH
      INTEGER(longint), DIMENSION(2,1023) :: MPQRS
      CHARACTER(LEN=120) :: NAME, OLDNAM
      INTEGER, DIMENSION(670) :: NBAS
      INTEGER, DIMENSION(200) :: NDTRF
      INTEGER, DIMENSION(20) :: NOB, NOB0, NOB0L, NOBC, NOBE, NOBL, 
     &                          NOBP, NOBV
      INTEGER, DIMENSION(6) :: NPFLG
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NRI, NSM
      LOGICAL :: QMOLN
      LOGICAL :: ukrmolp_ints
      INTEGER :: lcont !the number of the last continuum CSF after contraction
      LOGICAL :: write_ham_pars
C
C*** End of declarations rewritten by SPAG
C
      NAMELIST /input / nfte, nfti, lembf, ncont, scalem, thrhm, name, 
     &   ncorb, iodr, iexpc, icitg, icidg, nobc, idiag, ntgsym, nctgt, 
     &   notgt, mcont, gucont, iposit, numtgt, ntgtf, ntgts, nftg, 
     &   megul, scfuse, CRITES, CRITC, CRITR, ORTHO, maxiter, igh,
     &   lusme, qmoln, enh_factor, enh_factor_type, ukrmolp_ints, NOEX
C
      DATA NFTW/6/, MEGUL/13/, iexpc/0/, icitg/0/, icidg/1/
      DATA nfte/26/, nfti/16/, lembf/5000/, ncont/0/, scalem/1.0_wp/
      DATA THRHM/1.0E-15_wp/, ncorb/-1/, iodr/1/, nalm/0/, NMPRS/0/, 
     &     nelm/0/
c     set namelist defaults
      DATA ntgt/1/, nctgt/ntgtmx*1/, mcont/ntgtmx*0/, gucont/ntgtmx*0/, 
     &     ntgsym/1/, numtgt/ntgtmx*1/, nftg/0/, BLANK/'        '/,
     &     NTGTF/NTGTMX*0/, NTGTS/NTGTMX*0/, NOTGT/NTGTMX*0/, SCFUSE/0/, 
     &     nobc/20*0/, idiag/-1/, iposit/0/
      DATA CRITE/1D-12/, CRITC/1D-10/, CRITR/1D-8/, ORTHO/1D-7/, 
     &     igh/-99/, maxiter/-1/
      DATA nobe, nobp, nobv/20*0, 20*0, 20*0/, LUSME/0/
      DATA QMOLN/.FALSE./
      DATA ukrmolp_ints/.TRUE./
      DATA enh_factor_type/0/
      DATA enh_factor/1.0_wp/
      DATA NOEX/0/
C

      name=blank
      READ(5,input)
c
c   Trap some inconsistent inputs related to enhacement factors
      if (icitg.eq.0. and.enh_factor_type.ne.0) then
       stop "You should not be using an enhancement factor in target 
     & calculations"
      endif
c
      if (enh_factor_type.eq.0. and.enh_factor.ne.1.0_wp) then
       stop "If you want to use an enhancement factor you should also 
     & include enh_factor_type in your input"
      endif
      oldnam=blank
      CALL cwbopn(megul)
      READ(megul)oldnam, MGVN, S, SZ, R, PIN, NORB, NSRB, NOCSF, NELT, 
     &           lcdof, IDIAGT, NSYM, SYMTYP, lndof, npflg, thres, 
     &           nctarg, ntgcon
      WRITE(NFTW,14)'CONGEN', oldnam
C
      IF(name.EQ.blank)name=oldnam
      WRITE(NFTW,14)'SCATCI', NAME
 14   FORMAT(' ',A6,': ',A)
      IF(idiag.LT.0)idiag=idiagt
      WRITE(NFTW,16)MGVN, S, NELT, SZ, IDIAG, R, THRES, NSYM, MEGUL, 
     &              NOCSF, (NPFLG(I),I=1,6)
 16   FORMAT(' MGVN =',I10,5X,'S    =',F5.1,/,' NELT =',I10,5X,'SZ   =',
     &       F5.1,/,' IDIAG=',I10,5X,'R    =',F5.1,/,' THRES=',D10.1,5X,
     &       'NSYM =',I5,/,' MEGUL=',I3,', NOCSF=',I10,/,' NPFLG',I10,
     &       5I4)
      IF(nocsf.LE.0)THEN
         WRITE(nftw,116)
 116     FORMAT(' No CSFS selected: stop')
         GO TO 38
      END IF
c
      IF(ntgcon.GT.0)ntgsym=ntgcon
      ALLOCATE(kphz(nocsf),nodo(nocsf))
      CALL rdnfto(megul,nob,nob0,nobl,nob0l,nsym,ndtrf,nelt,nodo,nocsf,
     &            kphz,nctarg,nctgt,notemp,mcont,gucont,ntgsym,ntgcon,
     &            nobe,nobp,nobv,iposit,NOEX)
 
 
 
c     checking the NOB-values
c     NOBE(i) number of electronic orbitals
c             default nobe(i)=nob(i)
c     NOBP(i) number of positronic orbitals
c             default nobp(i)=0
c     NOBV(i) at the moment not used, but
c             printed on output
c             default nobv(i)=nob0(i)
      DO i=1, nsym
         IF(nobp(i).EQ.0)THEN
            nobe(i)=nob(i)
 
            WRITE(NFTW,*)'NOB(i)=', nob(i)
            WRITE(NFTW,*)'setting NOBE(i)=', nobe(i)
 
         END IF
         nobep=nobe(i)+nobp(i)
         IF(nobep.NE.nob(i))THEN
            WRITE(NFTW,*)'ERROR on input:'
            WRITE(NFTW,*)'not: NOB(i)=NOBE(i)+NOBP(i)'
            WRITE(NFTW,*)'i=', i
            WRITE(NFTW,*)'NOB(i)=', nob(i)
            WRITE(NFTW,*)'NOBE(i)=', nobe(i)
            WRITE(NFTW,*)'NOBP(i)=', nobp(i)
         END IF
         IF(nobv(i).EQ.0)THEN
            nobv(i)=nob0(i)
         END IF
         IF(nob0(i).GT.nob(i))THEN
            WRITE(NFTW,*)'ERROR on input:'
            WRITE(NFTW,*)'NOB0(i) > NOB(i)'
            WRITE(NFTW,*)'i=', i
            WRITE(NFTW,*)'NOB(i)=', nob(i)
            WRITE(NFTW,*)'NOB0(i)=', nob0(i)
         END IF
      END DO
 
 
 
 
      ntgcon=max(0,ntgcon)
      IF(ntgsym.GT.1)ntgcon=ntgsym
c
      nsym0=nsym
      IF(symtyp.EQ.1)THEN
         nsym0=nsym0+nsym
         IF(nobl(nsym0).LE.0)nsym0=nsym0-1
      END IF
      WRITE(NFTW,17)(NOBL(I),I=1,NSYM0)
      WRITE(NFTW,170)(NOB0L(I),I=1,NSYM0)
c     determine start of continuum orbital expansion if needed
      IF(iexpc.GT.0)THEN
         DO i=1, ntgsym
            IF(symtyp.NE.1)THEN
               isym=mcont(i)+1
            ELSE
               isym=mcont(i)+mcont(i)+1
               IF(mod(mcont(i),2).EQ.0 .AND. gucont(i).EQ.-1)isym=isym+1
               IF(mod(mcont(i),2).EQ.1 .AND. gucont(i).EQ.1)isym=isym+1
            END IF
            IF(nobc(isym).GT.0)CYCLE
            nobt=nob(isym)-notgt(i)
            IF(nobt.LT.nob0l(isym))THEN
               WRITE(nftw,1045)
               WRITE(nftw,1022)i, isym, nobt, nob0l(isym)
 1022          FORMAT(' For target symmetry type number',i3,
     &                ' and symmetry',i3,/
     &                ' Attempt to start continuum expansion at orbital'
     &                ,i3,/' But NOB0 =',i3,'.  Input data as follows:')
               WRITE(nftw,1010)ntgsym, (notgt(iu),iu=1,ntgsym)
               WRITE(nftw,1020)(nctgt(iu),iu=1,ntgsym)
               WRITE(nftw,1025)(numtgt(iu),iu=1,ntgsym)
               WRITE(nftw,1030)(mcont(iu),iu=1,ntgsym)
               IF(gucont(1).NE.0)WRITE(nftw,1040)
     &                                 (gucont(iu),iu=1,ntgsym)
               STOP
            END IF
            nobc(isym)=nobt
         END DO
         WRITE(NFTW,171)(NOBC(I),I=1,NSYM0)
      ELSE IF(nctarg.GT.0)THEN
         DO i=1, ntgsym
            notgt(i)=notemp(i)
         END DO
      END IF
c
      WRITE(NFTW,18)(NDTRF(I),I=1,NELT)
      WRITE(NFTW,102)NSRB, LNDOF, LCDOF
C
 17   FORMAT(' NOB  =',I10,20I4)
 170  FORMAT(' NOB0 =',I10,20I4)
 171  FORMAT(' NOBC =',I10,20I4)
 18   FORMAT(' NDTRF=',I10,19I4/(13X,20I4))
 102  FORMAT(' NSRB=',I5,4X,'LNDOF=',I10,4X,'LCDOF=',I10)
      WRITE(nftw,190)nfte, nfti, lembf, iodr
 190  FORMAT(' nfte =',i5,8x,'nfti =',i6,7x,'lembf=',i6,7x,'iodr =',i4)
C
      IF(SYMTYP.GE.2)THEN
         WRITE(NFTW,1910)SYMTYP
 1910    FORMAT(' MOLECULE SYMMETRY CASE,  SYMTYP =',I2)
         JUNK=MPROD(1,1,NPFLG(6),NFTW)
      END IF
      IF(iexpc.NE.0)WRITE(nftw,1005)
 1004 FORMAT(/' Error: CI target option can only run with IEXPC =1')
 1005 FORMAT(/' Prototype CSF option requested')
 1006 FORMAT(/' CI target option requested,'/
     &       ' Vectors to be read from unit NFTG =',i3)
 1007 FORMAT(/' CI target option requested,'/
     &       ' Vectors to be constructed using CONGEN input')
 1008 FORMAT(/' Error: inconsistant options:',/' NCORB =',i3,
     &       ', NFTE =,',i3)
 1009 FORMAT(' Phases for target vectors read in from CONGEN')
      IF(icitg.NE.0)THEN
         IF(nftg.GT.0)THEN
            WRITE(nftw,1006)nftg
            IF(nctarg.GT.0)WRITE(nftw,1009)
         END IF
         IF(nftg.LE.0)WRITE(nftw,1007)
         IF(ncorb.GE.0)THEN
            WRITE(nftw,1008)ncorb, nfte
            GO TO 38
         END IF
         IF(iexpc.EQ.0)THEN
            WRITE(6,1004)
            GO TO 38
         END IF
      END IF
c     print namelist input
      IF(iexpc.NE.0 .OR. icitg.NE.0)THEN
         WRITE(nftw,1010)ntgsym, (notgt(i),i=1,ntgsym)
 1010    FORMAT(/' Number of target symmetries in expansion,   NTGSYM ='
     &          ,i5/
     &          ' Number of continuum orbitals for each state, NOTGT =',
     &          20I5,/,(' ',20I5))
         WRITE(nftw,1020)(nctgt(i),i=1,ntgsym)
 1020    FORMAT(' Number of CI components for each state,      NCTGT =',
     &          20I10,/,('  ',20I5))
         WRITE(nftw,1025)(numtgt(i),i=1,ntgsym)
 1025    FORMAT(' Number of target states of each symmetry,   NUMTGT =',
     &          20I5,/,('  ',20I5))
         WRITE(nftw,1030)(mcont(i),i=1,ntgsym)
 1030    FORMAT(' Continuum M projection  for each state,      MCONT =',
     &          20I5,/,('  ',20I5))
         IF(gucont(1).NE.0)WRITE(nftw,1040)(gucont(i),i=1,ntgsym)
 1040    FORMAT(' Continuum G/U symmetry  for each state,     GUCONT =',
     &          20I5,/,('  ',20I5))
      END IF
      IF(iexpc.NE.0)THEN
c       error check
         DO i=1, ntgsym
            IF(nctgt(i).LE.0 .OR. notgt(i).LT.2)THEN
               WRITE(nftw,1045)
 1045          FORMAT(' Inconsistency detected in EXPCSF data')
               GO TO 38
            END IF
            IF(idiag.EQ.0)THEN
               WRITE(nftw,1046)
 1046          FORMAT(
     &               ' IDIAG=0 option not available with Prototype CSFs'
     &               )
               GO TO 38
            END IF
         END DO
      END IF
      IF(icidg.NE.0)WRITE(nftw,1055)
 1055 FORMAT(/' Diagonalisation of Hamiltonian matrix requested')
      IF(icidg.EQ.2)GO TO 5000
c
      ISS=S+S
      ISD=SZ+SZ
      IF(ISS.LT.ISD)THEN
         WRITE(NFTW,40)
 40      FORMAT('  S.LT.SZ')
         GO TO 38
      END IF
      nn0=0
      DO I=1, NSYM
         nn0=nn0+nob0l(i)
      END DO
      IF(((nftg.EQ.0 .AND. icitg.EQ.1) .OR. idiag.GE.2 .OR. ncorb.EQ.0)
     &   .AND. nn0.LE.0)THEN
         WRITE(nftw,42)idiag, ncorb, nftg
 42      FORMAT(//'  Parameter NOB0 must be set in CONGEN for:',/
     &         '  Option to neglect target configurations, IDIAG=2 or 3'
     &         ,/
     &       '  Phase corrected target CI option with NCORB=0 or NFTG=0'
     &       ,/'  Present values: IDIAG',i3,' NCORB',i3,' NFTG',i3/)
         GO TO 38
      END IF
      IF(IPOSIT.NE.0)THEN
         NORBP=NORB/2
      END IF
c
C-----------------------------------------------------------------------
      nummx=1
      ncitot=0
      ALLOCATE(kpt(nocsf))
c
      IF(iexpc.NE.0 .OR. icitg.NE.0)THEN
c
c     set up pointers for expanding energy expressions
c
         CALL MKPT(kpt,nobl,nobc,symtyp,notgt,nctgt,mcont,gucont,mocsf,
     &             nocsf,mxcsf,ncont,ncont2,ntgsym,ntgt,numtgt,ncitot,
     &             nummx,iexpc,icitg,lusme)
         WRITE(nftw,1050)nocsf, mxcsf, mocsf
 1050    FORMAT(' SCATCI will expand NOCSF =',i7,' prototype CSFs',/15x,
     &          'into MXCSF =',i7,' actual configurations',/15x,
     &          'and  MOCSF =',i7,' dimension final Hamiltonian')
      ELSE
         mocsf=nocsf
         mxcsf=nocsf
         ntgt=ntgsym
         IF(LUSME.NE.0)WRITE(lusme,*)mocsf
      END IF
      if (iposit .eq. 0) then
        WRITE(nftw,1060)ncont
      else 
        write(nftw,1060)(mocsf-(nocsf-ncont))
      endif
 1060 FORMAT(/' Number of last continuum CSF, NCONT =',i7)
C
C---- Set up the integral blocks table: Two different cases here
C                 either C-inf-v or D2h
C
      IF(SYMTYP.LE.1)THEN
         CALL TABLBA(NFTW,NSYM,NMPRS,LPQRS,NPQRS,MPQRS,SCFUSE)
      ELSE
         WRITE(NFTW,801)SYMTYP
         IF(npflg(6).GT.0)i=MPROD(1,1,NPflg(6),NFTW)
 801     FORMAT(5X,
     &          'ENERGY EXPRESSIONS COMPUTED FOR MOLECULE SYMMETRY CASE'
     &          ,' - SYMTYP=',I3)
         CALL TABLBM(NFTW,NSYM,NPQRS,MPQRS,NPFLG(6))
      END IF
C
C---- Generate the pointer arrays for blocks of one and two electron
C     integrals based on the table just evaluated
C
      CALL TABLA(NFTW,NSYM,NOBE,NPQRS,MPQRS,MBAS,NBASH,NBAS,NINT2e,MXRW,
     &           SYMTYP,NPFLG(6),lusme)
      I = sum(NOBE(1:NSYM)) !total number of orbitals
      allocate(NRI(I),NSM(I))
      NRI(:) = 0
      NSM(:) = 0
      CALL NOBID(NOBT,NSYM,NOBE,NRI,NSM,nobmx,LUSME)
c
C-----------------------------------------------------------------------
C
C     NIPAIR DEFINES THE MAX NUMBER OF I*(I-1)/2 VALUES
C        REQUIRED FOR SUBROUTINE INDEX / SEE COMMENTS AT BEGINNING OF
C        SUBROUTINE TABLA FOR DETAILS
C
      IF(symtyp.LE.1)nobmx=0
      NIPAIR=Max((NSYM*NSYM+NSYM)/2+1,mxrw,300,nobmx)
      ALLOCATE(kpair(nipair+1))
C
      CALL GENIPR(NIPAIR,kpair)
C
C     Lay out space to store integrals and (portions of) the final
C     Hamiltonian
C
      nint1e=NBASH(NSYM+1)-1
      IF(npflg(6).GT.0 .AND. SYMTYP.GT.1)
     &   CALL PRBLKS(MBAS,nint1e,NINT2e,10,nsym,NFTW)
C
C-----------------------------------------------------------------------
      WRITE(6,*)'  Read wavefunctions from CONGEN and load integrals'
C-----------------------------------------------------------------------
      ALLOCATE(knorb(norb),kmn(nsrb),kmm(nsrb),kmg(nsrb),kms(nsrb),
     &         kmcon(nsrb),kmcorb(norb),kindo(nocsf+1),kicdo(nocsf+1),
     &         kndo(lndof),knsrb(2*nsrb),cdo(lcdof))
C-----------------------------------------------------------------------
c    names of variables:
c     in EXDRVF  |  in MKORBS  | dimension
c     ---------- | ----------- | ----------
c     NSYM       |  NSYM       | -
c     KMN        |  MN         | NSRB
c     KMG        |  MG         | NSRB
c     KMM        |  MM         | NSRB
c     KMS        |  MS         | NSRB
c     NORB       |  NORB       | -
c     NSRB       |  NSRBD      | -
c     KNORB      |  MAP        | NORB
c     KNSRB      |  MPOS       | NSRB
c     KMCON      |  MCON       | NSRB
c     KMCORB     |  MCORB      | NORB
c     IPOSIT     |  IPOSIT     | -
c     NOBL       |  NOBL       | 20
c     NOB0L      |  NOB0L      | 20
c     SYMTYP     |  SYMTYP     | -
C-----------------------------------------------------------------------
c
      CALL RDWF(NOCSF+1,KICDO,KINDO,KNDO,LNDOF,CDO,LCDOF,megul)
C
      IF((iposit.NE.0) .AND. (symtyp.EQ.2))THEN
         CALL PMKORBS(nob,nobe,nob0,nsym,kmn,kmg,kmm,kms,kmcon,kmcorb,
     &                NORB,NSRB,knorb,knsrb,iposit,symtyp,LUSME,NFTW)
      ELSE
         CALL MKORBS(NSYM,KMN,KMG,KMM,KMS,NORB,NSRB,KNORB,KNSRB,kmcon,
     &               kmcorb,IPOSIT,NOBL,NOB0L,SYMTYP,LUSME,NFTW)
      END IF
C
      IF(NALM.NE.0)GO TO 38
C
      IF(NPFLG(1).NE.0)THEN
         WRITE(NFTW,104)
 104     FORMAT(/' INPUT FUNCTIONS IN PACKED FORM')
         CALL PTPWF(NFTW,NOCSF,NELT,NDTRF,NODO,KINDO,KICDO,KNDO,CDO)
      END IF
C-----------------------------------------------------------------------
C     ZM: Memory allocation for the 1p and 2p integrals occurs only
C         for the SWEDEN/ALCHEMY case since the UKRMol+ integrals are kept
C         in the interface module ukrmolp_interface.
C
      ALLOCATE(eig(ntgt),vecs(ncitot))
C
C     SWEDEN/ALCHEMY integrals on input
C
      if (.not.(ukrmolp_ints)) then
         WRITE(NFTW,'(/,10X,"Assuming SWEDEN integrals on input.")')
         ALLOCATE(x1e(2*nint1e),x2e(nint2e),stat=i)
         if (i .ne. 0) stop "Memory allocation error for integrals."
      else
         WRITE(NFTW,'(/,10X,"Assuming UKRMol+ integrals on input.")')
         call mpi_mod_start
      endif
C
c     load integrals into core
c
ccc      write(*,*) 'call INTSIN: load integrals into core'
      CALL INTSIN(nfte,nfti,lembf,x1e,nint1e,x2e,nint2e,mocsf,nftw,
     &            symtyp,nsym,nobe,iposit,scalem,name,nalm,qmoln,
     &            ukrmolp_ints)
ccc      write(*,*) 'INTSIN finished'
      IF(nalm.NE.0)GO TO 38
C-----------------------------------------------------------------------
      CALL MAKEMG(KMG,NSRB,NELT,NDTRF)
c
      IF(icitg.NE.0)THEN
         WRITE(6,*)'  Load CI target and construct matrix elements'
c
c     Target CI vectors and energies are needed
c
         IF(nftg.GT.0)THEN
c           read them from file
            CALL CIRMAT(NTGsym,NFTG,NTGTF,vecs,eig,NTGTS,NCTGT,NFTW,
     &                  numtgt,kphz,nctarg,symtyp)
         ELSE
c           construct them
            jump=2
            istart=1
            kvecs=1
            keigen=1
            idiagt=0
            IF(idiag.EQ.3)idiagt=1
            DO i=1, ntgsym
               IF(iexpc.EQ.0)jump=notgt(i)
               CALL ENRGMT(NORB,KMN,KMM,KMS,NELT,NDTRF,KMG,nctgt(i),
     &                     NODO,KINDO,KNDO,KICDO,CDO,THRES,NALM,NFTW,
     &                     IPOSIT,KNORB,KNSRB,kmcon,kmcorb,NSRB,SYMTYP,
     &                     idiagt,nelm,nfte,lembf,x1e,nint1e,x2e,thrhm,
     &                     IODR,NRI,NSM,NOBE,MBAS,NBAS,NBASH,kpair,
     &                     istart,jump,ukrmolp_ints,NOEX)
               IF(nalm.NE.0)GO TO 38
               CALL dgtarg(nfte,nctgt(i),numtgt(i),eig(keigen),
     &                     vecs(kvecs),i,nelm,CRITE,CRITC,CRITR,ORTHO,
     &                     maxiter,igh,symtyp)
               istart=istart+nctgt(i)*jump
               kvecs=kvecs+nctgt(i)*numtgt(i)
               keigen=keigen+numtgt(i)
            END DO
         END IF

! Do loop and write statements below commented as they are not needed
! at the moment. Left in case they are neeed for future development (JDG)
!     Required by cdenprop to determine arbitrary overall phase 
!     between target CI vectors generated here and in the target run
!     AlexH 2011

         k=0     
         do i=1, ntgsym
            do m=1,numtgt(i)
!               write(2000,*) ' Target state ', i,m 
               do j=1,nctgt(i)
                  k=k+1
!                  write(2000,*) vecs(k)
               end do
!               write(2000,*) ' '
            end do
         end do
!         write(2000,*) ' The target energies'
!         write(2000,*) eig
         write(2001) vecs(1:k) !Target CI vectors

      END IF
c      DEALLOCATE(kphz)
c
c*********************************************************************
c  A  deallocate statment for kphz is introduced in EXDVRF.
c  (see statment above). This provides that the array kphz is
c  deallocated in EXDVRF also for the case (nftg .eq. 0). This array
c  was deallocated above in this subroutine only "if (nftg .gt. 0)"
c  after the "call CIRMAT" statment. (N. Vinci 21/11/01).
c
c  Moved deallocation of kphz as we need it to be saved to the 
c  target CI vectors for input to CDENPROP. (A. Harvey 15/5/14)
C
c*********************************************************************
C
c     Now construction of the matrix elements
c
      idiag=min(idiag,2)
      IF(icitg.NE.0 .OR. ntgsym.GT.1)THEN

c    Enhancement factor to be used for the scattering calculations when
c    projectile is a positron
      if (enh_factor_type .eq. 1) then
        write(6,'(/,"Using an enhancement factor of ",f10.7)')enh_factor
      elseif (enh_factor_type .eq. 2) then
        enh_factor = 1.0
        write(6,'(/,"Computing MP2 enhancement factor")')
      elseif (enh_factor_type .eq. 3) then
        write(6,'(/,"MP2 enhancementfactor with ", f10.7)')
     & enh_factor
      else
         enh_factor = 1.0
       endif

         CALL ENRGMX(NORB,KMN,KMM,KMS,NELT,NDTRF,KMG,NOCSF,NODO,KINDO,
     &               KNDO,KICDO,CDO,THRES,NALM,IDIAG,nelm,NFTW,IPOSIT,
     &               KNORB,KNSRB,kmcon,kmcorb,NSRB,SYMTYP,ncont,ncont2,
     &               nfte,lembf,x1e,nint1e,x2e,thrhm,IODR,NRI,NSM,NOBE,
     &               MBAS,NBAS,NBASH,kpair,NCORB,kpt,mocsf,iexpc,icitg,
     &               vecs,eig,numtgt,ntgsym,nctgt,notgt,nummx,gucont,
     &               mcont,lusme,enh_factor,enh_factor_type,
     &               ukrmolp_ints,NOEX)
      ELSE
         CALL ENRGMS(NORB,KMN,KMM,KMS,NELT,NDTRF,KMG,NOCSF,NODO,KINDO,
     &               KNDO,KICDO,CDO,THRES,NALM,IDIAG,nelm,NFTW,IPOSIT,
     &               KNORB,KNSRB,kmcon,kmcorb,NSRB,SYMTYP,ncont,ncont2,
     &               nfte,lembf,x1e,nint1e,x2e,thrhm,IODR,NRI,NSM,NOBE,
     &               MBAS,NBAS,NBASH,kpair,NCORB,kpt,mocsf,iexpc,eig,
     &               notgt(1),nummx,LUSME,enh_factor,enh_factor_type,
     &               ukrmolp_ints,NOEX)
      END IF
c
      if (.not.(ukrmolp_ints)) deallocate(x1e,x2e)
      DEALLOCATE(nodo,kpt,kpair,knorb,kmn,kmm,kmg,kms,kmcon,kmcorb,
     &           kindo,kicdo,kndo,knsrb,cdo,eig,vecs)
C
      IF(NALM.NE.0)GO TO 38
C
      IF(npflg(4).GT.0)CALL prtem(nfte,npflg(4)-1,lembf,nftw)
c
c     diagonalise the Hamiltonian if requested
c
 5000 IF(icidg.NE.0)THEN
         CALL dgmain(ntgcon,numtgt,notgt,mcont,nfte,nelm,symtyp,s,sz,
     &               nelt,mgvn,nfti,kphz,size(kphz),ukrmolp_ints)
C Z. Masin: if not diagonalize then write the params needed to be able to perform the diagonalization and writing of fort.25 later
      ELSE
         write(nftw,'(/,"Diagonalization not selected")')

         lcont = mocsf-(nocsf-ncont) !the sequence number of the last continuum CSF; we use it to write out to fort.25 only the CI coefficients for the continuum CSFs
         write(nftw,'("Sequence number of the last continuum CSF:",i0)')
     &         lcont
         size_kphz = size(kphz)
         write_ham_pars = .true.
         CALL ham_pars_io(write_ham_pars,lcont,ntgcon,numtgt,notgt,mcont
     &                    ,nfte,nelm,symtyp,s,sz,nelt,mgvn,nfti
     &                    ,kphz,size_kphz,ukrmolp_ints)

         write(nftw,'(/,"Parameters written to: ham_data")')
      END IF

      DEALLOCATE(kphz)

      if (ukrmolp_ints) call mpi_mod_finalize

      RETURN
C
 38   WRITE(NFTW,46)
 46   FORMAT('  DUE TO ALARM CONDITION THIS RUN WAS TERMINATED')
      RETURN
      END SUBROUTINE EXDRVF
C
      SUBROUTINE ham_pars_io(write_ham_pars,lcont,ntgcon,numtgt,notgt,
     &                       mcont,nfte,nelm,symtyp,s,sz,nelt,mgvn,nfti,
     &                       kphz,size_kphz,ukrmolp_ints)
C     Reads or writes the file 'ham_data' containing the parameters needed on the call to dgmain
C     if write_ham_pars .eq. .true. then the file is written interpretting the arguments as input
C     and if write_ham_pars .eq. .false. then the file is read in interpretting the argumens as output
      USE precisn, ONLY : wp
      USE scatci_data, ONLY : NTGTMX   
      IMPLICIT NONE
C     Argument input/output
      LOGICAL, INTENT(IN) :: write_ham_pars
      INTEGER, INTENT(INOUT) :: NTGCON, NFTE, NELM, SYMTYP, NELT, MGVN,
     &                       NFTI, LCONT, size_kphz
      INTEGER, DIMENSION(ntgtmx), INTENT(INOUT) :: MCONT, NOTGT, NUMTGT
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: kphz(:)
      LOGICAL, INTENT(INOUT) :: ukrmolp_ints
      REAL(KIND=wp), INTENT(INOUT) :: S, SZ
C     Local variables
      INTEGER :: ihd, err
C     
      ihd = 444
      if (write_ham_pars) then
         open(unit=ihd,file='ham_data',status='replace',
     &        form='unformatted',iostat=err)
      else
         open(unit=ihd,file='ham_data',status='old',
     &        form='unformatted',iostat=err)
      endif

      if (err .ne. 0) then
         stop "Error opening/creating the file ham_data."
      endif

      if (write_ham_pars) then
         write(*,'(/,"Writing ham_data...")')
         write(ihd)lcont,ntgcon,numtgt,notgt,mcont,nfte,nelm,symtyp,s,sz
     &             ,nelt,mgvn,nfti,size_kphz,ukrmolp_ints
         if (.not.(allocated(kphz))) then
            stop "On input the array kphz has not been allocated."
         else
            write(ihd) kphz(1:size_kphz)
         endif
         write(*,'("Done")')

      else
         write(*,'(/,"Reading ham_data...")')
         read(ihd) lcont,ntgcon,numtgt,notgt,mcont,nfte,nelm,symtyp,s,sz
     &             ,nelt,mgvn,nfti,size_kphz,ukrmolp_ints

         if (allocated(kphz)) deallocate(kphz)
         allocate(kphz(1:size_kphz),stat=err)
         if (err .ne. 0) stop "Memory allocation failed."
         kphz(1:size_kphz) = 0

         read(ihd) kphz(1:size_kphz)

         write(*,'("Done")')
      endif

      close(ihd)

      return

      END SUBROUTINE ham_pars_io
!*==expcor.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE EXPCOR(jorb,ncont,nti,xnip0,xnip,nd)
c
c     EXPCOR expands off-diagonal matrix elements between a
c     continuum orbital and an L**2 CSF.
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints, pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: JORB, NCONT, ND, NTI
      INTEGER(longint), DIMENSION(2,nd,ncont) :: XNIP
      INTEGER(longint), DIMENSION(2,nti) :: XNIP0
      INTENT (IN) JORB, NCONT, ND, NTI
C
C Local variables
C
      INTEGER :: I, IT, LW
      INTEGER, DIMENSION(8) :: LWD
C
C*** End of declarations rewritten by SPAG
C
c     Expand each integral in turn
c
      DO it=1, nti
         CALL unpack8ints(xnip0(1,it),lwd)
c     test for occurances of the continuum orbital in the integral label
         lw=0
         DO i=1, 4
            IF(lwd(i).EQ.jorb)lw=i
         END DO
c     No occurances, integral same for all matrix elements
         IF(lw.EQ.0)THEN
            DO i=2, ncont
               xnip(1,it,i)=xnip0(1,it)
               xnip(2,it,i)=xnip0(2,it)
            END DO
c     One occurance, update index and repack
         ELSE
            DO i=2, ncont
               lwd(lw)=jorb+i-1
               CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,0,
     &                        xnip(1,it,i))
            END DO
         END IF
      END DO
c
      RETURN
      END SUBROUTINE EXPCOR
!*==expdg.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE EXPDG(jorb,ncont,nti,xnip0,xnip,nd)
c
c     EXPDG
c     (a) expands diagonal matrix elements for a continuum
c         orbital for each target state.
c     (b) expands off-diagonal matrix elements for the case where the
c         continuum orbitals are actually the same
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints, pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: JORB, NCONT, ND, NTI
      INTEGER(longint), DIMENSION(2,nd,*) :: XNIP
      INTEGER(longint), DIMENSION(2,nti) :: XNIP0
      INTENT (IN) JORB, NCONT, ND, NTI
C
C Local variables
C
      INTEGER :: I, IT, LW1, LW2
      INTEGER, DIMENSION(8) :: LWD
C
C*** End of declarations rewritten by SPAG
C
c     Expand each integral in turn
c
      DO it=1, nti
         CALL unpack8ints(xnip0(1,it),lwd)
c     test for occurances of the continuum orbital in the integral label
         lw1=0
         lw2=0
         DO i=1, 4
            IF(lwd(i).EQ.jorb)THEN
               IF(lw1.EQ.0)THEN
                  lw1=i
               ELSE
                  lw2=i
               END IF
            END IF
         END DO
c     No occurances, integral same for all matrix elements
         IF(lw1.EQ.0)THEN
            DO i=2, ncont
               xnip(1,it,i)=xnip0(1,it)
               xnip(2,it,i)=xnip0(2,it)
            END DO
            CYCLE
c     One occurance, update index and repack
         ELSE IF(lw2.EQ.0)THEN
            DO i=2, ncont
               lwd(lw1)=jorb+i-1
               CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,0,
     &                        xnip(1,it,i))
            END DO
c     General case: two occurances, update index and repack
         ELSE
            DO i=2, ncont
               lwd(lw1)=jorb+i-1
               lwd(lw2)=jorb+i-1
               CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,0,
     &                        xnip(1,it,i))
            END DO
         END IF
      END DO
c
      RETURN
      END SUBROUTINE EXPDG
!*==expndg.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE EXPNDG(jorba,jorbb,nconta,ncontb,leng,nti,xnip0,xnip,
     &                  nd,ipair)
c
c     EXPNDG expands off-diagonal (upper & lower triangle) matrix
c     elements between like continuum orbitals but different target
c     configurations
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints, pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: JORBA, JORBB, LENG, NCONTA, NCONTB, ND, NTI
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER(longint), DIMENSION(2,nd,*) :: XNIP
      INTEGER(longint), DIMENSION(2,nti) :: XNIP0
      INTENT (IN) IPAIR, JORBA, JORBB, LENG, NCONTA, NCONTB, ND, NTI
C
C Local variables
C
      INTEGER :: I, IA, IB, II, IT, LWA, LWB
      INTEGER, DIMENSION(8) :: LWD
      INTEGER(longint), DIMENSION(2) :: XIWD
C
C*** End of declarations rewritten by SPAG
C
c     Expand each integral in turn
c
      DO it=1, nti
         CALL unpack8ints(xnip0(1,it),lwd)
         IF(lwd(1).LT.lwd(2) .OR. lwd(3).LT.lwd(4) .OR. ipair(lwd(1))
     &      +lwd(2).LT.ipair(lwd(3))+lwd(4))WRITE(6,900)lwd
 900     FORMAT(' Unexpected orbital ordering in EXPNDG input:',
     &          /' LWD =',4I4)
c     test for occurances of the continuum orbital in the integral label
         lwa=0
         lwb=0
         DO i=1, 4
            IF(lwd(i).EQ.jorba)lwa=i
            IF(lwd(i).EQ.jorbb)lwb=i
         END DO
c     No occurances, integral same for all matrix elements
         IF(max(lwa,lwb).EQ.0)THEN
            DO i=2, leng
               xnip(1,it,i)=xnip0(1,it)
               xnip(2,it,i)=xnip0(2,it)
            END DO
            CYCLE
c     One occurance, update index and repack
         ELSE IF(lwb.EQ.0)THEN
            i=1
            ii=3
            DO ia=1, nconta
               lwd(lwa)=jorba+ia-1
               CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,0,
     &                        xiwd)
               DO ib=ii, ncontb
                  IF(ia.NE.ib)THEN
                     i=i+1
                     xnip(1,it,i)=xiwd(1)
                     xnip(2,it,i)=xiwd(2)
                  END IF
               END DO
               ii=1
            END DO
c     One occurance, update index and repack
         ELSE IF(lwa.EQ.0)THEN
            i=1
            ii=3
            DO ia=1, nconta
               DO ib=ii, ncontb
                  IF(ia.NE.ib)THEN
                     i=i+1
                     lwd(lwb)=jorbb+ib-2
                     CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),
     &                              0,0,0,xnip(1,it,i))
                  END IF
               END DO
               ii=1
            END DO
c     General case: two occurances, update index and repack
         ELSE
            i=1
            ii=3
            DO ia=1, nconta
               lwd(lwa)=jorba+ia-1
               DO ib=ii, ncontb
                  IF(ia.NE.ib)THEN
                     i=i+1
                     lwd(lwb)=jorbb+ib-2
                     IF(ipair(lwd(1))+lwd(2).GT.ipair(lwd(3))+lwd(4))
     &                  THEN
                        CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),
     &                                 lwd(5),0,0,0,xnip(1,it,i))
                     ELSE
                        CALL pack8ints(lwd(3),lwd(4),lwd(1),lwd(2),
     &                                 lwd(5),0,0,0,xnip(1,it,i))
                     END IF
                  END IF
               END DO
               ii=1
            END DO
         END IF
      END DO
c
      RETURN
      END SUBROUTINE EXPNDG
!*==expodl.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE EXPODL(jorba,jorbb,ncont,leng,nti,xnip0,xnip,nd)
c
c     EXPODL expands off-diagonal lower triangle matrix elements between
c     like continuum orbitals
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints, pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: JORBA, JORBB, LENG, NCONT, ND, NTI
      INTEGER(longint), DIMENSION(2,nd,*) :: XNIP
      INTEGER(longint), DIMENSION(2,nti) :: XNIP0
      INTENT (IN) JORBA, JORBB, LENG, NCONT, ND, NTI
C
C Local variables
C
      INTEGER :: I, IA, IB, IT, LWA, LWB
      INTEGER, DIMENSION(8) :: LWD
      INTEGER(longint), DIMENSION(2) :: XIWD
C
C*** End of declarations rewritten by SPAG
C
c     Expand each integral in turn
c
      DO it=1, nti
         CALL unpack8ints(xnip0(1,it),lwd)
c     test for occurences of the continuum orbital in the integral label
         lwa=0
         lwb=0
         DO i=1, 4
            IF(lwd(i).EQ.jorba)lwa=i
            IF(lwd(i).EQ.jorbb)lwb=i
         END DO
c     No occurances, integral same for all matrix elements
         IF(max(lwa,lwb).EQ.0)THEN
            DO i=2, leng
               xnip(1,it,i)=xnip0(1,it)
               xnip(2,it,i)=xnip0(2,it)
            END DO
            CYCLE
c     One occurance, update index and repack
         ELSE IF(lwb.EQ.0)THEN
            i=1
            DO ia=1, ncont-1
               lwd(lwa)=jorba+ia-1
               CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,0,
     &                        xiwd)
               DO ib=max(3,ia+1), ncont
                  i=i+1
                  xnip(1,it,i)=xiwd(1)
                  xnip(2,it,i)=xiwd(2)
               END DO
            END DO
c     One occurance, update index and repack
         ELSE IF(lwa.EQ.0)THEN
            i=1
            DO ia=1, ncont-1
               DO ib=max(3,ia+1), ncont
                  i=i+1
                  lwd(lwb)=jorbb+ib-2
                  CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,
     &                           0,xnip(1,it,i))
               END DO
            END DO
c     General case: two occurances, update index and repack
         ELSE
            i=1
            DO ia=1, ncont-1
               lwd(lwa)=jorba+ia-1
               DO ib=max(3,ia+1), ncont
                  i=i+1
                  lwd(lwb)=jorbb+ib-2
                  CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,
     &                           0,xnip(1,it,i))
               END DO
            END DO
         END IF
      END DO
c
      RETURN
      END SUBROUTINE EXPODL
!*==expoff.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE EXPOFF(jorba,jorbb,nconta,ncontb,leng,nti,xnip0,xnip,
     &                  nd)
c
c     EXPOFF expands a general off-diagonal matrix elements between the
c     continuum orbitals for different continua and different target sta
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints, pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: JORBA, JORBB, LENG, NCONTA, NCONTB, ND, NTI
      INTEGER(longint), DIMENSION(2,nd,*) :: XNIP
      INTEGER(longint), DIMENSION(2,nti) :: XNIP0
      INTENT (IN) JORBA, JORBB, LENG, NCONTA, NCONTB, ND, NTI
C
C Local variables
C
      INTEGER :: I, IA, IB, II, IT, LWA, LWB
      INTEGER, DIMENSION(8) :: LWD
      INTEGER(longint), DIMENSION(2) :: XIWD
C
C*** End of declarations rewritten by SPAG
C
c     Expand each integral in turn
c
      DO it=1, nti
         CALL unpack8ints(xnip0(1,it),lwd)
c     test for occurances of the continuum orbital in the integral label
         lwa=0
         lwb=0
         DO i=1, 4
            IF(lwd(i).EQ.jorba)lwa=i
            IF(lwd(i).EQ.jorbb)lwb=i
         END DO
c     No occurances, integral same for all matrix elements
         IF(max(lwa,lwb).EQ.0)THEN
            DO i=2, leng
               xnip(1,it,i)=xnip0(1,it)
               xnip(2,it,i)=xnip0(2,it)
            END DO
            CYCLE
c     One occurance, update index and repack
         ELSE IF(lwa.EQ.0)THEN
            i=1
            ii=2
            DO ia=1, nconta
               DO ib=ii, ncontb
                  lwd(lwb)=jorbb+ib-1
                  i=i+1
                  CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,
     &                           0,xnip(1,it,i))
               END DO
               ii=1
            END DO
c     One occurance, update index and repack
         ELSE IF(lwb.EQ.0)THEN
            i=1
            ii=2
            DO ia=1, nconta
               lwd(lwa)=jorba+ia-1
               CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,0,
     &                        xiwd)
               DO ib=ii, ncontb
                  i=i+1
                  xnip(1,it,i)=xiwd(1)
                  xnip(2,it,i)=xiwd(2)
               END DO
               ii=1
            END DO
c     General case: two occurances, update index and repack
         ELSE
            i=1
            ii=2
            DO ia=1, nconta
               lwd(lwa)=jorba+ia-1
               DO ib=ii, ncontb
                  lwd(lwb)=jorbb+ib-1
                  i=i+1
                  CALL pack8ints(lwd(1),lwd(2),lwd(3),lwd(4),lwd(5),0,0,
     &                           0,xnip(1,it,i))
               END DO
               ii=1
            END DO
         END IF
      END DO
c
      RETURN
      END SUBROUTINE EXPOFF
!*==genipr.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE GENIPR(N,MX)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: N
      INTEGER, DIMENSION(0:n) :: MX
      INTENT (IN) N
      INTENT (OUT) MX
C
C Local variables
C
      INTEGER :: IT, L
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     Generate the IPAIR table which is for the sequential index
C     along the diagonal of a lower half triangle:
C
C              |  0   \
C              |  1    \
C              |  2 3   \
C              |  4 5 6  \
C              |  7 8 9 10\
C              |   etc...  \
C              |   etc...  \
C              |            \
C              --------------
C     Thus,
C
C        MX = 0, 1, 3, 6, 10, ......
C
C***********************************************************************
C
      IT=0
      mx(0)=0
      DO L=1, N
         MX(L)=IT
         IT=IT+L
      END DO
C
      RETURN
      END SUBROUTINE GENIPR
!*==ijkpqrs.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      FUNCTION IJKPQRS(NA,MJK,NMAX)
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NMAX
      INTEGER :: IJKPQRS
      INTEGER(longint), DIMENSION(2,*) :: MJK
      INTEGER(longint), DIMENSION(2) :: NA
      INTENT (IN) MJK, NA, NMAX
C
C Local variables
C
      INTEGER :: N
C
C*** End of declarations rewritten by SPAG
C
c     check list of integral indices MJK for existence of new label NA
C
      DO N=1, NMAX
         IF(MJK(1,N).EQ.NA(1) .AND. MJK(2,N).EQ.NA(2))THEN
            ijkpqrs=n
            RETURN
         END IF
      END DO
C
      ijkpqrs=0
      RETURN
      END FUNCTION IJKPQRS
!*==intsin.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE INTSIN(nfte,nfti,lembf,xint1e,nint1e,xint2e,nint2e,
     &                  nocsf,nfta,isymtp,nsym1,nob1,iposit,scalem,name,
     &                  nalm,qmoln,ukrmolp_ints)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : ONE=>XONE, ZERO=>XZERO
      USE params, ONLY : ctrans1
      USE global_utils, ONLY : INTAPE
! the data_file_obj_id constant is used to identify the format of the UKRMol+ integral file
      USE ukrmol_interface_gbl, ONLY : read_ukrmolp_ints,
     &                                 data_file_obj_id, line_len,
     &                                 get_kinetic_energy_integral
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IPOSIT, ISYMTP, LEMBF, NALM, NFTA, NFTE, NFTI, NINT1E, 
     &           NINT2E, NOCSF, NSYM1
      CHARACTER(LEN=120) :: NAME
      REAL(KIND=wp) :: SCALEM
      INTEGER, DIMENSION(nsym1) :: NOB1
      REAL(KIND=wp), ALLOCATABLE :: XINT1E(:) !DIMENSION(2*nint1e)
      REAL(KIND=wp), ALLOCATABLE :: XINT2E(:)
      INTENT (IN) IPOSIT, ISYMTP, LEMBF, NAME, NFTE, NOB1, NOCSF, NSYM1
      INTENT (INOUT) NALM, XINT1E
      LOGICAL :: QMOLN
      INTENT (IN) QMOLN
      LOGICAL :: ukrmolp_ints
      INTENT (IN) ukrmolp_ints
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(20) :: CHARG, XNUC
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp) :: EN, SIGN, VSMALL, X1E
      INTEGER :: I, IA, IAT, IFAIL, IMAX, IMAXP, IMIN, IND, IONEIN, K, 
     &           L, LTRI, N, ND, NEND, NNUC, NSYM, NT
c     Alchemy header arrays
      CHARACTER(LEN=4), DIMENSION(4) :: LABEL
      CHARACTER(LEN=4), DIMENSION(33) :: NAM1
c     Sweden header arrays
      INTEGER, DIMENSION(8) :: NAO, NCORE
      INTEGER, DIMENSION(20) :: NHE, NOB
      REAL(KIND=wp), DIMENSION(NINT1E) :: XTEMPE, XTEMPP
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: XINTP
c     UKRMol+ header
      CHARACTER(LEN=line_len) :: ukrmolp_header
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     Taken from the Alchemy CI code. Reads the file of transformed
C     integrals from unit nfti and opens the file nfte used to store
C     the Hamiltonian matrix.
C
C     Note: This version is interface to Sweden integrals and can be
C           used for polyatomic calculations too.
C
C       ZM: Extended to work with UKRMol+ integrals. If the logical UKRMOLP_INTS
C           is .true. then the polyatomic integrals on input are assumed to be from UKRMol+. Otherwise from SWEDEN.
C
C***********************************************************************
C
      EQUIVALENCE(NHE(1),NOB(1))
      EQUIVALENCE(DTNUC(1),EN)
      EQUIVALENCE(DTNUC(22),XNUC(1))
      EQUIVALENCE(DTNUC(2),CHARG(1))
      DATA vsmall/1.E-10_wp/, nnuc/0/
      NINT1e=0
      NINT2e=0
C
C
C---- Read the header on the integrals written by the four index
C     Transformation code.
C
C               Either Alchemy or Sweden/UKRMol+ format
C                (ISYMTP.LT.2)    (ISYMTP.ge.2)
C
C     Note that Alchemy/UKRMol+ headers give info for computing the nuclear
C     interaction potential while Sweden gives the actual value.
C
      IF(ISYMTP.LT.2)THEN
         READ(NFTI)(NAM1(I),I=1,33), NSYM, NT, NNUC, ND, LTRI, 
     &             (NOB(I),I=1,NSYM), (ND,I=1,NT), (ND,I=1,NT), 
     &             (CHARG(I),I=1,NNUC), (XNUC(I),I=1,NNUC)
         WRITE(NFTA,1240)(NAM1(I),I=1,30), NOCSF, LEMBF, ISYMTP
         EN=ZERO
         DO N=2, NNUC
            IF(abs(charg(n)).LT.vsmall)CYCLE
            DO I=1, N-1
               IF(abs(charg(i)).LT.vsmall)CYCLE
               EN=EN+CHARG(N)*CHARG(I)/ABS(XNUC(N)-XNUC(I))
            END DO
         END DO
C
      ELSE
C
C        ZM:
C        Here we identify the format of the polyatomic transformed integrals file.
C
C        In case UKRMol+ integrals file is on input we transfer control to the routine READ_UKRMOLP_INTS which reads in the integrals
C        and performs operations equivalent to those performed in the SWEDISH branch of this routine.
C
C        If the integrals are SWEDEN format then we continue this routine without any modifications.
C
         if (ukrmolp_ints) then
C
            call READ_UKRMOLP_INTS(nfte,nfti,lembf,nint1e,nint2e,
     &              nocsf,nfta,isymtp,nsym1,nob1,iposit,scalem,name,
     &              nalm,qmoln)
C BC - write the kinetic energy integrals to fort file 413, same as Hemal's solution
C line 4798 onwards (requires updates from uk_rmol_interface.f90)
            if (qmoln) then
               ind=0
               DO k=1, nsym1
                  DO l=1, nob1(k)
                     ind=ind+1
                     WRITE(413,12) k-1, l, l,
     &                     get_kinetic_energy_integral(ind,ind)
                  END DO
               END DO
            end if
C End of kinetic energy integrals to file.
            RETURN
C
         ELSE
C
            if (size(xint1e) < 2*nint1e) then
               print *,size(xint1e),2*nint1e
               stop "The 1p integrals buffer is too small."
            endif
            if (size(xint2e) < nint2e) then
               print *,size(xint2e),nint2e
               stop "The 2p integrals buffer is too small."
            endif
C
         ENDIF
C
         READ(NFTI)(LABEL(I),I=1,4)
         WRITE(6,FMT='(/,5X,''Label on Sweden tape = '',4A)')
     &         (LABEL(I),I=1,4)
         READ(NFTI)NSYM, EN, (NAO(I),NOB(I),NCORE(I),I=1,NSYM)
C
         DO I=1, NSYM
            IF(NCORE(I).NE.0)THEN
               WRITE(NFTA,1591)I, NCORE(I)
               WRITE(NFTA,1592)EN
            END IF
         END DO
      END IF
C
C---- Write to printout
C
      WRITE(NFTA,1241)NSYM, (NOB(I),I=1,NSYM)
      WRITE(NFTA,1510)EN
      IF(ISYMTP.LT.2)THEN
         WRITE(NFTA,1242)NNUC, (CHARG(I),I=1,NNUC)
         WRITE(NFTA,1243)(XNUC(I),I=1,NNUC)
      ELSE
         WRITE(NFTA,1520)(NAO(I),NOB(I),NCORE(I),I=1,NSYM)
      END IF
 1240 FORMAT(/' Transformed integrals read:',/5X,30A4/' NOCSF=',I5,3X,
     &       3X,'LEMBF=',I5,3X,'SYMTYP=',I3)
 1241 FORMAT(' NSYM =',I5,3X,'NOB  =',20I5)
 1510 FORMAT(/,10x,'Nuclear repulsion energy = ',f15.7)
C
C---- Following are used only by the Alchemy (diatomic) option
C
 1242 FORMAT(' NNUC =',I5,3X,'CHARG=',10F10.0/21X,10F10.0)
 1243 FORMAT(15X,'XNUC =',10F10.5/21X,10F10.5)
C
C---- Following are used only by the Sweden (polyatomic) option
C
 1520 FORMAT(/10x,'  NAO  NMO  NCORE',/,(10x,2I5,i7))
 1591 FORMAT(10X,'Non-zero CORE in Sweden. Sym no. = ',i5,' Ncore = ',
     &       I5,/)
 1592 FORMAT(10X,'Core Energy = ',F15.7,' Hartrees ',/)
c
c     check for data input errors
c
      IF(NSYM.NE.NSYM1)THEN
         IF(NSYM.GT.NSYM1)THEN
            DO I=NSYM1+1, NSYM
               IF(NOB(I).NE.0)GO TO 450
            END DO
         ELSE
            DO I=NSYM+1, NSYM1
               IF(NOB1(I).NE.0)GO TO 450
            END DO
         END IF
      END IF
      nsym=min(nsym,nsym1)
C
      DO I=1, NSYM
         IF(NOB(I).NE.NOB1(I))GO TO 450
      END DO
C
C---- Write header on energy matrix file
C
      OPEN(UNIT=nfte,FORM='unformatted')
      WRITE(NFTE)NOCSF, lembf, 0, nocsf, 0, nsym, 0, 0, 0, 0, nnuc, 0, 
     &           NAME, NHE, DTNUC
c
      IF(isymtp.LE.1)THEN
c
c     Read Alchemy integrals
c
         ALLOCATE(xintp(nint1e))
         IA=1
         CALL RDINTS(NFTI,NSYM,NOB,LTRI,nint1e+1,IA,XINT1e,NALM)
         IF(NALM.NE.0)GO TO 901
         IA=1
         CALL RDINTS(NFTI,NSYM,NOB,LTRI,nint1e+1,IA,XINT1e,NALM)
         IF(NALM.NE.0)GO TO 901
         IAT=IA-1
         ia=1
         CALL RDINTS(NFTI,NSYM,NOB,LTRI,nint1e+1,IA,XINTp,NALM)
         IF(NALM.NE.0)GO TO 21
c
         IF(iat.NE.nint1e)GO TO 900
C
C---- ADD TOGETHER one electron CONTRIBUTIONS:
C     first term is K; second is V.
C     So scale K by mass of exotic particle and change sign for positron
         SIGN=ONE
         IF(ABS(IPOSIT).EQ.1)SIGN=-ONE
         DO I=1, IAT
            x1e=XINT1e(I)+XINTp(I)
CGDSTART: exotic particle
c        XINTp(I) = XINT1e(I)/SCALEM+SIGN*XINTp(I)
CGDEND:    electron
            XINT1e(I)=x1e
         END DO
C
 21      IMAX=0
         IMIN=0
         imaxp=0
         nend=-1
C
         ia=1
         CALL RDINT(IMIN,IMAX,NEND,imaxp,NFTI,IA,nint2e+1,XINT2e)
         IF(NEND.NE.1)GO TO 900
         DEALLOCATE(xintp)
      ELSE
C
C---- Read transformed Sweden integrals
C
C---- Number of one electron integrals per record now calculated.
C     N.B. All symmetries for each type on one record
C
         IONEIN=0
         DO I=1, NSYM
            IONEIN=IONEIN+(NOB(I)*(NOB(I)+1))/2
         END DO
C
C---- Pick up the one electron integrals for Hamiltonian. Currently
C     these are combined into the second record on the output file.
C     This needs changed for positrons as was done for Alchemy into
C     k.e. and nuclear integrals.
C
C     First locate header then skip the overlap integrals
C
         CALL SEARCH(NFTI,ctrans1,ifail)
         READ(NFTI)
         CALL INTAPE(NFTI,XINT1e,nint1e)
c        write(6,*) 'one electron integrals first call'
c        do 448 k=1,nint1e
c 448    write(6,*) k,xint1e(k)
C
C...START IMPLEMENTATION OF  ONE ELECTRON INTEGRALS (T-V) FOR POSITRON CASE
C
         IF(IPOSIT.EQ.0) GO TO 1918
 
         WRITE(126,*)'k+v integrals xtempe'
         DO I=1, NINT1E
            XTEMPE(I)=XINT1E(I)
         END DO
         CALL INTAPE(NFTI,XTEMPP,nint1e)
! Hemal Varambhia write statement 7th November 2008
         IF (QMOLN) THEN
            ind=0
            DO k=1, nsym
               DO l=1, nob(k)
                  ind=ind+l
                  WRITE(413,12)k-1, l, l, XTEMPP(ind)
               END DO
            END DO
         ENDIF
! End of Hemal Varambhia write statement 7th November 2008
         IF(iposit.NE.0)THEN
            WRITE(6,*)'k-v integrals xtempe'
            DO K=1, NINT1E
               XINT1E(K+NINT1E)=2.0_wp*XTEMPP(K)-XTEMPE(K)
            END DO
         END IF
               ! if statement added by Hemal Varambhia on 7th November 2008
 1918    CONTINUE
 
C
C..END IMPLEMENTATION OF ONE ELECTRON INTEGRALS (T-V) FOR POSITRON CASE
C
C     Read Sweden 2-electron integrals
C
         CALL CHN2E(NFTI,nfta,XINT2e,nint2e)
      END IF
c
c       write(6,*) 'one electron integrals after search are'
c       do 444 k=1,2*nint1e
c 444   write(6,*) K, xint1e(K)
 
c       write(6,*) 'two electron integrals are'
c       DO 445 k=1,nint2e
c 
c       write(6,*) k,xint2e(k)
c 445   continue
      WRITE(nfta,1600)nint1e, nint2e
 1600 FORMAT(' Integrals read successfully:',
     &       /' 1-electron integrals, NINT1e =',i10,
     &       /' 2-electron integrals, NINT2e =',i10)
      RETURN
c     error conditions
 450  NALM=1
      WRITE(NFTA,451)
      WRITE(NFTA,452)(NOB(I),I=1,NSYM)
      WRITE(NFTA,452)(NOB1(I),I=1,NSYM1)
 451  FORMAT(/' MISMATCH IN NOB BETWEEN INTEGRALS AND FORMULAE')
 452  FORMAT(' NOB=',20I5)
      RETURN
C
 900  NALM=1
 901  WRITE(NFTA,902)nalm, nend, iat, nint1e, nint2e
 902  FORMAT('  UNABLE TO READ INTEGRALS,  NALM =',i2,' NEND =',i3,
     &       /'  IAT =',i7,' NINT1e =',i7,' NINT2e =',i9)
 12   FORMAT(i6,3x,i6,3x,i6,5x,f12.7)
      RETURN
      END SUBROUTINE INTSIN
!*==makemg.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
 
      SUBROUTINE MAKEMG(MG,NSRB,NELT,NDTRF)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NELT, NSRB
      INTEGER, DIMENSION(*) :: MG
      INTEGER, DIMENSION(NELT) :: NDTRF
      INTENT (IN) NDTRF, NELT, NSRB
      INTENT (OUT) MG
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
      DO I=1, NSRB
         MG(I)=0
      END DO
      DO I=1, NELT
         MG(NDTRF(I))=I
      END DO
      RETURN
      END SUBROUTINE MAKEMG
!*==mkcfd.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE mkcfd(coef,ctgt,num,nctgt)
c
c     MKCFD makes a list of the target CI coefficients needed for
c     weighting the current integrals for a diagonal element
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NCTGT, NUM
      REAL(KIND=wp), DIMENSION(*) :: COEF
      REAL(KIND=wp), DIMENSION(nctgt,*) :: CTGT
      INTENT (IN) CTGT, NCTGT, NUM
      INTENT (OUT) COEF
C
C Local variables
C
      INTEGER :: I, IJ, J
C
C*** End of declarations rewritten by SPAG
C
      ij=0
      DO i=1, num
         DO j=i, num
            ij=ij+1
            coef(ij)=ctgt(1,i)*ctgt(1,j)
         END DO
      END DO
      RETURN
      END SUBROUTINE MKCFD
!*==mkcfg.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE mkcfg(coef,ctgt1,ictgt1,ctgt2,ictgt2,num1,num2,nctgt1,
     &                 nctgt2)
c
c     MKCFG makes a list of the target CI coefficients needed for
c     weight the current integrals: general off-diagonal case
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ICTGT1, ICTGT2, NCTGT1, NCTGT2, NUM1, NUM2
      REAL(KIND=wp), DIMENSION(*) :: COEF
      REAL(KIND=wp), DIMENSION(nctgt1,num1) :: CTGT1
      REAL(KIND=wp), DIMENSION(nctgt2,num2) :: CTGT2
      INTENT (IN) CTGT1, CTGT2, ICTGT1, ICTGT2, NCTGT1, NCTGT2, NUM1, 
     &            NUM2
      INTENT (OUT) COEF
C
C Local variables
C
      INTEGER :: I, IJ, J
C
C*** End of declarations rewritten by SPAG
C
      ij=0
      DO i=1, num1
         DO j=1, num2
            ij=ij+1
            coef(ij)=ctgt1(ictgt1,i)*ctgt2(ictgt2,j)
         END DO
      END DO
      RETURN
      END SUBROUTINE MKCFG
!*==mkcfl.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE mkcfl(coef,ctgt,num,nctgt)
c
c     MKCFL makes a list of the target CI coefficients needed for
c     weighting thecurrent integrals for a continuum -- L**2
c     off-diagonal matrix element
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NCTGT, NUM
      REAL(KIND=wp), DIMENSION(num) :: COEF
      REAL(KIND=wp), DIMENSION(nctgt,*) :: CTGT
      INTENT (IN) CTGT, NCTGT, NUM
      INTENT (OUT) COEF
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
      DO i=1, num
         coef(i)=ctgt(1,i)
      END DO
      RETURN
      END SUBROUTINE MKCFL
!*==mkcfo.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE mkcfo(coef,ctgt1,ctgt2,num,nctgt)
c
c     MKCFO makes a list of the target CI coefficients needed for
C     weighting the current integrals for an off diagonal element
c     within one target state.
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NCTGT, NUM
      REAL(KIND=wp), DIMENSION(*) :: COEF
      REAL(KIND=wp), DIMENSION(nctgt,*) :: CTGT1, CTGT2
      INTENT (IN) CTGT1, CTGT2, NCTGT, NUM
      INTENT (OUT) COEF
C
C Local variables
C
      INTEGER :: I, IJ, J
C
C*** End of declarations rewritten by SPAG
C
      ij=0
      DO i=1, num
         DO j=i, num
            ij=ij+1
            coef(ij)=(ctgt1(1,i)*ctgt2(1,j)+ctgt2(1,i)*ctgt1(1,j))
         END DO
      END DO
      RETURN
      END SUBROUTINE MKCFO
C
      SUBROUTINE MKCFT(COEF,CTGT1,CTGT2,NUM,NCTGT)
C
C     MKCFT is like MKCFO: it calculates the lower triangle of
C     the matrix of products of target CI expansion coefficients.
C     The difference with repect to MKCFO is that here no
C     explicit symmetrization is done.
C
      USE precisn, ONLY: wp
      IMPLICIT NONE
C
C Dummy arguments
C
      INTEGER :: NCTGT, NUM
      REAL(KIND=wp), DIMENSION(*) :: COEF
      REAL(KIND=wp), DIMENSION(NCTGT,*) :: CTGT1, CTGT2
      INTENT (IN) CTGT1, CTGT2, NCTGT, NUM
      INTENT (OUT) COEF
C
C Local variables
C
      INTEGER :: I, IJ, J
C
      IJ=0
      DO I=1, NUM
         DO J=I, NUM
            IJ=IJ+1
            COEF(IJ)=CTGT1(1,I)*CTGT2(1,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE MKCFT
!*==mkorbs.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE MKORBS(NSYM,MN,MG,MM,MS,NORB,NSRBD,MAP,MPOS,mcon,mcorb,
     &                  IPOSIT,NOBL,NOB0L,SYMTYP,LUSME,NFTW)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IPOSIT, LUSME, NORB, NSRBD, NSYM, SYMTYP, NFTW
      INTEGER, DIMENSION(NORB) :: MAP
      INTEGER, DIMENSION(nsrbd) :: MCON
      INTEGER, DIMENSION(*) :: MCORB, MG, MM, MN, MS
      INTEGER, DIMENSION(NSRBD) :: MPOS
      INTEGER, DIMENSION(20) :: NOB0L, NOBL
      INTENT (IN) IPOSIT, LUSME, NOB0L, NOBL, NORB, NSRBD, NSYM, SYMTYP
      INTENT (IN) NFTW
      INTENT (OUT) MCORB
      INTENT (INOUT) MAP, MCON, MG, MM, MN, MPOS, MS
C
C Local variables
C
      INTEGER, SAVE :: I, IC, ICON, IK, IKP, IPOS, IS, ISO, J, K, M, M1, 
     &                 MA, MB, N, NEP, NSRB
      INTEGER, DIMENSION(20), SAVE :: NOBLJ
C
C*** End of declarations rewritten by SPAG
C
c........................................
      DO is=1, nsym
         NOBLJ(is)=NOBL(is)
      END DO
 
      IF(IPOSIT.NE.0)THEN
         DO is=1, nsym
            NOBLJ(is)=NOBL(is)/2
         END DO
      END IF
 
      IF(SYMTYP.EQ.0)THEN
         IC=1
         ISO=4
      ELSE IF(SYMTYP.EQ.1)THEN
         IC=2
         ISO=4
      ELSE
         IC=1
         ISO=2
      END IF
      I=1
      MA=0
      DO J=1, IC
         M1=MA+1
         MA=MA+NOBLJ(J)
         IPOS=0
         DO N=M1, MA
            MAP(N)=N
 
            IF(n.LT.nob0l(j)+m1)THEN
               icon=1
            ELSE
               icon=0
            END IF
            mcorb(n)=icon
            MN(I)=N
            MG(I)=0
            MM(I)=0
            MS(I)=0
            MPOS(I)=IPOS
            mcon(i)=icon
            I=I+1
            MN(I)=N
            MG(I)=0
            MM(I)=0
            MS(I)=1
            MPOS(I)=IPOS
            mcon(i)=icon
            I=I+1
         END DO
      END DO
      K=MA+1
C
      DO M=IC+1, NSYM*IC
         MA=NOBLJ(M)
         MB=(M-1)/IC
         IPOS=0
         DO N=1, MA
            MAP(K)=K
 
            IF(n.LE.nob0l(m))THEN
               icon=1
            ELSE
               icon=0
            END IF
            mcorb(k)=icon
            DO J=1, ISO
               MN(I)=K
               MG(I)=0
               MM(I)=MB
               MS(I)=0
               MPOS(I)=IPOS
               mcon(i)=icon
               I=I+1
            END DO
            K=K+1
            MS(I-1)=1
            IF(SYMTYP.LE.1)THEN
               MM(I-1)=-MB
               MM(I-2)=-MB
               MS(I-3)=1
            END IF
         END DO
      END DO
      NSRB=I-1
 
      IF(IPOSIT.EQ.0)GO TO 955
C
C...................POSITRON LOOP BEGINS.......................
C
      WRITE(NFTW,*)' MAP : OLD, NEW'
      IK=K-1
      IKP=IK+1
C     NEP = total number of orbitals (electron + positron)
      NEP=IK+IK
      DO N=IKP, NEP
         MAP(N)=MAP(N-IK)
      END DO
      WRITE(NFTW,*)' MAP : OLD, NEW'
      DO N=1, NEP
         WRITE(NFTW,*)N, MAP(N)
      END DO
      IF(LUSME.NE.0)THEN
         WRITE(LUSME,*)NEP
         WRITE(LUSME,*)map
      END IF
      DO N=IKP, NEP
         MN(IK+IKP)=N
         IK=IK+1
         MN(IK+IKP)=N
         IK=IK+1
      END DO
      DO J=1, NEP
         MG(J+NEP)=MG(J)
         MM(J+NEP)=MM(J)
         MS(J+NEP)=MS(J)
         MPOS(J+NEP)=1
 
      END DO
      NSRB=2*NSRB
C
C...................POSITRON LOOP ENDS.......................
C
 955  CONTINUE
 
      WRITE(NFTW,'(19x,A,9x,A,9x,A,9x,A,9x,A,6x,A,6x,A)')
     & 'I', 'N', 'G', 'M', 'S', 'MPOS', 'MCON'
      DO J=1, NSRB
         WRITE(NFTW,3000) J,MN(J),MG(J),MM(J),MS(J),MPOS(J),MCON(J)
      END DO
 3000 FORMAT(10X,7I10)
      WRITE(NFTW,*)'GIVEN NSRBD=', NSRBD
      WRITE(NFTW,*)'CALCULATED NSRB=', NSRB
      IF(nsrb.NE.nsrbd)THEN
         WRITE(NFTW,960)nsrbd, nsrb
 960     FORMAT(' Inconsistency detected in subroutine MKORBS:',
     &          /' Given NSRBD =',i5,' calculated NSRB =',i5,' STOP')
         STOP
      END IF
      RETURN
      END SUBROUTINE MKORBS
!*==mkpt.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE mkpt(kpt,nob,nobc,symtyp,notgt,nctgt,mcont,gucont,
     &                mocsf,nocsf0,mxcsf,ncont,ncont2,ntgsym,ntgt,
     &                numtgt,ncitot,nummx,iexpc,icitg,lusme)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ICITG, IEXPC, LUSME, MOCSF, MXCSF, NCITOT, NCONT, 
     &           NCONT2, NOCSF0, NTGSYM, NTGT, NUMMX, SYMTYP
      INTEGER, DIMENSION(ntgsym) :: GUCONT, MCONT, NCTGT, NOTGT, NUMTGT
      INTEGER, DIMENSION(nocsf0) :: KPT
      INTEGER, DIMENSION(*) :: NOB, NOBC
      INTENT (IN) GUCONT, ICITG, IEXPC, LUSME, MCONT, NCTGT, NOB, NOBC, 
     &            NOCSF0, NOTGT, NTGSYM, SYMTYP
      INTENT (OUT) KPT
      INTENT (INOUT) MOCSF, MXCSF, NCITOT, NCONT, NCONT2, NTGT, NUMMX, 
     &               NUMTGT
C
C Local variables
C
      INTEGER :: I, IADD, ICSF, ISYM, JCSF, JORB, M, N
C
C*** End of declarations rewritten by SPAG
C
c     MKPT sets up pointers for mapping when expanding energy expresions
C     The pointers are set for MOCSF new CSFs expanded from NOCSF0
C     prototype CSFs generated in CONGEN as follows:
C
C     KPT orbital sequence number for the continuum orbital in this CSF
c
c     Also determine size parameters for CI target contraction
c
c
      mxcsf=nocsf0
      icsf=0
      jcsf=0
      ntgt=0
      ncont2=0
      ncont=0
      DO n=1, ntgsym
         IF(numtgt(n).LT.1)numtgt(n)=1
         nummx=max(nummx,numtgt(n))
         ntgt=ntgt+numtgt(n)
         ncitot=ncitot+nctgt(n)*numtgt(n)
         ncont2=ncont2+notgt(n)*numtgt(n)
         IF(iexpc.EQ.0)THEN
            ncont=ncont+nctgt(n)*notgt(n)
            CYCLE
         END IF
c     extra CSFs for nth target state
         iadd=(notgt(n)-2)
         mxcsf=mxcsf+iadd*nctgt(n)
c     symmetry of nth target state
         IF(symtyp.NE.1)THEN
            isym=mcont(n)+1
         ELSE
            isym=2*mcont(n)+1
            IF(mod(mcont(n),2).EQ.0 .AND. gucont(n).LT.0)isym=isym+1
            IF(mod(mcont(n),2).NE.0 .AND. gucont(n).GT.0)isym=isym+1
         END IF
c     orbital sequence number of first continuum orbital for nth target
         jorb=0
         DO i=1, isym-1
            jorb=jorb+nob(i)
         END DO
         jorb=jorb+nobc(isym)+1
         DO m=1, nctgt(n)
            icsf=icsf+1
            jcsf=jcsf+1
c     first  continuum orbital for each target state:
            kpt(icsf)=jorb
            icsf=icsf+1
            jcsf=jcsf+1
c     second continuum orbital for each target state:
            kpt(icsf)=jorb+1
            jcsf=jcsf+iadd
         END DO
      END DO
      IF(iexpc.NE.0)THEN
c        number of prototype continuum CSF
         ncont=icsf
c        pointers for L**2 CSFs
         DO i=icsf+1, nocsf0
            jcsf=jcsf+1
            kpt(i)=0
         END DO
      END IF
      IF(icitg.EQ.0)THEN
c     reset parameters when no CI target contraction requested
         DO n=1, ntgsym
            IF(nctgt(n).NE.1)THEN
               WRITE(6,99)
 99            FORMAT(/
     &              ' Attempt to do expand prototype CSFs for CI target'
     &              ,/' Option not implemented: STOP')
               STOP
            END IF
         END DO
      END IF
      mocsf=ncont2+nocsf0-ncont
      IF(LUSME.NE.0)WRITE(lusme,*)mocsf
      RETURN
      END SUBROUTINE MKPT
!*==movep.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE MOVEP(NFT,NTH,NALM,NPFLG,NFT1)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE scatci_data, ONLY : MEIG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NALM, NFT, NFT1, NPFLG, NTH
      INTENT (IN) NFT, NFT1, NPFLG
      INTENT (OUT) NALM
      INTENT (INOUT) NTH
C
C Local variables
C
      REAL(KIND=wp) :: E0, S, SZ
      REAL(KIND=wp), DIMENSION(MEIG) :: EIG
      CHARACTER(LEN=32) :: HEADER
      INTEGER :: I, IDUM, J, K, M, MGVN, NEIG, NELT, NNUC, NOCSF, NREC, 
     &           NSTAT, NT
      CHARACTER(LEN=120) :: NAME
      LOGICAL :: OPN
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     MOVEP LOCATES THE NTH-DATASET CONTAINED ON THE ALCHEMY CI
C           DUMPFILE (for polyatomics) ASSOCIATED WITH UNIT NFT
C
C           NPFLG = 0     NO PRINTOUT DURING SEARCH
C                 = 1     PRINT HEADER LABELS ENCOUNTERED
C           NALM  = 0     NO ERRORS DETECTED
C                 = 1     ERRORS DETECTED
C
C***********************************************************************
C
C
      NALM=0
      INQUIRE(UNIT=nft,OPENED=opn)
      IF(opn)THEN
         REWIND NFT
      ELSE
         OPEN(UNIT=nft,FORM='unformatted')
      END IF
      IF(NTH.EQ.1)RETURN
      M=NTH-1
      IF(NTH.EQ.0)M=2000
      DO I=1, M
         READ(NFT,END=50)header
         READ(NFT,END=50)NT, nrec, NAME, nnuc, nocsf, nstat, mgvn, s, 
     &                   sz, nelt, e0
         DO K=1, NNUC
            READ(NFT)
         END DO
C
         NEIG=MIN(MEIG,NSTAT)
         READ(NFT)(IDUM,J=1,nocsf), (EIG(J),J=1,NEIG)
C
         IF(NPFLG.GT.0)THEN
            WRITE(NFT1,104)header
            WRITE(NFT1,100)NT, name
            WRITE(nft1,103)mgvn, s, sz, nelt, nnuc
            WRITE(NFT1,101)nocsf, nstat, e0
            WRITE(NFT1,102)(EIG(J)+E0,J=1,NEIG)
         END IF
 104     FORMAT(/' Header:',4x,a)
 100     FORMAT(/' SET',I4,4X,A)
 101     FORMAT(/' NOCSF=',I5,4X,'NSTAT=',I5,4X,'EN   =',F20.10)
 103     FORMAT(/' MGVN =',I2,4x,'s =',f6.1,4x,'sz =',f6.1,4x,'NELT =',
     &          I5,4x,'NNUC =',I3)
 102     FORMAT(/' EIGEN-ENERGIES',/(16X,5F20.10))
c
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
c        BACKSPACE NFT
      END IF
      RETURN
      END SUBROUTINE MOVEP
!*==movew.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE MOVEW(NFT,NTH,NALM,NPFLG,NFT1)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE scatci_data, ONLY : MEIG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NALM, NFT, NFT1, NPFLG, NTH
      INTENT (IN) NFT, NPFLG
      INTENT (OUT) NALM
      INTENT (INOUT) NTH
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp), DIMENSION(MEIG) :: EIG
      INTEGER :: I, IDUM, J, K, M, NEIG, NHDIM, NSET, NSTAT, NT
      CHARACTER(LEN=120) :: NAME
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(20) :: NHE
      LOGICAL :: OPN
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     MOVEW LOCATES THE NTH-DATASET CONTAINED ON THE ALCHEMY CI
C           DUMPFILE ASSOCIATED WITH UNIT NFT
C
C           NPFLG = 0     NO PRINTOUT DURING SEARCH
C                 = 1     PRINT HEADER LABELS ENCOUNTERED
C           NALM  = 0     NO ERRORS DETECTED
C                 = 1     ERRORS DETECTED
C
C***********************************************************************
C
C
      NALM=0
      INQUIRE(UNIT=nft,OPENED=opn)
      IF(opn)THEN
         REWIND NFT
      ELSE
         OPEN(UNIT=nft,FORM='unformatted')
      END IF
      IF(NTH.EQ.1)RETURN
      M=NTH-1
      IF(NTH.EQ.0)M=2000
      DO I=1, M
         READ(NFT,END=50)NT, NHD, NAME, NHE, DTNUC
         NSET=NT
         NHDIM=NHD(8)
         NSTAT=NHD(3)
         NEIG=MIN(MEIG,NSTAT)
         READ(NFT)(IDUM,J=1,NHDIM), (EIG(J),J=1,NEIG)
C
         IF(NPFLG.GT.0)CALL PRTHD(NT,NHD,NAME,NHE,DTNUC,NEIG,EIG,NFT1)
         DO K=1, NSTAT+nhd(10)
            READ(NFT)
         END DO
C
      END DO
      RETURN
C
 50   IF(NTH.NE.0)THEN
         NALM=1
      ELSE
         NTH=NSET+1
         BACKSPACE NFT
      END IF
      RETURN
      END SUBROUTINE MOVEW
!*==mvdiag.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE MVDIAG(CH2,XMH2,CH1,XMH1,NH,NHC,THRES)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NH, NHC
      REAL(KIND=wp) :: THRES
      REAL(KIND=wp), DIMENSION(*) :: CH1, CH2
      INTEGER(longint), DIMENSION(2,*) :: XMH1, XMH2
      INTENT (IN) CH1, NH, THRES, XMH1
      INTENT (OUT) CH2, XMH2
      INTENT (INOUT) NHC
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
      DO I=1, NH
         IF(ABS(CH1(I)).LT.THRES)CYCLE
         NHC=NHC+1
         XMH2(1,NHC)=XMH1(1,I)
         XMH2(2,NHC)=XMH1(2,I)
         CH2(NHC)=CH1(I)
      END DO
c
      RETURN
      END SUBROUTINE MVDIAG
!*==nin2ea.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      FUNCTION nin2ea(lwd,nri,nsm,nob,ipair,mbas,nbas,map,iodr)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IODR
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER, DIMENSION(5) :: LWD
      INTEGER, DIMENSION(*) :: MAP, MBAS, NBAS, NOB, NRI, NSM
      INTEGER :: NIN2EA
      INTENT (IN) IODR, IPAIR, LWD, MAP, MBAS, NBAS, NOB, NRI, NSM
C
C Local variables
C
      INTEGER :: I, IBAS, IP, IPQRS, M, MD, MPQ, MPR, MRS, MV
      INTEGER, DIMENSION(4) :: MPP, NBP, NPP, NPQ
C
C*** End of declarations rewritten by SPAG
C
C    NIN2eA gives a pointer to a 2-electron integral for ALCHEMY
C
C
      DO I=1, 4
         IP=map(LWD(I))
         NPP(I)=NRI(IP)
         M=NSM(IP)
         MPP(I)=M
         NBP(I)=NOB(M+1)
      END DO
C
      IF(lwd(5).NE.0)THEN
         MV=MPP(1)+MPP(2)
      ELSE
         MV=ABS(MPP(1)-MPP(2))
      END IF
      IF(MV.NE.0)THEN
         MD=(MV-1)/2
      ELSE
         MD=-1
      END IF
C
      MPQ=MPP(1)-MD
      MRS=MPP(3)-MD
      MPR=IPAIR(MPQ)+MRS+MBAS(MV+1)
C
      IBAS=NBAS(MPR)-1
      DO I=1, 4, 2
         IF(MPP(I).NE.MPP(I+1))THEN
            NPQ(I)=NBP(I)*(NPP(I+1)-1)+NPP(I)
            NPQ(I+1)=NBP(I)*NBP(I+1)
         ELSE
            IF(NPP(I).LT.NPP(I+1))THEN
               NPQ(I)=IPAIR(NPP(I+1))+NPP(I)
            ELSE
               NPQ(I)=IPAIR(NPP(I))+NPP(I+1)
            END IF
            NPQ(I+1)=IPAIR(NBP(I)+1)
         END IF
      END DO
C
      IF(MPP(1).EQ.MPP(3) .AND. MPP(2).EQ.MPP(4))THEN
         IF(NPQ(1).LT.NPQ(3))THEN
            IPQRS=IPAIR(NPQ(3))+NPQ(1)
         ELSE
            IPQRS=IPAIR(NPQ(1))+NPQ(3)
         END IF
      ELSE
         IF(IODR.EQ.0)THEN
            IPQRS=NPQ(2)*(NPQ(3)-1)+NPQ(1)
         ELSE
            IPQRS=NPQ(4)*(NPQ(1)-1)+NPQ(3)
         END IF
      END IF
C
      nin2ea=IPQRS+IBAS
      RETURN
      END FUNCTION NIN2EA
!*==nin2em.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      FUNCTION nin2em(lwd,nri,nsm,nob,ipair,nbas,iodr,map)
C
C**********************************************************************
C
C     This routine (following M. YOSHIMINE'S routine, PINDEX,
C     for ALCHEMY) generates a single index for energy expression
C     integrals from and index pair or quadrupole; index is consistent
C     with SWEDEN integral block ordering i.e. D2h
C
C     In order to change this code to work with another quantum
C     chemistry integals system it is only necessary to modify
C     the routine which constructs the table of blocks.
C
C     NIN2eM calculates the index for 2-electron integral MP,MQ,MR,MS
C                  ASSUME MP GE MQ, MR GE MS, MPQ GE MRS
C                  NOTHING FURTHER
C
C     Input data:
C
C     IODR     ORDER OF TRANSFORMED INTEGRALS
C              SEE CODE FOR FURTHER NOTES
C     NRI(NOBT)    ORBITAL NUMBER WITHIN SYMMETRY
C     NSM(NOBT)    ORBITAL M-VALUE (NSYM-1)
C     NOB(NSYM)   NUMBER OF ORBS PER SYMMETRY
C     NBAS(NBLK+1) SERIAL NUMBER OF FIRST INTEGRAL IN (LMIN0) BLOCK
C                  OF TRANSFORMED INTEGRALS COUNT STARTS FROM NUMBER
C                  OF H-INTEGRALS+1
C                  NBAS(NBLK+1) IS TOTAL NUMBER OF INTS
C     IPAIR(MAX)   IPAIR(I)=(I*(I-1))/2
C                  MAX IS DETERMINED BY THREE CONDITIONS
C              1.  MAX NUMBER OF SYM PAIRS NSYM*(NSYM+1)/2
C              2.  MAX NOB(I)+1
C              3.  MAX ROW LENGTH FOR LM/LM OR LL/LL BLOCK
C                  THAT IS MAX OF NOB(I)*NOB(J) AND
C                        NOB(I)*(NOB(I)+1)/2
C                  THIS MAY BE REDUCED TO MAX OF (CONDITION 1. PLUS 1 A
C                  COND 3.)
C
C**********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: NOBTEST
      COMMON /NOBDEBUG/ NOBTEST
C
C Dummy arguments
C
      INTEGER :: IODR
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER, DIMENSION(4) :: LWD
      INTEGER, DIMENSION(*) :: MAP, NBAS, NOB, NRI, NSM
      INTEGER :: NIN2EM
      INTENT (IN) IODR, IPAIR, LWD, MAP, NBAS, NOB, NRI, NSM
C
C Local variables
C
      REAL(KIND=wp) :: AAA
      INTEGER :: I, IBAS, IP, IPQRS, M, MPQ, MPR, MRS, NOBPQ, NOBRS
      INTEGER, DIMENSION(4) :: MPP, NBP, NPP
C
C*** End of declarations rewritten by SPAG
C
      DO I=1, 4
         IP=map(LWD(I))
         NPP(I)=NRI(IP)
         M=NSM(IP)+1
         MPP(I)=M
         NBP(I)=NOB(M)
      END DO
C
C...... It is useful to remember that
C
C       NPP is the orbital number within its symmetry
C       MPP is the orbital symmetry
C       NBP is the number of orbitals total in the symmetry
C
C       Each of these is for the quadruplet formula element
C
 
      MPQ=IPAIR(MPP(1))+MPP(2)
      MRS=IPAIR(MPP(3))+MPP(4)
      MPR=IPAIR(MPQ)+MRS
 
c
c     after mapping, the positronic orbital
c     might have a lower orbital number
c     than the electronic orbital
c     only one case possible:
c     (AB|ij)
c     here A,B = positronic orbital
c          i,j = electronic orbital
c     after mapping this is
c     (ab|ij)
c     here a=map(A) and b=map(B)
c     this might be translated into
c     (ij|ab), if a<i or if a=i and b<j
c     the order a>b and i>j are conserved
c     integrals of the type (Ai|Bj) will not
c     occur.
c
      IF(MPQ.EQ.MRS .AND. 
     &   (npp(1).LT.npp(3) .OR. (npp(1).EQ.npp(3) .AND. npp(2).LT.npp(4)
     &   ) .OR. (npp(1).LT.npp(3) .AND. npp(2).LT.npp(4))))THEN
         aaa=npp(1)
         npp(1)=npp(3)
         npp(3)=aaa
         aaa=mpp(1)
         mpp(1)=mpp(3)
         mpp(3)=aaa
         aaa=npp(2)
         npp(2)=npp(4)
         npp(4)=aaa
         aaa=mpp(2)
         mpp(2)=mpp(4)
         mpp(4)=aaa
      END IF
 
      IBAS=NBAS(MPR)-1
      IF(mpp(1).LT.mpp(2) .OR. mpp(3).LT.mpp(4) .OR. mpq.LT.mrs .OR. 
     &   lwd(1).LT.lwd(2) .OR. lwd(3).LT.lwd(4) .OR. ipair(lwd(1))
     &   +lwd(2).LT.ipair(lwd(3))+lwd(4))THEN
         WRITE(6,900)lwd, mpp, npp
      END IF
C
C...... IBAS is the serial number less one of the storage for this
C       particular block
C
C       Now, for each quadruplet, find the relative storage position
C       within its symmetry block. There are one of four possibilities:
C
C                 a) IIII - all the same symmetry
C                 b) IJIJ - block pairs identical
C                 c) IIJJ - block pairs not identical
C                 d) IJKL - all different
C
C       Each is dealt with separately here. Note that the first three
C       involve triangularity on two levels. The index here must be
C       consistent with the ordering adopted by the integrals code.
C       The calculation here mirrors that adopted in the Sweden
C       ordering code and thus is for Sweden transformed integrals.
C
C       Note that for the fourth type of block - IJKL - where all
C       symmetries are different there are two kinds of ordering
C       in common use. In Alchemy these are defined as
C
C                    PQRS and RSPQ  (IODR=1 and IODR=0)
C
C       It appears that Sweden uses the format PQRS while Alchemy
C       classed RSPQ as the 'normal' format. The variable IODR
C       controls which type of ordering is used; the default is to
C       have IODR=1 but it can be changed in the namelist input.
C
C
      IF(MPP(1).EQ.MPP(2) .AND. MPP(3).EQ.MPP(4) .AND. MPP(1).EQ.MPP(3))
     &   THEN
         IPQRS=IPAIR(NPP(1))+NPP(2)
         IF(ipqrs.GT.nobtest)WRITE(6,910)ipqrs, nobtest
         IPQRS=IPAIR(IPQRS-1)+IPQRS-1+IPAIR(NPP(3))+NPP(4)
      ELSE IF(MPP(1).EQ.MPP(3) .AND. MPP(2).EQ.MPP(4))THEN
         IPQRS=(NPP(1)-1)*NBP(2)+NPP(2)-1
         IF(ipqrs.GT.nobtest)WRITE(6,910)ipqrs, nobtest
         IPQRS=IPAIR(IPQRS)+IPQRS+(NPP(3)-1)*NBP(4)+NPP(4)
      ELSE IF(MPP(1).EQ.MPP(2) .AND. MPP(3).EQ.MPP(4))THEN
         NOBRS=IPAIR(NBP(3))+NBP(3)
         IPQRS=(IPAIR(NPP(1))+NPP(2)-1)*NOBRS+IPAIR(NPP(3))+NPP(4)
      ELSE IF(IODR.EQ.0)THEN
         NOBPQ=NBP(1)*NBP(2)
         IPQRS=((NPP(3)-1)*NBP(4)+NPP(4)-1)*NOBPQ+(NPP(1)-1)*NBP(2)
     &         +NPP(2)
      ELSE IF(IODR.EQ.1)THEN
         NOBRS=NBP(3)*NBP(4)
         IPQRS=((NPP(1)-1)*NBP(2)+NPP(2)-1)*NOBRS+(NPP(3)-1)*NBP(4)
     &         +NPP(4)
      END IF
C
C...... Return the absolute storage pointer. Note we must
C       bias IPQRS by one to account for the fact that Sweden stores
C       the block descriptor as the first element of the block.
C
      nin2em=IPQRS+IBAS+1
C
C---- Subroutine return point
C
      RETURN
 900  FORMAT(' Possible integral indexing error in NIN2EM:',/' LWD =',
     &       4I4,/' MPP =',4I4,/' NPP =',4I4,
     &       /' Please contact Jonathan Tennyson')
 910  FORMAT(' Probable integral indexing error in NIN2EM:',/' IPQRS =',
     &       i8,' maximum allowable =',i8,
     &       /' Please contact Jonathan Tennyson')
C
      END FUNCTION NIN2EM
!*==nind1e.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      FUNCTION nind1e(m,np,nq,ipair,nbash)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: M, NP, NQ
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER, DIMENSION(*) :: NBASH
      INTEGER :: NIND1E
      INTENT (IN) IPAIR, M, NBASH, NP, NQ
C
C*** End of declarations rewritten by SPAG
C
C    NIND1e gives a pointer to a 1-electron integral
C    (works for both ALCHEMY & SWEDEN)
C
      IF(NP.LT.NQ)THEN
         nind1e=IPAIR(NQ)+NP+NBASH(M+1)-1
      ELSE
         nind1e=IPAIR(NP)+NQ+NBASH(M+1)-1
      END IF
      RETURN
      END FUNCTION NIND1E
!*==nobid.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE NOBID(NOBT,NSYM,NOB,NRI,NSM,nobmx,LUSME)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: NOBTEST
      COMMON /NOBDEBUG/ NOBTEST
C
C Dummy arguments
C
      INTEGER :: LUSME, NOBMX, NOBT, NSYM
      INTEGER, DIMENSION(*) :: NOB, NRI, NSM
      INTENT (IN) LUSME, NOB, NSYM
      INTENT (INOUT) NOBMX, NOBT, NRI, NSM
C
C Local variables
C
      INTEGER :: I, L
C
C*** End of declarations rewritten by SPAG
C
C     Generates an orbital table from the definition of NSYM and NOB.
C
      NOBT=0
      nobmx=0
      DO L=1, NSYM
         nobmx=max(nobmx,nob(l)**2)
         DO I=1, NOB(L)
            NOBT=NOBT+1
            NRI(NOBT)=I
            NSM(NOBT)=L-1
         END DO
      END DO
      nobtest=nobmx
c
      IF(LUSME.NE.0)THEN
         WRITE(lusme,*)NOBT, NOBMX
         WRITE(lusme,*)(NRI(I),I=1,NOBT), (NSM(I),I=1,NOBT)
      END IF
      RETURN
      END SUBROUTINE NOBID
!*==pindex.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
 
      SUBROUTINE PINDEX(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,KT,XNIP,
     &                  CFP,ndim,xint1e,nint1e,xint2e,elem,isymtp,map,
     &                  NCORB,mcorb,nci,LUSME,enh_factor,
     &                  enh_factor_type,ukrmolp_ints,NOEX)
C
C***********************************************************************
C
C     Computes a sequential index for each formula element which is used
C     to access the integral arrays XINT1e and XINT2e.
C     Loops over ELEM allow for case of CI target contraction with
C     more than one target state for a given symmetry.
C     This routine is for C-infinity-v symmetry (ALCHEMY) and SWEDEN.
C
C     ZM: modified to work with UKRMol+ integrals. This is controlled by 
C         the input flag ukrmolp_ints. Note that the routine can be rewritten
C         to take advantage of UKRMol+ capability to compute multiple inidices
C         at once. I.e. the orbital indices can be unpacked first and then passed
C         at once to the indexing routines which should result in some perf. gain
C         compared to the one-by-one calls.
C
C     KT      NUMBER OF INDICES PASSED
C     NIP(KT)  PACKED INDICES / ONE PER BYTE
C              H-INTEGRALS (0,NP,0,NQ)
C                  NO ASSUMPTION MADE ON ORDER OF NP AND NQ
C              2-INTEGRALS (NP,NQ,NR,NS)
C                  ASSUME MP GE MQ, MR GE MS, MPQ GE MRS
C                  NOTHING FURTHER
C
C     IODR = 0, PQRS AND =I, RSPQ (for ALCHEMY) (usual IODR=1)
C     IODR = 0 or 1 for SWEDEN
C
C     IPAIR(I)=I*(I-1)/2 ADDED TO M. YOSHIMINE'S ROUTINE
C
C**********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE consts, ONLY : ZERO=>XZERO
      USE integer_packing, ONLY : unpack8ints
      USE ukrmol_interface_gbl, ONLY : GET_INTEGRAL, GET_ORB_ENERGY
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IODR, ISYMTP, KT, LUSME, NCI, NCORB, NDIM, NINT1E, NOEX
      REAL(KIND=wp), DIMENSION(ndim,nci) :: CFP
      REAL(KIND=wp), DIMENSION(nci) :: ELEM
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER, DIMENSION(*) :: MAP, MBAS, MCORB, NBAS, NBASH, NOB, NRI, 
     &                         NSM
      REAL(KIND=wp), ALLOCATABLE :: XINT1E(:), XINT2E(:)
      INTEGER(longint), DIMENSION(2,kt) :: XNIP
      LOGICAL :: ukrmolp_ints
      INTENT (IN) CFP, ISYMTP, KT, LUSME, MCORB, NCI, NCORB, NDIM, 
     &            NINT1E, XINT1E, XINT2E, enh_factor, ukrmolp_ints,
     &            enh_factor_type
      INTENT (INOUT) ELEM
C
C Local variables
C
      INTEGER :: I, IA, IB, IDX, IQ, NK, orb_n, orb_num,
     &           belec_orb, velec_orb, posi_orb, enh_factor_type
      INTEGER, DIMENSION(8) :: LWD
      INTEGER, DIMENSION(4) :: mapped
      REAL(KIND=wp) :: X, enh_factor, orb_diff, oe,
     &                 mp2_corr, x1
      real(kind=wp), dimension(4) :: occu_orbe, virt_orbe
      logical :: with_pos
C
C*** End of declarations rewritten by SPAG
C
C---- Loop to 200 is over formula elements. These are assigned a
C     sequential index which is calculated from the orbital labels.
C
C      Two possibilities exist here:
C
C   1) Two electron integrals (treated differently by Alchemy & Sweden & UKRmol+)
C   2) One electron integrals (treated the same in Alchemy & Sweden but differently in UKRmol+)
C
c
      elem=zero
c
      DO NK=1, KT
C
C..... Unpack the orbital indices into array LWD; then scan the
C      quartet for the occurence of the dummy continuum orbital.
C      Should this be found, set the coefficient to zero and then
C      proceed to sorting as normal. In the case that no dummy
C      continuum is being used (NCORB=-1), skip the search and
C      start ordering immediately.
C
         CALL unpack8ints(XNIP(1,NK),LWD)
c
         IF(NCORB.GT.0)THEN
            DO IQ=1, 4
               IF(LWD(IQ).EQ.NCORB)GO TO 200
            END DO
         ELSE IF(NCORB.EQ.0)THEN
            DO IQ=1, 4
               IF(lwd(iq).LE.0)CYCLE
               IF(mcorb(LWD(IQ)).EQ.0)GO TO 200
            END DO
         END IF
C
         IF(LWD(1).EQ.0)THEN
c        1-electron case
            IA=map(lwd(2))
            IB=map(lwd(4))
            if (NOEX .eq. 1) lwd(5) = 0
            if (isymtp.ge.2 .and. ukrmolp_ints) then !branch to the UKRMol+ integral
               x=get_integral(ia,ib,0,0,lwd(5))
            else !SWEDEN/ALCHEMY
               idx=nind1e(nsm(ia),nri(ia),nri(ib),ipair,nbash)
c              positron case
               IF(lwd(5).NE.0)idx=idx+nint1e
               x=xint1e(idx)
            endif
         ELSE IF(isymtp.LE.1)THEN
c        2-electron  ALCHEMY case
            idx=nin2ea(lwd,nri,nsm,nob,ipair,mbas,nbas,map,iodr)
            x=xint2e(idx)
         ELSE IF (isymtp.ge.2 .and. (.not. ukrmolp_ints))THEN
c        2-electron SWEDEN case
            idx=nin2em(lwd,nri,nsm,nob,ipair,nbas,iodr,map)
            x=xint2e(idx)
         ELSE IF (isymtp.ge.2 .and. ukrmolp_ints)THEN
c        2-electron UKRMol+ case
            mapped(1) = map(LWD(1))
            mapped(2) = map(LWD(2))
            mapped(3) = map(LWD(3))
            mapped(4) = map(LWD(4))
            x=get_integral(mapped(1),mapped(2),mapped(3),mapped(4),
     &                     lwd(5))

            !~~~~~~~~~~~ENHANCEMENT FACTOR INCLUSION ~~~~~~~~~~~~~~~~
            if (enh_factor_type .ne. 0) then
             with_pos=.False.
             do i=1, 4
              if (mapped(i) .ne. lwd(i)) then
                with_pos=.True.
              endif
             enddo
! Now multiply elements depending on which enhancement factor has been
! chosen: MP2 produces an integral dependent factor.
            if (enh_factor_type .gt. 1) then
              if (with_pos) then
                belec_orb = 0
                velec_orb = 0
                posi_orb = 0

                do orb_n=1, 4
                  virt_orbe(orb_n) = 0.0
                  occu_orbe(orb_n) = 0.0

                  orb_num = mapped(orb_n) !orb index

                  oe = GET_ORB_ENERGY(orb_num)
                  if (oe .gt. 0.0) then !virtual or continuum orbital
                    virt_orbe(orb_n) = oe
                    if (mapped(orb_n) .ne. lwd(orb_n)) then
                      velec_orb = velec_orb + 1
                    else
                      posi_orb = posi_orb + 1
                    endif
                  else !bound orbital
                    occu_orbe(orb_n) = oe
                    if (mapped(orb_n) .ne. lwd(orb_n)) then
                      belec_orb = belec_orb + 1
                    else
                      posi_orb = posi_orb + 1
                    endif
                  endif
                enddo

                orb_diff = 0.0
                if (belec_orb .eq. 1 .and. velec_orb .eq. 1) then
                  !Energy diff between bound and virtual/continuum
                  !orbitals
                  orb_diff = sum(occu_orbe(1:4)) - sum(virt_orbe(1:4))!Will always be negative
                  x1 = x
                elseif (mapped(1) .eq. mapped(2) .and.
     &                  mapped(3) .eq. mapped(4) .and.
     &                  mapped(1) .eq. mapped(3) .and.
     &                  velec_orb .ne. 0) then
                  !Electron and positron occupying the same orbital -
                  !virtual Ps
                  orb_diff= 1.0 !Orbital energy that all particles are in.
                  x1 = 1
                endif

                if (orb_diff .ne. 0.0) then !include MP2 correction
                  if (x .ge. 0) then
                    mp2_corr = ((x*x1) / orb_diff) * enh_factor
                  else
                    mp2_corr = (-(x*x1) / orb_diff) * enh_factor
                  endif
                  x = x + mp2_corr
                endif
              endif
            else
              if (with_pos) then
                x = x * enh_factor !use provided average enhancement factor. Apply to all (positron) 2-particle integrals
              endif
            endif
            endif
            !~~~~~~~~~~~~~END ENHANCEMENT FACTOR INCLUSION~~~~~~~~~~~~~~~~~~~~

         END IF
         DO i=1, nci
            elem(i)=elem(i)+cfp(nk,i)*x
         END DO
 
         IF(LUSME.NE.0)THEN
            WRITE(lusme,*)lwd(1), lwd(2), lwd(3), lwd(4), lwd(5)
            WRITE(lusme,*)(cfp(nk,i),i=1,nci)
         END IF
 
 200  END DO
      RETURN
      END SUBROUTINE PINDEX
!*==pmkorbs.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE PMKORBS(nob,nobe,nob0,nsym,mn,mg,mm,ms,mcon,mcorb,NORB,
     &                   NSRBD,map,mpos,iposit,symtyp,LUSME,NFTW)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IPOSIT, LUSME, NORB, NSRBD, NSYM, SYMTYP, NFTW
      INTEGER, DIMENSION(norb) :: MAP, MCORB
      INTEGER, DIMENSION(nsrbd) :: MCON, MG, MM, MN, MS
      INTEGER, DIMENSION(2*nsrbd) :: MPOS
      INTEGER, DIMENSION(nsym) :: NOB, NOB0, NOBE
      INTENT (IN) LUSME, NOB, NOB0, NOBE, NORB, NSRBD, NSYM, SYMTYP
      INTENT (IN) NFTW
      INTENT (OUT) MCON, MCORB, MG, MM, MN, MPOS, MS
      INTENT (INOUT) MAP
C
C Local variables
C
      INTEGER :: AMO, EMO, ICON, IMO, IPOS, ISO, ISPIN, ISYM, JMO, 
     &           MAXMO, MAXSPIN, MINMO, NSRB
C
C*** End of declarations rewritten by SPAG
C
c     setting up the following arrays:
c       mn() = orbital number
c       mg() = ??? destinguish gerade and ungerade
c       mm() = ??? m-quantum number for degenerate MOs
c       ms() = spin (for alpha=0, for beta=1)
c       mpos() = flag for positron (for e-=0, for p+=1)
c     NOTE: only implemented for poly-atomic code (symtyp=2)
 
c.........................................
      IF(symtyp.NE.2)THEN
         WRITE(NFTW,*)' ERROR in PMKORBS: calculation with positrons'
         WRITE(NFTW,*)'                   only possible for SYMTYP=2'
         WRITE(NFTW,*)'                   (abelian groups).'
         WRITE(NFTW,*)' here: SYMTYP=', SYMTYP
         STOP
 
      END IF
 
      maxspin=2
 
c     imo = mo-number
c     iso = so-number
c     ipos = positron-flag
c          = 0 for 1..nobe(isym)
c          = 1 for nobe(isym)..nob(isym)
 
      imo=0
      iso=0
      emo=0
      DO isym=1, nsym
 
c     electronic MOs
 
         ipos=0
         icon=1
         maxmo=nobe(isym)
         amo=emo
         DO jmo=1, maxmo
            imo=imo+1
            emo=emo+1
            IF(jmo.GT.nob0(isym))THEN
               icon=0
            END IF
            mcorb(imo)=icon
            map(imo)=emo
            DO ispin=1, maxspin
               iso=iso+1
               mn(iso)=imo
               mg(iso)=0
               mm(iso)=isym-1
               ms(iso)=ispin-1
               mpos(iso)=ipos
               mcon(iso)=icon
ccc            write(6,*),iso,imo,mn(iso)
            END DO
         END DO
 
c     positronic MOs
 
         ipos=1
         minmo=nobe(isym)+1
         maxmo=nob(isym)
         icon=0
c        shift=ipos*nobe(isym)
         emo=amo
         DO jmo=minmo, maxmo
            imo=imo+1
            emo=emo+1
            mcorb(imo)=icon
c          map(imo)=imo-shift
            map(imo)=emo
            DO ispin=1, maxspin
               iso=iso+1
               mn(iso)=imo
               mg(iso)=0
               mm(iso)=isym-1
               ms(iso)=ispin-1
               mpos(iso)=ipos
               mcon(iso)=icon
            END DO
         END DO
 
      END DO
 
      nsrb=iso
 
      IF(LUSME.NE.0)THEN
         WRITE(LUSME,*)NORB
         WRITE(LUSME,*)MAP
      END IF
c     control number of spin orbitals
 
      WRITE(NFTW,*)'GIVEN NSRBD=', NSRBD
      WRITE(NFTW,*)'CALCULATED NSRB=', NSRB
 
      IF(NSRB.NE.NSRBD)THEN
         WRITE(NFTW,1010)NSRB, NSRBD
 1010    FORMAT(' HELP!!! MKORBS: NSRB, NSRBD = ',2I6)
         STOP
      END IF
      RETURN
      END SUBROUTINE PMKORBS
!*==posadd.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE POSADD(LPOSIT,NDTC,NORB,I1,I2,MPOS)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: I1, I2, LPOSIT, NORB
      INTEGER, DIMENSION(*) :: MPOS
      INTEGER, DIMENSION(4) :: NDTC
      INTENT (IN) I1, I2, MPOS, NDTC
      INTENT (OUT) LPOSIT
C
C Local variables
C
      INTEGER :: I, ICOUNT, IP
C
C*** End of declarations rewritten by SPAG
C
      ICOUNT=0
      lposit=0
      DO I=I1, I2
         IP=NDTC(I)
         IF(MPOS(IP).NE.0)ICOUNT=ICOUNT+1
      END DO
      IF(ICOUNT.EQ.0)RETURN
      IF(ICOUNT.EQ.2)THEN
         LPOSIT=1
         RETURN
      END IF
      WRITE(6,*)' ERROR IN POSADD - INVALID ICOUNT'
      WRITE(6,*)' ICOUNT, NDTC = ', ICOUNT, (NDTC(I),I=1,4)
      STOP
      END SUBROUTINE POSADD
!*==pqrs.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE PQRS(NP,NQ,NR,NS,MLA,MS,CFD,XMJK,CJK)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C COMMON variables
C
      INTEGER :: NALM, NJK, NJKMX
      COMMON /STPQ  / NJK, NALM, NJKMX
C
C Dummy arguments
C
      REAL(KIND=wp) :: CFD
      INTEGER :: MLA, MS, NP, NQ, NR, NS
      REAL(KIND=wp), DIMENSION(*) :: CJK
      INTEGER(longint), DIMENSION(2,*) :: XMJK
      INTENT (IN) CFD, MS, NP, NQ, NR, NS
      INTENT (INOUT) CJK
C
C Local variables
C
      INTEGER, SAVE :: IA, IB, IC, ID, IJK
      REAL(KIND=wp), SAVE :: X
      INTEGER(longint), DIMENSION(2), SAVE :: XNA
C
C*** End of declarations rewritten by SPAG
C
      IF(NP.LT.NQ)THEN
         IA=NQ
         IB=NP
      ELSE
         IA=NP
         IB=NQ
      END IF
      IF(NR.LT.NS)THEN
         IC=NS
         ID=NR
      ELSE
         IC=NR
         ID=NS
      END IF
C
      IF(IA.GT.IC .OR. (IA.EQ.IC .AND. IB.GT.ID))THEN
         CALL pack8ints(IA,IB,IC,ID,mla,0,0,0,XNA)
      ELSE
         CALL pack8ints(IC,ID,IA,IB,mla,0,0,0,XNA)
      END IF
C
c     IF (MLA .NE. 0) XNA=-XNA
      X=CFD
      IF(MS.NE.0)X=-X
      ijk=ijkpqrs(xna,xmjk,njk)
C
      IF(IJK.GT.0)THEN
         CJK(IJK)=CJK(IJK)+X
      ELSE
         NJK=NJK+1
         IF(NJK.LE.NJKMX)THEN
            XMJK(1,NJK)=XNA(1)
            XMJK(2,NJK)=XNA(2)
            CJK(NJK)=X
         ELSE
            NALM=1
            RETURN
         END IF
      END IF
      RETURN
      END SUBROUTINE PQRS
!*==prblks.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE PRBLKS(MBAS,nser1e,NSER2e,NPBXDX,nsym,NFTW)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NFTW, NPBXDX, NSER1E, NSER2E, NSYM
      INTEGER, DIMENSION(*) :: MBAS
      INTENT (IN) MBAS, NFTW, NPBXDX, NSER1E, NSER2E, NSYM
C
C Local variables
C
      INTEGER :: I, K, KK, KP, L, M1, M2, M3, M4, M4L
      INTEGER, DIMENSION(6) :: ISER
      INTEGER, DIMENSION(4,6) :: M
C
C*** End of declarations rewritten by SPAG
C
C**********************************************************************
C
C     Subroutine to print out the block labels and serial
C     starting values for ordering ALMLOFs integrals.
C
C     Input Data:
C           MBAS Array Containing serial number if block exists
C                                 zero if it does not
C         NSER?e Total number of integrals 1/2 e integrals
C         NPBXDX Print flag
C           NSYM Total number of symmetries less one
C
C**********************************************************************
C
      WRITE(NFTW,5)NSER1e, NSER2e
 5    FORMAT(/,'TOTAL NUMBER OF ORDERED 1e INTEGRALS = ',I10,/,
     &       'TOTAL NUMBER OF ORDERED 2e INTEGRALS = ',I10)
      IF(NPBXDX.LT.5)RETURN
C
      WRITE(NFTW,10)
 10   FORMAT(/,54X,'BLOCK DESCRIPTIONS',/,6('  M1 M2 M3 M4    NSER'))
C
      I=0
      K=1
      KP=1
C
C---- Descend into loop over integral blocks
C
      DO M1=1, nsym
         DO M2=1, M1
            DO M3=1, M1
               M4L=M3
               IF(M3.EQ.M1)M4L=M2
               DO M4=1, M4L
                  IF(MBAS(KP).NE.0)THEN
                     M(1,K)=M1
                     M(2,K)=M2
                     M(3,K)=M3
                     M(4,K)=M4
                     ISER(K)=MBAS(KP)
                     I=I+1
                     IF(K.EQ.6)THEN
                        WRITE(NFTW,105)((M(L,KK),L=1,4),ISER(KK),KK=1,6)
                        K=0
                        M(1,1)=0
                     END IF
                     K=K+1
                  END IF
                  KP=KP+1
               END DO
            END DO
         END DO
      END DO
C
C---- Empty any blocks not yet written out
C
      IF(K.EQ.1 .AND. M(1,1).EQ.0)THEN
         WRITE(NFTW,125)I
      ELSE
         WRITE(NFTW,105)((M(L,KK),L=1,4),ISER(KK),KK=1,K)
 105     FORMAT(6(I4,3I3,I8))
         WRITE(NFTW,125)I
 125     FORMAT(/10X,'TOTAL NUMBER OF BLOCKS = ',I10)
      END IF
C
      RETURN
      END SUBROUTINE PRBLKS
!*==prtem.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE PRTEM(NFTE,NPFLG,LEMBF,NFT)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LEMBF, NFT, NFTE, NPFLG
      INTENT (IN) NFT, NPFLG
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(LEMBF) :: EM
      INTEGER :: I, J, NT
      INTEGER, DIMENSION(2,LEMBF) :: IJ
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     PRTEM PRINTS OUT THE NONZERO ELEMENTS OF THE ENERGY MATRIX
C           STORED ON THE DISK FILE ATTACHED TO UNIT NFTE
C
C           NPFLG  LE  0  DIAGONAL ELEMENTS ONLY PRINTED
C                  GT  0  ALL NON-ZERO ELEMENTS PRINTED
C
C***********************************************************************
C
C
      REWIND NFTE
      READ(NFTE)
      CALL REMX(NFTE,NT,IJ,EM,LEMBF)
C
      IF(NPFLG.GT.0)THEN
C
         WRITE(NFT,210)
 100     CONTINUE
         WRITE(NFT,220)(IJ(1,I),IJ(2,I),EM(I),I=1,NT)
         CALL REMX(NFTE,NT,IJ,EM,LEMBF)
         IF(NT.NE.0)GO TO 100
      ELSE
C
         WRITE(NFT,230)
 130     DO I=1, NT
            J=I-1
            IF(IJ(1,I).NE.IJ(2,I))GO TO 250
         END DO
         J=0
         WRITE(NFT,220)(IJ(1,I),IJ(2,I),EM(I),I=1,NT)
         CALL REMX(NFTE,NT,IJ,EM,LEMBF)
         IF(NT.NE.0)GO TO 130
 250     IF(J.NE.0)WRITE(NFT,220)(IJ(1,I),IJ(2,I),EM(I),I=1,J)
      END IF
      REWIND NFTE
      RETURN
C
 210  FORMAT('  NON-ZERO ELEMENTS OF ENERGY MATRIX'//)
 220  FORMAT(1X,5(2I5,D16.8))
 230  FORMAT('  DIAGONAL ELEMENTS OF ENERGY MATRIX'//)
      END SUBROUTINE PRTEM
!*==prthd.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE PRTHD(NSET,NHD,HEAD,NHE,DTNUC,NEIG,EIG,NFT)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER(LEN=*) :: HEAD
      INTEGER :: NEIG, NFT, NSET
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp), DIMENSION(*) :: EIG
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
C***********************************************************************
C
C
      WRITE(NFT,100)NSET, HEAD
      E0=DTNUC(1)
      WRITE(NFT,101)NHD, E0
      WRITE(NFT,103)(NHE(I),I=1,NHD(4))
      ND=21+NHD(9)
      IF(nd.GT.200 .OR. nd.LE.0)nd=24
      WRITE(NFT,105)(DTNUC(I),I=22,ND)
      ND=1+NHD(9)
      IF(nd.GT.200 .OR. nd.LE.0)nd=4
      WRITE(NFT,104)(DTNUC(I),I=2,ND)
      WRITE(NFT,102)(EIG(I)+E0,I=1,NEIG)
C
      RETURN
C
 100  FORMAT(' SET',I4,4X,A)
 101  FORMAT(10X,'ICVC =',I5,4X,'NOCSF=',I5,4X,'NSTAT=',I5,4X,'NSYM =',
     &       I5,4X,'IDFLG=',2I5/10X,'IDEN =',I5,4X,'NHDIM=',I5,4X,
     &       'NNUC =',I5,4X,'LDEN =',I5,4X,'EN   =',F20.10)
 102  FORMAT(10X,'EIGEN-ENERGIES',/(16X,5F20.10))
 103  FORMAT(10X,'NOB  =',20I5)
 104  FORMAT(10X,'CHARG=',10F10.4/(16X,10F10.4))
 105  FORMAT(10X,'XNUC =',10F10.4/(16X,10F10.4))
C
      END SUBROUTINE PRTHD
!*==ptpwf.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
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
      INTEGER, DIMENSION(*) :: ICDI, INDI, NDI, NDTRF, NODI
      INTENT (IN) CDI, ICDI, INDI, NDI, NDTRF, NELT, NFTW, NOCSF, NODI
C
C Local variables
C
      INTEGER, SAVE :: I, K, MA, MB, MC, MD, N
C
C*** End of declarations rewritten by SPAG
C
      WRITE(NFTW,139)(NDTRF(I),I=1,NELT)
 139  FORMAT(' REFERENCE DETERMINANT'//(1X,20I5))
      WRITE(NFTW,137)
 137  FORMAT('  CSF',9X,'COEFFICIENT',2X,'NSO'/)
      DO N=1, NOCSF
         MA=NODI(N)
         MB=INDI(N)
         MC=ICDI(N)-1
         MD=NDI(MB)
         WRITE(NFTW,138)N, CDI(MC+1), MD, (NDI(MB+I),I=1,2*MD)
 138     FORMAT(1X,I4,D20.10,I5,2X,20I5/(32X,20I5))
         MB=MB+MD+MD+1
         DO K=2, MA
            MD=NDI(MB)
            WRITE(NFTW,140)CDI(MC+K), MD, (NDI(MB+I),I=1,2*MD)
 140        FORMAT(5X,D20.10,I5,2X,20I5/(32X,20I5))
            MB=MB+MD+MD+1
         END DO
      END DO
      RETURN
      END SUBROUTINE PTPWF
!*==pzero.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE PZERO(NBL,NA,NB,NZ,NBB)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NZ
      INTEGER, DIMENSION(*) :: NA, NB, NBB, NBL
      INTENT (IN) NA, NB, NBL
      INTENT (INOUT) NBB, NZ
C
C Local variables
C
      INTEGER, SAVE :: I, M, MA, MB
C
C*** End of declarations rewritten by SPAG
C
      MA=NA(1)
      MB=NB(1)
      DO I=2, MA+1
         M=NBL(NA(I))
         NBB(M)=NBB(M)-1
      END DO
      DO I=MA+2, MA+MA+1
         M=NBL(NA(I))
         NBB(M)=NBB(M)+1
      END DO
C
      DO I=2, MB+1
         M=NBL(NB(I))
         NBB(M)=NBB(M)+1
      END DO
      DO I=MB+2, MB+MB+1
         M=NBL(NB(I))
         NBB(M)=NBB(M)-1
      END DO
C
      NZ=0
      DO I=2, MA+MA+1
         M=NBL(NA(I))
         NZ=NZ+ABS(NBB(M))
         NBB(M)=0
      END DO
C
      DO I=2, MB+MB+1
         M=NBL(NB(I))
         NZ=NZ+ABS(NBB(M))
         NBB(M)=0
      END DO
      RETURN
C
      END SUBROUTINE PZERO
!*==rdins.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE RDINS(NFTA,XINT,NT,NALM)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NALM, NFTA, NT
      REAL(KIND=wp) :: XINT(NT)
      INTENT (IN) NFTA, NT
      INTENT (OUT) NALM, XINT
C
C Local variables
C
      INTEGER :: NA, NB, NC
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     Reads vector XINT of length NT from file NFTA along with
C     three other constants which are not returned. Used to read
C     transformed integrals for Hamiltonian construction
C
C***********************************************************************
C
      NALM=0
      READ(NFTA,END=200,ERR=200)NA, NB, NC, XINT
C
C---- Standard return
C
      RETURN
C
C---- Error condition
C
 200  NALM=2
      RETURN
C
      END SUBROUTINE RDINS
!*==rdint.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE RDINT(IMIN,IMAX,NEND,IMAXP,NFTA,IA,ILST,XINT)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IA, ILST, IMAX, IMAXP, IMIN, NEND, NFTA
      REAL(KIND=wp) :: XINT(:)
      INTENT (IN) ILST, IMAXP, NFTA
      INTENT (INOUT) IA, IMAX, IMIN, NEND, XINT
C
C Local variables
C
      INTEGER :: I, IB, IC, K, KC
      INTEGER, SAVE :: KA, KB
      INTEGER, DIMENSION(14), SAVE :: NPAR
C
C*** End of declarations rewritten by SPAG
C
c     Read Alchemy 2e integrals
C
C     NFTA          FILE NO. OF INTEGRAL TAPE
C     IA            INITIAL INDEX FOR XINT
C     ILST          DIMENSION OF XINT
C     IMIN          INITIAL INDEX OF THE CURRENT BLOCK
C     IMAX          LAST    INDEX OF THE CURRENT BLOCK
C     XINT(ILST)
C     NEND          =0,  MORE INTEGRALS
C                   =1,  END OF FILE
C                   =2,  ERROR
C
      IF(nend.EQ.-1)THEN
         nend=0
         GO TO 100
      END IF
C
      IB=IMAX-IMAXP
      IF(IB.GT.0)THEN
         IC=IMAXP-IMIN+1
         DO I=1, IB
            XINT(I)=XINT(IC+I)
         END DO
         IA=IB+1
      END IF
      IMAX=IMAXP
      IF(NEND.EQ.1)GO TO 300
C
 110  KC=KA
      IF(NPAR(2).NE.NPAR(4) .OR. NPAR(3).NE.NPAR(5))THEN
         DO K=KC, KB
            IB=IA+NPAR(10)-1
            IF(IB.GT.ILST)THEN
               KA=K
               GO TO 200
            ELSE
               READ(NFTA)(XINT(I),I=IA,IB)
            END IF
            IA=IB+1
         END DO
      ELSE
         DO K=KC, KB
            IB=IA+K-1
            IF(IB.GT.ILST)THEN
               KA=K
               GO TO 200
            ELSE
               READ(NFTA)(XINT(I),I=IA,IB)
            END IF
            IA=IB+1
         END DO
      END IF
C
 100  READ(NFTA,END=300,ERR=400)NPAR
      KA=1
      KB=NPAR(11)
      GO TO 110
C
 300  NEND=1
      REWIND NFTA
C
 200  IMIN=IMAX+1
      IMAX=IMIN+IA-2
      IA=1
      RETURN
C
 400  NEND=2
      RETURN
      END SUBROUTINE RDINT
!*==rdints.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE RDINTS(NFTA,NSYM,NOB,LTRIB,ILST,IA,XINT,NALM)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IA, ILST, LTRIB, NALM, NFTA, NSYM
      INTEGER, DIMENSION(*) :: NOB
      REAL(KIND=wp) :: XINT(:)
      INTENT (IN) ILST, LTRIB, NOB, NSYM
      INTENT (INOUT) IA, NALM
C
C Local variables
C
      INTEGER :: IB, L, LTRI, N, NINT, NT
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     Reads the one electron integrals for each symmetry.
C
C     Input:
C
C      NFTA  Logical unit holding the integrals
C      NSYM  Number of symmetries in the wavefunction
C       NOB  Number of orbitals per symmetry
C     LTRIB  The number of words per record on file NFTA
C      ILST  Highest storage positin to be used in array XINT
C        IA  Initial storage position to be used in array XINT
C
C     Output:
C
C       XINT  Vector holding the integrals themselves
C       NALM  Return code (0= Success ; 1= failure)
C
C***********************************************************************
C
C
C---- Descend into loop over symemtries. Note that NINT is the number
C     of one electron integrals per symmetry. Integrals are blocked
C     in units of LTRI so that more than one record is possible. See
C     routine RDINS to understand why LTRI=LTRIB-3.
C
      LTRI=LTRIB-3
      DO L=1, NSYM
         NINT=(NOB(L)*(NOB(L)+1))/2
         IF(NINT.EQ.0)CYCLE
C
C....... Loop over records of integrals and read each one. Uses
C        auxiliary routine for this.
C
         DO N=1, NINT, LTRI
            IB=N+LTRI-1
            IF(IB.GT.NINT)IB=NINT
            NT=IB-N+1
            IF(IA+NT.GT.ILST)THEN
C Error condition
               NALM=1
               RETURN
            END IF
            CALL RDINS(NFTA,XINT(IA:IA+NT),NT,NALM)
            IF(NALM.NE.0)RETURN
            IA=IA+NT
         END DO
      END DO
      RETURN
      END SUBROUTINE RDINTS
C
      SUBROUTINE SETPOSWF (IPOSIT, NCTARG, NTGCON, MEGUL)
C
C     Rewind and reposition CONGEN unit MEGUL to the beginning
C     of the wavefunction data. This needs to be compatible with
C     the subroutine RDNFTO below!
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: IPOSIT, NCTARG, NTGCON, MEGUL
C
      REWIND (MEGUL)
      READ (MEGUL)
      IF (IPOSIT == 0) THEN
         IF (NCTARG > 0) THEN
            READ (MEGUL)
            READ (MEGUL)
         END IF
      ELSE
         IF (NTGCON > 0) THEN
            READ (MEGUL)
            READ (MEGUL)
            READ (MEGUL)
            READ (MEGUL)
         END IF
         READ (MEGUL)
         READ (MEGUL)
         READ (MEGUL)
         READ (MEGUL)
         READ (MEGUL)
      END IF
C
      END SUBROUTINE SETPOSWF
!*==rdnfto.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE RDNFTO(megul,nob,nob0,nobl,nob0l,nsym,ndtrf,nelt,nodo,
     &                  nocsf,iphz,nctarg,nctgt,notgt,mcont,gucont,
     &                  ntgsym,ntgcon,nobe,nobp,nobv,iposit,NOEX)
c
c     read CSFs as generated and projected by CONGEN (part 1)
C            NODO(NOCSF) NO OF DTRS IN EACH INPUT WF
C            NDTRF(NELT) REFERENCE DTR IN EXPANDED FORM
c
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IPOSIT, MEGUL, NCTARG, NELT, NOCSF, NSYM, NTGCON, 
     &           NTGSYM, NOEX
      INTEGER, DIMENSION(ntgsym) :: GUCONT, MCONT, NCTGT, NOTGT
      INTEGER, DIMENSION(nctarg) :: IPHZ
      INTEGER, DIMENSION(NELT) :: NDTRF
      INTEGER, DIMENSION(NSYM) :: NOB
      INTEGER, DIMENSION(nsym) :: NOB0, NOBE, NOBP, NOBV
      INTEGER, DIMENSION(2*nsym) :: NOB0L, NOBL
      INTEGER, DIMENSION(NOCSF) :: NODO
      INTENT (IN) MEGUL, NCTARG, NELT, NOCSF, NSYM, NTGCON, NTGSYM
      INTENT (OUT) GUCONT, IPHZ, MCONT, NDTRF, NOB, NOB0, NOB0L, NOBE, 
     &             NOBL, NOBP, NOBV, NODO, NOTGT
      INTENT (INOUT) IPOSIT, NCTGT
C
C Local variables
C
      INTEGER :: I, ITG
C
C*** End of declarations rewritten by SPAG
C
      IF(iposit.EQ.0)THEN
         IF(nctarg.GT.0)THEN
            IF(ntgcon.GT.0)READ(megul)iphz, nctgt, notgt, mcont, gucont
            IF(ntgcon.LE.0)READ(megul)iphz
            READ(megul)NOB, NDTRF, NODO, IPOSIT, NOB0, nobl, nob0l
         END IF
      END IF
      IF(iposit.NE.0)THEN
         IF(ntgcon.GT.0)THEN
            READ(megul)iphz
            READ(megul)(nctgt(itg),itg=1,ntgcon)
            READ(megul)(notgt(itg),itg=1,ntgcon)
            READ(megul)(mcont(itg),itg=1,ntgcon)
            READ(megul)(gucont(itg),itg=1,ntgcon)
            WRITE(*,*)'in rdnfto-gt: nctgt=', nctgt
         END IF
         IF(ntgcon.LE.0)THEN
            READ(megul)iphz
            WRITE(*,*)'in rdnfto-le: nctgt=', nctgt
         END IF
         READ(megul)NOB, NDTRF, NODO, IPOSIT, NOB0, nobl, nob0l
         READ(megul)(NOBE(I),I=1,nsym)
         READ(megul)(NOBP(I),I=1,nsym)
         READ(megul)(NOBV(I),I=1,nsym)
      END IF
      RETURN
      END SUBROUTINE RDNFTO
!*==rdwf.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE RDWF(M,ICDO,INDO,NDO,LNDOF,CDO,LCDOF,megul)
c
c     read CSFs as generated and projected by CONGEN (part 2)
C            NDO(LNDOF)  DTRS IN PACKED FORM
C            CDO(LCDOF)  CORRESPONDING COEFFICIENTS
C            INDO(NOCSF+1) INDEX ON DTRS
C            ICDO(NOCSF+1) INDEX ON COEF.
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LCDOF, LNDOF, M, MEGUL
      REAL(KIND=wp), DIMENSION(LCDOF) :: CDO
      INTEGER, DIMENSION(M) :: ICDO, INDO
      INTEGER, DIMENSION(LNDOF) :: NDO
      INTENT (IN) LCDOF, LNDOF, M, MEGUL
      INTENT (OUT) CDO, ICDO, INDO, NDO
C
C*** End of declarations rewritten by SPAG
C
      READ(megul)ICDO, INDO
      READ(megul)NDO
      READ(megul)CDO
      RETURN
      END SUBROUTINE RDWF
!*==readcid.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
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
      CALL CIVIO(NFTW,1,NOCSF,NSTAT,EIG,VEC,NALM,iphase,vec)
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
!*==readcip.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE READCIP(nftw,nciset,nocsf,nstat1,nctgt,mgvn,s,sz,nelt,
     &                   EIG,VEC,iphase,NFT)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE scatci_data, ONLY : MEIG
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: MGVN, NCISET, NCTGT, NELT, NFT, NFTW, NOCSF, NSTAT1
      REAL(KIND=wp) :: S, SZ
      REAL(KIND=wp), DIMENSION(*) :: EIG, VEC
      INTEGER, DIMENSION(*) :: IPHASE
      INTENT (IN) NCTGT, NSTAT1
      INTENT (OUT) MGVN, NELT, S, SZ
      INTENT (INOUT) NOCSF
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CHARGE, XNUC, YNUC, 
     &                                            ZNUC
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: CNAME
      REAL(KIND=wp) :: E0
      INTEGER :: I, NALM, NEIG, NNUC, NREC, NSET, NSTAT
      CHARACTER(LEN=120) :: NAME
C
C*** End of declarations rewritten by SPAG
C
C**********************************************************************
C
C     READCIP reads CI data from unit NFTW in format used for
C     polyatomic targets
C
C**********************************************************************
C
c
c.... Find file
      CALL MOVEP(NFTW,nciset,NALM,0,NFT)
      IF(NALM.NE.0)GO TO 1800
C
      READ(nftw)
      READ(NFTW)Nset, nrec, NAME, nnuc, nocsf, nstat, mgvn, s, sz, nelt, 
     &          e0
      ALLOCATE(xnuc(nnuc),ynuc(nnuc),znuc(nnuc),charge(nnuc),
     &         CNAME(nnuc))
      DO i=1, nnuc
         READ(nftw)cname(i), xnuc(i), ynuc(i), znuc(i), charge(i)
      END DO
      IF(nset.NE.nciset)GO TO 1800
      IF(nstat1.GT.nstat)THEN
         WRITE(nft,2980)nstat1, nstat
 2980    FORMAT(/' Target state number',i4,' requested',
     &          /' Only               ',i4,' on file: STOP')
         STOP
      END IF
C
C     READ CI COEFFICIENTS, ZM changed nstat to nstat1 to ensure only
C     the states requested by CIRMAT are read-in.
C
      CALL CIVIO(NFTW,1,NOCSF,nstat1,EIG,VEC,NALM,iphase,vec)
      IF(NALM.NE.0)GO TO 2900
      IF(NOCSF.EQ.NCTGT)THEN
C
c....   Print summary of target data
         WRITE(NFT,100)NSET, name
         WRITE(NFT,101)nocsf, nstat, nnuc
         WRITE(NFT,120)(cname(i),xnuc(i),ynuc(i),znuc(i),charge(i),i=1,
     &                 nnuc)
         NEIG=MIN(MEIG,nstat1)
         WRITE(NFT,102)(EIG(I)+E0,I=1,NEIG)
      ELSE
c....   Error
         WRITE(NFT,2920)NOCSF, NCTGT
         STOP
      END IF
      DEALLOCATE(xnuc,ynuc,znuc,charge,cname)
      RETURN
c
 1800 WRITE(nft,1804)nftw
 1804 FORMAT(/' CIDATA NOT FOUND ON UNIT',I3)
      STOP
c
 2900 WRITE(NFT,2910)
 2910 FORMAT(/' UNABLE TO GET CI-TARGET VECTOR ')
      STOP
c
 100  FORMAT(' SET',I4,4X,A)
 101  FORMAT(10X,'NOCSF=',I5,4X,'NSTAT=',I5,4X,'NNUC =',I5)
 102  FORMAT(10X,'EIGEN-ENERGIES OF THE REQ. STATES',/(16X,5F20.10))
 120  FORMAT(/' Nuclear data     X         Y         Z       Charge'/
     &       (3x,a8,2x,4F10.6))
 2920 FORMAT(/' HAMILTONIAN TRANSFORMATION DATA INCONSISTENT ',2I10)
c
      END SUBROUTINE READCIP
!*==remx.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE REMX(MXU,NELM,IJ,EMX,NBUF)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: MXU, NBUF, NELM
      REAL(KIND=wp), DIMENSION(NBUF) :: EMX
      INTEGER, DIMENSION(2,NBUF) :: IJ
      INTENT (IN) MXU, NBUF
      INTENT (OUT) EMX, IJ, NELM
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     REMX READS THE FILE OF NONZERO HAMILTONIAN MATRIX ELEMENTS
C
C***********************************************************************
C
C
      READ(MXU)NELM, IJ, EMX
C
      RETURN
C
      END SUBROUTINE REMX
!*==search.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE SEARCH(IUNIT,A,ifail)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER(LEN=8) :: A
      INTEGER :: IFAIL, IUNIT
      INTENT (IN) A, IUNIT
      INTENT (OUT) IFAIL
C
C Local variables
C
      CHARACTER(LEN=8) :: A1
      CHARACTER(LEN=8), DIMENSION(4) :: B
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C     Utility to search a dataset IUNIT for header A where the dataset
C     is assumed to have MOLECULE-SWEDEN convention headers. The header
C     convention is
C
C     '********', '        ', '        ', 'ABCDEFGH'
C
C     where ABCDEFGH is a character string such as ONEELINT etc.
C
C***********************************************************************
C
C
      DATA A1/'********'/
C
 1    READ(IUNIT,END=3)B
      IF(B(1).NE.A1 .OR. B(4).NE.A)GO TO 1
C
      ifail=0
      RETURN
C
C---- Process error condition namely, header not found by end of file
C
 3    ifail=1
      RETURN
C
      END SUBROUTINE SEARCH
!*==stdiag.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE STDIAG(CH,CHC,XMHC,CHR,XMHR,NORB,NORBL,NHC,NHR,THRES,
     &                  IDIAG,MPOS)
C THIS ROUTINE REFLECTS TWO TECHNIQUES FOR LABELLING POSITRON ONE-E
C INTEGRALS: ARRAY ARGUMENT (FOR INTERNAL WORKING) AND NEGATIVE SIGN
C (FOR USE, VIA SORT, IN THE CI STEP)
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints, pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IDIAG, NHC, NHR, NORB, NORBL
      REAL(KIND=wp) :: THRES
      REAL(KIND=wp), DIMENSION(*) :: CH, CHC, CHR
      INTEGER, DIMENSION(2*norbl) :: MPOS
      INTEGER(longint), DIMENSION(2,*) :: XMHC, XMHR
      INTENT (IN) CHR, IDIAG, MPOS, NHR, NORB, NORBL, THRES
      INTENT (OUT) CHC
      INTENT (INOUT) CH, NHC
C
C Local variables
C
      INTEGER, SAVE :: I, N
      INTEGER, DIMENSION(8), SAVE :: ITEMP
C
C*** End of declarations rewritten by SPAG
C
      IF(IDIAG.EQ.0)THEN
         DO I=1, NHR
            CALL unpack8ints(XMHR(1,I),ITEMP)
            N=ITEMP(2)
            IF(itemp(5).NE.0)N=N+NORB
            CH(N)=CH(N)-CHR(I)
         END DO
      END IF
      NHC=0
      DO I=1, NORBL
         IF(ABS(CH(I)).LT.THRES)CYCLE
         NHC=NHC+1
         IF(mpos(2*I-1).NE.0)THEN
            CALL pack8ints(0,i,0,i,1,0,0,0,XMHC(1,NHC))
            CHC(NHC)=CH(I)
         ELSE
            CALL pack8ints(0,i,0,i,0,0,0,0,XMHC(1,NHC))
            CHC(NHC)=CH(I)
         END IF
      END DO
c
      RETURN
      END SUBROUTINE STDIAG
!*==store.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE STORE(NCH0,XNPQP,CPQP,NORB,THRES,CH,NPAM)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NCH0, NORB, NPAM
      REAL(KIND=wp) :: THRES
      REAL(KIND=wp), DIMENSION(*) :: CH, CPQP
      INTEGER(longint), DIMENSION(2,*) :: XNPQP
      INTENT (IN) CH, NORB, NPAM, THRES
      INTENT (OUT) CPQP
      INTENT (INOUT) NCH0
C
C Local variables
C
      INTEGER, SAVE :: I, J, M, MV, NCHP
C
C*** End of declarations rewritten by SPAG
C
      mv=0
      IF(NPAM.EQ.2 .OR. NPAM.EQ.4)mv=1
      M=0
      NCHP=NCH0
      IF(NPAM.LE.2)THEN
         DO I=1, NORB
            DO J=1, I
               M=M+1
               IF(ABS(CH(M)).GE.THRES)THEN
                  NCHP=NCHP+1
                  CPQP(NCHP)=CH(M)
                  CALL pack8ints(I,I,J,J,mv,0,0,0,XNPQP(1,NCHP))
c         NPQP(NCHP)=IPACK4(I,I,J,J)
               END IF
            END DO
         END DO
      ELSE
         DO I=1, NORB
            DO J=1, I
               M=M+1
               IF(ABS(CH(M)).GE.THRES)THEN
                  NCHP=NCHP+1
                  CPQP(NCHP)=CH(M)
                  CALL pack8ints(I,J,I,J,mv,0,0,0,XNPQP(1,NCHP))
c         NPQP(NCHP)=IPACK4(I,J,I,J)
               END IF
            END DO
         END DO
      END IF
      NCH0=NCHP
      RETURN
      END SUBROUTINE STORE
!*==tabla.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE TABLA(NFTW,NSYM,NOB,LPQRS,MPQRS,MBAS,NBASH,NBAS,ICUR,
     &                 MAXRW,SYMTYP,NPFLG,lusme)
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE global_utils, ONLY : INDFUNC
      USE integer_packing, ONLY : unpack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: ICUR, LPQRS, LUSME, MAXRW, NFTW, NPFLG, NSYM, SYMTYP
      INTEGER, DIMENSION(*) :: MBAS, NBAS, NBASH, NOB
      INTEGER(longint), DIMENSION(2,1023) :: MPQRS
      INTENT (IN) LPQRS, LUSME, MPQRS, NFTW, NOB, NPFLG, NSYM, SYMTYP
      INTENT (INOUT) ICUR, MAXRW, MBAS, NBAS, NBASH
C
C Local variables
C
      INTEGER :: I, M, MD, MM, MP, MPA, MPB, MVL, N, NAM1, NAM2, NPQ, 
     &           NPR, NRS
      INTEGER, DIMENSION(8) :: MS
      INTEGER(longint), DIMENSION(2) :: ND
      INTEGER, DIMENSION(4) :: NP
C
C*** End of declarations rewritten by SPAG
C
C*******************************************************************
C
C     This routine generates a set of integral block pointers for the
C     blocks of one and two electron integrals in the problem;
C     information which is then used in the sorting step to place each
C     integral into a box.
C
C     Note that this routine has been modified from the Alchemy I
C     version in order to interface to Sweden-Molecule transformed
C     integrals.
C
C     Input data:
C     SYMTYP(I*4)  .lt. 2 for ALCHEMY case
C                  .ge. 2 for MOLECULE case
C          NPFLG   PRINT SYM INT BLOCK DATA FOR MOLECULE CASES
C           NSYM   Number of symmetries in the problem
C            NOB   Number of orbitals per symmetry
C          LPQRS   Number of elements in the integral block table
C       MPQRS( )   Packed integral block table
C
C   Output data: (* See comments in code *)
C          MAXRW   maximum of NPQ and NRS sets
C       MBAS ( )   LENGTH=2*NSYM+1
C       NBASH( )   LENGTH=NSYM+1
C       NBAS ( )   LENGTH=(NSYM*(NSYM+1))/2+ (NSYM-I+2)*(NSYM-I+1),
C                     I=2,NSYM)
C           ICUR   TOTAL NO OF 2e INTEGRALS+1
C
C***********************************************************************
C
C---- Statement function definition
C
C JMC      ITRIAD(M)=(M*(M+1))/2
C JMC replacing this statement function with references to INDFUNC(M+1,0)
c (see the definition of INDFUNC in the global_utils module).
C
C---- Initialization
C
      MAXRW=0
      DO I=1, 4
         MS(I)=0
         NP(I)=0
      END DO
      MBAS(1)=0
C
C---- Establish how many blocks of two electron integrals there
C     is at most. For Sweden this is straightforward
C     triangulation on symmetry; for Alchemy must look at group
C     degeneracy.
C
      IF(SYMTYP.GE.2)THEN
C
C....... MOLECULE case
C
         MVL=1
         MBAS(2)=INDFUNC(NSYM+1,0)
         MBAS(2)=INDFUNC(MBAS(2)+1,0)
      ELSE
C
C....... ALCHEMY case
C
         MBAS(2)=INDFUNC(NSYM+1,0)
         MVL=2*NSYM-1
         DO M=2, MVL
            MD=(M-2)/2
            MM=NSYM-MD-1
            MBAS(M+1)=MBAS(M)+INDFUNC(MM+1,0)
         END DO
      END IF
C
C---- Establish the block pointers for one electron integrals in
C     the array NBASH. Note that we just triangulate on the number
C     of orbitals per symmetry as we have Hamiltonian type integrals
C     only here. It is useful to note that NBASH points to the first
C     storage location of the block
C
      NBASH(1)=1
      DO M=1, NSYM
         NBASH(M+1)=NBASH(M)+INDFUNC(NOB(M)+1,0)
      END DO
C
      DO M=1, MBAS(MVL+1)
         NBAS(M)=0
      END DO
C
C---- ICUR counts the current next 2e integral index. Here we set it
C
      ICUR=1
C
C---- Now loop over the different blocks of two electron integrals;
C     compute the number of integrals per block and store the pointer
C     index into array NBAS. At the end ICUR will be the total number
C     of integrals (one+two electron) for this orbital set. Note in
C     particular that each block of Molecule-Sweden integrals will
C     contain one extra number, a block descriptor, hence the +1 in
C     line number 170 below.
C
      IF(SYMTYP.GE.2)THEN
C
C....... MOLECULE case
C
         IF(npflg.GT.0)WRITE(6,9911)
 9911    FORMAT(//,5X,'Integral Storage Table follows : ',/)
C
         DO N=1, LPQRS
            ND(1)=MPQRS(1,N)
            ND(2)=MPQRS(2,N)
C
cmc         CALL UNPAK4(ND,MS)
            CALL unpack8ints(ND,MS)
C
            DO I=1, 4
               MS(I)=MS(I)+1
               NP(I)=NOB(MS(I))
            END DO
            IF(MS(1).EQ.MS(2))THEN
               NPQ=INDFUNC(NP(1)+1,0)
            ELSE
               NPQ=NP(1)*NP(2)
            END IF
            NAM1=INDFUNC(MAX(MS(1),MS(2))+1,0)-ABS(MS(1)-MS(2))
C
            IF(MS(3).EQ.MS(4))THEN
               NRS=INDFUNC(NP(3)+1,0)
            ELSE
               NRS=NP(3)*NP(4)
            END IF
            NAM2=INDFUNC(MAX(MS(3),MS(4))+1,0)-ABS(MS(3)-MS(4))
C
            MP=INDFUNC(MAX(NAM1,NAM2)+1,0)-ABS(NAM1-NAM2)
            NBAS(MP)=ICUR
C
            IF(npflg.GT.0)WRITE(6,9912)MS(1), MS(2), MS(3), MS(4), MP, 
     &                                 ICUR
 9912       FORMAT(3X,4I3,1X,I5,1X,I10)
            IF(NAM1.EQ.NAM2)THEN
               NPR=INDFUNC(NPQ+1,0)
            ELSE
               NPR=NPQ*NRS
            END IF
            MAXRW=MAX(MAXRW,NPQ,NRS)
            IF(npr.GT.0)ICUR=ICUR+NPR+1
         END DO
      ELSE
C
C....... ALCHEMY case
C
         DO N=1, LPQRS
            ND(1)=MPQRS(1,N)
            ND(2)=MPQRS(2,N)
            IF(ND(1).GT.0)THEN
               MD=0
            ELSE
               MD=1
               ND(1)=-ND(1)
            END IF
cmc         CALL UNPAK4(ND,MS)
            CALL unpack8ints(ND,MS)
            DO I=1, 4
               NP(I)=NOB(MS(I)+1)
            END DO
C
            IF(MD.EQ.0)THEN
               M=MS(1)+MS(2)
            ELSE
               M=ABS(MS(1)-MS(2))
            END IF
            IF(M.EQ.0)THEN
               MD=-1
            ELSE
               MD=(M-1)/2
            END IF
            MPA=MS(1)-MD
            MPB=MS(3)-MD
            IF(MPA.LT.MPB)THEN
               MP=(MPB*(MPB-1))/2+MPA+MBAS(M+1)
            ELSE
               MP=(MPA*(MPA-1))/2+MPB+MBAS(M+1)
            END IF
C
            NBAS(MP)=ICUR
C
            IF(MS(1).EQ.MS(2))THEN
               NPQ=INDFUNC(NP(1)+1,0)
            ELSE
               NPQ=NP(1)*NP(2)
            END IF
            IF(MS(3).EQ.MS(4))THEN
               NRS=INDFUNC(NP(3)+1,0)
            ELSE
               NRS=NP(3)*NP(4)
            END IF
            MAXRW=MAX(NPQ,NRS,MAXRW)
            IF(MS(1).EQ.MS(3) .AND. MS(2).EQ.MS(4))THEN
               ICUR=ICUR+INDFUNC(NPQ+1,0)
            ELSE
               ICUR=ICUR+NPQ*NRS
            END IF
C
         END DO
      END IF
C
C---- Presently ICUR in the pointer to the start of the next block.
C     Since we are in fact finished with all blocks then ICUR-1 is
C     the total number of integrals for the orbital set in question.
C
      ICUR=ICUR-1
C
C---- Write information to the printer about the blocks
C
      WRITE(NFTW,400)(NOB(I),I=1,NSYM)
      IF(LUSME.NE.0)WRITE(lusme,401)(NOB(I),I=1,NSYM)
 400  FORMAT(' NOBE =',10I10/(7X,10I10))
 401  FORMAT(10I10/(7X,10I10))
      WRITE(NFTW,410)(MBAS(I),I=1,MVL+1)
      IF(LUSME.NE.0)WRITE(lusme,402)(MBAS(I),I=1,MVL+1)
 402  FORMAT(10I10/(7X,10I10))
 410  FORMAT(' MBAS =',10I10/(7X,10I10))
      WRITE(NFTW,420)(NBASH(I),I=1,NSYM+1)
      IF(LUSME.NE.0)WRITE(lusme,403)(NBASH(I),I=1,NSYM+1)
 403  FORMAT(10I10/(7X,10I10))
 420  FORMAT(' NBASH=',10I10/(7X,10I10))
      IF(LUSME.NE.0)WRITE(lusme,440)(NBAS(I),I=1,MBAS(MVL+1)), ICUR
 440  FORMAT(10I10/(7X,10I10))
      IF(SYMTYP.LT.2)THEN
         WRITE(NFTW,430)(NBAS(I),I=1,MBAS(MVL+1))
 430     FORMAT(' NBAS =',10I10/(7X,10I10))
      END IF
C
      RETURN
      END SUBROUTINE TABLA
!*==tablba.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE TABLBA(NFTW,NSYM,NMPRS,LPQRS,NPQRS,MPQRS,SCFUSE)
C
C     ALCHEMY version
C     MODIFICATION OF M. YOSHIMINES ROUTINE TO CREATE A LIST OF
C        TRANSFORMAD INTEGRAL BLOCKS  (PSB 5/74)
C     MODIFICATIONS PERMIT CREATION OF A LIST THAT PLACES SCF
C                                BLOCKS FIRST
C
C     NSYM          NUMBER OF SYMMETRIES
C     NMPRS         LE 0 CREATE LIST
C                   GT 0 LIST NMPRS LONG AS INPUT
C                   HOWEVER IF SYMTYP GE 2 LIST IS ALWAYS CREATED
C     LPQRS(5,2)    INPUT LIST OR TEMP STORAGE FOR ALCHEMY
C                   NOT USED IN MOLECULE
C                   AN I*2 ARRAY
C     NPQRS         LENGTH OF BLOCK LIST
C     MPQRS(2)      PACKED BLOCK INDICES / MPQRS LT 0 THEN M SMALLER
C                   MPQRS(2) GE 0 THEN M GREATER OR ONLY M
C     SCFUSE        EQ 0 CREATE BLOCKS IN USUAL ORDER
C                   NE 0 CREATE BLOCK LIST WITH SCF BLOCKS FIRST
C
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : pack8ints
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NFTW, NMPRS, NPQRS, NSYM, SCFUSE
      INTEGER, DIMENSION(5,*) :: LPQRS
      INTEGER(longint), DIMENSION(2,1023) :: MPQRS
      INTENT (IN) NFTW, NMPRS, NSYM, SCFUSE
      INTENT (OUT) MPQRS
      INTENT (INOUT) LPQRS, NPQRS
C
C Local variables
C
      INTEGER :: I, IBGN, IF, IPASS, J, K, M, MSYM, N, NBGN
      INTEGER, DIMENSION(5) :: NN
      INTEGER(longint), DIMENSION(2) :: NOG
C
C*** End of declarations rewritten by SPAG
C
      NPQRS=NMPRS
      IF(NPQRS.GT.0)GO TO 80
      K=0
      MSYM=NSYM+NSYM-1
      IF(SCFUSE.EQ.0)THEN
C
C        PROCESS ALCHEMY ASSUMING STANDARD LIST
C
         DO M=1, MSYM
            NBGN=M/2+1
            NN(1)=M-1
            DO N=NBGN, NSYM
               NN(2)=N-1
               NN(3)=ABS(N-M)
               DO I=NBGN, N
                  NN(4)=I-1
                  NN(5)=ABS(I-M)
                  K=K+1
                  DO J=1, 5
                     LPQRS(J,K)=NN(J)
                  END DO
               END DO
            END DO
         END DO
      ELSE
C
C        PROCESS ALCHEMY ASSUMING SCF BLOCKS FIRST
C
         DO IPASS=1, 2
            DO M=IPASS, MSYM
               IBGN=M/2+1
               NN(1)=M-1
               DO N=M/2+IPASS, NSYM
                  NN(2)=N-1
                  NN(3)=ABS(N-M)
                  IF=N+1-IPASS
                  IF(M.NE.1 .AND. IPASS.EQ.1)IBGN=IF
                  DO I=IBGN, IF
                     NN(4)=I-1
                     NN(5)=ABS(I-M)
                     K=K+1
                     DO J=1, 5
                        LPQRS(J,K)=NN(J)
                     END DO
                  END DO
               END DO
            END DO
            MSYM=MSYM-2
         END DO
      END IF
C
      NPQRS=K
 80   DO I=1, NPQRS
cmc      NOG=IPACK4(LPQRS(2,I),LPQRS(3,I),LPQRS(4,I),LPQRS(5,I))
         CALL pack8ints(LPQRS(2,I),LPQRS(3,I),LPQRS(4,I),LPQRS(5,I),0,0,
     &                  0,0,NOG(1))
         IF(LPQRS(2,I)+LPQRS(3,I).NE.LPQRS(1,I))NOG(1)=-NOG(1)
         MPQRS(1,I)=NOG(1)
         MPQRS(2,I)=NOG(2)
      END DO
      WRITE(NFTW,100)NPQRS, ((LPQRS(I,J),I=1,5),J=1,NPQRS)
 100  FORMAT(' NMPRS=',I5//' LPQRS'/(2X,5(2X,I3)))
C
      RETURN
      END SUBROUTINE TABLBA
!*==tablbm.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE tablbm(iwrite,nsym,npqrs,mpqrs,npflg)
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE integer_packing, ONLY : unpack8ints, pack8ints
      USE global_utils, ONLY : MPROD
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IWRITE, NPFLG, NPQRS, NSYM
      INTEGER(longint), DIMENSION(2,1023) :: MPQRS
      INTENT (IN) IWRITE, NPFLG, NSYM
      INTENT (INOUT) MPQRS, NPQRS
C
C Local variables
C
      INTEGER :: I, J, JEND, K, L, LEND, MF12, MF34, MFTEST
      INTEGER, DIMENSION(8) :: NN
      INTEGER(longint), DIMENSION(2) :: NX
C
C*** End of declarations rewritten by SPAG
C
c***********************************************************************
c
c     For symmetry groups that are subgroups of D2h counts the number
c     of non-zero blocks of two electron integrals.
c
c     Description:
c
c     In the Abelian group D2h we consider pairs of symmetry types.
c     For N symmetries there are N*(N+1)/2 = NX unique pairs. For
c     two electron integrals (which are pairs of pairs) then there are
c
c                         NX*(NX+1)/2 unique
c
c     combinations. However not all combinations can be valid due to
c     the group theoretical restriction that the integrand must belong
c     to the overall symmetric irreducible representation, for example
c     for C2v (N=4) there are 55 possible integral pairs. This routine
c     examines all the possible pairs and calculates which are valid.
c
c     Note:
c
c     This routine performs the same function as routine TABLBA and
c     TABLBM which were part of the original Alchemy package. The
c     difference here is that the integral block ordering corresponds
c     to that produced by MOLECULE-SWEDEN instead of Alchemy. For
c     reasons of design compatibility this routine has the same argument
c     list as both the corresponding alchemy routines; this leads to
c     some redundant parameters being passed. Also note that in every
c     Sweden integral block there is an extra integral; this real number
c     is used to hold some block description data for that particular
c     block of integrals. Thus when we evaluate the number of integrals
c     per block in this routine we always add one.
c
c     Input data:
c         iwrite  logical unit for the printer
c           nsym  number of symmetries in the orbital set
c
c     Output data:
c           npqrs length of created block list
c           mpqrs packed block indices
c
c***********************************************************************
c
c---- Format statements
c
 500  FORMAT('1',/,10x,'D2h Two Electron Integral Box ','Information')
 510  FORMAT(/,10x,'No. of Boxes = ',i4,/)
 520  FORMAT(10x,'Box Descriptions: ',/)
 530  FORMAT(10x,4I5)
c
c---- Note that symmetry table stores sym number less one !
c
      npqrs=0
c
c---- All IIII blocks of integrals
      WRITE(*,*)'All IIII blocks of integrals'
c
      DO i=1, nsym
         j=i-1
         npqrs=npqrs+1
cmc         mpqrs(npqrs)=ipack4(j,j,j,j)
         CALL pack8ints(j,j,j,j,0,0,0,0,nx)
         mpqrs(1,npqrs)=nx(1)
         mpqrs(2,npqrs)=nx(2)
         WRITE(*,*)npqrs, j, j, j, j
      END DO
c
      IF(nsym.EQ.1)GO TO 150
c
c---- All IJIJ blocks
      WRITE(*,*)'All IJIJ blocks'
c
      DO i=2, nsym
         jend=i-1
         DO j=1, jend
            lend=j-1
            npqrs=npqrs+1
cmc            mpqrs(npqrs)=ipack4(jend,lend,jend,lend)
            CALL pack8ints(jend,lend,jend,lend,0,0,0,0,nx)
            mpqrs(1,npqrs)=nx(1)
            mpqrs(2,npqrs)=nx(2)
            WRITE(*,*)npqrs, jend, lend, jend, lend
         END DO
      END DO
c
c---- The IIJJ blocks of integrals
      WRITE(*,*)'The IIJJ blocks of integrals'
c
      DO i=2, nsym
         jend=i-1
         DO j=1, jend
            lend=j-1
            npqrs=npqrs+1
cmc            mpqrs(npqrs)=ipack4(jend,jend,lend,lend)
            CALL pack8ints(jend,jend,lend,lend,0,0,0,0,nx)
            mpqrs(1,npqrs)=nx(1)
            mpqrs(2,npqrs)=nx(2)
            WRITE(*,*)npqrs, jend, jend, lend, lend
         END DO
      END DO
c
      IF(nsym.LT.4)GO TO 150
c
c---- The IJKL blocks
      WRITE(*,*)'The IJKL blocks'
c
      DO i=2, nsym
         jend=i-1
         DO j=1, jend
            DO k=2, i
               lend=k-1
               IF(i.EQ.k)lend=j
               DO l=1, lend
                  IF(i.EQ.k .AND. j.EQ.l)CYCLE
                  mf12=mprod(i,j,0,iwrite)
                  mf34=mprod(k,l,0,iwrite)
                  mftest=mprod(mf12,mf34,0,iwrite)
                  IF(mftest.NE.1)CYCLE
                  npqrs=npqrs+1
cmc                  mpqrs(npqrs)=ipack4(i-1,j-1,k-1,l-1)
                  CALL pack8ints(i-1,j-1,k-1,l-1,0,0,0,0,nx)
                  mpqrs(1,npqrs)=nx(1)
                  mpqrs(2,npqrs)=nx(2)
                  WRITE(*,*)npqrs, i-1, j-1, k-1, l-1
               END DO
            END DO
         END DO
      END DO
c
c---- Write a summary now to the printer if requested
c
 150  IF(npflg.LE.0)RETURN
      WRITE(iwrite,500)
      WRITE(iwrite,510)npqrs
      WRITE(iwrite,520)
      DO i=1, npqrs
         nx(1)=mpqrs(1,i)
         nx(2)=mpqrs(2,i)
cmc         call unpak4(nx,nn)
         CALL unpack8ints(nx,nn)
         WRITE(iwrite,530)(nn(j)+1,j=1,4)
      END DO
c
      RETURN
      END SUBROUTINE TABLBM
!*==tpinde.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE TPINDE(IODR,NRI,NSM,NOB,MBAS,NBAS,NBASH,IPAIR,KT,XNIP,
     &                  CFP,ndim,xint1e,nint1e,xint2e,elem,isymtp,map,
     &                  NCORB,mcorb,nci,ukrmolp_ints)
C
C***********************************************************************
C
C     Computes a sequential index for each formula element which is used
C     to access the integral arrays XINT1e and XINT2e.
C     Loops over ELEM allow for case of CI target contraction with
C     more than one target state for a given symmetry.
C     This routine is for C-infinity-v symmetry (ALCHEMY) and SWEDEN.
C
C     ZM: modified to work with UKRMol+ integrals. This is controlled by 
C         the input flag ukrmolp_ints.
C
C     KT      NUMBER OF INDICES PASSED
C     NIP(KT)  PACKED INDICES / ONE PER BYTE
C              H-INTEGRALS (0,NP,0,NQ)
C                  NO ASSUMPTION MADE ON ORDER OF NP AND NQ
C              2-INTEGRALS (NP,NQ,NR,NS)
C                  ASSUME MP GE MQ, MR GE MS, MPQ GE MRS
C                  NOTHING FURTHER
C
C     IODR = 0, PQRS AND =I, RSPQ (for ALCHEMY) (usual IODR=1)
C     IODR = 0 or 1 for SWEDEN
C
C     IPAIR(I)=I*(I-1)/2 ADDED TO M. YOSHIMINE'S ROUTINE
C
C**********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE precisn, ONLY : longint ! for specifying 64-bit integers
      USE consts, ONLY : ZERO=>XZERO
      USE integer_packing, ONLY : unpack8ints
      USE ukrmol_interface_gbl, ONLY : GET_INTEGRAL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IODR, ISYMTP, KT, NCI, NCORB, NDIM, NINT1E
      REAL(KIND=wp), DIMENSION(ndim,nci) :: CFP
      REAL(KIND=wp), DIMENSION(nci) :: ELEM
      INTEGER, DIMENSION(0:*) :: IPAIR
      INTEGER, DIMENSION(*) :: MAP, MBAS, MCORB, NBAS, NBASH, NOB, NRI, 
     &                         NSM
      REAL(KIND=wp), ALLOCATABLE :: XINT1E(:), XINT2E(:)
      INTEGER(longint), DIMENSION(2,kt) :: XNIP
      LOGICAL :: ukrmolp_ints
      INTENT (IN) CFP, ISYMTP, KT, MCORB, NCI, NCORB, NDIM, NINT1E, 
     &            XINT1E, XINT2E, ukrmolp_ints
      INTENT (INOUT) ELEM
C
C Local variables
C
      INTEGER :: I, IA, IB, IDX, IQ, NK
      INTEGER, DIMENSION(4) :: mapped
      INTEGER, DIMENSION(8) :: LWD
      REAL(KIND=wp) :: THRESH, X
C
C*** End of declarations rewritten by SPAG
C
C---- Loop to 200 is over formula elements. These are assigned a
C     sequential index which is calculated from the orbital labels.
C
C      Two possibilities exist here:
C
C   1) Two electron integrals (treated differently by Alchemy & Sweden & UKRMol+)
C   2) One electron integrals (treated the same by Alchemy & Sweden but differently in UKRMol+)
C
c
ccc      write(6,*) '*****************START PINDEX***************'
      DO i=1, nci
         elem(i)=zero
      END DO
c
      DO NK=1, KT
C
C..... Unpack the orbital indices into array LWD; then scan the
C      quartet for the occurence of the dummy continuum orbital.
C      Should this be found, set the coefficient to zero and then
C      proceed to sorting as normal. In the case that no dummy
C      continuum is being used (NCORB=-1), skip the search and
C      start ordering immediately.
C
         CALL unpack8ints(XNIP(1,NK),LWD)
c
         IF(NCORB.GT.0)THEN
            DO IQ=1, 4
               IF(LWD(IQ).EQ.NCORB)GO TO 200
            END DO
         ELSE IF(NCORB.EQ.0)THEN
            DO IQ=1, 4
               IF(lwd(iq).LE.0)CYCLE
               IF(mcorb(LWD(IQ)).EQ.0)GO TO 200
            END DO
         END IF
C
         IF(LWD(1).EQ.0)THEN
c        1-electron case
            IA=map(lwd(2))
            IB=map(lwd(4))
            if (isymtp.ge.2 .and. ukrmolp_ints) then !branch to the UKRMol+ integral
               x=get_integral(ia,ib,0,0,lwd(5))
            else !SWEDEN/ALCHEMY
               idx=nind1e(nsm(ia),nri(ia),nri(ib),ipair,nbash)
c              positron case
               IF(lwd(5).NE.0)idx=idx+nint1e
               x=xint1e(idx)
            endif
         ELSE IF(isymtp.LE.1)THEN
c        2-electron  ALCHEMY case
            idx=nin2ea(lwd,nri,nsm,nob,ipair,mbas,nbas,map,iodr)
            x=xint2e(idx)
         ELSE IF (isymtp.ge.2 .and. (.not. ukrmolp_ints))THEN
c        2-electron SWEDEN case
            idx=nin2em(lwd,nri,nsm,nob,ipair,nbas,iodr,map)
            x=xint2e(idx)
         ELSE IF (isymtp.ge.2 .and. ukrmolp_ints)THEN
c        2-electron UKRMol+ case
            mapped(1) = map(LWD(1))
            mapped(2) = map(LWD(2))
            mapped(3) = map(LWD(3))
            mapped(4) = map(LWD(4))
            x=get_integral(mapped(1),mapped(2),mapped(3),mapped(4),
     &                     lwd(5))
         END IF
         THRESH=1.0E-8_wp
         DO i=1, nci
            IF(ABS(cfp(nk,i)).GT.THRESH)THEN
               elem(i)=elem(i)+cfp(nk,i)*x
            END IF
         END DO
 200  END DO
      RETURN
      END SUBROUTINE TPINDE
!*==wrtem.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
 
 
      SUBROUTINE WRTEM(NFTE,NT,LEMBF,NEMBF,EM)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LEMBF, NFTE, NT
      REAL(KIND=wp), DIMENSION(LEMBF) :: EM
      INTEGER, DIMENSION(2,LEMBF) :: NEMBF
      INTENT (IN) EM, LEMBF, NEMBF, NFTE, NT
C
C*** End of declarations rewritten by SPAG
C
C**********************************************************************
C
C     WRTEM WRITES THE ENERGY MATRIX EM TO THE DISK FILE ASSOCIATED
C          WITH UNIT NFTE
C
C**********************************************************************
C
C
      WRITE(NFTE)NT, NEMBF, EM
C
      RETURN
      END SUBROUTINE WRTEM
!*==write_svn_info.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
 
 
 
C
C     Write svn info for file to unit
C
      SUBROUTINE WRITE_SVN_INFO(unit)
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: UNIT
      INTENT (IN) UNIT
C
C*** End of declarations rewritten by SPAG
C
C
C     Please do not edit text between pair of $'s
C     subversion automatically edits these entries.
C
C     The double colon (::) specifies fixed width which
C     is important to work with fortran fixed format files.
C
C     The format statement has been carefully constructed to
C     give maximum width between $$'s
C
C     If edits are not occuring then enable svn keywords with:
C     svn propset svn:keywords "URL Author Date Rev Id" thisfile.f
C
 9876 FORMAT(/,1X,                                                     '
     &$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     &'/,1X,'$',64X,'$',/,
     &     1X,'$',16X,'Subversion Revision Information',17X,'$',/,
     &     1X,'$',64X,'$',/,1X,                                        '
     &$Id::                                                            $
     &',/,1X,'$',64X,'$',/,1X,                                         '
     &$URL::                                                           $
     &',/,1X,                                                          '
     &$LastChangedBy::                                                 $
     &',/,1X,                                                          '
     &$LastChangedDate::                                               $
     &',/,1X,                                                          '
     &$LastChangedRevision::                                           $
     &',/,1X,'$',64X,'$',/,1X,                                         '
     &$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     &',/)
C
      WRITE(unit,9876)
C
      END SUBROUTINE WRITE_SVN_INFO

      END MODULE SCATCI_ROUTINES
