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
!*==dgmain.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE dgmain(ntgsym,numtgt,notgs,mcont,nfte,nelm,symtyp,s,sz,
     &                  nelt,mgvn,nfti,kphz,size_kphz,ukrmolp_ints)
!***********************************************************************
!
!     Hamiltonian is read, energy shifted if required and
!     diagonalized to yield eigenvalues and eigenvectors.
!
!***********************************************************************
      USE precisn ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      USE params, ONLY : STAR=>C8STARS, CBLANK
      USE scatci_data, ONLY : NTGTMX
      USE SCATCI_ROUTINES, ONLY: movew, prthd, civio
#ifdef arpack
      use m_4_mkarp            ! pgg
#endif
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: MGVN, NELT, NFTE, NFTI, NTGSYM, SYMTYP, size_kphz
      REAL(KIND=wp) :: S, SZ
      INTEGER, DIMENSION(ntgsym) :: MCONT
      INTEGER, DIMENSION(NTGTMX) :: NOTGS, NUMTGT
      integer, dimension(size_kphz) :: kphz
      LOGICAL :: ukrmolp_ints
      INTENT (IN) NUMTGT, SYMTYP, ukrmolp_ints
!
! Local variables
!
      REAL(KIND=wp) :: CRITC, CRITE, CRITR, EN, EONE=XZERO, ORTHO, 
     &                THRHM, THRPRT, TOL, XREF
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: DGEM, EIG, HMT
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp), DIMENSION(ntgtmx) :: ESHIFT
      INTEGER :: I, IGH, ITGT, LARGE, LEMBF, M, MINDAV, MSHIFT, 
     &           N, NCIPFG, NCISET, NFTA, NFTW, NHDIM, NKEY, NNUC, 
     &           NOCSF, NOPVEC, NPCVC, NPVEC, NSHIFT, NSTAT, NSYM, NTGT

      INTEGER :: maxiter, IERR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPHZ
      INTEGER, DIMENSION(20) :: KEYCSF, NHE
      character(len=120) :: NAME
      character(len=4), DIMENSION(30) :: NAMP
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(NTGTMX) :: NOTGT
      INTEGER, DIMENSION(2) :: NPFLG
      integer :: i1, i2
      integer :: nelm                   !mal 10/06/2011 : Note that this variable 'nelm'
                                        !is not passed to mkarp. Therefore the value
                                        !for the number of non-zero elements in the N+1
                                        !Hamiltonian is obtained inside mkarp by reading
                                        !the fort.26 file. Inside mkarp 'nelm' declared 
                                        !as integer(kind=longint). This allows for large
                                        !values for the number of non-zero elements even
                                        !without compiling with the -i8 flag. 
!
!*** End of declarations rewritten by SPAG
!
      EQUIVALENCE(NHD(2),NOCSF)
      EQUIVALENCE(NHD(4),NSYM)
      EQUIVALENCE(NHD(8),NHDIM)
      EQUIVALENCE(NHD(9),NNUC)
      EQUIVALENCE(DTNUC(1),EN)
!
!     NAMELIST DEFINITION :
!
      NAMELIST /CINORN/ NFTA, NFTW, NCISET, NPFLG, NPCVC, NAME, NSTAT, 
     &   NOPVEC, NOTGT, NCIPFG, NKEY, KEYCSF, LARGE, THRPRT, NTGT, 
     &   ESHIFT, ITGT, CRITE, CRITC, CRITR, ORTHO, maxiter, igh
!
      DATA THRHM/1.0E-15_wp/, NFTA/6/, NSTAT/0/, mindav/750/
      DATA NFTW/25/, NCISET/1/, NPFLG/2*0/, NPCVC/1/, NOPVEC/0/, 
     &     ESHIFT/NTGTMX*XZERO/, NHD/10*0/, itgt/0/
      DATA TOL/1.E-12_wp/
      DATA NCIPFG/3/, NKEY/0/, LARGE/8/, THRPRT/.0025_wp/
      DATA CRITE/0.0_wp/, CRITC/1D-10/, CRITR/1D-8/, ORTHO/1D-7/, 
     &     igh/-99/, maxiter/-1/
!
!MAL 13/06/2011 : Default fault value of CRITE has been changed from 1D-12 based on recommendation
!from PGG who found 0.0_wp to be the best value for arpack.

      write(*,*) "==================================================="
      write(*,*) " inside subroutine --> 'dgmain' --> begins. " 
      write(*,*)

      NAME=CBLANK
      ntgt=0
      DO n=1, ntgsym
         DO m=1, numtgt(n)
            ntgt=ntgt+1
            notgt(ntgt)=notgs(n)
         END DO
      END DO
!
!---- Read namelist input...
!
      READ(5,CINORN,END=9998)
!
 9998 CONTINUE
      NKEY=ABS(NKEY)
      LARGE=ABS(LARGE)
      THRPRT=ABS(THRPRT)
!
      WRITE(NFTA,1239)(STAR,I=1,15)
      WRITE(NFTA,2010)NAME
 2010 FORMAT(/' MATRIX DIAGONALIZATION SECTION'/' (RUN)',5X,A)
!
      IF(NTGT.GT.NTGTMX)GO TO 2910
!
      WRITE(NFTA,124)NFTE, NFTW
 124  FORMAT('  DATA SET FILE NUMBERS AND LENGTH SPECIFICATIONS :',/,5X,
     &       ' NFTE =',I3,6X,'NFTW =',I7/)
!
      WRITE(NFTA,126)NPFLG, NCISET
 126  FORMAT(' NPFLG : ',2I6,/,' NCISET: ',I6)
!
      WRITE(NFTA,200)NPCVC, THRHM
 200  FORMAT(/,' NPCVC =',I5,/,' THRHM =',D9.2,/)
!
!     Initializes target energy shifting
      NSHIFT=0
      MSHIFT=0
      DO I=1, NTGT
         IF(ABS(ESHIFT(I)).GT.TOL)THEN
            IF(NOTGT(I).GT.0)THEN
               NSHIFT=NSHIFT+1
            ELSE IF(itgt.GT.0)THEN
               MSHIFT=MSHIFT+1
            ELSE
               WRITE(nfta,247)I, itgt
 247           FORMAT(/' WARNING: ESHIFT.NE.0 BUT NOTGT.EQ.0 FOR STATE',
     &                I3,2X,' and ITGT =',i3,
     &                /' HAMILTONIAN WILL NOT BE MODIFIED')
            END IF
         END IF
      END DO
!
!     READ HEADER ON ENERGY MATRIX FILE
!
      REWIND NFTE
      READ(NFTE)NOCSF, LEMBF, NHD, NAMP, NHE, DTNUC
!
      WRITE(NFTA,2106)NAMP
      NHDIM=Nocsf
      ALLOCATE(iphz(nocsf),dgem(nocsf))
      iphz=0
      IF(nstat.LE.0 .OR. nstat.GT.nocsf)nstat=nocsf
      nhd(3)=nstat
      nhd(5)=ntgt
      IF(nstat.LT.nocsf)THEN
         DO n=1, ntgt
            iphz(n)=notgt(n)
         END DO
      END IF
!
!     HAMILTONIAN CONSTRUCTION AND DIAGONALIZATION
!
      WRITE(NFTA,248)nstat, NHDIM
 248  FORMAT(/I6,' Eigenvectors to be computed for',/I6,
     &       ' dimensional final Hamiltonian matrix')
      IF(igh.GT.0)WRITE(nfta,249)
 249  FORMAT(/' Givens-Householder diagonalisation forced')
!
!     Which diagonalisation method is to be used?
!
      IF(igh.LT.-10)THEN
         IF(nocsf/nstat.LT.5 .OR. nocsf.LT.mindav)THEN
            igh=1  ! givens householder
         ELSE IF(nstat.LE.3)THEN
            igh=0 ! davidson
         ELSE
#ifdef arpack
            igh=-1  ! arpack
#else
            WRITE(NFTA,'("Warning: Arpack not available")')
            igh=1
#endif
         END IF
      END IF

      IF(igh.EQ.0)THEN

!     Use Davidson for larger matrices when less 3 eigenvectors required.
         ALLOCATE(hmt(nhdim*nstat),eig(nhdim))
         CALL MKDVM(hmt,eig,dgem,nhdim,nstat,nfta,nfte,lembf,nelm,CRITE,
     &              CRITC,CRITR,ORTHO,maxiter,IERR)
         IF(IERR.NE.0) STOP 1
#ifdef arpack
      ELSE IF(igh.LT.0)THEN

!     Use ARPACK for larger matrices when less than 20% of
!     eigenvectors required.

!====================================================================!
! pgg --> 1 --> begin.

! 1. Allocating ram space for eigenpairs.

         allocate ( hmt(nhdim*nstat) ) 
         allocate ( eig(nhdim) )

! 2. Calling the "mkarp" to diagonalize the produced
!    and written to the hard disk "fort.26" file.

         call mkarp(nhdim,nstat,maxiter,crite,nfte)

! 3. Opening the produced file "eigensystem_arpack_file.bin" file.
!    This file contains the produced eigenpairs and it is
!    written to the hard disk.

         open( unit=11,
     1         file="___tmp_eigensys", 
     1         status="old", 
     1         action="read",
     1         form="unformatted" )
         
! 4. Reading the eigenvalues.

         do i1 = 1, nstat
           read(unit=11) eig(i1)
         end do

! 5. Reading the eigenvectors.

         do i1 = 1, nstat
           do i2 = 1, nhdim
             read(unit=11) hmt( (i1-1)*nhdim + i2 )
           end do
         end do
  
! 6. Closing the unit and deleting the file 
!    "___tmp_eigensys".

         close( unit=11, status="delete" )

! pgg --> 1 --> end.
!====================================================================!

         REWIND nfte
         READ(nfte)
         CALL ldham(dgem,hmt,iphz,iphz,nocsf,nfte,lembf,0,nfta)
#endif
      ELSE
!        Or use in core Givens-Householder method
         ALLOCATE(hmt(nhdim*nhdim),eig(nhdim))
!
         CALL HMAT(hmt,NTGT,NOTGT,ESHIFT,nshift,NHDIM,NFTA,nfte,lembf,
     &             eone)
!
! The two writes below to files 12345 and 123456 may be used in the
! future
!        write(123456) NHDIM,hmt ! AH test write hamiltonian to file
         IF(npflg(1).NE.0)CALL PRTHAM(hmt,NHDIM,NFTA,npflg(1))
         CALL QLDIAG(NHDIM,HMT,EIG)
!        write(12345,*) EIG ! AH test write hamiltonian to file
         IF(nstat.LT.nocsf)THEN
            REWIND nfte
            READ(nfte)
            CALL ldham(dgem,hmt,iphz,iphz,nocsf,nfte,lembf,0,nfta)
         END IF
      END IF
!
!     PRINTOUT OF RESULTS FROM H-DIAGONALIZATION
!
!     REFERENCE ENERGY
      XREF=EONE+EN
!     Perform energy shifts on final eigenvalues if requested
      IF(mshift.GT.0)THEN
         WRITE(nfta,246)ntgt, (eshift(i),i=1,ntgt)
 246     FORMAT(/' Lowest',i3,' eigen energies shifted'/' ESHIFT',
     &          10F10.5/(7X,10F10.5))
         DO i=1, ntgt
            eig(i)=eig(i)+eshift(i)
         END DO
      END IF
!
      IF(NOPVEC.EQ.-1)THEN
         NPVEC=NSTAT
      ELSE
         NPVEC=NOPVEC
      END IF
      IF(NPVEC.GT.0)CALL PRCI(NOCSF,NPVEC,EIG,EN,XREF,HMT,NCIPFG,NKEY,
     &                        KEYCSF,LARGE,THRPRT,NFTA)
!
!     WRITE CI COEFFICIENTS TO LIBRARY FILE
!
!     A. Harvey 15.5.14 make sure we save the phases.
!     Z.M. iphz and kphz can have different dimensions so we have to make the assignment below correctly:
      i = min(size_kphz,nocsf)
      iphz(1:i)=kphz(1:i)

      IF(NCISET.GE.0)THEN
         IF(symtyp.GE.2)THEN
            CALL writcip(nftw,nciset,name,en,nocsf,nstat,mgvn,s,sz,nelt,
     &                   eig,hmt,iphz,dgem,ntgsym,mcont,notgs,nfti,nfta,
     &                   npflg(2),npcvc,ukrmolp_ints)
         ELSE
            CALL writcid(nftw,nciset,name,nhe,nhd,dtnuc,nocsf,nstat,eig,
     &                   hmt,iphz,dgem,nfta,npflg(2),npcvc)
         END IF
      END IF
      DEALLOCATE(iphz,dgem,hmt,eig)
!
      RETURN
 2910 WRITE(NFTA,2991)NTGT, NTGTMX
 2991 FORMAT(/' NTGT =',I3,' EXCEEDS MAXIMUM VALUE, NTGTMX =',I3)
!
      RETURN
!
 2106 FORMAT(/' (EMX)',5X,30A4)
 1239 FORMAT(/2X,15A8)
!

      write(*,*) " "
      write(*,*) " inside subroutine --> 'dgmain' --> ends. "
      write(*,*) "==================================================="
      write(*,*) " "


      END SUBROUTINE DGMAIN

!*==dgtarg.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE DGTARG(nfte,nhdim,nstat,eigs,vecs,itgt,nelm,CRITE,
     &                  CRITC,CRITR,ORTHO,maxiter,igh,symtyp)
C
C***********************************************************************
C
C     Target Hamiltonian is read and diagonalized
C     to yield nstat eigenvalues and eigenvectors.
C
C***********************************************************************
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE scatci_data, ONLY : MEIG
      USE SCATCI_ROUTINES, ONLY: prthd
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: NTDUM=1
C
C Dummy arguments
C
      REAL(KIND=wp) :: CRITC, CRITE, CRITR, ORTHO
      INTEGER :: IGH, ITGT, MAXITER, NELM, NFTE, NHDIM, NSTAT, SYMTYP
      REAL(KIND=wp), DIMENSION(nstat) :: EIGS
      REAL(KIND=wp), DIMENSION(nhdim,nstat) :: VECS
      INTENT (IN) IGH, SYMTYP
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: DGEM, EIG
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp) :: EN, EONE
      REAL(KIND=wp), DIMENSION(ntdum) :: ESHIFT1
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: HMT
      INTEGER :: I, LEMBF, MINDAV, NEIG, NFTA, NNUC, NOCSF, NOCSF1, 
     &           NPFLG, NSYM, NTGT, IERR
      character(len=120) :: NAMP
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(20) :: NHE
      INTEGER, DIMENSION(ntdum) :: NOTGT1
c      dimension eigs(nstat),vecs(nhdim,nstat),hmt(nhdim,nhdim),
c     * eig(nhdim)
C
C*** End of declarations rewritten by SPAG
C
      EQUIVALENCE(NHD(2),NOCSF)
      EQUIVALENCE(NHD(4),NSYM)
      EQUIVALENCE(NHD(9),NNUC)
      EQUIVALENCE(DTNUC(1),EN)
C
C
      DATA NPFLG/0/, NFTA/6/, mindav/750/
C
C     READ HEADER ON ENERGY MATRIX FILE
C
      ALLOCATE(hmt(nhdim,nhdim),eig(nhdim),dgem(nhdim))
      REWIND NFTE
      READ(NFTE)NOCSF1, LEMBF, nhd, namp, nhe, dtnuc
C
      nocsf=nhdim
C
C     HAMILTONIAN CONSTRUCTION AND DIAGONALIZATION
C
      WRITE(NFTA,248)nstat, NHDIM
 248  FORMAT(/I6,' Eigenvectors to be computed for',/I6,
     &       ' dimensional target Hamiltonian matrix')
      IF(igh.GT.0)WRITE(nfta,249)
 249  FORMAT(/' Givens-Householder diagonalisation forced')
C
c     Which diagonalisation method is to be used?
c
      IF(nocsf/nstat.GT.5 .AND. nocsf.GT.mindav .AND. igh.LE.0)THEN
c        Use Davidson for larger matrices when less than 20% of
c        eigenvectors required.
         CALL MKDVM(hmt,eig,dgem,nhdim,nstat,nfta,nfte,lembf,nelm,CRITE,
     &              CRITC,CRITR,ORTHO,maxiter,IERR)
         IF(IERR.NE.0) STOP 1
      ELSE
c
         ntgt=1
         CALL HMAT(hmt,NTGT,NOTGT1,ESHIFT1,0,NHDIM,NFTA,nfte,lembf,eone)
         IF(npflg.NE.0)CALL PRTHAM(hmt,NHDIM,NFTA,npflg)
         CALL QLDIAG(NHDIM,hmt,EIG)
      END IF
C
C     PRINTOUT/STORE  RESULTS FROM H-DIAGONALIZATION
C
      IF(symtyp.LE.1)THEN
         CALL PRTHD(ITGT,NHD,NAMP,NHE,DTNUC,NSTAT,EIG,NFTA)
      ELSE
c...  Print summary of output data for polyatomic case
         WRITE(NFTA,100)itgt, namp
         WRITE(NFTA,101)nocsf, nstat, en
         NEIG=MIN(MEIG,nstat)
         WRITE(NFTA,102)(EIG(I)+En,I=1,NEIG)
      END IF
c
 100  FORMAT(/' SET',I4,4X,a)
 101  FORMAT(/' NOCSF=',I5,4X,'NSTAT=',I5,4X,'EN   =',F20.10)
 102  FORMAT(/' EIGEN-ENERGIES',/(16X,5F20.10))
c
      CALL MVCIVC(NOCSF,NSTAT,eigs,vecs,EIG,hmt)
      WRITE(NFTA,2182)ITGT, NSTAT
 2182 FORMAT(' STATE',I3,':',I3,' transformation vectors selected')
c     reposition Hamiltonian scratch file
      OPEN(UNIT=nfte,FORM='unformatted')
      READ(NFTE)
C
      DEALLOCATE(eig,hmt,dgem)
      RETURN
      END SUBROUTINE DGTARG
!*==geth.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE GETH(H,ndim,emx,IJ2,nbuf,nfte,NFT)
C
C***********************************************************************
C
C     GETH READS THE NONZERO ENERGY MATRIX ELEMENTS AND CONSTRUCTS
C          THE LOWER TRIANGULAR THE HAMILTONIAN  MATRIX
C
C***********************************************************************
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE SCATCI_ROUTINES, ONLY: REMX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NBUF, NDIM, NFT, NFTE
      REAL(KIND=wp), DIMENSION(nbuf) :: EMX
      REAL(KIND=wp), DIMENSION(ndim,ndim) :: H
      INTEGER, DIMENSION(2,nbuf) :: IJ2
      INTENT (IN) NDIM, NFT
      INTENT (OUT) H
C
C Local variables
C
      INTEGER :: L, LREAD, NELM
C
C*** End of declarations rewritten by SPAG
C
      lread=0
      L=0
      nelm=0
 311  IF(L.EQ.NELM)THEN
         CALL REMX(nfte,NELM,IJ2,EMX,NBUF)
         lread=lread+nelm
         IF(NELM.EQ.0)GO TO 312
         L=0
      END IF
      L=L+1
      IF(max(ij2(1,l),ij2(2,l)).GT.ndim)WRITE(NFT,*)
     &   'HELP: attemp to create element', ij2(1,l), ij2(2,l)
      H(ij2(1,l),ij2(2,l))=EMX(L)
      GO TO 311
C
 312  CLOSE(nfte)
      WRITE(NFT,100)lread, NFTE
 100  FORMAT(/,I10,' Hamiltonian elements read from unit NFTE =',i3)
C
      RETURN
      END SUBROUTINE GETH
!*==hmat.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE HMAT(EMX,NTGT,NOTGT,ESHIFT,nshift,NHDIM,NFT,nfte,nbuf,
     &                eone)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : ZERO=>XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: EONE
      INTEGER :: NBUF, NFT, NFTE, NHDIM, NSHIFT, NTGT
      REAL(KIND=wp), DIMENSION(nhdim,nhdim) :: EMX
      REAL(KIND=wp), DIMENSION(ntgt) :: ESHIFT
      INTEGER, DIMENSION(ntgt) :: NOTGT
      INTENT (IN) ESHIFT, NOTGT, NSHIFT, NTGT
      INTENT (OUT) EONE
      INTENT (INOUT) EMX
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(NBUF) :: BUF
      REAL(KIND=wp) :: ES
      INTEGER :: I, J, K
      INTEGER, DIMENSION(2,NBUF) :: IJ
C
C*** End of declarations rewritten by SPAG
C
C
C*********************************************************************
C
C     HMAT CONSTRUCTS FULL HAMILTONIAN MATRIX... UPPER TRIANGLE BUT
C     INCLUDING ZERO ELEMENTS for Givens-Householder diagonalisation
C
C*********************************************************************
C
C     READ ENERGY MATRIX ELEMENTS
C
 
      emx=zero
 
      CALL geth(emx,nhdim,BUF,ij,nbuf,nfte,NFT)
C
      IF(NSHIFT.GT.0)THEN
         WRITE(NFT,246)NSHIFT, (ESHIFT(I),I=1,NTGT)
         K=0
         DO I=1, NTGT
            ES=ESHIFT(I)
            DO J=1, NOTGT(I)
               k=k+1
               emx(K,k)=emx(K,k)+ES
            END DO
         END DO
      END IF
      eone=emx(1,1)
C
      RETURN
C
 246  FORMAT(/' APPLYING THRESHOLD SHIFTS TO',I3,
     &       ' TARGET STATES'/' ESHIFT',10F10.5/(7X,10F10.5))
      END SUBROUTINE HMAT
!*==ldhama.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE LDHAMA(elem,lj1,lj2,ndim,nfte,nbuf,nfta)
C
C***********************************************************************
C
C     LDHAMA READS THE NONZERO ENERGY MATRIX ELEMENTS
C     and, if room (INCORE=1) stores elements in elem
c     in preparation for Arpack diagonalisation
C***********************************************************************
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE SCATCI_ROUTINES, ONLY: REMX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NBUF, NDIM, NFTA, NFTE
      REAL(KIND=wp), DIMENSION(*) :: ELEM
      INTEGER, DIMENSION(*) :: LJ1, LJ2
      INTENT (IN) NDIM, NFTA
      INTENT (OUT) ELEM, LJ1, LJ2
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(nbuf) :: EMX
      INTEGER, DIMENSION(2,nbuf) :: IJ2
      INTEGER :: L, LOFF, NELM
C
C*** End of declarations rewritten by SPAG
C
      L=0
      nelm=0
      loff=0
 311  IF(L.EQ.NELM)THEN
         CALL REMX(nfte,NELM,IJ2,EMX,NBUF)
         IF(NELM.EQ.0)GO TO 312
         L=0
      END IF
      L=L+1
      IF(max(ij2(1,l),ij2(2,l)).GT.ndim .OR. min(ij2(1,l),ij2(2,l))
     &   .LT.1)THEN
         WRITE(NFTA,*)' HELP: attempt to create element', ij2(1,l), 
     &                ij2(2,l), loff, l
         STOP
      END IF
      loff=loff+1
      elem(loff)=EMX(L)
      lj1(loff)=ij2(1,l)
      lj2(loff)=ij2(2,l)
      GO TO 311
C
 312  CLOSE(nfte)
      WRITE(NFTA,100)loff, NFTE
 100  FORMAT(/,I10,' Hamiltonian elements read from unit NFTE =',i3)
C
      RETURN
      END SUBROUTINE LDHAMA
!*==ldhamb.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE LDHAMB(ndim,nfte,nbuf,nfta,loff)
C
C***********************************************************************
c     LDHAMB counts the number of nonzero elements of the Hamiltonian
c     matrix. This is needed in case of a restart of the ARPACK
c     diagonalization as we skip over the Hamiltonian matrix contruction
c     phase. The Hamiltonian matrix recalculation can be avoided in case of
c     an ARPACK restart by setting ICIDG=2 in the SCATCI &INPUT namelist.
c     Date: 4/09/2008, Amar Dora.
 
C***********************************************************************
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE SCATCI_ROUTINES, ONLY: REMX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: LOFF, NBUF, NDIM, NFTA, NFTE
      INTENT (IN) NDIM, NFTA
      INTENT (INOUT) LOFF
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(nbuf) :: EMX
      INTEGER, DIMENSION(2,nbuf) :: IJ2
      INTEGER :: L, NELM
C
C*** End of declarations rewritten by SPAG
C
      L=0
      nelm=0
      loff=0
 311  IF(L.EQ.NELM)THEN
         CALL REMX(nfte,NELM,IJ2,EMX,NBUF)
         IF(NELM.EQ.0)GO TO 312
         L=0
      END IF
      L=L+1
      IF(max(ij2(1,l),ij2(2,l)).GT.ndim .OR. min(ij2(1,l),ij2(2,l))
     &   .LT.1)THEN
         WRITE(NFTA,*)' HELP: attempt to create element', ij2(1,l), 
     &                ij2(2,l), loff, l
         STOP
      END IF
      loff=loff+1
      GO TO 311
C
 312  CONTINUE
      REWIND nfte
      READ(nfte)
C
      RETURN
      END SUBROUTINE LDHAMB
!*==ldham.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE LDHAM(diag,offdg,lj1,lj2,ndim,nfte,nbuf,incore,nfta)
C
C***********************************************************************
C
C     LDHAM READS THE NONZERO ENERGY MATRIX ELEMENTS
C     Stores diagonals in DIAG
C     and, if room (INCORE=1) stores off-diagonals in OFFDG
c     in preparation for Davidson diagonalisation
C***********************************************************************
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE SCATCI_ROUTINES, ONLY: REMX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: INCORE, NBUF, NDIM, NFTA, NFTE
      REAL(KIND=wp), DIMENSION(*) :: DIAG, OFFDG
      INTEGER, DIMENSION(*) :: LJ1, LJ2
      INTENT (IN) INCORE, NDIM, NFTA
      INTENT (OUT) DIAG, LJ1, LJ2, OFFDG
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(nbuf) :: EMX
      INTEGER, DIMENSION(2,nbuf) :: IJ2
      INTEGER :: L, LOFF, LREAD, NELM
C
C*** End of declarations rewritten by SPAG
C
      IF(incore.EQ.0)lj2(1)=nfte
      lread=0
      L=0
      nelm=0
      loff=0
 311  IF(L.EQ.NELM)THEN
         CALL REMX(nfte,NELM,IJ2,EMX,NBUF)
         lread=lread+nelm
         IF(NELM.EQ.0)GO TO 312
         L=0
      END IF
      L=L+1
      IF(max(ij2(1,l),ij2(2,l)).GT.ndim .OR. min(ij2(1,l),ij2(2,l))
     &   .LT.1)THEN
         WRITE(NFTA,*)' HELP: attempt to create element', ij2(1,l), 
     &                ij2(2,l), loff, l
         STOP
      END IF
      IF(ij2(1,l).EQ.ij2(2,l))diag(ij2(1,l))=EMX(L)
      IF(incore.EQ.1)THEN
         loff=loff+1
         offdg(loff)=EMX(L)
         lj1(loff)=ij2(1,l)
         lj2(loff)=ij2(2,l)
      END IF
      GO TO 311
C
 312  CLOSE(nfte)
      WRITE(NFTA,100)lread, NFTE
 100  FORMAT(/,I10,' Hamiltonian elements read from unit NFTE =',i3)
C
      RETURN
      END SUBROUTINE LDHAM

!*==mkdvm.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE MKDVM(hmat,eig,dgem,nhdim,nume,nfta,nfte,lembf,nelm,
     &                 CRITE,CRITC,CRITR,ORTHO,maxiter,IERR)
C
C***********************************************************************
C
C     MKDVM sets up matrix elements for Davidson diagonalisation
C     and calls diagonaliser DVDSON
C     Parameter INCORE determines if matrix elements can be held in core
c     (=1) or not (=0)
C***********************************************************************
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: CRITC, CRITE, CRITR, ORTHO
      INTEGER :: LEMBF, MAXITER, NELM, NFTA, NFTE, NHDIM, NUME, IERR
      REAL(KIND=wp), DIMENSION(nhdim) :: DGEM
      REAL(KIND=wp), DIMENSION(nume) :: EIG
      REAL(KIND=wp), DIMENSION(nhdim*nume) :: HMAT
      INTENT (IN) NUME
      INTENT (OUT) EIG, HMAT, IERR
      INTENT (INOUT) MAXITER
C
C Local variables
C
      LOGICAL :: HIEND
      INTEGER :: I, IERROR, IHIGH, ILOW, INCORE, IRWSZ, ITEMP,
     &           LIM, MBLOCK, NDIM, NLOOPS, NMV
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISELE, LJ
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: OFEM, WORK
      EXTERNAL OPCORE, OPDISK
C
C*** End of declarations rewritten by SPAG
C
      IERR=0
      NLOOPS=0
      NMV=0
c
c *** Check that there is enough space to hold all elements in core
c     calculate workspace for DVDSON routine
      LIM=MIN(NUME+20,NHDIM)
      IRWSZ=LIM*(2*NHDIM+LIM+9)+LIM*(LIM+1)/2+nume
      ALLOCATE(work(irwsz),isele(lim))
      ALLOCATE(ofem(nelm),lj(2*nelm),stat=ierror)
      IF(ierror.EQ.0)THEN
         incore=1
         WRITE(NFTA,2481)
 2481    FORMAT(' Non-zero matrix elements held in core')
         itemp=nelm
      ELSE
         incore=0
         WRITE(NFTA,2482)
 2482    FORMAT(' Non-zero matrix elements to be read from disk',
     &          ' at each iteration')
         itemp=lembf
         ALLOCATE(ofem(lembf),lj(2*lembf))
      END IF
c
c     Load diagonal elements, and others if enough room
      CALL LDHAM(dgem,ofem,lj,lj(itemp+1),nhdim,nfte,lembf,incore,nfta)
      ilow=1
      ihigh=nume
      IF(maxiter.LE.0)MAXITER=MAX(NUME*40,500)
c      write(nfta,2843)  CRITE,CRITC,CRITR,ORTHO,maxiter
c 2843 format(' CRITE =',D8.2,' CRITC =',D8.2,' CRITR =',D8.2,
c     1    ' ORTHO =',D8.2,' MAXITER =',I5)
      IF(incore.EQ.1)THEN
         mblock=1
         CALL DVDSON(opcore,nhdim,LIM,dgem,ILOW,IHIGH,ISELE,0,MBLOCK,
     &               CRITE,CRITC,CRITR,ORTHO,MAXITER,work,IRWSZ,HIEND,
     &               NLOOPS,NMV,IERR,ofem,lj,nelm)
      ELSE
         mblock=nume
         CALL DVDSON(opdisk,nhdim,LIM,dgem,ILOW,IHIGH,ISELE,0,MBLOCK,
     &               CRITE,CRITC,CRITR,ORTHO,MAXITER,work,IRWSZ,HIEND,
     &               NLOOPS,NMV,IERR,ofem,lj,lembf)
      END IF
c
      ndim=nume*nhdim
      DO i=1, ndim
         hmat(i)=work(i)
      END DO
      DO i=1, nume
         eig(i)=work(i+ndim)
      END DO
c
      WRITE(nfta,2800)nloops, nmv, IERR
 2800 FORMAT(/' Davidson diagonalisation completed:'/I4,' iterations',
     &       i10,' matrix vector multiplies',' IERR =',i4)
c
      DEALLOCATE(work,isele,lj,ofem)
      RETURN
      END SUBROUTINE MKDVM
!*==mvcivc.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE MVCIVC(NOCSF,NVEC,eigs,vecs,EIG,vec)
c
c     transfer the target vectors required for later
c
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NOCSF, NVEC
      REAL(KIND=wp), DIMENSION(nocsf) :: EIG
      REAL(KIND=wp), DIMENSION(nvec) :: EIGS
      REAL(KIND=wp), DIMENSION(nocsf,nocsf) :: VEC
      REAL(KIND=wp), DIMENSION(nocsf,nvec) :: VECS
      INTENT (IN) EIG, NOCSF, NVEC, VEC
      INTENT (OUT) EIGS, VECS
C
C Local variables
C
      INTEGER :: I, J
C
C*** End of declarations rewritten by SPAG
C
      DO i=1, nvec
         eigs(i)=eig(i)
         DO j=1, nocsf
            vecs(j,i)=vec(j,i)
         END DO
      END DO
      RETURN
      END SUBROUTINE MVCIVC
!*==mvcore.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE MVCORE(N,X,Y,A,LJ1,LJ2,NELM)
 
*=======================================================================
*
*       Computes the product of matrix A with vector x
*                               y = A x
*       where A(NxN) is a Symmetric Sparse matrix. Only the nonzero
*       elements of either the upper or lower triangular matrix are
*       kept and managed by Indices.
*
*=======================================================================
*
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: N, NELM
      REAL(KIND=wp), DIMENSION(NELM) :: A
      INTEGER, DIMENSION(NELM) :: LJ1, LJ2
      REAL(KIND=wp), DIMENSION(N) :: X, Y
      INTENT (IN) A, LJ1, LJ2, N, NELM, X
      INTENT (INOUT) Y
C
C Local variables
C
      INTEGER :: IBRA, IELM, IKET
C
C*** End of declarations rewritten by SPAG
C
C
************************************************************************
*
*   on entry
*   --------
*   N           the order of the matrix A
*   X           the vector to multiply A with
*   A           Linear array keeping the the nonzero elements of
*               the matrix A.
*   LJ1         Indices of Bras of A
*   LJ2         Indices of Kets of A
*
*   on exit
*   -------
*
*   Y           the result of the multiplication (dim=N)
*
************************************************************************
C
      y=XZERO
      DO ielm=1, nelm
         ibra=lj1(ielm)
         iket=lj2(ielm)
         y(ibra)=y(ibra)+a(ielm)*x(iket)
         IF(ibra.NE.iket)y(iket)=y(iket)+a(ielm)*x(ibra)
      END DO
C
      RETURN
      END SUBROUTINE MVCORE
!*==mvdisk.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE MVDISK(N,x,y,emx,ij2,nfte,nbuf)
*=======================================================================
*
*       Computes the product of matrix A with vectors x
*                               y = A x
*       where A(NxN) is a Symmetric Sparse matrix. Only the nonzero
*       elements of either the upper or lower triangular matrix are
*       kept and managed by Indices.  A is read from disk.
*
*=======================================================================
*
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      USE SCATCI_ROUTINES, ONLY: REMX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: N, NBUF, NFTE
      REAL(KIND=wp), DIMENSION(nbuf) :: EMX
      INTEGER, DIMENSION(2,nbuf) :: IJ2
      REAL(KIND=wp), DIMENSION(N) :: X, Y
      INTENT (IN) N, X
      INTENT (INOUT) Y
C
C Local variables
C
      INTEGER :: IBRA, IKET, L, NELM
C
C*** End of declarations rewritten by SPAG
C
C
************************************************************************
*
*   on entry
*   --------
*   N           the order of the matrix A
*   x           the vector to multiply A with
*   A           Linear array keeping the the nonzero elements of
*               the matrix A.
*   LJ1         Indices of Bras of A
*   LJ2         Indices of Kets of A
*
*   on exit
*   -------
*
*   y           the result of the multiplication (dim=N)
*
************************************************************************
C
 
      y=XZERO
 
      REWIND nfte
      READ(nfte)
      L=0
      nelm=0
 311  IF(L.EQ.NELM)THEN
         CALL REMX(nfte,NELM,IJ2,EMX,NBUF)
         IF(NELM.EQ.0)RETURN
         L=0
      END IF
      L=L+1
C
      ibra=ij2(1,l)
      iket=ij2(2,l)
      y(ibra)=y(ibra)+emx(l)*x(iket)
      IF(ibra.NE.iket)y(iket)=y(iket)+emx(l)*x(ibra)
 
      GO TO 311
C
      END SUBROUTINE MVDISK
!*==prci.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE PRCI(NOCSF,NSTATE,EE,ET,XR,CR,NCIPFG,NKEY,KEYCSF,LARGE,
     &                THRPRT,NFT)
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : ZERO=>XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: ET, THRPRT, XR
      INTEGER :: LARGE, NCIPFG, NFT, NKEY, NOCSF, NSTATE
      REAL(KIND=wp), DIMENSION(nocsf,*) :: CR
      REAL(KIND=wp), DIMENSION(*) :: EE
      INTEGER, DIMENSION(2) :: KEYCSF
      INTENT (IN) KEYCSF, NCIPFG, NKEY, NSTATE, THRPRT
      INTENT (INOUT) LARGE
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(10) :: CBIG
      REAL(KIND=wp), DIMENSION(20) :: CKEY, CKEY2, CSUMC
      REAL(KIND=wp) :: CSUMR
      INTEGER :: I, I1, INC, J, KC, NKEYN, NKF, NS, NSF, NSI
      INTEGER, DIMENSION(10) :: IBIG
      character(len=4) :: IBL='    '
      INTEGER, DIMENSION(20) :: KCSFI
C
C*** End of declarations rewritten by SPAG
C
      NKF=MIN(NKEY,20)
      LARGE=MIN(LARGE,10)
      IF(NKF.LE.0)GO TO 35
C
      I1=0
      DO I=1, NKF
         KC=KEYCSF(I) ! JMC are the array bounds really OK here for keycsf???
         DO J=1, nocsf
            IF(KC.EQ.J)GO TO 20
         END DO
         CYCLE
 20      I1=I1+1
         KCSFI(I1)=KC
      END DO
      NKF=I1
C
 35   IF(NCIPFG.NE.1)GO TO 200
      IF(NKF.LE.0)GO TO 120
C
      DO I=1, NKF
         CSUMC(I)=ZERO
      END DO
C
      NSF=0
 50   WRITE(NFT,60)NOCSF, (KCSFI(I),I=1,NKF)
 60   FORMAT('1EIGENVECTOR RECAPITULATION  NOCSF =',I6,/,9X,'ENERGY',3X,
     &       10I9/18X,10I9)
      WRITE(NFT,65)
 65   FORMAT(' ')
C
      NSI=NSF+1
      INC=20
      IF(NKF.GE.10)INC=10
      NSF=MIN(NSF+INC,NSTATE)
C
      DO I=NSI, NSF
         CSUMR=ZERO
         DO J=1, NKF
            CKEY(J)=CR(KCSFI(J),I)
            CKEY2(J)=CKEY(J)*CKEY(J)
            CSUMC(J)=CSUMC(J)+CKEY2(J)
            CSUMR=CSUMR+CKEY2(J)
         END DO
         WRITE(NFT,85)I, EE(I)+ET, (CKEY(J),J=1,NKF)
 85      FORMAT(I4,F14.6,2X,10F9.4/22X,10F9.4)
         WRITE(NFT,105)(IBL,CKEY2(J),J=1,NKF), IBL, CSUMR
 105     FORMAT(21X,10(A1,'(',F6.3,')')/21X,11(A1,'(',F6.3,')'))
         WRITE(NFT,65)
      END DO
      IF(NSF.LT.NSTATE)GO TO 50
      WRITE(NFT,115)(CSUMC(J),J=1,NKF)
 115  FORMAT(13X,'SUMS',3X,10F9.4/28X,10F9.4)
C
C     PRINT LARGE COEFFICIENTS
C
 120  IF(LARGE.LE.0)GO TO 200
      WRITE(NFT,121)
 121  FORMAT(//)
      DO NS=1, NSTATE
         DO I=1, LARGE
            CBIG(I)=ZERO
            IBIG(I)=0
         END DO
         DO I=1, LARGE
            DO J=1, nocsf
C
               IF(ABS(CR(J,NS)).LE.ABS(CBIG(I)))CYCLE
               IF(I.GT.1)THEN
                  IF(ABS(CR(J,NS)).EQ.ABS(CBIG(I-1)))THEN
                     IF(J.LE.IBIG(I-1))CYCLE
                  ELSE IF(ABS(CR(J,NS)).GT.ABS(CBIG(I-1)))THEN
                     CYCLE
                  END IF
               END IF
               IBIG(I)=j
               CBIG(I)=CR(J,NS)
            END DO
         END DO
C
         IF(NKF.GT.0)THEN
            DO I=1, NKF
               CKEY(I)=CR(KCSFI(I),NS)
            END DO
         END IF
C
         WRITE(NFT,165)NS, EE(NS)+ET, (IBIG(I),CBIG(I),I=1,LARGE)
 165     FORMAT(I4,F14.6,3X,10(I4,F7.3))
         IF(NKF.GT.0)WRITE(NFT,170)(KCSFI(J),CKEY(J),J=1,NKF)
 170     FORMAT(21X,10(I4,F7.3))
         WRITE(NFT,65)
      END DO
 200  IF(NCIPFG.EQ.2 .AND. NKF.GT.0)THEN
         NKEYN=0
         DO NS=1, NSTATE
            CSUMR=ZERO
            DO I=1, NKF
               CSUMR=CSUMR+CR(KCSFI(I),NS)*CR(KCSFI(I),NS)
            END DO
            IF(CSUMR.LT.THRPRT)CYCLE
            NKEYN=NKEYN+1
C
            CALL PRCIVC(NOCSF,NS,EE(NS),ET,XR,CR(1,NS),NFT)
C
         END DO
      END IF
C
      IF(NCIPFG.EQ.3)THEN
         DO NS=1, NSTATE
            CALL PRCIVC(NOCSF,NS,EE(NS),ET,XR,CR(1,NS),NFT)
         END DO
      END IF
C
      RETURN
      END SUBROUTINE PRCI
!*==prcivc.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE PRCIVC(NOCSF,NSTAT,EE,EN,XR,CIV,NFT)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: EE, EN, XR
      INTEGER :: NFT, NOCSF, NSTAT
      REAL(KIND=wp), DIMENSION(50,*) :: CIV
      INTENT (IN) CIV, EE, EN, NFT, NOCSF, NSTAT, XR
C
C Local variables
C
      REAL(KIND=wp) :: ED, ET
      INTEGER :: I, IA, IMAX, J, JA, JB, JMAX
C
C*** End of declarations rewritten by SPAG
C
      ET=EE+EN
      ED=ET-XR
      WRITE(NFT,500)NOCSF, NSTAT
      WRITE(NFT,510)EE, ET, XR, ED
C
      JMAX=(nocsf-1)/50+1
      IMAX=nocsf-50*(JMAX-1)
      JA=1
C
 200  JB=JA+7
      IF(JB.GE.JMAX)GO TO 300
      WRITE(NFT,520)((50*(j-1)+i,CIV(I,J),J=JA,JB),I=1,50)
      JA=JB+1
      WRITE(NFT,530)
      GO TO 200
C
 300  IA=1
 310  DO I=IA, IMAX
         WRITE(NFT,520)(50*(j-1)+I,CIV(I,J),J=JA,JMAX)
      END DO
      IF(IMAX.EQ.50 .OR. JA.EQ.JMAX)RETURN
      JMAX=JMAX-1
      IA=IMAX+1
      IMAX=50
      GO TO 310
C
      RETURN
C
 500  FORMAT('1 NOCSF=',I6,8X,' STATE',I5,9X,//,5X,'E (ELECTRONIC)',11X,
     &       'E (TOTAL)',7X,'E (REFERENCE)',6X,'E (DIFFERENCE)')
 510  FORMAT(/F19.8,3F20.8,/)
 520  FORMAT((1X,8(I5,F11.6)))
 530  FORMAT('1'/'0'/'0'/'0')
      END SUBROUTINE PRCIVC
!*==prtham.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE PRTHAM(HAM,NHDIM,NFT,npflg)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NFT, NHDIM, NPFLG
      REAL(KIND=wp), DIMENSION(nhdim,nhdim) :: HAM
      INTENT (IN) HAM, NFT, NHDIM, NPFLG
C
C Local variables
C
      INTEGER :: I, J
C
C*** End of declarations rewritten by SPAG
C
C
C***********************************************************************
C
C     PRHAM PRINTS OUT THE ELEMENTS OF THE ENERGY MATRIX
C
C           NPFLG  LE  0  DIAGONAL ELEMENTS ONLY PRINTED
C                  GT  0  ALL NON-ZERO ELEMENTS PRINTED
C
C***********************************************************************
C
C
C
      IF(NPFLG.GT.0)THEN
C
         WRITE(NFT,210)
         DO i=1, nhdim
            WRITE(NFT,220)(i,j,ham(i,j),j=1,i)
         END DO
      ELSE IF(npflg.LT.0)THEN
C
         WRITE(NFT,230)
         WRITE(NFT,220)(i,i,ham(i,i),i=1,nhdim)
      END IF
      RETURN
C
 210  FORMAT('  ELEMENTS OF ENERGY MATRIX'//)
 220  FORMAT(3(2I5,D16.8))
 230  FORMAT('  DIAGONAL ELEMENTS OF ENERGY MATRIX'//)
      END SUBROUTINE PRTHAM
!*==qldiag.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE QLDIAG(MK,VEC,EIG)
      USE precisn, ONLY : wp ! for specifying the kind of reals
      USE blas_lapack_gbl, ONLY : blasint
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: MK
      REAL(KIND=wp), DIMENSION(MK) :: EIG
      REAL(KIND=wp), DIMENSION(MK,MK) :: VEC
C
C Local variables
C
      INTEGER :: ERR
      INTEGER(KIND=blasint) :: N, LIWORK, IFAIL, LX
      INTEGER(KIND=blasint), ALLOCATABLE, DIMENSION(:) :: IWORK
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: X
C
C*** End of declarations rewritten by SPAG
C
C
C*********************************************************************
C
C     QLDIAG DIAGONALIZES THE FULL HAMILTONIAN MATRIX
C
C*********************************************************************
C
C     HOUSEHOLDER - QL DIAGONALIZATION  ( MUST BE IN-CORE )
C

      ALLOCATE(IWORK(3+5*MK),X(1+6*MK+2*MK**2),STAT=ERR)
      IF (ERR /= 0) THEN
         PRINT *, 'allocation error with workspace for dsyevd'
         STOP
      END IF
C
      N=MK
      IFAIL=0
      LX = 1 + 6*MK + 2*MK**2
      LIWORK = 3 + 5*MK
C
      CALL DSYEVD('V','L',N,VEC,N,EIG,X,LX,IWORK,LIWORK,IFAIL)
C
      DEALLOCATE(IWORK,X,STAT=ERR)
      IF(ERR /= 0)THEN
         PRINT *, 'deallocation error with workspace for dsyevd'
         STOP
      END IF
C
      RETURN
      END SUBROUTINE QLDIAG
!*==writcid.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE WRITCID(nftw,nciset,NAME,NHE,NHD,DTNUC,NOCSF,NSTAT,EIG,
     &                   VEC,iphz,dgem,NFTA,npflg,npcvc)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE scatci_data, ONLY : MEIG
      USE SCATCI_ROUTINES, ONLY: movew, prthd, civio
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      character(len=120) :: NAME
      INTEGER :: NCISET, NFTA, NFTW, NOCSF, NPCVC, NPFLG, NSTAT
      REAL(KIND=wp), DIMENSION(nocsf) :: DGEM
      REAL(KIND=wp), DIMENSION(41) :: DTNUC
      REAL(KIND=wp), DIMENSION(*) :: EIG, VEC
      INTEGER, DIMENSION(*) :: IPHZ
      INTEGER, DIMENSION(10) :: NHD
      INTEGER, DIMENSION(20) :: NHE
      INTENT (IN) NCISET, NPFLG
      INTENT (INOUT) NHD
C
C Local variables
C
      INTEGER :: NALM, NEIG, NTH
C
C*** End of declarations rewritten by SPAG
C
C
C**********************************************************************
C
C     WRITCID writes CI data to unit NFTW in format used for diatomic
C     targets
C
C**********************************************************************
C
c
      NTH=NCISET
      WRITE(NFTA,3822)NTH
 3822 FORMAT(' VECTOR MATRIX STORED AS SET NTH =',I6,' LISTED LAST',/)
C
      CALL MOVEW(NFTW,NTH,NALM,NPCVC,NFTA)
C
      IF(NALM.NE.0)THEN
         WRITE(NFTA,2940)
 2940    FORMAT(' ERROR POSITIONING FOR OUTPUT OF CI COEFFICIENTS',//)
         STOP
      END IF
C
      IF(nstat.LT.nocsf)nhd(10)=1
      WRITE(NFTW,ERR=2950)NTH, NHD, NAME, NHE, DTNUC
C
      neig=nstat
      IF(npflg.LE.0 .AND. nstat.GT.meig)neig=meig
C
      CALL PRTHD(NTH,NHD,NAME,NHE,DTNUC,NEIG,EIG,NFTA)
      CALL CIVIO(NFTW,0,NOCSF,NSTAT,EIG,VEC,NALM,iphz,dgem)
      IF(NALM.NE.0)GO TO 2950
 
      REWIND NFTW
c
      RETURN
c
 2950 WRITE(NFTA,2960)
 2960 FORMAT(' ERROR WRITING CI OUTPUT FILE',//)
      STOP
c
      END SUBROUTINE WRITCID
!*==writcip.spg  processed by SPAG 6.56Rc at 10:10 on  8 Nov 2010
      SUBROUTINE WRITCIP(nftw,nset,NAME,e0,nocsf,nstat,mgvn,s,sz,nelt,
     &                   EIG,VEC,iphase,dgem,ntgsym,mcont,notgs,nfti,
     &                   NFT,npflg,npcvc,ukrmolp_ints)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE params, ONLY : ccidata, c8stars, cblank, cpoly
      USE scatci_data, ONLY : MEIG
      USE ukrmol_interface_gbl, ONLY : GET_GEOM
      USE SCATCI_ROUTINES, ONLY: movep, search, civio
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: E0, S, SZ
      INTEGER :: MGVN, NELT, NFT, NFTI, NFTW, NOCSF, NPCVC, NPFLG, NSET, 
     &           NSTAT, NTGSYM
      character(len=120) :: NAME
      REAL(KIND=wp), DIMENSION(*) :: DGEM, EIG, VEC
      INTEGER, DIMENSION(*) :: IPHASE
      INTEGER, DIMENSION(ntgsym) :: MCONT, NOTGS
      LOGICAL :: ukrmolp_ints
      INTENT (IN) E0, MCONT, MGVN, NAME, NELT, NOTGS, NPFLG, NTGSYM, S, 
     &            SZ, ukrmolp_ints
      INTENT (INOUT) NSET
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: CHARGE, XNUC, YNUC, 
     &                                            ZNUC
      character(len=8), ALLOCATABLE, DIMENSION(:) :: CNAME
      INTEGER :: I, IFAIL, II, NALM, NEIG, NNUC, NREC, NTH
C
C*** End of declarations rewritten by SPAG
C
C
C**********************************************************************
C
C     WRITCIP writes CI data to unit NFTW in format used for
C     polyatomic targets.  Data is always appended to the end of
C     existing files
C
C**********************************************************************
c
c.... Position file for next set
      NTH=NSET
C
      CALL MOVEP(NFTW,NTH,NALM,NPCVC,NFT)
      IF(NALM.NE.0)THEN
         WRITE(NFT,2940)
 2940    FORMAT(' ERROR POSITIONING FOR OUTPUT OF CI COEFFICIENTS',//)
         STOP
      END IF
C
C
C.... Read integral header for geometry and charge data
      if (ukrmolp_ints) then !Take the UKRmol+ data from the orbital_basis_data object in the interface module
         call GET_GEOM(nnuc,cname,xnuc,ynuc,znuc,charge)
      else !Read-in the data from the SWEDEN transformed integrals file
         REWIND nfti
         CALL SEARCH(nfti,cpoly,ifail)
         READ(nfti)nnuc
         ALLOCATE(xnuc(nnuc),ynuc(nnuc),znuc(nnuc),charge(nnuc),
     &            CNAME(nnuc))
         DO i=1, nnuc
            READ(nfti)cname(i), ii, xnuc(i), ynuc(i), znuc(i), charge(i)
         END DO
      endif
c
c.... Write header
      WRITE(nftw)c8stars, cblank, cblank, ccidata
      nset=nth
      nrec=nnuc+nstat+1
      WRITE(NFT,3822)Nset
 3822 FORMAT(' CI data stored as set number',I3)
      WRITE(NFTw)Nset, nrec, NAME, nnuc, nocsf, nstat, mgvn, s, sz, 
     &           nelt, e0, ntgsym, (mcont(i),notgs(i),i=1,ntgsym)
      DO i=1, nnuc
         WRITE(nftw)cname(i), xnuc(i), ynuc(i), znuc(i), charge(i)
      END DO
c
C     Write CI eigenvalues and eigenvectors
C
      CALL CIVIO(NFTW,0,NOCSF,nstat,EIG,VEC,NALM,iphase,dgem)
      IF(NALM.NE.0)GO TO 2900
C
c...  Print summary of output data
      WRITE(NFT,100)NSET, name
      WRITE(nft,103)mgvn, s, sz, nelt, nnuc
      WRITE(NFT,120)(cname(i),xnuc(i),ynuc(i),znuc(i),charge(i),i=1,
     &              nnuc)
      WRITE(NFT,101)nocsf, nstat, e0
      IF(ntgsym.GT.0)WRITE(nft,104)ntgsym, 
     &                             (i,mcont(i),notgs(i),i=1,ntgsym)
      neig=nstat
      IF(npflg.LE.0 .AND. nstat.GT.meig)neig=meig
      WRITE(NFT,102)(EIG(I)+E0,I=1,NEIG)
c
      DEALLOCATE(xnuc,ynuc,znuc,charge,cname)
      RETURN
c
 2900 WRITE(NFT,2960)
 2960 FORMAT(' ERROR WRITING CI OUTPUT FILE',//)
      STOP
c
 100  FORMAT(/' SET',I4,4X,A)
 101  FORMAT(/' NOCSF=',I5,4X,'NSTAT=',I5,4X,'EN   =',F20.10)
 103  FORMAT(/' MGVN =',I2,4x,'s =',f6.1,4x,'sz =',f6.1,4x,'NELT =',I5,
     &       4x,'NNUC =',I3)
 104  FORMAT(/' NTGSYM=',i4/'  I   MCONT  NOTGT'/(i3,2I7))
 102  FORMAT(/' EIGEN-ENERGIES',/(16X,5F20.10))
 120  FORMAT(/' Nuclear data     X         Y         Z       Charge'/
     &       (3x,a8,2x,4F10.6))
c
      END SUBROUTINE WRITCIP
!*==opcore.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE OPCORE(N,M,B,C,A,LJ1,LJ2,NELM)
C=======================================================================
C
C       Computes the product of matrix A with a block of vectors B(N,M)
C                               C=A B
C       where A(NxN) is a Symmetric Sparse matrix. Only the nonzero
C       elements of either the upper or lower triangular matrix are
C       kept and managed by Indices.
C
C=======================================================================
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      USE maths, ONLY : MATHS_CVECA_PLUS_VECB
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: M, N, NELM
      REAL(KIND=wp), DIMENSION(NELM) :: A
      REAL(KIND=wp), DIMENSION(N,M) :: B, C
      INTEGER, DIMENSION(NELM) :: LJ1, LJ2
      INTENT (IN) A, LJ1, LJ2, NELM
C
C Local variables
C
      REAL(KIND=wp) :: ELEM
      INTEGER :: IBRA, IELM, IKET
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C   on entry
C   --------
C   N           the order of the matrix A
C   B           the matrix (block of vectors) to multiply A with
C   A           Linear array keeping the the nonzero elements of
C               the matrix A.
C   LJ1         Indices of Bras of A
C   LJ2         Indices of Kets of A
C
C   on exit
C   -------
C
C   C           the result of the multiplication (dim=NxM)
C
C***********************************************************************
C
      c=XZERO
 
      DO ielm=1, nelm
         ibra=lj1(ielm)
         iket=lj2(ielm)
         elem=a(ielm)
         CALL MATHS_CVECA_PLUS_VECB(m,elem,b(iket,1),n,c(ibra,1),n)
         IF(ibra.NE.iket)
     &      CALL MATHS_CVECA_PLUS_VECB(m,elem,b(ibra,1),n,c(iket,1),n)
      END DO
C
      RETURN
      END SUBROUTINE OPCORE
!*==opdisk.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      SUBROUTINE OPDISK(N,M,B,C,emx,ij2,nfte,nbuf)
C=======================================================================
C
C       Computes the product of matrix A with a block of vectors B(N,M)
C                               C=A B
C       where A(NxN) is a Symmetric Sparse matrix. Only the nonzero
C       elements of either the upper or lower triangular matrix are
C       kept and managed by Indices.  A is read from disk.
C
C=======================================================================
C
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      USE consts, ONLY : XZERO
      USE maths, ONLY : MATHS_CVECA_PLUS_VECB
      USE SCATCI_ROUTINES, ONLY: REMX
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: M, N, NBUF, NFTE
      REAL(KIND=wp), DIMENSION(N,M) :: B, C
      REAL(KIND=wp), DIMENSION(nbuf) :: EMX
      INTEGER, DIMENSION(2,nbuf) :: IJ2
C
C Local variables
C
      REAL(KIND=wp) :: ELEM
      INTEGER :: IBRA, IKET, L, NELM
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C   on entry
C   --------
C   N           the order of the matrix A
C   B           the matrix (block of vectors) to multiply A with
C   A           Linear array keeping the the nonzero elements of
C               the matrix A.
C   LJ1         Indices of Bras of A
C   LJ2         Indices of Kets of A
C
C   on exit
C   -------
C
C   C           the result of the multiplication (dim=NxM)
C
C***********************************************************************
C
 
      c=XZERO
 
      REWIND nfte
      READ(nfte)
      L=0
      nelm=0
 311  IF(L.EQ.NELM)THEN
         CALL REMX(nfte,NELM,IJ2,EMX,NBUF)
         IF(NELM.EQ.0)RETURN
         L=0
      END IF
      L=L+1
C
      ibra=ij2(1,l)
      iket=ij2(2,l)
      elem=emx(l)
      CALL MATHS_CVECA_PLUS_VECB(m,elem,b(iket,1),n,c(ibra,1),n)
      IF(ibra.NE.iket)
     &   CALL MATHS_CVECA_PLUS_VECB(m,elem,b(ibra,1),n,c(iket,1),n)
 
      GO TO 311
C
      END SUBROUTINE OPDISK
