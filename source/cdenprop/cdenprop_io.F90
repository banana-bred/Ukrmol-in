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
      module cdenprop_io

      contains

      !Subroutines to read in the CSF file
      !***********************************
      !***********************************


      !subroutines taken from SCATCI
      !*****************************
      SUBROUTINE RDNFTO(megul,nob,nob0,nobl,nob0l,nsym,ndtrf,nelt,nodo,nocsf,iphz,nctarg,nctgt,notgt,mcont,gucont,ntgsym,ntgcon)
      !
      !
      !     read CSFs as generated and projected by CONGEN (part 1)
      !            NODO(NOCSF) NO OF DTRS IN EACH INPUT WF
      !            NDTRF(NELT) REFERENCE DTR IN EXPANDED FORM
      !IPOSIT is a flag for positron scattering
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NOB(NSYM),NDTRF(NELT),NODO(NOCSF),iphz(nctarg)
      DIMENSION nobl(2*nsym),nob0(nsym),nob0l(2*nsym)
      dimension nctgt(ntgsym),notgt(ntgsym),mcont(ntgsym)
      integer gucont(ntgsym)
      !
      if (nctarg.gt.0) then
         if (ntgcon .gt. 0) then
            read(megul) iphz,nctgt,notgt,mcont,gucont

         end if

         if (ntgcon .le. 0) then
            read(megul) iphz

         end if

      endif

      read(megul) NOB,NDTRF,NODO,IPOSIT,NOB0,nobl,nob0l

      RETURN
      END SUBROUTINE RDNFTO

      SUBROUTINE RDWF(M,ICDO,INDO,NDO,LNDOF,CDO,LCDOF,megul)
      !
      !     read CSFs as generated and projected by CONGEN (part 2)
      !            NDO(LNDOF)  DTRS IN PACKED FORM
      !            CDO(LCDOF)  CORRESPONDING COEFFICIENTS
      !            INDO(NOCSF+1) INDEX ON DTRS
      !            ICDO(NOCSF+1) INDEX ON COEF.
      !
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ICDO(M),INDO(M),NDO(LNDOF),CDO(LCDOF)
      !
      read(megul) ICDO,INDO
      read(megul) NDO
      read(megul) CDO

      RETURN

      END subroutine RDWF

      SUBROUTINE CWBOPN(LUNIT)
      !***********************************************************************
      !
      !
      !     OPENS UP A SEQUENTIAL UNIT
      !
      !***********************************************************************
      !
            LOGICAL OP
      !
            INQUIRE(UNIT=LUNIT,OPENED=OP)
      !
            IF (.NOT.OP) OPEN(UNIT=LUNIT,STATUS='UNKNOWN',FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
      !
            REWIND LUNIT
      !
            RETURN
      END SUBROUTINE CWBOPN

      SUBROUTINE WRNFTO(SNAME,MGVN,S,SZ,R,PIN,NORB,NSRB, &
          &                  NOCSF,NELT,IDIAG,NSYM,SYMTYP, &
          &                  NOB,NDTRF,NODO,M,ICDO,INDO, &
          &                  NDO,lndi,CDO,lcdi,NFTO,nobl,nx,&
          &                  npflg,thres,iposit,nob0,nob0l,nctarg,&
          &                  ntgsym,notgt,nctgt,mcont,gucont,iphz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER SYMTYP,gucont
      CHARACTER(LEN=120) :: NAME
      character(len=80) :: sname
      DIMENSION NOB(NSYM),NDTRF(NELT),NODO(NOCSF),ICDO(M),INDO(M),  &
      &          NDO(LNDI),CDO(LCDI),notgt(ntgsym),nctgt(ntgsym), &
      &          mcont(ntgsym),gucont(ntgsym),iphz(nctarg)
      DIMENSION NOBW(20),nobl(nx),nob0(nsym),nob0l(nx),npflg(6)
      name=sname
      NORBW = NORB
      NSRBW = NSRB
      DO I = 1, NSYM
        NOBW(I) = NOB(I)
      END DO
      IF (IPOSIT .LT. 0) THEN
        NORBW = 0
        DO I = 1, NSYM
          NOBW(I) = NOB(I) - NOB0(I)
          NORBW = NORBW + NOBW(I)
        END DO
        NSRBW = 2 * (2 * NORBW - NOBW(1))
      ENDIF
      REWIND NFTO
      WRITE(NFTO) NAME,MGVN,S,SZ,R,PIN,NORBW,NSRBW, &
      &            NOCSF,NELT,lcdi,IDIAG,NSYM,SYMTYP,lndi, &
      &            npflg,thres,nctarg,ntgsym
      if (ntgsym.gt.0) write(nfto) iphz,nctgt,notgt,mcont,gucont
      if (ntgsym.le.0) write(nfto) iphz
      WRITE(NFTO) (NOBW(I),I=1,NSYM),NDTRF,NODO, &
      &             IPOSIT,NOB0,nobl,nob0l
      WRITE(NFTO) ICDO,INDO
      WRITE(NFTO) NDO
      WRITE(NFTO) CDO
      RETURN
      END subroutine WRNFTO

      SUBROUTINE WRWF (NFT,N1,NODO,N2,CDO,N3,NDO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NODO(N1),NDO(N3),CDO(N2)

      REWIND NFT

      WRITE(NFT) N1,NODO
      WRITE(NFT) N2,CDO
      WRITE(NFT) N3,NDO

      REWIND NFT
      RETURN
      END subroutine WRWF

      SUBROUTINE PTPWF(NFTW,NOCSF,NELT,NDTRF,NODI,INDI,ICDI,NDI,CDI)
      !Prints packed? CSFs
      !*******************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER(I-N)
      save
      DIMENSION NDTRF(*),NODI(*),NDI(*),CDI(*)
      DIMENSION INDI(*),ICDI(*)
!
      WRITE(NFTW,139) (NDTRF(I),I=1,NELT)
  139 FORMAT(' REFERENCE DETERMINANT'//(1X,20I5))
      WRITE(NFTW,137)
  137 FORMAT('  CSF',9X,'COEFFICIENT',2X,'NSO'/)
      DO N=1,NOCSF

        MA=NODI(N)

        MB=INDI(N)

        MC=ICDI(N)-1

        MD=NDI(MB)

        WRITE(NFTW,138) N,CDI(MC+1),MD,(NDI(MB+I),I=1,2*MD)
  138   FORMAT(1X,I4,D20.10,I5,2X,20I5/(32X,20I5))
        MB=MB+MD+MD+1
        DO K=2,MA
          MD=NDI(MB)
          WRITE(NFTW,140) CDI(MC+K),MD,(NDI(MB+I),I=1,2*MD)
  140     FORMAT(5X,D20.10,I5,2X,20I5/(32X,20I5))
          MB=MB+MD+MD+1
        END DO
      END DO
      RETURN
      END SUBROUTINE PTPWF

      ! IO for the target CI vectors
      !*******************************************************************************
      !********************************************************************************

      SUBROUTINE CIVIO (NFT,NRW,NK,NS,NSRW,EI,CV,NALM,iphz,dg)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!***********************************************************************
!
!
!     CIVIO CONTROLS THE I/O OF CI COEFFICIENTS AND STATE DATA
!
!***********************************************************************
!
      DIMENSION EI(ns),CV(nk,nsrw),iphz(nk),dg(nk)
!
!     Read/Write ENERGY AND CSF SPECIFICATION FOR EACH CI STATE
!     annd then Read/Write COEFFICIENTS FOR EACH CI STATE
      NALM=0
      IF (NRW.eq.0) then
      WRITE (NFT,ERR=200) iphz,EI,dg
        do i=1,nsrw
          WRITE (NFT,ERR=200) i,(CV(j,i),j=1,nk)
        end do
      else
        READ  (NFT,ERR=200) iphz,EI
        do i=1,nsrw
          READ  (NFT,ERR=200) M,(CV(j,i),j=1,nk)
        end do
      endif
      return
!
200  nalm=1
      RETURN
      END SUBROUTINE CIVIO

!       SUBROUTINE MOVEP (NFT,NTH,NALM,NPFLG,NFT1)
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! !
! !***********************************************************************
! !
! !     MOVEP LOCATES THE NTH-DATASET CONTAINED ON THE ALCHEMY CI
! !           DUMPFILE (for polyatomics) ASSOCIATED WITH UNIT NFT
! !
! !           NPFLG = 0     NO PRINTOUT DURING SEARCH
! !                 = 1     PRINT HEADER LABELS ENCOUNTERED
! !           NALM  = 0     NO ERRORS DETECTED
! !                 = 1     ERRORS DETECTED
! !
! !***********************************************************************
! !
!       CHARACTER*120 NAME
!       character*32 header
!       PARAMETER (MEIG=20)
!       logical opn
!       DIMENSION EIG(MEIG)
! !
!       NALM = 0
!       inquire(unit=nft,opened=opn)
!       if (opn) then
!         REWIND NFT
!       else
!         open(unit=nft,form='unformatted')
!       endif
!       IF (NTH.EQ.1) RETURN
!       M = NTH-1
!       IF (NTH.EQ.0) M=2000
!       DO 30 I=1,M
!       READ(NFT,END=50) header
!       READ(NFT,END=50) NT,nrec,NAME,nnuc,nocsf,nstat,mgvn,s,sz,nelt,e0
!       DO 10 K=1,NNUC
!       READ (NFT)
! 10    continue
! !
!       NEIG = MIN(MEIG,NSTAT)
!       READ(NFT) (IDUM,J=1,nocsf),(EIG(J),J=1,NEIG)
! !
!       IF (NPFLG.GT.0) then
!          WRITE (NFT1,104)  header
!          WRITE (NFT1,100)  NT,name
!          write (nft1,103)  mgvn,s,sz,nelt,nnuc
!          WRITE (NFT1,101)  nocsf,nstat,e0
!          WRITE (NFT1,102) (EIG(J)+E0,J=1,NEIG)
!
!       END IF
!   104 format (/' Header:',4x,a)
!   100 FORMAT (/' SET',I4,4X,A)
!   101 FORMAT (/' NOCSF=',I5,4X,'NSTAT=',I5,4X,'EN   =',F20.10)
!   103 format (/' MGVN =',I2,4x,'s =',f6.1,4x,'sz =',f6.1,4x, &
!     &             'NELT =',I5,4x,'NNUC =',I3)
!   102 FORMAT (/' EIGEN-ENERGIES',/(16X,5F20.10))
! !
!       DO 20 K=1,NSTAT
!       READ (NFT)
! 20    continue
! !
! 30    CONTINUE
!
!       RETURN
! !
! 50    IF (NTH.NE.0) THEN
!         NALM = 1
!
!       ELSE
!         NTH = NT+1
! !        BACKSPACE NFT
!
!       ENDIF
!
!       RETURN
!
!       END SUBROUTINE MOVEP
      SUBROUTINE MOVEP(NFT,NTH,NALM,NPFLG,NFT1)
      use cdenprop_defs
!       USE precisn, ONLY : wp ! for specifying the kind of reals
!       USE scatci_data, ONLY : MEIG
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: NALM, NFT, NFT1, NPFLG, NTH
      INTENT (IN) NFT, NFT1, NPFLG
      INTENT (OUT) NALM
      INTENT (INOUT) NTH
!
! Local variables
!
      integer, parameter :: meig=20
      REAL(KIND=idp) :: E0, S, SZ
      REAL(KIND=idp), DIMENSION(MEIG) :: EIG
      CHARACTER(LEN=32) :: HEADER
      INTEGER :: I, IDUM, J, K, M, MGVN, NEIG, NELT, NNUC, NOCSF, NREC, &
     &           NSTAT, NT
      CHARACTER(LEN=120) :: NAME
      LOGICAL :: OPN
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     MOVEP LOCATES THE NTH-DATASET CONTAINED ON THE ALCHEMY CI
!           DUMPFILE (for polyatomics) ASSOCIATED WITH UNIT NFT
!
!           NPFLG = 0     NO PRINTOUT DURING SEARCH
!                 = 1     PRINT HEADER LABELS ENCOUNTERED
!           NALM  = 0     NO ERRORS DETECTED
!                 = 1     ERRORS DETECTED
!
!***********************************************************************
!
!
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
         READ(NFT,END=50)NT, nrec, NAME, nnuc, nocsf, nstat, mgvn, s,&
     &                   sz, nelt, e0
         DO K=1, NNUC
            READ(NFT)
         END DO
!
         NEIG=MIN(MEIG,NSTAT)
         READ(NFT)(IDUM,J=1,nocsf), (EIG(J),J=1,NEIG)
!
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
 103     FORMAT(/' MGVN =',I2,4x,'s =',f6.1,4x,'sz =',f6.1,4x,'NELT =',&
     &          I5,4x,'NNUC =',I3)
 102     FORMAT(/' EIGEN-ENERGIES',/(16X,5F20.10))
!
         DO K=1, NSTAT
            READ(NFT)
         END DO
!
      END DO
      RETURN
!
 50   IF(NTH.NE.0)THEN
         NALM=1
      ELSE
         NTH=NT+1
!        BACKSPACE NFT
      END IF
      RETURN
      END SUBROUTINE MOVEP


      SUBROUTINE FINDP(NFT,NTH,NALM,NPFLG,NFT1,mgvn_in, spin_in,gflag,itarg)
      use cdenprop_defs
!       USE precisn, ONLY : wp ! for specifying the kind of reals
!       USE scatci_data, ONLY : MEIG
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: NALM, NFT, NFT1, NPFLG, NTH,mgvn_in,gflag,itarg
      real(kind=idp) :: spin_in
      INTENT (IN) NFT, NFT1, NPFLG,mgvn_in,spin_in,gflag,itarg
      INTENT (OUT) NALM
      INTENT (INOUT) NTH
!
! Local variables
!
      integer, parameter :: meig=20
      REAL(KIND=idp) :: E0, S, SZ
      REAL(KIND=idp), DIMENSION(MEIG) :: EIG
      CHARACTER(LEN=32) :: HEADER
      INTEGER :: I, IDUM, J, K, M, MGVN, NEIG, NELT, NNUC, NOCSF, NREC, &
     &           NSTAT, NT
      CHARACTER(LEN=120) :: NAME
      LOGICAL :: OPN
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     MOVEP LOCATES THE NTH-DATASET CONTAINED ON THE ALCHEMY CI
!           DUMPFILE (for polyatomics) ASSOCIATED WITH UNIT NFT
!
!           NPFLG = 0     NO PRINTOUT DURING SEARCH
!                 = 1     PRINT HEADER LABELS ENCOUNTERED
!           NALM  = 0     NO ERRORS DETECTED
!                 = 1     ERRORS DETECTED
!
!***********************************************************************
!
!
      NALM=0
      INQUIRE(UNIT=nft,OPENED=opn)
      IF(opn)THEN
         REWIND NFT
      ELSE
         OPEN(UNIT=nft,FORM='unformatted')
      END IF
!       IF(NTH.EQ.1)RETURN
!       M=NTH-1
      NTH=0
      IF(NTH.EQ.0)M=2000
      DO I=1, M
         READ(NFT,END=50)header
         READ(NFT,END=50)NT, nrec, NAME, nnuc, nocsf, nstat, mgvn, s,&
     &                   sz, nelt, e0

         if (gflag == 0) then
             if ( nt .eq. itarg) then
                 NTH=I
                 return
             end if
         else if ((mgvn .eq. mgvn_in) .and. ((2*s+1) .eq. spin_in)) then
             NTH=I
             return
         end if

         DO K=1, NNUC
            READ(NFT)
         END DO
!
         NEIG=MIN(MEIG,NSTAT)
         READ(NFT)(IDUM,J=1,nocsf), (EIG(J),J=1,NEIG)
!
         DO K=1, NSTAT
            READ(NFT)
         END DO
!
      END DO
      RETURN
!
 50   NALM=1
      RETURN
      END SUBROUTINE FINDP

      subroutine read_all_ci_vectors(nftw,ntgsym,numtgt,itarget_symmetry_order,&
     &                               all_ntgsym,all_numtgt, all_ci_vec_target,&
     &                               iwrite,gflag)
      use cdenprop_defs
      implicit none
!
!     **********************************************************************
!
!     read_all_ci_energies reads CI vectors from unit NFTW in
!     format used for polyatomic targets.
!
!     **********************************************************************

!     Arguments
      integer :: nftw,ntgsym,iwrite
      type (CIvect), allocatable, dimension(:) :: all_ci_vec_target
      integer, dimension(:) :: numtgt
      integer, dimension(:,:) :: itarget_symmetry_order

!     Local variables
      integer, parameter :: meig=20
      integer :: nciset,all_ntgsym
      integer, allocatable :: all_numtgt(:)
      logical, allocatable :: not_in_np1_congen(:)
      real(kind=idp) :: e0, s, sz,spin_temp
      real(kind=idp), dimension(meig) :: eig
      character(len=32) :: header
      integer :: i, idum, j, k, m, mgvn, neig, nelt, nnuc, nocsf, nrec, &
     &           nstat, nt, itarg,ifail,gflag
      character(len=120) :: name
      logical :: opn
      type (CIvect) :: ci_vec_target_temp

      intent(in)  :: nftw,itarget_symmetry_order,ntgsym,numtgt,gflag
      intent(out) :: all_ntgsym,all_numtgt,all_ci_vec_target

      inquire(unit=nftw,opened=opn)
      if(opn)then
         rewind nftw
      else
         open(unit=nftw,form='unformatted')
      end if

!     First we determine the number of symmetries.
      nciset=0
      do
         read(nftw,end=50) header
         read(nftw,end=50) nt, nrec, name, nnuc, nocsf, nstat, mgvn, s,&
     &                       sz, nelt, e0
         do k=1, nnuc
            read(nftw)

         end do
         neig=min(meig,nstat)
         read(nftw)(idum,j=1,nocsf), (eig(j),j=1,neig)
!
         do k=1, nstat
            read(nftw)

         end do

         nciset=nciset+1

      end do

 50   rewind nftw
      all_ntgsym=nciset
      allocate(all_numtgt(all_ntgsym))
      all_numtgt=0
      allocate(not_in_np1_congen(all_ntgsym))
      not_in_np1_congen=.TRUE.

      allocate(all_ci_vec_target(all_ntgsym))

!     Read symmetries used in N+1 congen first (in the same order as N+1 congen).
      do itarg=1,ntgsym
         ifail=0
         nciset=0
         spin_temp=real(itarget_symmetry_order(itarg,2))
         CALL FINDP(nftw,nciset,ifail,0,6,itarget_symmetry_order(itarg,1)-1,spin_temp,gflag,itarg)
         if (ifail .ne. 0) then
            write(iwrite,'(" MGVN = ",i2," SPIN = ",f4.1, " NOT FOUND ON UNIT ", i8)') &
                itarget_symmetry_order(itarg,1)-1, spin_temp,nftw
            stop ' ERROR'
         end if

         call read_ci_vector(nftw,nciset,ci_vec_target_temp,0,iwrite)

!        Multiply in the phase correction factor
         do i=1, ci_vec_target_temp%nstat
           do j=1, ci_vec_target_temp%nocsf
             ci_vec_target_temp%CV(j,i)=ci_vec_target_temp%CV(j,i)*ci_vec_target_temp%iphz(j)
           end do
         end do

         all_ci_vec_target(itarg)=ci_vec_target_temp
         not_in_np1_congen(nciset)=.FALSE.
         all_numtgt(itarg)=ci_vec_target_temp%nstat
      end do

!     Now read the remaining sets if there are any. (NOTE: THESE ARE NOT PRESENTLY PHASE CORRECTED)
      if (ntgsym .lt. all_ntgsym) then
         itarg=ntgsym+1
         do nciset=1,all_ntgsym
            if (not_in_np1_congen(nciset)) then
               call read_ci_vector(nftw,nciset,ci_vec_target_temp,0,iwrite)
               all_ci_vec_target(itarg)=ci_vec_target_temp
               all_numtgt(itarg)=ci_vec_target_temp%nstat
               itarg=itarg+1

            end if

         end do

      end if

      return
60    stop ' ERROR READING TARGET VECTORS'
      end subroutine read_all_ci_vectors

      subroutine READCIP (nftw, nciset, name, nocsf, nstat, nctgt, mgvn, s, sz, nelt, e0, EIG, VEC, iphase, NFT, cname, charge, &
                          nnuc, xnuc, ynuc, znuc, nstat_read)
      use cdenprop_defs
      implicit none
!
!     **********************************************************************
!
!     READCIP reads CI data from unit NFTW in format used for
!     polyatomic targets
!
!     **********************************************************************
!
!     Input variables - AH
!
!     nftw = unit no. ; nciset = set no. ; NFT = output summary file
!     nctgt = no. csf (must be equal to nocsf)? ; nstat1= No of targ. states per sym.
!
!     Output Variables
!
!     nocsf = no. csf ; mgvn = symmetry q. no. ; s & sz = spin
!     nelt = no. elec. ; EIG = array of ?
!     VEC = array of ? ; iphase = array of phases
!     ***********************************************************************

!     Arguments

      integer, parameter ::MEIG=20,maxsym=8
      integer :: nftw,nciset,nocsf,nstat,nctgt,mgvn,nelt,NFT,nnuc,nstat_read
      real(kind=idp) :: s,sz
      character(len=8), dimension(:) :: cname(maxnuc)
      character(len=120) :: NAME
      real(kind=idp), dimension(:) :: xnuc(maxnuc),ynuc(maxnuc),znuc(maxnuc),charge(maxnuc)

      integer, allocatable, dimension(:) :: iphase
      real(kind=idp),allocatable, dimension(:) :: EIG
      real(kind=idp),allocatable, dimension(:,:) :: VEC
      real(kind=idp),allocatable :: temp(:)


!     Local variables
      integer :: i, nalm, nrec, nset, neig, nsr
      real(kind=idp) :: e0

!
!     Find file
      CALL MOVEP (NFTW,nciset,NALM,0,NFT)
      IF (NALM.NE.0) GO TO 1800
!
      read(nftw)
      READ(NFTW) Nset,nrec,NAME,nnuc,nocsf,nstat,mgvn,s,sz,nelt,e0
!       print *, 'MGVN', MGVN, nstat, nocsf

      do  i=1,nnuc
         read(nftw) cname(i),xnuc(i),ynuc(i),znuc(i),charge(i)
      end do

      if(nset.ne.nciset) go to 1800

!     READ CI COEFFICIENTS
!ZM read only those that I need
      nsr = nstat_read
      if (nsr > nstat) then
         print *,'incorrect nstat_read on input',nstat_read,nstat
         stop
      endif
      if (nsr .le. 0) nsr = nstat
!
      if (nsr .ne. nstat) then
         write(NFT,'(1X,"READING ONLY ",i10," EIGENPAIRS.")') nstat_read
      endif

      allocate(iphase(nocsf),eig(nstat),vec(nocsf,nsr))
      EIG=0d0; VEC=0d0; iphase=0

      CALL CIVIO (NFTW,1,NOCSF,nstat,nsr,EIG,VEC,NALM,iphase,vec)
      IF (NALM.NE.0) GO TO 2900
      NCTGT=nocsf !For N+1 calcualtions the below comparison is not true as CSF file is prototype basis
!     and CI vectors in terms of contracted basis. Consistancy checking to be done elsewhere.*******
      IF (NOCSF.EQ.NCTGT) THEN
!        Print summary of target data
         WRITE (NFT,100)  NSET,name
         WRITE (NFT,101)  nocsf,nstat,nnuc
         WRITE (NFT,120) (cname(i),xnuc(i),ynuc(i),znuc(i), &
    &                   charge(i),i=1,nnuc)
         NEIG = MIN(MEIG,nstat)
         WRITE (NFT,102) (EIG(I)+E0,I=1,NEIG)

      ELSE
!        Error
         WRITE(NFT,2920) NOCSF,NCTGT
         stop

      end if
!ZM read only those that I need
      nstat = nsr
!
      if (nstat .ne. nsr) then
         call move_alloc(eig,temp)
         allocate(eig(nsr))
         eig(1:nsr) = temp(1:nsr)
         deallocate(temp)
      endif

      write(NFT,'(/,1X,"READCIP finished",/)')

      return

!     Format statements

 1800 WRITE (nft,1804) nftw
 1804 FORMAT(/' CIDATA NOT FOUND ON UNIT',I3)
      STOP
!
 2900 WRITE(NFT,2910)
 2910 FORMAT(/' UNABLE TO GET CI-TARGET VECTOR ')
      stop
!
  100 FORMAT (' SET',I4,4X,A)
  101 FORMAT (10X,'NOCSF=',I10,4X,'NSTAT=',I10,4X,'NNUC =',I5)
  102 FORMAT (10X,'EIGEN-ENERGIES',/(16X,5F20.10))
  120 FORMAT(/' Nuclear data     X         Y         Z       Charge' &
    & /(3x,a8,2x,4f10.6))
 2920 FORMAT(/' HAMILTONIAN TRANSFORMATION DATA INCONSISTENT ',2I10)
!
      end subroutine READCIP

      subroutine read_property_integrals2(nftprop,pintegrals,iwrite,ukrmolp_ints)
!     Wrapper for read_property_integrals
      use cdenprop_defs
      use ukrmol_interface_gbl, only: construct_pintegrals
      implicit none
      integer :: nftprop,iwrite
      type (property_integrals) :: pintegrals
      logical, intent(in) :: ukrmolp_ints

      integer :: i,j

      if (ukrmolp_ints) then !UKRmol+ property integrals

         !Read all basis sets and integrals into memory and convert the UKRmol+ property integrals into the GAUSPROP storage format.
         call construct_pintegrals(nftprop,iwrite,pintegrals%no_of_integrals,pintegrals%no_of_properties,&
     &                             pintegrals%nilmq,pintegrals%lp,pintegrals%mp,pintegrals%qp,&
     &                             pintegrals%property_name,pintegrals%nob,pintegrals%mob,pintegrals%mpob)

         pintegrals%ukrmolp_ints = .true.

      else !SWEDEN property integrals

         call read_property_integrals(nftprop,pintegrals%no_of_integrals,pintegrals%no_of_properties,&
     &                                pintegrals%nilmq,pintegrals%lp,pintegrals%mp,pintegrals%qp,&
     &                                pintegrals%property_name,pintegrals%nob,pintegrals%mob,pintegrals%mpob,pintegrals%indexv,&
     &                                pintegrals%xintegrals,iwrite)

         allocate(pintegrals%inverted_indexv(pintegrals%no_of_integrals,pintegrals%no_of_properties))
         pintegrals%inverted_indexv=0
         pintegrals%ukrmolp_ints = .false.

         do i=1,pintegrals%no_of_properties
            do j=1,pintegrals%no_of_integrals
               if (pintegrals%indexv(j,i).ne. 0) then
                  pintegrals%inverted_indexv(pintegrals%indexv(j,i),i)=j

               end if

            end do

         end do
      endif

      end subroutine read_property_integrals2



      SUBROUTINE read_property_integrals(nftprop,NINTS,nlmq,nilmq,lp,mp,qp,pname,nob,mob,mpob,indexv,XBUF,iwrite)
!***********************************************************************
!
!     GPWRIT - GP WRITe  computed property integrals to disk
!
!     Input data:
!        NFTPROP  Logical unit for the disk file on which the integrals
!                 will be written
!         IWRITE  Logical unit for the printer
!          NINTS  Number of one electron integrals in the buffers
!         INDEXV  Index defining location of each integral in the
!                 lower triangle of ordered integrals
!
!     Output data:
!
!            None
!
!***********************************************************************
      USE cdenprop_defs
      IMPLICIT NONE

      INTEGER :: IWRITE, NFTPROP, NINTS, NLMQ, NNUC, NSYM, NTRANSF,lbuf
      CHARACTER(LEN=132) :: NAME
      REAL(KIND=idp), allocatable, DIMENSION(:) :: CHARG, XNUC, YNUC, ZNUC
      INTEGER, allocatable, DIMENSION(:,:) :: INDEXV
      INTEGER, allocatable, DIMENSION(:) :: LP, MP, NILMQ, QP
      INTEGER, allocatable, DIMENSION(:) :: NOB
      CHARACTER(LEN=8), allocatable, DIMENSION(:) :: PNAME
      REAL(KIND=idp), allocatable, DIMENSION(:,:) :: XBUF
      REAL(KIND=idp) :: EN
      INTEGER :: I, IFINISH, ISTART, J, K, L, LMQ,lmq1, M, N, NCODT, NOBT,&
     &           NPASS, NREC, NREM
      INTEGER, allocatable, DIMENSION(:) :: MOB, MPOB
!
      READ(nftprop)NAME, Nsym, nobt, NNUC, NINTS, lbuf

      allocate(CHARG(nnuc), XNUC(nnuc), YNUC(nnuc), ZNUC(nnuc))
      allocate(nob(nsym), mob(nobt), mpob(nobt))
      rewind(nftprop)
      READ(nftprop)NAME, Nsym, nobt, NNUC, NINTS, lbuf, &
     &              (Nob(I),I=1,Nsym), (mpob(I),I=1,Nobt),&
     &              (mob(I),I=1,Nobt), (CHARG(I),I=1,NNUC), EN,&
     &              (xnuc(i),I=1,NNUC), (ynuc(i),I=1,NNUC),&
     &              (znuc(i),I=1,NNUC)
!
      READ(NFTprop)Nlmq, ntransf


      allocate(LP(nlmq), MP(nlmq), NILMQ(nlmq), QP(nlmq),pname(nlmq))
      allocate(indexv(nints,nlmq),XBUF(nints,nlmq))

!
!---- Loop over all properties
      DO lmq=1, nlmq
!
!---- How many passes over the buffers will we have to make ?
!
         READ(nftprop)pname(lmq), lmq1, nrec, nilmq(lmq), lp(lmq),&
     &                 mp(lmq), qp(lmq)

         NPASS=Nilmq(lmq)/LBUF
         NREM=Nilmq(lmq)-NPASS*LBUF
         nrec=npass
         IF(nrem.GT.0) THEN
            nrec=nrec+1

         END IF

         ISTART=1
         DO I=1, NPASS
            IFINISH=ISTART+LBUF-1
            READ(nftprop)(Xbuf(J,lmq),J=ISTART,IFINISH),&
     &                    (indexv(j,lmq),J=ISTART,IFINISH), LBUF
            ISTART=IFINISH+1

         END DO
!
         IFINISH=ISTART+NREM-1
         IF(nrem.GT.0)READ(nftprop)(XBUF(J,lmq),J=ISTART,IFINISH),&
     &                              (indexv(j,lmq),J=ISTART,IFINISH),&
     &                              nrem
!
         WRITE(IWRITE,5050)pname(lmq), nilmq(lmq)
!
      END DO
!

      RETURN
!
 5000 FORMAT(4(i5,d15.6))
 5001 FORMAT(//1x,a8,6I5)
 5050 FORMAT(/' Number of property integrals of type ',a8,' read',&
     &       ' is',i10)
      END SUBROUTINE  read_property_integrals

      !In preparation for changes to the rmat I/O
      !******************************************

      subroutine read_civ (nftw, nciset, name, nocsf, nstat, nctgt, mgvn, s, sz, nelt, e0, EIG, VEC, iphase, NFT, cname, &
                           charge, nnuc, xnuc, ynuc, znuc, nstat_read)
      use cdenprop_defs
      implicit none
      integer, intent(in) :: NFT
      integer, intent(inout) :: nftw,nciset,nocsf,nctgt,mgvn,nelt,nnuc,nstat_read
      integer, intent(out) :: nstat
      character(len=120), intent(out) :: NAME
      real(kind=idp), intent(inout) :: s, sz
      integer,allocatable, dimension(:) :: iphase
      real(kind=idp), allocatable, dimension(:) :: EIG
      real(kind=idp),allocatable, dimension(:,:) :: VEC
      real(kind=idp), dimension(:) :: charge,xnuc,ynuc,znuc
      character(len=8), dimension(:) :: cname
      real(kind=idp) :: e0

      call READCIP (nftw, nciset, name, nocsf, nstat, nctgt, mgvn, s, sz, nelt, e0, EIG, VEC, iphase, NFT, cname, &
                    charge, nnuc, xnuc, ynuc, znuc, nstat_read)

      end subroutine read_civ


      subroutine read_ci_vector(nftw,nciset,ci_vec,nstat_read,iwrite)
      use cdenprop_defs
      implicit none
      type(CIvect) :: a, ci_vec
      integer :: nftw,nciset, nocsf, nstat,nstat_read,iwrite

      call ci_vec % dealloc
      call read_civ (nftw, nciset, ci_vec % name, ci_vec % nocsf, ci_vec % nstat, ci_vec % nocsf, ci_vec % mgvn, ci_vec % s, &
                     ci_vec % sz, ci_vec % nelt, ci_vec % e0, ci_vec % EI, ci_vec % CV, ci_vec % iphz, iwrite, ci_vec % cname, &
                     ci_vec % charge, ci_vec % nnuc, ci_vec % xnuc, ci_vec % ynuc, ci_vec % znuc, nstat_read)

      end subroutine read_ci_vector

!     ------------------------------------------------
!
!     READ BOUND STATE COEFFICIENTS FROM OUTER REGION
!     ROUTINE BOUND
!
!     ------------------------------------------------
      subroutine read_states_from_bound(lubnd,nbset,bform,states_from_bound,iwrite)
      use cdenprop_defs
      implicit none

!     Arguments
      integer :: lubnd, nbset, iwrite
      type(CIvect) :: states_from_bound
      character(len=11) :: bform

!     Local variables
      integer :: nchan
      real(kind=idp), allocatable, dimension(:):: vtemp, xvec

      call states_from_bound % dealloc
      call read_boundcoeff(lubnd,nbset,bform,states_from_bound%nocsf,nchan,states_from_bound%nstat,&
                           states_from_bound%ei,vtemp,states_from_bound%CV, states_from_bound%mgvn, xvec, iwrite)


      end subroutine read_states_from_bound

      subroutine read_boundcoeff(lubnd, nbset, bform, nstat, nchan, nbound, etot,vtemp,bcoef, mgvn, xvec, iwrit)
      use cdenprop_defs
      implicit none

!     Arguments
      integer :: lubnd, nbset, nstat, nchan, nbound, mgvn
      character(len=11) :: bform
      real(kind=idp), allocatable, dimension(:):: etot, vtemp, coef, xvec
      real(kind=idp),  allocatable, dimension(:,:) :: bcoef

!     Local variables
      integer :: stot, gutot, iprint, iwrit, ifail, iprnt
      real(kind=idp) :: RR
      integer :: i,j
      real(kind=idp) :: sumb

      intent(in) :: lubnd, nbset,bform
      intent(out) :: nstat, nchan,nbound
      intent(inout) :: etot,vtemp,bcoef, xvec

      iprnt=0

      call readbh(lubnd,nbset,nchan,mgvn,stot,gutot,nstat,nbound,rr,bform,iprnt,iwrit,ifail)
      if (ifail /= 0) stop "read_boundcoeff: error in readbh"

      allocate(etot(nbound),vtemp(nbound),coef(nstat*nbound),xvec(nchan*nbound),bcoef(nstat,nbound))

      call readbc(lubnd,nstat,etot,vtemp,coef,nbound,nchan,xvec,bform,iprnt,iwrit)

!     Unpack coefficients
      do i=1, nbound
        bcoef(:,i)=coef((i-1)*nstat+1:(i-1)*nstat+nstat)
      end do
!     Write state information to standard output
      write(iwrit,'(/,"  Bound states read from unit ",i6)') lubnd
      write(iwrit,'(  "  Symmetry =  ",i2,"   Spin = ",i2)') mgvn,stot
      write(iwrit,'(  "  Basis size =  ",i8)') nstat
      write(iwrit,'(  "  No. bound states =  ",i8,/)') nbound
      end subroutine read_boundcoeff

!     -------------------------------------------------
!
!     READBH AND GETSET TAKEN FROM OUTER REGION
!
!
!     -------------------------------------------------

      SUBROUTINE READBH(LUBND0,NBSET,NCHAN,MGVN,STOT,GUTOT,NSTAT,&
     &NBOUND,RR,BFORM0,IPRNT0,IWRIT0,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!***********************************************************************
!
!     READBH locates and reads bound state set NBSET on unit LUBND
!
!***********************************************************************

      INTEGER STOT,GUTOT
      CHARACTER(LEN=11) BFORM,BFORM0
      CHARACTER(LEN=80) TITLE
      DATA KEYBC/11/
      SAVE
!
      IFAIL=0
      BFORM = BFORM0
      LUBND = LUBND0
      IPRNT = IPRNT0
      IWRITE= IWRIT0

!---- Locate set number NTSET on unit LUBND
      NSET = NBSET
      CALL GETSET(LUBND,NSET,KEYBC,BFORM,IFAIL)
      IF(IFAIL.NE.0) GO TO 99

!---- Read header
      IF(BFORM.EQ.'FORMATTED') THEN
        READ(LUBND,10)  KEYBC,NSET
        READ(LUBND,13) TITLE
        READ(LUBND,10) NBOUND,MGVN,STOT,GUTOT,NSTAT,NCHAN
        READ(LUBND,11) RR
      ELSE
        READ(LUBND) KEYBC,NSET
        READ(LUBND) TITLE
        READ(LUBND) NBOUND,MGVN,STOT,GUTOT,NSTAT,NCHAN
        READ(LUBND) RR
      ENDIF

!---- Print header information
      IF(IPRNT.NE.0)  THEN
        WRITE(IWRITE,14)
        WRITE(IWRITE,100) KEYBC,NSET
        WRITE(IWRITE,130) TITLE
        WRITE(IWRITE,100) NBOUND,MGVN,STOT,GUTOT,NSTAT,NCHAN
        WRITE(IWRITE,110) RR
      ENDIF

      RETURN

 99   WRITE(IWRITE,98) NBSET,LUBND
      IFAIL = 1
      RETURN

 11   FORMAT(10F20.13)
 10   FORMAT(10I5)
 13   FORMAT(A80)
 14   FORMAT(/' Header on LUBND')
 100  FORMAT(1X,10I5)
 110  FORMAT(1X,4E20.13)
 130  FORMAT(1X,A80)
 98   FORMAT(/' UNABLE TO FIND BOUND STATE SET',I3,' ON UNIT',I3)
      END SUBROUTINE READBH

      SUBROUTINE READBC(LUBND,NSTAT,ETOT,VTEMP,COEF,NBOUND,NCHAN,XVEC,BFORM,IPRNT,IWRITE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(LEN=11) BFORM
      DIMENSION COEF(*),ETOT(*),VTEMP(*),XVEC(*)

      DO I=1,NBOUND
        IF(BFORM.EQ.'FORMATTED') THEN
          READ(LUBND,11) ETOT(I),VTEMP(I),(COEF((I-1)*NSTAT+J),&
     &J=1,NSTAT)
          READ(LUBND,11) (XVEC((I-1)*NCHAN+J),J=1,NCHAN)
        ELSE
         READ(LUBND) ETOT(I),VTEMP(I),(COEF((I-1)*NSTAT+J),J=1,NSTAT)
         READ(LUBND) (XVEC((I-1)*NCHAN+J),J=1,NCHAN)
        ENDIF
        IF(IPRNT.NE.0) THEN
           WRITE(IWRITE,110) ETOT(I),VTEMP(I),&
     &                       (COEF((I-1)*NSTAT+J),J=1,NSTAT)
           WRITE(IWRITE,110) (XVEC((I-1)*NCHAN+J),J=1,NCHAN)
        ENDIF
      END DO

      RETURN

 11   FORMAT(10F20.13)
 110  FORMAT(1X,4E20.13)
      END SUBROUTINE READBC

      SUBROUTINE GETSET(LUNIT,NSET,INKEY,FORM,IFAIL)

!***********************************************************************
!
!     GETSET locates set number NSET on unit LUNIT
!     If NSET = 0 file is positioned at end-of-information
!             = 1 the file is opened
!             = n file is positioned at the beginning of set number n
!     On return NSET = sequence number of current set
!
!***********************************************************************

      LOGICAL OP
      CHARACTER(LEN=11) FORM
      CHARACTER(LEN=80) LINE

      INQUIRE(UNIT=LUNIT,OPENED=OP)
      IF(.NOT.OP) OPEN(UNIT=LUNIT,ERR=99,FORM=FORM,status='unknown')
      IF(NSET.EQ.1) THEN
        REWIND(UNIT=LUNIT,ERR=100)
 100    RETURN
      ELSE
        REWIND LUNIT
        INSET = NSET
        ITIME = 1

!------ LOCATE SET NUMBER INSET
 5      IF(FORM.EQ.'FORMATTED') THEN
 15       READ(LUNIT,*,END=9) KEY,NSET,NREC
          IF(NSET.EQ.INSET.AND.KEY.EQ.INKEY) THEN
            BACKSPACE LUNIT
            RETURN
          ELSE
            DO I=1,NREC
              READ(LUNIT,*,END=99)
            END DO
          ENDIF
          GO TO 15
        ELSE
 25       READ(LUNIT,END=9) KEY,NSET,NREC
          IF(NSET.EQ.INSET.AND.KEY.EQ.INKEY) THEN
            BACKSPACE LUNIT
            RETURN
          ELSE
            DO I=1,NREC
              READ(LUNIT,END=99)
            END DO
          ENDIF
        GO TO 25
        ENDIF

!------ END OF FILE HAS BEEN REACHED
 9      IF(INSET.EQ.0) THEN
          BACKSPACE LUNIT
          RETURN
        ELSE IF(ITIME.EQ.1) THEN
          ITIME = 2
          REWIND LUNIT
          GO TO 5
        ENDIF
      ENDIF
!
 99   IFAIL = 1
      RETURN
      END SUBROUTINE GETSET

      subroutine read_csf_head(megul,csf_header,iwrite)
      use cdenprop_defs
      implicit none
      type (CSFheader) :: csf_header
      integer :: megul,iwrite
      integer :: i,m,d

      !Deallocate derived type allocatable arrays (should be a nicer way to do this.)
      !------------------------------------------------
      call csf_header % dealloc
      !---------------------------------------------------

      read(megul) csf_header%NAME,csf_header%MGVN,csf_header%S,csf_header%SZ,csf_header%R,csf_header%PIN,csf_header%NORB, &
                  & csf_header%NSRB, csf_header%NOCSF,csf_header%NELT,csf_header%lcdof,csf_header%IDIAG,csf_header%NSYM,&
                  & csf_header%SYMTYP,csf_header%lndof, csf_header%npflg,csf_header%thres,csf_header%nctarg,csf_header%ntgsym

      !ZM prevent wrong allocation and call to RDNFTO below if no target state symmetries are found on the CSF file.
      d = max(csf_header%ntgsym,1)

      !Write file header - This should become a function in the io module

      write(iwrite,"(' ',A6,': ',A)") 'SCATCI',csf_header%NAME
      !if (idiag.lt.0) idiag=idiagt
      write(iwrite,16) csf_header%MGVN,csf_header%S,csf_header%NELT,csf_header%SZ,csf_header%IDIAG,csf_header%R,csf_header%THRES,&
                  & csf_header%NSYM,MEGUL,csf_header%NOCSF,(csf_header%NPFLG(I),I=1,6)
      16 FORMAT(' MGVN =',I10,5X,'S    =',F5.1,/,' NELT =',I10,5X,'SZ   =' &
          &       ,F5.1,/,' IDIAGT=',I10,5X,'R    =',F5.1,/, &
          &       ' THRES=',D10.1,5X,'NSYM =',I5,/,' MEGUL=',I3,', NOCSF=' &
          &       ,I10,/,' NPFLG',I10,5I4)

      allocate(csf_header % NOB(csf_header % NSYM), &
               csf_header % NDTRF(csf_header % NELT), &
               csf_header % NODO(csf_header % NOCSF), &
               csf_header % iphz(csf_header % nctarg),&
               csf_header % nobl(2 * csf_header % nsym), &
               csf_header % nob0(csf_header % nsym), &
               csf_header % no_l2_virtuals(csf_header % nsym), &
               csf_header % nob0l(2 * csf_header % nsym), &
               csf_header % nctgt(d), &
               csf_header % notgt(d), &
               csf_header % mcont(d), &
               csf_header % gucont(d))

      csf_header%no_l2_virtuals=0


      call RDNFTO(megul,csf_header%nob,csf_header%nob0,csf_header%nobl,csf_header%nob0l,csf_header%nsym,csf_header%ndtrf, &
                  & csf_header%nelt,csf_header%nodo,csf_header%nocsf,csf_header%iphz,csf_header%nctarg,csf_header%nctgt, &
                  & csf_header%notgt,csf_header%mcont,csf_header%gucont,d,csf_header%ntgsym)

      end subroutine read_csf_head


    subroutine read_csf_body (megul, csf_header, csf_body, iwrite)

        use cdenprop_defs,  only: CSFheader, CSFbody
        use mpi_gbl,        only: shared_communicator, local_rank, master, shared_enabled
        use mpi_memory_gbl, only: mpi_memory_allocate_integer, mpi_memory_allocate_real, mpi_memory_synchronize

        type(CSFbody),   intent(inout) :: csf_body
        type(CSFheader), intent(in)    :: csf_header
        integer,         intent(in)    :: megul, iwrite

        call csf_body % dealloc

        ! allocate common memory area for all processes on a single node
        csf_body % indo_window = mpi_memory_allocate_integer(csf_body % indo, csf_header % nocsf + 1, shared_communicator)
        csf_body % icdo_window = mpi_memory_allocate_integer(csf_body % icdo, csf_header % nocsf + 1, shared_communicator)
        csf_body % ndo_window  = mpi_memory_allocate_integer(csf_body % ndo,  csf_header % lndof,     shared_communicator)
        csf_body % cdo_window  = mpi_memory_allocate_real   (csf_body % cdo,  csf_header % lcdof,     shared_communicator)

        ! let the local master process read the CSFs
        if (local_rank == master .or. .not. shared_enabled) then
            call rdwf(csf_header % nocsf + 1, &
                      csf_body % icdo,    &
                      csf_body % indo,    &
                      csf_body % ndo,     &
                      csf_header % lndof, &
                      csf_body % cdo,     &
                      csf_header % lcdof, &
                      megul)
        end if

        ! synchronize with other processes
        call mpi_memory_synchronize(csf_body % indo_window, shared_communicator)
        call mpi_memory_synchronize(csf_body % icdo_window, shared_communicator)
        call mpi_memory_synchronize(csf_body % ndo_window,  shared_communicator)
        call mpi_memory_synchronize(csf_body % cdo_window,  shared_communicator)

    end subroutine read_csf_body


      subroutine write_csf(unit,csf_header, csf_body)
      use cdenprop_defs
      implicit none
      type (CSFheader) :: csf_header
      type (CSFbody) :: csf_body
      integer :: unit, nfto,M,nx

      nfto=unit
      M=csf_header%NOCSF+1
      nx= size(csf_header%nobl)

      call WRNFTO (csf_header % NAME, csf_header % MGVN, csf_header % S, csf_header % SZ, csf_header % R, csf_header % PIN, &
                   csf_header % NORB, csf_header % NSRB, csf_header % NOCSF, csf_header % NELT, csf_header % IDIAG, &
                   csf_header % NSYM, csf_header % SYMTYP, csf_header % NOB, csf_header % NDTRF, csf_header % NODO, M, &
                   csf_body % ICDO, csf_body % INDO, csf_body % NDO, csf_header % lndof, csf_body % CDO, csf_header % lcdof, &
                   NFTO, csf_header % nobl, nx, csf_header % npflg, csf_header % thres, csf_header % iposit, csf_header % nob0, &
                   csf_header % nob0l, csf_header % nctarg, csf_header % ntgsym, csf_header % notgt, csf_header % nctgt, &
                   csf_header % mcont, csf_header % gucont, csf_header % iphz)

      ! nx, m, are potential problems

      end subroutine write_csf


    !> \brief  Read transition dipoles
    !> \author A Harvey, J Benda
    !> \date   2017 - 2021
    !>
    !> Read the complete target property file and select those properties relevant for the CI target states used in CDENPROP
    !> (via `numtgt`). Uses `molecular_properties_data`.
    !>
    !> \param[in]    iwrite  Output unit for text output.
    !> \param[in]    iutdip  Logical unit associated with the property file.
    !> \param[in]    ntset   Dataset number in the property file.
    !> \param[in]    nstat   Total number of CI target states.
    !> \param[in]    numtgt  Number of CI target states per spin-symmetry used in CDENPROP.
    !> \param[in]    itarget_symmetry_order  Irreducible representation and spin multiplicity per CI target spin-symmetry.
    !> \param[in]    nnuc    Number of nuclei.
    !> \param[in]    ismax   Highest angular momentum to read (not used).
    !> \param[inout] prop    On exit: allocated and populated matrix of properties (per l*l+l+m+1 property index).
    !> \param[out]   ifail   Success exit indicator.
    !>
    subroutine read_transdip2 (iwrite, iutdip, ntset, nstat, numtgt, itarget_symmetry_order, nnuc, ismax, prop, gflag, prop_order, &
                               ifail)

        use cdenprop_defs, only: idp, maxprop
        use class_molecular_properties_data, only: molecular_properties_data
        use mpi_gbl, only: mpi_xermsg

        implicit none

        integer, intent(in)  :: iwrite, iutdip, ntset, nstat, nnuc, numtgt(:), itarget_symmetry_order(:,:), ismax, gflag
        integer, intent(in)  :: prop_order(:)
        integer, intent(out) :: ifail

        real(kind=idp), allocatable, intent(inout) :: prop(:, :, :)

        type(molecular_properties_data) :: target_properties

        integer :: itgt, iidx, itgt_idx, itgt_irr, itgt_spn, itgt_sym
        integer :: jtgt, jidx, jtgt_idx, jtgt_irr, jtgt_spn, jtgt_sym
        integer :: i, j, l, m, lm

        ifail = 0

        write (iwrite, *) 'READING PROPERTIES FILE'

        call target_properties % read_properties(iutdip, ntset)   ! read all properties from the file

        if (target_properties % nuclear_prop_flag == 0) then
            call mpi_xermsg('cdenprop_io', 'read_transdip2', &
                            'Target properties read by CDENPROP seem to contain nuclear multipole (isw = 0)!' // new_line(' ') // &
                            ' *  But CDENPROP adds this contribution on its own, separately. Is the input correct?', 0, 0)
        end if

        allocate (prop(nstat, nstat, (maxprop + 1)**2))

        prop = 0

        ! loop over all (non-zero) properties and pick the needed ones
        do i = 1, target_properties % non_zero_properties

            ! get bra/ket state for this property as well as the l,m quantum numbers
            itgt = target_properties % properties_index(1, i)                  ! energy order of bra state
            jtgt = target_properties % properties_index(2, i)                  ! energy order of ket state
            l    = target_properties % properties_index(3, i)                  ! angular momentum
            m    = target_properties % properties_index(4, i)                  ! angular momentum component

            ! get quantum numbers of the bra state
            if (gflag == 0) then
                itgt_idx = prop_order(itgt)                                    ! index within spin-symmetry
            else
                itgt_idx = target_properties % energies_index(1, itgt)         ! index within spin-symmetry
            end if
            itgt_sym = target_properties % energies_index(2, itgt)             ! spin-symmetry index
            itgt_irr = target_properties % symmetries_index(1, itgt_sym) + 1   ! irreducible representation (mvgn + 1)
            itgt_spn = target_properties % symmetries_index(2, itgt_sym)       ! spin multiplicity

            ! get quantum numbers of the ket state
            if (gflag == 0) then
                jtgt_idx = prop_order(jtgt)                                    ! index within spin-symmetry
            else
                jtgt_idx = target_properties % energies_index(1, jtgt)         ! index within spin-symmetry
            end if
            jtgt_sym = target_properties % energies_index(2, jtgt)             ! spin-symmetry index
            jtgt_irr = target_properties % symmetries_index(1, jtgt_sym) + 1   ! irreducible representation (mgvn + 1)
            jtgt_spn = target_properties % symmetries_index(2, jtgt_sym)       ! spin multiplicity

            ! find out if both these bra/ket target states are needed by CDENPROP
            if (gflag == 0) then
                iidx = itgt_idx
                jidx = jtgt_idx
            else
                iidx = -1  ! CONGEN order of bra state
                jidx = -1  ! CONGEN order of ket state
                do j = 1, size(itarget_symmetry_order, 1)
                    if (itarget_symmetry_order(j, 1) == itgt_irr .and. &
                        itarget_symmetry_order(j, 2) == itgt_spn .and. itgt_idx <= numtgt(j)) then
                        iidx = sum(numtgt(1:j-1)) + itgt_idx
                    end if
                    if (itarget_symmetry_order(j, 1) == jtgt_irr .and. &
                        itarget_symmetry_order(j, 2) == jtgt_spn .and. jtgt_idx <= numtgt(j)) then
                        jidx = sum(numtgt(1:j-1)) + jtgt_idx
                    end if
                end do
            end if

            ! store (and symmetrize) the required property
            if (iidx > 0 .and. jidx > 0 .and. l <= maxprop) then
                lm = l*l + l + m + 1
                prop(iidx, jidx, lm) = target_properties % properties(i)
                prop(jidx, iidx, lm) = target_properties % properties(i)
            end if

        end do

    end subroutine read_transdip2


      subroutine read_transformation_matrix(lutransmat_l2contract,l2block_contract)
      USE cdenprop_defs
      implicit none

!     Argument variables
      integer :: lutransmat_l2contract
      real(kind=idp),allocatable,dimension(:,:) :: l2block_contract

!     Local variables
      integer :: noriginal_l2_basis, ncontracted_l2_basis

      noriginal_l2_basis=size(l2block_contract,1)
      ncontracted_l2_basis=size(l2block_contract,2)

      read(lutransmat_l2contract) noriginal_l2_basis,ncontracted_l2_basis
      rewind(lutransmat_l2contract)

      allocate(l2block_contract(noriginal_l2_basis,ncontracted_l2_basis))

      read(lutransmat_l2contract) noriginal_l2_basis,ncontracted_l2_basis,l2block_contract


      end subroutine  read_transformation_matrix


!     ******************************************************************
!
!     OUTPUT ROUTINES
!
!     ******************************************************************

      subroutine write_properties (lupropw, max_multipole, nstat, csf_head_i, csf_head_j, ci_vec_i, ci_vec_j, &
                                   threshold, properties, pintegrals)
      use cdenprop_defs
      implicit none

!     Arguments
      integer :: lupropw,max_multipole
      integer, dimension(:) :: nstat
      type (CIvect) :: ci_vec_i,ci_vec_j
      type (CSFheader) :: csf_head_i,csf_head_j
      real(kind=idp),dimension(:,:,:) :: properties
      real(kind=idp) :: threshold
      type (property_integrals) :: pintegrals

!     Local variables
      integer :: non_zero_properties,istate,jstate, i,j,k, ielement
      character(len=36) :: str_description
      character(len=5) :: str_istate


!     Determine the number of non-zero dipoles
      non_zero_properties = 0
      do istate=1,nstat(1)
         do jstate=1,nstat(2)
            do i=1,(max_multipole+1)**2
                if (abs(properties(istate,jstate,i)).gt. threshold) then
                   non_zero_properties = non_zero_properties + 1
                end if

            end do

         end do

      end do

      write(6,'(/," Writing properties to unit ",i6)') lupropw

!     Write header information to file
!     ------------------------------------------------------------------
!       if ((nstat(1)+nstat(2)) .lt. 1000) then
!           write(lupropw,'(I1,I3,I9,I3,I4,I9,I3,1x,3D20.12)') 6,1,non_zero_properties + nstat(1) + nstat(2) + ci_vec_i%nnuc,ci_vec_i%nnuc ,nstat(1) + nstat(2),non_zero_properties,0,0.d0,0.d0,0.d0
!       else
          write(lupropw,'(I1,I3,I9,I9,I9,I9,I3,1x,3D20.12)') &
            6, 1, non_zero_properties + nstat(1) + nstat(2) + ci_vec_i % nnuc, ci_vec_i % nnuc, nstat(1) + nstat(2), &
            non_zero_properties, 0, 0.d0, 0.d0, 0.d0
!       end if

      do i = 1,ci_vec_i%nnuc
         ! Determine the element index
         if ((ci_vec_i%cname(i)(1:4) .eq. 'Scat') .or. (ci_vec_i%cname(i)(1:4) .eq. 'Pseu')) then
            write(lupropw,'(I1,I3,a3,i3,f10.4,3F20.10)') &
                8, i, ci_vec_i % cname(i)(1:2), int(ci_vec_i % charge(i)), 0.0_idp, &
                ci_vec_i % xnuc(i), ci_vec_i % ynuc(i), ci_vec_i % znuc(i)

         else
            do j=1,103
               if (ci_vec_i%cname(i)(1:2) .eq. ASYMB(j)) then
                  ielement=j

               end if

            end do

            write(lupropw,'(I1,I3,a3,i3,f10.4,3F20.10)') &
                8, i, ci_vec_i % cname(i)(1:2), int(ci_vec_i % charge(i)), AMASS(ielement), &
                ci_vec_i % xnuc(i), ci_vec_i % ynuc(i), ci_vec_i % znuc(i)

         end if

      end do

!     Write states and energies to file
!     ------------------------------------------------------------------
!     TODO SOMETHING WRONG WITH SPIN

      do istate = 1,nstat(1)
         write(str_istate,'(i5)') istate
         if(csf_head_i%nsym .eq. 4) then
            str_description="State No. "//str_istate //" "//  spin_symb(int(2*(csf_head_i%S)+1)) //cc2v_symb(csf_head_i%MGVN+1)

         else
            str_description="State No. "//str_istate //" "//  spin_symb(int(2*(csf_head_i%S)+1)) //cd2h_symb(csf_head_i%MGVN+1)

         end if

         write(lupropw,'(I1,I8,I3,I8,4I3,D20.12,2X,A36)') &
            5, istate, 0, 0, csf_head_i % MGVN, int(2 * (csf_head_i % S) + 1), 0, 0, &
            ci_vec_i % ei(istate) + ci_vec_i % e0, str_description

      enddo

      do jstate = 1,nstat(2)
         write(str_istate,'(i5)') jstate
         if(csf_head_i%nsym .eq. 4) then
            str_description="State No. "//str_istate //" "// spin_symb(int(2*(csf_head_j%S)+1)) //cc2v_symb(csf_head_j%MGVN+1)

         else
            str_description="State No. "//str_istate //" "//  spin_symb(int(2*(csf_head_j%S)+1)) //cd2h_symb(csf_head_j%MGVN+1)

         end if

         write(lupropw,'(I1,I8,I3,I8,4I3,D20.12,2X,A36)') &
            5, jstate + nstat(1), 0, 0, csf_head_j % MGVN, int(2 * (csf_head_j % S) + 1), 0, 0, &
            ci_vec_j % ei(jstate) + ci_vec_j % e0, str_description

      end do

!     Write properties to file
!     ------------------------------------------------------------------

      i=0
      do istate=1,nstat(1)
         do jstate=1,nstat(2)
            do j=1,(max_multipole+1)**2
               if (abs(properties(istate,jstate,j)).gt.threshold) then
                  write(lupropw,'(I1,I8,I3,I8,4I3,D20.12)') &
                    1, istate, csf_head_i % MGVN, jstate + nstat(1), csf_head_j % MGVN, ci_vec_i % nnuc + 1, &
                    pintegrals % lp(j), pintegrals % mp(j) * pintegrals % qp(j), properties(istate, jstate, j)
                  i=i+1

               endif

            end do

         end do

      end do

      end subroutine write_properties

      subroutine write_dyson_gcsf_orbitals(ludyson_gcsf,icontains_continuum_orbitals,csf_head_i,csf_head_j,&
     &                                icontinuum_orbital_table_i,icontinuum_orbital_table_j,     &
     &                                dyson_orbitals_cl2,dyson_orbitals_l2c)
      use cdenprop_defs
      implicit none

      integer :: ludyson_gcsf,icontains_continuum_orbitals
      type (CSFheader) :: csf_head_i,csf_head_j
      integer,        allocatable, dimension(:,:)   :: icontinuum_orbital_table_i, icontinuum_orbital_table_j
      real(kind=idp), intent(in), dimension(:,:,:) :: dyson_orbitals_cl2, dyson_orbitals_l2c

      if (icontains_continuum_orbitals .eq. 3) then
         if(csf_head_j%l2nocsf .ne. 0) then
            write(ludyson_gcsf) size(dyson_orbitals_cl2,1),size(dyson_orbitals_cl2,2),size(dyson_orbitals_cl2,3)
            write(ludyson_gcsf) size(icontinuum_orbital_table_j,1),size(icontinuum_orbital_table_j,2)
            write(ludyson_gcsf) size(csf_head_i%nob)
            write(ludyson_gcsf) size(csf_head_i%nob0)
            write(ludyson_gcsf) dyson_orbitals_cl2
            write(ludyson_gcsf) icontinuum_orbital_table_j
            write(ludyson_gcsf) csf_head_i%nob
            write(ludyson_gcsf) csf_head_i%nob0

         end if

      else if (icontains_continuum_orbitals .eq. 1) then
         write(ludyson_gcsf) size(dyson_orbitals_cl2,1),size(dyson_orbitals_cl2,2),size(dyson_orbitals_cl2,3)
         write(ludyson_gcsf) size(icontinuum_orbital_table_i,1),size(icontinuum_orbital_table_i,2)
         write(ludyson_gcsf) size(csf_head_i%nob)
         write(ludyson_gcsf) size(csf_head_i%nob0)
         write(ludyson_gcsf) dyson_orbitals_cl2
         write(ludyson_gcsf) icontinuum_orbital_table_i
         write(ludyson_gcsf) csf_head_i%nob
         write(ludyson_gcsf) csf_head_i%nob0

      else if (icontains_continuum_orbitals .eq. 2) then
         write(ludyson_gcsf) size(dyson_orbitals_l2c,1),size(dyson_orbitals_l2c,2),size(dyson_orbitals_l2c,3)
         write(ludyson_gcsf) size(icontinuum_orbital_table_j,1),size(icontinuum_orbital_table_j,2)
         write(ludyson_gcsf) size(csf_head_j%nob)
         write(ludyson_gcsf) size(csf_head_j%nob0)
         write(ludyson_gcsf) dyson_orbitals_l2c
         write(ludyson_gcsf) icontinuum_orbital_table_j
         write(ludyson_gcsf) csf_head_j%nob
         write(ludyson_gcsf) csf_head_j%nob0

      end if

      end subroutine write_dyson_gcsf_orbitals

      subroutine write_dyson_orbitals (ludyson, icontains_continuum_orbitals, iuse_bound, ibra_or_ket_bound, csf_head_i, &
                                       csf_head_j, icontinuum_orbital_table_i, icontinuum_orbital_table_j,     &
                                       dyson_orbitals_cl2, dyson_orbitals_l2c, ci_vector_i, nstat_i, TRANSI, ci_vector_j, nstat_j, &
                                       TRANSJ, dyson_orbitals, idyson_orbital_irrep, itarget_irrep, itarget_spin, &
                                       dyson_orbital_norms, pmat, mpi_i_write, iwrite)
      use cdenprop_defs
      implicit none

#if defined(mpi) && defined(scalapack)
      external dgsum2d
#endif

!     Arguments
      integer :: ludyson,icontains_continuum_orbitals, iuse_bound, ibra_or_ket_bound,nstat_i,nstat_j,iwrite
      logical :: mpi_i_write
      type (CSFheader) :: csf_head_i,csf_head_j
      integer,        allocatable, dimension(:,:)   :: icontinuum_orbital_table_i, icontinuum_orbital_table_j
      real(kind=idp), allocatable, dimension(:,:,:) :: dyson_orbitals_cl2, dyson_orbitals_l2c
      real(kind=idp), allocatable, dimension(:,:)   :: ci_vector_i, ci_vector_j
      real(kind=idp), allocatable, dimension(:,:,:) :: dyson_orbitals
      integer,        allocatable, dimension(:)     :: idyson_orbital_irrep, itarget_irrep, itarget_spin
      real(kind=idp), allocatable, dimension(:,:)   :: dyson_orbital_norms
      character(len=1), intent(in) :: TRANSI,TRANSJ
      type(CIvect),     intent(in) :: pmat          ! only needed for information about current BLACS context

!     Local variables
      integer :: i, j, icount, no_orbitals, no_target_states, no_bound_states, istart_l2_csf,      &
     &           itarget, ibound,max_orbital_relative, icsf, igcsf, iorb, iorb_absolute
      integer(blasint) :: one = 1, lngth
      integer, allocatable, dimension(:) :: is_target_orb, map_absolute_to_relative_index
      integer, allocatable, dimension(:,:) :: map_relative_to_absolute_index
      real(kind=idp) :: phase

!     First figure some stuff out
!     **********************************************************************************************
!

!     We must determine whether I or J (bra or ket) is the bound state and whether the
!     bound state includes continuum orbitals. If the bound state contains continuum orbitals
!     then so do the dyson orbitals.
!     ----------------------------------------------------------------------------------------------

      select case(icontains_continuum_orbitals)
      case(1)
!        J is the bound state and contains no continuum orbitals.
         ibra_or_ket_bound=2

      case(2)
!        I is the bound state and contains no continuum orbitals.
         ibra_or_ket_bound=1

      case(3)
!        Either I or J could be the bound state and contains continuum orbitals
         select case(iuse_bound)
         case(1)
!           I is the bound state and contains continuum orbitals (BOUND used)
            ibra_or_ket_bound=1
         case(2)
!           J is the bound state and contains continuum orbitals (BOUND used)
            ibra_or_ket_bound=2
         case(0)
!           BOUND not used I or J could be the bound state
            select case(ibra_or_ket_bound)
            case(0)
!              We don't know which is the bound state: exit subroutine
               write(iwrite,*) " Warning: ibra_or_ket_bound must be set to output dyson orbitals if iuse_bound=0 &
                          &and both bra and ket contain continuum orbitals."
               return
            case(1)
!              I is the bound state and contains continuum orbitals
            case(2)
!              J is the bound state and contains continuum orbitals
            end select

         end select

      end select

!     We set up the relative (within a symmetry) orbital index to absolute orbital index
!     mapping for the target orbitals. Continuum orbitals are already held in absolute index
!     form in icontinuum_orbital_table_i/j.
!     ----------------------------------------------------------------------------------------------

      select case(ibra_or_ket_bound)
      case(1)
!        J at least contains continuum orbitals and so the info we require
         max_orbital_relative=maxval(csf_head_j%nob0+csf_head_j%no_l2_virtuals)
         allocate(map_relative_to_absolute_index(max_orbital_relative,csf_head_j%nsym))

         map_relative_to_absolute_index=0

         icount=1
         do i =1,csf_head_j%nsym
!           Bound orbitals
            do j=1, csf_head_j%nob0(i)+csf_head_j%no_l2_virtuals(i)
               map_relative_to_absolute_index(j,i)=icount
               icount=icount+1

            end do
!           Skip continuum and virtual orbitals (if virtuals are not grouped with the continuum)
            do j=1, csf_head_j%nob(i)-csf_head_j%nob0(i)-csf_head_j%no_l2_virtuals(i)
               icount=icount+1

            end do

         end do

         no_orbitals=sum(csf_head_j%nob)
         no_target_states=sum(csf_head_j%numtgt)
         no_bound_states=nstat_i

         allocate(map_absolute_to_relative_index(no_orbitals))
         iorb_absolute=1
         do i =1,csf_head_j%nsym
            do j=1, csf_head_j%nob(i)

               map_absolute_to_relative_index(iorb_absolute)=j
               iorb_absolute=iorb_absolute+1

            end do

         end do

         if (icontains_continuum_orbitals .eq. 3) then
            istart_l2_csf=size(icontinuum_orbital_table_i,1) + 1

         else
            istart_l2_csf=1
         end if

      case(2)
!        I at least contains continuum orbitals and so the info we require.
         max_orbital_relative=maxval(csf_head_i%nob0+csf_head_i%no_l2_virtuals)
         allocate(map_relative_to_absolute_index(max_orbital_relative,csf_head_i%nsym))
         map_relative_to_absolute_index=0

         icount=1
         do i =1,csf_head_i%nsym
!           Bound orbitals.
            do j=1, csf_head_i%nob0(i)+csf_head_i%no_l2_virtuals(i)
               map_relative_to_absolute_index(j,i)=icount
               icount=icount+1

            end do

!           Skip continuum and virtual orbitals (if virtuals are not grouped with the continuum).
            do j=1, csf_head_i%nob(i)-csf_head_i%nob0(i)-csf_head_i%no_l2_virtuals(i)
               icount=icount+1

            end do

         end do

         no_orbitals=sum(csf_head_i%nob)
         no_target_states=sum(csf_head_i%numtgt)
         no_bound_states=nstat_j

         allocate(map_absolute_to_relative_index(no_orbitals))
         iorb_absolute=1
         do i =1,csf_head_i%nsym
            do j=1, csf_head_i%nob(i)

               map_absolute_to_relative_index(iorb_absolute)=j
               iorb_absolute=iorb_absolute+1

            end do

         end do

         if (icontains_continuum_orbitals .eq. 3) then
            istart_l2_csf=size(icontinuum_orbital_table_j,1) + 1

         else
            istart_l2_csf=1
         end if

      end select

      allocate(is_target_orb(no_orbitals))
      do iorb_absolute =1, no_orbitals
         if (ANY(map_relative_to_absolute_index .eq. iorb_absolute)) then
            is_target_orb(iorb_absolute)=1
         else
            is_target_orb(iorb_absolute)=0
         end if
      end do



!     Figure out dyson orbital symmetry
!     ---------------------------------
      allocate(idyson_orbital_irrep(no_target_states), itarget_irrep(no_target_states), itarget_spin(no_target_states) )
      idyson_orbital_irrep=0

      select case (icontains_continuum_orbitals)
      case(3) ! We can get the irrep from the table
         select case (ibra_or_ket_bound)
         case(1)

            do i =1, size(icontinuum_orbital_table_i,1)
               idyson_orbital_irrep(icontinuum_orbital_table_i(i,1))=icontinuum_orbital_table_i(i,4)
               itarget_irrep(icontinuum_orbital_table_i(i,1)) = icontinuum_orbital_table_i(i,3)-1
               itarget_spin(icontinuum_orbital_table_i(i,1)) = csf_head_i%itarget_symmetry_order(icontinuum_orbital_table_i(i,6),2)
            end do

         case(2)
            do i =1, size(icontinuum_orbital_table_j,1)
               idyson_orbital_irrep(icontinuum_orbital_table_j(i,1))=icontinuum_orbital_table_j(i,4)
               itarget_irrep(icontinuum_orbital_table_j(i,1)) = icontinuum_orbital_table_j(i,3)-1
               itarget_spin(icontinuum_orbital_table_j(i,1)) = csf_head_j%itarget_symmetry_order(icontinuum_orbital_table_j(i,6),2)

            end do

         end select
      case default !We need to figure out the irrep from the target irrep
         select case (ibra_or_ket_bound)
         case(1)

            do i =1, size(icontinuum_orbital_table_j,1)

               idyson_orbital_irrep(icontinuum_orbital_table_j(i,1))=IPD2H( csf_head_i%mgvn+1, icontinuum_orbital_table_j(i,3) ) - 1
               itarget_irrep(icontinuum_orbital_table_j(i,1)) = icontinuum_orbital_table_j(i,3)-1
               itarget_spin(icontinuum_orbital_table_j(i,1)) = csf_head_j%itarget_symmetry_order(icontinuum_orbital_table_j(i,6),2)

            end do

         case(2)
            do i =1, size(icontinuum_orbital_table_i,1)
               idyson_orbital_irrep(icontinuum_orbital_table_i(i,1))=IPD2H( csf_head_j%mgvn+1, icontinuum_orbital_table_i(i,3) ) - 1
               itarget_irrep(icontinuum_orbital_table_i(i,1)) = icontinuum_orbital_table_i(i,3)-1
               itarget_spin(icontinuum_orbital_table_i(i,1)) = csf_head_i%itarget_symmetry_order(icontinuum_orbital_table_i(i,6),2)

            end do

         end select

      end select

!     Now multiply the GCSF dyson orbitals by the ci vectors to produce the dyson orbitals
!     **********************************************************************************************

      allocate( dyson_orbitals(no_orbitals, no_target_states, no_bound_states) )
      dyson_orbitals=0

      select case (ibra_or_ket_bound)
      case(1)
         do ibound=1,no_bound_states
            do itarget=1, no_target_states

!              The continuum CSFs (just continuum orbitals)
               if (icontains_continuum_orbitals .eq. 3) then
                  do igcsf=1, istart_l2_csf-1

                     if (icontinuum_orbital_table_i(igcsf,4) == idyson_orbital_irrep(itarget) .and. &
                         icontinuum_orbital_table_i(igcsf,1) == itarget) then
                        iorb_absolute= icontinuum_orbital_table_i(igcsf,5)
                        if (TRANSI .eq. 'T') then
                           dyson_orbitals(iorb_absolute, itarget, ibound) = &
                           dyson_orbitals(iorb_absolute, itarget, ibound) + ci_vector_i(igcsf,ibound)
                        else
                           dyson_orbitals(iorb_absolute, itarget, ibound) = &
                           dyson_orbitals(iorb_absolute, itarget, ibound) + ci_vector_i(ibound,igcsf)
                        endif
                     end if

                  end do

               end if
!              The L2 csfs (just target orbitals)
               if (allocated(dyson_orbitals_l2c)) then
                  do icsf=1, size(dyson_orbitals_l2c,2)
                     call pmat % global_to_local(icsf, 1, igcsf, j) ! CSFs distributed vertically within this BLACS context
                     do iorb=1,max_orbital_relative
            
                        iorb_absolute= map_relative_to_absolute_index(iorb, idyson_orbital_irrep(itarget)+1)

                        if (iorb_absolute .ne. 0) then
                           if (TRANSI .eq. 'T') then
                              dyson_orbitals(iorb_absolute, itarget, ibound) = &
                              dyson_orbitals(iorb_absolute, itarget, ibound) + &
                              dyson_orbitals_l2c(iorb, icsf, itarget) * ci_vector_i(igcsf,ibound)
                           else
                              dyson_orbitals(iorb_absolute, itarget, ibound) = &
                              dyson_orbitals(iorb_absolute, itarget, ibound) + &
                              dyson_orbitals_l2c(iorb, icsf, itarget) * ci_vector_i(ibound,igcsf)
                           endif
                        end if

                     end do
                  end do
#if defined(usempi) && defined(scalapack)
                   ! now combine the partial results obtained by individual members of this BLACS context
                   if (pmat % CV_is_scalapack) then
                      lngth = size(dyson_orbitals)
                      call dgsum2d(pmat % blacs_context, 'R', ' ', lngth, one, dyson_orbitals, lngth, -one, -one)
                   end if
#endif
               end if

            end do
         end do

      case(2)
         do ibound=1,no_bound_states
            do itarget=1, no_target_states

!              The continuum CSFs (just continuum orbitals)
               if (icontains_continuum_orbitals .eq. 3) then
                  do igcsf=1, istart_l2_csf-1

                     if (icontinuum_orbital_table_j(igcsf,4) == idyson_orbital_irrep(itarget) .and. &
                         icontinuum_orbital_table_j(igcsf,1) == itarget) then
                        iorb_absolute= icontinuum_orbital_table_j(igcsf,5)
                        if (TRANSJ .eq. 'T') then
                           dyson_orbitals(iorb_absolute, itarget, ibound) = &
                           dyson_orbitals(iorb_absolute, itarget, ibound) + ci_vector_j(ibound, igcsf)
                        else
                           dyson_orbitals(iorb_absolute, itarget, ibound) = &
                           dyson_orbitals(iorb_absolute, itarget, ibound) + ci_vector_j(igcsf, ibound)
                        endif
                     end if
                  end do

               end if
!              The L2 csfs (just target orbitals)
               if (allocated(dyson_orbitals_cl2)) then
                  do icsf=1, size(dyson_orbitals_cl2,2)
                     call pmat % local_to_global(1, icsf, j, igcsf)  ! CSFs distributed horizontaly within current BLACS context
                     do iorb=1,max_orbital_relative
            
                         iorb_absolute= map_relative_to_absolute_index(iorb, idyson_orbital_irrep(itarget)+1)

                         if (iorb_absolute .ne. 0) then
                            if (TRANSJ .eq. 'T') then
                               dyson_orbitals(iorb_absolute, itarget, ibound) = &
                               dyson_orbitals(iorb_absolute, itarget, ibound) + &
                               dyson_orbitals_cl2(iorb, icsf, itarget) * ci_vector_j(ibound, igcsf)
                            else
                               dyson_orbitals(iorb_absolute, itarget, ibound) = &
                               dyson_orbitals(iorb_absolute, itarget, ibound) + &
                               dyson_orbitals_cl2(iorb, icsf, itarget) * ci_vector_j(igcsf, ibound)
                            endif
                         end if

                      end do
                   end do
#if defined(usempi) && defined(scalapack)
                   ! now combine the partial results obtained by individual members of this BLACS context
                   if (pmat % CV_is_scalapack) then
                      lngth = size(dyson_orbitals)
                      call dgsum2d(pmat % blacs_context, 'C', ' ', lngth, one, dyson_orbitals, lngth, -one, -one)
                   end if
#endif
               end if

            end do
         end do

      end select

              allocate(dyson_orbital_norms( no_target_states, no_bound_states))

        !     Calculate the norm of each dyson orbital.
              do ibound=1,no_bound_states
                 do itarget=1, no_target_states
                   dyson_orbital_norms(itarget, ibound) = &
                    sqrt(dot_product(dyson_orbitals(:, itarget, ibound), dyson_orbitals(:, itarget, ibound)))
                 end do
              end do

        !     Finally write the dyson orbitals to file
   !     ZM: only master writes
   !     Note: to get dyson orbital norms and coefficients for specific target spin states
   !     the Dyson orbital norm and coefficient that are written to file must be multiplied by
   !     the absolute value of the appropriate Clebsch-Gordan coefficient.
   !     E.g. for doublet ionic states with singlet neutral state the factor is sqrt(1/2). For triplet
   !     ionic states with doublet neutral state the factor is sqrt(1/3) for ionic M=0  and sqrt(2/3)
   !     for ionic M_S= +/- 1.
   !     In general the factor is |< S_ion M_{S_ion} 1/2, m_s |S_neutral M_{S_neutral}>|
        !     ----------------------------------------------------------------------------------------------
              if (mpi_i_write) then
                 write(ludyson, '(" Dyson orbitals - Linear combination of the scattering run swedmos orbitals")')
                 write(ludyson, '(" bound state, target index, target spin, target IR, dyson orbital IR, is target orbital, &
                                   &relative orbital index, absolute orbital index, orbital coefficient")')
              endif

              do ibound=1,no_bound_states
                 do itarget=1, no_target_states
                    phase = csf_head_i%itarget_overall_phase(itarget)*csf_head_j%itarget_overall_phase(itarget)
                    if (mpi_i_write) write(ludyson, '(" Dyson orbital norm = ",d20.5)')  dyson_orbital_norms(itarget, ibound)
                    do iorb_absolute=1, no_orbitals

                       !Multiply in the relative phase between the
                       !target states to the continuum orbitals:
                       if (is_target_orb(iorb_absolute) .eq. 0) then
                          dyson_orbitals(iorb_absolute, itarget, ibound) = phase*dyson_orbitals(iorb_absolute, itarget, ibound)
                       endif

                       if (mpi_i_write) write(ludyson, '(8i5, d20.5)') ibound, itarget, itarget_spin(itarget), &
             &                                        itarget_irrep(itarget), idyson_orbital_irrep(itarget),&
     &                                        is_target_orb(iorb_absolute),map_absolute_to_relative_index(iorb_absolute),&
     &                                        iorb_absolute, dyson_orbitals(iorb_absolute, itarget, ibound)

            end do
            if (mpi_i_write) write(ludyson,*) ""
         end do
         if (mpi_i_write) write(ludyson,*) ""
      end do

      end subroutine write_dyson_orbitals

!
!     Write target phases to file for reading by outer region codes.
!     (note: in a future version this should go in the CI vectors file
!      the best way to do this would be to figure out the phase in scatci
!      instead of waiting until cdenprop)
!
      subroutine write_target_phases(ntgsym,numtgt,itarget_overall_phase)
      implicit none
!     Arguments
      integer :: ntgsym
      integer, dimension(:) :: numtgt,itarget_overall_phase
!     Local
      integer :: no_target_states, itarg, i, j

      no_target_states=sum(numtgt(1:ntgsym))

      OPEN(UNIT=55555, FILE='target.phases.data',status='replace')
      write(55555,'(i10)') no_target_states
      itarg=0
      do i=1, ntgsym
         do j=1,numtgt(i)
            itarg=itarg+1
            write(55555, '(3i10)') i, j, itarget_overall_phase(itarg)

         end do
      end do
      close(55555)

      end subroutine write_target_phases

      end module cdenprop_io

