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
!*==numcbas.spg  processed by SPAG 6.56Rc at 11:55 on  9 Mar 2010
      PROGRAM NUMCBAS
C
      USE PRECISN, ONLY : WP                        
      USE GLOBAL_UTILS, ONLY : PRINT_UKRMOL_HEADER
      USE consts, ONLY : ZERO => xzero
      USE NUMCBAS_DATA, ONLY : HRX, IRX, IRA, NIX, MAXRX, NFTA, RTOL
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL(KIND=wp) :: BTOL=0.2_wp, CHARGE=0.0_wp, ECMAX=10.0_wp, 
     &                 RLIM=10.0_wp, RMAT, TINY=1.0E-11_wp, TRANGE
      INTEGER :: I, IER, IRXO, LUNUMB=13, LVAL=0
      INTEGER, DIMENSION(3) :: IBUG=0
      CHARACTER(LEN=120) :: TITLE
C
C*** End of declarations rewritten by SPAG
C
C     MAIN DRIVING ROUTINE
C
C
      NAMELIST /INPUT / TITLE, LUNUMB, NIX, IRX, HRX, lval, IBUG, BTOL, 
     &   TINY, ECMAX, RLIM, CHARGE
C
      CALL PRINT_UKRMOL_HEADER(NFTA)
      WRITE(NFTA,1000)
C
      DO I=1, LEN(TITLE)
         TITLE(I:I)=' '
      END DO
C
C --- Set default values for NIX, HRX and IRX
C
      NIX=3
      HRX=(/ 1.E-02_wp, 2.E-02_wp, 2.605E-02_wp, (ZERO,I=4,MAXRX) /)
      IRX=(/ 30, 120, 500, (0,I=4,MAXRX) /)
C
C --- Read input data via NAMELIST
C
      READ(5,INPUT)
      RMAT=RLIM
      WRITE(NFTA,50)TITLE, LUNUMB
C
C --- Set up integration mesh and check its consistency
C --- Check that we cover the whole range [0;RLIM]. The number of
C     points in the last subrange is changed to ensure that we do
C
 
C --- Check for inappropriate or incomplete input data concerning
C     NIX,IRX,HRX
      IF(NIX>MAXRX)THEN
         WRITE(NFTA,80)NIX,MAXRX
         STOP
      END IF
      DO I=1, NIX
         IF(IRX(I)==0 .OR. HRX(I)==ZERO)THEN
            WRITE(NFTA,90)I
            STOP
         END IF
      END DO
 
      WRITE(NFTA,51)NIX
      IER=0
      IRXO=0
      TRANGE=ZERO
      DO I=1, NIX
         TRANGE=TRANGE+(IRX(I)-IRXO)*HRX(I)
         IF(I.EQ.NIX .AND. ABS(TRANGE-RMAT).GE.RTOL)THEN
            IRX(I)=INT((RMAT-TRANGE)/HRX(I))+IRX(I)
            IF(MOD(IRX(I),2).NE.0)IRX(I)=IRX(I)+1
            WRITE(NFTA,52)IRX(I)
         END IF
         IRXO=IRX(I)
         WRITE(NFTA,53)I, IRX(I), HRX(I)
C
C --  Check that the number of points in each subrange is even
C
         IF(MOD(IRX(I),2).NE.0)IER=1
      END DO
C
      IF(IER.NE.0)THEN
         WRITE(NFTA,70)
         STOP
      END IF
      IRA=IRX(NIX)
C
      WRITE(NFTA,68)lval, charge, RMAT
      WRITE(NFTA,69)ECMAX, BTOL, TINY
C
C --- CALCULATE THE NUMERICAL BASIS
C
      CALL BASIS(TITLE,ecmax,charge,lval,tiny,btol,lunumb,ibug)
C
 50   FORMAT(//4X,70('*')/4X,'*',68X,'*'/4X,'*',2X,A64,2X,'*'/4X,'*',
     &       68X,'*'/4X,70('*')///6X,'LUNUMB  =',I3,5X,
     &       'Output file for the basis')
 51   FORMAT(//11X,'INTEGRATION MESH INPUT DATA'//6X,'NIX   =',I3,5X,
     &       'No. of integration regions with different step-',
     &       'sizes'//10X,'I',7X,'IRX',15X,'HRX')
 52   FORMAT(//5X,'Number of points (IRX) in last subrange has',
     &       ' been changed to',5X,I5)
 53   FORMAT(1x,2(5X,I5),5X,D20.10,5X,I5)
 68   FORMAT(//11X,'NUMERICAL BASIS CALCULATION INPUT DATA'//6X,
     &       'LVAL   =',I8,5X,'Angular Momentum'/6x,'CHARGE =',F8.1,5X,
     &       'Effective charge'/6x,'RMAT   =',F8.1,5X,
     &       'R-matrix boundary radius')
 69   FORMAT(//11X,'SEARCHING PROCEDURE PARAMETERS'//6X,'ECMAX  =',F8.2,
     &       5X,'Maximum energy for the ','eigensolutions'/6X,
     &       'BTOL   =',D8.2,5X,'Iteration starting tolerance'/6X,
     &       'TINY   =',D8.2,5X,'Eigensolution convergence parameter')
 70   FORMAT(///11X,'***  INPUT DATA ERROR - JOB ABORTED  ***'//6X,
     &       'ERROR IN THE INTEGRATION MESH SPECIFICATION'/6X,
     &       'THE NUMBER OF STEPS IN ONE OR MORE REGIONS IS NOT ',
     &       'EXACTLY DIVISIBLE BY 2')
 1000 FORMAT(//11x,' Program NUMCBAS',//)
 80   FORMAT(///11X,'***  INPUT DATA ERROR - JOB ABORTED  ***'//6X,
     &       'ERROR IN THE INTEGRATION MESH SPECIFICATION'/6X,'NIX =',
     &       I3,'; MUST BE <= ',I3)
 90   FORMAT(///11X,'***  INPUT DATA ERROR - JOB ABORTED  ***'//6X,
     &       'ERROR IN IRX AND/OR HRX STARTING AT'/6X,'COMPONENT ',I2,
     &       ' - SOME VALUES ZERO')
C
      STOP
C
      END PROGRAM NUMCBAS
!*==basis.spg  processed by SPAG 6.56Rc at 11:55 on  9 Mar 2010
C
      SUBROUTINE BASIS(TITLE,ecmax,charge,lval,tiny,btol,lunumb,ibug)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, HALF=>XHALF, ONE=>XONE, TWO=>XTWO,
     &                   THREE=>XTHREE, IZERO
      USE NUMCBAS_DATA, ONLY : MAXORB, HRX, IRX, IRA, NIX, NFTA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: BTOL, CHARGE, ECMAX, TINY
      INTEGER :: LUNUMB, LVAL
      CHARACTER(LEN=*) :: TITLE
      INTEGER, DIMENSION(3) :: IBUG
      INTENT (IN) CHARGE, IBUG
C
C Local variables
C
      REAL(KIND=wp) :: A, A1, A2, FUN1, FUN2, H, RMAT, RMAT1, RR, X, Y
      REAL(KIND=wp), DIMENSION(maxorb) :: AK
      REAL(KIND=wp), DIMENSION(IRA+1,maxorb) :: DORB, ORB
      INTEGER :: I, IR, IR1, IR2, IRA1, J, JRR, L, NCO, NFTOT, NIX1
      REAL(KIND=wp), DIMENSION(IRA) :: R
      CHARACTER(LEN=120) :: TITLE1
      REAL(KIND=wp), DIMENSION(2*ira) :: VV
C
C*** End of declarations rewritten by SPAG
C
C     THIS SUBROUTINE PERFORMS THE CALCULATION OF THE NUMERICAL BASIS
C
 50   FORMAT(///21X,16('-')//21X,'SUBROUTINE BASIS'//21X,16('-')//11X,
     &       'NUMERICAL FUNCTIONS CALCULATION INPUT DATA')
 54   FORMAT(//6X,'POTENTIAL ON HALF MESH'//10x,'R',12x,'VV')
 55   FORMAT(i5,f10.5,1x,D16.8)
 65   FORMAT(7x,i5,5x,I5,5X,D20.10)
 57   FORMAT(///11X,'NUMERICAL BASIS INTEGRATION MESH'//11X,'NIX   =',
     &       I5,5X,'NUMBER OF INTEGRATION REGIONS'//11X,'I',7X,'IRX',
     &       15X,'HRX')
 58   FORMAT(//6X,'IRA   =',I5,10X,'TOTAL NUMBER OF INTEGRATION ',
     &       'STEPS'/6X,'RMAT  =',F9.4,6X,'R-MATRIX BOUNDARY RADIUS'/6X,
     &       'IRA1  =',I5,10X,'NUMBER OF OUTWARD INTEGRATION STEPS'/6X,
     &       'RMAT1 =',F9.4,6X,'MATCHING RADIUS')
 60   FORMAT(i5,f10.4,D16.8)
 61   FORMAT(6X,I5,5x,D16.8,4X,D16.8)
 71   FORMAT(///21X,'CHECK THE CONTENTS OF FILE LUNUMB'//1X,A//11X,
     &       'NTOT =',I5,9X,'L =',I3)
 74   FORMAT(//5X,'AMPLITUDES AT   R = 0.0',12X,'R =',F6.2/)
 75   FORMAT(//5X,28('-')/12X,'FINAL BASIS'/5X,28('-'))
 77   FORMAT(//11X,'FUNCTION NO.',I3/)
 81   FORMAT(//' POLE POSITIONS'/(5D15.6))
C
      IF(ibug(1).NE.0)WRITE(NFTA,50)
C
C     Set up MESH PARAMETERS PRIOR TO SOLVING DIFFERENTIAL EQUATIONS
C
      A1=TWO
      A2=TWO
      X=ZERO
      IR2=0
C
      DO I=1, NIX
         IR1=IR2+1
         IR2=IRX(I)
         H=HRX(I)
C jmc Y set but never used???        Y=H/THREE
C
         DO IR=IR1, IR2
            X=X+H
            A1=A1+A2
            A2=-A2
            R(IR)=X
         END DO
C
C jmc Y set but never used???         IF(I.LT.NIX)Y=(H+TWO*HRX(I+1))/THREE
      END DO
C
      RMAT=R(IRA)
      IF(NIX.LE.1)GO TO 7
      NIX1=NIX-1
      IF(IRA-IRX(NIX1).LE.100)GO TO 8
 7    NIX1=NIX
      NIX=NIX+1
      IRX(NIX)=IRA
      IRX(NIX1)=IRA-100
      HRX(NIX)=HRX(NIX1)
 8    CONTINUE
      IRA1=IRX(NIX1)
      RMAT1=R(IRA1)
      IF(ibug(1).NE.0)THEN
         WRITE(NFTA,57) NIX
         WRITE(NFTA,65) (I,IRX(I),HRX(I),I=1,NIX)
         WRITE(NFTA,58) IRA, RMAT, IRA1, RMAT1
      END IF
C
C --- Set up potential
C
      IF(IBUG(1).GT.0)WRITE(NFTA,54)
      A2=REAL(lval*(lval+1),kind=wp)
      IR2=0
      X=ZERO
      DO I=1, NIX
         IR1=IR2+1
         IR2=2*IRX(I)
         H=half*HRX(I)
         DO IR=IR1, IR2
            X=X+H
            VV(IR)=A2/(X*X)-two*charge/X
            IF(IBUG(1).NE.0)WRITE(NFTA,55) ir, x, VV(IR)
         END DO
      END DO
C
C     CALL  SUBROUTINE SEARCH TO CALCULATE NUMERICAL FUNCTIONS
C
      CALL SEARCH(NCO,AK,ORB,DORB,rmat,ecmax,vv,lval,tiny,btol,maxorb)
C
C     WRITE THE FUNCTIONS PRODUCED BY SEARCH ON LUNUMB
C
      CALL WRHEAD(NCO,lval,ira,r,TITLE,lunumb)
      IF(IBUG(2).NE.0)WRITE(NFTA,75)
      DO J=1, nco
C     NOTE RENUMBERING SO THAT ORB(2) IS NOW ORB(1)
         RR=ONE/R(1)
         FUN1=ORB(1,j)*RR
         RR=ONE/R(2)
         FUN2=RR*ORB(2,j)
         IR2=0
         JRR=0
         DO I=1, NIX
            IR1=IR2+1
            IR2=IRX(I)
            DO IR=IR1, IR2
               JRR=JRR+1
               X=ONE/R(IR)
               ORB(JRR,j)=X*ORB(IR,j)
            END DO
         END DO
C
C     IF L = 0 USE A TWO-POINTS LAGRANGE INTERPOLATION TO CALCULATE
C     THE VALUE AT THE ORIGIN - SET IT = 0.0 OTHERWISE
C
         IF(lval.EQ.0)THEN
            X=TWO*FUN1-FUN2
         ELSE
            X=ZERO
         END IF
C
         WRITE(LUNUMB)X, (ORB(IR,j),IR=1,IRA)
         IF(IBUG(2).NE.0)THEN
            WRITE(NFTA,77) J
            WRITE(NFTA,60) izero, zero, X, (ir,r(ir),ORB(IR,j),IR=1,IRA)
         END IF
      END DO
C
C     ADD POLE POSITIONS at end of LUNUMB
C
      WRITE(LUNUMB)(AK(I),I=1,NCO)
C
C     CHECK THE CONTENTS OF LUNUMB
C
      IF(ibug(3).NE.0)THEN
         REWIND LUNUMB
         READ(LUNUMB)TITLE1
         READ(LUNUMB)NFTOT, L, ira
         WRITE(NFTA,71) TITLE1, NFTOT, L
         READ(LUNUMB)A, (R(IR),IR=1,IRA)
         WRITE(NFTA,74) RMAT
         DO I=1, NFTOT
            READ(LUNUMB)A, (R(IR),IR=1,IRA)
            WRITE(NFTA,61) I, A, R(IRA)
         END DO
         READ(LUNUMB)(AK(I),I=1,NFTOT)
         WRITE(NFTA,81) (AK(I),I=1,NFTOT)
      END IF
C
      RETURN
C
      END SUBROUTINE BASIS
!*==search.spg  processed by SPAG 6.56Rc at 11:55 on  9 Mar 2010
C
      SUBROUTINE SEARCH(NCO,AK,ORB,DORB,rmat,ecmax,vv,lval,tiny,btol,
     &                  maxorb)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO
      USE NUMCBAS_DATA, ONLY : IRA, NFTA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: BTOL, ECMAX, RMAT, TINY
      INTEGER :: LVAL, MAXORB, NCO
      REAL(KIND=wp), DIMENSION(maxorb) :: AK
      REAL(KIND=wp), DIMENSION(IRA+1,maxorb) :: DORB, ORB
      REAL(KIND=wp), DIMENSION(2*IRA) :: VV
      INTENT (IN) ECMAX, MAXORB, TINY
      INTENT (OUT) NCO
      INTENT (INOUT) AK
C
C Local variables
C
      REAL(KIND=wp) :: BVALUE, E, EIGEN, ETA, ETRIAL
      INTEGER :: IR, J, NFIND, NODES, NODMAX
      REAL(KIND=wp), DIMENSION(ira+1) :: ORBL
C
C*** End of declarations rewritten by SPAG
C
C     Searches for eigen-solutions to the differential equations
C
 61   FORMAT(//21X,'SUMMARY TABLE'//6X,'Partial wave L =',I3,3X,
     &       'No. of eigensolutions  =',I4//18X,'Nodes',4X,
     &       'Eigenenergy (Ryd.)')
 62   FORMAT(1X,2(5X,I5),5X,D16.8)
 65   FORMAT(6X,'NODMAX =',I8,5X,'Last eigensolution within ',
     &       'the energy range')
 66   FORMAT(/' *** ERROR *** NODMAX < 0')
 
C
C     CALCULATES THE MAXIMUM NUMBER OF NODES OF THE
C     EIGENSOLUTIONS WITHIN THE SPECIFIED ENERGY RANGE
C
      E=ECMAX
      CALL BASFUN(lval,NODMAX,E,BVALUE,VV,ORB,DORB,rmat)
      IF(BVALUE.GT.ZERO)THEN
         NODMAX=NODMAX-1
         WRITE(NFTA,65) NODMAX
      END IF
 
C --- JMC apply test on nodmax; stop if <0
      IF(NODMAX<0)THEN
         WRITE(NFTA,66)
         STOP
      END IF
C
C --- Find eigen-solutions with nodes = 0 to nodmax
      ETA=TINY
      ETRIAL=ETA
C
      DO NODES=0, NODMAX
C
         CALL FINDER(lval,NODES,BTOL,ETA,ETRIAL,EIGEN,VV,ORBL,DORB,rmat)
         NFIND=NODES+1
C
C --- Store the solution
         DO IR=2, IRA+1
            orb(ir-1,nfind)=orbl(ir)
         END DO
C
C --- STORE THE EIGEN-ENERGY
         AK(NFIND)=EIGEN
C
      END DO
C
c --- Print summary
      WRITE(NFTA,61) lval, NFIND
      DO J=1, NFIND
         WRITE(NFTA,62) J, j-1, AK(J)
      END DO
C
      NCO=NFIND
C
      RETURN
      END SUBROUTINE SEARCH
!*==basfun.spg  processed by SPAG 6.56Rc at 11:55 on  9 Mar 2010
C
      SUBROUTINE BASFUN(LC,NODES,WINIT,BVALUE,POVALU,ORB,DORB,ra)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, HALF=>XHALF, ONE=>XONE, 
     &                   ONEPT5=>ONEP5, TWO=>XTWO, THREE=>XTHREE, 
     &                   FOUR=>XFOUR, TWELVE=>XTWELVE
      USE NUMCBAS_DATA, ONLY : TINYY, HRX, IRX, IRA, NIX, NFTA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: BVALUE, RA, WINIT
      INTEGER :: LC, NODES
      REAL(KIND=wp), DIMENSION(ira+1) :: DORB, ORB
      REAL(KIND=wp), DIMENSION(2*ira) :: POVALU
      INTENT (IN) LC, RA, WINIT
      INTENT (OUT) BVALUE
      INTENT (INOUT) DORB, NODES, ORB
C
C Local variables
C
      REAL(KIND=wp) :: B2, B3, DU, DUH, DYIN, DYOUT, EPSILON, FM, 
     &                FR, FRHVAL, FRM, H, H1, H12, HH12, HIMT, HS, 
     &                ORB1, ORB2, PV1, PV2, S, 
     &                TLC, TWOZ, U, WR, X, X1, X2, 
     &                YIN, YOUT, YR
      REAL(KIND=wp), DIMENSION(2) :: BNDRY, DBNDRY
      REAL(KIND=wp), DIMENSION(ira+1) :: DUS, US
      INTEGER :: I, I1, I2, I9, IMATCH, IMT, INT, IR, ITST, J, JR, K, 
     &           KM, LP1, LP2, LSWT, MMM, NIZ, NUTTY
C
C*** End of declarations rewritten by SPAG
C
C***********************************************************************
C
C                      INPUT
C
C***********************************************************************
C
C  1. THE USER MUST PROVIDE THE FOLLOWING INPUT DATA...
C
C     LC........ THE  ANGULAR MOMENTUM VALUE
C     RA........ THE BOUNDARY RADIUS
C     WINIT..... THE INITIAL ENERGY
C     NIX....... THE NUMBER OF CHANGES OF INTEGRATION STEP OR THE NUMBER
C                OF SUBINTERVALS INTO WHICH THE INTERVAL X=0 TO RA IS
C                DIVIDED
C     HRX(I),I=1,NIX.... THE INTEGRATION STEP IN EACH SUBINTERVAL
C     IRX(I),I=1,NIX.... THE TOTAL NUMBER OF INTEGRATION STEPS TO THE
C                        END OF THE ITH SUBINTERVAL
C
C
C  2. THE USER MUST PROVIDE THE POTENTIAL FUNCTION AND STORE IT IN THE
C      ARRAY POVALU
C
C     THE ODD ELEMENTS POVALU(2N-1),N=1,IRX(NIX), SHOULD
C     CONTAIN THE FUNCTION VALUES AT THE HALF-MESH POINTS.
C
C     THE EVEN ELEMENTS POVALU(2N),N=1,IRX(NIX), SHOULD
C     CONTAIN THE VALUES AT THE MESH POINTS.
C
C***********************************************************************
C
C                       OUTPUT RESULTS
C
C***********************************************************************
C
C     ORB(K) , K=1,IRX(NIX)+1 CONTAINS THE FINAL SOLUTION
C     NODES....  IS THE NUMBER OF NODES IN THE FUNCTION
C     BVALUE IS THE  LOGARITHMIC DERIVATIVE VALUE AT R=RA
C
C***********************************************************************
C
 1004 FORMAT(/' ERROR - IRX(',I1,')  IS AN ODD NUMBER IN BASFUN'/)
 1012 FORMAT(/' ERROR - NIX,RA,IRX AND HRX ARE INCOMPATIBLE IN BASFUN'/
     &       (I5,F10.6))
 1016 FORMAT(/' *** ERROR *** DIVISION BY ZERO VIA ORB(I9)'/
     &       ' ANGULAR MOMENTUM TOO LARGE? LVAL =',I5)
 1017 FORMAT(/' *** ERROR *** DIVISION BY ZERO VIA X1; STOPPING'/)
C
C      CHECK COMPATIBILITY OF NIX,HRX,IRX AND RA
C
      S=ZERO
      J=0
      DO I=1, NIX
         S=S+HRX(I)*REAL(IRX(I)-J,kind=wp)
         J=IRX(I)
      END DO
      IF(ABS(S-RA).GT.TINYY)THEN
         WRITE(NFTA,1012) NIX, RA, (IRX(K),HRX(K),K=1,NIX)
         RETURN
      END IF
C
C      CHECK THAT IRX(I) ARE EVEN INTEGERS
C
      DO I=1, NIX
         IF(MOD(IRX(I),2).NE.0)THEN
            WRITE(NFTA,1004) I
            STOP
         END IF
      END DO
C
C      EVALUATE AND INITIALIZE SOME COMMONLY USED PARAMETERS
C
      FM=ZERO
C *** PATCH
      IMATCH=IRX(NIX)-30
      ITST=1
C *** END PATCH
      I9=IRX(NIX)+1
      DO I=1, I9
         ORB(I)=ZERO
         DORB(I)=ZERO
      END DO
      TLC=REAL(LC*(LC+1),kind=wp)
      HIMT=HRX(1)
      WR=WINIT
      PV2=TLC/(HIMT*HIMT)-POVALU(2)
      TWOZ=ZERO
      IF(PV2.GE.TINYY)THEN
         PV1=FOUR*TLC/(HIMT*HIMT)-POVALU(1)
         NIZ=PV1/PV2-HALF  ! jmc is this real -> integer implicit conversion relied upon for the following?
         IF(NIZ.GE.0)THEN
            TWOZ=PV1*HIMT*HALF
         END IF
      END IF
C
C      INITIALIZATION OF FUNCTION AT X=0.0 WHICH IS USED AS THE FIRST
C      POINT IN THE SIMPSONS RULE NORMALIZATION
C
      X=ZERO
      U=ZERO
      DU=ZERO
      US(1)=ZERO
      DUS(1)=ZERO
C
C      EVALUATE THE FUNCTION AND  DERIVATIVE AT HINT AND STORE THE
C      FUNCTION
C
      LP1=LC+1
      LP2=LC+2
      X=HIMT
      H=HIMT
      U=X**LP1*(ONE-TWOZ*X/(TWO*REAL(LP1,kind=wp)))
      DU=HALF*X**LP1*(REAL(LP1,kind=wp)-TWOZ*REAL(LP2,kind=wp)*X/
     &   (TWO*REAL(LP1,kind=wp)))
      US(2)=U
      DUS(2)=DU*TWO/H
C
C      EVALUATE FR AT HINT
C
      MMM=1
C
      HH12=H*H/TWELVE
      FM=(POVALU(2)-WR)*HH12
      FR=FM*U
C
C      SET UP FRM AT HINT/2
C
      X=HALF*X
      FRM=U
      U=X**LP1*(ONE-TWOZ*X/(TWO*REAL(LP1,kind=wp)))
      ! jmc unnecessary bit of code here... mistake?
C
      HH12=H*H/THREE
      U=FRM
      FRM=(POVALU(1)-WR)*HH12*U
C
      KM=2
C
C      INTEGRATE OUT TO THE MATCHING POINT
C
      LSWT=1
      I1=2
      DO INT=1, NIX
         H1=H
         H=HRX(INT)
         HS=H*H/THREE
         I2=IRX(INT)
         IF(IMATCH.LT.I2)I2=IMATCH
         IF(INT.GT.1)THEN
            LSWT=2
            I1=IRX(INT-1)+1
            H12=H/H1
         END IF
C
C      INTEGRATE OVER A RANGE OF EQUAL INTERVALS
C
         DO IR=I1, I2
            JR=IR+1
            FRHVAL=POVALU(JR)-WR
            IF(FRHVAL*FM.GE.ZERO)THEN
               FM=FRHVAL ! jmc is this a mistake?  what is the point of resetting fm here???
            ELSE IF(ITST.EQ.0)THEN ! jmc which it never will... 
               IMATCH=JR
               ITST=1
            END IF
            CALL DEVGL(DU,FR,FRM,HS,U,YR,POVALU,WR,lswt,mmm,h12,km)
C
C      STORE THE FUNCTIONS AT EACH INTEGRATION
C
            US(JR)=U
            DUS(JR)=DU*TWO/H
C
            IF(IR.EQ.IMATCH)GO TO 42
         END DO
      END DO
C
C      STORE THE FUNCTIONS AND DERIVATIVES AT THE MATCHING POINT FOR THE
C      OUTWARD INTEGRATION
C
 42   YOUT=U
      DYOUT=DU*TWO/H
C
C      INITIALIZE ARRAYS FOR DEVOGELAERE INTEGRATION INWARDS
C
      MMM=-1
      ITST=1
      BNDRY(1)=ONE
      DBNDRY(1)=ZERO
      BNDRY(2)=ZERO
      DBNDRY(2)=ONE
C
      DO NUTTY=1, 2
C
C      EVALUATE THE FUNCTION AND DERIVATIVES AT RA
C
         H=HRX(NIX)
         HS=H*H/THREE
         U=BNDRY(NUTTY)
         DU=-DBNDRY(NUTTY)*H*HALF
C
C      EVALUATE FR AT RA
C
         KM=2*IRX(NIX)
         FR=(POVALU(KM)-WR)*HS/FOUR*U
C
C      EVALUATE FRM AT RA+H/2
C
         FRM=(U-DU+ONEPT5*FR)*(ONEPT5*POVALU(KM)-HALF*POVALU(KM-1)-WR)
     &       *HS
         JR=I9
C
C      STORE THE FUNCTION AT RA
C
         IF(NUTTY.EQ.1)THEN
            US(JR)=U
            DUS(JR)=DU*TWO/H
         ELSE
            ORB(JR)=U
            DORB(JR)=DU*TWO/H
         END IF
C
C      INTEGRATE IN TO THE MATCHING POINT
C
         LSWT=1
         DO INT=1, NIX
            H1=H
            IMT=NIX-INT+1
            H=-HRX(IMT)
            HS=H*H/THREE
            IF(IMT.EQ.1)THEN
               I1=0
            ELSE
               I1=IRX(IMT-1)
            END IF
            I2=IRX(IMT)
            IF(IMATCH.GT.I1)I1=IMATCH
            IF(INT.GT.1)THEN
               LSWT=2
               H12=H/H1
            END IF
C
C      INTEGRATE OVER A RANGE OF EQUAL INTEGRALS
C
            DO IR=I1+1, I2
               JR=JR-1
               CALL DEVGL(DU,FR,FRM,HS,U,YR,POVALU,WR,lswt,mmm,h12,km)
C
C      STORE THE FUNCTIONS AT EACH INTERATION
C
               IF(NUTTY.NE.2)THEN
                  US(JR)=U
                  DUS(JR)=DU*TWO/H
               ELSE
                  ORB(JR)=U
                  DORB(JR)=DU*TWO/H
               END IF
               IF(JR.EQ.IMATCH+1)GO TO 81
            END DO
         END DO
C
C      STORE THE FUNCTIONS AND DERIVATIVES AT THE MATCHING POINT FOR
C      THE INWARD INTEGRATION
C
 81      IF(NUTTY.EQ.1)THEN
            YIN=U
            DYIN=DU*TWO/H
         END IF
C
      END DO
C
      DUH=DU*TWO/H
      B2=YOUT*DUH-U*DYOUT
      B3=YIN*DYOUT-YOUT*DYIN
C
C      EVALUATE THE FINAL UNNORMALIZED FUNCTION AT THE MESH POINTS
C
      DO I=1, IMATCH
         ORB(I)=US(I)
         DORB(I)=DUS(I)
      END DO
      DO I=IMATCH+1, I9
         ORB(I)=ORB(I)*B3+US(I)*B2
         DORB(I)=DORB(I)*B3+DUS(I)*B2
      END DO
      EPSILON=tiny(1.0_wp)   ! tiny() gives the smallest positive number, in the F95 standard onwards
      IF(ABS(ORB(I9))<EPSILON)THEN
         WRITE(NFTA,1016) LC
         STOP
      END IF
      BVALUE=B3*RA/ORB(I9) ! JMC possible site for division by zero, e.g. if yout and dyout are zero
                           ! (they come from raising a number < 1 to a power that may be large)
                           ! See the test (against tiny(1.0_wp)) above to trap such a case
C
C      NORMALIZE THE SOLUTION
C
      X1=ZERO
      I1=1
      DO I=1, NIX
         H=HRX(I)
         X2=ORB(I1)*ORB(I1)+FOUR*ORB(I1+1)*ORB(I1+1)
         I1=I1+2
         I2=IRX(I)-1
         DO J=I1, I2, 2
            X2=X2+TWO*ORB(J)*ORB(J)+FOUR*ORB(J+1)*ORB(J+1)
         END DO
         I1=I2+2
         X1=X1+(X2+ORB(I1)*ORB(I1))*H
      END DO
C
      IF(X1<EPSILON)THEN
         WRITE(NFTA,1017)
         STOP
      END IF
      X2=SQRT(THREE/X1)
      NODES=0
      IF(ORB(2).LT.ZERO)X2=-X2
C
C      EVALUATE NODES, THE NUMBER OF NODES IN THE FINAL FUNCTION
C
      ORB1=ORB(1)
      DO I=2, I9
         ORB2=ORB(I)*X2
         ORB(I)=ORB2
         DORB(I)=DORB(I)*X2
         IF(ORB1*ORB2.LT.ZERO .OR. (ORB1.NE.ZERO .AND. ORB2.EQ.ZERO))
     &      NODES=NODES+1                                      ! jmc comparison of real with 0.0 exactly???
         ORB1=ORB2
      END DO
C
      RETURN
      END SUBROUTINE BASFUN
!*==devgl.spg  processed by SPAG 6.56Rc at 11:55 on  9 Mar 2010
C
      SUBROUTINE DEVGL(DY,FR,FRM,HS,Y,YR,POVALU,WR,lswt,mmm,h12,km)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : EIGHTH, FOURTH=>QUART, HALF=>XHALF, ONE=>XONE
      USE NUMCBAS_DATA, ONLY : IRA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: DY, FR, FRM, H12, HS, WR, Y, YR
      INTEGER :: KM, LSWT, MMM
      REAL(KIND=wp), DIMENSION(2*ira) :: POVALU
      INTENT (IN) H12, HS, MMM, POVALU, WR
      INTENT (INOUT) DY, FR, FRM, KM, LSWT, Y, YR
C
C Local variables
C
      REAL(KIND=wp) :: H12S
C
C*** End of declarations rewritten by SPAG
C
C     DE VOGELAERE INTEGRATION ROUTINE
C
      IF(LSWT.EQ.1)THEN
C
C     FIXED MESH
         DY=DY+FR
         YR=Y+DY
         Y=YR+FR-EIGHTH*FRM
      ELSE
C
C     VARIABLE MESH
         H12S=H12*H12
         LSWT=1
         DY=H12*DY+H12S*FR
         YR=Y+DY
         Y=YR+HALF*H12S*(FR*(H12+ONE)-FOURTH*H12*FRM)
      END IF
C
      KM=KM+MMM
      FRM=(POVALU(KM)-WR)*HS*Y
      DY=DY+FRM
      Y=YR+DY
C
      KM=KM+MMM
      FR=(POVALU(KM)-WR)*HS*FOURTH*Y
      DY=DY+FR
C
      RETURN
      END SUBROUTINE DEVGL
!*==finder.spg  processed by SPAG 6.56Rc at 11:55 on  9 Mar 2010
C
      SUBROUTINE FINDER(LC,NODES,BCON,ETA,ETRIAL,EIGEN,POVALU,ORB,DORB,
     &                  ra)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, FOURTH=>QUART, HALF=>XHALF, 
     &                   ONE=>XONE, TEN=>XTEN
      USE NUMCBAS_DATA, ONLY : ABSACC, IRX, IRA, NIX, NFTA
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: BCON, EIGEN, ETA, ETRIAL, RA
      INTEGER :: LC, NODES
      REAL(KIND=wp), DIMENSION(ira+1) :: DORB, ORB
      REAL(KIND=wp), DIMENSION(2*ira) :: POVALU
      INTENT (IN) BCON, ETA, NODES
      INTENT (OUT) EIGEN
      INTENT (INOUT) ETRIAL
C
C Local variables
C
      REAL(KIND=wp) :: B1, B2, BHIGH, BLOW, BVALUE, DEL, DELE, 
     &                E1, E2, E3, EHIGH, ELOW, EPA, ESTEP, FDASH, 
     &                PMIN, TENETA, TINY=1.E-6_wp
C JMC The above initialization of TINY replaces a DATA statement.
C In either case, the value of TINY is saved between different
C instances of this routine.
      INTEGER :: I, IHERE, NHIGH, NODE
C
C*** End of declarations rewritten by SPAG
C
C      THIS ROUTINE LOCATES THE EIGENVALUE WITH A GIVEN NUMBER OF NODES
C
C     LC........ THE  ANGULAR MOMENTUM VALUE
C     NODES....  IS THE NUMBER OF NODES IN THE FUNCTION
C     BCON.....  IS THE TOLERANCE ON THE LOGARITHMIC DERIVATIVE
C     ETA......  IS THE TOLERANCE ON THE EIGENVALUE
C     ETRIAL...  IS THE INITIAL ESTIMATE FOR THE EIGENVALUE
C
C      INITIALIZE DEL AND THE ENERGY VALUES
C
      PMIN=-ONE
      DO I=1, IRX(NIX)
         PMIN=MIN(PMIN,POVALU(I))
      END DO
      B1=ZERO
      IF(TINY.GT.ETA)TINY=ETA
      E1=ZERO
      TENETA=TEN*ETA
      E2=ZERO
      E3=ZERO
      IHERE=0
      DEL=ONE/RA ! JMC RA should not be zero as long as hrx(1:nix) is OK
 1    CALL BASFUN(LC,NODE,ETRIAL,BVALUE,POVALU,ORB,DORB,ra)
      E1=E2
      E2=E3
      E3=ETRIAL
      IF(ABS(E1-E3).LT.TINY)DEL=DEL*HALF
      IF(node.EQ.nodes)GO TO 4
      IF(NODE.LT.NODES)THEN
         ETRIAL=ETRIAL+DEL
      ELSE
         ETRIAL=ETRIAL-DEL
         IF(ETRIAL.LT.PMIN)GO TO 6
      END IF
      GO TO 1
C
C      WE HAVE THE CORRECT NUMBER OF NODES
C
 4    IF(ABS(BVALUE).LT.TINY)GO TO 22
      IF(BVALUE.LT.zero)GO TO 9
      ELOW=ETRIAL
      BHIGH=BVALUE
      IHERE=1
 5    ETRIAL=ETRIAL+DEL
      CALL BASFUN(LC,NODE,ETRIAL,BVALUE,POVALU,ORB,DORB,ra)
      IF(node.EQ.nodes)GO TO 4
      IF(node.LT.nodes)GO TO 6
      IF(NODE.GT.NODES+1)THEN
         ETRIAL=ETRIAL-DEL
         DEL=DEL*HALF
         IF(DEL.LT.ABS(ETRIAL*ABSACC))GO TO 6
         GO TO 5
      END IF
C
C      WE HAVE AN UPPER BOUND TO THE EIGENVALUE
C
 9    EHIGH=ETRIAL
      BLOW=BVALUE
      NHIGH=NODE
      IF(IHERE.EQ.1)GO TO 14
 10   ETRIAL=ETRIAL-DEL
      CALL BASFUN(LC,NODE,ETRIAL,BVALUE,POVALU,ORB,DORB,ra)
      IF(node.GT.nodes)GO TO 6
      IF(node.LT.nodes)GO TO 13
      IF(BVALUE.LT.zero)GO TO 9
      ELOW=ETRIAL
      BHIGH=BVALUE
      GO TO 14
 13   ETRIAL=ETRIAL+DEL
      DEL=DEL*HALF
      IF(DEL.LT.ABS(ETRIAL*ABSACC))GO TO 6
      GO TO 10
C
C      WE NOW HAVE HIGH AND LOW ENERGIES TO THE EIGENVALUE
C
 14   IF(NHIGH.EQ.NODES)GO TO 19
C
C      WE NEED TO FIND CLOSER ENERGY BOUNDS
C
 15   ESTEP=(EHIGH-ELOW)*FOURTH
      IF(ESTEP.LT.ABS(ELOW*ABSACC))GO TO 116
      ETRIAL=ELOW
 16   ETRIAL=ETRIAL+ESTEP
      CALL BASFUN(LC,NODE,ETRIAL,BVALUE,POVALU,ORB,DORB,ra)
      IF(NODE.EQ.NODES)GO TO 17
      EHIGH=ETRIAL
      BLOW=BVALUE
      GO TO 15
 17   IF(BVALUE.LT.zero)GO TO 18
      ELOW=ETRIAL
      BHIGH=BVALUE
      GO TO 16
C
C      WE NOW HAVE THE EIGENVALUE PINPOINTED ACCURATELY ENOUGH
C
 18   EHIGH=ETRIAL
      BLOW=BVALUE
 19   IF(ABS(BHIGH-BLOW).GT.BCON)GO TO 15
      ETRIAL=HALF*(EHIGH+ELOW)
 20   EPA=-ETA
 21   CALL BASFUN(LC,NODE,ETRIAL,BVALUE,POVALU,ORB,DORB,ra)
      B2=B1
      B1=BVALUE
      EPA=EPA+ETA
      ETRIAL=ETRIAL+ETA
      IF(ABS(EPA-ETA).GE.TINY)GO TO 21
      FDASH=(B1-B2)/ETA
      DELE=-B1/FDASH
      ETRIAL=ETRIAL-ETA+DELE
      IF(ABS(DELE).GT.TENETA)GO TO 20
      CALL BASFUN(LC,NODE,ETRIAL,BVALUE,POVALU,ORB,DORB,ra)
      IF(NODE.NE.NODES)GO TO 15
C
 22   EIGEN=ETRIAL
      RETURN
C
 116  WRITE(NFTA,98) NODES, ELOW, BHIGH
 98   FORMAT(/' *** WARNING *** SOLUTION WITH',I3,' NODES,  ENERGY =',
     &       F10.4,/16x,' HAS LOG DERIV',D10.2,
     &       ' AND SHOULD BE DISCARDED')
      ETRIAL=ELOW
      CALL BASFUN(LC,NODE,ETRIAL,BVALUE,POVALU,ORB,DORB,ra)
      EIGEN=ETRIAL
      RETURN
C
 6    WRITE(NFTA,99) NODE, NODES, ETRIAL, BVALUE
 99   FORMAT(/' FAILED TO FIND EIGENVALUE WITH',I3,' NODES',2X,I5,
     &       2D15.6)
      STOP
      END SUBROUTINE FINDER
!*==wrhead.spg  processed by SPAG 6.56Rc at 11:55 on  9 Mar 2010
C
      SUBROUTINE WRHEAD(NCO,lval,ira,r,TITLE,lunumb)
      USE PRECISN, ONLY : WP                        
      USE consts, ONLY : ZERO => xzero
      USE global_utils, ONLY : CWBOPN
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IRA, LUNUMB, LVAL, NCO
      CHARACTER(LEN=*) :: TITLE
      REAL(KIND=wp), DIMENSION(ira) :: R
      INTENT (IN) IRA, LVAL, NCO, R, TITLE
C
C Local variables
C
      INTEGER :: IR
C
C*** End of declarations rewritten by SPAG
C
C     THIS SUBROUTINE WRITES THE HEADER OF FILE LUNUMB WHICH WILL CONTAIN
C     THE NUMERICAL BASIS IN A FORM SUITABLE FOR GAUSBAS
C
C
C     WRITE HEADER ON FILE LUNUMB
C
      CALL CWBOPN(LUNUMB)
      WRITE(LUNUMB)TITLE
      WRITE(LUNUMB)NCO, lval, ira
      WRITE(LUNUMB)zero, (r(IR),IR=1,IRA)
c
      RETURN
      END SUBROUTINE WRHEAD
