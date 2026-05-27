      MODULE NUMCBAS_DATA
C JMC module for collecting various parameters used within numcbas
      USE precisn, ONLY : wp
      IMPLICIT NONE
      PUBLIC ! as it only contains data at present
      SAVE ! for the variables from the BASCON common block

      INTEGER, PARAMETER :: MAXORB=100

      REAL(KIND=wp), PARAMETER :: RTOL=5.E-02_wp

      REAL(KIND=wp), PARAMETER :: TINYY=1.E-6_wp

      REAL(KIND=wp), PARAMETER :: ABSACC=1.E-14_wp

      INTEGER, PARAMETER :: MAXRX=10

      INTEGER, PARAMETER :: NFTA=6 ! unit number for printing

C Variables from the BASCON common block
      REAL(KIND=wp), DIMENSION(MAXRX) :: HRX
      INTEGER :: IRA, NIX
      INTEGER, DIMENSION(MAXRX) :: IRX

      END MODULE NUMCBAS_DATA
