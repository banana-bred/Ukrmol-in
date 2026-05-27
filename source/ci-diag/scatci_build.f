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
!*==scatci.spg  processed by SPAG 6.56Rc at 17:35 on  9 Nov 2010
      PROGRAM SCATCI
C
C     MAIN PROGRAM
C     CFT77 version: all ENTRY points removed
C
      USE global_utils, ONLY : UTILS_DATE_TIME, PRINT_UKRMOL_HEADER
      USE SCATCI_ROUTINES, ONLY : EXDRVF
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      CHARACTER(LEN=8) :: CURDAT
      CHARACTER(LEN=10) :: TIM
C
C*** End of declarations rewritten by SPAG
C
c --- Date stamp output
      CALL UTILS_DATE_TIME(curdat,tim)
      CALL PRINT_UKRMOL_HEADER(6)
      WRITE(6,10)CURDAT, tim
 10   FORMAT(//' PROGRAM SCATCI (19 Oct 2009) : '//24X,' DATE = ',
     &       A8//24X,' TIME = ',a10,//)
C
C---- Call main routine
C
      CALL EXDRVF
C
      STOP
      END PROGRAM SCATCI

