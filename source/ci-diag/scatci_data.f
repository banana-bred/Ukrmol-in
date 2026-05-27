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
      MODULE SCATCI_DATA
C
C Module containing parameters / data required throughout the
C scatci program.
C
C Joanne Carr, November 2010
C
      USE precisn, ONLY : wp ! for specifying the kind of reals
      IMPLICIT NONE
      SAVE

c     NTGTMX defines the size os several arrays. It is the maximum
c     number of target electronic states/symmetries to be used.
      INTEGER, PARAMETER :: NTGTMX=1000

      INTEGER, PARAMETER :: MEIG=400

      END MODULE SCATCI_DATA
