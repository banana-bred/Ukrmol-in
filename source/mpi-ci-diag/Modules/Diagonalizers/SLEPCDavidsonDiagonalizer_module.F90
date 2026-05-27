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

!> \brief   Diagonalizer type using SLEPc (Davidson) backend
!> \authors A Al-Refaie
!> \date    2017
!>
!> Based on SLEPCDiagonalizer. This module is only included in the build if SLEPC_LIBRARIES are given
!> on the CMake command line.
!>
!> \todo Stub only, no actual functionality implemented.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module SLEPCDavidsonDiagonalizer_module

    use const_gbl,                only: stdout
    use mpi_gbl,                  only: mpi_xermsg
    use SLEPCDiagonalizer_module, only: SLEPCDiagonalizer
    use slepceps

    implicit none

#include <finclude/slepceps.h>

    public SLEPCDavidsonDiagonalizer

    private

    type,extends(SLEPCDiagonalizer) :: SLEPCDavidsonDiagonalizer
    contains
        procedure, public :: SelectEPS => Select_Davidson
    end type

contains

    subroutine Select_Davidson(this,eps)
        class(SLEPCDavidsonDiagonalizer) :: this
        EPS, intent(inout)               :: eps
        PetscErrorCode                   :: ierr

        write(stdout,"('--------Generalized Davidson used as Diagonalizer------')")

        !     ** Create eigensolver context
       !call PetscSetFPTrap(PETSC_FP_TRAP_ON, ierr)
        call EPSSetType(eps, EPSGD, ierr)
       !call EPSGDSetDoubleExpansion(eps, PETSC_TRUE)
        if (ierr /= 0) then
            call mpi_xermsg('SLEPCDavidsonDiagonalizer_module', 'Select_Davidson', 'Failed to set solver context', int(ierr), 1)
        end if

    end subroutine Select_Davidson

end module SLEPCDavidsonDiagonalizer_module
