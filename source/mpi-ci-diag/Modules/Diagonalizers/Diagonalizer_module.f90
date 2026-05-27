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

!> \brief   Base type for all diagonalizers
!> \authors A Al-Refaie
!> \date    2017
!>
!> Individual specializations are included in the build depending on the CMake configuration flags.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Diagonalizer_module

    implicit none

    type, abstract :: BaseDiagonalizer
    contains
        procedure(generic_diagonalize), deferred :: diagonalize
    end type BaseDiagonalizer

    abstract interface
        subroutine generic_diagonalize(this, matrix_elements, num_eigenpair, dresult, all_procs, option, integrals)
            use DiagonalizerResult_module, only: DiagonalizerResult
            use BaseIntegral_module,       only: BaseIntegral
            use BaseMatrix_module,         only: BaseMatrix
            use Options_module,            only: Options
            use precisn,                   only: wp

            import :: BaseDiagonalizer

            class(BaseDiagonalizer)         :: this
            class(DiagonalizerResult)       :: dresult
            class(BaseMatrix),   intent(in) :: matrix_elements
            class(BaseIntegral), intent(in) :: integrals
            type(Options),       intent(in) :: option
            integer,             intent(in) :: num_eigenpair
            logical,             intent(in) :: all_procs
        end subroutine generic_diagonalize
    end interface

end module Diagonalizer_module
