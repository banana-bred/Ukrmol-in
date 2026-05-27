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

!> \brief   ELPA matrix module
!> \authors J Benda
!> \date    2019
!>
!> Only compiled in when ELPA library is available during the configuration step.
!>
module ELPAMatrix_module

    use blas_lapack_gbl,        only: blasint
    use Parallelization_module, only: grid => process_grid
    use precisn_gbl,            only: wp
    use SCALAPACKMatrix_module, only: SCALAPACKMatrix
    use Timing_module,          only: master_timer

    implicit none

    private

    !> \brief   ELPA distributed matrix
    !> \authors J Benda
    !> \date    2019
    !>
    !> This matrix class is exactly the same as the SCALAPACKMatrix, except that the matrix finalization routine also
    !> takes care of explicit symmetrization; SCALAPACK is fine with just one triangle, but ELPA requires that both are set.
    !>
    type, public, extends(SCALAPACKMatrix) :: ELPAMatrix
    contains
        procedure, public :: finalize_matrix => finalize_ELPA
    end type

contains

    subroutine finalize_ELPA (this)

        class(ELPAMatrix) :: this
        integer(blasint)  :: one = 1, i, k, l, iprow, ipcol
        real(wp)          :: alpha = 1

        ! finalize parent
        call this % SCALAPACKMatrix % finalize_matrix

        ! perform symmetrization (add transposed triangle)
        call master_timer % start_timer("Symmetrize matrix")
        call pdtradd('U', 'T', this % mat_dimen, this % mat_dimen, &
                     alpha, this % a_local_matrix, one, one, this % descr_a_mat, &
                     alpha, this % a_local_matrix, one, one, this % descr_a_mat, this % a_local_matrix)

        ! scale the diagonal
        do i = 1, this % mat_dimen
            call infog2l(i, i, this % descr_a_mat, grid % gprows, grid % gpcols, grid % mygrow, grid % mygcol, k, l, iprow, ipcol)
            if (iprow == grid % mygrow .and. ipcol == grid % mygcol) then
                this % a_local_matrix(k, l) = this % a_local_matrix(k, l) / 2
            end if
        end do
        call master_timer % stop_timer("Symmetrize matrix")

    end subroutine finalize_ELPA

end module ELPAMatrix_module
