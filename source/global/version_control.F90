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
! Copyright 2018 
!
! Greg Armstrong, Jakub Benda, Andrew Brown, Daniel Clarke, Michael Lysaght,
! Zdenek Masin, Robert McGibbon, Laura Moore, Jonathan Parker, Martin Plummer,
! Ken Taylor, Hugo van der Hart, Jack Wragg                                    
!
! This file is part of RMT.
! 
!     RMT is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     RMT is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with RMT (in trunk/COPYING). Alternatively, you can also visit 
!     <https://www.gnu.org/licenses/>.
!

MODULE version_control

    IMPLICIT NONE

    PRIVATE

    PUBLIC print_git_revision

CONTAINS

    !> \brief    Prints Git revision info to standard output
    !> \authors  J Benda
    !> \date     2018
    !>
    !> The header contains information from the Git repository (if available). It is expected
    !> that the information is given to the compiler during translation of this unit. The following
    !> preprocessor macros are used: GIT_AUTH, GIT_HASH, GIT_DATE. They correspond to outputs
    !> of the following commands, in the same order:
    !> \verbatim
    !>     git log -1 --pretty=format:"%an"
    !>     git log -1 --pretty=format:"%h"
    !>     git log -1 --pretty=format:"%ad"
    !> \endverbatim
    !>
    !> \param lu  Output unit
    !>
    SUBROUTINE print_git_revision (lu)
        INTEGER, INTENT(IN) :: lu
#ifdef GIT_AUTH
        CHARACTER(LEN=57), PARAMETER :: str_auth = GIT_AUTH
#else
        CHARACTER(LEN=57), PARAMETER :: str_auth = 'unknown'
#endif
#ifdef GIT_HASH
        CHARACTER(LEN=57), PARAMETER :: str_hash = GIT_HASH
#else
        CHARACTER(LEN=57), PARAMETER :: str_hash = 'unknown'
#endif
#ifdef GIT_DATE
        CHARACTER(LEN=57), PARAMETER :: str_date = GIT_DATE
#else
        CHARACTER(LEN=57), PARAMETER :: str_date = 'unknown'
#endif
        WRITE(lu, '("+-------------------------------------------------------------------------+")')
        WRITE(lu, '("| Last revision: ",A57,                                                  "|")') str_hash
        WRITE(lu, '("| Last author:   ",A57,                                                  "|")') str_auth
        WRITE(lu, '("| Last date:     ",A57,                                                  "|")') str_date
        WRITE(lu, '("+-------------------------------------------------------------------------+")')

    END SUBROUTINE

END MODULE version_control
