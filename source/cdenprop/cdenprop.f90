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
program cdenprop

    use global_utils, only: print_ukrmol_header
    use cdenprop_procs
    use cdenprop_defs

    implicit none

    type (CIvect) :: ci_vec_i, ci_vec_j !The N+1 CI vectors

    call print_ukrmol_header(6)

    write(6,'(/,20X,"CDENPROP (2011) ",//,20X,"Transition dipoles program ",//)')

    call cdenprop_drv(ci_vec_i, ci_vec_j) !the ci_vec_* structures are only dummy structures for this cdenprop call

end program
