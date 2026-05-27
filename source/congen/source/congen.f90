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
!> \brief Program CONGEN
!>
!> This is the CONGEN entry point. Very little work is done here, actually, only retrieval of the current date
!> and time, write-out of the program header with SVN revision information, and of used real and integer byte
!> sizes. Otherwise, all work is done in the main subroutine \ref congen_driver::csfgen.
!>
!> Text written by the program before the call to \ref congen_driver::csfgen is printed to the standard output (or anywhere
!> the initial \ref congen_data::nftw points to). Only later output can be redirected to a file by resetting the value of
!> \ref congen_data::nftw in the namelist &state/.
!>
program congen

    use global_utils,  only : utils_get_bytesizes, utils_date_time, print_ukrmol_header
    use congen_data,   only : lratio, nftw
    use congen_driver, only : csfgen

    implicit none

    integer :: linteg, lreal
    character(len=8) :: curdat
    character(len=10) :: tim

    call utils_date_time (curdat, tim)
    call print_ukrmol_header (nftw)

    write(nftw,'("1",//," PROGRAM CONGEN :",//,24X," DATE = ",A8//,24X," TIME = ",A10,//)') curdat, tim
  
    call utils_get_bytesizes (linteg, lreal)
    lratio = lreal / linteg

    write(nftw,'(//"Using ",I3,"-byte integers"/"  and ",I3,"-byte reals"/)') linteg, lreal

    call csfgen

end program congen
