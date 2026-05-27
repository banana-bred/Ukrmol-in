! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
program test_molden_print
   use mpi_gbl
   use atomic_basis_gbl
   use molecular_basis_gbl
   use molden_gbl
   use const_gbl, only: line_len, sym_op_nam_len
   use precisn_gbl
   implicit none

   type(molden_input_obj) :: molden_input

   integer :: io, idat

   !namelist variables
   character(len=line_len) :: molden_file = ''
   integer :: nob(1:8) = -1, no_sym_op = 0
   character(len=sym_op_nam_len) :: sym_op(1:3) = (/'  ','  ','  '/)
   real(kind=cfp) :: a = -1.0_cfp

   namelist /target_data/ molden_file, nob, no_sym_op, sym_op, a

      call mpi_mod_start

      open(newunit=idat,file='./inp',form='formatted',status='unknown')

      read(idat,nml=target_data)

      close(idat)

      !open and analyze the Molden input file
      io = 1
      call molden_input%init(molden_file,io)

      call molden_input%print_energy_sorted_orbital_table

      !close the Molden file
      call molden_input%final()

      call mpi_mod_finalize

end program test_molden_print
