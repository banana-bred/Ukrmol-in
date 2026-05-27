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
program test_atomic_basis_mod
   use mpi_mod
   use atomic_basis_mod
   use molden_mod
   use symmetry
   use const, only: line_len
   use precisn
   implicit none

   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
   type(BTO_shell_data_obj) :: BTO_shell_data
   type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)

   type(geometry_obj) :: symmetry_inp
   type(molden_input_obj) :: molden_input
   character(len=line_len) :: molden_file
   integer :: err, n, io, i

      call mpi_mod_start

      !open and analyze the Molden input file
      molden_file = '../level2/common_obj/test/bin/pyrazine.cation.molden'
      io = 1
      call molden_input%init(molden_file,io)

      !read the target geometry and the number of nuclei
      call molden_input%read_geometry(symmetry_inp%nucleus,symmetry_inp%no_nuc)

      !set-up the symmetry input
      symmetry_inp%no_sym_op = 3 !no_sym_op
      symmetry_inp%sym_op(1:3) = (/'X','Y','Z'/) !sym_op(1:3)
      symmetry_inp%use_symmetry = .true.

      BTO_shell_data%l = 4
      BTO_shell_data%bspline_index = 5
      BTO_shell_data%number_of_functions = 2*BTO_shell_data%l+1
      BTO_shell_data%non_zero_at_boundary = .true.
      call BTO_shell_data%bspline_grid%init_grid(5.0_cfp,30.0_cfp,10,34)

      !read-in the GTO basis set and initialize the symmetry data
      call molden_input%read_basis(CGTO_shell_data)
      call molden_input%final()

      n = size(CGTO_shell_data)
      err = atomic_orbital_basis%init(n+1,symmetry_inp)

      do i=1,n
         call atomic_orbital_basis%add_shell(CGTO_shell_data(i))
      enddo

      call atomic_orbital_basis%add_shell(BTO_shell_data)

      call atomic_orbital_basis%write('./basis')
      call atomic_orbital_basis%read('./basis')

      call mpi_mod_finalize

end program test_atomic_basis_mod
