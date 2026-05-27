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
   use basis_data_generic_mod
   use atomic_basis_mod
   use molden_mod
   use symmetry
   use const, only: line_len
   implicit none

   type(atomic_orbital_basis_obj) :: atomic_orbital_basis
   type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)
   type(orbital_data_obj), allocatable :: orbital_data(:)
   type(geometry_obj) :: symmetry_inp
   type(molden_input_obj) :: molden_input
   character(len=line_len) :: molden_file
   integer :: err, n, io, nob(8), i

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

      n = molden_input%get_no_cgto_shells()
      err = atomic_orbital_basis%init(n,symmetry_inp)

      !read-in the GTO shell data
      call molden_input%read_basis(CGTO_shell_data)

      !add it to the basis
      do i=1,n
         call atomic_orbital_basis%add_shell(CGTO_shell_data(i))
      enddo

      nob = 0
      nob(1:8) = (/16,6,12,4,13,5,11,4/)

      !read-in the orbitals
      call molden_input%read_orbitals(CGTO_shell_data,nob,orbital_data)

      call molden_input%final()

      !open the Molden file for output
      molden_file = './molden.test'
      io = 2
      call molden_input%init(molden_file,io)

      call molden_input%write_geometry(atomic_orbital_basis%symmetry_data%nucleus)
      call molden_input%write_basis(CGTO_shell_data)
      call molden_input%write_orbitals(atomic_orbital_basis,orbital_data)

      call molden_input%final()

      !Set point group symmetry ID of all orbital sets
      i = atomic_orbital_basis%symmetry_data%get_pg()
      orbital_data(:)%point_group = i

      call atomic_orbital_basis%write('./basis')
      call atomic_orbital_basis%read('./basis')

      !open the Molden file for output
      molden_file = './molden.copy'
      io = 2
      call molden_input%init(molden_file,io)

      call molden_input%write_geometry(atomic_orbital_basis%symmetry_data%nucleus)
      call molden_input%write_basis(CGTO_shell_data)

      call mpi_mod_finalize

end program test_atomic_basis_mod
