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
   use symmetry
   use precisn
   implicit none

   type(atomic_orbital_basis_obj) :: atomic_orbital_basis
   type(CGTO_shell_data_obj) :: CGTO_shell_data
   type(geometry_obj) :: symmetry_inp
   integer :: err, i, n

      call mpi_mod_start

      CGTO_shell_data%l = 2
      CGTO_shell_data%number_of_functions = 2*CGTO_shell_data%l+1
      CGTO_shell_data%center(1:3) = (/-1.5_cfp,2.5_cfp,0.25_cfp/)
      CGTO_shell_data%number_of_primitives = 3
      allocate(CGTO_shell_data%exponents(CGTO_shell_data%number_of_primitives),CGTO_shell_data%contractions(CGTO_shell_data%number_of_primitives),CGTO_shell_data%norms(CGTO_shell_data%number_of_primitives))
      CGTO_shell_data%exponents(1:CGTO_shell_data%number_of_primitives) = (/0.25_cfp,400.0_cfp,15.0_cfp/)
      CGTO_shell_data%contractions(1:CGTO_shell_data%number_of_primitives) = (/-2.25_cfp,4.0_cfp,1.5_cfp/)
      CGTO_shell_data%non_zero_at_boundary = .false.

      !set-up the symmetry input
      symmetry_inp%no_sym_op = 0
      symmetry_inp%sym_op(1:3) = "" !sym_op(1:3)
      symmetry_inp%use_symmetry = .true.
      symmetry_inp%no_nuc = 1
      allocate(symmetry_inp%nucleus(symmetry_inp%no_nuc))
      symmetry_inp%nucleus(1)%center = CGTO_shell_data%center
      symmetry_inp%nucleus(1)%nuc = 1
      symmetry_inp%nucleus(1)%charge = 1.0_cfp
      symmetry_inp%nucleus(1)%name = "T"

      n = 3
      err = atomic_orbital_basis%init(n,symmetry_inp)

      do i=1,n
         call atomic_orbital_basis%add_shell(CGTO_shell_data)
      enddo !i

      call atomic_orbital_basis%write('./basis')
      call atomic_orbital_basis%read('./basis')

      call mpi_mod_finalize

end program test_atomic_basis_mod
