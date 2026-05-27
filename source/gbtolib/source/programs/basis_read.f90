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
program test
   use precisn_gbl
   use const_gbl, only: line_len, stdout, set_verbosity_level
   use mpi_gbl
   use atomic_basis_gbl
   use molecular_basis_gbl

   implicit none

   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
   type(molecular_orbital_basis_obj) :: molecular_orbital_basis

   character(len=line_len) :: path, val
   real(kind=cfp) :: A,B,delta_r,rmat_radius
   integer :: err, length
   logical :: calculate_radial_densities
   real(kind=cfp), allocatable :: amplitudes(:,:)

      call mpi_mod_start
      call set_verbosity_level(3)

      call get_command_argument(1,path,length,err)
      if (err > 0) then
         call mpi_xermsg('basis_read', 'main', &
                         'The path to the data file has not been given: it must be given as the first command line argument.', 1, 1)
      end if

      call atomic_orbital_basis%read(path)

      molecular_orbital_basis%ao_basis => atomic_orbital_basis
      call molecular_orbital_basis%read(path)

      call molecular_orbital_basis%print

      call molecular_orbital_basis%print_energy_sorted_orbital_table

      call molecular_orbital_basis%print_orbitals

      !Calculate radial charge densities of all orbitals.
      calculate_radial_densities = .true.

      call get_command_argument(2,val,length,err)
      if (err > 0) calculate_radial_densities = .false.
      if (calculate_radial_densities) read(val,*) rmat_radius

      call get_command_argument(3,val,length,err)
      if (err > 0) calculate_radial_densities = .false.
      if (calculate_radial_densities) read(val,*) A

      call get_command_argument(4,val,length,err)
      if (err > 0) calculate_radial_densities = .false.
      if (calculate_radial_densities) read(val,*) B

      call get_command_argument(5,val,length,err)
      if (err > 0) calculate_radial_densities = .false.
      if (calculate_radial_densities) read(val,*) delta_r

      if (calculate_radial_densities) then
         write(stdout,'("Radial charge densities will be calculated on the grid (A,B,delta_r): ",3e25.15)') A,B,delta_r
         write(stdout,'("For R-matrix radius: ",e25.15)') rmat_radius
         call molecular_orbital_basis%radial_charge_density(rmat_radius,A,B,delta_r,.true.,amplitudes)
      endif

!FINALIZE

      err = molecular_orbital_basis%final()
      if (err .ne. 0) call mpi_xermsg('basis_read','main','Error finalizing molecular_orbital_basis.',2,1)

      call mpi_mod_finalize

end program test
