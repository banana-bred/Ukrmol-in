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
program integral_interface
   use precisn_gbl
   use const_gbl, only: line_len, stdout, set_verbosity_level, one_elham, two_el_ints
   use mpi_gbl
   use atomic_basis_gbl
   use molecular_basis_gbl
   use parallel_arrays_gbl
   use integral_storage_gbl

   implicit none

   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
   type(molecular_orbital_basis_obj) :: molecular_orbital_basis
   type(p2d_array_obj), target :: one_electron_integrals, two_electron_integrals
   type(integral_options_obj) :: options

   character(len=line_len) :: path
   integer :: err, length, one_elham_column, two_el_column
   integer :: two_ind(2,1), four_ind(4,1), ind(1)
   real(kind=wp) :: x

      !Start MPI: this program has been tested so far only with one MPI task
      call mpi_mod_start

      !How much information should the library produce. 3 = most detailed output.
      call set_verbosity_level(3)

      !Command line argument is the path to the integral file produced by scatci_integrals
      call get_command_argument(1,path,length,err)
      if (err > 0) then
         call mpi_xermsg('get_integrals_template', 'main', &
                         'The path to the data file has not been given: it must be given as the first command line argument.', 1, 1)
      end if

      !Read the atomic and molecular basis sets and the one electron and two electron integrals
      err = molecular_orbital_basis%read_basis_and_integrals(path, &
                                                             atomic_orbital_basis, &
                                                             one_electron_integrals, &
                                                             two_electron_integrals, &
                                                             options)

      !Which column of the integal array contains the 1-electron Hamiltonian matrix elements
      one_elham_column = one_electron_integrals%find_column_matching_name(one_elham)

      !Column of the integal array containing the 2-electron matrix elements (always 1st one)
      two_el_column = 1

      call molecular_orbital_basis%print

      call molecular_orbital_basis%print_energy_sorted_orbital_table

      !Absolute indices of the orbitals for which we want the 1-electron integral
      two_ind(1,1) = 1
      two_ind(2,1) = 1

      !Compute the integral index for the 1-electron integral
      ind(1:1) = molecular_orbital_basis%integral_index(one_elham,two_ind,options%two_p_continuum)

      !Retrieve the integral from the storage. Even zero integrals are stored.
      x = one_electron_integrals%a(ind(1),one_elham_column)

      write(stdout,'("1-electron integral")')
      write(stdout,'(2i0,e25.15)') two_ind(1:2,1), x

      !The same for the 2-electron integral
      four_ind(1,1) = 1
      four_ind(2,1) = 1
      four_ind(3,1) = 1
      four_ind(4,1) = 1

      ind(1:1) = molecular_orbital_basis%integral_index(two_el_ints,four_ind,options%two_p_continuum)

      !If the index is not positive then the integral is zero 
      x = 0.0_wp
      if (ind(1) > 0) x = two_electron_integrals%a(ind(1), two_el_column)

      write(stdout,'("2-electron integral")')
      write(stdout,'(4i0,e25.15)') four_ind(1:4,1), x

!FINALIZE

      err = molecular_orbital_basis%final()
      if (err /= 0) call mpi_xermsg('get_integrals_template','main','Error finalizing molecular_orbital_basis.',2,1)

      err = one_electron_integrals%final()
      if (err /= 0) call mpi_xermsg('get_integrals_template','main','Error finalizing 1-el integrals.',3,1)

      err = two_electron_integrals%final()
      if (err /= 0) call mpi_xermsg('get_integrals_template','main','Error finalizing 2-el integrals.',4,1)

      call mpi_mod_finalize

end program
