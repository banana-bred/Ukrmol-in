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
   use molecular_basis_mod
   use molden_mod
   use symmetry
   use const, only: line_len
   use integral_storage_mod
   implicit none

   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
   type(molecular_orbital_basis_obj) :: molecular_orbital_basis

   type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)
   type(orbital_data_obj), allocatable :: orbital_data(:)
   type(integral_options_obj) :: integral_options
   type(integral_storage_obj) :: integral_storage
   type(p2d_array_obj), target :: ao_integrals
   type(geometry_obj) :: symmetry_inp
   type(sym_ortho_io) :: sym_ortho
   type(molden_input_obj) :: molden_input
   character(len=line_len) :: molden_file
   integer :: err, n, io, nob(8), i
   real(kind=cfp), allocatable :: overlap_matrix(:,:)

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

      !read-in the GTO basis set and initialize the symmetry data
      call molden_input%read_basis(CGTO_shell_data)

      n = size(CGTO_shell_data)
      err = atomic_orbital_basis%init(n,symmetry_inp)

      do i=1,n
         call atomic_orbital_basis%add_shell(CGTO_shell_data(i))
      enddo 

      nob = 0
      nob(1:8) = (/16,6,12,4,13,5,11,4/)

      !read-in the orbitals
      call molden_input%read_orbitals(CGTO_shell_data,nob,orbital_data)

      call molden_input%final()

      !Set point group symmetry ID of all orbital sets
      i = atomic_orbital_basis%symmetry_data%get_pg()
      orbital_data(:)%point_group = i

      molecular_orbital_basis%ao_basis => atomic_orbital_basis
      err = molecular_orbital_basis%init(8,symmetry_inp)
      if (err .ne. 0) stop "error initializing molecular orbital basis"

      do n=1,8
         call molecular_orbital_basis%add_shell(orbital_data(n))
      enddo

      call molecular_orbital_basis%print
      call molecular_orbital_basis%print_orbitals

      !describe where the integrals will be stored 
      err = integral_storage%init(memory=ao_integrals) !(disk='./oneel')
      if (err .ne. 0) then
         print *,err
         call mpi_xermsg('main','main','error initializing the target ao_integral_storage',1,1)
      endif

      !Calculate the atomic overlap integrals:
      integral_options%a = 10.0_cfp !-1.0_cfp !R-matrix radius?
      integral_options%two_p_continuum = .false. !Calculate integrals for 2p in the continuum?
      integral_options%calculate_overlap_ints = .true.
      integral_options%calculate_kinetic_energy_ints = .true.
      integral_options%calculate_nuclear_attraction_ints = .true.
      integral_options%calculate_one_el_hamiltonian_ints = .true.
      integral_options%calculate_property_ints = .true.
      integral_options%max_property_l = 2
      
      call atomic_orbital_basis%one_electron_integrals(integral_storage,integral_options)

      call atomic_orbital_basis%assemble_overlap_matrix(ao_integrals,overlap_matrix)

      do i=1,8
         !GS orthogonalize the target orbitals among themselves
         if (nob(i) > 0) call molecular_orbital_basis%orthogonalize(gramm_schmidt=.true.,overlap_matrix=overlap_matrix,symmetry=i,&
                             &active_start=1,active_end=nob(i),check_overlaps=.true.)
      enddo

      call atomic_orbital_basis%write('./basis')
      call molecular_orbital_basis%write('./basis')

      call atomic_orbital_basis%read('./basis')
      call molecular_orbital_basis%read('./basis')

      call molecular_orbital_basis%print
      call molecular_orbital_basis%print_orbitals

      sym_ortho%del_thrs = 1.0e-08_cfp
      do i=1,8
         !symmetrically orthogonalize the target orbitals among themselves
         if (nob(i) > 0) call molecular_orbital_basis%orthogonalize(symmetric=.true.,sym_ortho_data=sym_ortho,overlap_matrix=overlap_matrix,symmetry=i,&
                             &active_start=1,active_end=nob(i),check_overlaps=.true.)
      enddo

      call mpi_mod_finalize

end program test_atomic_basis_mod
