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
   use free_scattering_mod
   use molden_mod
   use symmetry
   use integral_storage_mod
   use const, only: line_len
   implicit none

   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
   type(molecular_orbital_basis_obj) :: molecular_orbital_basis
   type(BTO_shell_data_obj) :: BTO_shell_data
   type(orbital_data_obj), allocatable :: orbital_data(:)
   type(integral_options_obj) :: integral_options
   type(integral_storage_obj), target :: atomic_integral_storage, molecular_integral_storage
   type(p2d_array_obj), target :: ao_integrals, mo_integrals
   type(sym_ortho_io) :: sym_ortho

   type(geometry_obj) :: symmetry_inp
   integer :: err, n, i
   real(kind=cfp) :: min_energy = 0.0_cfp, max_energy = 3.0_cfp
   integer :: nE = 300
   logical, allocatable :: is_continuum(:)
   real(kind=cfp), allocatable :: overlap_matrix(:,:)

      call mpi_mod_start

      !set-up the symmetry input
      symmetry_inp%no_sym_op = 3 !no_sym_op
      symmetry_inp%sym_op(1:3) = (/'X','Y','Z'/) !sym_op(1:3)
      symmetry_inp%use_symmetry = .true.

      symmetry_inp%no_nuc = 1
      allocate(symmetry_inp%nucleus(symmetry_inp%no_nuc))

      symmetry_inp%nucleus(1)%name = "D"
      symmetry_inp%nucleus(1)%center = 0.0_cfp

      BTO_shell_data%l = 4
      BTO_shell_data%number_of_functions = 2*BTO_shell_data%l+1
      BTO_shell_data%non_zero_at_boundary = .true.
      call BTO_shell_data%bspline_grid%init_grid(0.0_cfp,30.0_cfp,10,34)

      n = BTO_shell_data%bspline_grid%n
      err = atomic_orbital_basis%init(n-2,symmetry_inp) !exclude those radial B-splines with a non-zero derivative at r=A

      do i=3,n
         BTO_shell_data%bspline_index = i
         call atomic_orbital_basis%add_shell(BTO_shell_data)
      enddo

      !describe where the AO integrals will be stored 
      err = atomic_integral_storage%init(memory=ao_integrals) !(disk='./oneel')
      if (err .ne. 0) then
         print *,err
         call mpi_xermsg('main','main','error initializing the target atomic_integral_storage',1,1)
      endif

      integral_options%a = BTO_shell_data%bspline_grid%B !-1.0_cfp !R-matrix radius?
      integral_options%calculate_overlap_ints = .true.
      integral_options%calculate_kinetic_energy_ints = .true.
      
      call atomic_orbital_basis%one_electron_integrals(atomic_integral_storage,integral_options)

      call ao_integrals%print(.true.)

      !Set point group symmetry ID of all orbital sets
      allocate(orbital_data(8))
      i = atomic_orbital_basis%symmetry_data%get_pg()
      orbital_data(:)%point_group = i

      molecular_orbital_basis%ao_basis => atomic_orbital_basis
      err = molecular_orbital_basis%init(8,symmetry_inp)
      if (err .ne. 0) stop "error initializing molecular orbital basis"

      do n=1,8
         call atomic_orbital_basis%get_continuum_flags(n,is_continuum) !mark those AOs with IRR=n which represent continuum functions centered on the CMS.
         orbital_data(n)%irr = n
         call orbital_data(n)%add_cms_continuum(is_continuum) !add the initial set of continuum functions for this symmetry to the orbital data set
         call molecular_orbital_basis%add_shell(orbital_data(n))
      enddo

      call atomic_orbital_basis%assemble_overlap_matrix(ao_integrals,overlap_matrix)

      do i=1,8
         !symmetrically orthogonalize the continuum orbitals among themselves
         if (orbital_data(i)%number_of_functions > 0) then
            sym_ortho%del_thrs = 1e-14_cfp !del_thrs
            call molecular_orbital_basis%orthogonalize(symmetric=.true.,sym_ortho_data=sym_ortho,overlap_matrix=overlap_matrix,symmetry=i,&
                                                      &active_start=1,active_end=orbital_data(i)%number_of_functions,check_overlaps=.true.)
            call molecular_orbital_basis%delete_orbitals(i,sym_ortho%to_delete)
         endif
      enddo !i

      !describe where the MO integrals will be stored 
      err = molecular_integral_storage%init(memory=mo_integrals) !(disk='./oneel')
      if (err .ne. 0) then
         print *,err
         call mpi_xermsg('main','main','error initializing the target molecular_integral_storage',1,1)
      endif

      molecular_orbital_basis%ao_integral_storage => atomic_integral_storage !point to the storage for the atomic integrals
      call molecular_orbital_basis%one_electron_integrals(molecular_integral_storage,integral_options)

      call free_scattering(molecular_integral_storage,molecular_orbital_basis,integral_options%a,min_energy,max_energy,nE)

      call mpi_mod_finalize

end program test_atomic_basis_mod
