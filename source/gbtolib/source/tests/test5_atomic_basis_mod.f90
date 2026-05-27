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
   use const, only: line_len, sym_op_nam_len
   use integral_storage_mod
   use precisn
   implicit none

   !These parameters are not a limitation on the maximum number of continuum exponents and maximum L but only on the size of the auxiliary arrays below.
   integer, parameter :: max_exp = 20
   integer, parameter :: max_cont_l = 15

   type(atomic_orbital_basis_obj) :: atomic_orbital_basis
   type(CGTO_shell_data_obj), allocatable :: CGTO_shell_data(:)
   type(geometry_obj) :: symmetry_inp
   type(molden_input_obj) :: molden_input
   type(integral_options_obj) :: integral_options
   type(integral_storage_obj) :: integral_storage
   type(p2d_array_obj), target :: ao_integrals

   integer :: err, n, io, i, j, number_of_continuum_shells, last_irr, no_cont_exps(0:max_cont_l), number_of_target_shells

   !namelist variables
   character(len=line_len) :: molden_file = '', mo_integrals_file_name = './moints', ao_integrals_file_name = './aoints'
   integer :: nob(1:8) = -1, no_sym_op = 0, max_l = -1, n_mpi_per_numa = -1, nE = -1
   character(len=sym_op_nam_len) :: sym_op(1:3) = (/'  ','  ','  '/)
   real(kind=cfp) :: del_thrs(1:8) = (/-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp,-1.0_cfp/)
   real(kind=cfp) :: exponents(1:max_exp,0:max_cont_l) = 0.0_cfp
   real(kind=cfp) :: a = -1.0_cfp, min_energy = 0.0_cfp, max_energy = 1.0_cfp
   logical :: run_free_scattering = .false., save_ao_integrals_to_disk = .false., do_two_particle_integrals = .true., two_p_continuum = .false.
   logical :: check_target_target_orbital_overlaps = .true., check_target_continuum_orbital_overlaps = .true., check_continuum_continuum_orbital_overlaps = .true.

   namelist /target_data/ molden_file, nob, no_sym_op, sym_op, a
   namelist /continuum_data/ del_thrs, exponents, max_l, run_free_scattering, min_energy, max_energy, nE
   namelist /process_control/ n_mpi_per_numa, ao_integrals_file_name, mo_integrals_file_name, save_ao_integrals_to_disk, do_two_particle_integrals, &
                             &check_target_target_orbital_overlaps, check_target_continuum_orbital_overlaps, check_continuum_continuum_orbital_overlaps, two_p_continuum

      call mpi_mod_start

      !this routine is contained in this program at the bottom: read-in the namelist variables and check for sanity of the input
      !todo the inp file should be passed as a command-line argument!!!
      call process_namelist(last_irr,no_cont_exps,number_of_continuum_shells)

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
      number_of_target_shells = size(CGTO_shell_data)

      n = number_of_target_shells+number_of_continuum_shells
      err = atomic_orbital_basis%init(n,symmetry_inp) !+1 = one extra shell

      do i=1,number_of_target_shells
         call atomic_orbital_basis%add_shell(CGTO_shell_data(i))
      enddo 

      do i=0,max_l

         CGTO_shell_data(1)%l = i
         CGTO_shell_data(1)%number_of_functions = 2*CGTO_shell_data(1)%l+1
         CGTO_shell_data(1)%center(1:3) = (/0.0_cfp,0.0_cfp,0.0_cfp/)
         CGTO_shell_data(1)%number_of_primitives = 1
         call CGTO_shell_data(1)%make_space(CGTO_shell_data(1)%number_of_primitives)

         do j=1,no_cont_exps(i)
            CGTO_shell_data(1)%exponents(1:CGTO_shell_data(1)%number_of_primitives) = exponents(j,i)
            CGTO_shell_data(1)%contractions(1:CGTO_shell_data(1)%number_of_primitives) = 1.0_cfp
            CGTO_shell_data(1)%non_zero_at_boundary = .true.
            call atomic_orbital_basis%add_shell(CGTO_shell_data(1))
         enddo !j

      enddo !i

      !describe where the integrals will be stored 
      err = integral_storage%init(memory=ao_integrals) !(disk='./oneel')
      if (err .ne. 0) then
         print *,err
         call mpi_xermsg('main','main','error initializing the target ao_integral_storage',1,1)
      endif

      integral_options%a = 10.0_cfp !-1.0_cfp !R-matrix radius?
      integral_options%two_p_continuum = .false. !Calculate integrals for 2p in the continuum?
      integral_options%calculate_overlap_ints = .true.
      integral_options%calculate_kinetic_energy_ints = .true.
      integral_options%calculate_nuclear_attraction_ints = .true.
      integral_options%calculate_one_el_hamiltonian_ints = .true.
      integral_options%calculate_property_ints = .true.
      integral_options%max_property_l = 2
      
      call atomic_orbital_basis%one_electron_integrals(integral_storage,integral_options)

      call ao_integrals%print(.true.)

      call atomic_orbital_basis%write('./basis')
      call atomic_orbital_basis%read('./basis')

      call mpi_mod_finalize

contains

   subroutine process_namelist(last_irr,no_cont_exps,number_of_continuum_shells)
      use precisn
!      use const, only: input_unit
      implicit none
      integer, intent(out) :: last_irr,no_cont_exps(0:),number_of_continuum_shells

      integer :: i, j, idat

         open(newunit=idat,file='./inp',form='formatted',status='unknown')
         read(idat,nml=target_data)

         if (molden_file .eq. '') call mpi_xermsg('test','process_namelist','No Molden input file specified.',1,1)
         if (no_sym_op < 0) call mpi_xermsg('test','process_namelist','The number of symmetry operations cannot be smaller than 0.',2,1)

         last_irr = 0
         i = 0
         do
            i = i + 1
            if (i > 8) then 
               last_irr = i-1
               exit
            endif
    
            if (nob(i) .le. 0) then 
               last_irr = i-1
               exit
            endif
         enddo

         if (last_irr .le. 0) call mpi_xermsg('test','process_namelist','Error in input NOB array: the index of the last IRR is .le. 0.',3,1)

         read(idat,nml=continuum_data)

         no_cont_exps(:) = 0
         number_of_continuum_shells = 0

         do i=0,max_l
            j = 0
            do
               j = j + 1
               if (j > max_exp) exit

               if (exponents(j,i) .le. 0.0_cfp) exit
            enddo
            j = j - 1
            
            if (j .le. 0) call mpi_xermsg('test','process_namelist','No exponents given in L shell < max_L.',4,1)
            no_cont_exps(i) = j
            number_of_continuum_shells = number_of_continuum_shells + j
         enddo

         read(idat,nml=process_control)

         if (ao_integrals_file_name == '') call mpi_xermsg('test','process_namelist','The input variable ao_integrals_file_name is empty but must not be empty.',5,1)
         if (mo_integrals_file_name == '') call mpi_xermsg('test','process_namelist','The input variable mo_integrals_file_name is empty but must not be empty.',6,1)
         
         if (ao_integrals_file_name .eq. mo_integrals_file_name) call mpi_xermsg('test','process_namelist','The input variables ao/mo_integrals_file_name must be different.',7,1)
         if (n_mpi_per_numa .le. 0) call mpi_xermsg('test','process_namelist','The number of MPI processes per NUMA region (n_mpi_per_numa) must be a positive number.',8,1)

         close(idat)

   end subroutine process_namelist

end program test_atomic_basis_mod
