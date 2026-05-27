! Copyright 2019
!
! Jan Senk with contributions from others (see the UK-AMOR website)                               
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
   use const_gbl, only: line_len, level1, set_verbosity_level, Fock_blocks_header, orbs_line
   use mpi_gbl
   use atomic_basis_gbl
   use molecular_basis_gbl
   use basis_data_generic_gbl
   use data_file_gbl

   implicit none

   type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
   type(molecular_orbital_basis_obj) :: molecular_orbital_basis
   type(orbital_data_obj) :: orbital_data
   type(data_file_obj) :: data_file

   integer, parameter :: verbosity = 3

   character(len=line_len) :: path, val
   integer :: err, length, irr, i, j, k, m, lunit, first_record, last_record

   real(kind=cfp), allocatable :: target_energies(:), continuum_energies(:), TT_block(:,:), TC_block(:,:), CC_block(:,:)
   integer :: n_target, n_continuum
   logical, allocatable :: is_continuum(:)
      
      !INITIALIZATION OF THE PROGRAM

      call mpi_mod_start

      call set_verbosity_level(verbosity)

      !READING THE PATH TO THE STORAGE FILE

      call get_command_argument(1,path,length,err)
      if (err > 0) then
         call mpi_xermsg('read_fock_blocks', 'main', &
                         'The path to the data file has not been given: it must be given as the first command line argument.', 1, 1)
      end if

      !READING THE ATOMIC AND MOLECULAR BASES DATA

      call atomic_orbital_basis%read(path)

      molecular_orbital_basis%ao_basis => atomic_orbital_basis
      call molecular_orbital_basis%read(path)

      call molecular_orbital_basis%print

      call molecular_orbital_basis%print_energy_sorted_orbital_table

      call molecular_orbital_basis%print_orbitals

      !OPENING THE FILE

      call data_file%open(path)
      lunit = data_file%get_unit_no()
      err = data_file%find_header(Fock_blocks_header,first_record,last_record)
      if (err /= 0) then
         call xermsg ('read_fock_blocks', 'main', &
                      'Error locating the header corresponding to the Fock matrix blocks.', err, 1)
      end if

      !The following may occur if writing of the data record has not been finished or the record is corrupt.
      if (first_record <= 0 .or. last_record <= 0) then
         call xermsg ('basis_data_generic_obj', 'read_bdg', &
                      'The requested data record is missing or not complete.', 2, 1)
      end if
      last_record = first_record

      !READING AND PRINTING THE ORBITAL ENERGIES AND THE FOCK MATRIX BLOCKS

      do irr = 1, molecular_orbital_basis%no_irr
         write(level1,'(/,3X,"Symmetry: ",i2,/)') irr
         call molecular_orbital_basis%get_shell_data(irr,orbital_data)

         read(lunit,pos=last_record,err=10) n_target, n_continuum
         inquire(lunit,pos=last_record)

         do i=1,n_target
            write(level1,'(5X,i0," Target,    energy: ",e25.15)') i, orbital_data%energy(i)
         enddo
         do i=n_target+1,n_target+n_continuum
            write(level1,'(5X,i0," Continuum, energy: ",e25.15)') i, orbital_data%energy(i)
         enddo

         if (n_target > 0) then
            allocate(TT_block(n_target,n_target))
            read(lunit,pos=last_record,err=10) TT_block
            inquire(lunit,pos=last_record)
            
            ! Writing the block into the log file
            write(level3,'(/,2X,"Target-target block of Fock matrix")')
            k = 0
            do i=1,n_target/orbs_line
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
               do j=1,n_target
                  write(level3,'(i0,50e25.15)') j,TT_block(j,k+1:k+orbs_line)
               enddo !j
               k = k + orbs_line
            enddo !i

            m = mod(n_target,orbs_line)
            if (m > 0) then
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
               do j=1,n_target
                  write(level3,'(i0,50e25.15)') j,TT_block(j,k+1:k+m)
               enddo !j
            endif
            ! End of writing the block into the log file

         endif
         if (n_target > 0 .and. n_continuum > 0) then
            allocate(TC_block(n_target,n_continuum))
            read(lunit,pos=last_record,err=10) TC_block
            inquire(lunit,pos=last_record)

            ! Writing the block into the log file
            write(level3,'(/,2X,"Continuum-target block of Fock matrix (a transpose of the saved block)")')
            k = 0
            do i=1,n_target/orbs_line
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+orbs_line)
               do j=1,n_continuum
                  write(level3,'(i0,50e25.15)') j+n_target,TC_block(k+1:k+orbs_line,j)
               enddo !j
               k = k + orbs_line
            enddo !i

            m = mod(n_target,orbs_line)
            if (m > 0) then
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1,k+m)
               do j=1,n_continuum
                  write(level3,'(i0,50e25.15)') j+n_target,TC_block(k+1:k+m,j)
               enddo !j
            endif
            ! End of writing the block into the log file

         endif
         if (n_continuum > 0) then
            allocate(CC_block(n_continuum,n_continuum))
            read(lunit,pos=last_record,err=10) CC_block
            inquire(lunit,pos=last_record)

            ! Writing the block into the log file
            write(level3,'(/,2X,"Continuum-continuum block of Fock matrix")')
            k = 0
            do i=1,n_continuum/orbs_line
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1+n_target,k+orbs_line+n_target)
               do j=1,n_continuum
                  write(level3,'(i0,50e25.15)') j+n_target,CC_block(k+1:k+orbs_line,j)
               enddo !j
               k = k + orbs_line
            enddo !i

            m = mod(n_continuum,orbs_line)
            if (m > 0) then
               write(level3,'(/,10X,50(i0,2X))') (j,j=k+1+n_target,k+m+n_target)
               do j=1,n_continuum
                  write(level3,'(i0,50e25.15)') j+n_target,CC_block(k+1:k+m,j)
               enddo !j
            endif
            ! End of writing the block into the log file

         endif

         deallocate(is_continuum)
         if (allocated(TT_block)) deallocate(TT_block)
         if (allocated(TC_block)) deallocate(TC_block)
         if (allocated(CC_block)) deallocate(CC_block)
      enddo !irr

      !FINALIZE

      err = molecular_orbital_basis%final()
      if (err .ne. 0) call mpi_xermsg('read_fock_blocks','main','Error finalizing molecular_orbital_basis.',3,1)

      call mpi_mod_finalize

10      call xermsg ('read_fock_blocks', 'main', 'Error reading the data from the file and position given.', 4, 1)

end program test
