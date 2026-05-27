! Copyright 2025
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

!> \brief   Storage for blocks in sparse transformation
!> \authors J Benda
!> \date    2025
!>
!> Auxiliary module for use in the sparse AO-to-MO transformation of two-electron integrals.
!> Contains particullarly the type "sparse_block_storage_obj", which is responsible for
!> abstraction of access to memory- or disk-stored integrals in the sparse transformation.
!>
module sparse_block_storage_gbl

   use file_mapping_gbl, only: file_mapping
   use iso_c_binding,    only: c_f_pointer
   use iso_fortran_env,  only: int64
   use mpi_gbl,          only: myrank
   use precisn_gbl,      only: cfp

   implicit none

   private

   public sparse_block_storage_obj

   !> \brief   Sparse two-electron integral block
   !> \authors J Benda
   !> \date    2025
   !>
   !> A simple class that holds a single block of two electron integrals in one-dimensional
   !> arrays. The array `Rv` holds the integral values, the array `Ri` holds the combined
   !> integral multi-index.
   !>
   type :: sparse_block_obj
      integer,   allocatable :: Ri(:)
      real(cfp), allocatable :: Rv(:)
   end type sparse_block_obj

   !> \brief   Storage of integral blocks
   !> \authors J Benda
   !> \date    2025
   !>
   !> A class that offers unified way of accessing sparse integral blocks from both a disk-based
   !> storage as well as from memory.
   !>
   type :: sparse_block_storage_obj
      logical                             :: in_memory = .true.
      integer                             :: Rn = 0, Ri_fd = 0, Rv_fd = 0
      type(file_mapping)                  :: Ri_mapping, Rv_mapping
      type(sparse_block_obj), allocatable :: blocks(:)
      integer, pointer                    :: Ri(:, :) => null()
      real(cfp), pointer                  :: Rv(:, :) => null()
   contains
      procedure :: reset
      procedure :: add_block
      procedure :: combine_blocks
      final     :: destruct
   end type sparse_block_storage_obj

contains

   !> \brief   Initialize block storage
   !> \authors J Benda
   !> \date    2025
   !>
   !> If non-empty `scratch` is given, this subroutine will create two binary files for the arrays `Ri` and `Rv`.
   !> Otherwise, a memory-only slot array is set up for later use.
   !>
   subroutine reset(this, nblocks, scratch)

      class(sparse_block_storage_obj), intent(inout) :: this

      integer,          intent(in) :: nblocks
      character(len=*), intent(in) :: scratch
      character(len=100)           :: rank

      this%Rn = 0
      this%in_memory = len(scratch) == 0

      if (this%in_memory) then
         if (allocated(this%blocks)) deallocate (this%blocks)
         allocate (this%blocks(nblocks))
      else
         write (rank, '(i0)') myrank
         open (newunit = this%Ri_fd, file = trim(scratch)//'Ri'//rank, access = 'stream', action = 'write', status = 'replace')
         open (newunit = this%Rv_fd, file = trim(scratch)//'Rv'//rank, access = 'stream', action = 'write', status = 'replace')
      end if

   end subroutine reset


   !> \brief   Store a single integral block
   !> \authors J Benda
   !> \date    2025
   !>
   !> Copy data from the provided array to the integral storage (to memory if no scratch files are used, or to disk).
   !> The writing to the files is protected by a OpenMP critical section. This might be a huge bottleneck, because in theory
   !> multiple threads should be able to write to the file at once. However, for that we would need to carefully fiddle with
   !> the file pointer. This is left for a future optimization.
   !>
   subroutine add_block(this, iblk, nnz, Wi, Wv)

      class(sparse_block_storage_obj), intent(inout) :: this

      integer,   intent(in) :: iblk, nnz, Wi(:)
      real(cfp), intent(in) :: Wv(:)

      !$omp atomic update
      this%Rn = this%Rn + nnz

      if (this%in_memory) then
         this%blocks(iblk)%Ri = Wi
         this%blocks(iblk)%Rv = Wv
      else
         !$omp critical(sparse_storage_add_block)
         write (this%Ri_fd) Wi
         write (this%Rv_fd) Wv
         !$omp end critical(sparse_storage_add_block)
      end if

   end subroutine add_block


   !> \brief   Merge blocks and update pointers
   !> \authors J Benda
   !> \date    2025
   !>
   !> Update the pointers this%Ri and this%Rv to point to a storage with consecutive data from the individual integral blocks
   !> provided before. This is trivial in case of the disk storage, where the pointers are simply set up to point to the
   !> memory mapping of the scratch files. In case of operation without the scratch files, an annonymous mapping is created
   !> (for consistent, unified treatment), which is not backed up by any file, and filled with the block data stored in memory.
   !>
   subroutine combine_blocks(this)

      class(sparse_block_storage_obj), intent(inout) :: this

      integer(int64)      :: offset, Ri_length, Rv_length, iblk, nnz
      character(len=1000) :: Ri_filename, Rv_filename

      offset = 0
      Ri_length = this%Rn * storage_size(1, int64)/8
      Rv_length = this%Rn * storage_size(1._cfp, int64)/8

      ! clean up from a possible previous use
      call this % Ri_mapping % finalize
      call this % Rv_mapping % finalize

      ! set up the mapping
      if (this%in_memory) then
         ! create anonymous (in-memory) mapping for unified treatment
         call this % Ri_mapping % init('', offset, Ri_length, .true.)
         call this % Rv_mapping % init('', offset, Rv_length, .true.)
      else
         ! recover names of the scratch files
         inquire (this%Ri_fd, name = Ri_filename)
         inquire (this%Rv_fd, name = Rv_filename)

         ! map the existing pair of scratch files from disk
         call this % Ri_mapping % init(Ri_filename, offset, Ri_length, .true.)
         call this % Rv_mapping % init(Rv_filename, offset, Rv_length, .true.)

         ! unlink the files (as of now they will persist only internally, through the mapping)
         close (this%Ri_fd, status = 'delete');  this%Ri_fd = 0
         close (this%Rv_fd, status = 'delete');  this%Rv_fd = 0
      end if

      ! set up Fortran pointers to the file mapping
      call c_f_pointer(this%Ri_mapping%ptr, this%Ri, [this%Rn, 1])
      call c_f_pointer(this%Rv_mapping%ptr, this%Rv, [this%Rn, 1])

      if (this%in_memory) then
         ! merge the blocks
         do iblk = 1, size(this%blocks)
            if (allocated(this%blocks(iblk)%Ri)) then
               nnz = size(this%blocks(iblk)%Ri)
               this%Ri(offset + 1 : offset + nnz, 1) = this%blocks(iblk)%Ri
               this%Rv(offset + 1 : offset + nnz, 1) = this%blocks(iblk)%Rv
               deallocate (this%blocks(iblk)%Ri)
               deallocate (this%blocks(iblk)%Rv)
               offset = offset + nnz
            end if
         end do
      else
         ! nothing to do for scratch storage - already combined by construction (albeit in unpredictable order)
      end if

   end subroutine combine_blocks


   !> \brief   Finalize the storage object
   !> \authors J Benda
   !> \date    2025
   !>
   !> Remove any remaining scratch files. The rest will be destructed automatically.
   !>
   subroutine destruct(this)

      type(sparse_block_storage_obj), intent(inout) :: this

      if (this%Ri_fd /= 0) close (this%Ri_fd, status = 'delete')
      if (this%Rv_fd /= 0) close (this%Rv_fd, status = 'delete')

      this%Ri_fd = 0
      this%Rv_fd = 0

   end subroutine destruct

end module sparse_block_storage_gbl
