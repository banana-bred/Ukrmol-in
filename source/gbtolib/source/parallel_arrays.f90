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
!> This module implements the p2d_array_obj which is used throughout the code to store and manipulate the calculated integrals.
module parallel_arrays_gbl
   use file_mapping_gbl, only: file_mapping, map
   use precisn_gbl
   use mpi_gbl
   use mpi_memory_gbl
   use const_gbl, only: stdout, level2, level3, Mib, line_len

   private

   public p2d_array_obj

   !> This is an array object which is used to store the atomic and molecular integrals in memory.
   !> The primary data structure is the 2D array this%a(1:this%d1,1:this%d2). The object may use the auxiliary integer array this%block_offset(1:this%no_blocks)
   !> as in the case of computation of the atomic integrals in the MPI regime. The use of the block_offset array is controled by the initialization routine init.
   !> Each column of this%a has a name which is contained in this%column_descriptor(1:this%d2). The intention is that the different types of integrals (e.g. overlap, kinetic energy, etc.) are kept in
   !> different columns of the array this%a.
   !>
   !> There are multiple possible locations where the main two-dimensional data array is stored, which are governed by features
   !> enabled during compilation of the library and by arguments passed to the `init` routine:
   !>
   !>    1. If a non-empty argument `scratch` is given in the call to `init` and mapping is enabled by `file_mapping_gbl::map`,
   !>       the data array will be stored on disk and accessed by file memory mapping. Upon destruction of this object, the file
   !>       will be deleted.
   !>    2. Otherwise, if MPI-3 shared memory is enabled via `mpi_gbl::shared_enabled`, the array will be allocated only once
   !>       per node. In this case it is the caller's responsibility to assure that all processes on a node expect the same
   !>       layout of the memory.
   !>    3. If neither of the above holds, the array is allocated by each process as specified by the dimensions in the call.
   !>
   !> Alternatively, the routine `init` can be skipped and the reading routine `read` called directly. Then the data will be either
   !> mapped to memory (if `map` is enabled; this time no `scratch` is needed), or shared on a node (if `shared_enabled` is on),
   !> or simply duplicated on each process. There is no parallel distribution of data on reading.
   !>
   type :: p2d_array_obj
      !> Set to .true. following a call to init.
      logical, private :: initialized = .false.
      !> Set to .true. by init if the object has been initialized with no_blocks > 0. This implies the use of the block_offset array.
      logical, private :: indexed = .false.
      !> Number of blocks.
      integer, private :: no_blocks = -1
      !> Local portion of the global array 'a': a_loc(my_start:my_end,1:d2). In the local mode my_start = 1 on all processes.
      !> The array can be shared on a single node in order to reduce memory usage in MPI-Scatci (Ahmed's implementation).
      !> \warning Never change the dimensions of this array by hand: if you need to do it then implement a method for this object otherwise you're making the d1,d2,my_start,etc. values inconsistent.
      !> \warning The type of this array variable cannot be ``real(kind=cfp), allocatable, target" since that is not allowed by Fortran 2003 standard. 
      !> \warning Therefore if MPI 3.0 is to be used this object must be rewritten to work with array a(:,:) as an external input declared with the ``target" attribute or using module-local variables.
      !> \warning I.e. the type-bound array a(:,:) should be replaced by ``real(kind=cfp), pointer :: a(:,:)" and the pointer associated when ``init" or ``read" is called.
      real(kind=cfp), pointer :: a(:,:) => null()
      !> Used in case of shared memory allocation of this%a. It stores the MPI window that this%a is stored in. If it is -1 then we are in the non-shared mode for this%a.
      integer, private :: shared_window = -1
      !> Dimensions of the array 'a'.
      integer, private :: d1 = -1, d2 = -1
      !> Pointers to start of blocks of data within the array 'a'. block_offset has dimension no_blocks.
      integer, allocatable :: block_offset(:)
      !> This array is used to name the data stored in the columns of 'a'.
      character(len=line_len), allocatable, private :: column_descriptor(:)
      !> Helper type for direct virtual memory file mapping.
      type(file_mapping), private :: mmap
   contains
      !> Initializes the object. 
      !> Arguments:
      !> \param[in] n
      !> n is integer: number of rows of this%a
      !> \param[in] d2_end
      !> d2_end is integer: the number of columns of this%a
      !> \param[in] no_blocks:
      !> no_blocks is integer: number of elements of this%block_offset
      !> \param[in] column_descriptor
      !> column_descriptor is a character array of dimension d2_end: contains descriptions of the values stored in the columns of this%a.
      procedure :: init
      !> Deallocates all arrays and restores the object to an uninitialized state. It is a blocking routine.
      procedure :: final
      !> If the array is initialized with no_blocks > 0 this routine can be used to copy the contents of an external integer array to the array this%block_offset which must have the same size.
      procedure :: set_block_offset
      !> Returns the values d1,d2,no_blocks.
      procedure :: get_array_dim
      !> Sorts the array 'a' using keys given as an array of integers of the same dimension as 'a'.
      !> It is a non-blocking routine.
      procedure :: sort_array_with_keys
      !> Prints out the non-zero contents of the array this%a (and this%block_offset, if this%no_blocks > 0). It is a non-blocking routine.
      !> If optional integer parameter out_unit is present, the contents are written into that unit instead of level2.
      procedure :: print
      !> Writes to the disk the arrays this%a,this%block_offset from the process specified. It is a blocking routine.
      procedure :: write
      !> The arrays this%a,this%block_offset are read-in from the disk by the master task and redistributed to other processes. 
      !> If shared-memory MPI is used then each NODE keeps only one copy of the arrays this%a,this%block_offset otherwise the array is kept by every MPI task. 
      !> It is a blocking routine.
      procedure :: read
      !> In-place reduction of the local arrays this%a on all processes using MPI_ALLREDUCE. It is a blocking routine.
      procedure :: reduce_a_local
      !> Reduction of the local arrays this%a on all processes to the array kept by the master. It is a blocking routine.
      procedure :: reduce_a_to_master
      !> Merges (gathers) indexed arrays with equal number of blocks from all tasks into a single indexed array on each task.
      !> The same array this % a is used for sends and receives: this assumes that the block of data of size sendcount sent by the process with rank j is
      !> located in the (j+1)st block of this % a. We assume that the
      !> block_offsets supplied on input by the process j have been adjusted before to point to the rank j's data in the merged array this % a.
      !> I.e. this routine doesn't perform any shifts of the pointers to the data in this % a.
      !> \param[in] sendcount  Number of elements in this % a(:, d2) that this process sends. This number may differ accross processes.
      !> \param[in] d2         Index in the 2nd dimenstion of this % a in which the merge will be done.
      !> \param[in] split_within_block If .true. then it is assumed that the values from each block have been redistributed
      !>                               among the tasks, i.e. if a non-zero value appears in within a block it should be kept by a single task only.
      procedure :: merge_indexed_array_from_all_tasks
      !> Logical function: does the object contain the array this%block_offset? It is a non-blocking routine.
      procedure :: have_offsets
      !> Returns the column index in the array this%a corresponding to the name of the column specified on input.
      procedure :: find_column_matching_name
      !> Returns the character array column_descriptor.
      procedure :: get_column_descriptor
      !> Combines two columns of this%a into a third one using a given operation: '+','-','*','/'.
      procedure :: combine_columns
      !> Multiplies a given column of this%a by a real number.
      procedure :: multiply_column
   end type p2d_array_obj

contains

   function have_offsets(this)
      implicit none
      class(p2d_array_obj) :: this
      logical :: have_offsets

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'have_offsets', 'The object has not been initialized.', 1, 1)
         end if

         have_offsets = this%indexed

   end function have_offsets

   subroutine sort_array_with_keys(this,n,d2,keys,sort_method)
      use sort_gbl, only: heap_sort_int_float, sort_int_float
      use omp_lib
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: n, d2, sort_method
      integer, intent(inout) :: keys(n,d2)

      real(kind=wp) :: start_t, end_t

         start_t = omp_get_wtime()

         write(level3,'("--------->","p2d_array_obj:sort_array_with_keys")')

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'sort_array_with_keys', 'The object has not been initialized.', 1, 1)
         end if

         if (n > this % d1 .or. n <= 0 .or. d2 <= 0 .or. d2 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'sort_array_with_keys', &
                             'On input n and/or d2 were out of allowed limits for the array this%a.', 2, 1)
         end if

         select case (sort_method)
            case (1)
               call heap_sort_int_float(n,keys(:,d2),this%a(:,d2))
            case (2)
               call sort_int_float(n,d2,keys,this%a)
            case default
               call mpi_xermsg ('parallel_arrays', 'sort_array_with_keys', &
                                'On input sort_method was out of the allowed range [1,2].', 3, 1)
         end select

         write(level3,'("<---------","done:p2d_array_obj:sort_array_with_keys")')

         end_t = omp_get_wtime()

         write(level2,'("sort_array_with_keys took [s]",f25.15)') end_t-start_t

   end subroutine sort_array_with_keys

   function init(this,n,d2_end,no_blocks,column_descriptor,scratch_opt)
      implicit none
      integer :: init
      class(p2d_array_obj) :: this
      integer, intent(in) :: n, d2_end,no_blocks
      character(len=line_len), intent(in) :: column_descriptor(:)
      character(len=*), intent(in), optional :: scratch_opt

      integer :: err, i, j
      integer(int64) :: offset, memory
      integer, allocatable :: nelem(:)
      character(len=:), allocatable :: scratch

         init = 0

         if (present(scratch_opt)) then
            if (scratch_opt /= '') then
               scratch = scratch_opt
            end if
         end if

         call check_mpi_running

         call mpi_mod_barrier(err)

         if (this%initialized) then
            err = this%final()
            if (err .ne. 0) call mpi_xermsg('parallel_arrays','init','Error finalizing the object.',1,1)
         endif

         if (d2_end .le. 0 .or. n .le. 0) then
            print *,d2_end,n
            call mpi_xermsg('parallel_arrays','init','d2_end .le. 0 .or. n .le. 0 but both must be .ge. 1.',2,1)
         endif

         if (size(column_descriptor) .ne. d2_end) then
            print *,size(column_descriptor),d2_end
            call mpi_xermsg ('parallel_arrays', 'init', &
                             'size of the input array column_descriptor .ne. d2_end but they must be equal.', 3, 1)
         endif

         if (no_blocks > 0) then
            write(level3,'("--------->","p2d_array_obj:init; indexed array")')
            this%indexed = .true.
         else
            write(level3,'("--------->","p2d_array_obj:init; unindexed array")')
            this%indexed = .false.
         endif

         !Transfer the column labels
         allocate(this%column_descriptor,source=column_descriptor,stat=err)
         if (err .ne. 0) call mpi_xermsg('parallel_arrays','init','Memory allocation 0 failed.',err,1)

         allocate(nelem(1:nprocs),stat=err)
         if (err .ne. 0) call mpi_xermsg('parallel_arrays','init','Memory allocation 1 failed.',err,1)

         !Check that all processes request the same type of array (shared or local)
         nelem(:) = 0
         call mpi_mod_allgather(no_blocks,nelem)

         do i=1,nprocs
            do j=1,i
               if (nelem(i) > 0 .and. nelem(j) <= 0) then
                    call mpi_xermsg ('parallel_arrays', 'init', &
                                     'At least one process requests shared array while another one requests local array.', 4, 1)
               end if
            enddo
         enddo

         !Check that the d2_end dimension is the same on all processes: choose arbitrarily the root's d2_start value as the reference.
         nelem(:) = 0
         call mpi_mod_allgather(d2_end,nelem)

         if (nelem(myrank+1) /= nelem(1)) then
            call mpi_xermsg ('parallel_arrays', 'init', 'd2_start value must be the same for all processes.', 5, 1)
         end if

         !Send the size of the local array 'a' to everyone so that each process knows exactly how many elements every other process keeps.
         nelem(:) = 0
         call mpi_mod_allgather(n,nelem)

         call mpi_mod_barrier(err) !putting the barrier here ensures that we catch any errors in the dimensioning before the allocations take place.

         this%d1 = n
         this%d2 = d2_end
         this%no_blocks = no_blocks

         write(level2,'(10X,"rank ",i0,", array dimensions: ",i0,1X,i10)') myrank, n, d2_end
         if (this%indexed) then
            write(level2,'(10X,"rank ",i0,", number of blocks: ",i0)') myrank, no_blocks
         endif

         !count the number of bytes
         offset = 0
         memory = n*d2_end*cfp_bytes
         write(level2,'("Memory to be allocated for the array [MiB]: ",f12.3)') memory/(Mib*1.0_cfp)

         ! map array from a scratch disk file or allocate in memory
         if (map .and. allocated(scratch)) then
            call this % mmap % init(scratch, offset, memory, create_new = .true., unique_name = .true., read_write = .true.)
            call c_f_pointer(this%mmap%ptr, this%a, [this%d1, this%d2])
         else
            this%shared_window = mpi_memory_allocate_real_2dim(this%a, this%d1, this%d2)
         end if

         this%a = 0.0_cfp

         if (this%indexed) then
            memory = no_blocks*bit_size(memory)/8
            write(level2,'("Memory to be allocated for the block offset [MiB]: ",f12.3)') memory/(Mib*1.0_cfp)
            allocate(this%block_offset(this%no_blocks),stat=err)
            this%block_offset = -1
         endif

         this%initialized = .true.

         call mpi_mod_barrier(err)

         write(level3,'("<---------","done:p2d_array_obj:init")')

   end function init

   function final(this)
      use iso_c_binding, only: c_associated
      implicit none
      integer :: final
      class(p2d_array_obj) :: this

      integer :: ierr

         final = 0

         call check_mpi_running

         call mpi_mod_barrier(ierr)

         write(level3,'("--------->","p2d_array_obj:final")')

         ierr = 0
         if (allocated(this%block_offset)) deallocate(this%block_offset,stat=ierr)
         if (ierr .ne. 0) call mpi_xermsg('parallel_arrays','final','Memory deallocation 3 failed.',ierr,1)
         ierr = 0
         if (allocated(this%column_descriptor)) deallocate(this%column_descriptor,stat=ierr)
         if (ierr .ne. 0) call mpi_xermsg('parallel_arrays','final','Memory deallocation 4 failed.',ierr,1)

         if (map .and. c_associated(this % mmap % ptr)) then
            call this % mmap % finalize
         else if (associated(this%a)) then
            call mpi_memory_deallocate_real_2dim(this%a,size(this%a),this%shared_window)
         endif

         nullify(this % a)

         this%initialized = .false.
         this%indexed = .false.
         this%d1 = -1
         this%d2 = -1
         this%no_blocks = -1
         this%shared_window = -1

         call mpi_mod_barrier(ierr)

         write(level3,'("<---------","done:p2d_array_obj:final")')

   end function final

   subroutine set_block_offset(this,block_offset)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: block_offset(:)

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'set_block_offset', 'The object has not been initialized.', 1, 1)
         end if

         if (this % no_blocks /= size(block_offset)) then
            call mpi_xermsg ('parallel_arrays', 'set_block_offset', &
                             'On input size of block_offset does not match the size for which the object was initialized.', 2, 1)
         end if

         this%block_offset(1:this%no_blocks) = block_offset(1:this%no_blocks)

   end subroutine set_block_offset

   subroutine get_array_dim(this,d1,d2,no_blocks)
      use omp_lib
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(out) :: d1,d2,no_blocks

         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','get_array_dim','The object has not been initialized.',1,1)
         

         d1 = this%d1
         d2 = this%d2
         no_blocks = this%no_blocks

   end subroutine get_array_dim

   subroutine print(this,only_non_zero,out_unit)
      implicit none
      class(p2d_array_obj) :: this
      logical, intent(in) :: only_non_zero

      integer :: i, j, u

      integer, intent(in), optional :: out_unit


         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','print','The object has not been initialized.',1,1)

         write(level3,'("--------->","p2d_array_obj:print")')

         if (present(out_unit)) then
            u = out_unit
            write(level3,'("Printing redirected into the file under unit: ",i0)') out_unit
         else
            u = level2
         endif

         do j=1,this%d2
            if (this%indexed) then
               write(u,'("Number of block offsets: ",2i10)') j, this%no_blocks
               if (only_non_zero) then
                  do i=1,this%no_blocks
                     if (this%block_offset(i) >= 0) write(u,'(i10,1X,i20)') i,this%block_offset(i)
                  enddo
               else
                  do i=1,this%no_blocks
                     write(u,'(i10,1X,i20)') i,this%block_offset(i)
                  enddo
               endif
            endif

            write(u,'("Column ",i5," descriptor: ",a)') j, this%column_descriptor(j)
            if (only_non_zero) then
               do i=1,this%d1
                  if (this%a(i,j) .ne. 0.0_cfp) write(u,'(i0,1X,e25.15)') i,this%a(i,j)
               enddo !i
            else
               do i=1,this%d1
                  write(u,'(i0,1x,e25.15)') i,this%a(i,j)
               enddo !i
            endif
         enddo !j

         write(level3,'("<---------","done:p2d_array_obj:print")')

   end subroutine print

   function find_column_matching_name(this,column_name)
      implicit none
      class(p2d_array_obj) :: this
      character(len=*), intent(in) :: column_name
      integer :: find_column_matching_name

      integer :: i
      integer :: a_size
      logical :: found

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'find_column_matching_name', 'The object has not been initialized.', 1, 1)
         end if
         a_size = size(this%a,2)
         found = .false.
         do i=1,a_size
            !write (stdout, '(5A)') '"', trim(this%column_descriptor(i)), '" "', trim(column_name), '"'
            if (trim(this%column_descriptor(i)) .eq. trim(column_name)) then
               find_column_matching_name = i
               found = .true.
               exit
            endif
         enddo !i

         if (.not.(found)) then
            write (stdout, '(A)') 'required column:'
            write (stdout, '(5x,3A)') '"', trim(column_name), '"'
            write (stdout, '(A)') 'available columns:'
            do i = 1, a_size
               write (stdout, '(5x,3A)') '"', trim(this%column_descriptor(i)), '"'
            end do
            call mpi_xermsg ('parallel_arrays', 'find_column_matching_name', &
                             'A column matching the input name has not been found in the array.', 2, 1)
         endif

   end function find_column_matching_name

   subroutine get_column_descriptor(this,column_descriptor)
      implicit none
      class(p2d_array_obj) :: this
      character(len=line_len), allocatable :: column_descriptor(:)

      integer :: err

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'get_column_descriptor', 'The object has not been initialized.', 1, 1)
         end if

         if (allocated(column_descriptor)) deallocate(column_descriptor)

         allocate(column_descriptor,source=this%column_descriptor,stat=err)
         if (err .ne. 0) call mpi_xermsg('parallel_arrays','get_column_descriptor','Memory allocation failed.',err,1)

   end subroutine get_column_descriptor

   !> \todo the 'rank' parameter should be allowed to be < 0 in which case all processes write equal portion of the integrals array thereby speeding up the write on Lustre system.
   subroutine write(this,lunit,record_start,position_after_write,rank)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: lunit, rank
      integer, intent(in) :: record_start
      integer, intent(out) :: position_after_write

      integer :: err

      real :: start_t, end_t

         call cpu_time(start_t)

         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','write','The object has not been initialized.',1,1)

         write(level3,'("--------->","p2d_array_obj:write")')

         call mpi_mod_barrier(err)

         if (rank < 0 .or. rank > nprocs-1) call mpi_xermsg('parallel_arrays','write','On input rank is out of range.',3,1)

         if (myrank .eq. rank) then
            write(lunit,pos=record_start,err=10) size(this%a,1),size(this%a,2)
            write(lunit) this%column_descriptor
            write(lunit,err=10) this%a

            write(lunit) this%no_blocks
            if (this%no_blocks > 0) then
               write(lunit) this%block_offset
            endif
            inquire(lunit,pos=position_after_write)
         endif

         !get the position_after_write on all processes
         call mpi_mod_bcast(position_after_write,int(rank,kind=mpiint))

         write(level3,'("<---------","done:p2d_array_obj:write")')

         call cpu_time(end_t)
         
         write(level2,'("write took [s] ",f25.15)') end_t-start_t

         return

 10      call mpi_xermsg ('parallel_arrays', 'write', &
                          'Error executing the write command while writing the array data to the disk.', 4, 1)

   end subroutine write

   subroutine read(this,file_name,lunit,record_start,position_after_read)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: lunit
      integer, intent(in) :: record_start
      character(len=*), intent(in) :: file_name
      integer, intent(out) :: position_after_read

      integer :: err, d2, i, j
      integer(longint) :: offset, length
      integer(mpiint)  :: local_master_array(nprocs)
      real :: start_t, end_t

         call cpu_time(start_t)
         write(level3,'("--------->","p2d_array_obj:read")')

         call mpi_mod_barrier(err)

         if (this%initialized) then
            err = this%final()
         endif

         if (myrank .eq. master) then
            read(lunit,pos=record_start,err=50) this%d1,this%d2
         endif
         !Everyone needs this before hand
         call mpi_mod_bcast(this%d1,master)
         call mpi_mod_bcast(this%d2,master)
         length = this%d1 * this%d2 * cfp_bytes

         !Allocate the shared memory   
         if (.not. map) this%shared_window = mpi_memory_allocate_real_2dim(this%a,this%d1,this%d2)
         !now we get process one in the local node to do this
         !Fence off the memory for it

         call mpi_memory_synchronize(this%shared_window)

         if (myrank .eq. master) then    
            allocate(this%column_descriptor(this%d2),stat=err)
            if (err .ne. 0) then
               print *,this%d1,this%d2
               call mpi_xermsg('parallel_arrays','read','Master memory allocation failed.',err,1)
            endif

            read(lunit) this%column_descriptor
            inquire(unit=lunit,pos=position_after_read)
            offset = position_after_read - 1
            if (.not. map) read(lunit) this%a
            position_after_read = position_after_read + length
            read(lunit, pos=position_after_read) this%no_blocks

            if (this%no_blocks > 0) then
               allocate(this%block_offset(this%no_blocks),stat=err)
               if (err .ne. 0) call mpi_xermsg('parallel_arrays','read','Memory allocation 1 failed.',err,1)
               read(lunit) this%block_offset
            endif

            inquire(unit=lunit,pos=position_after_read)
         endif

         !Update every other process with the memory
         call mpi_memory_synchronize(this%shared_window)

         !Master must ensure every process knows what the array's dimensions are.
         call mpi_mod_bcast(this%no_blocks,master)

         !Master must ensure every process knows where in the file the data record ends.
         call mpi_mod_bcast(position_after_read,master)

         !If memory mapping is enabled, let each process map the file independently
         if (map) then
            call mpi_mod_bcast(offset, master)
            call this % mmap % init(file_name, offset, length)
            call c_f_pointer(this % mmap % ptr, this % a, [this%d1, this%d2])
         end if

         !Allocate everthing else except for a
         if (myrank .ne. master) then
            allocate(this%column_descriptor(this%d2),stat=err)
            if (err .ne. 0) call mpi_xermsg('parallel_arrays','read','Memory allocation 2 failed.',err,1)
            if (this%no_blocks > 0) then
               allocate(this%block_offset(this%no_blocks),stat=err)
               if (err .ne. 0) call mpi_xermsg('parallel_arrays','read','Memory allocation 3 failed.',err,1)
            endif
         endif

         call mpi_mod_bcast(this%column_descriptor,master)

         if (shared_enabled) then !Every node has only one copy of the array this%a

            call mpi_memory_synchronize(this%shared_window)
            if (local_rank == master) then
               do j=1,this%d2
                  if (.not. map) then
                     call mpi_mod_bcast(this % a(1:this%d1, j), master, local_masters_communicator)
                  end if
                  if (this%no_blocks > 0) then
                     call mpi_mod_bcast(this % block_offset(1:this%no_blocks), master, local_masters_communicator)
                  end if
               end do
            end if
            call mpi_memory_synchronize(this%shared_window)

         else !every MPI task keeps a copy of the arrays this%a,this%block_offset

            if (this%no_blocks > 0) call mpi_mod_bcast(this%block_offset,master)

            !Broadcast the array one dimension at a time - this can be improved
            !to read the full 2d array at once but for most applications what
            !follows should be enough.
            do d2=1,this%d2
               if (.not. map) call mpi_mod_bcast(this%a(1:this%d1,d2),master)
            enddo

         endif

         this%initialized = .true.

         call mpi_mod_barrier(err)

         write(level3,'("<---------","done:p2d_array_obj:read")')

         call cpu_time(end_t)

         write(level2,'("read took [s] ",f25.15)') end_t-start_t

         return
   
 50      call mpi_xermsg ('parallel_arrays', 'read', &
                          'Error executing the read command while reading the array data from the disk.', 3, 1)
   
   end subroutine read
   
   subroutine reduce_a_local(this)
      implicit none
      class(p2d_array_obj) :: this

      integer :: err, i, nelem
      real :: start_t, end_t

         call cpu_time(start_t)

         if (.not.(this%initialized)) call mpi_xermsg('parallel_arrays','reduce_a_local','The object has not been initialized.',1,1)

         write(level3,'("--------->","p2d_array_obj:reduce_a_local")')

         call mpi_mod_barrier(err)

         !For the reduction to work all processes must keep the same number of elements of the array to be reduced. The index range for 'a' is guaranteed to be the same if the number of elements 
         !agree since the local array 'a' on each process starts with index 1. To be absolutely sure that the dimension is the same on all processes we arbitrarily select the master as the reference
         !process and use its value of size(this%a,1) to check that all processes's local arrays 'a' have the same dimension.
         nelem = size(this%a,1)
         call mpi_mod_bcast(nelem,master)
         
         if (nelem .ne. size(this%a,1)) then
            print *,myrank,nelem,size(this%a,1)
            call mpi_xermsg ('parallel_arrays', 'reduce_a_local', &
                             'The local dimensions of the array a to be reduced must be the same on all processes but isnt.', 3, 1)
         endif

         write(level2,'(10X,"Reducing (summing) ",i0," rows of the local array a accross all processes...")') nelem
       
         !This can be rewritten so that the reduction is on the whole 2D array but for most applications this should be good enough. 
         do i=1,size(this%a,2)
            call mpi_reduceall_inplace_sum_cfp(this%a(1:size(this%a,1),i), nelem)
         enddo

         write(level2,'(10X,"The local array on all processes has been replaced by its reduced version &
                            &obtained by summing contributions from all processes.")')

         write(level3,'("--------->","done:p2d_array_obj:reduce_a_local")')

         call cpu_time(end_t)

         write(level2,'("reduce_a_local took [s] ",f25.15)') end_t-start_t

   end subroutine reduce_a_local

   subroutine reduce_a_to_master(this)
      implicit none
      class(p2d_array_obj) :: this

      integer :: err, i, nelem
      real :: start_t, end_t

      integer, parameter :: tag = 1

         call cpu_time(start_t)

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'reduce_a_to_master', 'The object has not been initialized.', 1, 1)
         end if

         write(level3,'("--------->","p2d_array_obj:reduce_a_to_master")')

         call mpi_mod_barrier(err)

         !For the reduction to work all processes must keep the same number of elements of the array to be reduced. The index range for 'a' is guaranteed to be the same if the number of elements 
         !agree since the local array 'a' on each process starts with index 1. To be absolutely sure that the dimension is the same on all processes we arbitrarily select the master as the reference
         !process and use its value of size(this%a,1) to check that all processes's local arrays 'a' have the same dimension.
         nelem = size(this%a,1)
         call mpi_mod_bcast(nelem,master)
         
         if (nelem .ne. size(this%a,1)) then
            print *,myrank,nelem,size(this%a,1)
            call mpi_xermsg ('parallel_arrays', 'reduce_a_to_master', &
                             'The local dimensions of the array a to be reduced must be the same on all processes but isnt.', 3, 1)
         endif

         write(level2,'(10X,"Master is reducing (summing) ",i0," rows of the local array a using &
                            &contributions from all processes...")') nelem

         !This can be rewritten so that the reduction is on the whole 2D array but for most applications this should be good enough. 
         do i=1,size(this%a,2)
            call mpi_reduce_inplace_sum_cfp(this%a(1:size(this%a,1),i), nelem)
         enddo

         write(level2,'(10X,"The local array on the Master process has been replaced by its reduced version obtained &
                            &by summing contributions from all processes.")')

         write(level3,'("--------->","done:p2d_array_obj:reduce_a_to_master")')

         call cpu_time(end_t)

         write(level2,'("reduce_a_to_master took [s] ",f25.15)') end_t-start_t

   end subroutine reduce_a_to_master

   subroutine merge_indexed_array_from_all_tasks(this, sendcount, d2, split_within_block)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: sendcount, d2
      logical, intent(in) :: split_within_block

      integer :: err, myrank_offset, i, j, k, p, offset, jth_offset, displs(nprocs), block_size, ii, jj
      integer, allocatable :: block_offset_all(:)
      real :: start_t, end_t
      
         write(level3,'("--------->","p2d_array_obj:merge_indexed_array_from_all_tasks")')

         call cpu_time(start_t)

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'merge_indexed_array_from_all_tasks', 'The object has not been initialized.', 1, 1)
         end if

         if (.not.this % indexed) then
            call mpi_xermsg ('parallel_arrays', 'merge_indexed_array_from_all_tasks', &
                             'The array to be merged must be initialized as indexed.', 2, 1)
         endif

         call mpi_mod_allgather(sendcount, displs)

         if (sendcount < 0 .or. sum(displs) > this % d1) then
            print *,sendcount, sum(displs), this % d1
            call mpi_xermsg ('parallel_arrays', 'merge_indexed_array_from_all_tasks', &
                             'On input sendcount had invalid size.', 3, 1)
         endif

         if (d2 /= this % d2) then
            print *, d2, this % d2
            call mpi_xermsg ('parallel_arrays', 'merge_indexed_array_from_all_tasks', &
                             'On input d2 was incompatible with this % d2.', 5, 1)
         endif

         ! this array may be large
         allocate(block_offset_all(this % no_blocks * nprocs), stat=err)
         if (err /= 0) then
            call mpi_xermsg ('parallel_arrays', 'merge_indexed_array_from_all_tasks', &
                             'Memory allocation failed.', err, 1)
         endif
         block_offset_all = -1

         !The block offsets received from the jth rank will be placed in the jth section of block_offset_all.
         !Each section has length this % no_blocks so offset for block_offset from the rank j is: j * this % no_blocks
         displs(1) = 0
         do i=2,nprocs
            displs(i) = displs(i-1) + this % no_blocks
         enddo

         !A copy of the block_offset for this task into the final array is needed for IN PLACE gather
         myrank_offset = displs(myrank+1)
         block_offset_all(myrank_offset+1 : myrank_offset+this % no_blocks) = this % block_offset(1 : this % no_blocks)

         call mpi_mod_allgatherv_in_place(this % no_blocks, block_offset_all)

         write(level2,'("Allgather of block_offset completed")')
         flush(stdout)

         !We assume that the input data for each task is already in the
         !correct area where that process would receive its own contribution to
         !the receive buffer, i.e. in analogy to the case of block_offset_all above only in this case the sendcount may
         !be a different number for each task.
         call mpi_mod_allgatherv_in_place(sendcount, this % a(:, d2))

         write(level2,'("Allgather of this%a completed")')
         flush(stdout)

         this % block_offset = -1

         if (split_within_block) then

            ! In this case each block has the same size block_size and the values within the block are assumed to be split between all tasks
            ! so we have to combine the non-zero values from each block from all tasks.
            do j = 2, nprocs
               p = 0
               do i = 1, this % no_blocks

                  jth_offset = (j-1) * this % no_blocks

                  ! pointer to start of i-th block of data on the j-th task
                  offset = block_offset_all(jth_offset + i)

                  if (offset /= -1) then 

                     ! find out where the next block starts and its size
                     ii = i
                     jj = j
                     do
                        ii = ii + 1
                        if (ii > this % no_blocks) then
                           if (jj == nprocs) then
                              block_size = size(this%a, 1) - offset
                              exit
                           endif
                           ! move onto the blocks for the next task
                           jj = jj + 1
                           jth_offset = (jj-1) * this % no_blocks
                           ii = 1
                        endif
                        if (block_offset_all(jth_offset + ii) > 0) then
                           block_size = block_offset_all(jth_offset + ii) - offset
                           exit
                        endif
                     enddo

                     if (block_size < 0) then
                        call mpi_xermsg ('parallel_arrays', 'merge_indexed_array_from_all_tasks', &
                                         'Error determining block_size.', 6, 1)
                     endif

                     do k = 1, block_size
                        ! index of the k-th value in i-th block.
                        p = p + 1

                        if (this % a(offset + k, d2) /= 0.0_cfp) then
                           ! check that the data are correctly redistributed: only one task should keep data for the same value
                           if (this % a(p, d2) /= 0.0_cfp) then
                              write(stdout,'(i10,2e25.15)') p, this % a(p, d2),this % a(offset + k, d2)
                              call mpi_xermsg ('parallel_arrays', 'merge_indexed_array_from_all_tasks', &
                                               'Attempt to overwrite a non-zero value in the merged array.', 7, 1)
                           endif
                           this % a(p, d2) = this % a(offset + k, d2) 
                        endif

                     enddo
                  endif

               enddo
            enddo

            ! copy the block offsets from that for the first process
            do i = 1, this%no_blocks
               this % block_offset(i) = block_offset_all(i)
            enddo

         else
            !Merge the block offsets from all tasks: we assume that the block offsets from this % block_offset have been supplied
            !with the correct shifts for the global array so we don't adjust the actual values of the block_offset pointers.
            do j = 1, nprocs
               jth_offset = (j-1) * this % no_blocks
   
               do i = 1, this % no_blocks
                  if (block_offset_all(jth_offset + i) /= -1) then
                     this % block_offset(i) =  block_offset_all(jth_offset + i)
                  endif
               enddo
   
            enddo !j
         endif

         deallocate(block_offset_all)

         call cpu_time(end_t)

         write(level2,'("took: ",f25.15)') end_t-start_t

         write(level3,'("<---------","p2d_array_obj:merge_indexed_array_from_all_tasks")')

   end subroutine merge_indexed_array_from_all_tasks

   subroutine combine_columns(this,col1,op,col2,col3)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: col1,col2,col3
      character(len=1), intent(in) :: op

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'The object has not been initialized.', 1, 1)
         end if

         if (col1 <= 0 .or. col1 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'On input col1 was out of range.', 2, 1)
         end if
         if (col2 <= 0 .or. col2 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'On input col2 was out of range.', 3, 1)
         end if
         if (col3 <= 0 .or. col3 > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'combine_columns', 'On input col3 was out of range.', 4, 1)
         end if

         if (shared_enabled) then

            call mpi_memory_synchronize(this%shared_window)
            if (local_rank == local_master) then
               select case (op)
                  case ('+')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) + this%a(1:this%d1,col2)
                  case ('-')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) - this%a(1:this%d1,col2)
                  case ('*')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) * this%a(1:this%d1,col2)
                  case ('/')
                     this%a(1:this%d1,col3) = this%a(1:this%d1,col1) / this%a(1:this%d1,col2)
                  case default
                     call mpi_xermsg('parallel_arrays','combine_columns','On input op was not one of: +,-,*,/.',5,1)
               end select
            endif
            call mpi_memory_synchronize(this%shared_window)

          else !each MPI task has its own copy of this%a:

            select case (op)
               case ('+')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) + this%a(1:this%d1,col2)
               case ('-')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) - this%a(1:this%d1,col2)
               case ('*')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) * this%a(1:this%d1,col2)
               case ('/')
                  this%a(1:this%d1,col3) = this%a(1:this%d1,col1) / this%a(1:this%d1,col2)
               case default
                  call mpi_xermsg('parallel_arrays','combine_columns','On input op was not one of: +,-,*,/.',6,1)
            end select

         endif
   end subroutine combine_columns

   subroutine multiply_column(this,col,fac)
      implicit none
      class(p2d_array_obj) :: this
      integer, intent(in) :: col
      real(kind=cfp), intent(in) :: fac

         if (.not. this % initialized) then
            call mpi_xermsg ('parallel_arrays', 'multiply_column', 'The object has not been initialized.', 1, 1)
         end if
         if (col <= 0 .or. col > this % d2) then
            call mpi_xermsg ('parallel_arrays', 'multiply_column', 'On input col was out of range.', 2, 1)
         end if

         if (shared_enabled) then

            call mpi_memory_synchronize(this%shared_window)

            if (local_rank == local_master) then
               this%a(1:this%d1,col) = this%a(1:this%d1,col)*fac
            endif

            call mpi_memory_synchronize(this%shared_window)
         else
            this%a(1:this%d1,col) = this%a(1:this%d1,col)*fac
         endif

   end subroutine multiply_column

end module parallel_arrays_gbl
