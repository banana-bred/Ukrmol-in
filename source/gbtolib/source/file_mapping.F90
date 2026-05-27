! Copyright 2020
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
!> \brief   Virtual memory file mapping
!> \author  J Benda
!> \date    2020 - 2021
!>
!> Direct mapping of disk files into virtual memory. This allows direct access to data in files without actually
!> reading them to memory. The system performs all necessary I/O on the fly, to obtain the required values at the
!> requested position. For simplicity, this module is accompanied with a C source code, which calls the necessary
!> system API.
!>
module file_mapping_gbl

    use iso_c_binding,   only: c_associated, c_int, c_int64_t, c_loc, c_null_char, c_null_ptr, c_ptr, c_size_t, c_char
    use iso_fortran_env, only: int64

    implicit none

    !> Map integrals to virtual memory instead of fully reading them to the physical memory.
#ifdef usemapping
    logical, parameter :: map = .true.
#else
    logical, parameter :: map = .false.
#endif

    !> An auxiliary derived type which holds the pointer to the mapped data file window
    !> and performs the mapping and unmapping on call to `init` or `final`.
    type file_mapping
        type(c_ptr) :: addr = c_null_ptr    !< Internal page-aligned pointer needed for finalization
        type(c_ptr) :: ptr  = c_null_ptr    !< Pointer pointing to the requested location in file
        integer(c_size_t) :: length = 0     !< Length of the mapped window
    contains
        procedure :: init => init_file_mapping
        procedure :: finalize => finalize_file_mapping
        final     :: destruct_file_mapping
    end type file_mapping

    interface
        subroutine map_file_window (filename, offset, length, addr, ptr, rw, scratch, create) bind(C, name='map_file_window')
            import :: c_int, c_ptr, c_size_t, c_char
            type(c_ptr)              :: addr, ptr
            character(kind=c_char)   :: filename(*)
            integer(c_size_t), value :: offset, length
            integer(c_int), value    :: rw, scratch, create
        end subroutine map_file_window
    end interface

    interface
        subroutine unmap_file_window (addr, length) bind(C, name='unmap_file_window')
            import :: c_ptr, c_size_t
            type(c_ptr),        value :: addr
            integer(c_size_t),  value :: length
        end subroutine unmap_file_window
    end interface

contains

    !> \brief   Initialize file virtual memory mapping
    !> \author  J Benda
    !> \date    2020 - 2025
    !>
    !> Set up a data file mapping window. The file name, byte offset and byte length of the window within the file needs to be
    !> given. Internally, this opens the file using the POSIX function `open` (or `OpenFileMapping` in Windows), requests the
    !> mapping starting from the virtual memory page boundary just below the requested offset as needed by `mmap`, and constructs
    !> the pointer to the actually wanted window, which starts at the specified offset.
    !>
    !> If `create_new` is given and set to `.true.`, then the file is created first before mapping from it. If `read_write` is
    !> given and set to `.true.`, the mapped view will be writable. If `unique_name` is given and set to `.true.`, the provided
    !> file name is interpreted only as a prefix for a file path, which will be suffixed to generate a unique file path. (This
    !> makes sense only together with `create_new`.) In this last case the file will be considered temporary and removed once
    !> the file mapping ends.
    !>
    !> \param[in] this          A file mapping object.
    !> \param[in] filename      Name of file to open (or prefix if `unique_name` is given).
    !> \param[in] offset        Position in file where to start mapped segment.
    !> \param[in] length        Length of the mapped segment.
    !> \param[in] read_write    Whether to enable writing.
    !> \param[in] unique_name   Whether to suffix `filename` with characters in order to make a unique file name.
    !>
    subroutine init_file_mapping (this, filename, offset, length, read_write, unique_name, create_new)

        class(file_mapping), intent(inout) :: this
        character(len=*),    intent(in)    :: filename
        integer(int64),      intent(in)    :: offset, length
        logical, optional,   intent(in)    :: read_write, unique_name, create_new

        character(len=:, kind=c_char), allocatable :: cfilename
        integer(c_int) :: ccreate_new, cunique_name, cread_write

        cfilename = trim(filename) // c_null_char

        ccreate_new  = 0;  if (present(create_new))  ccreate_new  = merge(1, 0, create_new)
        cunique_name = 0;  if (present(unique_name)) cunique_name = merge(1, 0, unique_name)
        cread_write  = 0;  if (present(read_write))  cread_write  = merge(1, 0, read_write)

        this % length = length

        call map_file_window(cfilename,             &
                             int(offset, c_size_t), &
                             int(length, c_size_t), &
                             this % addr,           &
                             this % ptr,            &
                             cread_write,           &
                             cunique_name,          &
                             ccreate_new)

    end subroutine init_file_mapping


    !> \brief   Finalize file virtual memory mapping
    !> \author  J Benda
    !> \date    2020
    !>
    !> Cancel the file mapping. This simply calls the POSIX `munmap` function.
    !>
    subroutine finalize_file_mapping (this)

        class(file_mapping), intent(inout) :: this

        if (c_associated(this % addr)) then

            call unmap_file_window(this % addr, this % length)

            this % addr   = c_null_ptr
            this % ptr    = c_null_ptr
            this % length = 0

        end if

    end subroutine finalize_file_mapping


    !> \brief   Destructor of the type "file_mapping"
    !> \author  J Benda
    !> \date    2021
    !>
    !> Called when an instance of "file_mapping" is being destructed.
    !>
    subroutine destruct_file_mapping (this)

        type(file_mapping), intent(inout) :: this

        call this % finalize

    end subroutine destruct_file_mapping

end module file_mapping_gbl
