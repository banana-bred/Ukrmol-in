! Copyright 2019
!
! For a comprehensive list of the developers that contributed to these codes
! see the UK-AMOR website.
!
! This file is part of UKRmol-in (UKRmol+ suite).
!
!     UKRmol-in is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     UKRmol-in is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  UKRmol-in (in source/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.

!> \brief   Memory manager module
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module has an extremely basic memory management system in place, this will hopefully be upgraded in due time
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial revision.
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module MemoryManager_module

    use precisn,   only: longint, wp
    use const_gbl, only: stdout

    implicit none

    public master_memory

    private

    !> \brief   This is a simple class to handle memory management tracking
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    type :: MemoryManager
        private
        integer(longint) :: total_local_memory = 0      !< The total memory we have been assigned
        integer(longint) :: available_local_memory = 0  !< Memory we have left over
        integer(longint) :: cached_tracked_memory = 0
        logical          :: initialized = .false.
    contains
        private
        procedure, public :: construct => Memorymanager_ctor
        procedure, public :: init_memory  !< Initialize the memory with a new value
        procedure, public :: track_memory !< Tracks
        procedure, public :: free_memory
        procedure, public :: get_available_memory
        procedure, public :: get_scaled_available_memory
        procedure, public :: print_memory_report
    end type MemoryManager

    !> \brief   This is the global memory manager.
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> No other memory manager can exist
    !>
    type(MemoryManager) :: master_memory

contains

    !> \brief   Constructor, used to initialize total available memory
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this       Manager object to update.
    !> \param[in] total_memory  Total available memory for the process in bytes
    !>
    subroutine Memorymanager_ctor (this, total_memory)
        class(MemoryManager)         :: this
        integer(longint), intent(in) :: total_memory

        this % total_local_memory = total_memory
        this % available_local_memory = this % total_local_memory
        this % initialized  = .true.
        this % available_local_memory = this % available_local_memory - this % cached_tracked_memory

    end subroutine Memorymanager_ctor


    !> \brief   Same as constructor, used to reinitialize te memory again (not really used)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this       Manager object to update.
    !> \param[in] total_memory  Total available memory for the process in bytes
    !>
    subroutine init_memory (this, total_memory)
        class(MemoryManager)         :: this
        integer(longint), intent(in) :: total_memory

        this % total_local_memory = total_memory
        this % available_local_memory = this % total_local_memory
        this % available_local_memory = this % available_local_memory - this % cached_tracked_memory
        this % initialized = .true.

    end subroutine init_memory


    !> \brief   Tracks memory allocation, usually called after an allocation
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this    Manager object to update.
    !> \param[in] elem_size  Size in byte of each element
    !> \param[in] nelem      Number of elements of size elem_size
    !> \param[in] stat error Value after an allocate call
    !> \param[in] array_name Assigned name of array, used to debug where a memory allocation error occured.
    !>
    subroutine track_memory (this, elem_size, nelem, stat, array_name)
        class(MemoryManager)         :: this
        character(len=*), intent(in) :: array_name
        integer,          intent(in) :: elem_size,nelem,stat
        integer(longint)             :: total_mem

        !Check if allocation was successfull
        if (stat /= 0) then
            !If not print name and error stats
            write (stdout, "('Array : ',a)") array_name
            write (stdout, "('stat returned error ',i16,' when trying to allocate nelems ',i16,' of size ', i16)") &
                stat, nelem, elem_size
            !Halt the program
            stop "[Memory manager] Memory allocation error"
        end if

        ! calculate in bytes total memory used
        total_mem = elem_size * nelem

        !Remove from available memory
        !$OMP ATOMIC
        this % cached_tracked_memory = this % cached_tracked_memory + total_mem
        if (.not. this % initialized) return
        !$OMP ATOMIC
        this % available_local_memory = this % available_local_memory - total_mem

        !If we've hit negative then throw error
        if (this % available_local_memory < 0) then
            call this % print_memory_report()
            write (stdout, "('Array : ',a)") array_name
            write (stdout, "('Overrun allowed memory space when trying to allocate nelems ',i16,' of size ', 2i16,' vs',i16)") &
                nelem, elem_size, total_mem, this % get_available_memory()
            stop "[Memory manager] Overrun space"
        end if
    end subroutine track_memory


    !> \brief   Tracks memory deallocation, usually called before a deallocation
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \param[inout] this    Manager object to update.
    !> \param[in] elem_size  Size in byte of each element
    !> \param[in] nelem      Number of elements of size elem_size
    !>
    subroutine free_memory (this, elem_size, nelem)
        class(MemoryManager) :: this
        integer, intent(in)  :: elem_size, nelem
        integer(longint)     :: total_mem

        ! Calculate total deallocated
        total_mem = elem_size * nelem
        !$OMP ATOMIC
        this % cached_tracked_memory = this % cached_tracked_memory - total_mem
        !Add memory back to pool
        if (.not. this % initialized) return
        !$OMP ATOMIC
        this % available_local_memory = this % available_local_memory + total_mem

        !If we have more memory than available then flag as there is a mismatching of track and free
        if (this % available_local_memory > this % total_local_memory) then
            stop "Memory mismatch"
        end if

    end subroutine free_memory


    !> \brief   Get currently available memory
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> \result Total memory in bytes
    !>
    integer(longint) function get_available_memory (this)
        class(MemoryManager) :: this

        get_available_memory = this % available_local_memory

    end function get_available_memory


    !> \brief   Get currently available memory scaled.
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This is particularly useful when trying to determine how many elements we can fit into memory by giving a margin of safety
    !>
    !> \param[inout] this    Manager object to update.
    !> \param[in]    scale_  Scale valued of available memory (0.0 - 1.0).
    !>
    !> \result Total memory in bytes
    !>
    integer(longint) function get_scaled_available_memory (this, scale_)
        class(MemoryManager) :: this
        real(wp), intent(in) :: scale_
        real(wp)             :: s

        !Clamp scale between 0.0 ad 1.0
        s = min(abs(scale_), 1.0_wp)

        !Return scaled
        get_scaled_available_memory = int(real(this % available_local_memory) * s, longint)

    end function get_scaled_available_memory


    !> \brief   Prints a very simple memory report
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine print_memory_report (this)
        class(MemoryManager) :: this
        real(wp) :: mem_total, mem_avail, percentage

        !Compute total in gibibytes
        mem_total = real(this % total_local_memory) / real(1024**3)

        !Compute available in gibibytes
        mem_avail = real(this % available_local_memory) / real(1024**3)

        !Compute percentage available
        percentage = mem_avail * 100.0 / mem_total
        write (stdout, "('Currently we have ',f18.8,'GB / ',f18.8,'GiB available')") mem_avail, mem_total
        write (stdout, "('Memory availability at ',f6.1,'%')") percentage

    end subroutine print_memory_report

end module MemoryManager_module
