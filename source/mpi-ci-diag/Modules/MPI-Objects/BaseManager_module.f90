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

!> \brief   Base manager type
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module contains the BaseManager class that is used to track all memory allocations within a single process,
!> as of now it is fairly simple in implementation.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module BaseManager_module

    use precisn,        only: longint
    use BaseMPI_module, only: BaseMPI

    implicit none

    public BaseManager, memory_man

    private

    !>This is a simple class to handle memory management tracking
    type, extends(BaseMPI) :: BaseManager
        private
        integer(longint) :: total_local_memory     !< The total memory we have been assigned
        integer(longint) :: available_local_memory !< Memory we have left over
    contains
        private
        procedure :: basemanager_ctor
        procedure :: init_memory  !< Initialize the memory with a new value
        procedure :: track_memory !< Tracks
        procedure :: free_memory
        procedure :: get_available_memory
        procedure :: get_available_global_memory
        generic, public :: construct => basemanager_ctor
    end type BaseManager

    !Declare a static BaseManager for entire program scope
    type(BaseManager) :: memory_man

    !interface    BaseManager
    !    procedure    ::    BaseManager_ctor
    !end interface

contains

    !> \brief   Type constructor
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine basemanager_ctor (this, total_memory)
        class(BaseManager)           :: this
        integer(longint), intent(in) :: total_memory

        this % total_local_memory = total_memory
        this % available_local_memory = this % total_local_memory
        call this % construct

    end subroutine


    !> \brief   Initialize memory
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine init_memory (this, total_memory)
        class(BaseManager)           :: this
        integer(longint), intent(in) :: total_memory

        this % total_local_memory = total_memory
        this % available_local_memory = this % total_local_memory

    end subroutine init_memory


    !> \brief   Track memory
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine track_memory (this, alloc_memory, error)
        class(BaseManager)            :: this
        integer(longint), intent(in)  :: alloc_memory
        integer,          intent(out) :: error

        this % available_local_memory = this % available_local_memory - alloc_memory

        if (this % available_local_memory < 0) then
            error = -1
            call this % verbose_log_message('Run out of memory')
        else
            error = 0
        end if

    end subroutine track_memory


    !> \brief   Free memory
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine free_memory (this, alloc_memory, error)
        class(BaseManager)            :: this
        integer(longint), intent(in)  :: alloc_memory
        integer,          intent(out) :: error

        this % available_local_memory = this % available_local_memory + alloc_memory

        if (this % available_local_memory > this % total_local_memory) then
            error = -1
            call this % verbose_log_message('Memory mismatch')
        else
            error = 0
        end if

    end subroutine free_memory


    !> \brief   Obtain available local memory
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer(longint) function get_available_memory (this)
        class(BaseManager) :: this

        get_available_memory = this % available_local_memory

    end function get_available_memory


    !> \brief   Obtain available global memory
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer(longint) function get_available_global_memory (this)
        class(BaseManager) :: this
        integer(longint)   :: local(1), global(1)

        local = this % available_local_memory
        !call mpi_reduceall_inplace_sum_cfp(local,1)
        get_available_global_memory = local(1)

    end function get_available_global_memory

end module BaseManager_module
