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

!> \brief   Base MPI module
!> \authors A Al-Refaie
!> \date    2017
!>
!> Uses MPI module from GBTOlib.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module BaseMPI_module

    use mpi_mod, only: check_mpi_running, myrank, nprocs

    implicit none

    public BaseMPI

    private

    !> \brief Base MPI type
    !>
    !> This abstract class that most will inherit from to gain automatic information on rank and number of processes
    !>
    type, abstract :: BaseMPI
        private
        !> MPI Rank
        integer :: rank
        !> Number of processes
        integer :: nProcs
    contains
        !Constructor function
        procedure, pass(this) :: basempi_ctor
        !>Only allows process 0 to print
        procedure, public     :: log_message
        !> Any process prints
        procedure, public     :: verbose_log_message
        !> Which constructor is invoked
        generic, public       :: construct => basempi_ctor
    end type BaseMPI

contains

    !> \brief   Constructs the BaseMPI object by assigning rank and nProcs
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Checks if MPI is running; if it is then assign values.
    !>
    subroutine basempi_ctor (this)
        class(BaseMPI) :: this

        call check_mpi_running

        this % rank   = myrank
        this % nProcs = nprocs

    end subroutine basempi_ctor

    !> \brief Will log a message whilst attaching a number in the beginning
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine verbose_log_message (this, message)
        class(BaseMPI),   intent(in) :: this
        character(len=*), intent(in) :: message

        print *, '[', this % rank, ']', message

    end subroutine verbose_log_message


    !> \brief Will only log if the rank is zero
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine log_message (this, message)
        class(BaseMPI)               :: this
        character(len=*), intent(in) :: message

        if (this % rank == 0 ) then
            call this % verbose_log_message(message)
        end if

    end subroutine log_message

end module BaseMPI_module
