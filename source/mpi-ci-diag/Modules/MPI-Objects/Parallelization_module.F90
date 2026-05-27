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

!> \brief   Distribution of processes into a grid
!> \authors J Benda
!> \date    2019
!>
!> This module contains utility routines and types which aid with the "parallelization or parallelization", i.e.,
!> with concurrent distributed diagonalizations of Hamiltonians for multiple irreducible representations.
!>
module Parallelization_module

#ifdef usempi
    use mpi
#endif
    use blas_lapack_gbl, only: blasint
    use mpi_gbl,         only: mpiint, myrank, nprocs, local_rank, local_nprocs, shared_communicator

    implicit none

    !> \brief   MPI process grid layout
    !> \authors J Benda
    !> \date    2019
    !>
    type :: ProcessGrid
        integer(blasint) :: wcntxt   !< BLACS context containing all MPI processes in the world communicator
        integer(blasint) :: wprows   !< MPI world communicator grid row count
        integer(blasint) :: wpcols   !< MPI world communicator grid column count
        integer(blasint) :: mywrow   !< row position of this process within the MPI world communicator
        integer(blasint) :: mywcol   !< column position of this process within the MPI world communicator

        integer :: igroup   !< zero-based index of the MPI group this process belongs to
        integer :: ngroup   !< total number of MPI groups partitioning the world communicator

        integer :: gprocs   !< number of processes in the current MPI group
        integer :: grank    !< rank of this process within the MPI group
        integer :: lprocs   !< number of processes of the current MPI group localized on a single node
        integer :: lrank    !< rank of this process within the local MPI group

        integer(blasint) :: gcntxt   !< BLACS context containing all MPI processes in the MPI group communicator
        integer(blasint) :: gprows   !< MPI group communicator grid row count
        integer(blasint) :: gpcols   !< MPI group communicator grid column count
        integer(blasint) :: mygrow   !< row position of this process within the MPI group communicator
        integer(blasint) :: mygcol   !< column position of this process within the MPI group communicator

        integer(mpiint) :: gcomm    !< MPI group communicator
        integer(mpiint) :: lcomm    !< subset of the MPI group communicator localized on a single node

        integer, allocatable :: groupnprocs(:)  !< Number of processes per MPI group.

        logical :: sequential       !< Whether the diagonalizations will be done sequentially (one after another) or not
    contains
        procedure, public,  pass   :: setup => setup_process_grid
        procedure, public,  pass   :: is_my_group_work
        procedure, public,  pass   :: which_group_is_work
        procedure, public,  pass   :: group_master_world_rank
        procedure, public,  pass   :: summarize
        procedure, private, nopass :: square
    end type ProcessGrid

    ! global instance of the process grid
    type(ProcessGrid) :: process_grid

contains

    !> \brief   Initialize the process grid
    !> \authors J Benda
    !> \date    2019
    !>
    !> Splits the world communicator to the given number of MPI groups. Sets up all global group communicators and local group
    !> communicators (subset on one node). If there are more groups than processes, then all groups are equal to the MPI
    !> world. If the number of processes is not divisible by the number of groups, the leading mod(nprocs,ngroup) processes will
    !> contain 1 process more than the rest of the groups.
    !>
    !> \param this        Process grid to set up.
    !> \param ngroup      Number of MPI groups to create.
    !> \param sequential  If "true", all diagonalizations will be done in sequence (not concurrently, even if there are
    !>                    enough CPUs to create requested groups). This is needed to have the eigenvectors written to disk,
    !>                    which does not happen with concurrent diagonalizations.
    !>
    subroutine setup_process_grid (this, ngroup, sequential)

        use const_gbl, only: stdout

        class(ProcessGrid), intent(inout) :: this
        integer,            intent(in)    :: ngroup
        logical,            intent(in)    :: sequential

        integer, allocatable :: groupstarts(:), groupends(:)
        integer(mpiint)      :: n, color, ierr
        integer(blasint)     :: cntxt, rows, cols
        integer              :: i, j, k, rank
        integer(blasint), allocatable :: groupmap(:,:)

        allocate (this % groupnprocs(ngroup), groupstarts(ngroup), groupends(ngroup))

        ! no context splitting needed if too few processes, or if concurrent diagonalizations are not desired
        this % sequential = (sequential .or. ngroup == 1 .or. nprocs < ngroup)

        ! default initialization
        this % gcomm  = 0
        this % lcomm  = shared_communicator
        this % igroup = 0
        this % ngroup = merge(1, ngroup, this % sequential)
        this % gprocs = nprocs
        this % grank  = myrank
        this % lprocs = local_nprocs
        this % lrank  = local_rank
        this % groupnprocs = nprocs

#ifdef usempi
#   ifdef scalapack
        ! allocate BLACS context for the world
        call this % square(int(nprocs), this % wprows, this % wpcols)
        call blacs_get(-1_blasint, 0_blasint, this % wcntxt)
        call blacs_gridinit(this % wcntxt, 'r', this % wprows, this % wpcols)
        call blacs_gridinfo(this % wcntxt, this % wprows, this % wpcols, this % mywrow, this % mywcol)
#   endif

        if (this % sequential) then
            this % gcntxt = this % wcntxt
            this % gcomm  = MPI_COMM_WORLD
            this % groupnprocs = nprocs
            this % gprows = this % wprows
            this % gpcols = this % wpcols
            this % mygrow = this % mywrow
            this % mygcol = this % mywcol
            return
        end if

        ! spread processes among groups (leading groups will end with more processes on imbalance)
        this % groupnprocs(:) = 0
        do i = 1, nprocs
            j = 1 + mod(i - 1, this % ngroup)
            this % groupnprocs(j) = this % groupnprocs(j) + 1
        end do
        groupstarts(1)    = 0
        groupends(this % ngroup) = nprocs - 1
        do i = 1, this % ngroup - 1
            groupstarts(i + 1) = groupstarts(i) + this % groupnprocs(i)
            groupends(i) = groupstarts(i + 1) - 1
        end do

        ! find out which group this process belongs to (store zero-based index)
        do i = 1, this % ngroup
            if (groupstarts(i) <= myrank .and. myrank <= groupends(i)) then
                this % igroup = i - 1
                exit
            end if
        end do

        ! create a MPI group and split it to sub-groups on individual nodes
        color = this % igroup
        call MPI_Comm_split(MPI_COMM_WORLD, color, myrank, this % gcomm, ierr)
        call MPI_Comm_split_type(this % gcomm, MPI_COMM_TYPE_SHARED, myrank, MPI_INFO_NULL, this % lcomm, ierr)
        call MPI_Comm_rank(this % gcomm, n, ierr); this % grank  = n
        call MPI_Comm_size(this % gcomm, n, ierr); this % gprocs = n
        call MPI_Comm_rank(this % lcomm, n, ierr); this % lrank  = n
        call MPI_Comm_size(this % lcomm, n, ierr); this % lprocs = n

#   ifdef scalapack
        ! Allocate BLACS context for the group. This has to be done separately for each group,
        ! not only once with a group-dependent content of `groupmap`, because `blacs_gridmap` uses
        ! `MPI_Group_incl`, which does not support creation of several groups at one time.

        write (stdout, '(/," Creating BLACS groups")')
        write (stdout, '(  " ---------------------")')

        do i = 1, this % ngroup

            ! populate matrix of ranks that belong to the this group's grid
            call this % square(this % groupnprocs(i), rows, cols)
            allocate (groupmap(rows, cols))
            rank = 0
            do k = 1, cols
                do j = 1, rows
                    groupmap(j, k) = groupstarts(i) + rank
                    rank = rank + 1
                end do
            end do

            ! split off a new BLACS context in a collective call
            call blacs_get(-1_blasint, 0_blasint, cntxt)
            write (stdout, '(" Group #",I0," size ",I0,"x",I0," grid ",*(1x,I0))') i, rows, cols, groupmap
            call blacs_gridmap(cntxt, groupmap, rows, rows, cols)
            deallocate (groupmap)

            ! if this group is mine, store the group information
            if (i == this % igroup + 1) then
                this % gprows = rows
                this % gpcols = cols
                this % gcntxt = cntxt
            end if

        end do

        ! find out this process' position within its BLACS grid
        call blacs_gridinfo(this % gcntxt, this % gprows, this % gpcols, this % mygrow, this % mygcol)
#   endif
#endif

    end subroutine setup_process_grid


    !> \brief   Write current grid layout to stdout
    !> \authors J Benda
    !> \date    2019
    !>
    subroutine summarize (this)

        use const_gbl, only: stdout

        class(ProcessGrid), intent(in) :: this

        integer :: i

        write (stdout, *)
        write (stdout, '(1x,A)') 'Process grid'
        write (stdout, '(1x,A)') '------------'
        write (stdout, '(1x,A,L1)') 'Sequential diagonalizations: ', this % sequential
        write (stdout, '(1x,A,I0)') 'Number of groups: ', this % ngroup
        write (stdout, '(1x,A)', advance = 'no') 'Number of processes per group:'
        do i = 1, this % ngroup
            write (stdout, '(1x,I0)', advance = 'no') this % groupnprocs(i)
        end do
        write (stdout, *)
        write (stdout, '(1x,A,I0)') 'This process belongs to group: ', this % igroup
        write (stdout, '(1x,A,I0)') 'This group size: ', this % gprocs
        write (stdout, '(1x,A,I0)') 'This group rank: ', this % grank
        write (stdout, '(1x,A,I0)') 'Shared-memory sub-group size: ', this % lprocs
        write (stdout, '(1x,A,I0)') 'Shared-memory sub-group rank: ', this % lrank
        write (stdout, *)

    end subroutine summarize


    !> \brief   Check whether this work-item is to be processed by this process' group
    !> \authors J Benda
    !> \date    2019
    !>
    !> The work-item index is expected to be greater than or equal to 1. Work-items (indices) larger
    !> than the number of MPI groups wrap around.
    !>
    logical function is_my_group_work (this, i) result (verdict)

        class(ProcessGrid), intent(inout) :: this
        integer,            intent(in)    :: i

        verdict = (this % which_group_is_work(i) == this % igroup)

    end function is_my_group_work


    !> \brief   Find out which group workitem will be processed by
    !> \authors J Benda
    !> \date    2019
    !>
    !> The work-item index is expected to be greater than or equal to 1. Work-items (indices) larger
    !> than the number of MPI groups wrap around. Groups are numbered from 0.
    !>
    integer function which_group_is_work (this, i) result (igroup)

        class(ProcessGrid), intent(inout) :: this
        integer,            intent(in)    :: i

        igroup = mod(i - 1, this % ngroup)

    end function which_group_is_work


    !> \brief   Find out world rank of the master process of a given MPI group
    !> \authors J Benda
    !> \date    2019
    !>
    integer function group_master_world_rank (this, igroup) result (rank)

        class(ProcessGrid), intent(inout) :: this
        integer,            intent(in)    :: igroup

        rank = sum(this % groupnprocs(1:igroup))

    end function group_master_world_rank


    !> \brief   Given integer area of a box, calculate its edges
    !> \authors A Al-Refaie, J Benda
    !> \date    2017 - 2019
    !>
    !> Return positive `a` and `b` such that their product is exactly equal to `n`, and the difference
    !> between `a` and `b` is as small as possible. On return `a` is always less than or equal to `b`.
    !>
    subroutine square (n, a, b)

        integer, intent(in)  :: n
        integer(blasint), intent(out) :: a, b

        integer :: i

        do i = 1, int(sqrt(real(n)) + 1)
            if (mod(n, i) == 0) then
                a = i
                b = n / i
            end if
        end do

    end subroutine square

end module Parallelization_module
