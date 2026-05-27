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

!> \brief   Timing module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Timing_Module

    use precisn,        only: longint, wp
    use const_gbl,      only: stdout
    use Utility_module, only: get_real_time, string_hash

    implicit none

    public master_timer

    private

    integer, parameter :: TIMER_DEFAULT_SIZE = 1000
    integer, parameter :: NAME_LEN           =   40 ! Max length of timer name

    type time_data
        logical                 :: used      ! Slot used?
        logical                 :: active    ! Currently active?
        character(len=name_len) :: name      ! Timer name
        integer(longint)        :: calls     ! Number of times the timer was invoked

        ! All times below are in seconds
        real(wp)               :: real_time  ! Total real time on this timer
        real(wp)               :: real_kids  ! Real time spent in nested timers
        real(wp)               :: real_start ! For active timers, time of activation
        integer(longint)       :: stack_p    ! For active timers, position in the stack
    end type time_data

    type :: Timer
        integer                :: process_id
        type(time_data)        :: timers(TIMER_DEFAULT_SIZE)
        integer                :: nested_timers(TIMER_DEFAULT_SIZE)
        integer                :: order(TIMER_DEFAULT_SIZE)
        integer                :: timer_count  = 0         ! Number of defined timers
        integer                :: timer_active = 0         ! Number of currently active timers
        real(wp)               :: program_start
        logical                :: initialized = .false.
    contains
        procedure, public      :: initialize
        procedure, public      :: start_timer
        procedure, public      :: stop_timer
        procedure, public      :: report_timers
        procedure, private     :: insert_time
       !procedure, private     :: sort_timers
    end type Timer

    type(Timer) :: master_timer

contains

    !> \brief   Initialize timers
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine initialize (this)
        class(Timer) :: this

        this % timers(:) % active = .false.
        this % timers(:) % used = .false.
        this % program_start = get_real_time()

        this % initialized = .true.

    end subroutine initialize


    !> \brief   Start a new named timer
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine start_timer(this, name)
        class(Timer)                 :: this
        character(len=*), intent(in) :: name  ! Timer name
        integer                      :: timer_idx

        timer_idx = this % insert_time(name)

        if (this % timers(timer_idx) % active) then
            write (stdout, "('TimerStart: timer ',a,' is already active')") trim(name)
            stop 'TimerStart - nested timer'
        end if

        !  Push the new timer to the timer stack
        this % timer_active = this % timer_active + 1
        this % nested_timers(this % timer_active) = timer_idx

        this % timers(timer_idx) % active     = .true.
        this % timers(timer_idx) % stack_p    = this % timer_active
        this % timers(timer_idx) % calls      = this % timers(timer_idx) % calls + 1
        this % timers(timer_idx) % real_start = get_real_time()
       !this % timers(timer_idx) % cpu_start  = get_cpu_time ()

    end subroutine start_timer


    !> \brief   Stop a named timer
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine stop_timer (this, name)
        class(Timer)                 :: this
        character(len=*), intent(in) :: name  ! Timer name
        integer                      :: timer_idx
        real(wp)                     :: real_time
        type(time_data), pointer     :: pt_timer

        timer_idx = this % insert_time(name)
        real_time = get_real_time() - this % timers(timer_idx) % real_start

        this % timers(timer_idx) % real_time = this % timers(timer_idx) % real_time + real_time
        this % timers(timer_idx) % active = .false.

        this % timer_active = this % timer_active - 1

        if (this % timer_active > 0) then
            this % timers(this % nested_timers(this % timer_active)) % real_kids &
                = this % timers(this % nested_timers(this % timer_active)) % real_kids + real_time
        end if

    end subroutine stop_timer


    !> \brief   Insert a new named timer
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer function insert_time (this, name)
        class(Timer)                 :: this
        character(len=*), intent(in) :: name  ! Timer name

        insert_time = string_hash(name, TIMER_DEFAULT_SIZE)

        search: do
            if (.not. this % timers(insert_time) % used) then
                ! This is a new key, insert it
                this % timer_count = this % timer_count + 1
                if (this % timer_count >= TIMER_DEFAULT_SIZE / 5) then
                    write (stdout, "('Too many timers. Increase table_size in timer.f90 to at least ',i5)") this % timer_count * 5
                    stop 'timer%insert_item'
                end if
                this % order(this % timer_count)       = insert_time
                this % timers(insert_time) % used      = .true.
                this % timers(insert_time) % active    = .false.
                this % timers(insert_time) % name      = name
                this % timers(insert_time) % calls     = 0
                this % timers(insert_time) % real_time = 0
               !this % timers(insert_time) % cpu_time  = 0
                this % timers(insert_time) % real_kids = 0
               !this % timers(insert_time) % cpu_kids  = 0
                exit search
            end if

            if (this % timers(insert_time) % name == name) then
                ! This is an existing key, simply return the location
                exit search
            end if

            insert_time = 1 + modulo(insert_time - 2, TIMER_DEFAULT_SIZE)
        end do search

    end function insert_time


    !> \brief   Print a table of timers
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine report_timers (this)
        class(Timer) :: this
        real(wp)     :: real_now, real_time, cpu_time, real_kids, cpu_kids, real_threshold
       !real(wp)     :: cpu_now, cpu_threshold
        integer      :: ord, pos, kid_pos, omitted

        type(time_data), pointer :: t, k
        character(len=1)         :: active

        real_now = get_real_time()
        real_threshold = 0.01_wp * (real_now - this % program_start)

        write (stdout, "(t2,'     ',t38,'     ',t45,'Total time (seconds)',t67,'Self time (seconds)')")
        write (stdout, "(t2,'Timer',t38,'Calls',t45,'--------------------',t67,'-------------------')")
        write (stdout, "(t2,'-----',t38,'-----',t50,'Real',t61,'CPU',t72,'Real',t83,'CPU')")
        write (stdout, *)

        omitted = 0

        scan: do ord = 1, this % timer_count

            pos = this % order(ord)
            if (.not. this % timers(pos) % used) then
                write (stdout, "('Timer ',i4,' in slot ',i5,' is defined but unused?!')") ord, pos
                stop 'TimerReport - logic error'
            end if

            ! Calculate active-timer corrections
            real_time = 0 ; real_kids = 0 ;
            cpu_time  = 0 ; cpu_kids  = 0 ;
            active     = ' '
            if (this % timers(pos) % active) then
                real_time = real_now - this % timers(pos) % real_start
                !cpu_time  = cpu_now  - t%cpu_start
                if (this % timer_active /= this % timers(pos) % stack_p) then
                    ! If we are not at the top of the stack, adjust
                    ! cumulative children time.
                    kid_pos = this % nested_timers(this % timers(pos) % stack_p + 1)
                    real_kids = real_now - this % timers(kid_pos) % real_start
                   !cpu_kids  = cpu_now  - k%cpu_start
                end if
                active = '*'
            end if

            real_time = real_time + this % timers(pos) % real_time
            !cpu_time  = cpu_time  + t%cpu_time
            real_kids = real_kids + this % timers(pos) % real_kids
            !cpu_kids  = cpu_kids  + t%cpu_kids

            !  Output needed?
            if (real_time < real_threshold) then
                omitted = omitted + 1
                cycle scan
            end if

            !  Output
            write (stdout, "(t2,a30,t33,a1,t35,I8,t45,2(f9.1,1x,f9.1,3x))") &
                this % timers(pos) % name, active, this % timers(pos) % calls, real_time, cpu_time, &
                real_time - real_kids, cpu_time - cpu_kids

        end do scan

        if (omitted > 0) then
            write (stdout, "(/' (',i3,' timers contributing less than 1% are not shown)')") omitted
        end if

        write (stdout, *)

    end subroutine report_timers

end module Timing_Module
