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

!> \brief   MOLPRO core potentials
!> \authors A Al-Refaie
!> \date    2017
!>
!> Based on BaseCorePotential.
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module MOLPROCorePotential_module

    use precisn,                  only: wp
    use const,                    only: stdout
    use BaseCorePotential_module, only: BaseCorePotential
    use input
    use mpi_mod,                  only: myrank

    implicit none

    integer, parameter          :: MAX_MOLPRO_ENERGIES = 3
    character(len=*), parameter :: MOLPRO_ENERGY_LINE = 'TOTAL ENERGIES'
    character(len=*), parameter :: MOLPRO_NUCLEAR = 'NUCLEAR REPULSION ENERGY'
    character(len=*), parameter :: MOLPRO_RHF_ENERGY_LINE ='HF-SCF'

    type, extends(BaseCorePotential) :: MOLPROCorePotential
        integer, allocatable :: molpro_symmetry_map(:)
        real(wp)             :: nuclear_energy
    contains
        procedure, public    :: parse_ecp => parse_molpro_casscf_energies
        procedure, public    :: parse_input => parse_molpro_input
        procedure, public    :: modify_target_energies => molpro_energies
        procedure, private   :: search_state_energy
        procedure, private   :: search_rhf_state_energy
        procedure, private   :: read_state_energy
        procedure, private   :: search_core_energy
        procedure, private   :: detect_symmetry_ordering
    end type

contains

    !> \brief   Main build routine of the hamiltonian
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> All build must be done within this routine in order to be used by MPI-SCATCI.
    !>
    !> \param[inout] this    Vector object to update.
    !> \param[out] filename  Path to file with molpro output.
    !>
    subroutine parse_molpro_casscf_energies (this, filename)
        class(MOLPROCorePotential) :: this
        character(len=*) :: filename

        integer  :: targ_sym_do,current_target_sym
        integer  :: targ_stat_do,err
        integer  :: io
        real(wp) :: energy_val
        logical  :: eof
        logical  :: success
        integer  :: actual_targ_sym

        io = find_io(1) + myrank

        open (io, action = 'read', status = 'old', file = filename, iostat = err)

        if (err /= 0) then
            write (stdout, "('Error opening input file ',a)") filename
            stop "Could not open input file"
        end if

        success = this % search_core_energy(io)

        if (success == .false.) then
            stop "Could not find core energy in molpro"
        end if

        !success = this%detect_symmetry_ordering(io)
        !if(success ==.false.) then
        !    stop "Could not assign MOLPRO symmetries to SCATCI targets"
        !endif

        eof = this % search_rhf_state_energy(io, 1, 1, this % target_energies(:,1))
        if (.not. eof) stop "Could not find HF energy"
        !do targ_sym_do=1,this%num_target_sym
        !     ! current_target_sym = current_target_sym + 1
            !      !do  targ_stat_do=1,this%num_target_per_symmetries(targ_sym_do)
                 ! actual_targ_sym = this%molpro_symmetry_map(targ_sym_do)
            !           !energy_line = '!MCSCF STATE 1.1 Energy'
                  !eof = this%search_state_energy(io,targ_sym_do,this%num_target_per_symmetries(targ_sym_do),this%target_energies(:,actual_targ_sym))
             !            !length = len(trim(energy_line))
                     !
                   !  if(eof == .false.) then
              !               ! write(stdout,"(a,2x,f)") energy_line
            !      ! enddo
            !enddo

        close (io)

        this % target_energies(:,:) = this % target_energies(:,:) - this % nuclear_energy

    end subroutine parse_molpro_casscf_energies


    !> \brief   Main build routine of the hamiltonian
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> All build must be done within this routine in order to be used by MPI-SCATCI.
    !>
    subroutine parse_molpro_input(this)
        class(MOLPROCorePotential) :: this
    end subroutine parse_molpro_input


    subroutine molpro_energies (this, target_sym, target_energy)
        class(MOLPROCorePotential) :: this
        integer, intent(in) :: target_sym
        real(wp)            :: target_energy(:)
        integer             :: num_states

        num_states = this % num_target_per_symmetries(target_sym)

        if (target_sym == 1) then
            this % ground_energy = target_energy(1)
            !this%L2_ECP_deficit = (this%target_energies(1,1) -  target_energy(1))
            this % L2_ECP_deficit = 0.1
            print *, 'ground_energy ',this % ground_energy
            print *, 'target_ECP_energy ',this % target_energies(1,1)
            print *, 'L2 diff ', this % L2_ECP_deficit
        end if
        !target_energy(1:num_states) =target_energy(1:num_states) + 0.455207
        !find deficieny
        !target_energy(1:num_states) =target_energy(1:num_states) + this%target_energies(1,1)
    end subroutine molpro_energies


    logical function search_state_energy (this, io, target_sym, number_of_energies, target_energy)
        class(MOLPROCorePotential) :: this
        integer,  intent(in)       :: io, target_sym,number_of_energies
        real(wp), intent(out)      :: target_energy(:)
        integer                    :: reason
        character(len=line_len)    :: buffer
        character(len=line_len)    :: header_line
        character(len=4)           :: sym_char, stat_char

        search_state_energy = .true.

        write (sym_char, '(i4)') target_sym

        !CI vector for state symmetry 1
        !header_line = '!MCSCF STATE '//trim(adjustl(stat_char))//'.'//trim(adjustl(sym_char))//' Energy'

        header_line = 'CI vector for state symmetry ' // trim(adjustl(sym_char))

        do while (.true.)
            read (io, '(a)', iostat = reason) buffer
            !write(stdout,"(a)") trim(buffer(2:line_len))
            !write(stdout,"(a)") trim(header_line)
            if (reason == 0) then
                if (trim(buffer(2:line_len)) == trim(header_line)) then
                    search_state_energy = this % read_state_energy(io, number_of_energies, target_energy)
                    !write(stdout,*) buffer
                    return
                end if
            else
                search_state_energy = .false.
                exit
            end if
        end do
    end function search_state_energy


    logical function search_rhf_state_energy (this, io, target_sym, number_of_energies, target_energy)
        class(MOLPROCorePotential) :: this
        integer,  intent(in)       :: io, target_sym,number_of_energies
        real(wp), intent(out)      :: target_energy(:)
        integer                    :: reason
        character(len=line_len)    :: buffer
        character(len=line_len)    :: header_line
        character(len=4)           :: sym_char, stat_char

        search_rhf_state_energy = .true.

        !write( sym_char, '(i4)') target_sym
        !CI vector for state symmetry 1
        !header_line = '!MCSCF STATE '//trim(adjustl(stat_char))//'.'//trim(adjustl(sym_char))//' Energy'
        !header_line = 'CI vector for state symmetry '//trim(adjustl(sym_char))

        do while (.true.)
            read (io, '(a)', iostat = reason) buffer
            write (stdout, '(a)') trim(buffer(1:line_len))
            !write(stdout,"(a)") trim(header_line)
            if (reason == 0) then
                if (trim(buffer(9:14)) == trim(MOLPRO_RHF_ENERGY_LINE)) then
                    read (io, '(a)', iostat = reason) buffer
                    call parse(buffer)
                    call readf(target_energy(1))
                    write (stdout, *) buffer
                    return
                end if
            else
                search_rhf_state_energy = .false.
                exit
            end if
        end do
    end function search_rhf_state_energy


    logical function read_state_energy (this, io, number_of_energies, target_energy)
        class(MOLPROCorePotential) :: this
        integer,  intent(in)       :: io, number_of_energies
        real(wp), intent(out)      :: target_energy(:)
        character(len=line_len)    :: header_line
        integer                    :: num_found_energies, total_energies_line, curr_en, reason
        character(len=line_len)    :: buffer
        character(len=line_len)    :: fields

        num_found_energies = 0

        do while (.true.)
            read (io, '(a)', iostat = reason) buffer
            if (reason == 0) then
                if (trim(buffer(2:len(MOLPRO_ENERGY_LINE)+1)) == trim(MOLPRO_ENERGY_LINE)) then
                    call parse(buffer)
                    total_energies_line = nitems - 2
                    !write(stdout,*) total_energies_line
                    call readu(fields)
                    call readu(fields)
                    num_found_energies = 0
                    do while (total_energies_line > 0)
                        num_found_energies = num_found_energies + 1

                        if (num_found_energies > number_of_energies) exit

                        call readf(target_energy(num_found_energies))
                        !write(stdout,"(f12.6)") target_energy(num_found_energies)
                        total_energies_line = total_energies_line - 1

                        if (total_energies_line == 0) then
                            read (io, '(a)', iostat = reason) buffer
                            call parse(buffer)
                            if (nitems == 0) exit
                            total_energies_line = nitems - 2
                        end if
                    end do
                    read_state_energy = .true.
                    return
                end if
            else
                read_state_energy = .false.
                return
            end if
        end do

    end function read_state_energy


    logical function search_core_energy (this, io)
        class(MOLPROCorePotential) :: this
        integer, intent(in)        :: io
        integer                    :: reason
        character(len=line_len)    :: buffer
        character(len=line_len)    :: fields

        do while(.true.)
            read (io, '(a)', iostat = reason) buffer
            !write(stdout,"(a)") trim(buffer(2:line_len))
            !write(stdout,"(a)") trim(header_line)
            if (reason == 0) then
                if (trim(buffer(2:len(MOLPRO_NUCLEAR) + 1)) == trim(MOLPRO_NUCLEAR)) then
                    call parse(buffer)
                    call readu(fields)
                    call readu(fields)
                    call readu(fields)
                    call readf(this % nuclear_energy)
                    !write(stdout,"(f12.6)") this%nuclear_energy

                    search_core_energy = .true.

                    return
                end if
            else
                search_core_energy = .false.
                exit
            end if
        end do

    end function search_core_energy


    logical function detect_symmetry_ordering (this, io)
        class(MOLPROCorePotential) :: this
        integer, intent(in)        :: io
        integer                    ::    reason
        character(len=line_len)    ::    buffer
        character(len=line_len)    ::    fields
        character(len=line_len)    ::    header_line
        character(len=4)           ::    sym_char
        integer                    ::    target_sym, ido
        integer                    ::    detected_space, detected_spin

        allocate(this % molpro_symmetry_map(this % num_target_sym))

        do target_sym = 1, this % num_target_sym

            write (sym_char, '(i4)') target_sym
            header_line = 'State symmetry ' // trim(adjustl(sym_char))

            do while (.true.)

                read (io, '(a)', iostat = reason) buffer

                if (reason == 0) then
                    if (trim(buffer(2:line_len)) == trim(header_line)) then
                        !Empty_line
                        read(io,'(a)',iostat=reason) buffer
                        !Info
                        read(io,'(a)',iostat=reason) buffer

                        do ido = 1, len(buffer)
                            if (buffer(ido:ido) == "=") buffer(ido:ido) = " "
                        end do

                        call parse(buffer)
                        !write(stdout,*) buffer
                        !Number
                        call readu(fields)
                        !of
                        call readu(fields)
                        !electrons
                        call readu(fields)
                        !#
                        call readu(fields)

                        call readu(fields)

                        if (trim(fields) == 'SPIN') then
                            call readu(fields)
                            call readu(fields)
                            select case(trim(fields))
                                case('SINGLET')
                                    detected_spin = 1
                                case('DOUBLET')
                                    detected_spin = 2
                                case('TRIPLET')
                                    detected_spin = 3
                                case default
                                    write (stdout, *) fields
                                    stop "Invalid MOLPRO input on spin symmetry"
                            end select
                        else
                            write (stdout, *) fields
                            stop "Invalid MOLPRO input on spin expected first"
                        end if
                        call readu(fields)
                        if (trim(fields) == 'SPACE') then
                            call readu(fields)
                            call readi(detected_space)
                        else
                            write (stdout, *) fields
                            stop "Invalid MOLPRO input on space expected second"
                        end if

                        do ido = 1, this % num_target_sym
                            if (this % spatial_symmetry(ido) + 1 == detected_space .and. &
                                this % spin_symmetry(ido) == detected_spin) then
                                this % molpro_symmetry_map(target_sym) = ido
                                exit
                            else if (ido == this % num_target_sym) then
                                write (stdout, *) 'Could not find molpro symmetry mapping of Spatial=', &
                                    detected_space, ' and Spin=', detected_spin
                                stop "Invalid molpro space spin assignment"
                            end if
                        end do

                        exit
                    end if
                else
                    detect_symmetry_ordering = .false.
                    return
                end if
            end do
        end do

        detect_symmetry_ordering = .true.

        write (stdout, *) ' '
        write (stdout, *) '---------------------------------'
        write (stdout, *) 'MOLPRO -> SCATCI symmetry mapping'
        write (stdout, *) '---------------------------------'
        write (stdout, *) ' '

        do ido = 1, this % num_target_sym
            write (stdout, *) ido, '----->', this % molpro_symmetry_map(ido)
        end do

    end function detect_symmetry_ordering

end module MOLPROCorePotential_module
