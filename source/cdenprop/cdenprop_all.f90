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

program cdenprop_all

    use cdenprop_defs,                   only: CIvect, tsym_max
    use cdenprop_procs,                  only: cdenprop_drv
    use class_molecular_properties_data, only: molecular_properties_data
    use class_namelists,                 only: denprop_namelist, cdenprop_namelists
    use const_gbl,                       only: stdin, stdout, level1, set_verbosity_level
    use global_utils,                    only: print_ukrmol_header
    use mpi_gbl,                         only: mpi_mod_start, mpi_mod_finalize

    implicit none

    integer, parameter :: maxtgt = tsym_max
    integer, parameter :: maxpfg = 20
    integer, parameter :: max_multipole = 2  ! see construct_namelists_from_denprop_input
    logical, parameter :: master_writes_to_stdout = .true.

    type(denprop_namelist)          :: den_namelist
    type(cdenprop_namelists)        :: cden_namelist
    type(molecular_properties_data) :: mol_prop_data

    type (CIvect) :: ci_vec_i, ci_vec_j  ! dummy parameters for cdenprop call

    integer :: isym, jsym, no_symmetry_pairs, lupropw_final

    call mpi_mod_start(.true.)
    call print_ukrmol_header(stdout)

    write (stdout, '(/,20X,"CDENPROP ALL (2015) ",//,20X,"All states transition moments ",//)')

    den_namelist % zlast = .false.

    namelist_loop: do while (.not. den_namelist % zlast)

        call den_namelist % init(maxtgt, maxpfg)
        call den_namelist % read_namelist(stdin)
        call set_verbosity_level(den_namelist % verbosity)

        lupropw_final     = den_namelist % nftmt
        no_symmetry_pairs = ((den_namelist % ntgt + 1) * den_namelist % ntgt) / 2

        if ( den_namelist % ipol >= 0 ) then

            call mol_prop_data % preallocate_property_blocks(no_symmetry_pairs * (max_multipole + 1)**2)

            write (level1, '(/,"Number of target state symmetries: ",i5)') den_namelist % ntgt

            do isym = 1, den_namelist % ntgt
                do jsym = isym, den_namelist % ntgt
                    write (level1, '(/,"Properties for target state symmetries with indices: ",2I5)') isym, jsym
                    call cden_namelist % init(maxtgt)
                    call cden_namelist % construct_namelists_from_denprop_input(den_namelist, isym, jsym)
                    call cdenprop_drv(ci_vec_i, ci_vec_j, cden_namelist, mol_prop_data)
                    call cden_namelist % dealloc
                end do
            end do

            mol_prop_data % swintf_format = .true.

            call mol_prop_data % sort_by_energy
            call mol_prop_data % write_properties(lupropw_final, den_namelist % ukrmolp_ints)
            if (.not. mol_prop_data % contains_continuum) then
                call mol_prop_data % write_properties(stdout, den_namelist % ukrmolp_ints)
            end if

        end if

        if ( den_namelist % ipol /= 0 ) then
            call mol_prop_data % clean
            call mol_prop_data % read_properties(lupropw_final,1)
            call mol_prop_data % calculate_polarisability(stdout, den_namelist % qmoln)
        end if

    end do namelist_loop

    call mpi_mod_finalize

end program cdenprop_all
