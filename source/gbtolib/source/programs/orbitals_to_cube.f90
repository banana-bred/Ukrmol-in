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

!> \brief Write orbitals to Gaussian CUBE files
!>
!> The program reads the input namelist &cube from the standard input. A typical one can look like this:
!>
!> \verbatim
!>  &cube
!>
!>    filename = 'dyson_orbitals.ukrmolp',  ! or 'moints'
!>
!>    orbitals = 1,10,                      ! range of orbitals (default: 1,1)
!>    a = 15,                               ! inner region radius (default: 0)
!>
!>    nx = 201,  xmin = -10,  dx = 0.1,     ! x grid (default: 201, -10, 0.1)
!>    ny = 201,  ymin = -10,  dy = 0.1,     ! y grid (default: 201, -10, 0.1)
!>    nz = 201,  zmin = -10,  dz = 0.1,     ! z grid (default: 201, -10, 0.1)
!>
!>    verbosity = 2,                        ! verbosity level (default: 0)
!>
!>  /
!> \endverbatim
!>
!> The entry "filename" is a path to a GBTOlib data file containing an orbital basis, produced, e.g., by scatci_integrals (molecular
!> orbitals) or cdenprop (Dyson orbitals). The entry "orbitals" gives a range of orbitals to convert. These indices are global,
!> that is, they refer to orbitals in the concatenated array from all irreducible representations. The evaluated orbitals will be
!> written to files "orb_X.cube", where "X" is the orbital index. The remaining parameters define the Cartesian grid to use for
!> the CUBE files.
!>
!> When evaluating an orbital on grid, GBTOlib first evaluates all atomic orbitals on that grid (unless already done on the same
!> grid). This is a lot of data. To avoid excessive memory consumption, the evaluation is performed in cycles and only subsets
!> of the whole product grids are processed at one time.
!>
!> The molecular orbital expansion coefficients stored in the files are always given in terms of the Gaussians normalized to
!> the inner region. However, Gaussian orbitals stored in the files are stored in their original form. That requires normalizing
!> them prior to evaluation. For that, specification of the inner region radius is necessary. If this renormalization is not
!> needed (e.g. when there are no long-range Gaussians leaking out of the inner region), the radius can be omitted.
!>
program orbitals_to_cube

    use atomic_basis_gbl,     only: atomic_orbital_basis_obj
    use common_obj_gbl,       only: nucleus_type
    use const_gbl,            only: set_verbosity_level
    use molecular_basis_gbl,  only: molecular_orbital_basis_obj
    use mpi_gbl,              only: mpi_mod_start, mpi_mod_finalize
    use precisn_gbl,          only: cfp

    implicit none

    type(atomic_orbital_basis_obj), target :: atomic_orbital_basis
    type(molecular_orbital_basis_obj)      :: molecular_orbital_basis
    type(nucleus_type), pointer            :: atoms(:)

    real(cfp), allocatable :: r(:, :), orbital_at_r(:), evaluated_orbitals(:, :, :, :)
    integer,   allocatable :: sign_at_r(:)

    character(len=1000) :: filename

    real(cfp) :: a, x, y, z, xmin, ymin, zmin, dx, dy, dz
    integer   :: iorb, norbs, ipoint, npoints, ix, nx, iy, ny, iz, nz, iatom, natoms, u, orbitals(2), i, j, verbosity

    namelist /cube/ filename, orbitals, nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, a, verbosity

    ! -------------------------------------------------------------------- !
    ! 1. Initialize defaults and read the input namelist                   !
    ! -------------------------------------------------------------------- !

    a = 0

    nx = 201;  xmin = -10;  dx = 0.1
    ny = 201;  ymin = -10;  dy = 0.1
    nz = 201;  zmin = -10;  dz = 0.1

    orbitals = 1
    verbosity = 0

    read (*, nml=cube)

    call mpi_mod_start(.true.)
    call set_verbosity_level(verbosity)

    ! -------------------------------------------------------------------- !
    ! 2. Set up the orbital basis                                          !
    ! -------------------------------------------------------------------- !

    print '(3a,/)', 'Reading orbitals from file "', trim(filename), '"'

    molecular_orbital_basis % ao_basis => atomic_orbital_basis

    call atomic_orbital_basis % read(filename)
    call atomic_orbital_basis % normalize_continuum(a)
    call molecular_orbital_basis % read(filename)

    do i = 1, molecular_orbital_basis%no_irr
        norbs = molecular_orbital_basis % get_number_of_orbitals(i)
        print '(2x,a,i0,a,i0,a)', 'irr ', i, ' contains ', norbs, ' orbitals'
    end do

    atoms => atomic_orbital_basis % symmetry_data % nucleus
    natoms = count(atoms(:) % name /= 'sc')

    print '(/,a,/)', 'Nuclei'
    do iatom = 1, natoms
        print '(2x,a,3f8.3)', atoms(iatom) % name, atoms(iatom) % center(:)
    end do

    ! -------------------------------------------------------------------- !
    ! 3. Set up Cartesian evaluation grid                                  !
    ! -------------------------------------------------------------------- !

    ipoint = 0
    npoints = nx*ny*nz

    ! set up regular Cartesian grid for sampling the wave-function
    allocate (r(3, npoints))
    do ix = 1, nx
        x = xmin + (ix - 1)*dx
        do iy = 1, ny
            y = ymin + (iy - 1)*dy
            do iz = 1, nz
                z = zmin + (iz - 1)*dz
                ipoint = ipoint + 1
                r(:, ipoint) = [x, y, z]
            end do
        end do
    end do

    ! -------------------------------------------------------------------- !
    ! 4. Evaluate all orbitals on the Cartesian grid                       !
    ! -------------------------------------------------------------------- !

    allocate (evaluated_orbitals(nz, ny, nx, orbitals(1):orbitals(2)))

    print '(/,a,/)', 'Orbitals to evaluate'

    print '(2x,a10,a10,a10)', 'absolute', 'symmetry', 'relative'
    do iorb = orbitals(1), orbitals(2)
        print '(2x,i10,i10,i10)', iorb, &
            molecular_orbital_basis%get_orbital_symmetry(iorb), &
            molecular_orbital_basis%get_index_within_symmetry(iorb)
    end do

    print '(/,a,/)', 'Evaluate orbitals on grid'

    do ix = 1, nx
        ipoint = (ix - 1)*ny*nz
        print '(2x,*(a,i0))', 'evaluating at points ', ipoint + 1, '...', ipoint + ny*nz, ' of ', npoints
        do iorb = orbitals(1), orbitals(2)
            call molecular_orbital_basis%eval_orbital(iorb, r(:, ipoint + 1:ipoint + ny*nz), ny*nz, orbital_at_r, sign_at_r)
            evaluated_orbitals(:, :, ix, iorb) = reshape(orbital_at_r, [nz, ny])
        end do
    end do

    print '(/,a,/)', 'Write to files'

    do iorb = orbitals(1), orbitals(2)
        write (filename, '(a,i0,a)') 'orb_', iorb, '.cube'
        open (newunit = u, file = trim(filename), action = 'write', form = 'formatted')
        write (u, '(a)') 'Gaussian CUBE file'
        write (u, '(a,i0)') 'orb', iorb
        write (u, '(i5,3f13.6)') natoms, xmin, ymin, zmin
        write (u, '(i5,3f13.6)') nx, dx, 0._cfp, 0._cfp
        write (u, '(i5,3f13.6)') ny, 0._cfp, dy, 0._cfp
        write (u, '(i5,3f13.6)') ny, 0._cfp, 0._cfp, dz
        do iatom = 1, natoms
            write (u, '(i5,4f13.6)') nint(atoms(iatom)%charge), 0._cfp, atoms(iatom)%center(:)
        end do
        do ix = 1, nx
            do iy = 1, ny
                write (u, '(6e13.5)') evaluated_orbitals(:, iy, ix, iorb)
            end do
        end do
    end do

    ! -------------------------------------------------------------------- !
    ! 5. Finalize                                                          !
    ! -------------------------------------------------------------------- !

    call mpi_mod_finalize

end program orbitals_to_cube
