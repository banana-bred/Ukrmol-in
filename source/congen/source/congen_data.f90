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
!> \brief   Module containing parameters / data required throughout the CONGEN program.
!> \authors J Carr
!> \date    2010
!>
!> This module defines many variables that for instance set hard-coded limits on various arrays used in the program
!> that are not expected to need dynamic allocation (e.g. maximal number of irreducible representation), or other
!> frequently used specific numbers (e.g. id of standard input unit or general math constants).
!>
!> Besides these parameters a lot of variables shared between CONGEN's routines are defined (and documented) here to
!> avoid extensive argument lists.
!>
module congen_data

    use precisn, only : wp
    use consts,  only : xhalf

    save

    integer, parameter :: lg = 75       !< ??? used to dimension arrays below; also in CPLE[A,M].
    integer, parameter :: nu = 20       !< Maximal number of irreducible representation (= max dimension of NOB etc.).
    integer, parameter :: ny = 150      !< Used to dimension REFGU and REFORB in csfgen.  Limit on input variable NREFO.
    integer, parameter :: nz = 150      !< Used to dimension REFDET in \ref congen_driver::csfgen; maximum number of electrons.
    integer, parameter :: jx = 75       !< Used to dimension many arrays in \ref congen_driver::csfgen; limit on the sum of elements in the NSHLP array.
    integer, parameter :: jy = 10       !< Used to dimension various arrays in csfgen including NELECP and NSHLP.
                                        !< Limit on input variable NDPROD.
    integer, parameter :: jz = 10       !< Maximum number of constraints upon CSFs accepted (limit on input variable NTCON).
    integer, parameter :: ky = 75       !< Used to dimension REFGUG and REFORG in congen_driver::csfgen.  Limit on input variable NREFOG?
    integer, parameter :: kz = 500      !< Used to dimension REFDTG in csfgen; limit on input variable NELECG?
    integer, parameter :: mxtarg = 300   !< Maximum number of target state symmetries.

    real(kind=wp), parameter :: root2 = 0.7071067811865475_wp   !< Note that this is actually inverse square root!
    real(kind=wp), parameter :: thresh1 = 1.E-10_wp             !< Threshold used in \ref congen_distribution::cgcoef.
    real(kind=wp), parameter :: thresh2 = 1.0E-30_wp            !< Threshold used in \ref congen_projection::wfgntr.

    integer, parameter :: nitem = 9     !< Number of interesting things printed per line in the PRINT* routines.
    integer, parameter :: nftr = 5      !< Unit number for reading standard input.
    integer            :: nftw = 6      !< Unit number for printing; may be changed via the STATE namelist so not a parameter.

    character(len=4), dimension(3), parameter :: blnk43 = &
        (/ '    ', '    ', '    ' /)
    character(len=8), dimension(7), parameter :: header = &
        (/ 'SHELL NO', 'OCCSHL  ', 'SYM     ', 'GUSHL   ', 'PQN     ', 'QNSHL   ', 'CUP     ' /)
    character(len=4), dimension(7), parameter :: rhead  = &
        (/ '    ', 'MSHL', 'PQN ', 'NES ', 'SO1 ', 'SO2 ', 'G/U ' /)

    character(len=3), parameter :: lp = ' ( '
    character(len=1), parameter :: star = '*'
 
    character(len=1), dimension(132) :: head    !< String with output page number and user-specific header text.
    integer :: lppmax                           !< Pagination limit - how many lines between two headers.
    integer :: lppr                             !< Global output line counter -- number of lines on current page.
    integer :: npage                            !< Page counter, currently limited to 3 digits (see \ref congen_pagecontrol::addl).

    integer :: lratio                           !< Ratio of real size to integer size. Used to manage workspace data.
    integer :: icdi                             !< Position in workspace where the determinant coefficients start (one coeff per one determinant).
    integer :: iexcon                           !< Position in workspace where EX constraints start.
    integer :: indi                             !< Position in workspace where packed determinants start (each as size + replacements).
    integer :: inodi                            !< Position in workspace where CSFs start (as number of determinants per CSF).
    integer :: iqnstr                           !< Starting position in workspace for per-shell quantum numbers (e.g. \ref congen_distribution::cplem).
    integer :: irfcon                           !< 
    integer :: ndel                             !< Workspace position used for the CSFs read from input.
    integer :: ndel1                            !< Workspace position used for reading CSFs from the input.
    integer :: ndel2                            !< Workspace position used for reading CSFs from the input.
    integer :: gutot                            !< Total gerade (= +1) / ungerade (= -1) quantum number.
    integer :: nisz                             !< Total S_z, set by CSFGEN, defaults to first element of \ref qntot.
    integer :: nndel                            !< Number of workspace elements used for reservation of space needed when reading CSFs from input.
    integer :: nnlecg                           !< The total number of electrons which are movable within the current wfn group.
    integer :: symtyp                           !< Symmetry type, 0 = C_infv, 1 = D_infh, 2 = {D_2h, C_2v, C_s, E}.

    integer, dimension(3,lg)   :: cup           !< Coupling scheme data. ???
    integer, dimension(3,lg)   :: pqnst         !< 
    integer, dimension(lg)     :: gushl         !< Gerade/ungerade per shell.
    integer, dimension(lg)     :: mshl          !< Symmetry number from zero to n-1 (mvalue).
    integer, dimension(lg)     :: occshl        !< Compressed shell occupations, all zeros deleted and pseudo shells expanded.
    integer, dimension(lg)     :: qnshlr        !< Total quantum numbers for all shells in the set of orbitals being processed.
    integer, dimension(lg)     :: spnmin        !<
    integer, dimension(lg)     :: sshlst        !<
    integer, dimension(2,lg)   :: kslim         !<
    integer, dimension(2,lg)   :: mclim         !<
    integer, dimension(3,2*lg) :: qnshl         !<
    integer, dimension(3)      :: qntar         !< Coupling control for current wfn group. ???
    integer, dimension(3)      :: qntot         !< Total quantum numbers.
    integer, dimension(nu)     :: shlmx1        !< Maximal occupancy of orbitals per symmetry.

    integer, parameter :: jsmax = 20                                    !< Maximal order of precomputed binomial coefficients.
    real(kind=wp), dimension((jsmax*jsmax + 3*jsmax + 4)/2) :: binom    !< Table of precomputed Pascal triangle.
    integer, dimension(jsmax+2) :: ind                                  !< Pascal triangle row pointers into \ref congen_data::binom.

    integer :: ift      !< Used to signal to \ref congen_distribution::assign to initialize some variables.
    integer :: nobt     !< Total number of all orbitals considered in the calculation, essentially equal to sum(nob).
    integer :: ntcon    !< 
    integer :: nsym     !< Number of symmetries given in input file (right-trimmed of those with zero orbitals).

    integer, dimension(nu)    :: nobi           !< Running index of the first orbital in each symmetry.
    integer, dimension(nu)    :: nshsym         !< 
    integer, dimension(jz)    :: nrcon          !< 
    integer, dimension(jz)    :: test           !< Determines which of the two types of constraints will be used.
    integer, dimension(nu,lg) :: nst            !< Pointer to shell index and order count; only used in \ref congen_distribution::assign.

    integer :: confpf   !< Print flag for the amount of print given of configurations and prototypes.
    integer :: ncsf     !< Total number of CSFs generated.
    integer :: ne       !< 
    integer :: nshl     !< Number of shells in the set of orbitals currently being processed.
    integer :: ntyp     !< Prototype number, updated in \ref congen_distribution::distrb.
    integer :: cdimx    !< Maximal number of determinants.
    integer :: ndist    !< Number of distributions generated from set of shells (set in \ref congen_distribution::assign).
    integer :: nstate   !< Number of couplings computed in \ref congen_distribution::cplea or \ref congen_distribution::cplem.
    integer :: megul    !< File unit for binary output of generated configurations.
    integer :: ncall    !< Used to signal to "wfn" to initialize some variables.
    integer :: ndimx    !< Maximal number of workspace elements usable for packed determinant data.
    integer :: nncsf    !< Number of CSFs generated by \ref congen_distribution::wfn (total).
    integer :: nodimx   !< Maximal number of workspace elements usable for CSF information.
    integer :: ntso     !< Total number of spin-orbitals, given the orbitals and their maximal occupacy (from group properties).

    integer, allocatable   :: exdet(:)          !< \todo Local to \ref congen_distribution::wfn - move it there!
    integer, allocatable   :: exref(:)          !< Population of spin-orbitals (0 = not populated, 1 = populated).
    integer, dimension(lg) :: nonew             !< 
    integer, dimension(lg) :: norep             !< 
    integer, dimension(nu) :: nsoi              !< Running index of the first spin-orbital within given total symmetry.

    integer :: iidis2                           !< 
    integer :: lcdi                             !< Total number of determinants (including those already flushed to disk).
    integer :: lndi                             !< Total number of integers forming packed determinants.
    integer :: ni                               !< Number of determinants currently in memory.
    integer :: nid                              !< Length of the integer array containing all packed determinants.
    integer :: noi                              !< Number of determinants per CSF.
    integer :: next                             !< Position in the workspace indicating the first unused element.
    integer :: nx                               !< Full or remaining workspace size (depending on context).
    integer :: idcp                             !< 
    integer :: idop                             !< 
    integer :: ieltp                            !< 

end module congen_data
