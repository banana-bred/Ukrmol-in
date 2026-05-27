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

!> \brief   MPI SCATCI Constants module
!> \authors A Al-Refaie
!> \date    2017
!>
!> This module provides constants used by MPI SCATCI.
!>
!> \note 30/01/2017 - Ahmed Al-Refaie: Initial revision.
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module consts_mpi_ci

    use precisn, only: shortint, longint, wp

    implicit none

    public

    !> Maximal length of integer arrays in input namelists
    integer, parameter :: mxint = 1000

    !> Number of long integers used to store up to 8 (packed) integral indices
    integer, parameter :: NIDX = 2

    !> Default limit on the number of eigenvalues printed
    integer, parameter :: MAX_NEIG = 400

    !> A default memory value (in GiB), we use the default given for per proc in archer
    real(wp), parameter :: DEFAULT_ARCHER_MEMORY = 2.5

    !> The maximum length of a name
    integer, parameter :: NAME_LEN_MAX = 120

    !------------------SYMMETRIES--------------------------!
    !>This describes C_inf_v symmetries
    integer, parameter :: SYMTYPE_CINFV = 0
    !>This describes D_inf_h symmetries
    integer, parameter :: SYMTYPE_DINFH = 1
    !>This describes D_2_h symmetries
    integer, parameter :: SYMTYPE_D2H = 2

    !-------CI TARGET_FLAG------------------------!
    !> No Ci target contraction (target only run)
    integer, parameter :: NO_CI_TARGET = 0
    !>Perform CI target contraction (target+scattering run)
    integer, parameter :: DO_CI_TARGET_CONTRACT = 1

    !------IEXPC FLAG
    !> Use configuration state functions as is
    integer, parameter :: NORMAL_CSF = 0
    !> The configuration state functions are prototypes and therefore require expansion
    integer, parameter :: PROTOTYPE_CSF = 1
    !> Generate target wavefunction parameter
    integer, parameter :: GENERATE_TARGET_WF = 0

    !-----MATRIX_EVALUATION---------------------!
    !> Matrix is evaluated as a difference from the first element
    integer, parameter :: MATRIX_EVALUATE_DIFF = 0
    !> Matrix is evaluated as full
    integer, parameter :: MATRIX_EVALUATE_FULL = 1
    !> Same as diff but dont evalute integrals belonging to the target state
    integer, parameter :: NO_PURE_TARGET_INTEGRAL_DIFF = 2
    !> Same as full but dont evalute integrals belonging to the target state
    integer, parameter :: NO_PURE_TARGET_INTEGRAL = 3

    !-----DIAGONALIZATION-----------------------!
    !> No diagonalization
    integer, parameter :: NO_DIAGONALIZATION = 0
    !> Diagonalize
    integer, parameter :: DIAGONALIZE = 1
    !> Diagonalize but with no restart
    integer, parameter :: DIAGONALIZE_NO_RESTART = 2
    !> Do not diagonalize, just read from file
    integer, parameter :: READ_FROM_FILE = 3

    !----PHASE_CORRECTION------------------------!
    !> No phase correction
    integer, parameter :: NORMAL_PHASE = -1
    !> Phase correct by negelcting continuum orbitals
    integer, parameter :: NEGLECT_CONTINUUM = 0

    !--------DIAGONALIZER CHOICE--------------------------!
    !> SCATCI chooses the diagonalizer
    integer, parameter :: SCATCI_DECIDE = -10
    !> Force a dense diagonalizer (all eigenvalues)
    integer, parameter :: DENSE_DIAG = 1
    !> Force a Davidson diagonalizer
    integer, parameter :: DAVIDSON_DIAG = 0
    !> Force an iterative (like ARPACK) diagonalizer
    integer, parameter :: ITERATIVE_DIAG = -1

    !> Omit L2 section of the eigenvectors (only keep continuum channels)
    integer, parameter :: SAVE_CONTINUUM_COEFFS = 1  ! b'001'
    !> Eigenpairs will be saved in memory and passed to CDENPROP and outer interface
    integer, parameter :: PASS_TO_CDENPROP = 2  ! b'010'
    !> Keep only L2 section of the eigenvectors (omit continuum channels)
    integer, parameter :: SAVE_L2_COEFFS = 4  ! b'100'
    !> Do not discard any coefficients
    integer, parameter :: SAVE_ALL_COEFFS = SAVE_CONTINUUM_COEFFS + SAVE_L2_COEFFS  ! = b'101' = 5

    !--------------------MATRIX CONSTANTS ---------------------------------!
    !> Matrix is dense
    integer, parameter  :: MAT_DENSE  = 0
    !> Matrix is sparse
    integer, parameter  :: MAT_SPARSE = 1
    !> Matrix element storage threshold
    real(wp), parameter :: DEFAULT_MATELEM_THRESHOLD = 1e-15
    !> Matrix uses C indexing
    integer, parameter  :: MATRIX_INDEX_C = 0
    !> Matrix uses FORTRAN indexing
    integer, parameter  :: MATRIX_INDEX_FORTRAN = 1

    !> Sort via IJ ordering
    integer, parameter  :: MATRIX_ORDER_IJ = 1
    !> Sort via JI ordering
    integer, parameter  :: MATRIX_ORDER_JI = 2
    !> Default expansion size for the cache
    integer, parameter  :: DEFAULT_EXPAND_SIZE = 100000

    !--------------------HAMILTONIAN CONSTANTS ---------------------------------!
    !> Hamiltonain is a scattering one
    integer, parameter  :: MAIN_HAMILTONIAN  = 0
    !> Hamiltonain is target one
    integer, parameter  :: TARGET_HAMILTONIAN = 1

    !> Default threshold to store integrals
    real(wp), parameter :: DEFAULT_INTEGRAL_THRESHOLD = 0

    !----------ECP CONSTANTS----------------------------------------------!
    integer, parameter :: ECP_TYPE_NULL = 0
    integer, parameter :: ECP_TYPE_MOLPRO = 1

    integer, parameter :: SLEPC_MAX_EIGENVALUES = 100 !This is for the spectrum slicing method

end module consts_mpi_ci
