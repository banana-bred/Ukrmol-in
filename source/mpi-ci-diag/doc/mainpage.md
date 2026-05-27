\authors Ahmed Al-Refaie (University College London), Jakub Benda (Open University)
\tableofcontents

[TOC]

MPI-SCATCI serves for construction and diagonalization of the molecular Hamiltonian.

The program source is implemented in \ref mpi_scatci and uses a lot of modules that implement various molecular integral
formats and diagonalization approaches.

# Parallelization

MPI-SCATCI can be run as a serial application in the similar way as SCATCI,

    mpi-scatci < scatci.inp 1> scatci.out 2> scatci.err

but is also parallelized using MPI,

    mpiexec -n 4 mpi-scatci scatci.inp 1> scatci.out 2> scatci.err

For this purpose it uses the MPI module present in the UKRmol+ GBTO library. In the MPI mode, MPI-SCATCI can optionally use the
MPI-3 shared memory for storage of molecular integrals. This is automatically enabled when GBTOlib is compiled with MPI-3
support and it helps to reduce memory requirements of large calculations.


# Program operation

The program starts by initialization of the MPI module and redirection of non-master process outputs to log files
(named `log_file.X`, where X stands for the rank of the process). The output of the master process is not redirected
and goes to the standard output.

The input can be passed to MPI-SCATCI in the original SCATCI namelist format, or as a concatenation of several SCATCI inputs.
The module \ref Options_module then also reads the data from the binary output of CONGEN.

With the input read, all orbitals are constructed using the object \ref Orbital_module::OrbitalTable (based on original
SCATCI subroutines). These are then used for construction of all configurations by \ref CSF_module::CSFObject. First,
configurations are read from CONGEN binary output file. Then they are transformed to an internal bit field representation
used by MPI-SCATCI. In MPI-3 shared memory mode there is just one array with the CSF data, populated by the master process.
In non-shared memory mode all processes populate their own copies. Determinants are stored as bit fields and processed by
bit operations, which reduces the memory needed and generally offers better performance.

Next stage is the loading of molecular integrals. If we are dealing with non-abelian group (for diatomics), ALCHEMY integrals
are necessary (\ref ALCHEMY_module). Otherwise, the code uses either SWEDEN or UKRmol+ integrals, depending on the switch
used in the input configuration (\ref SWEDEN_module, \ref UKRMOL_module). Both ALCHEMY and SWEDEN reading routines are
implemented directly in MPI-SCATCI. UKRmol+ integrals, on the contrary, are read using the interface subroutines of the
GBTO library.

Two different Hamiltonian matrices can be diagonalized: uncontracted and contracted. Uncontracted Hamiltonian is a regular
Hamiltonian, typically used for calculation of eigenstates of the target (see \ref Uncontracted_Hamiltonian_module).
Contracted Hamiltonian (see \ref Contracted_Hamiltonian_module) is such Hamiltonian, where the so-called "CI target"
states are employed. This is used in "scattering", or (N + 1)-electron, runs. It is also possible to avoid any diagonalization
at all, forcing MPI-SCATCI to just write out the Hamiltonian to disk, together with information file *ham_data*
(\ref Options_module::ham_pars_io_wrapper) that can be used for further processing of the Hamiltonian elements.

If the contracted ((N+1)-electron) Hamiltonian diagonalization is requested, MPI-SCATCI first prepares the (N-electron)
target states (see \ref Target_RMat_CI_module). Again, depending on the user setup, these are either loaded from disk
(as a result of previous target-only, uncontracted diagonalization), or obtained by diagonalization on the fly;
the same method is used for this diagonalization as for the subsequent (N+1)-electron one. Target eigenstates
are then use to constrain the contracted Hamiltonian so that it only uses these as the "target part" of the eigenstates
containing continuum orbitals; for details see [1,2].

Once the diagonalization is done, eigen-energies and eigen-values are written to disk in the SCATCI `fort.25` format
using the object from \ref SolutionHandler_module. Most of the writing is done from inside the diagonalizer type
and the calculated eigenvectors are deallocated immediately, unless the calculation of derived quantities is requested.

Various input parameters allow MPI-SCATCI to perform further tasks, like exclude a specific row and column from
the Hamiltonian, pass the eigenstates to CDENPROP for calculation of transition dipole moments, or produce the outer channels
and boundary amplitudes for UKRmol-out. See \ref Postprocessing_module. Additionally, the outer interface has a basic automatic
support for the contracted virtual orbitals (if `auto_nvo` is set to `.true.` in the input namelist `outer_interface`):
When the number of continuum orbitals specified in GBTOlib (scatci_integrals) is
different than the number of CI target continuum orbitals specified in the input for MPI-SCATCI, the difference (i.e. target
orbitals included as fake continuum) are automatically recognized as contracted virtual orbitals. This is equivalent to the NVO
infrastructure in SWINTERF. If more fine-grained control over calculation of NVO is desired, use the legacy SCATCI together
with SWINTERF instead.


# Diagonalization methods

The available diagonalization methods are listed in the following table:

| Method                   | Available in serial | Available in parallel | Dependency    | Module                                |
| ------                   | ------------------- | --------------------- | ----------    | ------                                |
| direct (dense)           | yes                 |                       | LAPACK        | \ref LapackDiagonalizer_module        |
| direct (dense)           |                     | yes                   | ScaLAPACK     | \ref SCALAPACKDiagonalizer_module     |
| iterative (Davidson)     | yes                 |                       |               | \ref DavidsonDiagonalizer_module      |
| iterative (Davidson)     |                     | (unfinished)          | SLEPc + PetSc | \ref SLEPCDavidsonDiagonalizer_module |
| iterative (Arnoldi)      | yes                 |                       | Arpack        | \ref ARPACKDiagonalizer_module        |
| iterative (Krylov-Schur) |                     | yes                   | SLEPc + PetSc | \ref SLEPCDiagonalizer_module         |

As indicated in the table, at the present moment the parallel version of Davidson method is not available. MPI-SCATCI
uses Krylov-Schur method instead.


# Input namelist options

The following table lists all members of the namelist `&input`:

| Option         | Default   | Description                                                               |
| -------------- | --------- | ------------------------------------------------------------------------- |
| `icidg`        | `1`       | Diagonalize: `0` = no, `1` = yes, `2` = yes (no restart), `3` = read only.|
| `icitg`        | `0`       | Use contracted target: `0` = no (target run), `1` = yes (scattering run). |
| `iexpc`        | `0`       | Assume that CSFs contains continuum prototypes: `0` = no, `1` = yes.      |
| `idiag`        | `-1`      | Hamiltonian matrix evaluation flag. See below for details.                |
| `nfti`         | `16`      | File unit with molecular integrals.                                       |
| `nfte`         | `26`      | Hamiltonian scratch file unit.                                            |
| `lembf`        | `5000`    | Number of matrix elements per record in the Hamiltonian scratch file.     |
| `megul`        | `13`      | File unit with CSFs from CONGEN.                                          |
| `nftg`         | `0`       | File unit with target eigenstates for phase consistency.                  |
| `ntgsym`       | `0`       | Number of different spin-symmetries for target eigenstates.               |
| `numtgt`       | `0,0,...` | Number of target eigenstates per spin-symmetry.                           |
| `notgt`        | `0,0,...` | Number of continuum orbitals used with given target spin-symmetry.        |
| `nctgt`        | `1,1,...` | Number of CI components per target spin-symmetry. (Not used in UKRmol+.)  |
| `ntgtf`        | `0,0,...` | Input dataset in `nftg` per target eigenstate.                            |
| `ntgts`        | `0,0,...` | Position of the eigenstate in the input dataset.                          |
| `mcont`        | `0,0,...` | Continuum M-value (spatial symmetry) per target eigenstate.               |
| `ncorb`        | `-1`      | Phase correction by neglecting continuum: `-1` = no, `0` = yes.           |
| `iord`         | `1`       | ALCHEMY molecular integrals ordering flag. (Not used in UKRmol+.)         |
| `name`         |           | Textual description of the run.                                           |
| `scalem`       | `0`       | Exotic mass flag. (Not used in UKRmol+.)                                  |
| `scfuse`       | `0`       | ALCHEMY molecular integrals control flag. (Not used in UKRmol+.)          |
| `qmoln`        | `.false.` | Write additional data for use in Quantemol-EC.                            |
| `ukrmolp_ints` | `.true.`  | Assume UKRmol+ molecular integral format.                                 |
| `iposit`       | `0`       | Scattering particle flag. `0` = electron, `-1` = positron                 |
| `enh_factor`   | `1.0D0`   | Positron enhancement factor. Rectifies missing polarisation description.  |

The Hamiltonian matrix evaluation flag `idiag` can have the following values:
 - 0: Matrix is evaluated as a difference from the first element.
 - 1: Matrix is evaluated as full.
 - 2: Same as 0 but don't evaluate integrals belonging to the target state.
 - 3: Same as 1 but don't evaluate integrals belonging to the target state.

Typically, one uses `idiag = 0` (default) for target runs, but `idiag = 1` (or `idiag = 2`) for scattering.

When `icidg = 3` is specified, the Hamiltonian is neither constructed nor diagonalized, but an already existing eigensolution
is read from the file and dataset sepcified in `&cinorn`. This is particularly useful when re-processing existing eigensolutions
with the built-in outer interface (see below). The configuration specified by the input must be identical to the one used
for calculation of the results to be read.

The following table lists all members of the namelist `&cinorn`:

| Option      | Default     | Description                                                                 |
| ----------- | ----------- | --------------------------------------------------------------------------- |
| `nftw`      | `25`        | Output unit for eigenpairs.                                                 |
| `nciset`    | `1`         | Dataset to write in the output file `nftw`.                                 |
| `npcvc`     | `1`         | Print dataset header when reading through the `nftw` file: 0 = no, 1 = yes. |
| `name`      |             | Textual description of the run.                                             |
| `nstat`     | `0`         | Number of eigenstates to obtain.                                            |
| `eshift`    | `0.,0.,...` | Energy shifts to be applied to target state (in atomic units)               |
| `igh`       | `-10`       | Algorithm to use: `0` = Davidson, `-1` = Arpack, `1` = dense.               |
| `critc`     | `1e-10`     | Arpack control tolerance.                                                   |
| `critr`     | `1e-8`      | Arpack control tolerance.                                                   |
| `ortho`     | `1e-7`      | Arpack control tolerance.                                                   |
| `crite`     | `1e-12`     | Arpack control tolerance.                                                   |
| `maxiter`   | `-1`        | Limit for iterative diagonalizers (Davidson, Arpack, SLEPc).                |
| `memp`      | `2.5`       | Memory limit per process in GiB.                                            |
| `forse`     | `0`         | Force serial diagonalizer even in parallel mode. (For debugging.)           |
| `exrc`      | `-1`        | Exclude row & column. (For debugging.)                                      |
| `vecstore`  | `5`         | Eigenvector storage and processing flag, see below.                         |
| `ecp`       | `0`         | Effective core potential to use. (Not fully implemented.)                   |

The keyword `vecstore` conveys more information than a single true/false flag. Its value expected
to be the sum of a subset of the following flags:

| Flag                   | Value |
| ----                   | ----- |
| store continuum coeffs | 1     |
| process in CDENPROP    | 2     |
| store L2 coeffs        | 4     |

The default value (5 = 4 + 1) tells MPI-SCATCI to store all coefficients (as done by SCATCI) and to
not include that particular spin-symmetry in CDENPROP post-processing.

The keyword `eshift` accepts a list of energy shifts which are applied to the eigenvalues in numerical order
starting from the lowest (most negative) state. The shifts should be given in atomic units.

For example, if a calculation obtains energies (in a given symmetry e.g., A1)

    E_calc = -0.490,-0.220,-0.075

But the reference data is

    E_ref  = -0.500,-0.250,-0.125

then `eshift` could be set as follows:

    eshift =  0.000,-0.020,-0.040

Note that the ground state is not shifted. Only the vertical excitation energies are shifted. Shifting the ground state would
change the vertical excitation thresholds for all of the other states.

The following table lists all members of the namelist `&outer_interface`:

| Option      | Default   | Description                                                         |
| ----------- | --------- | ------------------------------------------------------------------- |
| `write_amp` | `.true.`  | Write outer channel information and raw boundary amplitudes.        |
| `write_dip` | `.false.` | Calculate and write transition multipole moments.                   |
| `write_rmt` | `.false.` | Compose the input file for RMT; see [3] and [4].                    |
| `all_props` | `.true.`  | Compute also IRR-with-itself properties and sort states by energy.  |
| `auto_nvo`  | `.false.` | Automatically determine contracted virtual orbitals (`nvo`).        |
| `luprop`    | `-1`      | Property file for CDENPROP if different than `nfti` (for all IRRs). |
| `lupropw`   | `624`     | Output unit for calculated properties if `write_dip = .true.`.      |
| `lutdip`    | `24`      | Target transition dipoles for CDENPROP (typically with ISW = 1).    |
| `lutarg`    | `24`      | Target properties for SWINTERF (typically with ISW = 0).            |
| `itarget_spin` | `2,...`| Spin multiplicities of target spin-symmetries used in CDENPROP.     |
| `isw`       | `0`       | Produce electronic-nuclear props (= 0) or purely electronic (= 1).  |
| `nfdm`      | `18`      | Number of in-sphere evaluation radii for RMT data.                  |
| `delta_r`   | `0.08`    | Spacing of in-sphere evaluation radii for RMT data.                 |
| `rmatr`     | `-1.0`    | R-matrix radius.                                                    |
| `ntarg`     | `0`       | When non-zero, the number of targets to use in the outer region.    |
| `idtarg`    | `0,...`   | Comma-separated list of target IDs to use in the outer regions.     |
| `luchan`    | `10`      | Output unit for outer channel information.                          |
| `lurmt`     | `21`      | Output unit for raw boundary amplitudes.                            |
| `cform`     | '`F`'     | Formatted/unformatted output flag for the channel info file.        |
| `rform`     | '`U`'     | Formatted/unformatted output flag for the boundary amplitudes file. |

When not specified, the options `ntarg` and `idtarg` are automatically obtained from the the other inputs and from the target
properties file.

As always, namelist keywords are case-insensitive.


# CDENPROP post-processing

When evaluation of properties (multipole moments) are desired between states of some irreducible representations,
one needs to add the switch

    &cinorn
        ...
        vecstore = 2,  ! or 7 (= 1 + 2 + 4) to also write all eigenvectors
        ...
    /

The property file (in DENPROP format) will be written to *fort.624*. However, there are some important **caveats** to point out:
In case of a scattering diagonalization, CDENPROP assumes that the target eigenvectors are present in file *fort.26*
and---presently---there is no option to change this from MPI-SCATCI. However, *fort.26* is also the default Hamiltonian output
file name. So, to avoid collision, in a scattering diagonalization one needs to use a different unit for the Hamiltonian
by defining the `nfte` keyword to some available number, e.g.

    &input
        ...
        nfte = 80,
        ...
    /

Also, in the scattering mode, CDENPROP requires the target properties (*fort.24*), so these need to be present, too.
When `all_props` is set to `.true.` (default), properties are calculated between all eigenvector sets, including
elements between the eigenvectors of the same Hamiltonian. This is necessary for RMT or when using CDENPROP stage
in place of DENPROP. However, when running a photoionization calculation, only matrix elements between eigenstates
of distinct Hamiltonians are required. This can be achieved by setting `all_props` to `.false.`, which will also
inhibit automatic sorting of all eigenstates by energy when writing them to the property file.


# Advanced tuning options

MPI-SCATCI distributes parallel matrices using the standard ScaLAPACK two-dimensional block-cyclic model. The individual
block size is by default 64 by 64. While optimal in most cases, this may result in unnecessarily fine granularity when
dealing with huge Hamiltonians (100k by 100k and above), manifesting by increased communication or increased memory
consumption in parallel IO. The default value of 64 can be then increased by the environment variable
`CDENPROP_SCALAPACK_BLOCK_SIZE`, as follows:

    export CDENPROP_SCALAPACK_BLOCK_SIZE=256


# References

[1] Tennyson J, *A new algorithm for Hamiltonian matrix construction in electron–molecule collision calculations*,
J. Phys. B: At. Mol. Opt. Phys. **29** (1996) 1817–1828.

[2] Al-Refaie A F, Tennyson J, *A parallel algorithm for Hamiltonian matrix construction in electron–molecule collision calculations: MPI-SCATCI*,
Comput. Phys. Commun. **221** (2017) 53–62.

[3] Brown A C et al., *RMT: R-matrix with time-dependence. Solving the semi-relativistic, time-dependent Schrodinger equation for general,
multi-electron atoms and molecules in intense, ultrashort, arbitrarily polarized laser pulses*, Comput. Phys. Commun. **250** (2020) 107062.

[4] Mašín Z et al., *UKRmol+: a suite for modelling of electronic processes in molecules interacting with electrons, positrons and photons
using the R-matrix method*, Comput. Phys. Commun. **249** (2020) 107092.

