UKRmol+ (inner)
===============

UKRmol-in is licensed with the Gnu Public License.  For a comprehensive list of the developers 
that contributed to these codes see [Zenodo](https://doi.org/10.5281/zenodo.2630495).

Before compiling and using the code, check [Zenodo](https://doi.org/10.5281/zenodo.2630495) for the latest 
version of UKRmol-in, that may contain bug fixes, new functionality, etc.

Table of contents
-----------------

 * Downloading
 * Building
     * Conda
     * Docker
     * Building with CMake
     * Compatible libraries
     * List of available preprocessor directives
     * Testing
 * Notes for developers
     * Git-flow
     * Working on your own developments


Downloading
-----------

You can obtain a release version from Zenodo following the link above. If you are planning to develop the code,
you can request access and download the current development version from the [GitLab repository](https://gitlab.com/Uk-amor/UKRMol/UKRmol-in/).

UKRmol+ uses Git as its version control system. To download it, request permission from the
project managers and then just issue the following command:

    git clone --recurse-submodules https://gitlab.com/UK-AMOR/UKRMol/UKRmol-in.git

You will be prompted for your user name and password, and then the master branch of the repository 
will be downloaded to your computer. The option `--recurse-submodules` is required, because
the project uses several other projects using the Git submodules feature.

Some versions of Git will not prompt you for your user name and automatically infer it from
the environment (username). In such cases you may want to force your username. This can
be done by inserting `username@` between `https://` and `gitlab.com`, where `username` is your
Gitlab username.

If you forget to use the option `--recurse-submodules`, you can retrieve the submodule later by the sequence
of commands

    git submodule init
    git submodule update


Building
--------

### Building options:
* Conda/Mamba: portable, reproducible toolchain. Two solver setups are supported:
  - ARPACK + PETSc/SLEPc
  - ScaLAPACK + ELPA
* Docker: containerized build & run environment (Ubuntu-based)
* CMake: manual build using your own toolchain and compatible libraries

### Conda/Mamba/Docker
Detailed instructions for building UKRmol+ using conda/mamba/docker can be found in the
[conda/README.md](conda/README.md) file.

### Building with CMake (manual option)

UKRmol+ uses the freely and open source tool [CMake](https://cmake.org/download/) for
configuration, compilation, testing and installation. The newest CMake can always be
downloaded as an unpack-and-go binary package.

For building and testing, it is recommended to create a separate build directory (anywhere) and carry out
the compilation in that directory. (This way, one can have different builds with various combinations
of precision, integer width, compilers etc.) When using the Intel oneAPI 2025 programming suite,
the build process can be condensed to the following commands

    export CC=$(which icx)                    # serial C compiler is enough
    export FC=$(which mpiifx)                 # only parallel Fortran compiler is needed
    export FFLAGS="-i8 -qmkl-ilp64=cluster"   # link to nontreaded MKL and promote integers to 8 bytes (ILP64 mode)
    export BLA_VENDOR=Intel10_64ilp           # choose the ILP64 MKL BLAS and LAPACK
    
    cmake ..
    cmake --build .

with the two dots replaced by the path to the UKRmol-in directory. This will compile UKRmol+ with Intel MPI
and with Intel MKL (including BLAS, LAPACK and ScaLAPACK). Note particularly the compiler
flag `-i8` that forces use of long integers in UKRmol+. This is vital in all but the tiniest
calculations. Alternatively, the environment variables can be passed directly to CMake on the command line:

    cmake -D CMAKE_C_COMPILER=$(which icx) \
          -D CMAKE_Fortran_COMPILER=$(which mpiifx) \
          -D CMAKE_Fortran_FLAGS="-i8 -qmkl-ilp64=cluster" \
          -D BLA_VENDOR=Intel10_64ilp \
          ..

    cmake --build .

UKRmol+ supports many compilers and their versions. See the file [BUILDING.md](BUILDING.md) for more details
about other options and caveats.

### Compatible libraries

While UKRmol+ itself has to be compiled with 8-byte integers (i.e., in ILP64 mode), the external libraries can
use either 4-byte or 8-byte integers. However, to support very large calculations, using 8-byte integer
linear algebra libraries (BLAS, LAPACK, ScaLAPACK, Arpack, PetSc, SLEPc) may be necessary.

Here are details on the individual libraries used by UKRmol+. The needed ones are:

* **BLAS**, **LAPACK**: Intel MKL provides ILP64 version of both. When using open-source tools,
  it is possible to use OpenBLAS instead, compiled with option `INTERFACE64=1`. In spite
  of its name, that library provides also LAPACK.

The optional ones are:

* **ARPACK**: Can be compiled in 64-bit mode if `-DINTERFACE64=ON` option is passed to CMake
  when compiling ARPACK. But then also BLAS and LAPACK (which are needed by ARPACK) must
  use ILP64 interface. The Arpack library can be downloaded from https://github.com/opencollab/arpack-ng/tags .

* **SLEPc**: Used by MPI-SCATCI only; the suite will compile without them but MPI-SCATCI will not
  have its full functionality. To enable SLEPc, `-D SLEPC_INCLUDE_DIRS=...` should be set to a list
  of PetSc and SLEPc include directories with Fortran modules and `-D SLEPC_LIBRARIES=...` should be set to
  a list of PetSc and SLEPc libraries. See [BUILDING.md](BUILDING.md) for further details.

* **ScaLAPACK**: Intel MKL provides ILP64 version of this library. When using
  open-source tools, one can use the Reference ScaLAPACK (https://github.com/Reference-ScaLAPACK/scalapack),
  which also supports compilation in ILP64 mode. Some compilers (e.g. Cray) link ScaLAPACK automatically
  and no `SCALAPACK_LIBRARIES` need to be speficied. In such cases, add `-D WITH_SCALAPACK=ON` to your
  CMake command line instead, so that the configuration script knows that it should build ScaLAPACK-dependent
  code even when there are no `SCALAPACK_LIBRARIES`.

* **ELPA**: The "Eigenvalue SoLvers for Petaflop Applications" package (https://elpa.mpcdf.mpg.de/)
  can be optionally used by
  MPI-SCATCI when the module include directory and the ELPA library are provided using the CMake
  options `ELPA_INCLUDE_DIRS` and `ELPA_LIBRARIES`, respectively. It typically results in faster
  diagonalizations.

* **MPI**: Can be both LP64 and ILP64; this is automatically detected by the CMake configuration
  script. In contrast to the linear algebra libraries, 4-byte integer MPI library interface is not
  limiting for UKRmol+ and can be used even for very large calculations.

It is also possible to compile the code without the need for MPI at all. To achieve this: 
(i) avoid the compiler option `-Dusempi`; (ii) use the plain Fortran compiler, without the MPI
wrapper suggested in the above example (e.g. use ifx instead of mpiifx); (iii) add the CMake
option `-D WITH_MPI=OFF`.

This release of UKRmol+ is known to build and run correctly with the following software:

 - Arpack-ng 3.9.1
 - CMake 4.2.0
 - ELPA 2025.01.002
 - GNU Compiler Collection 15.2
 - Intel oneAPI 2025.3 (Intel IFX, Intel MKL, Intel MPI)
 - OpenBLAS 0.3.30
 - Open MPI 5.0.9
 - PetSc 3.24.2, SLEPc 3.24.1
 - Reference ScaLAPACK 2.2.2

See the file [BUILDING.md](BUILDING.md) for specifics of individual libraries used by UKRmol+.


### List of available preprocessor directives

GBTOlib (the library used by scatci_integrals) can use either double precision (default) or quadruple 
precision real numbers. To use the latter, add the following preprocessor definition among your compiler options (`CMAKE_Fortran_FLAGS`):

* `-Dusequadprec`  
  Enable calculation of molecular integrals in quadruple precision. This will make the
  operation notably slower, but will allow you to extend the Gaussian continuum basis to
  larger distances.
* `-Dquadpreckind=X`  
  If `-Dusequadprec` is used, this additional switch may be used to define the compiler-dependent Fortran
  kind of quadruple-precision real numbers to use. Default is true quadruple precision, often represented
  by the kind "16". Some compilers also recognize the kind "10" (or equivalent), that corresponds to 10-byte
  extended precision numbers.

There are several other options that can be added, all of them are related to capabilities of
the MPI library used. The CMake script will normally determine which of them to use on its
own, so you need to worry about them only in case that the automatic analysis fails.
The automatically included options are added to the CMake variable `GBTOlib_Fortran_FLAGS`, which
is printed out during the configuration step. The configuration script compiles a few trivial MPI
programs and runs them using the MPI launcher, which is expected to be named `mpiexec`. If you
use a different launcher (e.g. `mpirun` or `aprun`), you need to provide that name to CMake
using the option `-D MPIEXEC_EXECUTABLE=$(which mpirun)`.

* `-Dusempi`: include to compile with MPI library.

* `-Dquadreduceworks`: include if you're compiling with quad precision and your MPI library correctly handles
 quad precision arithmetics. Ignored if you're not using MPI or double precision is used.

* `-Dsplitreduce`: include if MPI_REDUCE from your library does not work correctly for data structures
 larger than ~2GB (i.e. limit of 32bit integer address). Ignored if you're not using MPI.

* `-Dmpithree`: include if you're using MPI 3 standard. This feature only affects the module
 mpi_memory_mod.F90 which is used by the UKRmol+ program MPI-SCATCI.


### Testing

After the compilation (`make`) finishes, the binaries can be tested using (within the build directory)

    ctest -R serial
    mv ./Testing Testing_serial
    ctest -R parallel
    mv ./Testing Testing_parallel

This will run all preset test calculations distributed with the suite and save the outputs of each type of test
to a separate directory for a further inspection. All the parallel runs use 2 MPI processes.
To provide additional flags for your MPI launcher in these runs, give them as a semicolon-separated list
`MPIEXEC_PREFLAGS` when invoking CMake before executing ctest.
For instance, with Open MPI you may want to restrict the number of OpenMP threads and prevent automatic binding 
of the processes like this

    cmake [... all options as before ...] -D MPIEXEC_PREFLAGS="-x;OMP_NUM_THREADS=2;--bind-to;none"

You can also execute only a subset of the test suite. For instance, this command will select only tests whose name contains
the string "target_CAS" and "_serial":

    ctest -R 'target_CAS.*_serial'

A full list of test names can be obtained by `ctest -N`. To also enable  UKRmol-out-related tests, one needs to build
UKRmol-out jointly with UKRmol-in, which is done by passing the path to the UKRmol-out package to CMake using the
option `-D UKRMOL_OUT_DIR=`.

For more details see the test suite [README](tests/suite/README.md).


Notes for the developers
------------------------

A commmit.template is available to ensure that developers do not forget to provide sufficient
information on thei work they've done. We ask that you do:

    git config commit.template commit.template

This way, when you commit some changes, the commit log will be pre-filled with the information in the
template. Alternatively, you can supply it on every commit by

    git commit -t commit.template


### Git-flow

From version 1.1 the UKRmol-in repository has switched to using  an adapted version of
[Git-flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).
Git-flow defines a set of branches and prefixes which have a well-defined meaning.

In the case of UKRmol-in the set of permanent branches and the naming convention is the following
(the flow chart [here](https://gitlab.com/Uk-amor/UKRMol/UKRmol-in/wikis/home) helps understanding
the structure):

    master : main development branch from which feature branches are created
    release: branch that contains code ready for (or already) release
    release-X.Y.Z: hotfix of the tagged version X.Y of the release branch

* Release is the branch which contains versions of the code that have been thoroughly tested. It usually
  correspond to a code that has undergone a major upgrade and is ready for release. Each version will
  be tagged, either by going to the Tags area of the GitLab interface, or automatically when a release file
  is generated using the Release facility of GitLab.

* Fixes of important bugs in the release versions (i.e. hotfixes) should be implemented in a branch
  originating from `release` and equipped with an incremented version number `release-X.Y.Z`, where X and Y
  correspond to the latest release tag. Implementation of each hotfix must be followed by its merge into
  the `master` branch.

* The Master branch is an integration branch for the Feature branches and it is based on the latest
  `release` version. We can see that it has the meaning of the `develop` branch in the usual
  Git-flow system, i.e. the `master` branch is the one into which all features are merged. This branch will
  eventually be merged into the `release` branch to create a new release.

* New features and your own developments should be implemented branching from the `master`. The names of
  such branches must be functional, i.e. referring to a particular feature that you're developing.
  If you wish you can prefix the name of your development branch with the string `feature-`.
  If it is possible to split the development of a major feature into a series of smaller upgrades
  (branches) then please do that since it helps to visualize the progress of the development and to
  see what was done by who and when. Implementation of each feature is followed by its merge into
  the `master` branch.


### Working on your own developments ###

To create a branch for your own developments use the following commands:

    git checkout master
    git branch name_of_feature

You'll then checkout the branch, work on it and commit it when needed. Once you're satisfied
with your implementation you should create a merge request by using the facility in the GitLab
interface. One of the managers/maintaners of the project will check it and merge it to `master`.

(When the development of one or several new features has finished and the code in the `master`
branch been thoroughly tested, the codes will be merged into the  `release` branch  and tagged.)
