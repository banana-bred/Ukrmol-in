GBTOlib
=======


Zdeněk Mašín, Jimena Gorfinkiel, 2012 - 2019

Jakub Benda, 2018 - 2019

GBTOlib logo: Copyright (C) 2018  Katarzyna Krzyzanowska

Table of contents
-----------------

 * Capability
 * Downloading
 * Building 
    * List of available preprocessor directives
 * Generating input for GBTOlib
 * Running `scatci_integrals`
 * Interfacing GBTOlib with other codes
 * Notes for developers
    * Git-flow
    * Working on your own developments
    * Other


Capability
----------

* Evaluation of molecular integrals in the basis of atom-centered Gaussian orbitals
  and center-of-mass centered B-splines and/or Gaussians.

* The center-of-mass (i.e. continuum) basis can be built from either Gaussians, B-splines or
  a combination of the two.

* An arbitrary angular momentum for the continuum basis.

* Transformation of the atomic integrals into basis of molecular integrals.

* Orbital orthogonalization using Gramm-Schmidt and/or Symmetric orthogonalization.

* Flexible configuration of the library allowing it to be ran on various machines ranging
  from single-node workstations to massively-parallel HPC architectures.

* Possibility to configure the library to use quadruple precision arithmetics.

* Input of molecular geometry, Gaussian basis and orbitals via the standard 
  [MOLDEN](http://cheminf.cmbi.ru.nl/molden/molden_format.html) file.


Downloading
-----------

GBTOlib uses Git as its version control system. To download it, request permission from the
project managers and then just issue the following command:

    git clone https://gitlab.com/UK-AMOR/UKRMol/GBTOlib.git

You will be prompted for your user name and password, and then the master branch of the repository 
will be downloaded to your computer.

Some versions of Git will not prompt you for your user name and automatically infer it from
the environment (username). In such cases you may need to force your username. This can
be done by inserting `username@` between `https://` and `gitlab.com`, where `username` is your
Gitlab username.


Building
--------

### Building options

**Conda/Mamba:** portable, reproducible toolchain.  
**Docker:** containerized build & run environment (Ubuntu-based).  
**CMake:** manual build using your own toolchain and compatible libraries.

---

### Conda/Mamba/Docker

Detailed instructions for building **GBTOlib** using **conda**, **mamba**, or **docker** can be found in the  
[`conda/README.md`](conda/README.md) file.

---

### Building manually (advanced option)

GBTOlib uses CMake configuration / build system. It is recommended to create a subdirectory
"build" in the same directory as this "README.md" file and carry out the compilation in that
subdirectory. (This way, one can have different builds with various combinations of precision,
integer width, compilers etc.) In an ideal world, you would just type the following commands:


    mkdir build; cd build
    cmake ..
    make

However, this would work only when your system environment was perfectly set up. In reality,
in most cases you will want to customize the build by specifying which compiler to use, and
where the needed external libraries are located. For that purpose, customize and use the sample
build scripts included below.  These example use the Intel Fortran compiler, with
Intel MPI and the BLAS and LAPACK routines from Intel MKL. You will need to adjust
the paths to your system (see below for an explanation of the libraries and flags used).For
CMake version older than 3.13:


    MKL_BLAS_LAPACK="$MKLROOT/lib/intel64/libmkl_intel_ilp64.so"
    MKL_BLAS_LAPACK="$MKLROOT/lib/intel64/libmkl_sequential.so;$MKL_BLAS_LAPACK"
    MKL_BLAS_LAPACK="$MKLROOT/lib/intel64/libmkl_core.so;$MKL_BLAS_LAPACK"
    
    cmake -D CMAKE_C_COMPILER=$(which icx) \
          -D CMAKE_Fortran_COMPILER=$(which mpiifx) \
          -D CMAKE_Fortran_FLAGS="-i8" \
          -D BLAS_LIBRARIES=$MKL_BLAS_LAPACK \
          -D LAPACK_LIBRARIES=$MKL_BLAS_LAPACK \
          ..
    
    make


Or, if you have CMake 3.13 or newer (recommended), you can use:


    export BLA_VENDOR=Intel10_64ilp

    cmake -D CMAKE_C_COMPILER=$(which icx) \
          -D CMAKE_Fortran_COMPILER=$(which mpiifx) \
          -D CMAKE_Fortran_FLAGS="-qmkl -i8" \
          ..

    make


(BLA_VENDOR is required when using the Intel compiler; for other compilers it is normally safe not to specify this.
The flag `-qmkl` ensures linking to the BLAS, LAPACK, BLAS95 and LAPACK95 libraries).

The scripts will compile all components, placing the resulting programs into the "bin" subdirectory
of your "build" directory. (It is possible to compile scatci_integrals using 'make scatci_integrals' 
provided  cmake has already been executed and the GBTOlib has been compiled). 

In its current state,  GBTOlib has to be compiled in ILP64 mode, i.e. with the use of 8-byte 
integers, bacause number of basis elements processed by the GBTO library can very easily reach 
ten-digit amounts. In the case of GNU Fortran, you need to use `-fdefault-integer-8`, in the case 
of Intel Fortran  use `-i8` as a compiler option to achieve the correct integer width. The external
libraries (BLAS, LAPACK, MPI) also may---but need not to---have 8-byte integer interface. GBTOlib
determines the actual interface during the CMake configuration step.

* **BLAS**, **LAPACK**: Intel MKL provides ILP64 version of both. When using open-source tools,
  it is possible to use OpenBLAS instead, compiled with option `INTERFACE64=1`. In spite
  of its name, that library provides also LAPACK.

* **MPI**: Can be both LP64 and ILP64. When using open-source tools, you can use OpenMPI,
  which by default compiles to LP64 mode, but can be compiled with ILP64 Fortran interface
  if `FCFLAGS=-fdefault-integer-8` is used (or `FCFLAGS=-i8` when compiling with/for
  Intel suite).

It is also possible to compile the code without the need for MPI at all. To achieve that: (i)
add `-D WITH_MPI=OFF` to the cmake command line; (ii)  use the plain Fortran compiler, without
the MPI wrapper suggested in the above example.

GBTOlib strongly benefits from automatic vectorization. The best performance is achieved by compiling
for the actual architecture where the code will run, so that the maximal supported version of AVX is
employed. This is achieved by including the corresponding flags in `CMAKE_Fortran_FLAGS`
(Intel: `-xHOST`, GNU: `-march=native -ffast-math`). Particularly with GNU Fortran this can lead even
to severalfold speed-up.



### List of available preprocessor directives ###

GBTOlib can use either double precision (default) or quadruple precision real numbers. To use the
latter, add the following preprocessor definition among your compiler options (`CMAKE_Fortran_FLAGS`):

* `-Dusequadprec`  
  Enable quadruple precision floating point operation. This will make the
  operation notably slower, but will allow you to extend the Gaussian continuum basis to
  larger distances.

* `-Dquadpreckind=N`  
  When used together with `-Dusequadprec`, overrides the floating-point "kind" used to represent
  the quadruple precision numbers to the provided compiler-specific value `N`. The default GBTOlib
  configuration uses 16-byte numbers for quadruple precision. However, some compilers may offer
  also other extended-precision floating-point numbers that provide lower precision at higher speed.
  For example, gfortran will use hardware-accelerated 10-byte floating point numbers (providing
  three extra digits of precision compared to the standard double precision) when `N=10`.

There are several other options that can be added, all of them are related to capabilities of
the MPI library used. The CMake script will normally determine which of them to use on its
own, so you need to worry about them only in case that the automatic analysis fails.
The automatically included options are added to the CMake variable `GBTOlib_Fortran_FLAGS`, which
is printed out during the configuration step. The configuration script compiles a few trivial MPI
programs and runs them using the MPI launcher, which is expected to be named `mpiexec`. If you
use a different launcher (e.g. `mpirun` or `aprun`), you need to provide that name to CMake
using the option `-D MPIEXEC_EXECUTABLE=$(which mpirun)`.

If you want or have to avoid the automatic analysis altogether (for example on a cluster frontend
that blocks execution of parallel applications), you can provide the MPI capabilities yourself
on the CMake command line using the `GBTOlib_Fortran_FLAGS` variable, for instance:

    -D GBTOlib_Fortran_FLAGS="-Dusempi;-Dsplitreduce;-Dmpi64bitinteger"

The full list of these MPI-related options follows:

* `-Dusempi`  
 Include if you want to compile with MPI library. Without this flag the library will
 be parallelized using OpenMP only.

* `-Dquadreduceworks`  
 Include if you're compiling with quad precision and your MPI library correctly handles 
 quad precision arithmetics. This directive is ignored if you're not using MPI or double 
 precision is used.

* `-Dsplitreduce`  
 Include if MPI_REDUCE from your library does not work correctly for data structures
 larger than ~2GB (i.e. limit of 32bit integer address). This directive is ignored if
 you're not using MPI.

* `-Dmpithree`  
 Include if you're using MPI 3.0 standard. This feature only affects the module 
 mpi_memory_mod.F90 which is used by the UKRmol+ program MPI-SCATCI. Note that if GBTOlib
 is linked with MPI-SCATCI the integer precision of GBTOlib must be -i8 since that is explicitly
 required by the MPI-SCATCI module CSF_module.F90 which calls functions from mpi_memory_mod.F90.

When starting up, GBTOlib normally reads all molecular basis and integral data from disk into memory.
However, for applications that require only a subset of the integrals, or in environments with little
memory available, this can lead to unnecessary memory consumption. To overcome this, GBTOlib can be
compiled in a "memory mapping" mode, where large arrays are directly mapped from disk to (virtual)
memory:

* `-Dusemapping`  
  Replace reading of some arrays from disk by direct mapping of the data to memory. Implementation
  of mapping is provided for POSIX-compatible systems and for MS Windows. When the auxiliary
  environment variable `GBL_DEBUG_MMAP` is defined during execution of a program, any attempt to
  map memory by GBTOlib will result in additional diagnostic output.

The compilation itself, as well as operation of the code have been tested with the following
toolsets:

* Doxygen 1.8.14, 1.8.15 (note that version 1.6.1 does NOT work)

* CMake 3.2.3, 3.6.2, 3.12.0, 3.13.0

* GCC 9.0, OpenBLAS 0.3.1 ILP64, Open MPI 3.1.1 ILP64 

* Intel Parallel Studio XE 2017.1, MKL 2017.1, Intel MPI 5.

* Intel Parallel Studio XE 2017.3.191, MKL 2017.2, Intel MPI 5. 

* Intel Parallel Studio XE 2018.3, MKL 2018.3, Open MPI 3.1.1 ILP64.

On the contrary, the following software is known to be incompatible:

* Cray Fortran 8.5 (bug in repositioning stream files, fixed in 8.7.1)

Test suite
----------

* When the build has been completed tests can be run to determine if the suite is running correctly 
  in both the serial and the parallel mode.

* To run the test suite type

    make test

* However, these tests only execute the integral calculation but don't check correctness
  of the integrals! For a complete check of the calculation the test suite for UKRmol-in
  must be run.


Generating input for GBTOlib
----------------------------

* The atomic integral calculation can be executed when the molecular geometry and the target
  atom-centered Gaussian basis set has been specified. Evaluation of the molecular integrals requires
  also the coefficients for the molecular orbitals. This information can be passed to the
  library either manually following the steps implemented in the program `scatci_integrals` or
  with the help of the object `molden_input_obj` which allows to read this information directly from a
  file in the [MOLDEN](http://cheminf.cmbi.ru.nl/molden/molden_format.html) format. This is also
  the strategy implemented in `scatci_integrals`.

* The MOLDEN file can be generated by a range of Quantum Chemistry software. GBTOlib library has been
  written to work mostly with MOLDEN files generated by [MOLPRO](http://www.molpro.net/) and the 
  open-source software [PSI4](http://www.psicode.org/).

* Note that the format of the MOLDEN file has a rather loose specification so typically you'll find
  differences in the format of the file produced by different software. Therefore don't be surprised
  if the MOLDEN file produced by software other than MOLPRO and PSI4 does not automatically work. In 
  most cases this can be solved by simple manual tweaking of the produced files. If this is not acceptable
  for you then either contact the developers or contribute to the code by extending the capabilities of
  the `molden_input_obj`.

Running `scatci_integrals`
--------------------------

* The program `scatci_integrals` generates 1- and 2-electron atomic and molecular integrals for the basis
  specified using the externally generated Molden file and the custom continuum basis specified directly
  in the input file.

* See [the manual](doc/scatci_integrals.md) for a list of known input file options.

* Examples of the input files can be found in the 'test' directory. By default the program takes the input
  from the file 'inp' which must be placed in the same directory where `scatci_integrals` is launched. 
  Therefore to run manually one of the inputs go to the main directory for the test and copy the input file
  from the folder 'inputs', e.g.:
  ```
    cd tests/D2h_minimal_scattering_integrals
    cp ./inputs/target.integrals.inp ./inp
    mpirun -n 2 ./scatci_integrals
  ```

  This will launch `scatci_integrals` using two MPI tasks.

* Alternatively the input file can be supplied as a command-line argument thus removing the need to copy
  the input file, e.g.:
  ```
    cd tests/D2h_minimal_scattering_integrals
    mpirun -n 2 ./scatci_integrals ./inputs/target.integrals.inp
  ```

Interfacing GBTOlib with other codes
------------------------------------

* You can use the program `scatci_integrals` to generate the `moints` file containing the
  atomic and molecular bases and your selected integrals.

* The program `get_integrals` from `source/programs/get_integrals_template.f90` contains a
  simple example how the contents of the `moints` file can be loaded into your program and the
  individual 1- and 2-electron integrals retrieved.

Notes for developers
--------------------

A commmit.template is available to ensure that developers do not forget to provide sufficient
information on thei work they've done. We ask that you do:

    git config commit.template commit.template

This way, when you commit some changes, the commit log will be pre-filled with the information in the
template.


 ### Git-flow ###

From version 1.0 the GBTOlib repository has switched to using an adapted version of
[Git-flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).
Git-flow defines a set of branches and prefixes which have a well-defined meaning.

In the case of GBTOlib the set of permanent branches and the naming convention is the following
(the flow chart [here](https://gitlab.com/Uk-amor/UKRMol/GBTOLib/wikis/home) helps understanding
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


Working on your own developments
--------------------------------

To create a branch for your own developments use the following commands:

    git checkout master
    git branch name_of_feature

You'll then checkout the branch, work on it and commit it when needed. Once you're satisfied
with your implementation, you should create a merge request by using the facility in the GitLab
interface. One of the managers/maintaners of the project will check it and merge it to `master`.

(When the development of one or several new features has finished and the code in the `master`
branch been thoroughly tested, the codes will be merged into the  `release` branch  and tagged.)



### Other ###

* The library code adheres to using 132-character-long lines. If you're using VIM
  to edit the source code include the following line into your `.vimrc` file:
  `set wrap; set textwidth=132`.

* If you use IFORT as your default compiler then please do make sure that you compile and test
  your code also with GFORTRAN. The IFORT compiler contains a number of non-standard extensions
  which are not supported by GFORTRAN. Using these extensions therefore kills code portability.

* OpenMP sections sometimes use `DEFAULT(SHARED)`. This approach is necessary in case the
  list of `SHARED()` variables includes objects whose type-bound variables are accessed
  in the parallel section. GFORTRAN by default does not automatically make shared some
  internal variables needed to correctly access the type-bound variables of shared objects.
  Apparently, IFORT does not suffer from this problem. Due to this problem please pay extra
  attention to testing new OpenMP sections by compiling and running the library with GFORTRAN.

* The library code is predominantly F2003 standard. The only F2008 feature used is the 
  `newunit` parameter in the `open` statement. Avoid using F2008 standard as much as possible
  unless its use brings significant gains in performance or code clarity. This rule is in place 
  due to a still unsatisfactory reliability of the F2008 standard in the GFORTRAN and IFORT compilers.

* The binary representation of LOGICALs is heavily compiler dependent (for instance gfortran uses 0/1,
  while ifort 0/-1 for false/true). In contrast, INTEGER types are de facto standard (at least on
  the most common x86 architecture). For this reason, GBTOlib uses INTEGERs 0/1 rather than LOGICALs
  false/true when writing the "moints" file to achieve maximal portability across compilers.

* The GBTOlib container can be built on GitLab by committing a message containing the keyword "DEPLOY".

### Known issues ###

* When compiled with Cray Fortran compiler 8.7.7 and quad precision (i.e. using `-Dusequadprec`) 
  `scatci_integrals` crashes when reading the Molden file. This is a compiler issue which will be 
  hopefully resolved in newer versions. For details see: 
  [GitLab](https://gitlab.com/Uk-amor/UKRMol/GBTOLib/issues/35).
