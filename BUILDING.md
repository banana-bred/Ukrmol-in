Building UKRmol+
================

UKRmol+ uses the free and open source tool [CMake](https://cmake.org/download/) for configuration, compilation, testing and
installation. The newest CMake can always be downloaded as an unpack-and-go binary package. This document assumes CMake 4.2.0,
but a wide range of versions works equivalently.

The build process generally consists of the following stages:

 1. Configuration
 2. Compilation (command "`cmake --build .`", or often just "`make`")
 3. Installation (optional; command "`cmake --install .`", or often just "`make install`")
 4. Testing (optional; command "`ctest`", on often just "`make test`")

The first step is the only complex one and it depends on the compilers used and on the required functionality of the UKRmol+
program suite. In the examples below it is assumed that the configuration and building of UKRmol-in is done in a dedicated
subdirectory (named e.g. "build") in the root folder of the unpacked UKRmol-in source package. The first run of CMake is then
"`cmake [options discussed below] ..`", where the two dots refer to the parent directory containing UKRmol-in files. Equivalently,
a full path to that directory can be used, too. The most important CMake options that can be specified in the command are explained
below.

During the configuration step, CMake combines the UKRmol+ project setup with the user-provided details and produces instructions for
a lower-level build tool. In Unix-like systems it is by default the program "`make`". However, this automatic choice can be
overriden to a different generator using the flag `-G`. For example, "`cmake -G Ninja`" will switch to the `ninja` build
manager (needs to be installed separately); the command "`cmake --build .`" then becomes equivalent to simply typing "`ninja`".

There are many CMake configuration flags that affect how UKRmol+ will be built. A frequently used subset is listed below.
Some of them can be also specified by exporting a corresponding environment variable. However, a parameter passed directly
to CMake has always precedence.

 - `CMAKE_C_COMPILER` (environment variable `CC`): Which C compiler executable to use.
 - `CMAKE_Fortran_COMPILER` (environment variable `FC`): Which Fortran compiler executable to use.
 - `CMAKE_C_FLAGS` (environment variable `CFLAGS`): Command line options to pass to the C compiler.
 - `CMAKE_Fortran_FLAGS` (environment variable `FFLAGS`): Command line options to pass to the Fortran compiler.
 - `BLA_VENDOR` (environment variable `BLA_VENDOR`): Which BLAS/LAPACK implementation to use. See "BLAS/LAPACK details" below.
 - `GBTOlib_Fortran_FLAGS`: MPI-related flags to use when compiling the integral library GBTOlib.
 - `BUILD_DOC`: Build the automatically generated (Doxygen) documentation.
 - `BLAS_LIBRARIES`: List of BLAS libraries to use if `BLA_VENDOR` is not used.
 - `LAPACK_LIBRARIES`: List of LAPACK libraries to use if `BLA_VENDOR` is not used.
 - `SCALAPACK_LIBRARIES`: List of paths to ScaLAPACK libraries to link to.
 - `ARPACK_LIBRARIES`: Path to the Arpack library to link to.
 - `SLEPC_INCLUDE_DIRS`: List of paths containing PETSc and SLEPc Fortran include directories.
 - `SLEPC_LIBRARIES`: List of PETSc and SLEPc libraries to link to.
 - `UKRMOL_OUT_DIR`: Path to the UKRmol-out directory if joint compilation (and testing) of both is desired.
 - `CLBLAST_INCLUDE_DIRS`: Path to CLBlast C include directory (only used by UKRmol-out).
 - `CLBLAS_LIBRARIES`: Path to the CLBlast library to link to (only used by UKRmol-out).
 - `MPIEXEC_EXECUTABLE`: Path to the MPI launcher (`mpiexec`, `aprun`, etc).
 - `MPIEXEC_PREFLAGS`: Extra flags for the MPI launcher (e.g., processor affinity setup).
 - `WITH_MPI`: Set to `OFF` to explicitly avoid enabling MPI support even if MPI is available.

Note that if you ever link UKRmol+ to shared libraries, these have to be visible for the operating system both during compilation
and during later execution: Either they need to be installed in a standard system location (/usr/lib, /usr/local/lib, etc.) or
their directory has to be published by appending its path to the environment variable `LD_LIBRARY_PATH` in the current shell
session. Some packages (e.g. Intel oneAPI) are shipped with environment scripts that do this automatically upon sourcing them
in the current terminal shell. When compiling an external library it is often possible to control whether a shared (.so, .dll)
or static library (.a) will be produced by the configuration options `-D BUILD_SHARED_LIBS=ON` (in CMake based projects)
and `--enable-shared --disable-static` (in GNU Autotools based projects) and similar.

The compiler preprocessor flags affecting UKRmol+ through configuration of GBTOlib (`-Dusequadprec`, `-Dquadpreckind=X`, `-Dusempi`,
`-Dquadreduceworks`, `-Dsplitreduce`, and `-Dmpithree`) are explained in the main [README.md](README.md).


### Sample basic configurations

Below we provide configuration step examples for selected compiler suites. If none of them suits your needs, details on the
individual components are discussed in the next sections. The examples assume a Unix-like operating system, but analogical
instructions are valid for Microsoft Windows and other platforms supported by CMake.

##### Intel oneAPI with ILP64 BLAS, LAPACK, MPI, and ScaLAPACK

The easiest configuration is using the Intel-specific option `-qmkl-ilp64=cluster`. This automatically links the code also to
Intel implementation of ScaLAPACK, which is used by MPI-SCATCI. The configuration can be specified using environment variables

    export CC=$(which icc)
    export FC=$(which mpiifx)
    export FFLAGS="-i8 -qmkl-ilp64=cluster"
    export BLA_VENDOR=Intel10_64ilp_seq

    cmake ..

or using CMake command line arguments

    cmake -D CMAKE_C_COMPILER=$(which icc) \
          -D CMAKE_Fortran_COMPILER=$(which mpiifx) \
          -D CMAKE_Fortran_FLAGS="-i8 -qmkl-ilp64=cluster" \
          -D BLA_VENDOR=Intel10_64ilp_seq \
          ..

However, `-qmkl-ilp64=cluster` results in linking to non-threaded BLAS/LAPACK libraries. (The value `Intel10_64ilp_seq` is set
consistently.) Single-threaded BLAS/LAPACK might be suboptimal for some use cases of UKRmol+. If both threaded BLAS/LAPACK
as well as linking to ScaLAPACK is required, the ScaLAPACK libraries need to be specified by hand. (`BLA_VENDOR` can be then
set consistently to `Intel10_64ilp`.)

    cmake -D CMAKE_C_COMPILER=$(which icc) \
          -D CMAKE_Fortran_COMPILER=$(which mpiifx) \
          -D CMAKE_Fortran_FLAGS="-i8" \
          -D BLA_VENDOR=Intel10_64ilp \
          -D SCALAPACK_LIBRARIES="-L$MKLROOT/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64" \
          ..

The list of libraries in `SCALAPACK_LIBRARIES` can be specified using compiler flags like above or as a list of absolute paths:

    SCALAPACK_LIBRARIES="$MKLROOT/lib/intel64/libmkl_scalapack_ilp64.so;$MKLROOT/lib/intel64/libmkl_blacs_intelmpi_ilp64.so"

If, for some reason, the automatic discovery of BLAS and LAPACK does not work, it is possible to specify the paths on the
command line manually (again via flags or paths):

    cmake -D CMAKE_C_COMPILER=$(which icc) \
          -D CMAKE_Fortran_COMPILER=$(which mpiifx) \
          -D CMAKE_Fortran_FLAGS="-i8" \
          -D BLAS_LIBRARIES="-L$MKLROOT/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core" \
          -D LAPACK_LIBRARIES="-L$MKLROOT/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core" \
          -D SCALAPACK_LIBRARIES="-L$MKLROOT/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64" \
          ..

##### System GNU compilers without MPI

GNU compilers are often preinstalled in common GNU/Linux distributions, or available in system repositories. They do not provide
parallel MPI compilers by themselves, neither BLAS or LAPACK implementation. However, the latter is often available in the system
independently, in the form of reference BLAS and LAPACK (libblas, liblapack), OpenBLAS (libopenblas) or another. CMake will
automatically discover the best usable installed option. Only the switching to 8-byte integers is then needed to be specified
manually:

    cmake -D CMAKE_Fortran_FLAGS="-fdefault-integer-8" ..

Forcing a specific BLAS/LAPACK implementation other than the automatic pick by CMake can still be done using the `BLA_VENDOR`
variable or the `BLAS_LIBRARIES` and `LAPACK_LIBRARIES` options.

If a system installation of MPI is avaiable in some standard location and its environment activated, the above command will
enable MPI support in UKRmol+ automatically.


### Compiler details

UKRmol+ has been successfully used with Intel and GNU compilers. Other compilers (most notably Cray) may work as well; however,
compatibility with other compilers is not actively maintained and minor tweaks in the code might be required in some places
for a successful compilation. Whenever there are multiple compilers installed in a system, it might be necessary to specify the one
to use with CMake using the environment variables `CC` and `FC` or the CMake command line arguments `CMAKE_C_COMPILER` and
`CMAKE_Fortran_COMPILER`.

##### Intel oneAPI

As of this UKRmol+ release, Intel oneAPI 2025 is the current iteration of Intel compilers distribution. Its serial C and Fortran
compilers are called `icx` and `ifx`. The parallel compilers are `mpiicx` and `mpiifx`. Note that Intel compilers are generally
aggressively utilising stack memory to achieve maximal speed, so it may be necessary to force unlimited stack when using a UKRmol+
program compiled with Intel. This can be done by issuing the following command prior to executing a UKRmol+ program:

    ulimit -s unlimited

Outputs printed to the terminal by programs compiled with Intel compilers are wrapped after 80 characters, which might interfere
with the original intended output layout. This breaking of printed lines can be suppressed by adding `-no-wrap-margin` among the
Fortran compiler flags. The fastest programs can be achieved by adding `-xHOST`, which allows the compiler to use native
hardware-specific processor instructions. However, such binaries may not work correctly on other computers.

Use the `-i8` Fortran compiler option to promote default Fortan integers in UKRmol+ to long (8-byte) integers.

##### Intel classic

The "classic" compilers `icc` and `ifort` are only present in older distributions of the Intel development platform.
The corresponding MPI compilers are `mpiicc` and `mpiifort`.

Use the `-i8` Fortran compiler option to promote default Fortan integers in UKRmol+ to long (8-byte) integers.

##### GNU Compiler Collection (GCC)

The free GNU compiler suite can be used to build UKRmol+ since GCC 8.2, though more performance on modern CPUs can be achieved
with newer versions. The C compiler is called `gcc`, the Fortran one `gfortran`. When compiling the code on a computer for use
on the same computer, it is often useful to include `-march=native` among the Fortran compiler flags (`FFLAGS` or
`CMAKE_Fortran_FLAGS`). This will allow the compiler to use hardware-specific instructions instead of a slower generic machine code.
Additionally, it is often useful to include the flags `-fno-signed-zeros -fno-trapping-math -fassociative-math`, which relax
strict requirements on reproduciblity of floating point operations and which enable automatic use of vector instructions in
certain places. Finally, `gfortran` allows delegating matrix multiplications to a dedicated BLAS library by use of the compiler
option `-fexternal-blas`, which speeds up generic Fortran matrix multiplication. Note however that this is only possible when
using the default 4-byte integer BLAS interface.

Use the `-fdefault-integer-8` Fortran compiler option to promote default Fortan integers in UKRmol+ to long (8-byte) integers.

##### Cray

Older versions of Cray compilers (< 8.7.1) had a bug that rendered stream I/O with Fortran unusable. All tested Cray compilers seem
to struggle with sourced allocation of arrays of quadruple precision floating numbers.

Use the `-default64` Fortran compiler option to promote default Fortan integers in UKRmol+ to long (8-byte) integers.


### BLAS/LAPACK details

BLAS and LAPACK are the basic linear algebra libraries used thoughout the code. Intel ships its implementation as part
of the Math Kernel Library (MKL), which is a component of Intel oneAPI. Alternatives to Intel MKL include OpenBLAS, reference
BLAS/LAPACK and others. UKRmol+ supports 4-byte and 8-byte integer interfaces (LP64 and ILP64 model), which are detected
automatically by CMake.

##### Intel MKL BLAS/LAPACK

Intel MKL contains BLAS and LAPACK in both 4-byte and 8-byte integer interface versions. MKL can be used both with Intel compilers
as well as with GNU compilers. The library consists of several independent layers and choosing the correct combination depends
on the compiler used, the integer interface and threading model required. See the official
[link line advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html)
for all possibilities. In CMake, the required libraries can be listed explicitly in `BLAS_LIBRARIES` and `LAPACK_LIBRARIES`, or
the automatic discovery can be utilised. This is done by setting the CMake option or environment variable `BLA_VENDOR` to one
of the following:

| `BLA_VENDOR`        | threading       | integer interface |
| ---                 | ---             | ---               |
| `Intel10_64lp_seq`  | single-threaded | 4-byte            |
| `Intel10_64lp`      | multi-threaded  | 4-byte            |
| `Intel10_64ilp_seq` | single-threaded | 8-byte            |
| `Intel10_64ilp`     | multi-threaded  | 8-byte            |

Note that if the compiler option `-qmkl` or `-qmkl-ilp64` is used (see "ScaLAPACK details" below), the CMake variable `BLA_VENDOR`
should be set consistently.

##### OpenBLAS

OpenBLAS can be downloaded [from Github]() and compiled locally. The default configuration results in single-threaded, LP64 library.
To compile a multi-threaded UKRmol+-compatible version with ILP64 interface, additional options need to be given either on
command line (as below) or by overwriting the defaults in the file "Makefile.rule".

    cd OpenBLAS-0.3.30
    make PREFIX=$HOME/Software/openblas-0.3.30-ilp64 USE_THREAD=1 USE_OPENMP=1 INTERFACE64=1
    make PREFIX=$HOME/Software/openblas-0.3.30-ilp64 USE_THREAD=1 USE_OPENMP=1 INTERFACE64=1 install

Use of OpenBLAS can be forced by setting `BLA_VENDOR=OpenBLAS`. However, if it is available in some standard system location,
this is typically not needed, as CMake finds it automatically. Alternatively, path to the OpenBLAS library can be specified using
`BLAS_LIBRARIES` and `LAPACK_LIBRARIES`.


### MPI details

The MPI library allows parallel processing of the calculation on multiple cores and nodes. It can be utilised by SCATCI_INTEGRALS,
MPI-SCATCI and some UKRmol-out programs. When using Intel compiler suite, it is natural to use the bundled Intel MPI, which
comes in both LP64 and ILP64 variants. Other HPC environments may provide different MPI implementations; UKRmol+ was tested
also with Open MPI 5.0.9 and Microsoft MPI 10.1.3.

CMake is typically able to find MPI by itself, automatically, provided it is installed in a standard system location or its
environment is activated. In case of doubts, compiler wrappers (`mpiicx`, `mpicc`, `mpiifx`, `mpiifort`, `mpifort`, etc.) shipped
with the specific MPI implementation can be used instead of the plain compilers (via `CMAKE_C_COMPILER` and
`CMAKE_Fortran_COMPILER`). If the autodetection still fails to pick the required installation for some reason (e.g. if there are
more installation of MPI implementations), specifying the full path to `mpiexec` or `mpirun` using the CMake parameter
`MPIEXEC_EXECUTABLE` can be used for disambiguation:

    cmake -D MPIEXEC_EXECUTABLE=$(which mpiexec) \
          ... and so on ...

##### Intel MPI

When using the ILP64 variant of Intel MPI, writing large datasets using MPI-IO in the MPI-SCATCI program will fail, see
[the bug report](https://community.intel.com/t5/Intel-MPI-Library/Wrong-limit-on-disp-parameter-in-ILP64-version-of-MPI-File-set/m-p/1507790).
This is affecting only writing of the RMT data files. If use of this functionality is envisioned, LP64 MPI interface should be used
instead. Intel automatically selects the ILP64 version of MPI whenever the `-i8` flag is used together with the Intel MPI Fortran
compiler wrapper `mpiifx` or `mpiifort`. To keep the UKRmol+ 8-byte integers but opt out of the ILP64 MPI interface, the
combination `-i8 -no_ilp64` has to be used:

    cmake -D CMAKE_C_COMPILER=$(which icx) \
          -D CMAKE_Fortran_COMPILER=$(which mpiifx) \
          -D CMAKE_Fortran_FLAGS="-i8 -no_ilp64 -qmkl-ilp64=cluster" \
          ... and so on ...

Alternatively, one can avoid using `mpiifx` (or `mpiifort`); then CMake links to the LP64 version of Intel MPI. However, this also
requires dropping `-qmkl-ilp64=cluster` and hence specifying ScaLAPACK manually via `SCALAPACK_LIBRARIES`:

    cmake -D CMAKE_C_COMPILER=$(which icx) \
          -D CMAKE_Fortran_COMPILER=$(which ifx) \
          -D CMAKE_Fortran_FLAGS="-i8" \
          -D SCALAPACK_LIBRARIES="-L$MKLROOT/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64" \
          ... and so on ...

There is yet another complication with Intel MPI in the ILP64 version: When UKRmol+ is configured in the default way, it will
be linked both to the shared library "libmpifort" that implements the standard LP64 MPI interface and to "libmpi_ilp64"
which overrides the interface to the ILP64 version. On some platforms, when UKRmol+ programs are being executed, the order
of loading shared libraries can be unpredictable. To avoid "libmpifort" accidentally being loaded and used first, one should
use "`mpiexec -ilp64`" to execute UKRmol+ programs rather than "`mpiexec`" only. In practice this issue does not appear on GNU
systems. It can also be fully avoided by linking to static MPI libraries, where the library order is solved once and for all during
compilation. (The compiler flag `-static-intel` can be used for this purpose.) However, it is simply best to avoid the ILP64
layer altogether as it is only needed for applications that do not distinguish MPI and internal integers, which UKRmol+ does.

The integer size in MPI interfaces is in principle independent of BLAS/LAPACK integers as well as of UKRmol+ integers;
there is no requirement that any of the three are the same.

##### Open MPI

The current version of Open MPI (5.0.9) can be compiled in both LP64 and ILP64 mode. However, these are fully functionally
identical; there is no point in using the non-default ILP64 Fotran layer, as the calls are always forwarded to the same
standard LP64 MPI C routines.

The default configuration of Open MPI does not support arithmetic on quadruple precision floating point numbers. If required,
this can be enabled during custom compilation of Open MPI from source by the compiler flag `-mlong-double-128`. For example:

    cd openmpi-5.0.9
    mkdir build
    cd build
    ../configure \
        --prefix=$HOME/Software/openmpi-5.0.9 \
        --enable-mpi1-compatibility \
        --enable-mpi-fortran=usempi \
        CFLAGS="-mlong-double-128"
    make
    make install

The above-suggested option `--enable-mpi1-compatibility` is not needed by UKRmol+, but allows Open MPI to be used together
with Intel MKL ScaLAPACK (see below).

##### Microsoft MPI

The current version of MSMPI (10.1.3) does not support quadruple precision. Additionally, when used together with the current
GNU Compiler Collection (15.2), in-place parallel communication in SCATCI_INTEGRALS does not work correctly due to incorrect import
of `MPI_IN_PLACE` and other MPI parameters from "msmpi.dll". Consequently, SCATCI_INTEGRALS can be run in the serial mode only.
Additionally, compiler parameters `-DMPI_TYPECLASS_REAL=1` and `-DMPI_TYPECLASS_INTEGER=2` need to be added among the Fortran flags
to provide definitions missing in the MSMPI Fortran module.

Microsoft MPI does not provide all MPI-3 features and identifies only as MPI-2-compatible. This results in particular in suboptimal
memory consumption of MPI-SCATCI.


### ScaLAPACK details

ScaLAPACK allows operations on distributed matrices. If enabled, it will be used in MPI-SCATCI to perform diagonalization
in parallel runs. ScaLAPACK is included in Intel MKL. Alternatively, a reference implementation of ScaLAPACK can be downloaded
[from Github](https://github.com/Reference-ScaLAPACK/scalapack/releases) and compiled locally. ScaLAPACK is expected to use the
same integer interface (4- or 8-byte) as the BLAS and LAPACK libraries.

##### Intel MKL

ScaLAPACK in Intel MKL consists of two libraries: ScaLAPACK proper (computation part), and BLACS (parallel communication part).
The communication library BLACS can be backed either by Intel MPI (by linking to `libmkl_blacs_intelmpi_[i]lp64`), or some other
supported MPI library (e.g. the MPI1-enabled Open MPI, see above, by linking to `libmkl_blacs_openmpi_[i]lp64`). A standard way
of inclusion MKL ScaLAPACK is by setting the corresponding CMake option:

    -D SCALAPACK_LIBRARIES="-L$MKLROOT/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64"

Alternatively, UKRmol+ also supports the Intel-specific way of linking to ScaLAPACK by means of the `-qmkl` or `-qmkl-ilp64`
compiler option (in `FFLAGS` or `CMAKE_Fortran_FLAGS`). When using this variable, explicit specification of `SCALAPACK_LIBRARIES`
is not required. However, care should be taken to set `BLA_VENDOR` consistently (see BLAS/LAPACK above). The table below lists
the corresponding pairs.

| compiler flags                | corresponding `BLA_VENDOR`    |
| ---                           | ---                           |
| `-i8 -qmkl-ilp64`             | `Intel10_64ilp`               |
| `-i8 -qmkl-ilp64=parallel`    | `Intel10_64ilp`               |
| `-i8 -qmkl-ilp64=sequential`  | `Intel10_64ilp_seq`           |
| `-i8 -qmkl-ilp64=cluster`     | `Intel10_64ilp_seq`           |

Note that when using `-i8`, then `-qmkl` has the same meaning as `-qmkl-ilp64`. Without `-i8`, the flag
`-qmkl` links to LP64 libraries instead (`Intel10_64lp` etc).

##### Reference ScaLAPACK

An example of configuration and testing of the Reference ScaLAPACK 2.2.2 with OpenBLAS 0.3.30, Open MPI 5.0.9
and GNU Compiler Collection 15.2 in LP64 mode:

    cd scalapack-2.2.2
    mkdir build
    cd build
    cmake -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_C_COMPILER=$(which mpicc) \
          -D CMAKE_C_FLAGS="-std=gnu89" \
          -D CMAKE_CXX_COMPILER=$(which mpicxx) \
          -D CMAKE_Fortran_COMPILER=$(which mpifort) \
          -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch" \
          -D BLAS_LIBRARIES="$HOME/Software/openblas-0.3.30/lib/libopenblas.so" \
          ..
    cmake --build .
    cmake --install .
    export OPENBLAS_NUM_THREADS=1       # <-- avoid multithreading in MPI tests
    ctest

If 8-byte integer variant is desired, this can be achieved by additional compiler options and linking to appropriate BLAS variant:

    cmake -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_C_COMPILER=$(which mpicc) \
          -D CMAKE_C_FLAGS="-std=gnu89 -DInt=long" \
          -D CMAKE_CXX_COMPILER=$(which mpicxx) \
          -D CMAKE_Fortran_COMPILER=$(which mpifort) \
          -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fdefault-integer-8" \
          -D BLAS_LIBRARIES="$HOME/Software/openblas-0.3.30-ilp64/lib/libopenblas_ilp64.so" \
          ..


### Arpack details

[Arpack](https://github.com/opencollab/arpack-ng/tags) can be used by SCATCI or MPI-SCATCI to obtain a few eigenstates
of the Hamiltonian with the lowest energy. In UKRmol+
it is only available in serial calculations; it cannot be used in parallel calculations in MPI-SCATCI. It it expected to use the
same integers (4- or 8-byte) as the BLAS and LAPACK libraries. The default are 4-byte integers; switching to 8 bytes requires
adding the switch `-DINTERFACE64=ON` when configuring Arpack. For example:

    cd arpack-ng-3.9.1
    mkdir build
    cd build
    cmake -D CMAKE_INSTALL_PREFIX=$HOME/Software/arpack-ng-3.9.1-ilp64 \
          -D BLA_VENDOR=Intel10_64ilp \
          -D INTERFACE64=ON \
          ..
    cmake --build .
    cmake --install .

The compiled library is then included in UKRmol+ by providing the path during configuration of UKRmol+:

    cmake -D ARPACK_LIBRARIES=$HOME/Software/arpack-ng-3.9.1-ilp64/lib/libarpack.a \
          ... and so on ...


### PETSc and SLEPc details

SLEPc and its prerequisite PETSc can be downloaded from http://slepc.upv.es/download/ and
https://petsc.org/release/install/download/ .

SLEPc is a distributed eigenvalue solver that can be used to calculate a subset of Hamiltonian eigenstates in parallel runs
of MPI-SCATCI. It inherits integer bitness, as well as all other configuration, from PETSc. PETSc and SLEPc integers are
independent of all other libraries and, in the current version of PETSc (3.24.2), they can be switched from the default 4 bytes
to the larger 8 bytes by means of the configuration option `--with-64-bit-indices` if desired. Additionally,
PETSc links to BLAS and LAPACK and it has to use compatible integers to call routines in these two libraries. If your BLAS and
LAPACK use 8-byte integers, this can be communicated to PETSc during its configuration step by means of the option
`--with-64-bit-blas-indices=1`. The BLAS and LAPACK library has to be specified using `--with-blaslapack-lib`. For example:

    cd petsc-3.24.2
    ./configure \
        --prefix=$HOME/Software/petsc-3.24.2-ilp64 \
        --with-64-bit-blas-indices=1 \
        --with-64-bit-indices=1 \
        --with-blaslapack-lib="-L$MKLROOT/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5"

When configuring UKRmol+, `SLEPC_INCLUDE_DIRS` need to include paths to the Fortran modules in both PETSc and SLEPc, as well 
as to directory containing the "finclude" subdirectory, all separated by semicolons. Likewise, `SLEPC_LIBRARIES` need full paths
to both the PETSc and SLEPc library, again separated by semicolons. For example

    SLEPC_INCLUDE_DIRS="$HOME/petsc/include/petsc;$HOME/petsc/include;$HOME/slepc/include/slepc;$HOME/slepc/include"
    SLEPC_LIBRARIES="$HOME/slepc/lib/libslepc.so;$HOME/petsc/lib/libpetsc.so"

When using pre-4.2.0 CMake to configure UKRmol+, the subsequent compilation with a recent PETSc (3.24.2) with the `make` driver
results in a stall. Using the Ninja generator (`cmake -G Ninja`) instead or upgrading CMake to 4.2.0+ solves this issue.


### ELPA details

[ELPA](https://elpa.mpcdf.mpg.de/software/tarball-archive/ELPA_TARBALL_ARCHIVE.html) is yet another distributed eigensolver that
can efficiently calculate all eigenstates of a Hamiltonian. If enabled, it is used
in parallel runs of MPI-SCATCI. While ELPA does advertise support for ILP64 mode, it's usability hasn't been confirmed yet by
developers of UKRmol+ and is not recommended. In the default LP64 mode, though, ELPA works very well and generally beats ScaLAPACK
in speed. Configuration, compilation, installation and testing of ELPA might look like this:

    cd elpa-2025.01.002
    mkdir build
    cd build
    ./configure \
        --prefix=$HOME/Software/elpa-2025.01.002 \
        --disable-avx512-kernels \
        CC=$(which mpiicx) \
        FC=$(which mpiifx) \
        CFLAGS="-O3 -march=native" \
        FFLAGS="-O3 -march=native" \
        LDFLAGS="$HOME/scalapack-2.2.2/lib/libscalapack.so $HOME/openblas-0.3.30/lib/libopenblas.so"
    make
    make install
    make check

If ELPA is being compiled on an HPC system with AVX-512 instruction support, the flag `--disable-avx512-kernels` can be typically
omitted.
