Test suite for UKRmol+ and GBTOlib
==================================

The tests contained here are integration tests which execute a range of calculations using UKRmol+
and GBTOlib.

There are several types of tests:

1. Target tests: no particle in the continuum, so these do not include the continuum basis
   - Hartree-Fock (HF)
   - Complete Active Space (CAS)

   For each molecule three tests have been constructed, each corresponding to a different orientation 
   of the molecule with respect to the three Cartesian axes.

2. Scattering tests: these do include the continuum basis
   - Static Exchange plus Polarization (SEP) based on HF target representation
   - Close-Coupling (CC) based on CAS target representation
   - R-matrix with pseudostates based on the inclusion of a pseudocontinuum basis
   - Positron scattering tests, with and without use of pseudostates (currently disabled)
   - For each molecule a test for only one orientation is included. It is assumed that the effects
     dependent on orientation affect predominantly the target part of the calculation.

3. Photoionization test: these include the B-spline continuum basis
   - Close-Coupling (CC) based on CAS target representation
   - At the moment only one test (for H2 in D2h point group) is available

4. RMT input file production: these include the B-spline continuum basis
   - Close-Coupling (CC) based on CAS target representation
   - At the moment only one test (for H2 in D2h point group) is available

The electron scattering test have been constructed for 8 point groups: C_{2v}, D_{2h}, C_s, C_{2h}, 
C_2, C_1, C_i, D_2. Each (main) directory for all tests is named after the point group the target 
molecule belongs to. 

Executing the test suite
------------------------

The test suite is integrated with CMake. When only GBTOlib and UKRmol-in testing is required,
no extra parameters need to be passed to CMake in order to enable testing. If one wishes to 
run UKRmol-out as well, the following additional parameter needs to be passed to CMake:

    -D UKRMOL_OUT_DIR="...path to UKRmol-out directory..."

In the latter case, all UKRmol-out programs will also be compiled together with UKRmol-in. Once the
compilation (`make`) is finished, one can execute the full test suite still from the build directory.
Depending on the way the test suite was configured, either serial only, or both serial and parallel 
tests should be run. If MPI was enabled, most programs will be executed via the MPI launcher (mpirun,
mpiexec, aprun, etc.), even though only the truly parallel ones (scatci_integrals, mpi-scatci, mpi_rsolve) 
will use 2 processes; the rest will run in single process only. To control what launcher will be used 
and what additional parameters will be passed to it besides the number of processes, the following CMake 
options can be taken advantage of:

    -D MPIEXEC_EXECUTABLE=$(which aprun)              # example: use "aprun" as the MPI launcher
    -D MPIEXEC_PREFLAGS="-genv;OMP_NUM_THREADS;4"     # example: restrict Intel MPI tasks to 4 cores each

(so, for example, will run the tests using: ctest -D MPIEXEC_EXECUTABLE=$(which aprun) -D MPIEXEC_PREFLAGS="-genv;OMP_NUM_THREADS;4")

The inputs for the integral part of the tests are taken from the GBTOlib repository while inputs for UKRmol+
programmes are taken from the UKRmol+ repository. GBTOlib also contains the input Molden files, which
can be used as provided, or generated as part of the test run using Molpro or Psi4. Molpro or Psi4 can be
enabled by the options

    -D WITH_MOLPRO=ON
    -D WITH_PSI4=ON

In such a case the programs need to be present in $PATH. Note that different quantum chemistry programs
produce slightly numerically different Molden files, so using the test suite in this way may lead to
slightly different results than with the provided Molden files.

The tests have been prepared with UKRmol+ and GBTOlib compiled in double precision and 64bit integers.
To run all tests for the suite execute, from the build directory, the commands:

    ctest -R serial
    ctest -R parallel

This will execute all target and scattering tests. Running the full suite of tests may require a whole day
depending on the performance of your machine. The test suite generates the following set of files:

    output.target.properties.dat    # contains matrix elements calculated by denprop
    output.scattering.eigs.dat      # contains R-matrix poles from scattering calculations
    output.scattering.eigph.dat     # contains eigenphase sums from scattering calculations
    output.scattering.ixsecs.dat    # contains cross sections from scatterig calculations
    output.molecular_data.dat       # contains subset of RMT input (but only when Python 3 is available)

These files can be compared, after running in turn  the serial and/or paralel set of tests, with the equivalent ones from 
the directory:

    test_benchmark_outputs

Some of the programmes of the suite allocate rather large arrays on the stack so we recommend that your machine
does not impose limits on the stack size using `ulimit -s unlimited`.

The reference results provided were calculated on a machine with Intel Xeon E5-2630v4 CPU, with 10 OpenMP threads
per MPI process. The serial part of the test suite took 1 hour 47 minutes. The memory required by the test
suite is partly proportional to the number of threads employed, but a workstation with 32 GiB should be enough.

Note that one can display a list of all tests in the suite using the command `ctest -N`. The naming of most tests
works like this:

    <stage>_<group>_<class>_<model><orientation>_<mode>

where

  - `stage` is one of: `integrals`, `ukrmolin`, `ukrmolout`
  - `group` is one of: `C1`, `C2`, `Cs`, `Ci`, `C2h`, `C2v`, `D2h`
  - `class` is one of: `target`, `scattering`, `photoionization`, `rmt_data`, `scattering_PCO`, 
    `scattering_Positron`,`scattering_Positron_PCO`, `scattering_ECP`, etc.
  - `model` is one of: `HF` & `CAS` (for `target`), `SEP` & `CC` (for others)
  - `orientation` is one of: `1`, `2`, `3` (only used for `target`)
  - `mode` is one of: `serial`, `parallel`

It is possible to run any subset of the tests whose names match a simple regular expression. Here are some examples:

    ctest -R 'integrals'
    ctest -R 'ukrmol'
    ctest -R 'parallel'
    ctest -R 'ukrmolout_.*_scattering_.*_parallel'
    ctest -R 'PCO'

Checking correctness of the results
-----------------------------------

The test suite can automatically compare obtained results to the benchmark ones using the freely available
program [numdiff](https://www.nongnu.org/numdiff/). If use of numdiff is desired, it can be enabled by
the CMake option

    -D WITH_NUMDIFF=ON

In such a case numdiff will be searched in $PATH. Note however, that the results are often quite dependent on the
hardware, libraries and optimization level used, so that blind machine comparison results in many false positives
(and thus failed tests). It is thus often more desirable to compare the results manually.

We recommend using a graphical diff tool, for instance [xxdiff](http://furius.ca/xxdiff/) or [meld](https://meldmerge.org/),
to compare your test files with the reference ones from the directory `test_benchmark_outputs`. Here is a brief
guidance for estimation of quality of the results:

* The target energies should be accurate to at least 9 digits.

* The magnitude of the target properties should agree to at least 3 digits, although normally will be 
  around 6. (The signs of the transition property matrix elements can differ).

* The energies of the R-matrix poles should be accurate to least 9 digits for the lowest few ones. You
  can observe decrease in precision for the higher energy ones but these should still agree to within
  approx. 6 digits.

* The eigenphase sums should agree to at least 3 digits. The agreement can be better than this in some
  energy ranges. You can ignore the difference in the time-stamp on the line just below `EIGENPHASE SUMS`

* In case of the tests run in the parallel mode you can expect a reduction in relative precision of the 
  energies of the R-matrix poles and observables in comparison to the benchmark. When in doubt about correctness
  of the outputs plot the two sets of cross sections together and check that they visually agree.
