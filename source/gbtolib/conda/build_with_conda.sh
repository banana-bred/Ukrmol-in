#!/usr/bin/env bash
# ===============================================================
#  Build script for GBTOlib using GNU + OpenMPI toolchain (via conda)
# ===============================================================
#  This script:
#   1. Creates or updates the conda environment defined in environment.yml
#   2. Activates it
#   3. Configures and builds GBTOlib twice:
#        - in double precision  (default)
#        - in quadruple precision (via -Dusequadprec)
#   4. Installs both builds into:
#        ./install/double
#        ./install/quad
# ===============================================================


# Set the script to exit immediately if any command fails ('-e')
# and to treat the exit status of a pipeline as the status of the last command
# that failed, not the very last command ('-o pipefail'). This is a "strict mode".
set -eo pipefail

# ---------------------------------------------------------------
# Step 0. Detect Mamba or Conda
# ---------------------------------------------------------------
# Check if the 'mamba' command is available in the system's PATH.
if command -v mamba &> /dev/null; then
    # If mamba is found, assign it as the package manager. Mamba is a much faster alternative to conda.
    MAMBA_OR_CONDA="mamba"
    echo "==> Using Mamba for environment management."
else
    # If mamba is not found, fall back to using the standard 'conda' command.
    MAMBA_OR_CONDA="conda"
    echo "==> Mamba not found, falling back to Conda."
fi

# ---------------------------------------------------------------
# Step 1. Create or update the environment
# ---------------------------------------------------------------
# Initialize conda's shell functions (like 'conda activate') for the current bash session.
eval "$(conda shell.bash hook)"

# Check if the conda environment named 'gbtolib_env' already exists.
# `conda env list` lists all environments, and `grep -q` quietly searches for the name.
if ! $MAMBA_OR_CONDA env list | grep -q "gbtolib_env"; then
    # If the environment does not exist, create it from the specified YAML file.
    echo "==> Creating new environment 'gbtolib_env'..."
    # `dirname "$0"` gets the directory where the script itself is located.
    $MAMBA_OR_CONDA env create -f "$(dirname "$0")/environment.yml"
else
    # If the environment already exists, update it with the packages from the YAML file.
    # This ensures the environment has all the latest required packages.
    echo "==> Updating existing environment 'gbtolib_env'..."
    $MAMBA_OR_CONDA env update -f "$(dirname "$0")/environment.yml"
fi

# A confirmation message that the environment is set up.
echo "==> Environment 'gbtolib_env' is ready."
# Activate the conda environment. The `|| source activate` is for older conda versions.
conda activate gbtolib_env || source activate gbtolib_env
# Print a confirmation that the environment is active.
echo "==> Activated environment 'gbtolib_env'."
# Print a blank line for better formatting.
echo
# List all packages installed in the current environment for verification.
echo "==> Installed packages:"
conda list


# ---------------------------------------------------------------
# Step 2. Toolchain and compiler flags
# ---------------------------------------------------------------
# Set and export compiler environment variables so CMake and other build tools will use them.
# `which` finds the full path to the MPI compiler wrappers provided by the conda environment.
export CC=$(which mpicc)      # Use the MPI C compiler wrapper.
export CXX=$(which mpicxx)    # Use the MPI C++ compiler wrapper.
export FC=$(which mpifort)    # Use the MPI Fortran compiler wrapper.
# Set and export Fortran compiler flags that will be used for the build.
export FFLAGS="-fdefault-integer-8 -std=f2018 -O3 -fopenmp"
#  -fdefault-integer-8 : Promotes default INTEGER types to 8 bytes. Essential for large-scale calculations.
#  -std=f2018          : Enforces the modern Fortran 2018 standard, which can help avoid bugs and warnings.
#  -O3                 : Enables a high level of compiler optimization for performance.
#  -fopenmp            : Enables support for OpenMP, allowing for shared-memory parallelism (threading).

# Find the MPI execution command, trying 'mpirun' first, then 'mpiexec'.
MPIEXEC=$(command -v mpirun || command -v mpiexec || echo mpiexec)

# Define preprocessor flags specific to GBTOlib to control how it uses MPI.
# Semicolon-separated for use in CMake's `-D` flag.
GBTOMPI_FLAGS="-Dusempi;-Dsplitreduce;-Dmpithree"
#  -Dusempi        : A switch passed to the preprocessor to enable MPI parallelization code paths.
#  -Dsplitreduce   : Enables a specific algorithm for MPI reduction operations, often a workaround for large data sets.
#  -Dmpithree      : Informs the code to assume at least the MPI version 3.0 API is available.

# ---------------------------------------------------------------
# Step 3. Output directories for installation
# ---------------------------------------------------------------
# Define the base directory where the final compiled libraries and binaries will be installed.
PREFIX_BASE="$(pwd)/install"
# Create the installation directories for both double and quad precision builds.
# The curly braces `{}` expand to create both directories with one command.
mkdir -p "$PREFIX_BASE"/{double,quad}
# Define a short variable for the double precision installation path.
PD="$PREFIX_BASE/double"
# Define a short variable for the quad precision installation path.
PQ="$PREFIX_BASE/quad"


# Resolve the path to the BLAS/LAPACK shared library from the active conda environment.
# On conda-forge, OpenBLAS is commonly used and provides both BLAS and LAPACK functions.
BLAS_SO="$CONDA_PREFIX/lib/libopenblas.so"
# If libopenblas.so doesn't exist, check for a versioned file like libopenblas.so.0.
[[ -e "$BLAS_SO" ]] || BLAS_SO="$CONDA_PREFIX/lib/libopenblas.so.0"
# If OpenBLAS is not found, fall back to the generic libblas.so, which might be provided by a different package.
[[ -e "$BLAS_SO" ]] || BLAS_SO="$CONDA_PREFIX/lib/libblas.so"

# ---------------------------------------------------------------
# Step 4. Configure, build, and install (double precision)
# ---------------------------------------------------------------
echo "==> Building GBTOlib [double precision] ..."

# Configure the project using CMake. The options are explained below:
#   -S . -B b_gbto_d             : Use the current directory '.' as the source and 'b_gbto_d' as the build directory.
#   -D CMAKE_BUILD_TYPE=Release  : Compile with full optimizations and no debug symbols.
#   -D CMAKE_INSTALL_PREFIX="$PD": Set the directory where `cmake --install` will place the files.
#   -D CMAKE_Fortran_FLAGS       : Pass the global Fortran flags (ILP64, std, OpenMP, etc.) to the compiler.
#   -D GBTOlib_Fortran_FLAGS     : Pass the GBTOlib-specific preprocessor flags for MPI control.
#   -D CMAKE_PREFIX_PATH         : Tell CMake to search for libraries and headers in the conda environment.
#   -D BLAS_LIBRARIES            : Explicitly tell CMake where to find the BLAS library.
#   -D LAPACK_LIBRARIES          : Explicitly tell CMake where to find the LAPACK library (same as BLAS for OpenBLAS).
#   -D MPIEXEC_EXECUTABLE        : Tell CMake's testing system (CTest) which MPI launcher to use.
cmake -S . -B b_gbto_d \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$PD" \
  -D CMAKE_Fortran_FLAGS="$FFLAGS" \
  -D GBTOlib_Fortran_FLAGS="$GBTOMPI_FLAGS" \
  -D CMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -D BLAS_LIBRARIES="$BLAS_SO" \
  -D LAPACK_LIBRARIES="$BLAS_SO" \
  -D MPIEXEC_EXECUTABLE="$MPIEXEC"

# Build the project using the configuration in the `b_gbto_d` directory.
# `-j` will use all available processor cores to build in parallel, speeding up compilation.
cmake --build b_gbto_d -j
# Install the compiled files (libraries, executables, modules) into the prefix directory `$PD`.
cmake --install b_gbto_d

# ---------------------------------------------------------------
# Step 5. Configure, build, and install (quad precision)
# ---------------------------------------------------------------
echo "==> Building GBTOlib [quad precision] ..."

# Configure the project again for quad precision in a new build directory 'b_gbto_q'.
# The options are mostly the same as the double precision build, with two key changes:
#   -D GBTOlib_Fortran_FLAGS: We add the '-Dusequadprec' flag to enable 128-bit floating point arithmetic.
#   -D CMAKE_INSTALL_PREFIX : The installation path is changed to the quad-precision directory '$PQ'.
cmake -S . -B b_gbto_q \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$PQ" \
  -D CMAKE_Fortran_FLAGS="$FFLAGS" \
  -D GBTOlib_Fortran_FLAGS="$GBTOMPI_FLAGS;-Dusequadprec" \
  -D CMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -D BLAS_LIBRARIES="$BLAS_SO" \
  -D LAPACK_LIBRARIES="$BLAS_SO" \
  -D MPIEXEC_EXECUTABLE="$MPIEXEC"

# Build the quad precision version in parallel.
cmake --build b_gbto_q -j
# Install the quad precision version into the prefix directory `$PQ`.
cmake --install b_gbto_q


# Deactivate the conda environment, returning the shell to its previous state.
conda deactivate
echo "==> Deactivated conda environment 'gbtolib_env'."

# ---------------------------------------------------------------
# Step 6. Summary
# ---------------------------------------------------------------
echo
# Print a final success message summarizing where the builds were installed.
echo "GBTOlib successfully built and installed to:"
echo "  • $PD  (double precision)"
echo "  • $PQ  (quad precision)"