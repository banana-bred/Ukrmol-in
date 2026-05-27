#!/usr/bin/env bash
# =====================================================================
# UKRmol-in / GBTOlib — ScaLAPACK (source) + ELPA (generic kernels only)
# - ScaLAPACK 2.2.x shared (no tests)
# - ELPA generic kernels only, OpenMP + mpif.h, installed in $CONDA_PREFIX/opt
# - Uses bash as CONFIG_SHELL to avoid 'FCFLAGS: command not found'
# - Wires ELPA includes+rpaths into UKRmol (double & quad)
# - Forces default INTEGER*8 in UKRmol+ only; deps stay LP64
#   (BLAS/LAPACK/ScaLAPACK integers are integer(blasint),
#    MPI integers are integer(mpiint); kinds auto-detected by CMake)
#
# Env toggles:
#   UKRMOL_ENV_YML   (default: ./environment.ukrmol.scalapack.yml)
#   UKRMOL_ENV_NAME  (default: ukrmol-scalapack-env)
#   SCALAPACK_VER    (default: 2.2.2)
#   ELPA_VER         (default: 2025.06.001)
#   NPROC            (default: autodetect)
# =====================================================================

# 'set -e' will cause the script to exit immediately if a command fails.
# 'set -E' ensures that the ERR trap is inherited by functions.
# 'set -o pipefail' causes a pipeline to return the exit status of the last command that exited with a non-zero status.
set -Eeo pipefail
# Set a trap for the ERR signal. If any command fails, it will print an error message with the line number and the command.
trap 'echo "ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR
# If the DEBUG environment variable is set to "1", enable command tracing ('set -x').
[[ "${DEBUG:-0}" == "1" ]] && set -x

# ---------- helpers ----------
# Function to check if a file is a valid, non-empty gzipped tarball.
valid_tgz() { [[ -s "$1" ]] && tar tzf "$1" >/dev/null 2>&1; }
# Function to download a gzipped tarball from a URL.
fetch_tgz() {
  # Assign the first argument to 'url' and the second to 'out'.
  local url="$1" out="$2"
  # Print the URL being downloaded.
  echo "    -> $url"
  # Remove the output file if it already exists to ensure a fresh download.
  rm -f "$out"
  # Check if 'curl' is available.
  if command -v curl >/dev/null 2>&1; then
    # Use curl to download the file, following redirects (-L), failing on server errors (-f), and retrying up to 4 times.
    curl -fL --retry 4 -o "$out" "$url"
  else
    # If curl is not available, use wget.
    wget -O "$out" "$url"
  fi
  # Validate the downloaded file.
  valid_tgz "$out"
}

# ---------- conda env ----------
# Check if mamba is available; if so, use it as the package manager ($PKG), otherwise use conda. Mamba is faster.
if command -v mamba >/dev/null 2>&1; then PKG=mamba; else PKG=conda; fi
# Configure the package manager to use strict channel priority. '|| true' prevents the script from exiting if the command fails.
$PKG config --set channel_priority strict >/dev/null 2>&1 || true
# Check if conda is installed; if not, print an error and exit.
command -v conda >/dev/null 2>&1 || { echo "ERROR: 'conda' not found"; exit 1; }
# Initialize the conda shell functions for the current bash session.
eval "$(conda shell.bash hook)"

# Inform the user which package manager is being used.
echo "Using $PKG for environment management"
# Get the absolute path of the directory where the script is located.
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
# Get the absolute path of the parent directory (the repository root).
REPO_ROOT="$(cd "$BASE_DIR/.." && pwd)"

# Set the path to the conda environment YAML file. Use the environment variable if set, otherwise use a default path.
UKRMOL_ENV_YML="${UKRMOL_ENV_YML:-$BASE_DIR/environment.ukrmol.scalapack.yml}"
# Set the name of the conda environment. Use the environment variable if set, otherwise use a default name.
UKRMOL_ENV_NAME="${UKRMOL_ENV_NAME:-ukrmol-scalapack-env}"

# Check if the environment YAML file exists.
if [[ -f "$UKRMOL_ENV_YML" ]]; then
  # Check if the conda environment already exists.
  if conda env list | awk '{print $1}' | grep -qx "$UKRMOL_ENV_NAME"; then
    # If it exists, update it using the YAML file.
    echo "==> Updating env '$UKRMOL_ENV_NAME' from $UKRMOL_ENV_YML …"
    # The '--prune' option removes any packages in the environment that are not in the YAML file.
    $PKG env update --prune -y -f "$UKRMOL_ENV_YML"
    echo "==> Updated env '$UKRMOL_ENV_NAME'"
  else
    # If it does not exist, create it from the YAML file.
    echo "==> Creating env '$UKRMOL_ENV_NAME' from $UKRMOL_ENV_YML …"
    $PKG env create -y -f "$UKRMOL_ENV_YML"
    echo "==> Created env '$UKRMOL_ENV_NAME'"
  fi
  # Print the list of installed packages in the environment for verification.
  echo installed packages:
  $PKG list -n "$UKRMOL_ENV_NAME"
  # Activate the conda environment.
  echo "==> Activating env '$UKRMOL_ENV_NAME' …"
  conda activate "$UKRMOL_ENV_NAME"
else
  # If the YAML file is not found, check if the environment exists anyway.
  if conda env list | awk '{print $1}' | grep -qx "$UKRMOL_ENV_NAME"; then
    # If it exists, activate it.
    conda activate "$UKRMOL_ENV_NAME"
  else
    # If neither the YAML nor the environment exists, use the currently active environment.
    echo "==> No env YAML provided; using CURRENT env: ${CONDA_PREFIX:-<none>}"
  fi
fi

# Print the path of the active conda environment.
echo "CONDA_PREFIX = ${CONDA_PREFIX:-<none>}"
# Exit if there is no active conda environment.
[[ -n "$CONDA_PREFIX" ]] || { echo "ERROR: no active conda env"; exit 2; }

# ---------- compilers & flags ----------
# Unset common compiler environment variables to ensure a clean state.
unset CC CXX FC F77 F90
# Set the C compiler to the MPI C compiler wrapper. 'export' makes it available to child processes.
export CC="$(command -v mpicc)"
# Set the C++ compiler to the MPI C++ compiler wrapper.
export CXX="$(command -v mpicxx)"
# Set the Fortran compiler to the MPI Fortran compiler wrapper.
export FC="$(command -v mpifort)"
# Check that all MPI compilers were found and are executable.
[[ -x "$CC" && -x "$CXX" && -x "$FC" ]] || { echo "ERROR: mpicc/mpicxx/mpifort not found"; exit 3; }
# Clear the shell's command hash table to ensure it finds the newly defined compilers from the conda environment.
hash -r

# Prepend the conda environment path to CMAKE_PREFIX_PATH to help CMake find libraries and headers.
export CMAKE_PREFIX_PATH="$CONDA_PREFIX ${CMAKE_PREFIX_PATH:-}"
# Prepend the conda pkg-config path to PKG_CONFIG_PATH.
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
# Add the conda include directory to the C preprocessor flags.
export CPPFLAGS="-I$CONDA_PREFIX/include ${CPPFLAGS:-}"
# Add the conda library path to linker flags (-L) and set the runtime path (-rpath) for the linker.
export LDFLAGS="-Wl,-rpath,${CONDA_PREFIX}/lib -L$CONDA_PREFIX/lib ${LDFLAGS:-}"

# Define paths to the BLAS and LAPACK libraries within the conda environment.
BLAS_LIBRARIES="$CONDA_PREFIX/lib/libblas.so"
LAPACK_LIBRARIES="$CONDA_PREFIX/lib/liblapack.so"
# Check if the BLAS and LAPACK shared libraries exist.
[[ -f "$BLAS_LIBRARIES" && -f "$LAPACK_LIBRARIES" ]] || { echo "ERROR: BLAS/LAPACK not found"; exit 4; }

# Define custom Fortran flags for the UKRmol+ build.
# -fdefault-integer-8: Promotes default INTEGERs to 8 bytes (important for large calculations).
# -std=f2018: Use the Fortran 2018 standard.
# -O3: High level of compiler optimization.
# -fopenmp: Enable OpenMP for multi-threading.
UKRMOL_FFLAGS="-fdefault-integer-8 -std=f2018 -O3 -fopenmp"

# Find the MPI execution command, trying 'mpirun' first, then 'mpiexec'.
MPIEXEC="$(command -v mpirun || command -v mpiexec || echo mpiexec)"
# Set pre-flags for mpiexec, needed for some cluster/container environments. Semicolon-separated for CMake.
MPI_PREFLAGS="--oversubscribe;--allow-run-as-root"
# Determine the number of processors to use for parallel builds ('make -j'). Defaults to 4 if detection fails.
NPROC="${NPROC:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 4)}"

# Define a working directory for third-party source code and builds.
WORKDIR="$REPO_ROOT/_thirdparty"
# Create the working directory if it doesn't exist.
mkdir -p "$WORKDIR"

# ---------- ScaLAPACK (shared; tests OFF to silence warnings) ----------
# Set the ScaLAPACK version to use. Can be overridden by an environment variable.
SCALAPACK_VER="${SCALAPACK_VER:-2.2.2}"
# Construct the download URL for the ScaLAPACK source code.
SCA_URL="https://codeload.github.com/Reference-ScaLAPACK/scalapack/tar.gz/refs/tags/v${SCALAPACK_VER}"

# Change to the working directory.
cd "$WORKDIR"
# Announce the download.
echo "==> Fetching ScaLAPACK ${SCALAPACK_VER} …"
# Define the filename for the downloaded tarball.
SCA_TGZ="scalapack-${SCALAPACK_VER}.tar.gz"
# Call the fetch_tgz helper function and exit if the download or validation fails.
fetch_tgz "$SCA_URL" "$SCA_TGZ" || { echo "ERROR: invalid ScaLAPACK tarball"; exit 6; }

# Create a clean directory for the ScaLAPACK source.
SCA_SRC="$WORKDIR/scalapack-src"; rm -rf "$SCA_SRC"; mkdir -p "$SCA_SRC"
# Extract the tarball into the source directory.
tar xzf "$SCA_TGZ" -C "$SCA_SRC"
# Find the name of the top-level directory inside the extracted tarball.
SCA_TOP="$(find "$SCA_SRC" -mindepth 1 -maxdepth 1 -type d -print -quit)"
# Exit if the top-level directory could not be found.
[[ -n "$SCA_TOP" ]] || { echo "ERROR: cannot find ScaLAPACK top directory"; exit 7; }

# Path to the BLACS CMake file, which is a dependency of ScaLAPACK.
BLACS_INSTALL_DIR="$SCA_TOP/BLACS/INSTALL"
# Check if the BLACS CMake file specifies an old, unsupported CMake version (e.g., 2.x).
if [[ -f "$BLACS_INSTALL_DIR/CMakeLists.txt" ]] && grep -q -E 'cmake_minimum_required\(VERSION *2\.' "$BLACS_INSTALL_DIR/CMakeLists.txt"; then
  # If so, patch the file to require a more modern version (3.5) to ensure compatibility.
  echo "==> Patching BLACS/INSTALL CMakeLists.txt (min version -> 3.5)"
  sed -i 's/cmake_minimum_required(VERSION [0-9]\+\.[0-9]\+)/cmake_minimum_required(VERSION 3.5)/' "$BLACS_INSTALL_DIR/CMakeLists.txt"
fi

# Announce the start of the build process.
echo "==> Building ScaLAPACK (${SCALAPACK_VER}) …"
# Create a build directory and change into it (out-of-source build).
mkdir -p "$SCA_TOP/build" && cd "$SCA_TOP/build"
# Configure the build using CMake.
cmake -S .. -B . \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_SHARED_LIBS=ON \
  -D BUILD_TESTING=OFF \
  -D CMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -D CMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -D CMAKE_C_COMPILER="$CC" \
  -D CMAKE_CXX_COMPILER="$CXX" \
  -D CMAKE_Fortran_COMPILER="$FC" \
  -D BLAS_LIBRARIES="$BLAS_LIBRARIES" \
  -D LAPACK_LIBRARIES="$LAPACK_LIBRARIES" \
  -D CMAKE_C_STANDARD=90 \
  -D CMAKE_C_EXTENSIONS=ON \
  -D CMAKE_C_FLAGS="-std=gnu89 -Wno-implicit-function-declaration"

# Build the project in parallel using $NPROC cores.
cmake --build . -j "$NPROC"
# Install the compiled files into the CONDA_PREFIX. '|| true' prevents exit on non-critical install errors.
cmake --install . || true

# Define the expected path of the installed ScaLAPACK shared library.
SCALAPACK_SO="$CONDA_PREFIX/lib/libscalapack.so"
# Check if the library was installed successfully.
[[ -f "$SCALAPACK_SO" ]] || { echo "ERROR: libscalapack.so not installed"; exit 8; }
# Confirm successful installation.
echo "==> ScaLAPACK installed: $SCALAPACK_SO"

# ---------- ELPA (generic kernels ONLY) ----------
# Set the ELPA version. Can be overridden by an environment variable.
ELPA_VER="${ELPA_VER:-2025.06.001}"

# Define a function to build ELPA.
ELPA_build_generic() {
  # Assign the version argument to a local variable.
  local ver="$1"
  # Construct the download URL for the ELPA source code.
  local url="https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/${ver}/elpa-${ver}.tar.gz"
  # Define local variables for file and directory names.
  local tgz="elpa-${ver}.tar.gz"
  local srcd="elpa-${ver}"
  # Define the installation prefix inside the conda env (e.g., /path/to/env/opt/elpa-...).
  local prefix="$CONDA_PREFIX/opt/elpa-${ver}"

  # Create a source directory for ELPA and change into it.
  mkdir -p "$WORKDIR/elpa-src" && cd "$WORKDIR/elpa-src"
  # Announce the download.
  echo "==> Fetching ELPA ${ver} …"
  # Download and validate the ELPA tarball.
  fetch_tgz "$url" "$tgz" || { echo "ERROR: invalid ELPA tarball for ${ver}"; return 9; }

  # Remove any old source directory and extract the new one.
  rm -rf "$srcd" && tar xzf "$tgz"
  cd "$srcd"
  # Create a clean build directory and change into it.
  rm -rf build && mkdir build && cd build

  # Force bash for the configure script to avoid shell compatibility issues.
  export CONFIG_SHELL=/bin/bash
  export SHELL=/bin/bash
  # Unset Fortran flags to avoid them being picked up incorrectly by the configure script.
  unset FCFLAGS FFLAGS

  # Define compiler wrappers that include the necessary flags.
  local FC_WRAPPED="$(command -v mpifort) -O3 -fopenmp -fPIC -fallow-argument-mismatch"
  local CC_WRAPPED="$(command -v mpicc) -O3 -fopenmp -fPIC"
  local CXX_WRAPPED="$(command -v mpicxx)"
  # Define linker flags needed to link against ScaLAPACK, BLAS, and LAPACK.
  local SCALAPACK_LDFLAGS="-L${CONDA_PREFIX}/lib -lscalapack -llapack -lblas -Wl,-rpath,${CONDA_PREFIX}/lib"

  # Define flags to build only the generic, portable CPU kernels (no AVX, SSE, etc.).
  local KFLAGS="--enable-generic-kernels \
                --disable-sse --disable-sse-assembly \
                --disable-avx-kernels --disable-avx2-kernels --disable-avx512-kernels"

  # Announce the configuration step.
  echo "==> Configuring ELPA (${ver}) — GENERIC kernels only …"
  # Run the ELPA configure script with all options.
  bash ../configure \
    --prefix="${prefix}" \
    --libdir="${prefix}/lib" \
    CC="$CC_WRAPPED" CXX="$CXX_WRAPPED" FC="$FC_WRAPPED" \
    --with-mpi=yes --enable-openmp --enable-shared \
    --disable-mpi-module \
    --enable-allow-thread-limiting \
    --without-threading-support-check-during-build \
    --enable-option-checking=fatal \
    ${KFLAGS} \
    SCALAPACK_FCFLAGS="" \
    SCALAPACK_LDFLAGS="${SCALAPACK_LDFLAGS}"

  # Build in parallel. If it fails, retry with a single process which can help debug.
  make -j "${NPROC:-4}" || make
  # Install ELPA to the specified prefix.
  make install

  # Find the path to the installed ELPA shared library.
  local libdir="${prefix}/lib"
  local elpaso=
  for f in "$libdir"/libelpa_openmp.so "$libdir"/libelpa.so; do
    [[ -f "$f" ]] && { elpaso="$f"; break; }
  done
  # Exit if the library was not found.
  [[ -n "$elpaso" ]] || { echo "ERROR: ELPA lib not found in ${libdir}"; return 10; }

  # Find the paths to the installed ELPA Fortran module directories.
  local moddirs=()
  for d in "${prefix}"/include/elpa_openmp-*/modules "${prefix}"/include/elpa-*/modules \
           "${prefix}"/include/elpa_openmp/modules "${prefix}"/include/elpa/modules; do
    [[ -d "$d" ]] && moddirs+=("$d")
  done
  # Exit if the module directories were not found.
  [[ ${#moddirs[@]} -gt 0 ]] || { echo "ERROR: ELPA module dirs missing (${prefix})"; return 11; }

  # Set global variables with the results.
  ELPA_PREFIX="${prefix}"
  ELPA_SO="${elpaso}"
  ELPA_MOD_DIRS=("${moddirs[@]}")
  return 0
}

# Initialize variables for the ELPA build results.
ELPA_SO=""
ELPA_MOD_DIRS=()
ELPA_PREFIX=""
# Call the ELPA build function and exit if it fails.
ELPA_build_generic "$ELPA_VER" || { echo "ERROR: ELPA (generic) build failed for ${ELPA_VER}"; exit 12; }

# Create linker flags to add the ELPA library directory to the runtime path.
ELPA_RPATH_FLAGS="-Wl,-rpath,${ELPA_PREFIX}/lib"
# Create compiler include flags (e.g., "-I/path/to/mod1 -I/path/to/mod2") for the command line.
ELPA_I_FLAGS="$(printf ' -I%s' "${ELPA_MOD_DIRS[@]}")"
# Create a semicolon-separated list of include paths for use in CMake list variables.
ELPA_INCLUDE_DIRS="$(IFS=';'; echo "${ELPA_MOD_DIRS[*]}")"

# Print a summary of the ELPA installation.
echo "==> ELPA installed at: ${ELPA_PREFIX}"
echo "    ELPA library:       ${ELPA_SO}"
printf '    ELPA modules:      %s\n' "${ELPA_MOD_DIRS[@]}"

# ---------- UKRmol (double & quad) ----------
# Change back to the repository root directory.
cd "$REPO_ROOT"
# Define the base directory for the final UKRmol+ installations.
PREFIX_BASE="$REPO_ROOT/install"
# Define specific installation prefixes for the double and quad precision builds.
PDS="$PREFIX_BASE/double_SCALAPACK"
PQS="$PREFIX_BASE/quad_SCALAPACK"
# Create the installation directories.
mkdir -p "$PDS" "$PQS"

# Define preprocessor flags for the GBTOlib library to enable MPI features.
GBTO_MPI_FLAGS="-Dusempi;-Dsplitreduce;-Dmpithree"
# Combine the rpath flags for both ELPA and the conda environment.
RPATH_LINK_FLAGS="${ELPA_RPATH_FLAGS} -Wl,-rpath,${CONDA_PREFIX}/lib"

# Announce the configuration of the double precision build.
echo "==> Configuring UKRmol+ (double) …"
# Configure the double precision build with CMake.
cmake -S . -B build_double_SCALAPACK \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$PDS" \
  -D CMAKE_Fortran_FLAGS="$UKRMOL_FFLAGS ${ELPA_I_FLAGS}" \
  -D GBTOlib_Fortran_FLAGS="$GBTO_MPI_FLAGS" \
  -D UKRMOL_OUT_DIR="../UKRmol-out" \
  -D BLAS_LIBRARIES="$BLAS_LIBRARIES" \
  -D LAPACK_LIBRARIES="$LAPACK_LIBRARIES" \
  -D MPIEXEC_EXECUTABLE="$MPIEXEC" \
  -D MPIEXEC_PREFLAGS="$MPI_PREFLAGS" \
  -D WITH_SCALAPACK=ON \
  -D SCALAPACK_LIBRARIES="$SCALAPACK_SO" \
  -D WITH_ELPA=ON \
  -D ELPA_LIBRARIES="$ELPA_SO" \
  -D ELPA_INCLUDE_DIRS="$ELPA_INCLUDE_DIRS" \
  -D CMAKE_EXE_LINKER_FLAGS="$RPATH_LINK_FLAGS" \
  -D CMAKE_SHARED_LINKER_FLAGS="$RPATH_LINK_FLAGS"

# Build the double precision version in parallel.
cmake --build build_double_SCALAPACK -j "$NPROC"
# Install the double precision version.
cmake --install build_double_SCALAPACK

# Announce the configuration of the quad precision build.
echo "==> Configuring UKRmol+ (quad) …"
# Configure the quad precision build with CMake.
cmake -S . -B build_quad_SCALAPACK \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$PQS" \
  -D CMAKE_Fortran_FLAGS="$UKRMOL_FFLAGS ${ELPA_I_FLAGS}" \
  -D GBTOlib_Fortran_FLAGS="${GBTO_MPI_FLAGS};-Dusequadprec" \
  -D UKRMOL_OUT_DIR="../UKRmol-out" \
  -D BLAS_LIBRARIES="$BLAS_LIBRARIES" \
  -D LAPACK_LIBRARIES="$LAPACK_LIBRARIES" \
  -D MPIEXEC_EXECUTABLE="$MPIEXEC" \
  -D MPIEXEC_PREFLAGS="$MPI_PREFLAGS" \
  -D WITH_SCALAPACK=ON \
  -D SCALAPACK_LIBRARIES="$SCALAPACK_SO" \
  -D WITH_ELPA=ON \
  -D ELPA_LIBRARIES="$ELPA_SO" \
  -D ELPA_INCLUDE_DIRS="$ELPA_INCLUDE_DIRS" \
  -D CMAKE_EXE_LINKER_FLAGS="$RPATH_LINK_FLAGS" \
  -D CMAKE_SHARED_LINKER_FLAGS="$RPATH_LINK_FLAGS"

# Build the quad precision version in parallel.
cmake --build build_quad_SCALAPACK -j "$NPROC"
# Install the quad precision version.
cmake --install build_quad_SCALAPACK

# test psi4
# Announce the test of the Psi4 installation.
echo "==> Testing PSI4 installation …"
# Execute a short Python script to test the Psi4 import.
python - <<'PY'
import sys
try:
    import psi4
    print("Psi4 OK:", psi4.__version__)
except Exception as e:
    print("Psi4 import failed:", e, file=sys.stderr)
    sys.exit(1)
PY

# Deactivate the conda environment. '|| true' prevents an error if it's already inactive.
conda deactivate || true
# Print a blank line for readability.
echo
# Print a summary of the completed build and installation paths.
echo "BUILD COMPLETED (generic ELPA)"
echo "   • ScaLAPACK: ${SCALAPACK_SO}"
echo "   • ELPA     : ${ELPA_SO}"
echo "   • Double   : ${PDS}"
echo "   • Quad     : ${PQS}"
echo "   * Psi4     : Tested successfully"