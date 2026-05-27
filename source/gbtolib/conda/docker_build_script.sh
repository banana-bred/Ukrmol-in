#!/bin/bash
set -e

# ---------------------------------------------------------------------------- #
# GBTOlib Build (Standard LP64 / 4-byte Integers)                              #
# ---------------------------------------------------------------------------- #

# Configuration
INSTDIR=/opt/gbtolib
NPROC=$(nproc)
ARCH=x86-64-v3
CURL="curl -Ls"

# Versions (Matches UKRmol+)
zlib_version=1.3.1
libevent_version=2.1.12
hwloc_version=2.9.3
openmpi_version=5.0.9
openmpi_ver=5.0
openblas_version=0.3.26

export CC=gcc
export CXX=g++
export FC=gfortran

# ============================================================================ #
# 1. zlib
# ============================================================================ #
echo "--- Building zlib $zlib_version ---"
$CURL https://zlib.net/zlib-$zlib_version.tar.gz | tar xz
pushd zlib-$zlib_version
./configure --prefix=$INSTDIR --libdir=$INSTDIR/lib --enable-shared
make -j $NPROC
make install
popd
rm -rf zlib-$zlib_version

# ============================================================================ #
# 2. libevent
# ============================================================================ #
echo "--- Building libevent $libevent_version ---"
$CURL https://github.com/libevent/libevent/releases/download/release-$libevent_version-stable/libevent-$libevent_version-stable.tar.gz | tar xz
pushd libevent-$libevent_version-stable
./configure \
    --prefix=$INSTDIR \
    --libdir=$INSTDIR/lib \
    --enable-shared \
    --disable-openssl \
    --disable-static \
    CFLAGS="-march=$ARCH"
make -j $NPROC
make install
popd
rm -rf libevent-$libevent_version-stable

# ============================================================================ #
# 3. hwloc
# ============================================================================ #
echo "--- Building hwloc $hwloc_version ---"
$CURL https://download.open-mpi.org/release/hwloc/v${hwloc_version%.*}/hwloc-$hwloc_version.tar.gz | tar xz
pushd hwloc-$hwloc_version
./configure \
    --prefix=$INSTDIR \
    --libdir=$INSTDIR/lib \
    --disable-cairo \
    --disable-io \
    --disable-libxml2 \
    CFLAGS="-march=$ARCH"
make -j $NPROC
make install
popd
rm -rf hwloc-$hwloc_version

# ============================================================================ #
# 4. Open MPI 5.0.9 (Standard / LP64 Mode)
# ============================================================================ #
echo "--- Building OpenMPI $openmpi_version (Standard) ---"
$CURL https://download.open-mpi.org/release/open-mpi/v$openmpi_ver/openmpi-$openmpi_version.tar.gz | tar xz
pushd openmpi-$openmpi_version

./configure \
    --prefix=$INSTDIR \
    --libdir=$INSTDIR/lib \
    --with-hwloc=$INSTDIR \
    --with-libevent=$INSTDIR \
    --with-zlib=$INSTDIR \
    --disable-oshmem \
    --disable-static \
    --with-slurm \
    CFLAGS="-mlong-double-128 -march=$ARCH" \
    CXXFLAGS="-march=$ARCH" \
    FCFLAGS="-march=$ARCH"

make -j $NPROC
make install
popd
rm -rf openmpi-$openmpi_version

# ============================================================================ #
# 5. OpenBLAS
# ============================================================================ #

echo "--- Building OpenBLAS $openblas_version (Standard) ---"
$CURL https://github.com/OpenMathLib/OpenBLAS/releases/download/v$openblas_version/OpenBLAS-$openblas_version.tar.gz | tar xz
pushd OpenBLAS-$openblas_version

make -j $NPROC \
    DYNAMIC_ARCH=1 \
    USE_THREAD=1 \
    USE_OPENMP=1 \
    NUM_THREADS=64 \
    NO_STATIC=1

make PREFIX=$INSTDIR install
popd
rm -rf OpenBLAS-$openblas_version

# ============================================================================ #
# 6. GBTOlib Setup & Build
# ============================================================================ #
export PATH=$INSTDIR/bin:$PATH
export LD_LIBRARY_PATH=$INSTDIR/lib:$LD_LIBRARY_PATH
export CC=$INSTDIR/bin/mpicc
export CXX=$INSTDIR/bin/mpicxx
export FC=$INSTDIR/bin/mpifort

export FFLAGS="-fdefault-integer-8 -O3 -fopenmp -march=$ARCH"

echo "--- Building GBTOlib ---"
BASE_GBTO_FLAGS="-Dusempi;-Dsplitreduce;-Dmpithree"

for prec in double quad; do
    echo "--> Configuring $prec precision..."
    mkdir -p gbtolib-git/build.$prec
    pushd gbtolib-git/build.$prec

    if [ "$prec" == "quad" ]; then
        CURRENT_FLAGS="${BASE_GBTO_FLAGS};-Dusequadprec"
    else
        CURRENT_FLAGS="${BASE_GBTO_FLAGS}"
    fi

    cmake -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_Fortran_FLAGS="$FFLAGS" \
        -D GBTOlib_Fortran_FLAGS="$CURRENT_FLAGS" \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_BINDIR=bin.$prec \
        -D CMAKE_INSTALL_LIBDIR=lib.$prec \
        -D CMAKE_INSTALL_RPATH="\$ORIGIN/../lib:\$ORIGIN/../lib.$prec" \
        -D BUILD_SHARED_LIBS=ON \
        -D MPI_ROOT=$INSTDIR \
        -D BLAS_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D LAPACK_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D MPIEXEC_EXECUTABLE=$INSTDIR/bin/mpiexec \
        -D MPIEXEC_PREFLAGS="--allow-run-as-root;--oversubscribe" \
        ..

    cmake --build . --parallel $NPROC
    cmake --install .
    popd

    # Activation Script
    ENV_FILE="$INSTDIR/bin.$prec/env.bash"
    mkdir -p "$(dirname "$ENV_FILE")"
    cat <<EOF > "$ENV_FILE"
#!/bin/bash
export PATH=\$(echo \$PATH | sed -e "s|$INSTDIR/bin.double:||" -e "s|$INSTDIR/bin.quad:||")
export LD_LIBRARY_PATH=\$(echo \$LD_LIBRARY_PATH | sed -e "s|$INSTDIR/lib.double:||" -e "s|$INSTDIR/lib.quad:||")
export PATH="$INSTDIR/bin.$prec:\$PATH"
export LD_LIBRARY_PATH="$INSTDIR/lib.$prec:\$LD_LIBRARY_PATH"
echo -e "\033[1;32mActivated GBTOlib ($prec)\033[0m"
EOF
    chmod +x "$ENV_FILE"
done

# ============================================================================ #
# 7. Generate Version Information File
# ============================================================================ #
VERSION_FILE="$INSTDIR/versions.env"
echo "Generating version file at $VERSION_FILE..."

cat <<EOF > "$VERSION_FILE"
export GB_MPI_VER="$openmpi_version"
export GB_BLAS_VER="$openblas_version"
export GB_ZLIB_VER="$zlib_version"
export GB_HWLOC_VER="$hwloc_version"
export GB_LIBEVENT_VER="$libevent_version"
EOF

echo "GBTOlib Full Stack Build (OpenMPI 5.0 Standard) Complete."