#!/bin/bash
# ---------------------------------------------------------------------------- #
# Builds a portable binary distribution of UKRmol+                             #
# ---------------------------------------------------------------------------- #
# This script will download, patch, compile and install UKRmol+ and its deps.  #
# The compilation relies on some external development packages (C/C++/Fortran  #
# compilers, Autotools, CMake), which are expected to be present. The          #
# resulting directory 'ukrmolp-release' will be self-contained and portable.   #
# ---------------------------------------------------------------------------- #

# parameters
DIRNAME=ukrmolp-release             # name of the directory with the software
AUXDIR=$(pwd)/auxiliary-software    # SW only used for building
INSTDIR=$(pwd)/$DIRNAME             # SW that is part of the release
ARCH=x86-64-v3                      # target instruction set (AVX2)
NPROC=8                             # number of cores for compilation and tests
RUN_TESTS=false                     # execute UKRmol+ tests after compilation
CURL="env -i curl"                  # avoid this environment when downloading

# package versions
arpack_version=3.9.1
cmake_version=4.2.0
cython_version=3.2.0
eigen_version=3.4.1
elpa_version=2025.01.002
gcc_version=15.2.0
gmp_version=6.3.0
hwloc_version=2.12.2
libevent_version=2.1.12
libffi_version=3.5.2
mpc_version=1.3.1
mpfr_version=4.2.2
numpy_version=2.3.4
ninja_version=1.13.1
openblas_version=0.3.30
openclicd_version=2025.07.22
openmpi_version=5.0.9
openmpi_ver=5.0
openssl_version=3.5.4
patchelf_version=0.15.5
petsc_version=3.24.1
psi4_version=1.10
pybind11_version=2.13.6
python_version=3.11.4
python_ver=3.11
scalapack_version=2.2.2
setuptools_version=80.9.0
slepc_version=3.24.1
zlib_version=1.3.1

# 0. Compiler
# -----------
# A recent version of the GNU compiler suite to build the release

echo "--- Building UkRmol+ container ---"
echo "The script is using cashed files if they are present in the 'ukrmol-cache' folder."
echo "To force recompilation of a specific package, please delete the relevant files from the 'ukrmol-cache' folder."

if [ ! -f $AUXDIR/bin/gcc ]
then
    # preparations
    mkdir -p $AUXDIR/bin
    $CURL -OLs https://github.com/ninja-build/ninja/releases/download/v$ninja_version/ninja-linux.zip || exit 1
    unzip ninja-linux.zip || exit 1
    mv -v ninja $AUXDIR/bin/
    $CURL -Ls https://github.com/Kitware/CMake/releases/download/v$cmake_version/cmake-$cmake_version-linux-x86_64.tar.gz | tar xz || exit 1
    cp -rv cmake-$cmake_version-linux-x86_64/* $AUXDIR/
    # patchelf removed as it is not needed for fixed-path Docker builds
    
    # download all components
    $CURL -Ls https://ftp.gwdg.de/pub/misc/gcc/releases/gcc-$gcc_version/gcc-$gcc_version.tar.xz | tar xJ || exit 1
    pushd gcc-$gcc_version

    
    # GMP
    if [ ! -d gmp-$gmp_version ]; then
        echo "--- Downloading GMP $gmp_version ---"
        $CURL -fL --retry 5 --retry-delay 5 \
            https://gmplib.org/download/gmp/gmp-$gmp_version.tar.xz \
            -o gmp-$gmp_version.tar.xz || exit 1

        xz -t gmp-$gmp_version.tar.xz || exit 1
        tar xJf gmp-$gmp_version.tar.xz || exit 1
    else 
        echo "GMP found in cache -> SKIPING DOWNLOAD"
    fi
    rm -rf gmp && mv gmp-$gmp_version gmp

    # MPFR
    if [ ! -d mpfr-$mpfr_version ]; then
        echo "--- Downloading MPFR $mpfr_version ---"
        $CURL -fL --retry 5 --retry-delay 5 \
            https://www.mpfr.org/mpfr-$mpfr_version/mpfr-$mpfr_version.tar.xz \
            -o mpfr-$mpfr_version.tar.xz || exit 1

        xz -t mpfr-$mpfr_version.tar.xz || exit 1
        tar xJf mpfr-$mpfr_version.tar.xz || exit 1
    else 
        echo "MPFR found in cache -> SKIPING DOWNLOAD"
    fi
    rm -rf mpfr && mv mpfr-$mpfr_version mpfr


    $CURL -Ls https://ftp.gnu.org/gnu/mpc/mpc-$mpc_version.tar.gz | tar xz || exit 1
    rm -rf mpc && mv mpc-$mpc_version mpc
    popd
    # build the compiler suite
    mkdir -p gcc-$gcc_version/build
    pushd gcc-$gcc_version/build
    ../configure \
        --prefix=$AUXDIR \
        --disable-multilib \
        --enable-languages=fortran \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else
    echo "GCC found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

export PATH="$AUXDIR/bin:$INSTDIR/bin:$PATH"
export LD_LIBRARY_PATH="$AUXDIR/lib64:$INSTDIR/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$INSTDIR/lib/pkgconfig"
export CC="$(which gcc)"
export CXX="$(which g++)"
export FC="$(which gfortran)"

# copy compiler runtime libraries to the staging directory

echo "--- Copying compiler runtime libraries ---"
mkdir -p $INSTDIR/lib
cp -Lv \
    $AUXDIR/lib64/libgomp.so.1 \
    $AUXDIR/lib64/libstdc++.so.6 \
    $AUXDIR/lib64/libgcc_s.so.1 \
    $AUXDIR/lib64/libgfortran.so.5 \
    $AUXDIR/lib64/libquadmath.so.0 \
    $INSTDIR/lib

# 1. zlib
# -------
# Basic compression library (dependency of Open MPI)

if [ ! -f $INSTDIR/lib/libz.so ]
then
    $CURL -Ls https://zlib.net/zlib-$zlib_version.tar.gz | tar xz || exit 1
    mkdir -p zlib-$zlib_version/build
    pushd zlib-$zlib_version/build
    ../configure \
        --prefix=$INSTDIR \
        --libdir=$INSTDIR/lib \
        --enable-shared \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else
    echo "zlib found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 2. libevent
# -----------
# Event library (dependency of Open MPI)

if [ ! -f $INSTDIR/lib/libevent.so ]
then
    $CURL -Ls https://github.com/libevent/libevent/releases/download/release-$libevent_version-stable/libevent-$libevent_version-stable.tar.gz | tar xz || exit 1
    mkdir -p libevent-$libevent_version-stable/build
    pushd libevent-$libevent_version-stable/build
    ../configure \
        --prefix=$INSTDIR \
        --libdir=$INSTDIR/lib \
        --enable-shared \
        --disable-openssl \
        --disable-static \
        CFLAGS="-march=$ARCH" \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else 
    echo "libevent found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 3. hwloc
# --------
# Hardware locality library (dependency of Open MPI)

if [ ! -f $INSTDIR/lib/libhwloc.so ]
then
    $CURL -Ls https://download.open-mpi.org/release/hwloc/v2.12/hwloc-$hwloc_version.tar.gz | tar xz || exit 1
    mkdir -p hwloc-$hwloc_version/build
    pushd hwloc-$hwloc_version/build
    ../configure \
        --prefix=$INSTDIR \
        --libdir=$INSTDIR/lib \
        --disable-cairo \
        --disable-io \
        --disable-libxml2 \
        CFLAGS="-march=$ARCH" \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else 
    echo "hwloc found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 4. Open MPI
# -----------
# Implementation of MPI 3.1

if [ ! -f $INSTDIR/lib/libmpi.so ]
then
    $CURL -Ls https://download.open-mpi.org/release/open-mpi/v$openmpi_ver/openmpi-$openmpi_version.tar.gz | tar xz || exit 1
    mkdir -p openmpi-$openmpi_version/build
    pushd openmpi-$openmpi_version/build
    ../configure \
        --prefix=$INSTDIR \
        --libdir='${exec_prefix}/lib' \
        --with-ucx=no \
        --with-hcoll=no \
        --with-hwloc=$INSTDIR \
        --with-libevent=$INSTDIR \
        --with-ofi=no \
        --with-psm2=no \
        --with-zlib=$INSTDIR \
        --enable-mpi-fortran=usempi \
        --disable-oshmem \
        --disable-prte-prefix-by-default \
        --disable-sphinx \
        CFLAGS="-mlong-double-128 -march=$ARCH" \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else
    echo "Open MPI found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 5. OpenBLAS
# -----------
# BLAS/LAPACK implementation

if [ ! -f $INSTDIR/lib/libopenblas.so ]
then
    $CURL -Ls https://github.com/OpenMathLib/OpenBLAS/releases/download/v$openblas_version/OpenBLAS-$openblas_version.tar.gz | tar xz || exit 1
    pushd OpenBLAS-$openblas_version
    make PREFIX=$INSTDIR DYNAMIC_ARCH=1 USE_THREAD=1 USE_OPENMP=1 NUM_THREADS=128 NO_STATIC=1 || exit 1
    make PREFIX=$INSTDIR DYNAMIC_ARCH=1 USE_THREAD=1 USE_OPENMP=1 NUM_THREADS=128 NO_STATIC=1 install || exit 1
    popd
else
    echo "OpenBLAS found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 6. ScaLAPACK
# ------------
# Reference ScaLAPACK implementation

if [ ! -f $INSTDIR/lib/libscalapack.so ]
then
    $CURL -Ls https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v$scalapack_version.tar.gz | tar xz || exit 1
    sed -i 's/PILAENV = 32/PILAENV = 128/g' $scalapack_pkg/PBLAS/SRC/pilaenv.f
    mkdir -p scalapack-$scalapack_version/build
    pushd scalapack-$scalapack_version/build
    CMAKE_POLICY_VERSION_MINIMUM=3.5 cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_LIBDIR=lib \
        -D CMAKE_POLICY_DEFAULT_CMP0074=NEW \
        -D CMAKE_C_FLAGS="-std=gnu89 -march=$ARCH" \
        -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -march=$ARCH" \
        -D BLAS_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D LAPACK_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D MPI_ROOT=$INSTDIR \
        -D BUILD_SHARED_LIBS=ON \
        .. || exit 1
    cmake --build . || exit 1
    cmake --install . || exit 1
    popd
else
    echo "ScaLAPACK found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 7. ELPA
# -------
# High-performance parallel dense diagonalizer

if [ ! -f $INSTDIR/lib/libelpa.so ]
then
    $CURL -Ls https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/$elpa_version/elpa-$elpa_version.tar.gz | tar xz || exit 1
    mkdir -p elpa-$elpa_version/build
    pushd elpa-$elpa_version/build
    ../configure \
        --prefix=$INSTDIR \
        --libdir=$INSTDIR/lib \
        --disable-avx512-kernels \
        --disable-static \
        --with-test-programs=no \
        CC=$INSTDIR/bin/mpicc \
        FC=$INSTDIR/bin/mpifort \
        CFLAGS="-march=$ARCH" \
        LDFLAGS="$INSTDIR/lib/libscalapack.so $INSTDIR/lib/libopenblas.so" \
        || exit 1
    make -j 1 && make install || exit 1
    popd
else
    echo "ELPA found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 8. Arpack-NG
# ------------
# Iterative eigenvalue solver

if [ ! -f $INSTDIR/lib/libarpack.so ]
then
    $CURL -Ls https://github.com/opencollab/arpack-ng/archive/refs/tags/$arpack_version.tar.gz | tar xz || exit 1
    mkdir -p arpack-ng-$arpack_version/build
    pushd arpack-ng-$arpack_version/build
    cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_LIBDIR=lib \
        -D BLAS_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D LAPACK_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D BUILD_SHARED_LIBS=ON \
        .. || exit 1
    cmake --build . || exit 1
    cmake --install . || exit 1
    popd
else
    echo "Arpack-NG found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 11. libffi
# ----------
# Foreign function interfaces (dependency of Python)

if [ ! -f $INSTDIR/lib/libffi.so ]
then
    $CURL -Ls https://github.com/libffi/libffi/releases/download/v$libffi_version/libffi-$libffi_version.tar.gz | tar xz || exit 1
    mkdir -p libffi-$libffi_version/build
    pushd libffi-$libffi_version/build
    ../configure \
        --prefix=$INSTDIR \
        --libdir=$INSTDIR/lib \
        --enable-shared \
        --disable-multi-os-directory \
        --disable-static \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else
    echo "libffi found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 12. OpenSSL
# -----------
# Secure sockets library (dependency of Python)

if [ ! -f $INSTDIR/lib/libssl.so ]
then
    $CURL -Ls https://github.com/openssl/openssl/releases/download/openssl-$openssl_version/openssl-$openssl_version.tar.gz | tar xz || exit 1
    mkdir -p openssl-$openssl_version/build
    pushd openssl-$openssl_version/build
    ../Configure \
        --prefix=$INSTDIR \
        --libdir=lib \
        --openssldir=$INSTDIR/etc/ssl \
        --with-zlib-include=$INSTDIR/include \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else
    echo "OpenSSL found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 12. Python
# ----------
# The reference Python interpreter (a dependency of Psi4)

if [ ! -f $INSTDIR/lib/libpython3.so ]
then
    $CURL -Ls https://www.python.org/ftp/python/$python_version/Python-$python_version.tgz | tar xz || exit 1
    mkdir -p Python-$python_version/build
    pushd Python-$python_version/build
    ../configure \
        --prefix=$INSTDIR \
        --libdir=$INSTDIR/lib \
        --enable-shared \
        CFLAGS="-march=$ARCH" \
        || exit 1
    make -j $NPROC && make install || exit 1
    popd
else
    echo "Python found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 13. Setuptools
# --------------
# Python package management (a dependency of Cython)

if [ ! -d $INSTDIR/lib/python$python_ver/site-packages/setuptools ]
then
    $CURL -Ls https://github.com/pypa/setuptools/archive/refs/tags/v$setuptools_version.tar.gz | tar xz || exit 1
    pushd setuptools-$setuptools_version
    python3 setup.py build || exit 1
    python3 setup.py install || exit 1
    popd
else
    echo "Setuptools found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 14. Eigen
# ---------
# Linear algebra library (a dependency of Psi4)

if [ ! -d $INSTDIR/include/eigen3 ]
then
    $CURL -Ls https://gitlab.com/libeigen/eigen/-/archive/$eigen_version/eigen-$eigen_version.tar.gz | tar xz || exit 1
    mkdir -p eigen-$eigen_version/build
    pushd eigen-$eigen_version/build
    cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_LIBDIR=lib \
        -D CMAKE_CXX_FLAGS="-march=$ARCH" \
        -D BUILD_SHARED_LIBS=ON \
        -D EIGEN_BUILD_BLAS=OFF \
        -D EIGEN_BUILD_LAPACK=OFF \
        .. || exit 1
    cmake --build . || exit 1
    cmake --install . || exit 1
    popd
else
    echo "Eigen found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 15. Cython
# ----------
# Python-to-C compiler (a dependency of NumPy)

if [ ! -f $INSTDIR/bin/cython ]
then
    $CURL -Ls https://github.com/cython/cython/releases/download/$cython_version/cython-$cython_version.tar.gz | tar xz || exit 1
    pushd cython-$cython_version
    python3 setup.py build || exit 1
    python3 setup.py install || exit 1
    popd
else
    echo "Cython found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 16. NumPy
# ---------
# Python numerical library (a dependency of Psi4)

if [ ! -d $INSTDIR/lib/python$python_ver/site-packages/numpy ]
then
    $CURL -Ls https://github.com/numpy/numpy/releases/download/v$numpy_version/numpy-$numpy_version.tar.gz | tar xz || exit 1
    pushd numpy-$numpy_version
    python3 vendored-meson/meson/meson.py setup --prefix $INSTDIR --libdir $INSTDIR/lib build || exit 1
    python3 vendored-meson/meson/meson.py compile -C build || exit 1
    python3 vendored-meson/meson/meson.py install -C build || exit 1
    popd
else
    echo "NumPy found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 17. Pybind11
# ------------
# Python-C++ interface (a dependency of Psi4)

if [ ! -d $INSTDIR/include/pybind11 ]
then
    $CURL -Ls https://github.com/pybind/pybind11/archive/refs/tags/v$pybind11_version.tar.gz | tar xz || exit 1
    mkdir -p pybind11-$pybind11_version/build
    pushd pybind11-$pybind11_version/build
    cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_LIBDIR=lib \
        .. || exit 1
    cmake --build . || exit 1
    cmake --install . || exit 1
    popd
else
    echo "Pybind11 found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 18. Pydantic
# ------------
# (a run-time dependency of Psi4)

if [ ! -d $INSTDIR/lib/python$python_ver/site-packages/pydantic ]
then
    pip3 install pydantic || exit 1

else
    echo "Pydantic found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 19. Pint
# --------
# (a runtime dependency of Psi4)

if [ ! -d $INSTDIR/lib/python$python_ver/site-packages/pint ]
then
    pip3 install pint || exit 1

else
    echo "Pint found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 9. PetSc
# --------
# General linear algebra library

if [ ! -f $INSTDIR/lib/libpetsc.so ]
then
    $CURL -Ls https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$petsc_version.tar.gz | tar xz || exit 1
    pushd petsc-$petsc_version
    ./configure \
        --prefix=$INSTDIR \
        --with-blaslapack-lib="$INSTDIR/lib/libopenblas.so" \
        --with-debugging=0 \
        --with-mpi-dir=$INSTDIR \
        --with-x=0 \
        --COPTFLAGS="-O3 -march=$ARCH" \
        --CXXOPTFLAGS="-O3 -march=$ARCH" \
        --FOPTFLAGS="-O3 -march=$ARCH" \
        || exit 1
    make && make install || exit 1
    popd
else
    echo "PetSc found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 10. SLEPc
# ---------
# Eigenvalue solver based on PetSc

if [ ! -f $INSTDIR/lib/libslepc.so ]
then
    $CURL -Ls https://slepc.upv.es/download/distrib/slepc-$slepc_version.tar.gz | tar xz || exit 1
    pushd slepc-$slepc_version
    PETSC_DIR=$INSTDIR ./configure --prefix=$INSTDIR || exit 1
    PETSC_DIR=$INSTDIR make -j $NPROC || exit 1
    PETSC_DIR=$INSTDIR make install || exit 1
    popd
else
    echo "SLEPc found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 20. Psi4
# --------
# Quantum chemistry software

if [ ! -f $INSTDIR/bin/psi4 ]
then
    $CURL -Ls https://github.com/psi4/psi4/archive/refs/tags/v$psi4_version.tar.gz | tar xz || exit 1
    mkdir -p psi4-$psi4_version/build
    pushd psi4-$psi4_version/build
    cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_C_FLAGS="-march=$ARCH" \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_LIBDIR=lib \
        -D BLAS_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D LAPACK_LIBRARIES=$INSTDIR/lib/libopenblas.so \
        -D CMAKE_DISABLE_FIND_PACKAGE_pybynd11=ON \
        -D BUILD_SHARED_LIBS=ON \
        -D pybind11_DIR=$INSTDIR \
        -D Eigen3_DIR=$INSTDIR \
        -D ENABLE_XHOST=OFF \
        -D ENABLE_ecpint=ON \
        .. || exit 1
    # Psi4 cannot find Eigen3 (unlike libint2), we have to use CPLUS_INCLUDE_PATH
    CPLUS_INCLUDE_PATH=$INSTDIR/include cmake --build . || exit 1
    cmake --install . || exit 1
    # Psi4 hardcodes the interpreter path, we have to generalize it for relocatability
    sed -i "s;#!$INSTDIR/bin/python$python_ver;#!/usr/bin/env python$python_ver;g" $INSTDIR/bin/psi4
    popd
else
    echo "Psi4 found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 21. OpenCL headers
# ------------------
# Headers for the OpenCL C API

if [ ! -d $INSTDIR/include/CL ]
then
    $CURL -Ls https://github.com/KhronosGroup/OpenCL-Headers/archive/refs/tags/v2025.07.22.tar.gz | tar xz || exit 1
    mkdir -p OpenCL-Headers-$openclicd_version/build
    pushd OpenCL-Headers-$openclicd_version/build
    cmake \
        -G Ninja \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        .. || exit 1
    cmake --build . || exit 1
    cmake --install . || exit 1
    popd
else
    echo "OpenCL headers found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 22. OpenCL ICD
# --------------
# Generic OpenCL installable client driver

if [ ! -f $INSTDIR/lib/libOpenCL.so ]
then
    $CURL -Ls https://github.com/KhronosGroup/OpenCL-ICD-Loader/archive/refs/tags/v2025.07.22.tar.gz | tar xz || exit 1
    mkdir -p OpenCL-ICD-Loader-$openclicd_version/build
    pushd OpenCL-ICD-Loader-$openclicd_version/build
    cmake \
        -G Ninja \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_LIBDIR=lib \
        .. || exit 1
    cmake --build . || exit 1
    cmake --install . || exit 1
    popd
else
    echo "OpenCL ICD found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 23. CLBlast
# -----------
# OpenCL-based BLAS library

if [ ! -f $INSTDIR/lib/libclblast.so ]
then
    git clone https://github.com/CNugteren/CLBlast.git clblast-git
    mkdir -p clblast-git/build
    pushd clblast-git/build
    cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_INSTALL_PREFIX=$INSTDIR \
        -D CMAKE_INSTALL_LIBDIR=lib \
        -D TUNERS=OFF \
        .. || exit 1
    cmake --build . || exit 1
    cmake --install . || exit 1
    popd
else
    echo "CLBlast found in cache -> SKIPING DOWNLOAD AND BUILD"
fi

# 24. UKRmol+
# -----------
# Molecular R-matrix codes

if [ ! -d "ukrmol-in-git" ] || [ ! -d "ukrmol-out-git" ]; then
    echo "Error: UKRmol directories not found."
    exit 1
fi

# Allow Root User for OpenMPI 5.0+
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export PRTE_ALLOW_RUN_AS_ROOT=1
export PRTE_ALLOW_RUN_AS_ROOT_CONFIRM=1

# GBTO Flags explicit definition (matches Conda/CI logic)
GBTO_FLAGS="-Dusempi;-Dsplitreduce;-Dmpithree"

declare -A flags
flags[double]="-fdefault-integer-8 -march=$ARCH -Dusemapping -fno-signed-zeros -fno-trapping-math -fassociative-math"
flags[quad]="-fdefault-integer-8 -march=$ARCH -Dusequadprec"

for prec in double quad
do
    # ---------------------------------------------------------
    # 1. CONDITIONAL COMPILATION
    # Only runs if the directory does NOT exist
    # ---------------------------------------------------------
    if [ ! -d $INSTDIR/bin.$prec ]
    then
        mkdir -p ukrmol-in-git/build.$prec
        pushd ukrmol-in-git/build.$prec
        
        # We also define MPIEXEC explicitly to be safe
        MPI_EXEC_BIN=$(which mpiexec) 

        cmake \
            -G Ninja \
            -D CMAKE_BUILD_TYPE=Release \
            -D CMAKE_Fortran_FLAGS="${flags[$prec]}" \
            -D GBTOlib_Fortran_FLAGS="$GBTO_FLAGS" \
            -D CMAKE_INSTALL_BINDIR=bin.$prec \
            -D CMAKE_INSTALL_INCLUDEDIR=include/UKRmol+ \
            -D CMAKE_INSTALL_LIBDIR=lib.$prec \
            -D CMAKE_INSTALL_PREFIX=$INSTDIR \
            -D CMAKE_INSTALL_RPATH="\$ORIGIN/../lib.$prec:\$ORIGIN/../lib" \
            -D CMAKE_POLICY_DEFAULT_CMP0074=NEW \
            -D BUILD_DOC=ON \
            -D BUILD_SHARED_LIBS=ON \
            -D WITH_CLBLAST=ON \
            -D MPI_ROOT=$INSTDIR \
            -D ARPACK_LIBRARIES=$INSTDIR/lib/libarpack.so \
            -D BLAS_LIBRARIES=$INSTDIR/lib/libopenblas.so \
            -D LAPACK_LIBRARIES=$INSTDIR/lib/libopenblas.so \
            -D SCALAPACK_LIBRARIES=$INSTDIR/lib/libscalapack.so \
            -D ELPA_INCLUDE_DIRS=$INSTDIR/include/elpa-$elpa_version/modules \
            -D ELPA_LIBRARIES=$INSTDIR/lib/libelpa.so \
            -D SLEPC_INCLUDE_DIRS="$INSTDIR/include;$INSTDIR/include/petsc;$INSTDIR/include/petsc/finclude;$INSTDIR/include/slepc;$INSTDIR/include/slepc/finclude" \
            -D SLEPC_LIBRARIES="$INSTDIR/lib/libslepc.so;$INSTDIR/lib/libpetsc.so" \
            -D CLBLAST_INCLUDE_DIRS=$INSTDIR/include \
            -D CLBLAST_LIBRARIES=$INSTDIR/lib/libclblast.so \
            -D UKRMOL_OUT_DIR="../ukrmol-out-git" \
            -D MPIEXEC_EXECUTABLE=$MPI_EXEC_BIN \
            -D MPIEXEC_PREFLAGS="--allow-run-as-root;--oversubscribe" \
            .. || exit 1
            
        cmake --build . || exit 1
        cmake --install . || exit 1
        popd
    else
        echo "UKRmol+ ($prec) found in cache -> SKIPING BUILD"
        echo "To force recompilation, please delete the directory '$INSTDIR/bin.$prec' in ukrmol-cache and re-run this script."
    fi
    # END OF COMPILATION CHECK
    
    # ---------------------------------------------------------
    # 2. ALWAYS GENERATE ENV.BASH
    # This runs every time, even if files were cached
    # ---------------------------------------------------------
    
    # Ensure dir exists (in case it was manually deleted but build dir remains, though unlikely)
    mkdir -p "$INSTDIR/bin.$prec"

    echo "Generating env.bash for $prec..."
    ENV_FILE="$INSTDIR/bin.$prec/env.bash"
    
    cat <<EOF > "$ENV_FILE"
#!/bin/bash
DIR="\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )"
# lib.double or lib.quad
LIB_DIR="\$(dirname "\$DIR")/lib.$prec"
# Common lib for deps
COMMON_LIB="\$(dirname "\$DIR")/lib"

# Reset PATH/LD_LIBRARY_PATH (Remove old double/quad entries)
export PATH=\$(echo \$PATH | sed -e "s|/opt/ukrmolp/bin.double:||" -e "s|/opt/ukrmolp/bin.quad:||")
export LD_LIBRARY_PATH=\$(echo \$LD_LIBRARY_PATH | sed -e "s|/opt/ukrmolp/lib.double:||" -e "s|/opt/ukrmolp/lib.quad:||")

# Add current
export PATH="\$DIR:\$PATH"
export LD_LIBRARY_PATH="\$LIB_DIR:\$COMMON_LIB:\$LD_LIBRARY_PATH"

echo "Activated UKRmol+ ($prec)"
EOF
    chmod +x "$ENV_FILE"

done

if $RUN_TESTS
then
    pushd ukrmol-in-git/build.double
    OMP_NUM_THREADS=$(($NPROC / 2)) \
    PRTE_MCA_prte_hwloc_default_binding_policy=none \
        ctest -R parallel || exit 1
    popd
fi

# ============================================================================ #
# 26. Generate Version Information File
# ============================================================================ #
# We save the variables defined at the top of this script to a file
# so the .bashrc can read them at runtime.

VERSION_FILE="$INSTDIR/versions.env"
echo "Generating version file at $VERSION_FILE..."

cat <<EOF > "$VERSION_FILE"
export UK_GCC_VER="$gcc_version"
export UK_MPI_VER="$openmpi_version"
export UK_BLAS_VER="$openblas_version"
export UK_SCALAPACK_VER="$scalapack_version"
export UK_ELPA_VER="$elpa_version"
export UK_ARPACK_VER="$arpack_version"
export UK_PETSC_VER="$petsc_version"
export UK_SLEPC_VER="$slepc_version"
export UK_PSI4_VER="$psi4_version"
export UK_PYTHON_VER="$python_version"
export UK_NUMPY_VER="$numpy_version"
EOF

echo
echo "BUILD COMPLETED SUCCESSFULLY"
echo "The self-contained installation is ready in the './$DIRNAME' directory."
