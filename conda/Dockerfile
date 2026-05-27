# WHAT IT DOES:
#   Creates a portable, self-contained environment for UKRmol+ by compiling 
#   the entire toolchain from source. It ensures 100% reproducibility by not 
#   relying on system package managers for scientific libraries.
#
# WHAT IS DOWNLOADED & COMPILED:
#   1. Compiler Suite: GCC 15.2.0 (gfortran, g++, gcc)
#   2. MPI Stack:      OpenMPI 5.0.9 + hwloc + libevent + zlib
#   3. Math Libs:      OpenBLAS 0.3.30, ScaLAPACK 2.2.2, ELPA 2025.01, 
#                      PETSc 3.24.1, SLEPc 3.24.1, ARPACK-NG 3.9.1
#   4. Chemistry:      Psi4 1.10 (plus Python 3.11, NumPy 2.3.4, PCMSolver)
#   5. UKRmol+:        Compiled in both Double and Quadruple precision.
#
# HOW TO RUN (CRITICAL):
#   You must run the build from the PARENT directory containing both 
#   'ukrmol-in' and 'ukrmol-out' repositories.
#
#   Structure:
#     my_project/
#      ├── ukrmol-in/   <-- Contains this Dockerfile
#      └── ukrmol-out/
#
#   Command:
#     cd my_project
#     docker build -f ukrmol-in/conda/Dockerfile -t ukrmol-plus:latest .
# ============================================================================ #

# ============================================================================ #
# STAGE 1: BUILDER
# ============================================================================ #
# We use a separate 'builder' stage to compile all software. 
# This allows us to discard all source code, compilers (GCC), and intermediate 
# build artifacts in the final stage, keeping the production image lightweight.
# ============================================================================ #
FROM ubuntu:22.04 AS builder

LABEL stage="builder"
# Prevents apt-get from asking for user input (timezone, keyboard layout, etc.)
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install Bootstrap Dependencies
# These are the minimum tools required to run the 'docker_build_script.sh'.
# We need curl/git to download sources and make/g++ to bootstrap the build tools.
RUN echo "--- Installing bootstrap dependencies ---" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    g++ make unzip pkg-config doxygen graphviz curl git file ca-certificates xz-utils \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# 2. Copy Source Repositories
# These are copied explicitly to ensure the build script can find them.
COPY ukrmol-in ./ukrmol-in-git
COPY ukrmol-out ./ukrmol-out-git
COPY ukrmol-in/conda/docker_build_script.sh .

# 3. CACHE INJECTION (CRITICAL STEP)
# This is where the magic happens. We copy the contents of the local 'ukrmol-cache' 
# folder (on your host machine) into the Docker build context.
#
# - If 'ukrmol-cache' contains tarballs (e.g., gmp-6.3.0.tar.xz), the script skips downloading them.
# - If it contains compiled folders (e.g., auxiliary-software/), the script skips compiling them.
#
# This dramatically speeds up re-builds if you use the 'cache_ukrmol_docker_build.sh' helper.
COPY ukrmol-cache/ ./

# 4. Run the Master Build Script
# The script checks the presence of files copied in Step 3. 
# If found, it skips steps; otherwise, it downloads and compiles from scratch.
RUN chmod +x docker_build_script.sh && \
    ./docker_build_script.sh

# ============================================================================ #
# STAGE 2: FINAL IMAGE
# ============================================================================ #
# This stage creates the actual distributable image. It starts fresh (empty) 
# and we only COPY the final binaries from the 'builder' stage.
# ============================================================================ #
FROM ubuntu:22.04 AS final

LABEL description="UKRmol+ suite"

# 1. Install Runtime Dependencies
# These are libraries required to RUN the binaries (OpenMP, Fortran runtime, etc.),
# distinct from the tools needed to BUILD them.
RUN echo "--- Installing runtime dependencies ---" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    libgomp1 libquadmath0 vim nano less openssh-client \
    && rm -rf /var/lib/apt/lists/*

# 2. Copy Compiled Software
# We pull only the 'ukrmolp-release' (the suite) and 'auxiliary-software' (GCC libs)
# from the builder stage. All source code and build temp files are left behind.
COPY --from=builder /build/ukrmolp-release /opt/ukrmolp
COPY --from=builder /build/auxiliary-software /opt/compiler

# 3. Setup Shell Environment
# We install the custom .bashrc that handles the 'use_double' / 'use_quad' aliases.
COPY ukrmol-in/conda/ukrmol.bashrc /opt/ukrmolp/ukrmol.bashrc
RUN echo "source /opt/ukrmolp/ukrmol.bashrc" >> /etc/bash.bashrc && \
    chmod 644 /opt/ukrmolp/ukrmol.bashrc

# 4. Configure Entrypoint
# Ensures the environment is loaded correctly whether running interactively or not.
RUN echo '#!/bin/bash' > /opt/ukrmolp/entrypoint.sh && \
    echo 'source /opt/ukrmolp/ukrmol.bashrc' >> /opt/ukrmolp/entrypoint.sh && \
    echo 'exec "$@"' >> /opt/ukrmolp/entrypoint.sh && \
    chmod +x /opt/ukrmolp/entrypoint.sh

# 5. Create Non-Root User
# Good practice for HPC usage (Apptainer compatibility).
RUN useradd --create-home --shell /bin/bash ukrmol
USER ukrmol
WORKDIR /home/ukrmol

# 6. Define Environment Variables
# - OMPI_*: Allow OpenMPI to run as root (sometimes needed in containers).
# - PATH/LD_LIBRARY_PATH: Set up to find our custom-built compilers and libraries.
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    PRTE_ALLOW_RUN_AS_ROOT=1 \
    PRTE_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OPAL_PREFIX="/opt/ukrmolp" \
    CC="/opt/compiler/bin/gcc" \
    CXX="/opt/compiler/bin/g++" \
    FC="/opt/compiler/bin/gfortran" \
    PATH="/opt/ukrmolp/bin.double:/opt/ukrmolp/bin:/opt/compiler/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/ukrmolp/lib.double:/opt/ukrmolp/lib:/opt/ukrmolp/lib64:/opt/compiler/lib64:/opt/compiler/lib:${LD_LIBRARY_PATH}" \
    PYTHONPATH="/opt/ukrmolp/lib/python3.11/site-packages:/opt/ukrmolp/lib64/python3.11/site-packages:${PYTHONPATH}"

ENTRYPOINT ["/opt/ukrmolp/entrypoint.sh"]
CMD ["/bin/bash"]