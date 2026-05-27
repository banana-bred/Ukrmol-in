# ============================================================================ #
# GBTOlib Container (Debian Bookworm / OpenMPI 5.0.9 / Standard Integer)       #
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# STAGE 1: BUILDER
# ---------------------------------------------------------------------------- #
FROM debian:bookworm AS builder

LABEL stage="builder"
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install Build Dependencies
RUN echo "--- Installing build dependencies ---" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    g++ gfortran make cmake ninja-build \
    curl ca-certificates file git python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# 2. Copy Source and Scripts
COPY . ./gbtolib-git
COPY conda/docker_build_script.sh .

# 3. Run Build Script (Fixing Windows line endings on the fly)
RUN sed -i 's/\r$//' docker_build_script.sh && \
    chmod +x docker_build_script.sh && \
    ./docker_build_script.sh

# ---------------------------------------------------------------------------- #
# STAGE 2: FINAL IMAGE
# ---------------------------------------------------------------------------- #
FROM debian:bookworm-slim AS final

LABEL description="GBTOlib Standard (OpenMPI 5.0.9 / 4-byte Ints)"

# 1. Install Runtime Dependencies
RUN echo "--- Installing runtime dependencies ---" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    libgomp1 libquadmath0 libgfortran5 \
    openssh-client vim nano less \
    && rm -rf /var/lib/apt/lists/*

# 2. Copy Compiled Software from Builder
COPY --from=builder /opt/gbtolib /opt/gbtolib

# 3. Setup Shell Configuration (Bashrc)
COPY conda/gbtolib.bashrc /opt/gbtolib/gbtolib.bashrc
RUN echo "source /opt/gbtolib/gbtolib.bashrc" >> /etc/bash.bashrc && \
    chmod 644 /opt/gbtolib/gbtolib.bashrc

# 4. Create Entrypoint Script
# This ensures environment variables are loaded for ANY command run in the container
RUN echo '#!/bin/bash' > /opt/gbtolib/entrypoint.sh && \
    echo 'source /opt/gbtolib/gbtolib.bashrc' >> /opt/gbtolib/entrypoint.sh && \
    echo 'exec "$@"' >> /opt/gbtolib/entrypoint.sh && \
    chmod +x /opt/gbtolib/entrypoint.sh

# 5. User Setup
RUN useradd --create-home --shell /bin/bash gbto
USER gbto
WORKDIR /home/gbto

# 6. Environment Variables
# Note: PATH is explicitly set here, but .bashrc aliases provide the switching logic
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    PRTE_ALLOW_RUN_AS_ROOT=1 \
    PRTE_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OPAL_PREFIX="/opt/gbtolib" \
    PATH="/opt/gbtolib/bin.double:/opt/gbtolib/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/gbtolib/lib.double:/opt/gbtolib/lib:${LD_LIBRARY_PATH:-}"

# 7. Entrypoint & Command
ENTRYPOINT ["/opt/gbtolib/entrypoint.sh"]
CMD ["/bin/bash"]