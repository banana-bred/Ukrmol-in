# Conda/Mamba Build Guide

This document explains how to build and use **GBTOlib** - the general bound–continuum integral library used by UKRmol+ - with a **fully reproducible Conda/Mamba toolchain**, and optionally run it in a **Docker container** for cluster environments.

GBTOlib is built automatically in **two precision modes**:
- **Double precision** (default)
- **Quadruple precision** (via `-Dusequadprec`)

---

- [Conda/Mamba Build Guide](#condamamba-build-guide)
  - [What are Conda & Mamba?](#what-are-conda--mamba)
    - [Install Conda or Mamba](#install-conda-or-mamba)
    - [Useful links — Mamba installation & docs](#useful-links--mamba-installation--docs)
  - [Building GBTOlib with Conda/Mamba](#building-gbtolib-with-condamamba)
  - [Environment Overview](#environment-overview)
  - [Results Layout (after build)](#results-layout-after-build)
  - [Tips](#tips)
  - [Running GBTOlib in Docker & Apptainer](#running-gbtolib-in-docker--apptainer)
  - [Notes](#notes)

---

## What are Conda & Mamba?

* **Conda** creates isolated environments with compilers, MPI and scientific libraries; it won’t affect your system install.
* **Mamba** is a faster, drop‑in replacement for conda (same commands). If installed, the scripts prefer it automatically.
* We use the **conda‑forge** channel for portable builds on Linux, macOS and Windows (via WSL).

### Install Conda or Mamba
Choose the correct architecture:
- Most desktops/servers: `x86_64`
- Newer ARM boards/servers (e.g., some cloud instances, Raspberry Pi 4/5 64-bit): `aarch64`

```bash
# Download the latest Miniconda for Linux x86_64
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

**aarch64 example (bash):**
```bash
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh
```

```bash
bash Miniconda3-latest-Linux-x86_64.sh
# or for ARM:
# bash Miniconda3-latest-Linux-aarch64.sh
```

Follow the prompts:
- Press **Enter** to review the license, then **yes** to accept.
- Choose the install location (default: `~/miniconda3` is fine).
- Allow the installer to run `conda init` (recommended).

#### Activate Conda in your shell
If you allowed `conda init`, start a new shell or source your rc file:

```bash
# For bash:
source ~/.bashrc
# For zsh:
source ~/.zshrc
```

Update conda to the latest version:
```bash
conda update -n base -c defaults conda
```
#### Bonus: Install Mamba (faster solver)

Mamba is a near drop-in replacement for `conda` that dramatically speeds up environment solves.
```bash
# Using conda-forge (recommended)
conda install -n base -c conda-forge mamba
```

You can now use `mamba` instead of `conda` for faster installs, e.g.,
```bash
# Create env with mamba
mamba create -n fastenv python=3.11

# Install packages
mamba install numpy pandas scikit-learn ...
```

---


#### Useful links — Mamba installation & docs

- Mamba project (overview, source):  
  https://github.com/mamba-org/mamba

- Mamba user guide (installation & usage):  
  https://mamba.readthedocs.io/en/latest/user_guide/mamba.html

- Micromamba (single static binary, no base env; great for containers/CI):  
  https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html

- Mambaforge (Conda + Mamba + conda-forge as default):  
  https://github.com/conda-forge/miniforge#miniforge3-and-mambaforge

---


## Building GBTOlib with Conda/Mamba

The Conda/Mamba setup provides a portable, reproducible build of GBTOlib using the **GNU + OpenMPI + OpenBLAS** toolchain.  
All compilers, MPI wrappers, and libraries come from `conda-forge`, ensuring consistent ILP64 compatibility.

### Quick start

```bash
bash conda/build_with_conda.sh
```

**What happens**

1. Creates or updates the `gbtolib_env` environment from [`environment.yml`](./environment.yml).  
2. Activates it and configures compilers (`mpicc`, `mpifort`, `mpicxx`).  
3. Applies Fortran flags for **ILP64 integers**, **Fortran 2018**, **OpenMP**, and **optimization (-O3)**.  
4. Builds and installs GBTOlib twice:
   - **Double precision** → `install/double/`
   - **Quad precision** → `install/quad/`
5. Prints versions of compilers and linked libraries for reproducibility.

---

## Environment Overview

| Package | Purpose |
|----------|----------|
| `compilers` | GNU toolchain (gcc, g++, gfortran) |
| `openmpi` | MPI compiler wrappers and runtime |
| `openblas`, `libblas`, `liblapack` | Linear algebra libraries |
| `cmake`, `make` | Build system |
| `coreutils` | Utility tools (e.g., `env`, `timeout`) |

To create or activate the environment manually:

```bash
conda env create -f conda/environment.yml
conda activate gbtolib_env
```

---

## Results Layout (after build)

```
conda/
 ├─ environment.yml         # Conda environment definition
 ├─ build_with_conda.sh     # Main build script (double + quad precision)
 └─ README.md               # This file

install/
 ├─ double/                 # Installed double-precision binaries and libs
 └─ quad/                   # Installed quad-precision binaries and libs
```

---


## Tips

- Always activate the environment before building or running.  
- If `mpirun` isn’t found, the environment isn’t active (`which mpirun` should point inside Conda).

---
## Running GBTOlib in Docker & Apptainer

For HPC environments or local development where strict reproducibility is required, we provide a **Ubuntu 22.04-based Dockerfile** that builds the entire toolchain (**OpenMPI 5**, **OpenBLAS**, **zlib**, **hwloc**) **from source**.

### 1. Build the Docker Image

You must build the image from the root of the repository so that the build context includes the source files, but point to the Dockerfile inside the `conda` directory.

```bash
docker build -f conda/Dockerfile -t gbtolib:latest .
```

### 1b. Alternative: Pull Pre-built Image

Instead of building the image from source, you can pull the latest version directly from the GitLab Container Registry.

**1. Log in to the Registry**
(Required because the repository is private)
```bash
docker login registry.gitlab.com
# Enter your GitLab username and Personal Access Token when prompted
```
**2. Pull the Image**
```bash
docker pull registry.gitlab.com/uk-amor/ukrmol/gbtolib:latest
```

### 2. Run Locally (Docker)

To run interactively on your laptop or workstation, mount your current directory (where your input files are) to `/work` inside the container.

```bash
docker run --rm -it -v "$(pwd):/work"  gbtolib:latest
```

Once inside:
```bash
# Check status
gbtolib_help

# Run calculation 
mpirun -np 2 --oversubscribe scatci_integrals target.scatci_integrals.inp
```

### 3. Run on HPC (Apptainer)

On a cluster, you typically cannot run Docker directly. You must convert the image to Apptainer (Singularity) format.

The setup below is an example for a particular cluster running the SLURM scheduler. The arguments for the `salloc` command and paths for `scratch` will be different on your cluster.
Please consult your cluster admin if you're not sure how to modify this example.

**Step A: Export Image (On Local Machine)**
```bash
# Save and Compress
docker save gbtolib:latest | gzip > gbtolib.tar.gz
```
*Transfer `gbtolib.tar.gz` to your cluster.*

**Step B: Interactive Job on HPC**
Inside an interactive allocation (`salloc`). It is recommended to set Apptainer temporary directories to scratch to avoid filling up your home quota (Step 3 below).

```bash
# 1. Allocate resources (Adjust partition/mem/cpu as needed)
salloc -p edu --mem 10G --cpus-per-task 2 -t 2:0:0

# 2. Unzip archive
gunzip gbtolib.tar.gz

# 3. Configure Apptainer temp directories
export APPTAINER_CACHEDIR=/scratch/tmp/$USER/apptainer_tmp
export APPTAINER_TMPDIR=/scratch/tmp/$USER/apptainer_tmp
mkdir -p $APPTAINER_TMPDIR

# 4. Build the SIF image
apptainer build gbtolib.sif docker-archive://gbtolib.tar

# 5. Run the container
# --bind mounts your cluster work directory into the container
apptainer run --bind /work/$USER gbtolib.sif
```

### 4. Executing GBTOlib Inside the Container (HPC)

Once inside the container, you must prepare the environment to prevent the internal OpenMPI from conflicting with the cluster's Slurm scheduler.

**CRITICAL STEP: Detach from Slurm**
Because you are inside a Slurm job (`salloc`), OpenMPI 5 detects the `SLURM_` environment variables and tries to use the cluster's `srun` command, which does not exist inside the container. You must unset these variables before running:

```bash
# 1. Hide Slurm from the container's internal MPI
unset SLURM_JOBID SLURM_JOB_ID SLURM_NODELIST

# 2. Run the calculation
# -np should match your salloc --cpus-per-task
# --oversubscribe ensures it runs even if it detects a core mismatch
mpirun -np 2 --oversubscribe scatci_integrals target.scatci_integrals.inp
```

### Helper Commands

The container includes aliases to help manage precision modes:

| Command | Description |
| :--- | :--- |
| `gbtolib_help` | Shows installed versions (OpenMPI, OpenBLAS), paths, and current precision. |
| `use_double` | Switches environment to **Double Precision** (Default). |
| `use_quad` | Switches environment to **Quadruple Precision**. |

---
