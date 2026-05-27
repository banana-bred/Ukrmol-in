# Conda/Mamba and Docker Build Guide

This document explains how to build and run **UKRmol+** using two robust, reproducible methods:
1.  **Conda/Mamba:** Builds ScaLAPACK and ELPA from source against a modern Conda toolchain. This is a reliable method for direct use on a Linux workstation or cluster.
2.  **Docker (Recommended for Portability):** Creates a completely self-contained, portable container by building the **entire dependency stack from source**, including the compilers themselves. This is the most robust and reproducible method, guaranteeing that the software can run anywhere Docker is installed.

This document also covers the optional Psi4 analysis environment.
## Table of Contents

- [Conda/Mamba and Docker Build Guide](#condamamba-and-docker-build-guide)
  - [What are Conda & Mamba?](#what-are-conda--mamba)
    - [Install Conda or Mamba](#install-conda-or-mamba)
    - [Useful links — Mamba installation & docs](#useful-links--mamba-installation--docs)
  - [Getting Started: First-Time Setup](#getting-started-first-time-setup)
  - [Build Option: ScaLAPACK + ELPA](#build-option-scalapack--elpa)
  - [Optional: Psi4 analysis environment](#optional-psi4-analysis-environment)
  - [Results layout (after build)](#results-layout-after-build)
  - [Tips](#tips)
  - [User guide — working with the Conda environment](#user-guide--working-with-the-conda-environment)
    - [Activate the environment](#activate-the-environment)
    - [Running jobs](#running-jobs)
    - [Running the test suite](#running-the-test-suite)
    - [Common checks & tips](#common-checks--tips)
  - [**Docker & Apptainer Build Guide (Recommended)**](#docker--apptainer-build-guide-recommended)
    - [Directory Setup & Layout](#1-directory-setup--layout)
    - [Build the Docker Image](#2-build-the-docker-image)
    - [Run Locally (Docker)](#3-run-locally-docker)
    - [Run on HPC (Apptainer)](#4-run-on-hpc-apptainer)
    - [Helper Commands](#helper-commands)
---

## What are Conda & Mamba?

*   **Conda** creates isolated computing environments with compilers, MPI, and scientific libraries; it won’t affect your system installation.
*   **Mamba** is a faster, drop‑in replacement for conda (same commands). If installed, our scripts prefer it automatically.
*   We use the **conda‑forge** channel for portable builds on Linux, macOS, and Windows (via WSL).

### Install Conda or Mamba
Choose the correct architecture for your system.

```bash
# Download the latest Miniconda for Linux x86_64
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the prompts: accept the license, choose the install location (default is fine), and **allow the installer to run `conda init`**.

#### Activate Conda in your shell
If you allowed `conda init`, start a new shell or source your rc file:

```bash
# For bash:
source ~/.bashrc
```

#### Bonus: Install Mamba (faster solver)
Mamba dramatically speeds up environment creation and package installation.
```bash
conda install -n base -c conda-forge mamba
```
---

#### Useful links — Mamba installation & docs

- Mamba user guide (installation & usage):  
  https://mamba.readthedocs.io/en/latest/user_guide/mamba.html
- Mambaforge (Conda + Mamba + conda-forge as default):  
  https://github.com/conda-forge/miniforge#miniforge3-and-mambaforge

---

## Getting Started: First-Time Setup

If this is your first time building UKRmol+, follow these initial steps to prepare your workspace.

### 1. Create a Project Directory

First, create a single parent folder for your project.

```bash
mkdir my_ukrmol_project
cd my_ukrmol_project
```

### 2. Download the Required Repositories

The build system requires two repositories, `UKRmol-in` and `UKRmol-out`, to be located side-by-side.  

```bash
git clone https://gitlab.com/Uk-amor/UKRMol/UKRmol-in.git
git clone https://gitlab.com/Uk-amor/UKRMol/UKRmol-out.git
```

### 3. Verify the Directory Structure
Your directory structure should look like this:

```
my_ukrmol_project/
 ├─ UKRmol-in/      # UKRmol-in repository
 └─ UKRmol-out/     # UKRmol-out repository
```

### 4. Navigate into the Source Directory

**This is a critical step.** All build and run commands must be executed from *within* the `UKRmol-in` directory.

```bash
cd UKRmol-in
```

Your terminal is now in the correct location (`my_ukrmol_project/UKRmol-in`). You are ready to proceed.

---

## Build Option: ScaLAPACK + ELPA

This is the standard build for production runs, large eigenproblems, and robust distributed linear algebra. It uses a Conda environment to provide a modern toolchain, then builds ScaLAPACK and ELPA from source. The script installs UKRmol+ in both **double** and **quad** precision.

**Build UKRmol+**

```bash
bash conda/build_with_conda_sca_elpa.sh
```

**Run an example**
```bash
mpirun -n 4 ./install/double_SCALAPACK/bin/mpi-scatci input.in
```

---

## Optional: Psi4 analysis environment

If your workflow includes **Psi4** to generate molecular orbitals, the build script for the ScaLAPACK/ELPA environment automatically installs Psi4.

To use it, activate the environment:
```bash
conda activate ukrmol-scalapack-env
python -c "import psi4; print('Psi4', psi4.__version__)"
```
---

## Results layout (after build)

```
my_ukrmol_project/
 ├─ UKRmol-in/
 │  ├─ build_double_SCALAPACK/
 │  └─ build_quad_SCALAPACK/
 │
 │  ├─ conda/
 │  │  ├─ build_with_conda_sca_elpa.sh
 │  │  └─ environment.ukrmol.scalapack.yml
 │
 │  └─ install/
 │     ├─ double_SCALAPACK/       # ScaLAPACK+ELPA build (double)
 │     └─ quad_SCALAPACK/         # ScaLAPACK+ELPA build (quad)
 │
 └─ UKRmol-out/
```
---

## Tips
- Always `conda activate` the environment before building or running.
- If `mpirun` isn’t found, you’re likely outside the environment (`which mpirun` should point inside it).

---

## User guide — working with the Conda environment

This section assumes the ScaLAPACK+ELPA environment and installs already exist.

### Activate the environment
```bash
conda activate ukrmol-scalapack-env
```

### Running jobs
From your working directory containing your input files:

```bash
# Run a double precision calculation
mpirun -n 4 /path/to/my_ukrmol_project/UKRmol-in/install/double_SCALAPACK/bin/mpi-scatci input.in
```

*   Use `install/quad_SCALAPACK/bin/...` for quadruple precision.

### Running the test suite
1.  Navigate to the relevant build directory, e.g., `cd build_double_SCALAPACK`.
2.  Follow the instructions from `UKRmol-in/README.md` to run the test suite, e.g., `ctest -R serial`.

### Common checks & tips
*   **Is the environment active?** The command `which mpirun` should point to a path inside your conda environment.
*   **Performance:** The default ELPA build uses generic kernels for portability. For maximum performance on a specific x86 CPU, you can enable optimized kernels (AVX2, etc.) in the `build_with_conda_sca_elpa.sh` script.

---
## Docker & Apptainer Build Guide (Recommended)

For HPC environments or local development where strict reproducibility is required, we provide a **Ubuntu 22.04-based Dockerfile**. This build compiles the entire toolchain (**OpenMPI 5**, **OpenBLAS**, **ScaLAPACK**, **ELPA**, **Psi4**...) **from source**.

### 1. Directory Setup & Layout

**Crucial Step:** The Docker build process requires **both** the `ukrmol-in` and `ukrmol-out` repositories to be present in the build context. You must run the build command from the **Parent Directory** containing both repositories.

**Required Structure:**
```text
my_ukrmol_project/           <-- RUN DOCKER BUILD HERE
 ├── ukrmol-in/              <-- Clone of UKRmol-in
 │   └── conda/Dockerfile
 └── ukrmol-out/             <-- Clone of UKRmol-out
```

### 2. Build the Docker Image

Navigate to the parent project directory and run the build command, pointing to the Dockerfile inside the `ukrmol-in/conda` folder.

```bash
cd my_ukrmol_project

docker build -f ukrmol-in/conda/Dockerfile -t ukrmol-plus:latest .
```
*Note: This process compiles everything from source and may take several hours.*

### 3. Run Locally (Docker)

To run interactively on your laptop or workstation, mount your current directory (where your input files are) to `/work` inside the container.

```bash
cd /path/to/my/calculations
docker run --rm -it -v "$(pwd):/work" ukrmol-plus:latest
```

Once inside:
```bash
# Check status and installed versions
ukrmol_help

# Run calculation (Double Precision is default)
mpirun -n 4 --oversubscribe mpi-scatci input.in
```

### 4. Run on HPC (Apptainer)

On a cluster, you typically cannot run Docker directly. You must convert the image to Apptainer (Singularity) format.

The setup below is an example for a particular cluster running the SLURM scheduler. The arguments for the `salloc` command and paths for `scratch` will be different on your cluster.
Please consult your cluster admin if you're not sure how to modify this example.

**Step A: Export Image (On Local Machine)**
```bash
# Save and Compress
docker save ukrmol-plus:latest | gzip > ukrmol-plus.tar.gz
```
*Transfer `ukrmol-plus.tar.gz` to your cluster.*

**Step B: Interactive Job on HPC**
Inside an interactive allocation (`salloc`). It is recommended to set Apptainer temporary directories to scratch to avoid filling up your home quota (Step 3 below).

```bash
# 1. Allocate resources (Adjust partition/mem/cpu as needed)
salloc -p edu --mem 10G --cpus-per-task 4 -t 2:0:0

# 2. Unzip archive
gunzip ukrmol-plus.tar.gz

# 3. Configure Apptainer temp directories
export APPTAINER_CACHEDIR=/scratch/tmp/$USER/apptainer_tmp
export APPTAINER_TMPDIR=/scratch/tmp/$USER/apptainer_tmp
mkdir -p $APPTAINER_TMPDIR

# 4. Build the SIF image
apptainer build ukrmol-plus.sif docker-archive://ukrmol-plus.tar

# 5. Run the container
# --bind mounts your cluster work directory into the container
apptainer run --bind /work/$USER ukrmol-plus.sif
```

**Step C: Execution (Known Issue)**

Work in Progress: Running mpirun inside the container on HPC clusters (which use external schedulers like Slurm) is currently experiencing compatibility issues.

Status: A fix for HPC mpirun execution is currently being developed.

Current State: Jobs may fail to launch or communicate correctly across nodes until this patch is released.

### Helper Commands

The container includes aliases to help manage precision modes:

| Command | Description |
| :--- | :--- |
| `ukrmol_help` | Shows installed versions (OpenMPI, ScaLAPACK, Psi4, etc.), paths, and status. |
| `use_double` | Switches environment to **Double Precision** (Default). |
| `use_quad` | Switches environment to **Quadruple Precision**. |

---
