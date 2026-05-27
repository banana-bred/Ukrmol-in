# ============================================================================ #
# UKRmol+ Custom Environment Configuration (Apptainer Friendly)                #
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# 1. Core Logic (Runs in BOTH Interactive and Non-Interactive modes)
# ---------------------------------------------------------------------------- #

# Load Version Information
if [ -f /opt/ukrmolp/versions.env ]; then
    source /opt/ukrmolp/versions.env
fi

# Aliases for switching precision
alias use_double='source /opt/ukrmolp/bin.double/env.bash'
alias use_quad='source /opt/ukrmolp/bin.quad/env.bash'

# Help Function
ukrmol_help() {
  # Detect which bin directory is currently in PATH
  FULL_PATH=$(command -v mpi-scatci 2>/dev/null)
  
  if [ -z "$FULL_PATH" ]; then
    CURRENT_BIN="None"
    PREC="None"
  else
    CURRENT_BIN=$(dirname "$FULL_PATH")
    if [[ "$CURRENT_BIN" == *"quad"* ]]; then 
        PREC="QUADRUPLE"
    else 
        PREC="DOUBLE"
    fi
  fi
  
  echo -e '\n\033[1;34m╔══════════════════════════════════════════════════════════════╗\033[0m'
  echo -e '\033[1;34m║\033[0m              \033[1;32mUKRmol+ Environment Status\033[0m                      \033[1;34m║\033[0m'
  echo -e '\033[1;34m╚══════════════════════════════════════════════════════════════╝\033[0m'
  
  echo -e '\033[1;34m:: Configuration\033[0m'
  echo -e "  Precision:  \033[1;33m$PREC\033[0m"
  echo -n "  OS:         "; cat /etc/os-release | grep PRETTY_NAME | cut -d'"' -f2
  echo -e "  Compiler:   GCC $UK_GCC_VER"
  
  echo -e '\n\033[1;34m:: Linked Libraries (Source Build)\033[0m'
  echo -e "  OpenMPI:    $UK_MPI_VER"
  echo -e "  OpenBLAS:   $UK_BLAS_VER"
  echo -e "  ScaLAPACK:  $UK_SCALAPACK_VER"
  echo -e "  ELPA:       $UK_ELPA_VER"
  echo -e "  ARPACK-NG:  $UK_ARPACK_VER"
  echo -e "  PETSc:      $UK_PETSC_VER"
  echo -e "  SLEPc:      $UK_SLEPC_VER"

  echo -e '\n\033[1;34m:: Quantum Chemistry & Python\033[0m'
  echo -e "  Psi4:       $UK_PSI4_VER"
  echo -e "  Python:     $UK_PYTHON_VER"
  echo -e "  NumPy:      $UK_NUMPY_VER"

  echo -e '\n\033[1;34m:: Environment Commands\033[0m'
  echo -e "  \033[1;36muse_double\033[0m    Switch to \033[1;33mDouble Precision\033[0m"
  echo -e "  \033[1;36muse_quad\033[0m      Switch to \033[1;35mQuadruple Precision\033[0m"
  echo -e "  \033[1;36mukrmol_help\033[0m   Show this help message"

  echo -e '\n\033[1;34m:: Installed Programs ('$CURRENT_BIN')\033[0m'
  if [ -d "$CURRENT_BIN" ]; then
      ls --color=auto -w 90 "$CURRENT_BIN"
  else
      echo "No binaries currently in PATH (Run use_double or use_quad)"
  fi
  echo -e ''
}
export -f ukrmol_help

# ---------------------------------------------------------------------------- #
# 2. Interactive Visuals (Run ONLY if user is typing in a terminal)
# ---------------------------------------------------------------------------- #

# Check if shell is interactive (contains 'i')
if [[ $- == *i* ]]; then

    # Determine Container Type for Prompt
    if [ -f /.dockerenv ]; then
        CT_TYPE="Docker"
    else
        CT_TYPE="Apptainer"
    fi

    # Setup Prompt
    export PS1="\[\033[01;32m\]$CT_TYPE\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ "

    # Startup Banner
    echo -e '\033[1;34m╔════════════════════════════════════════════════════════════════╗\033[0m'
    echo -e '\033[1;34m║\033[0m\033[1;32m             Welcome to the UKRmol+ Container Suite             \033[0m\033[1;34m║\033[0m'
    echo -e '\033[1;34m╚════════════════════════════════════════════════════════════════╝\033[0m'
    echo -e ''
    echo -e '\033[1;32mCurrent Environment:\033[0m  \033[1;33mDOUBLE PRECISION\033[0m \033[1;32m(Default)\033[0m'
    echo -e ''
    echo -e '\033[1;32mEnvironment Commands:\033[0m'
    echo -e '   \033[1;36muse_quad\033[0m     \033[0;90m::\033[0m \033[1;32mSwitch to\033[0m \033[1;31mQuadruple Precision\033[0m'
    echo -e '   \033[1;36muse_double\033[0m   \033[0;90m::\033[0m \033[1;32mSwitch to\033[0m \033[1;33mDouble Precision\033[0m'
    echo -e '   \033[1;36mukrmol_help\033[0m  \033[0;90m::\033[0m \033[1;32mShow detailed status and libraries\033[0m'
    echo -e ''

fi