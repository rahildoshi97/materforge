#!/bin/bash

################################################################################
# LUMI Setup Script for MaterForge
# Post-January 2026 Maintenance (ROCm 6.3.4, CCE 19.0.0, LUMI/25.03)
#
# USAGE:
#   source ./lumi.sh cpu      # Load CPU environment (default)
#   source ./lumi.sh gpu      # Load GPU environment
#
# This script must be SOURCED, not executed, to persist environment variables
################################################################################

set -o pipefail

lumi_cpu() {
    # CPU partition with GNU compiler
    module --force purge
    module load LUMI/25.03 partition/C PrgEnv-gnu
    module load cray-python
    
    # Disable GPU MPI for CPU builds
    unset MPICH_GPU_SUPPORT_ENABLED

    # Activate MaterForge virtual environment
    if [[ -f /projappl/project_465002382/venvs/materforge/bin/activate ]]; then
        source /projappl/project_465002382/venvs/materforge/bin/activate
    fi
    
    # Print environment info
    echo "✅ LUMI CPU environment loaded!"
    echo "   Partition: CPU (partition/C)"
    echo "   Python: $(python3 --version)"
    echo "   Compiler: $(cc --version | head -1)"
    echo "   CMake: $(cmake --version | head -1)"
}

lumi_gpu() {
    # GPU partition with Cray compiler and ROCm
    module --force purge
    
    # Clear any lingering HIP/ROCm variables that might conflict
    unset HIP_FLAGS HIP_COMPILE_FLAGS HIP_LINK_FLAGS CMAKE_HIP_FLAGS
    unset CRAY_ROCM_INCLUDE_OPTS HIPCXX_FLAGS
    
    # Load GPU stack
    module load LUMI/25.03 partition/G PrgEnv-cray
    module load rocm/6.3.4        # Explicitly use 6.3.4 (actual version post-maintenance)
    module load cray-mpich/8.1.32 # Explicit version for reproducibility
    module load cray-python
    
    # Enable GPU support for MPI
    export MPICH_GPU_SUPPORT_ENABLED=1
    
    # ROCm 6.3.4 environment paths
    export ROCM_PATH=/opt/rocm-6.3.4
    export HIP_PATH=/opt/rocm-6.3.4
    export DEVICE_LIB_PATH=/opt/rocm-6.3.4/amdgcn/bitcode
    export HIP_DEVICE_LIB_PATH=/opt/rocm-6.3.4/amdgcn/bitcode
    export CMAKE_PREFIX_PATH="/opt/rocm-6.3.4:/opt/rocm-6.3.4/lib/cmake"
    export HIP_PLATFORM=amd
    export CMAKE_HIP_FLAGS=""
    
    # Make MPI headers visible to HIP compiler
    MPI_INCLUDE_PATH="/opt/cray/pe/mpich/8.1.32/ofi/crayclang/19.0/include"
    export CPATH="${MPI_INCLUDE_PATH}:${CPATH}"
    export C_INCLUDE_PATH="${MPI_INCLUDE_PATH}:${C_INCLUDE_PATH}"
    export CPLUS_INCLUDE_PATH="${MPI_INCLUDE_PATH}:${CPLUS_INCLUDE_PATH}"
    
    # Verify GPU support libraries are available
    if [[ -n "$PE_MPICH_GTL_DIR_amd_gfx90a" ]]; then
        export PE_MPICH_GTL_DIR_amd_gfx90a
        export PE_MPICH_GTL_LIBS_amd_gfx90a
    fi
    
    # Activate MaterForge virtual environment
    if [[ -f /projappl/project_465002382/venvs/materforge/bin/activate ]]; then
        source /projappl/project_465002382/venvs/materforge/bin/activate
    fi
    
    # Print environment info
    echo "✅ ROCm 6.3.4 GPU environment loaded!"
    echo "   Partition: GPU (partition/G)"
    echo "   Python: $(python3 --version)"
    echo "   Compiler: $(CC --version | head -1)"
    echo "   ROCM_PATH: $ROCM_PATH"
    echo "   HIP_PATH: $HIP_PATH"
    echo "   MPI: cray-mpich/8.1.32"
    
    # Debug: Show GTL status
    if [[ -n "$PE_MPICH_GTL_DIR_amd_gfx90a" ]]; then
        echo "   GTL: Enabled (PE_MPICH_GTL_DIR_amd_gfx90a=$PE_MPICH_GTL_DIR_amd_gfx90a)"
    fi
}

################################################################################
# Main Entry Point
################################################################################

LUMI_MODE="${1:-cpu}"  # Default to CPU if no argument provided

case "$LUMI_MODE" in
    gpu)
        lumi_gpu
        ;;
    cpu|*)
        lumi_cpu
        ;;
esac

# Ensure sourcing, not execution
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "⚠️  WARNING: This script should be SOURCED, not executed!"
    echo "Usage: source ./lumi.sh [cpu|gpu]"
    exit 1
fi
