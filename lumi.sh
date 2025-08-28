#!/bin/bash

# LUMI setup for materforge project

# CPU Setup - Default execution
lumi_cpu() {
    module purge
    module --force purge
    module load LUMI/24.03 partition/C
    module load PrgEnv-gnu
    module load buildtools/24.03 cray-python/3.11.7
    
    # Activate virtual environment
    if [[ -f /project/project_465001284/venvs/materforge/bin/activate ]]; then
        source /project/project_465001284/venvs/materforge/bin/activate
    fi
    
    echo "✅ LUMI CPU environment loaded!"
    echo "Partition: CPU (partition/C)"
    echo "Compiler: $(cc --version | head -1)"
}

# GPU Setup (ROCm)
lumi_gpu() {
    module purge  
    module --force purge
    
    # Clear ALL potentially problematic environment variables
    unset HIP_FLAGS
    unset HIPCXX_FLAGS
    unset CMAKE_HIP_FLAGS
    unset CRAY_ROCM_INCLUDE_OPTS
    unset HIP_COMPILE_FLAGS
    unset HIP_LINK_FLAGS
    
    module load LUMI/24.03 partition/G PrgEnv-cray buildtools/24.03 craype-accel-amd-gfx90a rocm/6.0.3
    module load cray-mpich
    
    # CRITICAL: Enable GPU support for MPI GTL
    export MPICH_GPU_SUPPORT_ENABLED=1

    # Force override the problematic variable AFTER module loading
    export CRAY_ROCM_INCLUDE_OPTS="-I/opt/rocm-6.0.3/include -I/opt/rocm-6.0.3/include/rocprofiler -I/opt/rocm-6.0.3/include/roctracer -I/opt/rocm-6.0.3/include/hip"
    unset CRAY_ROCM_INCLUDE_OPTS
    
    # Set ROCm environment variables
    export ROCM_PATH=/opt/rocm-6.0.3
    export HIP_PATH=/opt/rocm-6.0.3  
    export DEVICE_LIB_PATH=/opt/rocm-6.0.3/amdgcn/bitcode
    export HIP_DEVICE_LIB_PATH=/opt/rocm-6.0.3/amdgcn/bitcode
    export CMAKE_PREFIX_PATH="/opt/rocm-6.0.3:/opt/rocm-6.0.3/lib/cmake"
    export HIP_PLATFORM=amd
    
    # Additional safeguards - clear any cached CMake HIP flags
    export CMAKE_HIP_FLAGS=""
    
    # CRITICAL: Set environment variables to make MPI headers visible to HIP compiler
    MPI_INCLUDE_PATH="/opt/cray/pe/mpich/8.1.29/ofi/crayclang/17.0/include"
    export CPATH="${MPI_INCLUDE_PATH}:$CPATH"
    export C_INCLUDE_PATH="${MPI_INCLUDE_PATH}:$C_INCLUDE_PATH"  
    export CPLUS_INCLUDE_PATH="${MPI_INCLUDE_PATH}:$CPLUS_INCLUDE_PATH"
    
    # Debug: Check if GTL libraries are available
    echo "Checking GTL environment variables:"
    echo "PE_MPICH_GTL_DIR_amd_gfx90a: $PE_MPICH_GTL_DIR_amd_gfx90a"
    echo "PE_MPICH_GTL_LIBS_amd_gfx90a: $PE_MPICH_GTL_LIBS_amd_gfx90a"

    # Activate virtual environment
    if [[ -f /project/project_465001284/venvs/materforge/bin/activate ]]; then
        source /project/project_465001284/venvs/materforge/bin/activate
    fi
    
    echo "✅ ROCm 6.0.3 environment configured for materforge"
    echo "Partition: GPU (partition/G)"
    echo "Compiler: $(CC --version | head -1)"
    echo "MPI: $(module list cray-mpich 2>&1 | grep cray-mpich || echo 'Not loaded')"
    echo "MPI Include Path: $MPI_INCLUDE_PATH"
}

# Default behavior: setup CPU environment
if [[ $1 == "cpu" ]]; then
    lumi_cpu
elif [[ $1 == "gpu" ]]; then
    lumi_gpu
else
    # Default to CPU setup
    lumi_cpu
fi
