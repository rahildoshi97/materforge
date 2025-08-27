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
    module load LUMI/24.03 partition/G PrgEnv-cray buildtools/24.03 craype-accel-amd-gfx90a rocm/6.0.3

    # Set ROCm environment variables
    export ROCM_PATH=/opt/rocm-6.0.3
    export HIP_PATH=/opt/rocm-6.0.3  
    export DEVICE_LIB_PATH=/opt/rocm-6.0.3/amdgcn/bitcode
    export HIP_DEVICE_LIB_PATH=/opt/rocm-6.0.3/amdgcn/bitcode
    export CMAKE_PREFIX_PATH="/opt/rocm-6.0.3:/opt/rocm-6.0.3/lib/cmake"

    # Activate virtual environment
    if [[ -f /project/project_465001284/venvs/materforge/bin/activate ]]; then
        source /project/project_465001284/venvs/materforge/bin/activate
    fi

    echo "✅ ROCm 6.0.3 environment configured for materforge"
    echo "Partition: GPU (partition/G)"
    echo "Compiler: $(CC --version | head -1)"
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
