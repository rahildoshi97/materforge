#!/bin/bash
# ROCm setup for materforge project

# Load required modules
module purge  
module --force purge
module load LUMI/24.03 partition/G PrgEnv-cray buildtools/24.03 craype-accel-amd-gfx90a rocm/6.0.3

# Set ROCm environment variables
export ROCM_PATH=/opt/rocm-6.0.3
export HIP_PATH=/opt/rocm-6.0.3  
export DEVICE_LIB_PATH=/opt/rocm-6.0.3/amdgcn/bitcode
export HIP_DEVICE_LIB_PATH=/opt/rocm-6.0.3/amdgcn/bitcode
export CMAKE_PREFIX_PATH="/opt/rocm-6.0.3:/opt/rocm-6.0.3/lib/cmake"

echo "ROCm 6.0.3 environment configured for materforge"
