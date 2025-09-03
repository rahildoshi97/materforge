#!/bin/bash -l
#SBATCH --job-name=benchmark_materforge_lumi-g_single_node
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard-g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=0-00:30:00
#SBATCH --account=project_465001980

# Enable export of environment from this script to srun
unset SLURM_EXPORT_ENV

# Default values
BUILD_DIR=""
PROBLEM_SIZE=768  # Default cells per GPU for weak scaling

usage() {
    echo "Usage: $0 build_dir [problem_size]"
    echo "  build_dir    The build directory (required)"
    echo "  problem_size Cells per GPU (optional, default: 768)"
    exit 1
}

# Argument parsing
if [ $# -eq 0 ]; then
    echo "Missing build_dir argument." >&2
    usage
elif [ $# -gt 2 ]; then
    echo "Too many arguments provided." >&2
    usage
else
    BUILD_DIR=$1
    PROBLEM_SIZE=${2:-768}  # Default to 768 cells per GPU
fi

# Get build directory argument
BUILD_DIR=$(realpath ${BUILD_DIR})
SCRIPT_PATH="$(realpath $0)"
TIMESTAMP=$(date +%s)

echo "Using BUILD_DIR: ${BUILD_DIR}"
echo "Using PROBLEM_SIZE: ${PROBLEM_SIZE}"

# Print script for reproducibility
echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
cat "${SCRIPT_PATH}"
echo "---"

# Application configuration
BINARY="CodegenHeatEquationGPUScaling"

# Create job directory
JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR} || usage
cd ${JOB_DIR}

# Copy build logs
cp ${BUILD_DIR}/build_*.log . 2>/dev/null || echo "No build logs found"

# Load necessary modules
#module load LUMI/24.03 \
#    partition/G \
#    cpeCray/24.03 \
#    craype-x86-trento \
#    craype-accel-amd-gfx90a \
#    rocm/6.0.3 \
#    cray-mpich/8.1.29 \
#    craype-network-ofi \
#    cray-python/3.11.7 \
#    buildtools/24.03
module load LUMI/24.03 partition/G buildtools/24.03 rocm/6.0.3 craype-accel-amd-gfx90a PrgEnv-cray
module load cray-mpich

module list

set -x

# LUMI-G specific configuration
CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export MPICH_GPU_SUPPORT_ENABLED=1

CMD="${BUILD_DIR}/${BINARY}"

# Verify executable exists
if [ ! -f "${CMD}" ]; then
    echo "Error: Executable ${CMD} not found!"
    ls -la ${BUILD_DIR}/
    exit 1
fi

# Single node benchmark with weak scaling (optimal for single node)
echo "========================================="
echo "Single Node Weak Scaling Benchmark"
echo "Problem size per GPU: ${PROBLEM_SIZE}^3"
echo "Total GPUs: 8"
echo "========================================="

# Execute single node test with fixed total problem size
srun --cpu-bind=${CPU_BIND} ${CMD} single ${PROBLEM_SIZE}
