#!/bin/bash -l

#SBATCH --job-name=benchmark_materforge_lumi-c_single_node
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=0-00:30:00
#SBATCH --account=project_465001284

unset SLURM_EXPORT_ENV

# Default values
BUILD_DIR=""
PROBLEM_SIZE=384  # Default cells per CPU core for weak scaling

usage() {
    echo "Usage: $0 build_dir [problem_size]"
    echo "  build_dir    The build directory (required)"
    echo "  problem_size Cells per dimension per core (optional, default: 384)"
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
    PROBLEM_SIZE=${2:-384}
fi

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
BINARY="CodegenHeatEquationCPUScaling"  # Assuming CPU version exists

# Create job directory
JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR}
cd ${JOB_DIR}

# Copy build logs
cp ${BUILD_DIR}/build_*.log . 2>/dev/null || echo "No build logs found"

# Load LUMI-C modules
module load LUMI/24.03 partition/C
module load PrgEnv-gnu
module load buildtools/24.03 cray-python/3.11.7

module list

set -x

CMD="${BUILD_DIR}/${BINARY}"

echo "========================================="
echo "Single Node CPU Weak Scaling Benchmark"
echo "Problem size per core: ${PROBLEM_SIZE}^3"
echo "Total CPU cores: 128"
echo "========================================="

# Execute single node test with fixed total problem size
srun --cpu-freq=2200000-2200000 ${CMD} single ${PROBLEM_SIZE}
