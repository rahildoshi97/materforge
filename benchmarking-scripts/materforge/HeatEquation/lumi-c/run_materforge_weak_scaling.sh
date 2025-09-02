#!/bin/bash -l

#SBATCH --job-name=benchmark_materforge_lumi-c_weak_scaling
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=128
#SBATCH --time=0-02:00:00
#SBATCH --account=project_465001284

unset SLURM_EXPORT_ENV

BUILD_DIR=""
PROBLEM_SIZE=256  # Cells per CPU core

usage() {
    echo "Usage: $0 build_dir [cells_per_core]"
    echo "  build_dir     The build directory (required)"
    echo "  cells_per_core Cells per CPU core (optional, default: 256)"
    exit 1
}

if [ $# -eq 0 ]; then
    echo "Missing build_dir argument." >&2
    usage
elif [ $# -gt 2 ]; then
    echo "Too many arguments provided." >&2
    usage
else
    BUILD_DIR=$1
    PROBLEM_SIZE=${2:-256}
fi

BUILD_DIR=$(realpath ${BUILD_DIR})
SCRIPT_PATH="$(realpath $0)"

echo "Using BUILD_DIR: ${BUILD_DIR}"
echo "Using CELLS_PER_CORE: ${PROBLEM_SIZE}^3"

echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
cat "${SCRIPT_PATH}"
echo "---"

BINARY="CodegenHeatEquationCPUScaling"

JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR}
cd ${JOB_DIR}
cp ${BUILD_DIR}/build_*.log . 2>/dev/null

module load LUMI/24.03 partition/C
module load PrgEnv-gnu
module load buildtools/24.03 cray-python/3.11.7

module list

set -x

CMD="${BUILD_DIR}/${BINARY}"

# Weak scaling test: constant work per CPU core
for n in 1 2 4 8 16 32 64 128 256 512; do
    echo "========================================="
    echo "Weak Scaling: $n nodes"
    echo "Cells per core: ${PROBLEM_SIZE}^3"
    echo "CPU cores per node: 128, Total cores: $((n*128))"
    echo "Total cells: $((n*128*PROBLEM_SIZE*PROBLEM_SIZE*PROBLEM_SIZE))"
    echo "========================================="
    
    srun --cpu-freq=2200000-2200000 --nodes=$n ${CMD} weak ${PROBLEM_SIZE}
done
