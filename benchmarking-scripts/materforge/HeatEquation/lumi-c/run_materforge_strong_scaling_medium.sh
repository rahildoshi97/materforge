#!/bin/bash -l

#SBATCH --job-name=benchmark_materforge_lumi-c_strong_scaling_medium
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=128
#SBATCH --time=0-01:30:00
#SBATCH --account=project_465001284

unset SLURM_EXPORT_ENV

BUILD_DIR=""
PROBLEM_SIZE=768  # Medium fixed total problem size

usage() {
    echo "Usage: $0 build_dir [total_problem_size]"
    echo "  build_dir           The build directory (required)"
    echo "  total_problem_size  Total cells per dimension (optional, default: 768)"
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
    PROBLEM_SIZE=${2:-768}
fi

BUILD_DIR=$(realpath ${BUILD_DIR})
SCRIPT_PATH="$(realpath $0)"

echo "Using BUILD_DIR: ${BUILD_DIR}"
echo "Using TOTAL_PROBLEM_SIZE: ${PROBLEM_SIZE}^3"

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

# Medium strong scaling test: 8 to 64 nodes
for n in 8 16 32 64; do
    echo "========================================="
    echo "Medium Strong Scaling: $n nodes"
    echo "Total problem size: ${PROBLEM_SIZE}^3"
    echo "CPU cores per node: 128, Total cores: $((n*128))"
    echo "========================================="
    
    srun --cpu-freq=2200000-2200000 --nodes=$n ${CMD} strong ${PROBLEM_SIZE}
done
