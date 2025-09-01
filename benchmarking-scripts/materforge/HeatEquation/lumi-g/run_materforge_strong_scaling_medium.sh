#!/bin/bash -l

#SBATCH --job-name=benchmark_materforge_lumi-g_strong_scaling_medium
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard-g
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=0-01:30:00
#SBATCH --account=project_465001134

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

BINARY="CodegenHeatEquationGPUScaling"

JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR}
cd ${JOB_DIR}
cp ${BUILD_DIR}/build_*.log . 2>/dev/null

module load LUMI/24.03 partition/G buildtools/24.03 rocm/6.0.3 craype-accel-amd-gfx90a PrgEnv-cray
module load cray-mpich
module list

set -x

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export MPICH_GPU_SUPPORT_ENABLED=1

CMD="${BUILD_DIR}/${BINARY}"

# Medium strong scaling test: 8 to 64 nodes
for n in 8 16 32 64; do
    echo "========================================="
    echo "Medium Strong Scaling: $n nodes"
    echo "Total problem size: ${PROBLEM_SIZE}^3"
    echo "GPUs per node: 8, Total GPUs: $((n*8))"
    echo "========================================="
    
    srun --cpu-bind=${CPU_BIND} --nodes=$n ${CMD} strong ${PROBLEM_SIZE}
done
