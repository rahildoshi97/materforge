#!/bin/bash -l

#SBATCH --job-name=benchmark_materforge_lumi-g_weak_scaling
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard-g
#SBATCH --nodes=128
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=0-01:00:00
#SBATCH --account=project_465001980

unset SLURM_EXPORT_ENV

BUILD_DIR=""
PROBLEM_SIZE=512  # Cells per GPU

usage() {
    echo "Usage: $0 build_dir [cells_per_gpu]"
    echo "  build_dir     The build directory (required)"
    echo "  cells_per_gpu Cells per GPU (optional, default: 512)"
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
    PROBLEM_SIZE=${2:-512}
fi

BUILD_DIR=$(realpath ${BUILD_DIR})
SCRIPT_PATH="$(realpath $0)"

echo "Using BUILD_DIR: ${BUILD_DIR}"
echo "Using CELLS_PER_GPU: ${PROBLEM_SIZE}^3"

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

# Weak scaling test: constant work per GPU
for n in 1 2 4 8 16 32 64 128; do
    echo "========================================="
    echo "Weak Scaling: $n nodes"
    echo "Cells per GPU: ${PROBLEM_SIZE}^3"
    echo "GPUs per node: 8, Total GPUs: $((n*8))"
    echo "Total cells: $((n*8*PROBLEM_SIZE*PROBLEM_SIZE*PROBLEM_SIZE))"
    echo "========================================="
    
    srun --cpu-bind=${CPU_BIND} --nodes=$n ${CMD} weak ${PROBLEM_SIZE}
done
