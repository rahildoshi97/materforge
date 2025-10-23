#!/bin/bash -l
#SBATCH --job-name=couette_gpu_weak_scaling
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard-g
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=0-01:00:00
#SBATCH --account=project_465001980

unset SLURM_EXPORT_ENV

BUILD_DIR=""
PROBLEM_SIZE=128

usage() {
    echo "Usage: $0 build_dir [problem_size]"
    exit 1
}

[ $# -eq 0 ] && usage
BUILD_DIR=$(realpath $1)
[ $# -ge 2 ] && PROBLEM_SIZE=$2

JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR}
cd ${JOB_DIR}

cp ${BUILD_DIR}/build_*.log . 2>/dev/null || true
cp ${BUILD_DIR}/*.prm . 2>/dev/null || true

module load LUMI/24.03 partition/G buildtools/24.03 rocm/6.0.3 craype-accel-amd-gfx90a PrgEnv-cray
module load cray-mpich

set -x

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export MPICH_GPU_SUPPORT_ENABLED=1
export ROCR_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

CMD="${BUILD_DIR}/CouetteFlowGPUScaling"
PRM_FILE="${BUILD_DIR}/CouetteScaling.prm"

[ ! -f "${CMD}" ] && { echo "Error: Executable not found"; exit 1; }

echo "========================================="
echo "Weak Scaling Couette Flow Benchmark"
echo "Cells per GPU: ${PROBLEM_SIZE}^3"
echo "Total nodes: ${SLURM_NNODES}"
echo "Total GPUs: $((SLURM_NNODES * 8))"
echo "========================================="

srun -n $((SLURM_NNODES * 8)) --cpu-bind=${CPU_BIND} \
    ${CMD} ${PRM_FILE} weak ${PROBLEM_SIZE}
