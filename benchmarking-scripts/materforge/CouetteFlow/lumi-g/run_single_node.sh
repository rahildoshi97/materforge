#!/bin/bash -l
#SBATCH --job-name=couette_gpu_single_node
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard-g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=0-00:30:00
#SBATCH --account=project_465001980

unset SLURM_EXPORT_ENV

BUILD_DIR=""
PRM_FILE="CouetteScaling.prm"

usage() {
    echo "Usage: $0 build_dir [prm_file]"
    echo "  build_dir  The build directory (required)"
    echo "  prm_file   Configuration file (optional, default: CouetteScaling.prm)"
    exit 1
}

[ $# -eq 0 ] && usage
BUILD_DIR=$(realpath $1)
[ $# -ge 2 ] && PRM_FILE=$2

SCRIPT_PATH="$(realpath $0)"
JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR}
cd ${JOB_DIR}

cp ${BUILD_DIR}/build_*.log . 2>/dev/null || true
cp ${BUILD_DIR}/${PRM_FILE} . 2>/dev/null || true

module load LUMI/24.03 partition/G buildtools/24.03 rocm/6.0.3 craype-accel-amd-gfx90a PrgEnv-cray
module load cray-mpich

set -x

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export MPICH_GPU_SUPPORT_ENABLED=1
export ROCR_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

CMD="${BUILD_DIR}/CouetteFlowGPUScaling"

[ ! -f "${CMD}" ] && { echo "Error: Executable not found"; exit 1; }
[ ! -f "${BUILD_DIR}/${PRM_FILE}" ] && { echo "Error: .prm file not found"; exit 1; }

echo "========================================="
echo "Single Node Couette Flow Benchmark"
echo "Configuration: ${PRM_FILE}"
echo "Total GPUs: 8"
echo "========================================="

srun -n 8 --cpu-bind=${CPU_BIND} ${CMD} ${PRM_FILE}
