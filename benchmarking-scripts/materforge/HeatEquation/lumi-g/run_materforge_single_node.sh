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

# Get the build directory argument
if [ $# -eq 0 ]; then
    echo "Missing build_dir argument." >&2
    usage
elif [ $# -gt 1 ]; then
    echo "Too many arguments provided." >&2
    usage
else
    BUILD_DIR=$1
fi

# Get build directory argument
echo "BUILD_DIR = $BUILD_DIR"
BUILD_DIR=$(realpath ${BUILD_DIR})
SCRIPT_PATH="$(realpath $0)"
TIMESTAMP=$(date +%s)

echo "Using BUILD_DIR: ${BUILD_DIR}"
echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
cat "${SCRIPT_PATH}"
echo "---"

# Application configuration
APP_PATH=""
BINARY="CodegenHeatEquationGPUScaling"
ARGUMENT=""

# Create job directory
JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR} || usage
cd ${JOB_DIR}

# Copy build logs
cp ${BUILD_DIR}/build_*.log .

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

# Execute with srun
srun --cpu-bind=${CPU_BIND} ${CMD}
