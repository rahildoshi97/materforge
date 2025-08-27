#!/bin/bash -l
#SBATCH --job-name=benchmark_materforge_lumi-c_single_node_256
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=0-01:00:00
#SBATCH --account=project_465001284

unset SLURM_EXPORT_ENV
BUILD_DIR=$1

usage() {
    echo "Usage: $0 build_dir"
    echo "  build_dir  The build directory (required)"
    exit 1
}

if [ $# -eq 0 ]; then
    echo "Missing build_dir argument." >&2
    usage
else
    BUILD_DIR=$(realpath ${BUILD_DIR})
fi

SCRIPT_PATH="$(realpath $0)"
TIMESTAMP=$(date +%s)
date -d @${TIMESTAMP}
echo "Using BUILD_DIR: ${BUILD_DIR}"
echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
cat "${SCRIPT_PATH}"
echo "---"

BINARY="CodegenHeatEquationCPUScaling"
JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR} || exit 1
cd ${JOB_DIR}
cp ${BUILD_DIR}/build_*.log .
set -x

CMD="/project/project_465001284/repos/materforge/apps/cmake-build-lumi-release-cpu/${BINARY}"

# Load environment
module load LUMI/24.03 partition/C PrgEnv-gnu buildtools/24.03 cray-python/3.11.7
source /project/project_465001284/venvs/materforge/bin/activate

echo "=== Single Node Performance Tests - 256 Cells ==="

# Test different thread counts with 256 cells per process
export OMP_NUM_THREADS=1
echo "--- 1 Thread, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

export OMP_NUM_THREADS=2
echo "--- 2 Threads, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

export OMP_NUM_THREADS=4
echo "--- 4 Threads, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

export OMP_NUM_THREADS=8
echo "--- 8 Threads, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

export OMP_NUM_THREADS=16
echo "--- 16 Threads, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

export OMP_NUM_THREADS=32
echo "--- 32 Threads, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

export OMP_NUM_THREADS=64
echo "--- 64 Threads, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

export OMP_NUM_THREADS=128
echo "--- 128 Threads, 1 MPI rank ---"
srun --cpu-freq=2200000-2200000 --ntasks=1 --cpus-per-task=128 ${CMD} weak 256

# Test hybrid MPI+OpenMP with 256 cells
export OMP_NUM_THREADS=16
echo "--- 8 MPI ranks, 16 threads each ---"
srun --cpu-freq=2200000-2200000 --ntasks=8 --cpus-per-task=16 ${CMD} weak 256

export OMP_NUM_THREADS=8
echo "--- 16 MPI ranks, 8 threads each ---"
srun --cpu-freq=2200000-2200000 --ntasks=16 --cpus-per-task=8 ${CMD} weak 256

export OMP_NUM_THREADS=4
echo "--- 32 MPI ranks, 4 threads each ---"
srun --cpu-freq=2200000-2200000 --ntasks=32 --cpus-per-task=4 ${CMD} weak 256

export OMP_NUM_THREADS=2
echo "--- 64 MPI ranks, 2 threads each ---"
srun --cpu-freq=2200000-2200000 --ntasks=64 --cpus-per-task=2 ${CMD} weak 256

export OMP_NUM_THREADS=1
echo "--- 128 MPI ranks, 1 thread each ---"
srun --cpu-freq=2200000-2200000 --ntasks=128 --cpus-per-task=1 ${CMD} weak 256
