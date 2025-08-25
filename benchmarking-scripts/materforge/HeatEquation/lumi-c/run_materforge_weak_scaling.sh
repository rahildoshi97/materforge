#!/bin/bash -l
#SBATCH --job-name=benchmark_materforge_lumi-c_weak_scaling
#SBATCH --output=%x.%j.out
#SBATCH --partition=standard
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=128
#SBATCH --time=0-02:00:00
#SBATCH --account=project_465001284

unset SLURM_EXPORT_ENV

BUILD_DIR=$1
BUILD_DIR=$(realpath ${BUILD_DIR})

BINARY="CodegenHeatEquationCPUScaling"
JOB_DIR="$HOME/lss-rdm/jobs/$SLURM_JOBID/"
mkdir -p ${JOB_DIR} || exit 1
cd ${JOB_DIR}

CMD="${HOME}/rahil/repos/materforge/apps/cmake-build-lumi-release-cpu/${BINARY}"

module load LUMI/24.03 partition/C PrgEnv-gnu buildtools/24.03 cray-python/3.11.7
source ~/rahil/venvs/materforge/bin/activate

echo "=== Weak Scaling Test (Constant Work Per Process) ==="

# Constant problem size per process - true weak scaling
PROBLEM_SIZE=64  # Each process handles 64³ cells
export OMP_NUM_THREADS=1

for n in 1 2 4 8 16 32 64 128 256 512; do
    total_ranks=$((n * 128))
    
    echo "-------------------------------------------------------- Number of nodes: $n (${total_ranks} total ranks, ${PROBLEM_SIZE}³ cells per process) --------------------------------------------------------"
    
    srun --cpu-freq=2200000-2200000 --nodes=$n --ntasks-per-node=128 --cpus-per-task=1 ${CMD} weak ${PROBLEM_SIZE}
done
