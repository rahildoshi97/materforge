#!/bin/bash -l
#SBATCH --job-name=benchmark_materforge_lumi-c_strong_scaling_large
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

echo "=== Strong Scaling Test - Large (Fixed Total Size: 1024³) ==="

TOTAL_SIZE=1024
export OMP_NUM_THREADS=1

for n in 64 128 256 512; do
    echo "-------------------------------------------------------- Number of nodes: $n --------------------------------------------------------"
    
    srun --cpu-freq=2200000-2200000 --nodes=$n --ntasks-per-node=128 --cpus-per-task=1 ${CMD} strong ${TOTAL_SIZE}
done
