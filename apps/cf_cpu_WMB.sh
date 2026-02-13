#!/bin/bash -l
#SBATCH --job-name=cf_profile   # Job name
#SBATCH --partition=standard    # partition name
#SBATCH --nodes=16               # Total number of nodes 
#SBATCH --ntasks=2048           # Total number of mpi tasks
#SBATCH --time=0-00:30:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_465002382   # Project for billing
#SBATCH --output=cf_cpu_perftools.o%j # Name of stdout output file
#SBATCH --error=cf_cpu_perftools.e%j  # Name of stderr error file

#module --force purge
#module load LUMI/25.03 partition/C PrgEnv-cray cray-python
#module load perftools-base
#module load perftools-lite

#echo "=== MODULE LIST ==="
#module list

# ./build/lumi-release-cpu/CouetteFlowScaling \
# ./build/lumi-release-cpu/CouetteFlowScaling.prm \

srun --cpu-freq=2200000 \
  /projappl/project_465002382/repos/materforge/apps/build/lumi-release-cpu/CouetteFlowScaling \
  /projappl/project_465002382/repos/materforge/apps/build/lumi-release-cpu/CouetteFlowScaling.prm \
  -DomainSetup.blocks=\<16,16,8\> \
  -DomainSetup.cellsPerBlock=\<128,64,128\> \
  &> out/cf_cpu_tempdep_inter_strong_16N_2048P_b16x16x8_cpb128x64x128_WMB_3.txt
