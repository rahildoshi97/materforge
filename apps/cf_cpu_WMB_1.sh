#!/bin/bash -l
#SBATCH --job-name=cf_profile   # Job name
#SBATCH --partition=standard    # partition name
#SBATCH --nodes=64               # Total number of nodes 
#SBATCH --ntasks=8192           # Total number of mpi tasks
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
  /projappl/project_465002382/repos/materforge/apps/build/lumi-release-cpu-1/CouetteFlowScaling \
  /projappl/project_465002382/repos/materforge/apps/build/lumi-release-cpu-1/CouetteFlowScaling.prm \
  -DomainSetup.blocks=\<32,16,16\> \
  -DomainSetup.cellsPerBlock=\<64,64,64\> \
  &> out/cf_cpu_tempdep_inter_strong_64N_8192P_b32x16x16_cpb64x64x64_WMB_3.txt
