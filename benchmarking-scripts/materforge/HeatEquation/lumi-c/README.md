# Materforge Heat Equation Benchmarking Scripts for LUMI-C

## Overview
These scripts are designed to benchmark the materforge heat equation solver on LUMI-C CPU nodes.

## Scripts
- `build_materforge.sh` - Build script for the materforge application
- `run_materforge_single_node.sh` - Single node performance analysis
- `run_materforge_weak_scaling.sh` - Weak scaling study (1-512 nodes)
- `run_materforge_strong_scaling_small.sh` - Strong scaling (1-8 nodes, 256³ problem)
- `run_materforge_strong_scaling_medium.sh` - Strong scaling (8-64 nodes, 512³ problem)
- `run_materforge_strong_scaling_large.sh` - Strong scaling (64-512 nodes, 1024³ problem)

## Usage

### 1. Build the application
```bash
./build_materforge.sh
```

### 2. Submit benchmarking jobs
```bash
# Get build directory from build output
BUILD_DIR=/path/to/build_materforge_timestamp

# Submit jobs
sbatch run_materforge_single_node.sh $BUILD_DIR
sbatch run_materforge_weak_scaling.sh $BUILD_DIR
sbatch run_materforge_strong_scaling_small.sh $BUILD_DIR
sbatch run_materforge_strong_scaling_medium.sh $BUILD_DIR
sbatch run_materforge_strong_scaling_large.sh $BUILD_DIR
```

### 3. Monitor jobs
```bash
squeue -u $USER
```

## Expected Results
- **Single Node**: Optimal MPI/OpenMP configuration identification
- **Weak Scaling**: Constant time per timestep as nodes increase
- **Strong Scaling**: Speedup measurement for fixed problem sizes

## Notes
- Uses account `project_465001284` (available CPU hours)
- All jobs use `partition=standard` for multi-node runs
- CPU frequency pinned to 2.2GHz for consistent results
- Results saved in `$HOME/lss-rdm/jobs/$SLURM_JOBID/`
