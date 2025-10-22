#!/bin/bash -l

# ---
# type: buildscript
# cluster: lumi-g
# application: materforge-couette
# ---

timestamp=$(date +%s)
SCRIPT_PATH="$(realpath $0)"
NAME="couette"

MATERFORGE_REPO="https://i10git.cs.fau.de/rahil.doshi/materforge.git"
MATERFORGE_BRANCH="bm"
MATERFORGE_COMMIT=""

# Configurable paths
SRC_DIR="${MATERFORGE_SRC_DIR:-$(pwd)/materforge_src}"
BUILD_DIR="$(pwd)/build_${NAME}_${timestamp}"
LOGFILE="${BUILD_DIR}/build_${timestamp}.log"
VENV_PATH="${MATERFORGE_VENV:-/project/project_465001284/venvs/materforge}"

# Colors
NC='\033[0m'
ERROR='\033[1;31m'
SUCCESS='\033[1;32m'
INFO='\033[1;36m'
PROGRESS='\033[1;35m'

mkdir -p "${BUILD_DIR}" || { echo -e "${ERROR}Failed to create build directory: ${BUILD_DIR}${NC}"; exit 1; }

{
date -d @${timestamp}
echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
cat "${SCRIPT_PATH}"
echo "---"

# Clone or use existing repository with submodules
if [ ! -d "${SRC_DIR}" ]; then
    echo -e "${PROGRESS}Cloning into ${SRC_DIR}${NC}"
    git clone --recursive -b ${MATERFORGE_BRANCH} ${MATERFORGE_REPO} ${SRC_DIR}
else
    echo -e "${PROGRESS}Using existing repository at ${SRC_DIR}${NC}"
fi

cd "${SRC_DIR}" || exit 1

echo "Current directory: $(pwd)"
echo "Git Info:"
echo "Git url = $(git remote get-url origin 2>/dev/null || echo 'No remote')"
echo "Git branch = $(git branch --show-current 2>/dev/null || echo 'No branch')"
echo "Git commit = $(git rev-parse HEAD 2>/dev/null || echo 'No commit')"
echo

git checkout .
git fetch -a
git checkout $MATERFORGE_BRANCH
# Force reset to remote (discard any local commits)
git reset --hard origin/$MATERFORGE_BRANCH
[[ -n "$MATERFORGE_COMMIT" ]] && git checkout $MATERFORGE_COMMIT
# No git pull needed - already at origin state
# git pull

# Initialize and update submodules
echo -e "${PROGRESS}Updating submodules${NC}"
git submodule update --init --recursive || exit 1

# Load LUMI-G modules
module load LUMI/24.03 partition/G buildtools/24.03 rocm/6.0.3 craype-accel-amd-gfx90a PrgEnv-cray
module load cray-mpich
module list

# Activate virtual environment
if [[ -f "${VENV_PATH}/bin/activate" ]]; then
    source "${VENV_PATH}/bin/activate"
    echo -e "${INFO}Activated virtual environment: ${VENV_PATH}${NC}"
else
    echo -e "${ERROR}Virtual environment not found at ${VENV_PATH}${NC}"
    exit 1
fi

# Build in apps directory
cd ${SRC_DIR}/apps || exit 1

# Clean any existing builds
rm -rf cmake-build-lumi-release-gpu

echo -e "${PROGRESS}Building Couette Flow for GPU using cmake presets${NC}"
set -x

# Use existing cmake presets
cmake --preset lumi-release-gpu || { echo -e "${ERROR}CMake configure failed${NC}"; exit 1; }
cmake --build --preset lumi-release-gpu-build -j 8 || { echo -e "${ERROR}Build failed${NC}"; exit 1; }

set +x

echo "Contents of build directory:"
ls -lah cmake-build-lumi-release-gpu/

# Copy executables and .prm files to build directory
if ls cmake-build-lumi-release-gpu/CouetteFlowGPUScaling 1> /dev/null 2>&1; then
    cp cmake-build-lumi-release-gpu/CouetteFlowGPUScaling ${BUILD_DIR}/
    echo -e "${SUCCESS}CouetteFlowGPUScaling executable copied${NC}"
else
    echo -e "${ERROR}Error: CouetteFlowGPUScaling executable not found${NC}"
    exit 1
fi

# Copy .prm files if they exist
if ls cmake-build-lumi-release-gpu/*.prm 1> /dev/null 2>&1; then
    cp cmake-build-lumi-release-gpu/*.prm ${BUILD_DIR}/
    echo -e "${INFO}.prm files copied${NC}"
fi

echo
echo -e "${SUCCESS}Couette Flow GPU build completed in ${BUILD_DIR}${NC}"
echo "Available files:"
ls -lah ${BUILD_DIR}/

} 2>&1 | tee ${LOGFILE}
