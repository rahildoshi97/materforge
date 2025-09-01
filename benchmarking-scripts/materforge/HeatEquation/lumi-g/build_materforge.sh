#!/bin/bash -l
# ---
# type: buildscript
# cluster: lumi-g
# application: materforge
# ---

timestamp=$(date +%s)

SCRIPT_PATH="$(realpath $0)"
NAME="materforge"

MATERFORGE_REPO="https://i10git.cs.fau.de/rahil.doshi/materforge.git"
MATERFORGE_BRANCH="bm"
MATERFORGE_COMMIT=""

SRC_DIR="${MATERFORGE_SRC_DIR:-$(pwd)/materforge_src}"
BUILD_DIR="$(pwd)/build_${NAME}_${timestamp}"
LOGFILE="${BUILD_DIR}/build_${timestamp}.log"
VENV_PATH="${MATERFORGE_VENV:-/project/project_465001284/venvs/materforge}"

mkdir -p "${BUILD_DIR}" || { echo -e "${ERROR}Failed to create build directory: ${BUILD_DIR}${NC}"; exit 1; }

{
    date -d @${timestamp}
    echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
    cat "${SCRIPT_PATH}"
    echo "---"
    
    # Clone or use existing repository
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
    [[ -n "$MATERFORGE_COMMIT" ]] && git checkout $MATERFORGE_COMMIT
    git pull

    # Initialize and update submodules
    echo -e "${PROGRESS}Updating submodules${NC}"
    git submodule update --init --recursive || exit 1

    # Load LUMI-G modules
    module load LUMI/24.03 partition/G buildtools/24.03 rocm/6.0.3 craype-accel-amd-gfx90a PrgEnv-cray
    module load cray-mpich
    module list
    
    # python setup
    VENV_PATH="/project/project_465001284/venvs/materforge"
    if [[ -f "${VENV_PATH}/bin/activate" ]]; then
        source "${VENV_PATH}/bin/activate"
        echo -e "${INFO}Activated virtual environment: ${VENV_PATH}${NC}"
    else
        echo -e "${ERROR}Virtual environment not found at ${VENV_PATH}${NC}"
        exit 1
    fi

    # Activate virtual environment
    if [ -f "${VENV_PATH}/bin/activate" ]; then
        source "${VENV_PATH}/bin/activate"
    else
        echo "Warning: Virtual environment not found at ${VENV_PATH}"
    fi

    # Build in apps directory
    cd ${SRC_DIR}/apps || exit 1

    # Clean any existing builds
    rm -rf cmake-build-lumi-release-gpu

    echo -e "${PROGRESS}Building materforge using cmake system${NC}"
    set -x
	
    # Configure with proper flags for GTL library linking
    # CMAKE_CXX_FLAGS="--offload-arch=gfx90a -I$MPICH_DIR/include -L$MPICH_DIR/lib -lmpi $PE_MPICH_GTL_DIR_amd_gfx90a $PE_MPICH_GTL_LIBS_amd_gfx90a"

    # Use existing cmake presets but with corrected configuration
    # cmake --preset lumi-release-gpu \
        # -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS"

    # Use existing cmake presets
    cmake --preset lumi-release-gpu || { echo -e "${ERROR}CMake configure failed${NC}"; exit 1; }
    cmake --build --preset lumi-release-gpu-build || { echo -e "${ERROR}Build failed${NC}"; exit 1; }        

    set +x

    # Copy executables to build directory
    mkdir -p ${BUILD_DIR}/apps
    if ls cmake-build-lumi-release-gpu/CodegenHeatEquation* 1> /dev/null 2>&1; then
        cp cmake-build-lumi-release-gpu/CodegenHeatEquation* ${BUILD_DIR}/
        echo "Executables copied successfully"
    else
        echo "Warning: No CodegenHeatEquation executables found"
    fi

    echo
    echo -e "${SUCCESS}materforge build completed in ${BUILD_DIR}${NC}"
    echo "Available executables:"
    ls -la ${BUILD_DIR}/CodegenHeatEquation* 2>/dev/null || echo "No executables built"
    
} 2>&1 | tee ${LOGFILE}
