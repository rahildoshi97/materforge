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

# Updated paths
SRC_DIR="/project/project_465001284/repos/materforge"
BUILD_DIR="$(pwd)/build_${NAME}_${timestamp}"
LOGFILE="${BUILD_DIR}/build_${timestamp}.log"

# Print the values of the variables
echo
echo -e "${PROGRESS}--- Inside the materforge build script ---${NC}"
echo
echo -e "${INFO}BENCHMARK_NAME = $NAME"
echo -e "${INFO}TIME_LABEL = $TIME_LABEL"
echo -e "${INFO}BASE_DIR = $BASE_DIR"
echo -e "${INFO}SRC_DIR = $SRC_DIR"
echo -e "${INFO}BUILD_DIR = $BUILD_DIR"
echo -e "${INFO}BUILD_LOG = $LOGFILE"
echo -e "${INFO}SCRIPT_DIR = $(dirname "$SCRIPT_PATH")${NC}"
echo

mkdir -p "${BUILD_DIR}" || { echo -e "${ERROR}Failed to create build directory: ${BUILD_DIR}${NC}"; exit 1; }

{
    date -d @${timestamp}
    echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
    cat "${SCRIPT_PATH}"
    echo "---"
    
    # Clone or update the materforge source code
    if [ ! -d "${SRC_DIR}" ]; then
        echo "Cloning into ${SRC_DIR}"
        git clone -b ${MATERFORGE_BRANCH} ${MATERFORGE_REPO} ${SRC_DIR}
    else
        echo "Using existing repository at ${SRC_DIR}"
    fi

	echo "Git Info:"
	echo "Git url = $(git remote get-url origin 2>/dev/null || echo 'No remote')"
	echo "Git branch = $(git branch --show-current 2>/dev/null || echo 'No branch')"
	echo "Git commit = $(git rev-parse HEAD 2>/dev/null || echo 'No commit')"
	echo
    
    cd "${SRC_DIR}" || exit 1
    git checkout .
	git fetch -a
	git checkout $MATERFORGE_BRANCH
	[[ -n "$MATERFORGE_COMMIT" ]] && git checkout $MATERFORGE_COMMIT
	git pull

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

    PYTHON_VERSION=$(python --version 2>&1 | awk '{print $2}')
    REQUIRED_VERSION="3.10"
    if [[ "$(printf '%s\n' "$REQUIRED_VERSION" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED_VERSION" ]]; then
        echo -e "${ERROR}Python version $PYTHON_VERSION is less than required version $REQUIRED_VERSION${NC}"
        exit 1
    fi
    echo -e "${INFO}Python version $PYTHON_VERSION is compatible (>= $REQUIRED_VERSION)${NC}"

    # Use existing cmake build from apps directory
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
    echo "Building materforge using cmake presets"
    cmake --preset lumi-release-gpu || { echo "CMake confiCMakegure failed"; exit 1; }
    cmake --build --preset lumi-release-gpu-build || { echo " build failed"; exit 1; }
        
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
