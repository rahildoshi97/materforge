#!/bin/bash -l
# ---
# type: buildscript
# cluster: lumi-c
# application: apps/materforge
# ---

timestamp=$(date +%s)
SCRIPT_PATH="$(realpath $0)"

NAME="materforge"
MATERFORGE_REPO="https://i10git.cs.fau.de/rahil.doshi/materforge.git"
MATERFORGE_BRANCH="bm"
MATERFORGE_COMMIT=""

SRC_DIR="${HOME}/rahil/repos/materforge"
BUILD_DIR="$(pwd)/build_${NAME}_${timestamp}"
LOGFILE="${BUILD_DIR}/build_${timestamp}.log"

mkdir -p ${BUILD_DIR} || exit 1
{
	date -d @${timestamp}

	echo "--- currently executed script: $(basename ${SCRIPT_PATH})"
	cat "${SCRIPT_PATH}"
	echo "---"

	# Use existing repository instead of cloning
	echo "Using existing repository at ${SRC_DIR}"

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

	# Load LUMI environment
	module load LUMI/24.03 partition/C
	module load PrgEnv-gnu
	module load buildtools/24.03 cray-python/3.11.7

	module list

	# Python virtual environment setup
	source ~/rahil/venvs/materforge/bin/activate

	# Use existing cmake build from apps directory
	cd ${SRC_DIR}/apps
	
	# Clean any existing builds
	rm -rf cmake-build-lumi-*

	echo "Building materforge using existing cmake system"

	# Use existing cmake presets
	cmake --preset lumi-release-cpu
	cmake --build --preset lumi-release-cpu-build

	# Copy built executables to benchmark build directory
	mkdir -p ${BUILD_DIR}
	cp cmake-build-lumi-release-cpu/CodegenHeatEquation* ${BUILD_DIR}/ 2>/dev/null || true

	echo
	echo "materforge build completed in ${BUILD_DIR}"
	echo "Available executables:"
	ls -la ${BUILD_DIR}/CodegenHeatEquation* 2>/dev/null || echo "No executables built"
	
	# Also show where the actual build is
	echo "Actual build location: ${SRC_DIR}/apps/cmake-build-lumi-release-cpu"
	echo "Available executables there:"
	ls -la ${SRC_DIR}/apps/cmake-build-lumi-release-cpu/CodegenHeatEquation* 2>/dev/null || echo "No executables found"

} 2>&1 | tee ${LOGFILE}
