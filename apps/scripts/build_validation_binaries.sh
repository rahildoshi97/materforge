#!/bin/bash
# Build all constant-viscosity validation binaries plus tempdep.
#
# Run on the login node AFTER activating the materforge venv:
#   source ~/.venvs/materforge/bin/activate
#   bash apps/scripts/build_validation_binaries.sh
#
# Each const binary is compiled with the matching CONST_NU value so that the
# baked-in omega = 2/(6*nu+1) is physically correct.  The tempdep binary is
# rebuilt last (restoring USE_MATERFORGE=ON as the active cmake config).
#
# Output binaries (in apps/build/woody-release-cpu/):
#   CouetteFlowScaling_const_0.04  ...  CouetteFlowScaling_const_1.0
#   CouetteFlowScaling_tempdep

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"      # apps/scripts/
APPS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"                       # apps/
BUILD_DIR="${APPS_DIR}/build/woody-release-cpu"
LOG_DIR="${APPS_DIR}/logs/build"
mkdir -p "${LOG_DIR}"

eval "$(MODULESHOME=/apps/modules /usr/bin/modulecmd bash load cmake/3.30.5)"

# ── Constant-viscosity cases ─────────────────────────────────────────────────
NUS=(0.04 0.06 0.08 0.1 0.2 0.4 0.6 0.8 1.0)

for nu in "${NUS[@]}"; do
    echo ""
    echo "══════════════════════════════════════════════"
    echo " Building: const nu=${nu}"
    echo "══════════════════════════════════════════════"

    cmake --preset woody-release-cpu \
          -DUSE_MATERFORGE=OFF -DCONST_NU="${nu}" \
          -S "${APPS_DIR}" -B "${BUILD_DIR}" \
          > "${LOG_DIR}/configure_const_${nu}.log" 2>&1

    make -C "${BUILD_DIR}" -j8 CouetteFlowScaling \
          > "${LOG_DIR}/build_const_${nu}.log" 2>&1

    cp "${BUILD_DIR}/CouetteFlowScaling" \
       "${BUILD_DIR}/CouetteFlowScaling_const_${nu}"

    echo "  -> CouetteFlowScaling_const_${nu}  ($(date '+%H:%M:%S'))"
done

# ── Temperature-dependent case ───────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo " Building: tempdep (MaterForge)"
echo "══════════════════════════════════════════════"

cmake --preset woody-release-cpu \
      -DUSE_MATERFORGE=ON \
      -S "${APPS_DIR}" -B "${BUILD_DIR}" \
      > "${LOG_DIR}/configure_tempdep.log" 2>&1

make -C "${BUILD_DIR}" -j8 CouetteFlowScaling \
      > "${LOG_DIR}/build_tempdep.log" 2>&1

cp "${BUILD_DIR}/CouetteFlowScaling" \
   "${BUILD_DIR}/CouetteFlowScaling_tempdep"

echo "  -> CouetteFlowScaling_tempdep  ($(date '+%H:%M:%S'))"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo " All binaries built:"
ls -lh "${BUILD_DIR}"/CouetteFlowScaling_* | awk '{print "  "$NF, $5}'
echo "══════════════════════════════════════════════"
echo " Next step:"
echo "   sbatch apps/scripts/run_validation_array.sh"
