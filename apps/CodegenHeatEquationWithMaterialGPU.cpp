//======================================================================================================================
//
//! \file CodegenHeatEquationWithMaterialGPU.cpp
//! \author Rahil Doshi
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"
#include "pde/boundary/Neumann.h"
#include "stencil/D3Q19.h"
#include "timeloop/SweepTimeloop.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/timing/TimingPool.h"

// GPU-specific includes
#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/FieldCopy.h"
#include "gpu/GPUWrapper.h"
#include "gpu/HostFieldAllocator.h"
//#include "gpu/ParallelStreams.h"
#include "gpu/communication/MemcpyPackInfo.h"
#include "gpu/communication/UniformGPUScheme.h"

#include "gen/HeatEquationKernelWithMaterialGPU.hpp"

namespace walberla {

///////////////////////
/// Typedef Aliases ///
///////////////////////

typedef GhostLayerField<real_t, 1> ScalarField;
typedef gpu::GPUField<real_t> GPUScalarField;

void swapFields(StructuredBlockForest& blocks, BlockDataID uID, BlockDataID uTmpID) {
    for (auto block = blocks.begin(); block != blocks.end(); ++block) {
        GPUScalarField* u = block->getData<GPUScalarField>(uID);
        GPUScalarField* u_tmp = block->getData<GPUScalarField>(uTmpID);
        u->swapDataPointers(u_tmp);
    }
}

void initDirichletBoundaryNorth(const shared_ptr<StructuredBlockForest>& blocks,
                               BlockDataID uId, BlockDataID uTmpId,
                               BlockDataID uCpuId, BlockDataID uTmpCpuId) {
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        // Initialize on CPU first
        if (blocks->atDomainYMaxBorder(*block)) {
            ScalarField* u = block->getData<ScalarField>(uCpuId);
            ScalarField* u_tmp = block->getData<ScalarField>(uTmpCpuId);
            CellInterval xyz = u->xyzSizeWithGhostLayer();
            xyz.yMin() = xyz.yMax();
            for (auto cell = xyz.begin(); cell != xyz.end(); ++cell) {
                const Vector3<real_t> p = blocks->getBlockLocalCellCenter(*block, *cell);
                // Set constant temperature boundary condition (same as CPU)
                real_t v = real_c(3800.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
    }
    // Copy to GPU
    //gpu::fieldCpy<GPUScalarField, ScalarField>(blocks, uId, uCpuId);
    //gpu::fieldCpy<GPUScalarField, ScalarField>(blocks, uTmpId, uTmpCpuId);
}

void initDirichletBoundariesAllSides(const shared_ptr<StructuredBlockForest>& blocks,
                                    BlockDataID uId, BlockDataID uTmpId,
                                    BlockDataID uCpuId, BlockDataID uTmpCpuId) {
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        ScalarField* u = block->getData<ScalarField>(uCpuId);
        ScalarField* u_tmp = block->getData<ScalarField>(uTmpCpuId);
        // North boundary (Y max)
        if (blocks->atDomainYMaxBorder(*block)) {
            CellInterval north = u->xyzSizeWithGhostLayer();
            north.yMin() = north.yMax();
            for (auto cell = north.begin(); cell != north.end(); ++cell) {
                real_t v = real_c(3800.0); // Hot boundary
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
        // South boundary (Y min)
        if (blocks->atDomainYMinBorder(*block)) {
            CellInterval south = u->xyzSizeWithGhostLayer();
            south.yMax() = south.yMin();
            for (auto cell = south.begin(); cell != south.end(); ++cell) {
                real_t v = real_c(300.0); // Cold boundary
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
        // East boundary (X max)
        if (blocks->atDomainXMaxBorder(*block)) {
            CellInterval east = u->xyzSizeWithGhostLayer();
            east.xMin() = east.xMax();
            for (auto cell = east.begin(); cell != east.end(); ++cell) {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
        // West boundary (X min)
        if (blocks->atDomainXMinBorder(*block)) {
            CellInterval west = u->xyzSizeWithGhostLayer();
            west.xMax() = west.xMin();
            for (auto cell = west.begin(); cell != west.end(); ++cell) {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
        // Top boundary (Z max)
        if (blocks->atDomainZMaxBorder(*block)) {
            CellInterval top = u->xyzSizeWithGhostLayer();
            top.zMin() = top.zMax();
            for (auto cell = top.begin(); cell != top.end(); ++cell) {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
        // Bottom boundary (Z min)
        if (blocks->atDomainZMinBorder(*block)) {
            CellInterval bottom = u->xyzSizeWithGhostLayer();
            bottom.zMax() = bottom.zMin();
            for (auto cell = bottom.begin(); cell != bottom.end(); ++cell) {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
    }
    // Copy to GPU
    gpu::fieldCpy<GPUScalarField, ScalarField>(blocks, uId, uCpuId);
    gpu::fieldCpy<GPUScalarField, ScalarField>(blocks, uTmpId, uTmpCpuId);
}

int main(int argc, char** argv)
{
    mpi::Environment env(argc, argv);

    // GPU device selection
    gpu::selectDeviceBasedOnMpiRank();

    /////////////////////////////
    /// SIMULATION PARAMETERS ///
    /////////////////////////////

    // Use same parameters as CPU version
    constexpr uint_t x = 128;
    // cells in X,Y,Z direction per block
    const uint_t xCells = uint_c(x);
    const uint_t yCells = uint_c(x);
    const uint_t zCells = uint_c(x);
    // Physical domain size in X,Y,Z
    const real_t xSize = real_c(1.0);
    const real_t ySize = real_c(1.0);  // 2
    const real_t zSize = real_c(1.0);
    // Number of blocks in X,Y,Z direction
    const uint_t xBlocks = uint_c(1);
    const uint_t yBlocks = uint_c(1);  // 2
    const uint_t zBlocks = uint_c(1);

    const uint_t processes = uint_c(MPIManager::instance()->numProcesses());
    if (processes != xBlocks * yBlocks * zBlocks)
    { WALBERLA_ABORT("The number of processes must be equal to the number of blocks!"); }

    const real_t dx = xSize / real_c(xBlocks * xCells + uint_c(1));
    const real_t dy = ySize / real_c(yBlocks * yCells + uint_c(1));
    const real_t dz = zSize / real_c(zBlocks * zCells + uint_c(1));

    // WALBERLA_CHECK_FLOAT_EQUAL(dx, dy);

    const real_t dt = real_c(1);
    uint_t timeSteps = uint_c(2e4);
    uint_t vtkWriteFrequency = uint_c(0);

    ///////////////////////////
    /// BLOCK STORAGE SETUP ///
    ///////////////////////////

    auto aabb = math::AABB(real_c(0.5) * dx, real_c(0.5) * dy, real_c(0.5) * dz,
                          xSize - real_c(0.5) * dx, ySize - real_c(0.5) * dy, zSize - real_c(0.5) * dz);

    shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGrid(
        aabb, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, true, false, false, false);

    //////////////
    /// FIELDS ///
    //////////////

    // CPU fields for initialization and VTK output
    auto allocator = make_shared<gpu::HostFieldAllocator<real_t>>(); // use pinned memory for faster transfers
    BlockDataID uFieldCpuId = field::addToStorage<ScalarField>(blocks, "u_cpu", real_c(300.0), field::fzyx, uint_c(1), allocator);
    BlockDataID uTmpFieldCpuId = field::addToStorage<ScalarField>(blocks, "u_tmp_cpu", real_c(300.0), field::fzyx, uint_c(1), allocator);
    BlockDataID alphaFieldCpuId = field::addToStorage<ScalarField>(blocks, "thermal_diffusivity_cpu", real_c(0.0), field::fzyx, uint_c(1), allocator);

    // GPU fields for computation
    BlockDataID uFieldId = gpu::addGPUFieldToStorage<ScalarField>(blocks, uFieldCpuId, "u", true);
    BlockDataID uTmpFieldId = gpu::addGPUFieldToStorage<ScalarField>(blocks, uTmpFieldCpuId, "u_tmp", true);
    BlockDataID alphaFieldId = gpu::addGPUFieldToStorage<ScalarField>(blocks, alphaFieldCpuId, "thermal_diffusivity", true);

    /////////////////////
    /// COMMUNICATION ///
    /////////////////////

    constexpr bool cudaEnabledMPI = false; // Set to true if CUDA-aware MPI is available
    gpu::communication::UniformGPUScheme<stencil::D3Q19> commScheme(blocks, cudaEnabledMPI);
    auto packInfo = make_shared<gpu::communication::MemcpyPackInfo<GPUScalarField>>(uFieldId);
    commScheme.addPackInfo(packInfo);

    ////////////////////////////
    /// FIELD INITIALIZATION ///
    ////////////////////////////

    // Initialize alpha field properly
    //initAlphaField(blocks, alphaFieldId, alphaFieldCpuId);

    //////////////////////////
    /// DIRICHLET BOUNDARY ///
    //////////////////////////

    // Initialize boundary conditions
    //initDirichletBoundaryNorth(blocks, uFieldId, uTmpFieldId, uFieldCpuId, uTmpFieldCpuId);  // Neumann
    initDirichletBoundariesAllSides(blocks, uFieldId, uTmpFieldId, uFieldCpuId, uTmpFieldCpuId);  // Dirichlet

    ////////////////////////
    /// NEUMANN BOUNDARY ///
    ////////////////////////

    //pde::NeumannDomainBoundary<ScalarField> neumann(*blocks, uFieldCpuId);  // Neumann
    //neumann.excludeBoundary(stencil::N);  // Neumann

    ////////////////
    /// TIMELOOP ///
    ////////////////

    SweepTimeloop timeloop(blocks, timeSteps);

    timeloop.add() << BeforeFunction([&]() {
        // Apply Neumann boundaries on CPU
        //neumann();  // Neumann
        //initDirichletBoundariesAllSides(blocks, uFieldId, uTmpFieldId, uFieldCpuId, uTmpFieldCpuId);  // Dirichlet

        // Copy updated boundary values to GPU
        //gpu::fieldCpy<GPUScalarField, ScalarField>(blocks, uFieldId, uFieldCpuId);

        // Also copy boundary conditions to u_tmp field
        //gpu::fieldCpy<GPUScalarField, ScalarField>(blocks, uTmpFieldId, uTmpFieldCpuId);

        // GPU communication
        commScheme.getCommunicateFunctor();
    }, "GPU Communication")
    << Sweep(HeatEquationKernelWithMaterialGPU(alphaFieldId, uFieldId, uTmpFieldId, dt, dx), "HeatEquationKernelWithMaterialGPU")
    << AfterFunction([&]() {
        swapFields(*blocks, uFieldId, uTmpFieldId);

        // Copy results back to CPU for boundary condition updates
        //gpu::fieldCpy<ScalarField, GPUScalarField>(blocks, uFieldCpuId, uFieldId);
        //gpu::fieldCpy<ScalarField, GPUScalarField>(blocks, uTmpFieldCpuId, uTmpFieldId);

        // Reapply Dirichlet boundaries on CPU
        // initDirichletBoundaryNorth(blocks, uFieldId, uTmpFieldId, uFieldCpuId, uTmpFieldCpuId);
    }, "GPU Swap");

    if (vtkWriteFrequency > 0)
    {
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtkGPU", vtkWriteFrequency, 0, false, "vtk_out_gpu",
                                                       "simulation_step", false, true, true, false, 0);

        auto tempWriter = make_shared<field::VTKWriter<ScalarField>>(uFieldCpuId, "temperature");
        vtkOutput->addCellDataWriter(tempWriter);

        auto alphaWriter = make_shared<field::VTKWriter<ScalarField>>(alphaFieldCpuId, "thermal_diffusivity");
        vtkOutput->addCellDataWriter(alphaWriter);

        vtkOutput->addBeforeFunction([&]() {
            // Copy GPU data to CPU for VTK output
            // Only copy when actually writing VTK (every 400 steps)
            gpu::fieldCpy<ScalarField, GPUScalarField>(blocks, uFieldCpuId, uFieldId);
            gpu::fieldCpy<ScalarField, GPUScalarField>(blocks, alphaFieldCpuId, alphaFieldId);
        });

        timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output GPU");
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// BENCHMARK ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const uint_t warmupSteps = uint_c(5);
    const uint_t outerIterations = uint_c(1);

    for (uint_t i = 0; i < warmupSteps; ++i)
        timeloop.singleStep();

    auto remainingTimeLoggerFrequency = real_c(-1.0); // in seconds
    if (remainingTimeLoggerFrequency > 0)
    {
        auto logger = timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps() * outerIterations,
                                                 remainingTimeLoggerFrequency);
        timeloop.addFuncAfterTimeStep(logger, "remaining time logger");
    }

    for (uint_t outerIteration = 0; outerIteration < outerIterations; ++outerIteration)
    {
        WALBERLA_GPU_CHECK(gpuPeekAtLastError())
        timeloop.setCurrentTimeStepToZero();
        WcTimingPool timeloopTiming;
        WcTimer simTimer;

        WALBERLA_MPI_WORLD_BARRIER()
        WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
        WALBERLA_GPU_CHECK(gpuPeekAtLastError())

        WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timeSteps << " time steps")
        simTimer.start();
        timeloop.run(timeloopTiming);
        WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
        simTimer.end();
        WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")

        double simTime = simTimer.max();
        WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); }

        const auto reducedTimeloopTiming = timeloopTiming.getReduced();
        WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)

        uint_t mlups = timeSteps * xCells * yCells * zCells / (uint_c(simTime * 1000000.0));
        WALBERLA_LOG_RESULT_ON_ROOT("mlups:\t" << mlups)
    }

    return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
