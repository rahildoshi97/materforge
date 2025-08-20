//======================================================================================================================
//! \file CodegenHeatEquationGPUScaling.cpp
//! \author Rahil Doshi
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "core/Environment.h"
#include "field/AddToStorage.h"
#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/FieldCopy.h"
#include "gpu/HostFieldAllocator.h"
#include "gpu/communication/UniformGPUScheme.h"
#include "gpu/communication/MemcpyPackInfo.h"
#include "stencil/D3Q19.h"
#include "field/vtk/VTKWriter.h"
#include "timeloop/SweepTimeloop.h"
#include "core/timing/TimingPool.h"
#include "gen/HeatEquationKernelWithMaterialGPU.hpp"

namespace walberla {

typedef GhostLayerField<real_t, 1> ScalarField;
typedef gpu::GPUField<real_t> GPUScalarField;

void swapFields(StructuredBlockForest& blocks, BlockDataID uID, BlockDataID uTmpID) {
    for (auto block = blocks.begin(); block != blocks.end(); ++block) {
        GPUScalarField* u = block->getData<GPUScalarField>(uID);
        GPUScalarField* u_tmp = block->getData<GPUScalarField>(uTmpID);
        u->swapDataPointers(u_tmp);
    }
}

void initDirichletBoundaries(const shared_ptr<StructuredBlockForest>& blocks,
                            BlockDataID uId, BlockDataID uTmpId,
                            BlockDataID uCpuId, BlockDataID uTmpCpuId) {
    // Initialize on CPU, then copy to GPU
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        ScalarField* u = block->getData<ScalarField>(uCpuId);
        ScalarField* u_tmp = block->getData<ScalarField>(uTmpCpuId);
        // Apply boundaries on all domain borders
        // North boundary (Y max)
        if (blocks->atDomainYMaxBorder(*block)) {
            CellInterval north = u->xyzSizeWithGhostLayer();
            north.yMin() = north.yMax();
            for (auto cell = north.begin(); cell != north.end(); ++cell) {
                u->get(*cell) = real_c(3800.0);
                u_tmp->get(*cell) = real_c(3800.0);
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

int main(int argc, char** argv) {
    mpi::Environment env(argc, argv);
    gpu::selectDeviceBasedOnMpiRank();

    // Command line argument parsing
    bool weakScaling = false;
    uint_t baseCells = 256;  // Default cells per GPU

    if (argc > 1) {
        std::string scalingType(argv[1]);
        if (scalingType == "weak") weakScaling = true;
    }
    if (argc > 2) {
        baseCells = uint_c(std::stoi(argv[2]));
    }

    const uint_t numGPUs = uint_c(MPIManager::instance()->numProcesses());

    // Calculate domain decomposition for GPUs
    uint_t gpus_x = 1, gpus_y = 1, gpus_z = 1;
    if (numGPUs == 1) {
        gpus_x = 1; gpus_y = 1; gpus_z = 1;
    } else if (numGPUs == 2) {
        gpus_x = 2; gpus_y = 1; gpus_z = 1;
    } else if (numGPUs == 4) {
        gpus_x = 2; gpus_y = 2; gpus_z = 1;
    } else if (numGPUs == 8) {
        gpus_x = 2; gpus_y = 2; gpus_z = 2;
    } else if (numGPUs == 16) {
        gpus_x = 2; gpus_y = 2; gpus_z = 4;
    } else if (numGPUs == 32) {
        gpus_x = 2; gpus_y = 4; gpus_z = 4;
    } else if (numGPUs == 64) {
        gpus_x = 4; gpus_y = 4; gpus_z = 4;
    } else {
        WALBERLA_ABORT("Unsupported number of GPUs: " << numGPUs);
    }

    uint_t xCells, yCells, zCells;
    real_t xSize, ySize, zSize;

    if (weakScaling) {
        // Weak scaling: constant cells per GPU/process
        xCells = baseCells;
        yCells = baseCells;
        zCells = baseCells;
        // Domain size scales with GPU/processes count
        xSize = real_c(gpus_x * 1.0);  // or procs_x for CPU
        ySize = real_c(gpus_y * 1.0);
        zSize = real_c(gpus_z * 1.0);
    } else {
        // Strong scaling: fixed total domain size
        const uint_t totalSize = baseCells;
        xCells = totalSize / gpus_x;
        yCells = totalSize / gpus_y;
        zCells = totalSize / gpus_z;
        // Physical domain remains constant
        xSize = real_c(1.0);
        ySize = real_c(1.0);
        zSize = real_c(1.0);
    }

    const uint_t xBlocks = gpus_x;
    const uint_t yBlocks = gpus_y;
    const uint_t zBlocks = gpus_z;

    WALBERLA_LOG_INFO_ON_ROOT("GPU Scaling Type: " << (weakScaling ? "Weak" : "Strong"));
    WALBERLA_LOG_INFO_ON_ROOT("GPUs: " << numGPUs);
    WALBERLA_LOG_INFO_ON_ROOT("Cells per GPU: " << xCells << "×" << yCells << "×" << zCells);
    WALBERLA_LOG_INFO_ON_ROOT("Memory per GPU: ~" << (3 * xCells * yCells * zCells * 8 / 1024 / 1024) << " MB");

    if (numGPUs != xBlocks * yBlocks * zBlocks) {
        WALBERLA_ABORT("GPU count mismatch!");
    }

    const real_t dx = xSize / real_c(xBlocks * xCells + uint_c(1));
    const real_t dt = real_c(1);
    const uint_t timeSteps = uint_c(2e4);
    const uint_t vtkWriteFrequency = uint_c(400);

    // Block storage setup
    auto aabb = math::AABB(real_c(0.5) * dx, real_c(0.5) * dx, real_c(0.5) * dx,
                          xSize - real_c(0.5) * dx, ySize - real_c(0.5) * dx, zSize - real_c(0.5) * dx);

    shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGrid(
        aabb, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, true, false, false, false);

    // Field initialization with pinned memory
    auto allocator = make_shared<gpu::HostFieldAllocator<real_t>>();
    BlockDataID uFieldCpuId = field::addToStorage<ScalarField>(blocks, "u_cpu", real_c(300.0), field::fzyx, uint_c(1), allocator);
    BlockDataID uTmpFieldCpuId = field::addToStorage<ScalarField>(blocks, "u_tmp_cpu", real_c(300.0), field::fzyx, uint_c(1), allocator);
    BlockDataID alphaFieldCpuId = field::addToStorage<ScalarField>(blocks, "alpha_cpu", real_c(0.0), field::fzyx, uint_c(1), allocator);

    // GPU fields
    BlockDataID uFieldId = gpu::addGPUFieldToStorage<ScalarField>(blocks, uFieldCpuId, "u", true);
    BlockDataID uTmpFieldId = gpu::addGPUFieldToStorage<ScalarField>(blocks, uTmpFieldCpuId, "u_tmp", true);
    BlockDataID alphaFieldId = gpu::addGPUFieldToStorage<ScalarField>(blocks, alphaFieldCpuId, "alpha", true);

    // GPU Communication
    constexpr bool cudaEnabledMPI = false;
    gpu::communication::UniformGPUScheme<stencil::D3Q19> commScheme(blocks, cudaEnabledMPI);
    auto packInfo = make_shared<gpu::communication::MemcpyPackInfo<GPUScalarField>>(uFieldId);
    commScheme.addPackInfo(packInfo);

    // Initialize boundaries
    initDirichletBoundaries(blocks, uFieldId, uTmpFieldId, uFieldCpuId, uTmpFieldCpuId);

    // GPU Timeloop
    SweepTimeloop timeloop(blocks, timeSteps);
    //<< BeforeFunction([&]() {commScheme.getCommunicateFunctor();}, "GPU Communication")
    timeloop.add() << BeforeFunction(commScheme.getCommunicateFunctor(), "GPU Communication")
                   << Sweep(HeatEquationKernelWithMaterialGPU(alphaFieldId, uFieldId, uTmpFieldId, dt, dx), "HeatSolverGPU")
                   << AfterFunction([&]() {swapFields(*blocks, uFieldId, uTmpFieldId);}, "GPU Swap");

    if (vtkWriteFrequency > 0) {
        std::string scalingType = weakScaling ? "weak" : "strong";
        std::string vtkFilename = "vtk_GPU_" + scalingType + "_" + std::to_string(baseCells) + 
                                "cells_" + std::to_string(numGPUs) + "gpu(s)";
        std::string vtkDirectory = "vtk_out_gpu_" + scalingType + "_" + std::to_string(baseCells) + 
                                "cells_" + std::to_string(numGPUs) + "gpu(s)";
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, vtkFilename, vtkWriteFrequency, 0, false, vtkDirectory,
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

    // Benchmark
    const uint_t warmupSteps = uint_c(5);
    for (uint_t i = 0; i < warmupSteps; ++i)
        timeloop.singleStep();

    WcTimingPool timeloopTiming;
    WcTimer simTimer;

    WALBERLA_MPI_WORLD_BARRIER()
    WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
    WALBERLA_LOG_INFO_ON_ROOT("Starting GPU simulation with " << timeSteps << " time steps")

    simTimer.start();
    timeloop.run(timeloopTiming);
    WALBERLA_GPU_CHECK(gpuDeviceSynchronize())
    simTimer.end();

    double simTime = simTimer.max();
    WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); }

    const auto reducedTimeloopTiming = timeloopTiming.getReduced();
    WALBERLA_LOG_RESULT_ON_ROOT("GPU Time loop timing:\n" << *reducedTimeloopTiming)

    uint_t totalCellUpdates = timeSteps * xCells * yCells * zCells;
    uint_t mlups = totalCellUpdates / uint_c(simTime * 1000000.0);
    WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per GPU: " << mlups)

    uint_t globalCells = xCells * yCells * zCells * numGPUs;
    uint_t globalMLUPS = timeSteps * globalCells / uint_c(simTime * 1000000.0);
    WALBERLA_LOG_RESULT_ON_ROOT("Global GPU MLUPS: " << globalMLUPS)

    return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
