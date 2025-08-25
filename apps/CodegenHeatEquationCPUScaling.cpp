//======================================================================================================================
//! \file CodegenHeatEquationCPUScaling.cpp
//! \author Rahil Doshi
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "stencil/D3Q19.h"
#include "field/vtk/VTKWriter.h"
#include "timeloop/SweepTimeloop.h"
#include "core/timing/TimingPool.h"
#include "gen/HeatEquationKernelWithMaterialCPU.hpp"

namespace walberla {

typedef GhostLayerField<real_t, 1> ScalarField;

void swapFields(StructuredBlockForest& blocks, BlockDataID uID, BlockDataID uTmpID) {
    for (auto block = blocks.begin(); block != blocks.end(); ++block) {
        ScalarField* u = block->getData<ScalarField>(uID);
        ScalarField* u_tmp = block->getData<ScalarField>(uTmpID);
        u->swapDataPointers(u_tmp);
    }
}


// Automatic process decomposition for any number of processes
std::tuple<uint_t, uint_t, uint_t> calculateProcessDecomposition(uint_t numProcesses) {
    // Find best factorization for 3D decomposition
    uint_t procs_x = 1, procs_y = 1, procs_z = 1;
    // Start with cube root as base
    uint_t base = uint_c(std::round(std::cbrt(real_c(numProcesses))));
    // Find factors close to cube root
    for (uint_t z = base; z >= 1; --z) {
        if (numProcesses % z == 0) {
            uint_t remaining = numProcesses / z;
            uint_t base_xy = uint_c(std::round(std::sqrt(real_c(remaining))));     
            for (uint_t y = base_xy; y >= 1; --y) {
                if (remaining % y == 0) {
                    uint_t x = remaining / y;
                    procs_x = x; procs_y = y; procs_z = z;
                    goto found;
                }
            }
        }
    }
    found:   
    return std::make_tuple(procs_x, procs_y, procs_z);
}

void initDirichletBoundaries(const shared_ptr<StructuredBlockForest>& blocks,
                            BlockDataID uId, BlockDataID uTmpId) {
    // Simplified boundary initialization - hot north (3800K), cold elsewhere (300K)
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        ScalarField* u = block->getData<ScalarField>(uId);
        ScalarField* u_tmp = block->getData<ScalarField>(uTmpId);
        // Apply boundaries on all domain borders
        // North boundary (Y max)
        if (blocks->atDomainYMaxBorder(*block)) {
            CellInterval north = u->xyzSizeWithGhostLayer();
            north.yMin() = north.yMax();
            for (auto cell = north.begin(); cell != north.end(); ++cell) {
                u->get(*cell) = real_c(3800.0);  // Hot boundary
                u_tmp->get(*cell) = real_c(3800.0);
            }
        }
        // South boundary (Y min)
        if (blocks->atDomainYMinBorder(*block)) {
            CellInterval south = u->xyzSizeWithGhostLayer();
            south.yMax() = south.yMin();
            for (auto cell = south.begin(); cell != south.end(); ++cell) {
                u->get(*cell) = real_c(300.0); // Cold boundary
                u_tmp->get(*cell) = real_c(300.0);
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
}

int main(int argc, char** argv) {
    mpi::Environment env(argc, argv);

    // Command line argument parsing
    bool weakScaling = false;
    uint_t problemSize = 256;  // Default cells per block/rank

    if (argc > 1) {
        std::string scalingType(argv[1]);
        if (scalingType == "weak") weakScaling = true;
        else if (scalingType == "strong") weakScaling = false;
        else {WALBERLA_LOG_WARNING_ON_ROOT("Unknown scaling type: " << scalingType << ". Use 'weak' or 'strong'. Defaulting to strong.");}
    }
    if (argc > 2) {
        problemSize = uint_c(std::stoi(argv[2]));
    }

    const uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());

    // Automatic process decomposition
    uint_t procs_x, procs_y, procs_z;
    std::tie(procs_x, procs_y, procs_z) = calculateProcessDecomposition(numProcesses);

    uint_t xCells, yCells, zCells;
    real_t xSize, ySize, zSize;

    if (weakScaling) {
        // Weak scaling: constant cells per process
        xCells = problemSize;
        yCells = problemSize;
        zCells = problemSize;
        // Domain size scales with process count
        xSize = real_c(procs_x * 1.0);  // or gpus_x for GPU
        ySize = real_c(procs_y * 1.0);
        zSize = real_c(procs_z * 1.0);
        WALBERLA_LOG_INFO_ON_ROOT("Weak Scaling: " << problemSize << "^3 cells per process");
    } else {
        // Strong scaling: fixed total domain size
        const uint_t totalCellsPerDim = problemSize;
        xCells = totalCellsPerDim / procs_x;
        yCells = totalCellsPerDim / procs_y;
        zCells = totalCellsPerDim / procs_z;

        // Ensure cells are evenly divisible
        if (totalCellsPerDim % procs_x != 0 || totalCellsPerDim % procs_y != 0 || totalCellsPerDim % procs_z != 0) {
            WALBERLA_ABORT("Problem size " << problemSize << " not evenly divisible by process grid " 
                          << procs_x << "x" << procs_y << "x" << procs_z);
        }

        // Physical domain remains constant
        xSize = real_c(1.0);
        ySize = real_c(1.0);
        zSize = real_c(1.0);

        WALBERLA_LOG_INFO_ON_ROOT("Strong Scaling: " << problemSize << "^3 total cells");
    }

    const uint_t xBlocks = procs_x;
    const uint_t yBlocks = procs_y;
    const uint_t zBlocks = procs_z;

    WALBERLA_LOG_INFO_ON_ROOT("=== Scaling Configuration ===");
    WALBERLA_LOG_INFO_ON_ROOT("Scaling Type: " << (weakScaling ? "Weak" : "Strong"));
    WALBERLA_LOG_INFO_ON_ROOT("Processes: " << numProcesses);
    WALBERLA_LOG_INFO_ON_ROOT("Process Grid: " << procs_x << "x" << procs_y << "x" << procs_z);
    WALBERLA_LOG_INFO_ON_ROOT("Cells per process: " << xCells << "x" << yCells << "x" << zCells);
    WALBERLA_LOG_INFO_ON_ROOT("Total cells: " << (xCells*xBlocks) << "x" << (yCells*yBlocks) << "x" << (zCells*zBlocks));
    WALBERLA_LOG_INFO_ON_ROOT("Total cells count: " << (xCells*yCells*zCells*numProcesses));

    // Verification
    if (numProcesses != xBlocks * yBlocks * zBlocks) {
        WALBERLA_ABORT("Process decomposition failed: " << numProcesses << " != " << xBlocks*yBlocks*zBlocks);
    }

    const real_t dx = xSize / real_c(xBlocks * xCells);
    const real_t dt = real_c(1e-3);
    const uint_t timeSteps = uint_c(1e3);
    const uint_t vtkWriteFrequency = uint_c(0);

    // Block storage setup
    auto aabb = math::AABB(real_c(0.5) * dx, real_c(0.5) * dx, real_c(0.5) * dx,
                          xSize - real_c(0.5) * dx, ySize - real_c(0.5) * dx, zSize - real_c(0.5) * dx);

    shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGrid(
        aabb, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, true, false, false, false);

    // Field initialization
    BlockDataID uFieldId = field::addToStorage<ScalarField>(blocks, "u", real_c(300.0), field::fzyx, uint_c(1));
    BlockDataID uTmpFieldId = field::addToStorage<ScalarField>(blocks, "u_tmp", real_c(300.0), field::fzyx, uint_c(1));
    BlockDataID alphaFieldId = field::addToStorage<ScalarField>(blocks, "alpha", real_c(0.0), field::fzyx, uint_c(1));

    // Communication
    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commScheme(blocks);
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(uFieldId));

    // Boundary conditions
    initDirichletBoundaries(blocks, uFieldId, uTmpFieldId);

    // Timeloop
    SweepTimeloop timeloop(blocks, timeSteps);
    timeloop.add() << BeforeFunction(commScheme, "Communication")
                   << Sweep(HeatEquationKernelWithMaterialCPU(alphaFieldId, uFieldId, uTmpFieldId, dt, dx), "HeatSolverCPU")
                   << AfterFunction([blocks, uFieldId, uTmpFieldId]() {swapFields(*blocks, uFieldId, uTmpFieldId);}, "Swap");

    if (vtkWriteFrequency > 0) {
        std::string scalingType = weakScaling ? "weak" : "strong";
        std::string vtkFilename = "vtk_CPU_" + scalingType + "_" + std::to_string(problemSize) + 
                                "cells_" + std::to_string(numProcesses) + "proc(s)";
        std::string vtkDirectory = "vtk_out_cpu_" + scalingType + "_" + std::to_string(problemSize) + 
                                "cells_" + std::to_string(numProcesses) + "proc(s)";
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, vtkFilename, vtkWriteFrequency, 0, false, vtkDirectory,
                                                       "simulation_step", false, true, true, false, 0);
        auto tempWriter = make_shared<field::VTKWriter<ScalarField>>(uFieldId, "temperature");
        vtkOutput->addCellDataWriter(tempWriter);
        auto alphaWriter = make_shared<field::VTKWriter<ScalarField>>(alphaFieldId, "thermal_diffusivity");
        vtkOutput->addCellDataWriter(alphaWriter);
        timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output CPU");
    }

    // Benchmark
    const uint_t warmupSteps = uint_c(5);
    for (uint_t i = 0; i < warmupSteps; ++i)
        timeloop.singleStep();

    WcTimingPool timeloopTiming;
    WcTimer simTimer;

    WALBERLA_MPI_WORLD_BARRIER()
    WALBERLA_LOG_INFO_ON_ROOT("Starting CPU simulation with " << timeSteps << " time steps")

    simTimer.start();
    timeloop.run(timeloopTiming);
    simTimer.end();

    double simTime = simTimer.max();
    WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); }

    const auto reducedTimeloopTiming = timeloopTiming.getReduced();
    WALBERLA_LOG_RESULT_ON_ROOT("=== Performance Results ===")
    WALBERLA_LOG_RESULT_ON_ROOT("Total simulation time: " << simTime << " seconds")

    // Calculate performance metrics
    uint_t cellsPerProcess = xCells * yCells * zCells;
    uint_t totalCells = cellsPerProcess * numProcesses;
    uint_t totalCellUpdates = timeSteps * totalCells;
    
    double mlupsPerProcess = (timeSteps * cellsPerProcess) / (simTime * 1e6);
    double totalMLUPS = totalCellUpdates / (simTime * 1e6);
    
    WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per process: " << mlupsPerProcess)
    WALBERLA_LOG_RESULT_ON_ROOT("Total MLUPS: " << totalMLUPS)
    WALBERLA_LOG_RESULT_ON_ROOT("Time per timestep: " << (simTime / timeSteps * 1000) << " ms")
    WALBERLA_LOG_RESULT_ON_ROOT("Cells per process: " << cellsPerProcess)
    WALBERLA_LOG_RESULT_ON_ROOT("Total cells: " << totalCells)
    
    // Timing breakdown
    WALBERLA_LOG_RESULT_ON_ROOT("Detailed timing:\n" << *reducedTimeloopTiming)

    return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
