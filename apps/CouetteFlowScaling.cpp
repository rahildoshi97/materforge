//======================================================================================================================
//! \file CouetteFlowScaling.cpp
//! \author Rahil Doshi
//! \brief Unified CPU/GPU scaling benchmark for Couette flow
//!
//! Usage: ./CouetteFlowScaling config.prm
//!
//! This is a SINGLE executable that runs on both CPU and GPU based on build configuration.
//! The sweep code generation automatically produces CPU or GPU code depending on TARGET_PLATFORM.
//!
//! Build options:
//!   CPU:      cmake --preset lumi-cpu && cmake --build --preset lumi-cpu-build
//!   GPU-HIP:  cmake --preset lumi-gpu && cmake --build --preset lumi-gpu-build
//!   GPU-CUDA: cmake --preset lumi-cuda && cmake --build --preset lumi-cuda-build
//======================================================================================================================

#include <string>
#include <cmath>
#include "core/all.h"
#include "blockforest/all.h"
#include "stencil/all.h"
#include "field/all.h"
#include "timeloop/all.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/communication/StencilRestrictedPackInfo.h"
#include "domain_decomposition/SharedSweep.h"
#include "walberla/experimental/Sweep.hpp"
#include "core/timing/TimingPool.h"
#include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h"

// GPU support (if available)
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/FieldCopy.h"
#include "gpu/HostFieldAllocator.h"
#include "gpu/communication/UniformGPUScheme.h"
#include "gpu/communication/MemcpyPackInfo.h"
#endif
#include "gen/CouetteFlowSweeps.hpp"

namespace CouetteFlow
{

#define use_materForge 1

using namespace walberla;

// Type definitions
using ScalarField_T = field::GhostLayerField<real_t, 1>;
using VectorField_T = field::GhostLayerField<real_t, 3>;

using LbStencil = stencil::D3Q19;
using PdfField_T = field::GhostLayerField<real_t, LbStencil::Q>;

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
using CommScheme = gpu::communication::UniformGPUScheme<LbStencil>;
using GPUField = gpu::GPUField<real_t>;
#else
using CommScheme = blockforest::communication::UniformBufferedScheme<LbStencil>;
#endif
using PdfsPackInfo = field::communication::StencilRestrictedPackInfo<PdfField_T, LbStencil>;

void mysleep(IBlock*){
    sleep(5);
}

//======================================================================
// MAIN FUNCTION
//======================================================================
void run(int argc, char **argv)
{
    Environment env{argc, argv};
    auto cfgFile = env.config();
    
    if (!cfgFile) {
        WALBERLA_ABORT("Usage: " << argv[0] << " path-to-configuration-file\n");
    }
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    gpu::selectDeviceBasedOnMpiRank();
    WALBERLA_LOG_INFO_ON_ROOT("GPU support enabled");
#else
    WALBERLA_LOG_INFO_ON_ROOT("CPU-only build");
#endif
    
    // Log build information
    WALBERLA_LOG_INFO_ON_ROOT(*cfgFile);
    
    const uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
    
    //======================================================================
    // READ CONFIGURATION - All parameters from .prm file
    // Bash scripts override via: -DomainSetup.blocks <x,y,z>
    //======================================================================
    
    Config::BlockHandle domainSetup = cfgFile->getBlock("DomainSetup");
    const Vector3<uint_t> numBlocks = domainSetup.getParameter<Vector3<uint_t>>("blocks");
    const Vector3<uint_t> cellsPerBlock = domainSetup.getParameter<Vector3<uint_t>>("cellsPerBlock");
    const Vector3<bool> periodic = domainSetup.getParameter<Vector3<bool>>("periodic");
    
    Config::BlockHandle simParams = cfgFile->getBlock("Parameters");
    const real_t latticeViscosity = simParams.getParameter<real_t>("nu");
    const real_t channelVelocity = simParams.getParameter<real_t>("u_max");
    const real_t errorThreshold = simParams.getParameter<real_t>("errorThreshold");
    uint_t numTimesteps = simParams.getParameter<uint_t>("timesteps");
    const real_t T_bottom = simParams.getParameter<real_t>("T_bottom");
    const real_t T_top = simParams.getParameter<real_t>("T_top");
    
    Config::BlockHandle outputParams = cfgFile->getBlock("Output");
    const uint_t vtkWriteFrequency = outputParams.getParameter<uint_t>("vtkWriteFrequency", 0);
    
    Config::BlockHandle perfParams = cfgFile->getBlock("Performance");
    const uint_t performanceLogFrequency = perfParams.getParameter<uint_t>("performanceLogFrequency", 0);
    
    //======================================================================
    // VALIDATE CONFIGURATION
    //======================================================================

    WALBERLA_CHECK_EQUAL(numBlocks[0] * numBlocks[1] * numBlocks[2], numProcesses,
                        "Number of blocks (" << numBlocks[0] << "x" << numBlocks[1] << "x" << numBlocks[2]
                        << " = " << numBlocks[0] * numBlocks[1] * numBlocks[2]
                        << ") must equal number of MPI processes (" << numProcesses << ")");
    
    const uint_t totalCellsX = numBlocks[0] * cellsPerBlock[0];
    const uint_t totalCellsY = numBlocks[1] * cellsPerBlock[1];
    const uint_t totalCellsZ = numBlocks[2] * cellsPerBlock[2];
    const uint_t totalCells = totalCellsX * totalCellsY * totalCellsZ;
    const uint_t cellsPerProcess = cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2];
    
    //======================================================================
    // LOG CONFIGURATION
    //======================================================================
    
    WALBERLA_LOG_INFO_ON_ROOT("=== Couette Flow Scaling Configuration ===");
    WALBERLA_LOG_INFO_ON_ROOT("Number of processes: " << numProcesses);
    WALBERLA_LOG_INFO_ON_ROOT("Block grid: " << numBlocks[0] << " x " << numBlocks[1] << " x " << numBlocks[2]);
    WALBERLA_LOG_INFO_ON_ROOT("Cells per block: " << cellsPerBlock[0] << " x " << cellsPerBlock[1] << " x " << cellsPerBlock[2]);
    WALBERLA_LOG_INFO_ON_ROOT("Cells per process: " << cellsPerProcess);
    WALBERLA_LOG_INFO_ON_ROOT("Total domain: " << totalCellsX << " x " << totalCellsY << " x " << totalCellsZ << " = " << totalCells);
    WALBERLA_LOG_INFO_ON_ROOT("Timesteps: " << numTimesteps);
    WALBERLA_LOG_INFO_ON_ROOT("Lattice viscosity: " << latticeViscosity);
    WALBERLA_LOG_INFO_ON_ROOT("Channel velocity: " << channelVelocity);
    
    //======================================================================
    // CREATE BLOCK STRUCTURE
    //======================================================================
    
    auto blocks = blockforest::createUniformBlockGrid(
        numBlocks[0], numBlocks[1], numBlocks[2],               // numberOfBlocks
        cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],   // numberOfCellsPerBlock
        real_t(1.0),                                            // dx
        true,                                                   // oneBlockPerProcess
        periodic[0], periodic[1], periodic[2],                  // periodicity
        false                                                   // keepGlobalBlockInformation       
    );
    
    //======================================================================
    // FIELD SETUP
    //======================================================================
    
    BlockDataID pdfsId = field::addToStorage<PdfField_T>(blocks, "pdfs", real_c(0.0), field::fzyx, 1);
    BlockDataID rhoId = field::addToStorage<ScalarField_T>(blocks, "rho", real_c(1.0), field::fzyx, 0);
    BlockDataID uId = field::addToStorage<VectorField_T>(blocks, "u", real_c(0.0), field::fzyx, 0);
    BlockDataID tempId = field::addToStorage<ScalarField_T>(blocks, "temperature", real_c(300.0), field::fzyx, 0);
    BlockDataID viscId = field::addToStorage<ScalarField_T>(blocks, "viscosity", latticeViscosity, field::fzyx, 0);
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    // GPU field setup
    BlockDataID gpuPdfsId = gpu::addGPUFieldToStorage<PdfField_T>(blocks, pdfsId, "pdfs_gpu", true);
    BlockDataID gpuRhoId = gpu::addGPUFieldToStorage<ScalarField_T>(blocks, rhoId, "rho_gpu", true);
    BlockDataID gpuUId = gpu::addGPUFieldToStorage<VectorField_T>(blocks, uId, "u_gpu", true);
    BlockDataID gpuTempId = gpu::addGPUFieldToStorage<ScalarField_T>(blocks, tempId, "temp_gpu", true);
    BlockDataID gpuViscId = gpu::addGPUFieldToStorage<ScalarField_T>(blocks, viscId, "visc_gpu", true);
#endif
    
    //======================================================================
    // INITIALIZE FIELDS
    //======================================================================
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    // GPU: Initialize with GPU field IDs
    auto zeroVelocity = gen::Couette::SetZeroVelocity{gpuRhoId, gpuUId};
    auto initTemperature = gen::Couette::InitializeTemperature{blocks, gpuTempId, T_bottom, T_top};
    auto initializePdfs = gen::Couette::InitPdfs{gpuRhoId, gpuPdfsId, gpuUId};
    
    WALBERLA_LOG_INFO_ON_ROOT("Initializing fields on GPU...");
    for (auto &b : *blocks) {
        zeroVelocity(&b);
        initTemperature(&b);
        initializePdfs(&b);
    }
#else
    // CPU: Initialize with CPU field IDs
    auto zeroVelocity = gen::Couette::SetZeroVelocity{rhoId, uId};
    auto initTemperature = gen::Couette::InitializeTemperature{blocks, tempId, T_bottom, T_top};
    auto initializePdfs = gen::Couette::InitPdfs{rhoId, pdfsId, uId};
    
    WALBERLA_LOG_INFO_ON_ROOT("Initializing fields on CPU...");
    for (auto &b : *blocks) {
        zeroVelocity(&b);
        initTemperature(&b);
        initializePdfs(&b);
    }
#endif
    
    //======================================================================
    // CREATE STREAM-COLLIDE SWEEP
    //======================================================================
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    std::shared_ptr<gen::Couette::StreamCollide> streamCollidePtr;
    #if use_materForge
        WALBERLA_LOG_INFO_ON_ROOT("Using temperature-dependent viscosity (MaterForge) - GPU");
        streamCollidePtr = std::make_shared<gen::Couette::StreamCollide>(
            gpuRhoId, gpuPdfsId, gpuTempId, gpuUId, gpuViscId
        );
    #else
        WALBERLA_LOG_INFO_ON_ROOT("Using constant viscosity - GPU");
        streamCollidePtr = std::make_shared<gen::Couette::StreamCollide>(
            gpuRhoId, gpuPdfsId, gpuUId, gpuViscId
        );
    #endif
    auto streamCollide = makeSharedSweep(streamCollidePtr);
#else
    std::shared_ptr<gen::Couette::StreamCollide> streamCollidePtr;
    #if use_materForge
        WALBERLA_LOG_INFO_ON_ROOT("Using temperature-dependent viscosity (MaterForge) - CPU");
        streamCollidePtr = std::make_shared<gen::Couette::StreamCollide>(
            rhoId, pdfsId, tempId, uId, viscId
        );
    #else
        WALBERLA_LOG_INFO_ON_ROOT("Using constant viscosity - CPU");
        streamCollidePtr = std::make_shared<gen::Couette::StreamCollide>(
            rhoId, pdfsId, uId, viscId
        );
    #endif
    auto streamCollide = makeSharedSweep(streamCollidePtr);
#endif
    
    //======================================================================
    // COMMUNICATION SETUP
    //======================================================================
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    gpu::communication::UniformGPUScheme<LbStencil> communication(blocks, true);
    communication.addPackInfo(std::make_shared<gpu::communication::MemcpyPackInfo<GPUField>>(gpuPdfsId));
#else
    blockforest::communication::UniformBufferedScheme<LbStencil> communication(blocks);
    communication.addPackInfo(std::make_shared<PdfsPackInfo>(pdfsId));
#endif
    
    //======================================================================
    // BOUNDARY CONDITIONS
    //======================================================================
    
    auto intersectsUpperWall = [&](auto link) -> bool {
        blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
        return link.wallCell.z() > blocks->getDomainCellBB().zMax();
    };
    
    auto intersectsLowerWall = [&](auto link) -> bool {
        blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
        return link.wallCell.z() < blocks->getDomainCellBB().zMin();
    };
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    auto noSlip = gen::NoSlipFactory{blocks, gpuPdfsId}.fromLinks([&](auto link) {
        return intersectsLowerWall(link);
    });
    
    auto ubb = gen::UBBFactory{blocks, gpuPdfsId, channelVelocity}.fromLinks([&](auto link){
        return intersectsUpperWall(link);
    });
#else
    auto noSlip = gen::NoSlipFactory{blocks, pdfsId}.fromLinks([&](auto link) {
        return intersectsLowerWall(link);
    });
    
    auto ubb = gen::UBBFactory{blocks, pdfsId, channelVelocity}.fromLinks([&](auto link){
        return intersectsUpperWall(link);
    });
#endif
    
    //======================================================================
    // TIMELOOP SETUP
    //======================================================================
    
    using walberla::lbm_mesapd_coupling::psm::gpu::deviceSyncWrapper;
    SweepTimeloop loop{blocks->getBlockStorage(), numTimesteps};
    /*if (numProcesses > 1) {
        loop.addFuncBeforeTimeStep(communication.getCommunicateFunctor(), "LBM Communication");
    }*/
    loop.addFuncBeforeTimeStep([&]() { WALBERLA_MPI_BARRIER(); }, "Barrier before Communication");
    loop.addFuncBeforeTimeStep(communication.getCommunicateFunctor(), "LBM Communication");
    //loop.add() << Sweep(deviceSyncWrapper(streamCollide), "StreamCollide");
    loop.add() << Sweep(deviceSyncWrapper([](IBlock*){sleep(5);}), "StreamCollide");
    loop.add() << Sweep(deviceSyncWrapper(noSlip), "NoSlip");
    loop.add() << Sweep(deviceSyncWrapper(ubb), "UBB");
    
    //======================================================================
    // GPU SYNCHRONIZATION HELPER FOR VTK OUTPUT
    //======================================================================
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    auto syncGPUFieldsToCPU = [&]() {
        // Copy GPU fields back to CPU for VTK output
        gpu::fieldCpy<ScalarField_T, GPUField>(blocks, rhoId, gpuRhoId);
        gpu::fieldCpy<VectorField_T, GPUField>(blocks, uId, gpuUId);
        gpu::fieldCpy<ScalarField_T, GPUField>(blocks, tempId, gpuTempId);
        gpu::fieldCpy<ScalarField_T, GPUField>(blocks, viscId, gpuViscId);
    };
#endif
    
    //======================================================================
    // VTK OUTPUT (Optional)
    //======================================================================
    
    if (vtkWriteFrequency > 0) {
        std::string vtkName = "couette_scaling";
    #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
        vtkName += "_gpu";
    #else
        vtkName += "_cpu";
    #endif
    #if use_materForge
        vtkName += "_tempdep";
    #else
        vtkName += "_const";
    #endif
        
        std::string vtkOutputDir = "testvtk";
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, vtkName, vtkWriteFrequency, 0, 
                                                        false, vtkOutputDir, "simulation_step", 
                                                        false, true, true, false, 0);
        
        auto densityWriter = make_shared<field::VTKWriter<ScalarField_T>>(rhoId, "density");
        vtkOutput->addCellDataWriter(densityWriter);
        auto velWriter = make_shared<field::VTKWriter<VectorField_T>>(uId, "velocity");
        vtkOutput->addCellDataWriter(velWriter);
        auto tempWriter = make_shared<field::VTKWriter<ScalarField_T>>(tempId, "temperature");
        vtkOutput->addCellDataWriter(tempWriter);
        auto viscWriter = make_shared<field::VTKWriter<ScalarField_T>>(viscId, "viscosity");
        vtkOutput->addCellDataWriter(viscWriter);
        
    #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
        // Sync GPU and copy data to CPU for VTK to read
        vtkOutput->addBeforeFunction([blocks, uId, gpuUId, rhoId, gpuRhoId, tempId, gpuTempId, viscId, gpuViscId]() {
            WALBERLA_GPU_CHECK(gpuDeviceSynchronize());
            gpu::fieldCpy<VectorField_T, gpu::GPUField<real_t>>(blocks, uId, gpuUId);
            gpu::fieldCpy<ScalarField_T, gpu::GPUField<real_t>>(blocks, rhoId, gpuRhoId);
            gpu::fieldCpy<ScalarField_T, gpu::GPUField<real_t>>(blocks, tempId, gpuTempId);
            gpu::fieldCpy<ScalarField_T, gpu::GPUField<real_t>>(blocks, viscId, gpuViscId);
        });
    #endif
        
        loop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
    }
    
    //======================================================================
    // WARMUP
    //======================================================================
    
    const uint_t warmupSteps = uint_c(5);
    WALBERLA_LOG_INFO_ON_ROOT("Running " << warmupSteps << " warmup iterations...");
    for (uint_t i = 0; i < warmupSteps; ++i)
        loop.singleStep();
    
    //======================================================================
    // ADD PERFORMANCE LOGGER
    //======================================================================
    
    WcTimingPool timeloopTiming;
    
    if (performanceLogFrequency > 0) {
        loop.addFuncAfterTimeStep([&timeloopTiming, performanceLogFrequency, numTimesteps, totalCells]() {
            static uint_t step = 0;
            step++;
            if (step % performanceLogFrequency == 0) {
                const auto reducedTiming = timeloopTiming.getReduced();
                WALBERLA_LOG_RESULT_ON_ROOT("Step " << step << "/" << numTimesteps);
            }
        }, "Performance Logger");
    }
    
    // Add remaining time logger
    RemainingTimeLogger logger{numTimesteps};
    loop.addFuncAfterTimeStep(logger, "Remaining Time Logger");
    
    //======================================================================
    // TIMED RUN
    //======================================================================
    
    WcTimer simTimer;
    
    WALBERLA_MPI_WORLD_BARRIER()
    WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << numTimesteps << " timesteps")
    
    simTimer.start();
    loop.run(timeloopTiming);
    simTimer.end();
    
    WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
    
    //======================================================================
    // PERFORMANCE METRICS
    //======================================================================
    
    double simTime = simTimer.max();
    WALBERLA_MPI_SECTION() { 
        walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); 
    }
    
    const auto reducedTimeloopTiming = timeloopTiming.getReduced();
    
    const uint_t totalLatticeUpdates = numTimesteps * totalCells;
    const double mlupsPerProcess = (numTimesteps * cellsPerProcess) / (simTime * 1e6);
    const double totalMLUPS = totalLatticeUpdates / (simTime * 1e6);
    
    WALBERLA_LOG_RESULT_ON_ROOT("=== Performance Results ===")
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    WALBERLA_LOG_RESULT_ON_ROOT("Backend: GPU")
#else
    WALBERLA_LOG_RESULT_ON_ROOT("Backend: CPU")
#endif
    WALBERLA_LOG_RESULT_ON_ROOT("Total simulation time: " << simTime << " seconds")
    WALBERLA_LOG_RESULT_ON_ROOT("Domain size: " << totalCellsX << " x " << totalCellsY << " x " << totalCellsZ)
    WALBERLA_LOG_RESULT_ON_ROOT("Total cells: " << totalCells)
    WALBERLA_LOG_RESULT_ON_ROOT("Cells per process: " << cellsPerProcess)
    WALBERLA_LOG_RESULT_ON_ROOT("Number of processes: " << numProcesses)
    WALBERLA_LOG_RESULT_ON_ROOT("Timesteps: " << numTimesteps)
    WALBERLA_LOG_RESULT_ON_ROOT("Total lattice updates: " << totalLatticeUpdates)
    WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per process: " << mlupsPerProcess)
    WALBERLA_LOG_RESULT_ON_ROOT("Total MLUPS: " << totalMLUPS)
    WALBERLA_LOG_RESULT_ON_ROOT("Time per timestep: " << (simTime / numTimesteps * 1000.0) << " ms")
    WALBERLA_LOG_RESULT_ON_ROOT("")
    WALBERLA_LOG_RESULT_ON_ROOT("=== Detailed Timing ===")
    WALBERLA_LOG_RESULT_ON_ROOT(*reducedTimeloopTiming)
    
    WALBERLA_LOG_INFO_ON_ROOT("Simulation completed successfully!");
}

} // namespace CouetteFlow

int main(int argc, char **argv)
{
    CouetteFlow::run(argc, argv);
    return EXIT_SUCCESS;
}
