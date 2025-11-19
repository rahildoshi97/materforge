//======================================================================================================================
//! \file CouetteFlowScaling.cpp
//! \author Rahil Doshi
//! \brief Unified scaling benchmark for Couette flow (CPU/GPU)
//!
//! Usage: ./CouetteFlowScaling config.prm
//!
//! Bash scripts override parameters via command line:
//!   -DomainSetup.blocks <x,y,z>
//!   -DomainSetup.cellsPerBlock <x,y,z>
//!   -Parameters.timesteps <N>
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

// GPU support (if available)
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/FieldCopy.h"
#include "gpu/HostFieldAllocator.h"
#include "gpu/communication/UniformGPUScheme.h"
#include "gpu/communication/MemcpyPackInfo.h"
#include "gen/CouetteFlowSweepsGPU.hpp"
#else
#include "gen/CouetteFlowSweeps.hpp"
#endif

namespace CouetteFlow
{

constexpr bool use_materForge = true;

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

//======================================================================
// MAIN FUNCTION
//======================================================================
void run(int argc, char **argv)
{
    Environment env{argc, argv};
    auto cfgFile = env.config();
    
    if (!cfgFile) {
        WALBERLA_ABORT("Usage: " << argv << " path-to-configuration-file\n");
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
                        "Number of blocks (" << numBlocks << "x" << numBlocks << "x" << numBlocks
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
    auto gpuPdfsId = gpu::addGPUFieldToStorage<GPUField>(pdfsId, blocks, "pdfs_gpu", true);
    auto gpuRhoId = gpu::addGPUFieldToStorage<GPUField>(rhoId, blocks, "rho_gpu", true);
    auto gpuUId = gpu::addGPUFieldToStorage<GPUField>(uId, blocks, "u_gpu", true);
    auto gpuTempId = gpu::addGPUFieldToStorage<GPUField>(tempId, blocks, "temp_gpu", true);
    auto gpuViscId = gpu::addGPUFieldToStorage<GPUField>(viscId, blocks, "visc_gpu", true);
#endif
    
    //======================================================================
    // INITIALIZE FIELDS
    //======================================================================
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    auto zeroVelocity = gen::CouetteGPU::SetZeroVelocity{rhoId, uId};
    auto initTemperature = gen::CouetteGPU::InitializeTemperature{blocks, tempId, T_bottom, T_top};
    auto initializePdfs = gen::CouetteGPU::InitPdfs{rhoId, pdfsId, uId};
#else
    auto zeroVelocity = gen::Couette::SetZeroVelocity{rhoId, uId};
    auto initTemperature = gen::Couette::InitializeTemperature{blocks, tempId, T_bottom, T_top};
    auto initializePdfs = gen::Couette::InitPdfs{rhoId, pdfsId, uId};
#endif
    
    WALBERLA_LOG_INFO_ON_ROOT("Initializing fields...");
    for (auto &b : *blocks) {
        zeroVelocity(&b);
        initTemperature(&b);
        initializePdfs(&b);
    }
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    // Copy to GPU
    gpu::fieldCpy<PdfField_T, GPUField>(blocks, pdfsId, gpuPdfsId);
    gpu::fieldCpy<ScalarField_T, GPUField>(blocks, rhoId, gpuRhoId);
    gpu::fieldCpy<VectorField_T, GPUField>(blocks, uId, gpuUId);
    gpu::fieldCpy<ScalarField_T, GPUField>(blocks, tempId, gpuTempId);
    gpu::fieldCpy<ScalarField_T, GPUField>(blocks, viscId, gpuViscId);
#endif
    
    //======================================================================
    // CREATE STREAM-COLLIDE SWEEP
    //======================================================================
    
    auto streamCollide = [&]() {
        if constexpr (use_materForge) {
            WALBERLA_LOG_INFO_ON_ROOT("Using temperature-dependent viscosity (MaterForge)");
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            return makeSharedSweep(std::make_shared<gen::CouetteGPU::StreamCollide>(
                gpuRhoId, gpuPdfsId, gpuTempId, gpuUId, gpuViscId
            ));
#else
            return makeSharedSweep(std::make_shared<gen::Couette::StreamCollide>(
                rhoId, pdfsId, tempId, uId, viscId
            ));
#endif
        } else {
            WALBERLA_LOG_INFO_ON_ROOT("Using constant viscosity");
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            return makeSharedSweep(std::make_shared<gen::CouetteGPU::StreamCollide>(
                gpuRhoId, gpuPdfsId, gpuUId, gpuViscId
            ));
#else
            return makeSharedSweep(std::make_shared<gen::Couette::StreamCollide>(
                rhoId, pdfsId, uId, viscId
            ));
#endif
        }
    }();
    
    //======================================================================
    // COMMUNICATION SETUP
    //======================================================================
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    CommScheme comm{blocks, true};  // GPU-aware MPI
    comm.addPackInfo(std::make_shared<gpu::communication::MemcpyPackInfo<GPUField>>(gpuPdfsId));
#else
    CommScheme comm{blocks};
    comm.addPackInfo(std::make_shared<PdfsPackInfo>(pdfsId));
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
    auto noSlip = gen::NoSlipGPUFactory{blocks, gpuPdfsId}.fromLinks([&](auto link) {
        return intersectsLowerWall(link);
    });
    
    auto ubb = gen::UBBGPUFactory{blocks, gpuPdfsId, channelVelocity}.fromLinks([&](auto link){
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
    
    SweepTimeloop loop{blocks->getBlockStorage(), numTimesteps};
    loop.add() << Sweep(streamCollide, "StreamCollide") << AfterFunction(comm, "Communication");
    loop.add() << Sweep(noSlip, "NoSlip");
    loop.add() << Sweep(ubb, "UBB");
    
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
        if constexpr (use_materForge) {
            vtkName += "_tempdep";
        } else {
            vtkName += "_const";
        }
        
        std::string vtkOutputDir = vtkName + "_out";
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
    
    //======================================================================
    // CONVERGENCE CHECK
    //======================================================================
    
    WALBERLA_LOG_INFO_ON_ROOT("Checking for convergence...");
    
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
    // Copy velocity field back to CPU
    gpu::fieldCpy<GPUField, VectorField_T>(blocks, gpuUId, uId);
    
    auto computeVelocityError = gen::CouetteGPU::VelocityErrorLmax{
        blocks, uId, channelVelocity
    };
#else
    auto velocityErrorLmax = std::make_unique<real_t>(-std::numeric_limits<real_t>::infinity());
    auto computeVelocityError = gen::Couette::VelocityErrorLmax{
        blocks, uId, velocityErrorLmax.get(), channelVelocity
    };
#endif
    
    for(auto& b: *blocks) {
        computeVelocityError(&b);
    }
    
#ifndef WALBERLA_BUILD_WITH_GPU_SUPPORT
    mpi::reduceInplace(*velocityErrorLmax, mpi::MAX);
    
    WALBERLA_LOG_INFO_ON_ROOT("L-infinity error of x-velocity: " << *velocityErrorLmax);
    
    WALBERLA_ROOT_SECTION() {
        if (*velocityErrorLmax < errorThreshold) {
            WALBERLA_LOG_INFO_ON_ROOT("Convergence criterion met (error < " << errorThreshold << ")");
        } else {
            WALBERLA_LOG_WARNING_ON_ROOT("Solution has not converged to threshold " << errorThreshold);
        }
    }
#endif
    
    WALBERLA_LOG_INFO_ON_ROOT("Simulation completed successfully!");
}

} // namespace CouetteFlow

int main(int argc, char **argv)
{
    CouetteFlow::run(argc, argv);
    return EXIT_SUCCESS;
}
