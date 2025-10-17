//======================================================================================================================
//! \file CouetteFlowCPUScaling.cpp
//! \author Rahil Doshi
//! \brief CPU scaling benchmark for Couette flow using .prm configuration
//======================================================================================================================

#include <iostream>
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

#include "gen/CouetteFlowSweeps.hpp"

namespace CouetteFlow
{

using namespace walberla;

// Type definitions
using ScalarField_T = field::GhostLayerField<real_t, 1>;
using VectorField_T = field::GhostLayerField<real_t, 3>;

using LbStencil = stencil::D3Q19;
using PdfField_T = field::GhostLayerField<real_t, LbStencil::Q>;

using CommScheme = blockforest::communication::UniformBufferedScheme<LbStencil>;
using PdfsPackInfo = field::communication::StencilRestrictedPackInfo<PdfField_T, LbStencil>;

// Automatic process decomposition
std::tuple<uint_t, uint_t, uint_t> calculateProcessDecomposition(uint_t numProcesses) {
    uint_t procs_x = 1, procs_y = 1, procs_z = 1;
    uint_t base = uint_c(std::round(std::cbrt(real_c(numProcesses))));
    
    for (uint_t z = base; z >= 1; --z) {
        if (numProcesses % z == 0) {
            uint_t remaining = numProcesses / z;
            uint_t base_xy = uint_c(std::round(std::sqrt(real_c(remaining))));
            for (uint_t y = base_xy; y >= 1; --y) {
                if (remaining % y == 0) {
                    procs_x = remaining / y;
                    procs_y = y;
                    procs_z = z;
                    goto found;
                }
            }
        }
    }
    found:
    return std::make_tuple(procs_x, procs_y, procs_z);
}

void run(int argc, char **argv)
{
   Environment env{argc, argv};
   auto config = env.config();
   
   const uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   
   // Read parameters from config
   Config::BlockHandle simParams = config->getBlock("Parameters");
   const std::string scalingTest = simParams.getParameter<std::string>("scalingTest", "none");
   const uint_t problemSize = simParams.getParameter<uint_t>("problemSize", 64);
   const real_t channelVelocity = simParams.getParameter<real_t>("u_max");
   const real_t errorThreshold = simParams.getParameter<real_t>("errorThreshold");
   const uint_t numTimesteps = simParams.getParameter<uint_t>("timesteps");
   const real_t T_bottom = simParams.getParameter<real_t>("T_bottom");
   const real_t T_top = simParams.getParameter<real_t>("T_top");
   const real_t latticeViscosity = simParams.getParameter<real_t>("nu");
   
   // Read base domain configuration
   Config::BlockHandle domainSetup = config->getBlock("DomainSetup");
   Vector3<uint_t> baseCellsPerBlock = domainSetup.getParameter<Vector3<uint_t>>("cellsPerBlock");
   Vector3<bool> periodic = domainSetup.getParameter<Vector3<bool>>("periodic");
   
   // Determine if weak or strong scaling
   bool weakScaling = (scalingTest == "weak");
   bool strongScaling = (scalingTest == "strong");
   bool noScaling = (scalingTest == "none");
   
   uint_t xCells, yCells, zCells;
   uint_t xBlocks, yBlocks, zBlocks;
   
   // Automatic process decomposition
   uint_t procs_x, procs_y, procs_z;
   std::tie(procs_x, procs_y, procs_z) = calculateProcessDecomposition(numProcesses);
   
   if (noScaling) {
       // Use configuration from .prm file directly
       Vector3<uint_t> numBlocks = domainSetup.getParameter<Vector3<uint_t>>("blocks");
       xBlocks = numBlocks[0];
       yBlocks = numBlocks[1];
       zBlocks = numBlocks[2];
       
       xCells = baseCellsPerBlock[0];
       yCells = baseCellsPerBlock[1];
       zCells = baseCellsPerBlock[2];
       
       WALBERLA_LOG_INFO_ON_ROOT("Standard run (no scaling test)");
       
   } else if (weakScaling) {
       // Weak scaling: constant cells per process, blocks = process grid
       xCells = problemSize;
       yCells = problemSize / 2;  // Maintain 2:1:1 aspect ratio
       zCells = problemSize / 2;
       
       xBlocks = procs_x;
       yBlocks = procs_y;
       zBlocks = procs_z;
       
       WALBERLA_LOG_INFO_ON_ROOT("Weak Scaling Test");
       
   } else if (strongScaling) {
       // Strong scaling: fixed total domain, distribute across processes
       const uint_t totalCellsX = problemSize;
       const uint_t totalCellsY = problemSize / 2;
       const uint_t totalCellsZ = problemSize / 2;
       
       xBlocks = procs_x;
       yBlocks = procs_y;
       zBlocks = procs_z;
       
       xCells = totalCellsX / procs_x;
       yCells = totalCellsY / procs_y;
       zCells = totalCellsZ / procs_z;
       
       if (totalCellsX % procs_x != 0 || totalCellsY % procs_y != 0 || 
           totalCellsZ % procs_z != 0) {
           WALBERLA_ABORT("Problem size " << problemSize << 
               " not evenly divisible by process grid " << 
               procs_x << "x" << procs_y << "x" << procs_z);
       }
       
       WALBERLA_LOG_INFO_ON_ROOT("Strong Scaling Test");
   }
   
   const uint_t totalCellsX = xBlocks * xCells;
   const uint_t totalCellsY = yBlocks * yCells;
   const uint_t totalCellsZ = zBlocks * zCells;
   const uint_t totalCells = totalCellsX * totalCellsY * totalCellsZ;
   const uint_t cellsPerProcess = xCells * yCells * zCells;
   
   WALBERLA_LOG_INFO_ON_ROOT("=== CPU Scaling Configuration ===");
   WALBERLA_LOG_INFO_ON_ROOT("Test Type: " << scalingTest);
   WALBERLA_LOG_INFO_ON_ROOT("Processes: " << numProcesses);
   WALBERLA_LOG_INFO_ON_ROOT("Process Grid: " << procs_x << "x" << procs_y << "x" << procs_z);
   WALBERLA_LOG_INFO_ON_ROOT("Blocks: " << xBlocks << "x" << yBlocks << "x" << zBlocks);
   WALBERLA_LOG_INFO_ON_ROOT("Cells per block: " << xCells << "x" << yCells << "x" << zCells);
   WALBERLA_LOG_INFO_ON_ROOT("Cells per process: " << cellsPerProcess);
   WALBERLA_LOG_INFO_ON_ROOT("Total domain: " << totalCellsX << "x" << totalCellsY << "x" << totalCellsZ << " = " << totalCells);
   
   // Block storage setup
   const real_t dx = real_c(1.0);
   auto blocks = blockforest::createUniformBlockGrid(
       xBlocks, yBlocks, zBlocks,
       xCells, yCells, zCells,
       dx,
       true,            // oneBlockPerProcess
       periodic[0],     // xPeriodic (from .prm)
       periodic[1],     // yPeriodic (from .prm)
       periodic[2],     // zPeriodic (from .prm)
       false            // keepGlobalBlockInformation
   );
   
   // Field setup
   BlockDataID pdfsId = field::addToStorage<PdfField_T>(blocks, "pdfs", real_c(0.0), field::fzyx, 1);
   BlockDataID rhoId = field::addToStorage<ScalarField_T>(blocks, "rho", real_c(1.0), field::fzyx, 0);
   BlockDataID uId = field::addToStorage<VectorField_T>(blocks, "u", real_c(0.0), field::fzyx, 0);
   BlockDataID tempId = field::addToStorage<ScalarField_T>(blocks, "temperature", real_c(300.0), field::fzyx, 0);
   BlockDataID viscId = field::addToStorage<ScalarField_T>(blocks, "viscosity", latticeViscosity, field::fzyx, 0);
   
   // Initialize fields
   auto setAnalytical = gen::Couette::SetAnalyticalSolution{blocks, rhoId, uId, channelVelocity};
   auto initTemperature = gen::Couette::InitializeTemperature{blocks, tempId, T_bottom, T_top};
   auto initializePdfs = gen::Couette::InitPdfs{rhoId, pdfsId, uId};
   
   WALBERLA_LOG_INFO_ON_ROOT("Initializing fields...");
   for (auto &b : *blocks) {
      setAnalytical(&b);
      initTemperature(&b);
      initializePdfs(&b);
   }
   
   // Stream-collide sweep
   auto streamCollide = makeSharedSweep(
      std::make_shared<gen::Couette::StreamCollide>(
         rhoId, pdfsId, uId, viscId  // constant viscosity
         //rhoId, pdfsId, tempId, uId, viscId  // temperature-dependent viscosity
      )
   );
   
   // Communication
   CommScheme comm{blocks};
   comm.addPackInfo(std::make_shared<PdfsPackInfo>(pdfsId));
   
   // Boundary conditions
   auto intersectsUpperWall = [&](auto link) -> bool {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      return link.wallCell.z() > blocks->getDomainCellBB().zMax();
   };
   
   auto intersectsLowerWall = [&](auto link) -> bool {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      return link.wallCell.z() < blocks->getDomainCellBB().zMin();
   };
   
   auto noSlip = gen::NoSlipFactory{blocks, pdfsId}.fromLinks([&](auto link) {
      return intersectsLowerWall(link);
   });
   
   auto ubb = gen::UBBFactory{blocks, pdfsId, channelVelocity}.fromLinks([&](auto link){
      return intersectsUpperWall(link);
   });
   
   // Timeloop
   SweepTimeloop loop{blocks->getBlockStorage(), numTimesteps};
   loop.add() << Sweep(streamCollide, "StreamCollide") << AfterFunction(comm, "Communication");
   loop.add() << Sweep(noSlip, "NoSlip");
   loop.add() << Sweep(ubb, "UBB");
   
   // VTK output (optional, from config)
   Config::BlockHandle outputParams = config->getBlock("Output");
   const uint_t vtkWriteFrequency = outputParams.getParameter<uint_t>("vtkWriteFrequency", 0);
   
   if (vtkWriteFrequency > 0) {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, 
                                                       "vtk_out_couette_cpu", "simulation_step", 
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
   
   // Warmup
   const uint_t warmupSteps = uint_c(5);
   WALBERLA_LOG_INFO_ON_ROOT("Running " << warmupSteps << " warmup iterations...");
   for (uint_t i = 0; i < warmupSteps; ++i)
      loop.singleStep();
   
   // Add remaining time logger
   RemainingTimeLogger logger{numTimesteps};
   loop.addFuncAfterTimeStep(logger, "Remaining Time Logger");
   
   // Timed run
   WcTimingPool timeloopTiming;
   WcTimer simTimer;
   
   WALBERLA_MPI_WORLD_BARRIER()
   WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << numTimesteps << " timesteps")
   
   simTimer.start();
   loop.run(timeloopTiming);
   simTimer.end();
   
   WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
   
   // Get timing results
   double simTime = simTimer.max();
   WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); }
   
   const auto reducedTimeloopTiming = timeloopTiming.getReduced();
   
   // Performance metrics
   const uint_t totalLatticeUpdates = numTimesteps * totalCells;
   const double mlupsPerProcess = (numTimesteps * cellsPerProcess) / (simTime * 1e6);
   const double totalMLUPS = totalLatticeUpdates / (simTime * 1e6);
   
   WALBERLA_LOG_RESULT_ON_ROOT("=== Performance Results ===")
   WALBERLA_LOG_RESULT_ON_ROOT("Test type: " << scalingTest)
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
   
   // Check convergence
   WALBERLA_LOG_INFO_ON_ROOT("Checking for convergence...");
   auto velocityErrorLmax = std::make_unique<real_t>(-std::numeric_limits<real_t>::infinity());
   
   auto computeVelocityError = gen::Couette::VelocityErrorLmax{
      blocks, uId, velocityErrorLmax.get(), channelVelocity
   };
   
   for(auto& b: *blocks) {
      computeVelocityError(&b);
   }
   
   mpi::reduceInplace(*velocityErrorLmax, mpi::MAX);
   WALBERLA_LOG_INFO_ON_ROOT("L-infinity error of x-velocity: " << *velocityErrorLmax);
   
   WALBERLA_ROOT_SECTION() {
      if (*velocityErrorLmax < errorThreshold) {
         WALBERLA_LOG_INFO_ON_ROOT("Convergence criterion met (error < " << errorThreshold << ")");
      } else {
         WALBERLA_LOG_WARNING_ON_ROOT("Solution has not converged to threshold " << errorThreshold);
      }
   }
   
   WALBERLA_LOG_INFO_ON_ROOT("Simulation completed successfully!");
}

} // namespace CouetteFlow

int main(int argc, char **argv)
{
   CouetteFlow::run(argc, argv);
   return EXIT_SUCCESS;
}
