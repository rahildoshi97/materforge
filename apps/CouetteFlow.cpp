//======================================================================================================================
//! \file CouetteFlow.cpp
//! \author Rahil Doshi
//======================================================================================================================

#include <iostream>

#include "core/all.h"
#include "blockforest/all.h"
#include "stencil/all.h"
#include "field/all.h"
#include "timeloop/all.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/communication/StencilRestrictedPackInfo.h"
//#include "field/communication/PackInfo.h"
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

void run(int argc, char **argv)
{
   Environment env{argc, argv};
   auto config = env.config();
   auto blocks = blockforest::createUniformBlockGridFromConfig(config);

   // Read domain configuration for performance metrics
   Config::BlockHandle domainSetup = config->getBlock("DomainSetup");
   Vector3<uint_t> numBlocks = domainSetup.getParameter<Vector3<uint_t>>("blocks");
   Vector3<uint_t> cellsPerBlock = domainSetup.getParameter<Vector3<uint_t>>("cellsPerBlock");
   
   // Calculate total domain size
   const uint_t xCells = numBlocks[0] * cellsPerBlock[0];
   const uint_t yCells = numBlocks[1] * cellsPerBlock[1];
   const uint_t zCells = numBlocks[2] * cellsPerBlock[2];
   const uint_t totalCells = xCells * yCells * zCells;
   const uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());

   // Field setup: PDFs, density, velocity, temperature and viscosity
   BlockDataID pdfsId = field::addToStorage<PdfField_T>(blocks, "pdfs", real_c(0.0), field::fzyx, 1);
   BlockDataID rhoId = field::addToStorage<ScalarField_T>(blocks, "rho", real_c(1.0), field::fzyx, 0);
   BlockDataID uId = field::addToStorage<VectorField_T>(blocks, "u", real_c(0.0), field::fzyx, 0);
   BlockDataID tempId = field::addToStorage<ScalarField_T>(blocks, "temperature", real_c(300.0), field::fzyx, 0);
   BlockDataID viscId = field::addToStorage<ScalarField_T>(blocks, "viscosity", real_c(1.0/6.0), field::fzyx, 0);

   // Read parameters from the Parameters block
   Config::BlockHandle simParams = config->getBlock("Parameters");
   const real_t channelVelocity = simParams.getParameter<real_t>("u_max");
   const real_t errorThreshold = simParams.getParameter<real_t>("errorThreshold");
   const uint_t numTimesteps = simParams.getParameter<uint_t>("timesteps");
   // const bool useMaterForge = simParams.getParameter<bool>("useMaterForge", false);
   const real_t T_bottom = simParams.getParameter<real_t>("T_bottom");
   const real_t T_top = simParams.getParameter<real_t>("T_top");

   WALBERLA_LOG_INFO_ON_ROOT("=== Domain Configuration ===")
   WALBERLA_LOG_INFO_ON_ROOT("Blocks: " << numBlocks[0] << " x " << numBlocks[1] << " x " << numBlocks[2])
   WALBERLA_LOG_INFO_ON_ROOT("Cells per block: " << cellsPerBlock[0] << " x " << cellsPerBlock[1] << " x " << cellsPerBlock[2])
   WALBERLA_LOG_INFO_ON_ROOT("Total domain: " << xCells << " x " << yCells << " x " << zCells << " = " << totalCells << " cells")
   WALBERLA_LOG_INFO_ON_ROOT("Number of processes: " << numProcesses)
   WALBERLA_LOG_INFO_ON_ROOT("Cells per process: " << totalCells / numProcesses)

   // Prepare sweep functors using code-generated classes
   auto setAnalytical = gen::Couette::SetAnalyticalSolution{blocks, rhoId, uId, channelVelocity};
   auto initTemperature = gen::Couette::InitializeTemperature{blocks, tempId, T_bottom, T_top};
   auto initializePdfs = gen::Couette::InitPdfs{rhoId, pdfsId, uId};

   auto streamCollide = makeSharedSweep(
      std::make_shared<gen::Couette::StreamCollide>(
         rhoId, pdfsId, uId, viscId  // constant viscosity
         //rhoId, pdfsId, tempId, uId, viscId  // temperature-dependent viscosity
      )
   );

   // Set up initial state
   WALBERLA_LOG_INFO_ON_ROOT("Initializing fields...");
   for (auto &b : *blocks)
   {
      setAnalytical(&b);
      initTemperature(&b);
      initializePdfs(&b);
   }

   // Set up ghost layer communication
   CommScheme comm{blocks};
   comm.addPackInfo(std::make_shared<PdfsPackInfo>(pdfsId));
   //blockforest::communication::UniformBufferedScheme<stencil::D3Q19> comm(blocks);
   //comm.addPackInfo(make_shared<field::communication::PackInfo<PdfField_T>>(pdfsId));

   // Set up boundary conditions using code-generated factories
   auto intersectsUpperWall = [&](auto link) -> bool {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      return link.wallCell.z() > blocks->getDomainCellBB().zMax();
   };

   auto intersectsLowerWall = [&](auto link) -> bool {
      blocks->transformBlockLocalToGlobalCell(link.wallCell, link.block);
      return link.wallCell.z() < blocks->getDomainCellBB().zMin();
   };

   // Lower wall: no-slip boundary condition
   auto noSlip = gen::NoSlipFactory{blocks, pdfsId}.fromLinks([&](auto link) {
      return intersectsLowerWall(link);
   });

   // Upper wall: moving wall with velocity u_max (UBB - velocity bounce back)
   auto ubb = gen::UBBFactory{blocks, pdfsId, channelVelocity}.fromLinks([&](auto link){
      return intersectsUpperWall(link);
   });

   // Timeloop setup
   SweepTimeloop loop{blocks->getBlockStorage(), numTimesteps};
   loop.add() << Sweep(streamCollide, "StreamCollide Sweep") << AfterFunction(comm, "Communication");
   loop.add() << Sweep(noSlip, "NoSlip Sweep");
   loop.add() << Sweep(ubb, "UBB Sweep");

   // RemainingTimeLogger logger{numTimesteps};
   // loop.addFuncAfterTimeStep(logger);

   // VTK Output
   Config::BlockHandle outputParams = config->getBlock("Output");
   const uint_t vtkWriteFrequency = outputParams.getParameter<uint_t>("vtkWriteFrequency", 0);

   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, 
                                                       "vtk_out_couette", "simulation_step", 
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

   // ==================== BENCHMARK SECTION ====================
   
   // Warmup iterations
   const uint_t warmupSteps = uint_c(5);
   WALBERLA_LOG_INFO_ON_ROOT("Running " << warmupSteps << " warmup iterations...");
   for (uint_t i = 0; i < warmupSteps; ++i)
      loop.singleStep();

   // Add remaining time logger
   RemainingTimeLogger logger{numTimesteps};
   loop.addFuncAfterTimeStep(logger, "Remaining Time Logger");

   // Timed simulation run
   WcTimingPool timeloopTiming;
   WcTimer simTimer;

   WALBERLA_MPI_WORLD_BARRIER()
   WALBERLA_LOG_INFO_ON_ROOT("Starting Couette flow simulation with " << numTimesteps << " timesteps")
   
   simTimer.start();
   loop.run(timeloopTiming);
   simTimer.end();
   
   WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")

   // Get timing results
   double simTime = simTimer.max();
   WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); }

   const auto reducedTimeloopTiming = timeloopTiming.getReduced();
   
   // Calculate MLUPS (Million Lattice Updates Per Second)
   const uint_t cellsPerProcess = totalCells / numProcesses;
   const uint_t totalLatticeUpdates = numTimesteps * totalCells;
   const double mlupsPerProcess = (numTimesteps * cellsPerProcess) / (simTime * 1e6);
   const double totalMLUPS = totalLatticeUpdates / (simTime * 1e6);
   
   WALBERLA_LOG_RESULT_ON_ROOT("=== Performance Results ===")
   WALBERLA_LOG_RESULT_ON_ROOT("Total simulation time: " << simTime << " seconds")
   WALBERLA_LOG_RESULT_ON_ROOT("Domain size: " << xCells << " x " << yCells << " x " << zCells)
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

   // ==================== END BENCHMARK SECTION ====================

   // Check solution convergence
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
      //WALBERLA_CHECK_LESS(*velocityErrorLmax, errorThreshold);
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
