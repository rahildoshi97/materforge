//======================================================================================================================
//
//! \file CouetteFlow.cpp
//! \author Rahil Doshi
//
//======================================================================================================================

#include <iostream>

#include "core/all.h"
#include "blockforest/all.h"
#include "stencil/all.h"
#include "field/all.h"
#include "timeloop/all.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/communication/StencilRestrictedPackInfo.h"
#include "domain_decomposition/SharedSweep.h"
#include "walberla/experimental/Sweep.hpp"

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
   auto noSlip = gen::NoSlipFactory{blocks, pdfsId}.selectLinks([&](auto link) {
      return intersectsLowerWall(link);
   });

   // Upper wall: moving wall with velocity u_max (UBB - velocity bounce back)
   auto ubb = gen::UBBFactory{blocks, pdfsId, channelVelocity}.selectLinks([&](auto link){
      return intersectsUpperWall(link);
   });

   // Timeloop setup
   SweepTimeloop loop{blocks->getBlockStorage(), numTimesteps};
   loop.add() << Sweep(streamCollide) << AfterFunction(comm);
   loop.add() << Sweep(noSlip);
   loop.add() << Sweep(ubb);

   RemainingTimeLogger logger{numTimesteps};
   loop.addFuncAfterTimeStep(logger);

   // VTK Output
   Config::BlockHandle outputParams = config->getBlock("Output");
   const uint_t vtkWriteFrequency = outputParams.getParameter<uint_t>("vtkWriteFrequency", 0);

   if (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out_couette",
                                                       "simulation_step", false, true, true, false, 0);
      
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

   // Run the simulation
   WALBERLA_LOG_INFO_ON_ROOT("Commencing Couette flow simulation with " << numTimesteps << " timesteps");
   loop.run();

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
   }

   WALBERLA_LOG_INFO_ON_ROOT("Simulation completed successfully!");
}

} // namespace CouetteFlow

int main(int argc, char **argv)
{
   CouetteFlow::run(argc, argv);
   return EXIT_SUCCESS;
}
