//======================================================================================================================
//
//! \file CodegenHeatEquationWithMaterial.cpp
//! \author Rahil Doshi <rahil.doshi@fau.de>
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

#include "stencil/D2Q9.h"

#include "timeloop/SweepTimeloop.h"

#include "gen/HeatEquationKernelWithMaterial.hpp"

#include "lbm_generated/evaluation/PerformanceEvaluation.h"

namespace walberla
{
typedef GhostLayerField< real_t, 1 > ScalarField;

void swapFields(StructuredBlockForest& blocks, BlockDataID uID, BlockDataID uTmpID)
{
   for (auto block = blocks.begin(); block != blocks.end(); ++block)
   {
      ScalarField* u     = block->getData< ScalarField >(uID);
      ScalarField* u_tmp = block->getData< ScalarField >(uTmpID);

      u->swapDataPointers(u_tmp);
   }
}

void initDirichletBoundaryNorth(shared_ptr< StructuredBlockForest > blocks, BlockDataID uId, BlockDataID uTmpId)
{
   for (auto block = blocks->begin(); block != blocks->end(); ++block)
   {
      if (blocks->atDomainYMaxBorder(*block))
      {
         ScalarField* u     = block->getData< ScalarField >(uId);
         ScalarField* u_tmp = block->getData< ScalarField >(uTmpId);

         CellInterval xyz = u->xyzSizeWithGhostLayer();

         xyz.yMin() = xyz.yMax();
         for (auto cell = xyz.begin(); cell != xyz.end(); ++cell)
         {
            const Vector3< real_t > p = blocks->getBlockLocalCellCenter(*block, *cell);
            //  Set the dirichlet boundary to f(x) = 1 + sin(x) * x^2
            real_t v          = real_c(2000.0 + 1500 * std::sin(2 * math::pi * p[0]) * p[0] * p[0]);
            u->get(*cell)     = v;
            u_tmp->get(*cell) = v;
         }
      }
   }
}

int main(int argc, char** argv)
{
   mpi::Environment env(argc, argv);

   /////////////////////////////
   /// SIMULATION PARAMETERS ///
   /////////////////////////////

   //  Ensure matching aspect ratios of cells and domain.
   const uint_t xCells = uint_c(25);
   const uint_t yCells = uint_c(25);
   const uint_t zCells = uint_c(1);

   const real_t xSize = real_c(1.0);
   const real_t ySize = real_c(1.0);

   const uint_t xBlocks = uint_c(1);
   const uint_t yBlocks = uint_c(1);
   const uint_t zBlocks = uint_c(1);

   const uint_t processes = uint_c(MPIManager::instance()->numProcesses());

   if (processes != xBlocks * yBlocks)
   { WALBERLA_ABORT("The number of processes must be equal to the number of blocks!"); }

   const real_t dx = xSize / real_c(xBlocks * xCells + uint_c(1));
   const real_t dy = ySize / real_c(yBlocks * yCells + uint_c(1));

   WALBERLA_CHECK_FLOAT_EQUAL(dx, dy);

   const real_t dt    = real_c(1);
   // const real_t thermal_diffusivity = real_c(1.0);
   uint_t timeSteps = uint_c(2e4);
   uint_t vtkWriteFrequency = uint_c(0);

   ///////////////////////////
   /// BLOCK STORAGE SETUP ///
   ///////////////////////////

   auto aabb = math::AABB(real_c(0.5) * dx, real_c(0.5) * dy, real_c(0.0), xSize - real_c(0.5) * dx,
                          ySize - real_c(0.5) * dy, dx);

   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
      aabb, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, true, false, false, false);

   //////////////
   /// FIELDS ///
   //////////////

   BlockDataID uFieldId    = field::addToStorage< ScalarField >(blocks, "u", real_c(300.0), field::fzyx, uint_c(1));
   BlockDataID uTmpFieldId = field::addToStorage< ScalarField >(blocks, "u_tmp", real_c(300.0), field::fzyx, uint_c(1));
   BlockDataID kFieldId = field::addToStorage< ScalarField >(blocks, "heat_conductivity", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID rhoFieldId = field::addToStorage< ScalarField >(blocks, "density", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID cpFieldId = field::addToStorage< ScalarField >(blocks, "heat_capacity", real_c(0.0), field::fzyx, uint_c(1));
   BlockDataID alphaFieldId = field::addToStorage< ScalarField >(blocks, "thermal_diffusivity", real_c(0.0), field::fzyx, uint_c(1));

   /////////////////////
   /// COMMUNICATION ///
   /////////////////////

   blockforest::communication::UniformBufferedScheme< stencil::D2Q9 > commScheme(blocks);
   commScheme.addPackInfo(make_shared< field::communication::PackInfo< ScalarField > >(uFieldId));

   //////////////////////////
   /// DIRICHLET BOUNDARY ///
   //////////////////////////

   initDirichletBoundaryNorth(blocks, uFieldId, uTmpFieldId);

   ////////////////////////
   /// NEUMANN BOUNDARY ///
   ////////////////////////

   pde::NeumannDomainBoundary< ScalarField > neumann(*blocks, uFieldId);

   neumann.excludeBoundary(stencil::N);
   neumann.excludeBoundary(stencil::B);
   neumann.excludeBoundary(stencil::T);

   ////////////////
   /// TIMELOOP ///
   ////////////////

   SweepTimeloop timeloop(blocks, timeSteps);

   timeloop.add() // << BeforeFunction(commScheme, "Communication")
                  << BeforeFunction(neumann, "Neumann Boundaries")
                  << Sweep(HeatEquationKernelWithMaterial(alphaFieldId, uFieldId, uTmpFieldId, dt, dx), "HeatEquationKernelWithMaterial")
                  // << Sweep(HeatEquationKernelWithMaterial(alphaFieldId, cpFieldId, kFieldId, rhoFieldId, uFieldId, uTmpFieldId, dt, dx), "HeatEquationKernelWithMaterial")
                  << AfterFunction([blocks, uFieldId, uTmpFieldId]() { swapFields(*blocks, uFieldId, uTmpFieldId); }, "Swap");

   if constexpr (vtkWriteFrequency > 0)
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                      "simulation_step", false, true, true, false, 0);

      auto tempWriter = make_shared< field::VTKWriter< ScalarField > >(uFieldId, "temperature");
      vtkOutput->addCellDataWriter(tempWriter);

      auto kWriter = make_shared< field::VTKWriter< ScalarField > >(kFieldId, "heat_conductivity");
      vtkOutput->addCellDataWriter(kWriter);

      auto rhoWriter = make_shared< field::VTKWriter< ScalarField > >(rhoFieldId, "density");
      vtkOutput->addCellDataWriter(rhoWriter);

      auto cpWriter = make_shared< field::VTKWriter< ScalarField > >(cpFieldId, "heat_capacity");
      vtkOutput->addCellDataWriter(cpWriter);

      auto alphaWriter = make_shared< field::VTKWriter< ScalarField > >(alphaFieldId, "thermal_diffusivity");
      vtkOutput->addCellDataWriter(alphaWriter);

      timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   // lbm_generated::PerformanceEvaluation<BlockDataID> const performance(blocks, uFieldId, uFieldId);  //TODO

   WcTimingPool timeloopTiming;
   WcTimer simTimer;
   simTimer.start();
   timeloop.run(timeloopTiming);
   // WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
   simTimer.end();
   double simTime = simTimer.max();
   WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); }
   // performance.logResultOnRoot(timeSteps, simTime);
   const auto reducedTimeloopTiming = timeloopTiming.getReduced();
   WALBERLA_LOG_RESULT_ON_ROOT("Time loop timing:\n" << *reducedTimeloopTiming)
   WALBERLA_LOG_RESULT_ON_ROOT("Elapsed time:\t" << simTime)

   uint_t mlups = timeSteps * xCells * yCells / (uint_c(simTime * 1000000.0));
   WALBERLA_LOG_RESULT_ON_ROOT("mlups:\t" << mlups)

   return EXIT_SUCCESS;
}
}

int main(int argc, char** argv) { walberla::main(argc, argv); }
