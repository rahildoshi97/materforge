//======================================================================================================================
//
//! \file CodegenHeatEquationWithMaterialCPU.cpp
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

#include "gen/HeatEquationKernelWithMaterialCPU.hpp"

namespace walberla {

typedef GhostLayerField<real_t, 1> ScalarField;

void swapFields(StructuredBlockForest& blocks, BlockDataID uID, BlockDataID uTmpID)
{
    for (auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        ScalarField* u = block->getData<ScalarField>(uID);
        ScalarField* u_tmp = block->getData<ScalarField>(uTmpID);
        u->swapDataPointers(u_tmp);
    }
}

void initDirichletBoundaryNorth(const shared_ptr<StructuredBlockForest>& blocks,
                               BlockDataID uId, BlockDataID uTmpId)
{
    for (auto block = blocks->begin(); block != blocks->end(); ++block)
    {
        if (blocks->atDomainYMaxBorder(*block))
        {
            ScalarField* u = block->getData<ScalarField>(uId);
            ScalarField* u_tmp = block->getData<ScalarField>(uTmpId);
            CellInterval xyz = u->xyzSizeWithGhostLayer();
            xyz.yMin() = xyz.yMax();

            for (auto cell = xyz.begin(); cell != xyz.end(); ++cell)
            {
                const Vector3<real_t> p = blocks->getBlockLocalCellCenter(*block, *cell);
                // Set constant temperature boundary condition
                real_t v = real_c(3800.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }
    }
}

void initDirichletBoundariesAllSides(const shared_ptr<StructuredBlockForest>& blocks,
                                    BlockDataID uId, BlockDataID uTmpId)
{
    for (auto block = blocks->begin(); block != blocks->end(); ++block)
    {
        ScalarField* u = block->getData<ScalarField>(uId);
        ScalarField* u_tmp = block->getData<ScalarField>(uTmpId);

        // North boundary (Y max)
        if (blocks->atDomainYMaxBorder(*block))
        {
            CellInterval north = u->xyzSizeWithGhostLayer();
            north.yMin() = north.yMax();

            for (auto cell = north.begin(); cell != north.end(); ++cell)
            {
                real_t v = real_c(3800.0); // Hot boundary
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }

        // South boundary (Y min)
        if (blocks->atDomainYMinBorder(*block))
        {
            CellInterval south = u->xyzSizeWithGhostLayer();
            south.yMax() = south.yMin();

            for (auto cell = south.begin(); cell != south.end(); ++cell)
            {
                real_t v = real_c(300.0); // Cold boundary
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }

        // East boundary (X max)
        if (blocks->atDomainXMaxBorder(*block))
        {
            CellInterval east = u->xyzSizeWithGhostLayer();
            east.xMin() = east.xMax();

            for (auto cell = east.begin(); cell != east.end(); ++cell)
            {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }

        // West boundary (X min)
        if (blocks->atDomainXMinBorder(*block))
        {
            CellInterval west = u->xyzSizeWithGhostLayer();
            west.xMax() = west.xMin();

            for (auto cell = west.begin(); cell != west.end(); ++cell)
            {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }

        // Top boundary (Z max)
        if (blocks->atDomainZMaxBorder(*block))
        {
            CellInterval top = u->xyzSizeWithGhostLayer();
            top.zMin() = top.zMax();

            for (auto cell = top.begin(); cell != top.end(); ++cell)
            {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
                u_tmp->get(*cell) = v;
            }
        }

        // Bottom boundary (Z min)
        if (blocks->atDomainZMinBorder(*block))
        {
            CellInterval bottom = u->xyzSizeWithGhostLayer();
            bottom.zMax() = bottom.zMin();

            for (auto cell = bottom.begin(); cell != bottom.end(); ++cell)
            {
                real_t v = real_c(300.0);
                u->get(*cell) = v;
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

    // Ensure matching aspect ratios of cells and domain.
    constexpr uint_t x = 128;
    const uint_t xCells = uint_c(x);
    const uint_t yCells = uint_c(x);
    const uint_t zCells = uint_c(x);
    const real_t xSize = real_c(1.0);
    const real_t ySize = real_c(1.0);  // 2
    const real_t zSize = real_c(1.0);
    const uint_t xBlocks = uint_c(1);
    const uint_t yBlocks = uint_c(1);  // 2
    const uint_t zBlocks = uint_c(1);

    const uint_t processes = uint_c(MPIManager::instance()->numProcesses());
    if (processes != xBlocks * yBlocks * zBlocks)
    { WALBERLA_ABORT("The number of processes must be equal to the number of blocks!"); }

    const real_t dx = xSize / real_c(xBlocks * xCells + uint_c(1));
    const real_t dy = ySize / real_c(yBlocks * yCells + uint_c(1));
    const real_t dz = zSize / real_c(zBlocks * zCells + uint_c(1));

    WALBERLA_CHECK_FLOAT_EQUAL(dx, dy);

    const real_t dt = real_c(1);
    uint_t timeSteps = uint_c(2e4);
    uint_t vtkWriteFrequency = uint_c(400);

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

    BlockDataID uFieldId = field::addToStorage<ScalarField>(blocks, "u", real_c(300.0), field::fzyx, uint_c(1));
    BlockDataID uTmpFieldId = field::addToStorage<ScalarField>(blocks, "u_tmp", real_c(300.0), field::fzyx, uint_c(1));
    BlockDataID alphaFieldId = field::addToStorage<ScalarField>(blocks, "thermal_diffusivity", real_c(0.0), field::fzyx, uint_c(1));

    /////////////////////
    /// COMMUNICATION ///
    /////////////////////

    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commScheme(blocks);
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(uFieldId));

    //////////////////////////
    /// DIRICHLET BOUNDARY ///
    //////////////////////////

    //initDirichletBoundaryNorth(blocks, uFieldId, uTmpFieldId);  // Neumann
    initDirichletBoundariesAllSides(blocks, uFieldId, uTmpFieldId);  // Dirichlet

    ////////////////////////
    /// NEUMANN BOUNDARY ///
    ////////////////////////

    //pde::NeumannDomainBoundary<ScalarField> neumann(*blocks, uFieldId);  // Neumann
    //neumann.excludeBoundary(stencil::N);  // Neumann

    ////////////////
    /// TIMELOOP ///
    ////////////////

    SweepTimeloop timeloop(blocks, timeSteps);

    timeloop.add() << BeforeFunction(commScheme, "Communication")
                   //<< BeforeFunction(neumann, "Neumann Boundaries")  // Neumann
                   //<< BeforeFunction([&]() {initDirichletBoundariesAllSides(blocks, uFieldId, uTmpFieldId);}, "Dirichlet Boundaries")
                   << Sweep(HeatEquationKernelWithMaterialCPU(alphaFieldId, uFieldId, uTmpFieldId, dt, dx), "HeatEquationKernelWithMaterialCPU")
                   << AfterFunction([blocks, uFieldId, uTmpFieldId]() {swapFields(*blocks, uFieldId, uTmpFieldId);}, "Swap");

    if (vtkWriteFrequency > 0)
    {
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtkCPU", vtkWriteFrequency, 0, false, "vtk_out_cpu",
                                                       "simulation_step", false, true, true, false, 0);

        auto tempWriter = make_shared<field::VTKWriter<ScalarField>>(uFieldId, "temperature");
        vtkOutput->addCellDataWriter(tempWriter);

        auto alphaWriter = make_shared<field::VTKWriter<ScalarField>>(alphaFieldId, "thermal_diffusivity");
        vtkOutput->addCellDataWriter(alphaWriter);

        timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output CPU");
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
        timeloop.setCurrentTimeStepToZero();
        WcTimingPool timeloopTiming;
        WcTimer simTimer;

        WALBERLA_MPI_WORLD_BARRIER()

        WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timeSteps << " time steps")
        simTimer.start();
        timeloop.run(timeloopTiming);
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
