//======================================================================================================================
//! \file CouetteFlowSimulation.cpp - FINAL WORKING VERSION WITH DIRECT BOUNDARY CONDITIONS
//! \author Rahil Doshi  
//! \brief Couette Flow Simulation with temperature-dependent viscosity using waLBerla and materforge
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

// Include generated headers
#include "gen/CouetteFlowKernel.hpp"

namespace walberla {
namespace couette_flow {

using namespace walberla;

// Type definitions
typedef GhostLayerField<real_t, 1> ScalarField;
typedef GhostLayerField<real_t, 3> VectorField;
typedef GhostLayerField<real_t, 19> PdfField;

// Initialize fields
void initializeFields(const shared_ptr<StructuredBlockForest>& blocks,
                     BlockDataID densityId, BlockDataID forceId, BlockDataID pdfFieldId,
                     BlockDataID temperatureId, BlockDataID velocityId, BlockDataID viscosityId,
                     real_t hotTemperature, real_t coldTemperature, real_t wallVelocity,
                     real_t refTemp, real_t refViscosity) {
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        // Get field pointers
        VectorField* velocity = block->getData<VectorField>(velocityId);
        ScalarField* temperature = block->getData<ScalarField>(temperatureId);
        ScalarField* density = block->getData<ScalarField>(densityId);
        VectorField* force = block->getData<VectorField>(forceId);
        ScalarField* viscosity = block->getData<ScalarField>(viscosityId);

        // Get domain size
        const auto& domain = blocks->getDomain();
        const real_t domainHeight = domain.ySize();
        
        // Initialize fields
        for (auto cell = velocity->beginXYZ(); cell != velocity->end(); ++cell) {
            Cell localCell = cell.cell();
            Cell globalCell = localCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, *block);
            
            // Get normalized position
            Vector3<real_t> cellCenter = blocks->getCellCenter(globalCell);
            real_t y_normalized = std::max(0.0, std::min(1.0, 
                (cellCenter[1] - domain.yMin()) / domainHeight));
            
            // Linear temperature profile
            real_t temp = coldTemperature + (hotTemperature - coldTemperature) * y_normalized;
            temperature->get(localCell) = temp;
            
            // INITIAL COUETTE PROFILE for faster convergence
            real_t u_x = wallVelocity * y_normalized;
            velocity->get(localCell, 0) = u_x;
            velocity->get(localCell, 1) = 0.0;
            velocity->get(localCell, 2) = 0.0;

            // LBM normalized density
            density->get(localCell) = 1.0;

            // Zero force
            force->get(localCell, 0) = 0.0;
            force->get(localCell, 1) = 0.0;
            force->get(localCell, 2) = 0.0;

            // Viscosity
            viscosity->get(localCell) = refViscosity;
        }

        // Initialize PDFs
        CouetteFlowInit init(densityId, forceId, pdfFieldId, velocityId);
        init(&*block);

        WALBERLA_LOG_DEVEL_ON_ROOT("Initialized block with Couette profile and temperature range: " 
                                  << coldTemperature << " - " << hotTemperature << " K");
    }
}

// Main simulation function
int main(int argc, char** argv) {
    mpi::Environment env(argc, argv);
    
    // Simulation parameters
    const uint_t timesteps = 10000;
    const uint_t vtkWriteFrequency = 500;
    
    // Physical parameters - ADJUSTED for LBM units
    const real_t wallVelocity = 0.05;       // Reduced for stability
    const real_t hotTemperature = 600.0;    // K
    const real_t coldTemperature = 300.0;   // K
    
    // Domain size
    const uint_t xBlocks = 1, yBlocks = 1, zBlocks = 1;
    const uint_t xCells = 50, yCells = 25, zCells = 10;
    const real_t xSize = 1e-3;  // 1mm
    const real_t ySize = 0.5e-3; // 0.5mm
    const real_t zSize = 0.2e-3; // 0.2mm
    
    // Material properties
    const real_t refViscosity = 0.001;
    const real_t refTemperature = 300.0;
    
    WALBERLA_LOG_INFO_ON_ROOT("=== Couette Flow Simulation with Temperature-Dependent Viscosity ===");
    WALBERLA_LOG_INFO_ON_ROOT("Domain: " << xSize*1000 << " x " << ySize*1000 << " x " << zSize*1000 << " mm");
    WALBERLA_LOG_INFO_ON_ROOT("Grid: " << xCells << " x " << yCells << " x " << zCells << " cells");
    WALBERLA_LOG_INFO_ON_ROOT("Wall velocity: " << wallVelocity << " m/s");
    WALBERLA_LOG_INFO_ON_ROOT("Temperature range: " << coldTemperature << " - " << hotTemperature << " K");
    
    // Create block forest
    auto aabb = math::AABB(0, 0, 0, xSize, ySize, zSize);
    shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGrid(
        aabb, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, 
        true, false, true, false // periodic in x,z; walls in y
    );
    
    // Add fields
    BlockDataID densityId = field::addToStorage<ScalarField>(blocks, "density", 1.0, field::fzyx, 1);
    BlockDataID forceId = field::addToStorage<VectorField>(blocks, "force", real_t(0.0), field::fzyx, 1);
    BlockDataID pdfFieldId = field::addToStorage<PdfField>(blocks, "pdfs", real_t(0.0), field::fzyx, 1);
    BlockDataID temperatureId = field::addToStorage<ScalarField>(blocks, "temperature", refTemperature, field::fzyx, 1);
    BlockDataID velocityId = field::addToStorage<VectorField>(blocks, "velocity", real_t(0.0), field::fzyx, 1);
    BlockDataID viscosityId = field::addToStorage<ScalarField>(blocks, "viscosity", refViscosity, field::fzyx, 1);
    
    // Initialize fields
    initializeFields(blocks, densityId, forceId, pdfFieldId, temperatureId, velocityId, viscosityId,
                    hotTemperature, coldTemperature, wallVelocity, refTemperature, refViscosity);
    
    // Communication schemes
    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commVector(blocks);
    commVector.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(velocityId));
    commVector.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(forceId));
    
    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commScalar(blocks);
    commScalar.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(temperatureId));
    commScalar.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(viscosityId));
    commScalar.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(densityId));
    
    // Create sweep objects
    auto couetteFlowSweep = make_shared<CouetteFlowSweep>(densityId, forceId, pdfFieldId, temperatureId, velocityId, viscosityId);
    
    // Time loop - CORRECTED: Direct boundary condition implementation
    SweepTimeloop timeloop(blocks, timesteps);
    
    timeloop.add() 
        << BeforeFunction(commVector, "Vector Communication")
        << BeforeFunction(commScalar, "Scalar Communication")
        << Sweep([couetteFlowSweep](IBlock* block) { 
            (*couetteFlowSweep)(block); 
        }, "CouetteFlowSweep")
        << AfterFunction([&]() {
            // DIRECT boundary condition implementation - CORRECTED
            for (auto block = blocks->begin(); block != blocks->end(); ++block) {
                VectorField* velocity = block->getData<VectorField>(velocityId);
                ScalarField* temperature = block->getData<ScalarField>(temperatureId);
                
                // Get the size of the velocity field (including ghost layers)
                CellInterval fieldSize = velocity->xyzSizeWithGhostLayer();
                
                // Apply boundary conditions directly to ALL cells in boundary layers
                for (auto cell = fieldSize.begin(); cell != fieldSize.end(); ++cell) {
                    Cell localCell = *cell;
                    
                    // Top wall boundary (y = yMax) - moving wall
                    if (blocks->atDomainYMaxBorder(*block) && localCell.y() >= fieldSize.yMax() - 1) {
                        velocity->get(localCell, 0) = wallVelocity;  // x-component = wall velocity
                        velocity->get(localCell, 1) = 0.0;          // y-component = 0
                        velocity->get(localCell, 2) = 0.0;          // z-component = 0
                        temperature->get(localCell) = hotTemperature;
                        
                        WALBERLA_LOG_DEVEL_VAR_ON_ROOT(localCell);
                        WALBERLA_LOG_DEVEL_VAR_ON_ROOT(velocity->get(localCell, 0));
                    }
                    
                    // Bottom wall boundary (y = yMin) - stationary wall
                    if (blocks->atDomainYMinBorder(*block) && localCell.y() <= fieldSize.yMin() + 1) {
                        velocity->get(localCell, 0) = 0.0;          // x-component = 0
                        velocity->get(localCell, 1) = 0.0;          // y-component = 0
                        velocity->get(localCell, 2) = 0.0;          // z-component = 0
                        temperature->get(localCell) = coldTemperature;
                    }
                }
            }
        }, "Direct Boundary Conditions");
    
    // VTK output
    if (vtkWriteFrequency > 0) {
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "CouetteFlow_TempVisc", vtkWriteFrequency, 0, 
                                                      false, "vtk_out", "simulation_step", false, true, true, false, 0);
        
        auto velocityWriter = make_shared<field::VTKWriter<VectorField>>(velocityId, "velocity");
        auto temperatureWriter = make_shared<field::VTKWriter<ScalarField>>(temperatureId, "temperature");
        auto viscosityWriter = make_shared<field::VTKWriter<ScalarField>>(viscosityId, "viscosity");
        auto densityWriter = make_shared<field::VTKWriter<ScalarField>>(densityId, "density");
        
        vtkOutput->addCellDataWriter(velocityWriter);
        vtkOutput->addCellDataWriter(temperatureWriter);
        vtkOutput->addCellDataWriter(viscosityWriter);
        vtkOutput->addCellDataWriter(densityWriter);
        
        timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
        
        // Write initial state
        vtk::writeFiles(vtkOutput)();
    }
    
    // Run simulation
    WcTimingPool timeloopTiming;
    WcTimer simTimer;
    
    WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with " << timesteps << " time steps...");
    simTimer.start();
    timeloop.run(timeloopTiming);
    simTimer.end();
    
    // Performance results
    real_t simTime = real_c(simTimer.max());
    WALBERLA_MPI_SECTION() { 
        walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); 
    }
    
    const auto reducedTimeloopTiming = timeloopTiming.getReduced();
    
    WALBERLA_LOG_RESULT_ON_ROOT("=== Simulation Completed Successfully ===");
    WALBERLA_LOG_RESULT_ON_ROOT("Total simulation time: " << simTime << " seconds");
    WALBERLA_LOG_RESULT_ON_ROOT("Time per timestep: " << (simTime / timesteps * 1000) << " ms");
    WALBERLA_LOG_RESULT_ON_ROOT("Performance timing:\n" << *reducedTimeloopTiming);
    
    // Sample final results
    WALBERLA_ROOT_SECTION() {
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            ScalarField* temperature = block->getData<ScalarField>(temperatureId);
            ScalarField* viscosity = block->getData<ScalarField>(viscosityId);
            VectorField* velocity = block->getData<VectorField>(velocityId);
            
            Cell middleCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells/2), static_cast<cell_idx_t>(zCells/2));
            Cell topCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells-1), static_cast<cell_idx_t>(zCells/2));
            Cell bottomCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(0), static_cast<cell_idx_t>(zCells/2));
            
            real_t finalTempMid = temperature->get(middleCell);
            real_t finalViscMid = viscosity->get(middleCell);
            real_t finalVelXMid = velocity->get(middleCell, 0);
            real_t finalVelXTop = velocity->get(topCell, 0);
            real_t finalVelXBot = velocity->get(bottomCell, 0);
            
            WALBERLA_LOG_RESULT_ON_ROOT("Final center temperature: " << finalTempMid << " K");
            WALBERLA_LOG_RESULT_ON_ROOT("Final center viscosity: " << finalViscMid*1000 << " mPa·s");
            WALBERLA_LOG_RESULT_ON_ROOT("Final center velocity (x): " << finalVelXMid << " m/s");
            WALBERLA_LOG_RESULT_ON_ROOT("Final top velocity (x): " << finalVelXTop << " m/s");
            WALBERLA_LOG_RESULT_ON_ROOT("Final bottom velocity (x): " << finalVelXBot << " m/s");
            break;
        }
    }
    
    return EXIT_SUCCESS;
}

} // namespace couette_flow
} // namespace walberla

int main(int argc, char** argv) { 
    return walberla::couette_flow::main(argc, argv); 
}
