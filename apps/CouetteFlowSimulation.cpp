//======================================================================================================================
//! \file CouetteFlowThermalSimulation.cpp
//! \author Rahil Doshi
//! \brief Thermal Couette flow simulation using LBM with direct velocity BCs
//! \details This simulation models thermal Couette flow between two parallel plates
//!          using the Lattice Boltzmann Method (LBM). The top plate moves at
//!          a constant velocity while the bottom plate is stationary. A linear
//!          temperature gradient is imposed between the plates. The simulation
//!          uses direct velocity boundary conditions at the walls and a small
//!          pressure gradient force to sustain the flow.
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
namespace couette_flow_thermal {

using namespace walberla;

// Type definitions
typedef GhostLayerField<real_t, 1> ScalarField;
typedef GhostLayerField<real_t, 3> VectorField;
typedef GhostLayerField<real_t, 19> PdfField;

// Initialize fields with thermal Couette profile
void initializeFields(const shared_ptr<StructuredBlockForest>& blocks,
                     BlockDataID densityId, BlockDataID forceId, BlockDataID pdfFieldId,
                     BlockDataID temperatureId, BlockDataID velocityId, BlockDataID viscosityId,
                     real_t hotTemperature, real_t coldTemperature, real_t wallVelocity,
                     uint_t yCells, real_t refViscosity) {
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        VectorField* velocity = block->getData<VectorField>(velocityId);
        ScalarField* temperature = block->getData<ScalarField>(temperatureId);
        ScalarField* density = block->getData<ScalarField>(densityId);
        VectorField* force = block->getData<VectorField>(forceId);
        ScalarField* viscosity = block->getData<ScalarField>(viscosityId);
        
        for (auto cell = velocity->beginXYZ(); cell != velocity->end(); ++cell) {
            Cell localCell = cell.cell();
            Cell globalCell = localCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, *block);
            
            real_t y_normalized = real_t(globalCell.y()) / real_t(yCells - 1);
            y_normalized = std::max(0.0, std::min(1.0, y_normalized));
            
            // Linear temperature profile
            real_t temp = coldTemperature + (hotTemperature - coldTemperature) * y_normalized;
            temperature->get(localCell) = temp;
            
            // Initialize linear velocity profile
            velocity->get(localCell, 0) = wallVelocity * y_normalized;
            velocity->get(localCell, 1) = 0.0;
            velocity->get(localCell, 2) = 0.0;

            density->get(localCell) = 1.0;
            
            // Initialize force
            force->get(localCell, 0) = 0.0;
            force->get(localCell, 1) = 0.0;
            force->get(localCell, 2) = 0.0;

            viscosity->get(localCell) = refViscosity;
        }

        // Initialize PDFs
        CouetteFlowInit init(densityId, forceId, pdfFieldId, velocityId);
        init(&*block);
    }
    
    WALBERLA_LOG_INFO_ON_ROOT("Field initialization complete");
}

// Boundary condition application
void applySimpleCouetteBoundaries(const shared_ptr<StructuredBlockForest>& blocks,
                                  BlockDataID velocityId, real_t wallVelocity) {
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        VectorField* velocity = block->getData<VectorField>(velocityId);
        CellInterval domain = velocity->xyzSize();
        
        // TOP WALL (Moving wall)
        if (blocks->atDomainYMaxBorder(*block)) {
            int topY = domain.yMax();
            
            for (int x = domain.xMin(); x <= domain.xMax(); ++x) {
                for (int z = domain.zMin(); z <= domain.zMax(); ++z) {
                    Cell topCell(x, topY, z);
                    velocity->get(topCell, 0) = wallVelocity;
                    velocity->get(topCell, 1) = 0.0;
                    velocity->get(topCell, 2) = 0.0;
                }
            }
        }
        
        // BOTTOM WALL (Stationary wall)
        if (blocks->atDomainYMinBorder(*block)) {
            int bottomY = domain.yMin();
            
            for (int x = domain.xMin(); x <= domain.xMax(); ++x) {
                for (int z = domain.zMin(); z <= domain.zMax(); ++z) {
                    Cell bottomCell(x, bottomY, z);
                    velocity->get(bottomCell, 0) = 0.0;
                    velocity->get(bottomCell, 1) = 0.0;
                    velocity->get(bottomCell, 2) = 0.0;
                }
            }
        }
    }
}

/* ============================================================================
 * This function applies the external driving force required for LBM Couette
 * flow simulation. Without this force, direct velocity boundary conditions
 * create unphysical negative velocities in the interior domain.
 * ============================================================================ */
// Grid-independent pressure gradient force
void applyPressureGradientForce(const shared_ptr<StructuredBlockForest>& blocks,
                               BlockDataID forceId, real_t pressureGradient) {
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        VectorField* force = block->getData<VectorField>(forceId);
        
        for (auto cell = force->beginXYZ(); cell != force->end(); ++cell) {
            Cell localCell = cell.cell();
            
            // Constant pressure gradient force in flow direction (x) (grid-independent)
            // This provides the energy needed to sustain Couette flow
            // against LBM's natural tendency toward zero-velocity equilibrium
            force->get(localCell, 0) = pressureGradient;  // Force per unit mass F_x = -dP/dx / rho
            force->get(localCell, 1) = 0.0;               // No Y-force (wall normal)
            force->get(localCell, 2) = 0.0;               // No Z-force (spanwise)
        }
    }
}

int main(int argc, char** argv) {
    mpi::Environment env(argc, argv);
    
    // Parameters
    const uint_t timesteps = 5000;
    const uint_t vtkWriteFrequency = 50;
    const real_t wallVelocity = 0.05;
    const real_t hotTemperature = 600.0;
    const real_t coldTemperature = 300.0;
    const uint_t xBlocks = 1, yBlocks = 1, zBlocks = 1;
    const uint_t xCells = 128, yCells = 32, zCells = 32;
    const real_t refViscosity = 0.0925;
    const real_t pressureGradient = 0.0001;  // Small driving force

    WALBERLA_LOG_INFO_ON_ROOT("=== THERMAL COUETTE FLOW SIMULATION ===");
    WALBERLA_LOG_INFO_ON_ROOT("Grid: " << xCells << " x " << yCells << " x " << zCells);
    WALBERLA_LOG_INFO_ON_ROOT("Wall velocity: " << wallVelocity);
    WALBERLA_LOG_INFO_ON_ROOT("Temperature range: " << coldTemperature << " - " << hotTemperature);
    WALBERLA_LOG_INFO_ON_ROOT("Reference viscosity: " << refViscosity);
    WALBERLA_LOG_INFO_ON_ROOT("Pressure gradient: " << pressureGradient);
    WALBERLA_LOG_INFO_ON_ROOT("Timesteps: " << timesteps);

    // Create block forest
    auto aabb = math::AABB(0, 0, 0, xBlocks*xCells, yBlocks*yCells, zBlocks*zCells);
    shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGrid(
        aabb, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, 
        true,   // oneBlockPerProcess  
        true,   // xPeriodic = true (periodic in x)
        false,  // yPeriodic = false (walls in y)
        true,   // zPeriodic = true (periodic in z)
        false   // keepGlobalBlockInformation = false
    );

    // Add fields
    BlockDataID densityId = field::addToStorage<ScalarField>(blocks, "density", 1.0, field::fzyx, 1);
    BlockDataID forceId = field::addToStorage<VectorField>(blocks, "force", real_t(0.0), field::fzyx, 1);
    BlockDataID pdfFieldId = field::addToStorage<PdfField>(blocks, "pdfs", real_t(0.0), field::fzyx, 1);
    BlockDataID temperatureId = field::addToStorage<ScalarField>(blocks, "temperature", 300.0, field::fzyx, 1);
    BlockDataID velocityId = field::addToStorage<VectorField>(blocks, "velocity", real_t(0.0), field::fzyx, 1);
    BlockDataID viscosityId = field::addToStorage<ScalarField>(blocks, "viscosity", refViscosity, field::fzyx, 1);
    
    // Initialize fields
    initializeFields(blocks, densityId, forceId, pdfFieldId, temperatureId, 
                    velocityId, viscosityId, hotTemperature, coldTemperature, wallVelocity, yCells, refViscosity);
    
    // Apply driving force
    // applyPressureGradientForce(blocks, forceId, pressureGradient);

    // Communication scheme
    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commScheme(blocks);
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(velocityId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(temperatureId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(viscosityId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(densityId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<PdfField>>(pdfFieldId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(forceId));

    // Create sweeps
    auto couetteFlowSweep = make_shared<CouetteFlowSweep>(
        densityId, forceId, pdfFieldId, temperatureId, velocityId, viscosityId);

    // TIME LOOP
    SweepTimeloop timeloop(blocks, timesteps);
    
    timeloop.add() 
        << BeforeFunction(commScheme, "Communication")
        << Sweep([couetteFlowSweep](IBlock* block) { 
            (*couetteFlowSweep)(block);
        }, "CouetteFlowSweep")
        << AfterFunction([&]() {
            applySimpleCouetteBoundaries(blocks, velocityId, wallVelocity);
        }, "Boundary Conditions");
    
    // VTK output
    if (vtkWriteFrequency > 0) {
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "ThermalCouetteFlow", 
                                                      vtkWriteFrequency, 0, false, "vtk_out", 
                                                      "simulation_step", false, true, true, false, 0);
        
        auto velocityWriter = make_shared<field::VTKWriter<VectorField>>(velocityId, "velocity");
        auto temperatureWriter = make_shared<field::VTKWriter<ScalarField>>(temperatureId, "temperature");
        auto viscosityWriter = make_shared<field::VTKWriter<ScalarField>>(viscosityId, "viscosity");
        auto densityWriter = make_shared<field::VTKWriter<ScalarField>>(densityId, "density");
        // auto forceWriter = make_shared<field::VTKWriter<VectorField>>(forceId, "force");
        
        vtkOutput->addCellDataWriter(velocityWriter);
        vtkOutput->addCellDataWriter(temperatureWriter);
        vtkOutput->addCellDataWriter(viscosityWriter);
        vtkOutput->addCellDataWriter(densityWriter);
        // vtkOutput->addCellDataWriter(forceWriter);
        
        timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
        vtk::writeFiles(vtkOutput)();
    }
    
    // Run simulation
    WcTimingPool timeloopTiming;
    WcTimer simTimer;

    WALBERLA_LOG_INFO_ON_ROOT("Starting thermal Couette simulation with " << timesteps << " timesteps...");
    simTimer.start();
    timeloop.run(timeloopTiming);
    simTimer.end();
    
    // Results
    real_t simTime = real_c(simTimer.max());
    WALBERLA_MPI_SECTION() { 
        walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); 
    }
    
    const auto reducedTimeloopTiming = timeloopTiming.getReduced();
    
    WALBERLA_LOG_RESULT_ON_ROOT("=== Thermal Couette Results ===");
    WALBERLA_LOG_RESULT_ON_ROOT("Total time: " << simTime << " seconds");
    WALBERLA_LOG_RESULT_ON_ROOT("Time per timestep: " << (simTime / timesteps * 1000) << " ms");
    WALBERLA_LOG_RESULT_ON_ROOT("Performance timing:\n" << *reducedTimeloopTiming);
    
    // Final velocity analysis
    WALBERLA_ROOT_SECTION() {
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            ScalarField* temperature = block->getData<ScalarField>(temperatureId);
            ScalarField* viscosity = block->getData<ScalarField>(viscosityId);
            VectorField* velocity = block->getData<VectorField>(velocityId);
            
            Cell quarterCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells/4), static_cast<cell_idx_t>(zCells/2));
            Cell middleCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells/2), static_cast<cell_idx_t>(zCells/2));
            Cell threeQuarterCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(3*yCells/4), static_cast<cell_idx_t>(zCells/2));
            Cell topCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells-1), static_cast<cell_idx_t>(zCells/2));
            Cell bottomCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(0), static_cast<cell_idx_t>(zCells/2));
            
            WALBERLA_LOG_RESULT_ON_ROOT("=== FINAL VELOCITY PROFILE ===");
            WALBERLA_LOG_RESULT_ON_ROOT("Bottom (y=0): " << velocity->get(bottomCell, 0) << " (expected: 0.000)");
            WALBERLA_LOG_RESULT_ON_ROOT("Quarter (y=1/4): " << velocity->get(quarterCell, 0) << " (expected: " << wallVelocity*0.25 << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("Middle (y=1/2): " << velocity->get(middleCell, 0) << " (expected: " << wallVelocity*0.5 << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("3/4 (y=3/4): " << velocity->get(threeQuarterCell, 0) << " (expected: " << wallVelocity*0.75 << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("Top (y=1): " << velocity->get(topCell, 0) << " (expected: " << wallVelocity << ")");
            
            // Viscosity validation
            WALBERLA_LOG_RESULT_ON_ROOT("=== VISCOSITY PROFILE ===");
            real_t bottomVisc = viscosity->get(bottomCell);
            real_t middleVisc = viscosity->get(middleCell);
            real_t topVisc = viscosity->get(topCell);
            real_t bottomTemp = temperature->get(bottomCell);
            real_t middleTemp = temperature->get(middleCell);
            real_t topTemp = temperature->get(topCell);
            
            // Expected from materforge: 0.1 * exp(-0.0005 * (T - 300))
            real_t expectedBottomVisc = 0.1 * std::exp(-0.0005 * (bottomTemp - 300.0));
            real_t expectedMiddleVisc = 0.1 * std::exp(-0.0005 * (middleTemp - 300.0));
            real_t expectedTopVisc = 0.1 * std::exp(-0.0005 * (topTemp - 300.0));
            
            WALBERLA_LOG_RESULT_ON_ROOT("Bottom: visc=" << bottomVisc << " (expected=" << expectedBottomVisc << ", temp=" << bottomTemp << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("Middle: visc=" << middleVisc << " (expected=" << expectedMiddleVisc << ", temp=" << middleTemp << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("Top: visc=" << topVisc << " (expected=" << expectedTopVisc << ", temp=" << topTemp << ")");
            
            // Validation
            real_t viscError = std::abs(middleVisc - expectedMiddleVisc) / expectedMiddleVisc * 100.0;
            if (viscError < 1.0) {
                WALBERLA_LOG_RESULT_ON_ROOT("VISCOSITY VALIDATION: EXCELLENT (error < 1%)");
            } else {
                WALBERLA_LOG_RESULT_ON_ROOT("*** VISCOSITY VALIDATION: ERROR = " << viscError << "% ***");
            }

            break;
        }
    }
    
    return EXIT_SUCCESS;
}

} // namespace couette_flow_thermal
} // namespace walberla

int main(int argc, char** argv) { 
    return walberla::couette_flow_thermal::main(argc, argv); 
}
