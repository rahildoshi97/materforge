//======================================================================================================================
//! \file CouetteFlowThermalSimulation.cpp - USING GENERATED BOUNDARIES
//! \author Rahil Doshi
//! \brief Thermal Couette flow using generated boundary sweeps
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
            
            // SIMPLE initial velocity - let LBM develop the profile
            velocity->get(localCell, 0) = 0.0;  // Start from rest
            velocity->get(localCell, 1) = 0.0;
            velocity->get(localCell, 2) = 0.0;

            density->get(localCell) = 1.0;
            
            force->get(localCell, 0) = 0.0;
            force->get(localCell, 1) = 0.0;
            force->get(localCell, 2) = 0.0;

            viscosity->get(localCell) = refViscosity;
        }

        // Initialize PDFs
        CouetteFlowInit init(densityId, forceId, pdfFieldId, velocityId);
        init(&*block);
    }
}

// Function to apply boundary conditions
void applyBoundaryConditions(const shared_ptr<StructuredBlockForest>& blocks,
                            const shared_ptr<TopWallBC>& topWallBC,
                            const shared_ptr<BottomWallBC>& bottomWallBC) {
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        // Apply top wall boundary
        if (blocks->atDomainYMaxBorder(*block)) {
            // Get the boundary region - just the top layer
            VectorField* velocity = block->getData<VectorField>(topWallBC->f_velocityId());
            CellInterval domain = velocity->xyzSize();
            
            CellInterval topRegion(domain.xMin(), domain.yMax(), domain.zMin(),
                                  domain.xMax(), domain.yMax(), domain.zMax());
            
            topWallBC->runOnCellInterval(&(*block), topRegion);
        }
        
        // Apply bottom wall boundary
        if (blocks->atDomainYMinBorder(*block)) {
            // Get the boundary region - just the bottom layer
            VectorField* velocity = block->getData<VectorField>(bottomWallBC->f_velocityId());
            CellInterval domain = velocity->xyzSize();
            
            CellInterval bottomRegion(domain.xMin(), domain.yMin(), domain.zMin(),
                                     domain.xMax(), domain.yMin(), domain.zMax());
            
            bottomWallBC->runOnCellInterval(&(*block), bottomRegion);
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
    const uint_t timesteps = 10000;
    const uint_t vtkWriteFrequency = 50;
    const real_t wallVelocity = 0.02;
    const real_t hotTemperature = 600.0;
    const real_t coldTemperature = 300.0;
    const uint_t xBlocks = 1, yBlocks = 1, zBlocks = 1;
    const uint_t xCells = 128, yCells = 64, zCells = 32;
    const real_t refViscosity = 0.05;  // REDUCED viscosity for stability
    // Grid-independent pressure gradient
    const real_t pressureGradient = 0.0005;  // Force per unit mass (grid-independent)

    WALBERLA_LOG_INFO_ON_ROOT("=== Thermal Couette Flow Simulation ===");
    WALBERLA_LOG_INFO_ON_ROOT("Grid: " << xCells << " x " << yCells << " x " << zCells);
    WALBERLA_LOG_INFO_ON_ROOT("Wall velocity: " << wallVelocity);
    WALBERLA_LOG_INFO_ON_ROOT("Temperature range: " << coldTemperature << " - " << hotTemperature);
    WALBERLA_LOG_INFO_ON_ROOT("Reference viscosity: " << refViscosity);
    WALBERLA_LOG_INFO_ON_ROOT("Pressure gradient: " << pressureGradient);

    // Create block forest - CORRECTED periodicity
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
    BlockDataID pdfFieldTmpId = field::addToStorage<PdfField>(blocks, "pdfs_tmp", real_t(0.0), field::fzyx, 1);
    BlockDataID temperatureId = field::addToStorage<ScalarField>(blocks, "temperature", 300.0, field::fzyx, 1);
    BlockDataID velocityId = field::addToStorage<VectorField>(blocks, "velocity", real_t(0.0), field::fzyx, 1);
    BlockDataID viscosityId = field::addToStorage<ScalarField>(blocks, "viscosity", refViscosity, field::fzyx, 1);
    
    // Initialize fields
    initializeFields(blocks, densityId, forceId, pdfFieldId, temperatureId, 
                    velocityId, viscosityId, hotTemperature, coldTemperature, wallVelocity, yCells, refViscosity);
    
    // Communication
    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commVector(blocks);
    commVector.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(velocityId));
    commVector.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(forceId));
    
    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commScalar(blocks);
    commScalar.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(temperatureId));
    commScalar.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(viscosityId));
    commScalar.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(densityId));
    
    // Create sweeps using GENERATED boundary conditions
    auto couetteFlowSweep = make_shared<CouetteFlowSweep>(
        densityId, forceId, pdfFieldId, temperatureId, velocityId, viscosityId);
        // densityId, forceId, pdfFieldId, pdfFieldTmpId, temperatureId, velocityId, viscosityId);
    
    // Create generated boundary sweeps
    auto topWallBC = make_shared<TopWallBC>(velocityId, wallVelocity);
    auto bottomWallBC = make_shared<BottomWallBC>(velocityId);

    // Apply pressure gradient force (once, during initialization)
    applyPressureGradientForce(blocks, forceId, pressureGradient);

    // TIME LOOP - USING GENERATED BOUNDARIES
    SweepTimeloop timeloop(blocks, timesteps);
    
    timeloop.add() 
        /*<< BeforeFunction([&]() {
        applyPressureGradientForce(blocks, forceId, drivingForce);
        }, "Driving Force")*/
        << BeforeFunction(commVector, "Vector Communication")
        << BeforeFunction(commScalar, "Scalar Communication")
        << Sweep([couetteFlowSweep](IBlock* block) { (*couetteFlowSweep)(block);}, "CouetteFlowSweep")
        << AfterFunction([&]() {
            // Apply generated boundary conditions
            applyBoundaryConditions(blocks, topWallBC, bottomWallBC);
        }, "Generated Boundary Conditions")
        /*<< AfterFunction([&]() {
            // Manual field swapping
            for (auto block = blocks->begin(); block != blocks->end(); ++block) {
                PdfField* pdf = block->getData<PdfField>(pdfFieldId);
                PdfField* pdf_tmp = block->getData<PdfField>(pdfFieldTmpId);
                pdf->swapDataPointers(pdf_tmp);
            }
        }, "Field Swapping")*/;
    
    // VTK output
    if (vtkWriteFrequency > 0) {
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "ThermalCouetteFlow", 
                                                      vtkWriteFrequency, 0, false, "vtk_out", 
                                                      "simulation_step", false, true, true, false, 0);
        
        auto velocityWriter = make_shared<field::VTKWriter<VectorField>>(velocityId, "velocity");
        auto temperatureWriter = make_shared<field::VTKWriter<ScalarField>>(temperatureId, "temperature");
        auto viscosityWriter = make_shared<field::VTKWriter<ScalarField>>(viscosityId, "viscosity");
        auto densityWriter = make_shared<field::VTKWriter<ScalarField>>(densityId, "density");
        auto forceWriter = make_shared<field::VTKWriter<VectorField>>(forceId, "force");
        
        vtkOutput->addCellDataWriter(velocityWriter);
        vtkOutput->addCellDataWriter(temperatureWriter);
        vtkOutput->addCellDataWriter(viscosityWriter);
        vtkOutput->addCellDataWriter(densityWriter);
        vtkOutput->addCellDataWriter(forceWriter);
        
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
    
    // Sample final results - MORE DETAILED
    WALBERLA_ROOT_SECTION() {
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            ScalarField* temperature = block->getData<ScalarField>(temperatureId);
            ScalarField* viscosity = block->getData<ScalarField>(viscosityId);
            VectorField* velocity = block->getData<VectorField>(velocityId);
            
            // Sample multiple points along y-direction
            Cell quarterCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells/4), static_cast<cell_idx_t>(zCells/2));
            Cell middleCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells/2), static_cast<cell_idx_t>(zCells/2));
            Cell threeQuarterCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(3*yCells/4), static_cast<cell_idx_t>(zCells/2));
            Cell topCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells-1), static_cast<cell_idx_t>(zCells/2));
            Cell bottomCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(0), static_cast<cell_idx_t>(zCells/2));
            
            WALBERLA_LOG_RESULT_ON_ROOT("=== VELOCITY PROFILE ===");
            WALBERLA_LOG_RESULT_ON_ROOT("Bottom (y=0): " << velocity->get(bottomCell, 0));
            WALBERLA_LOG_RESULT_ON_ROOT("Quarter (y=1/4): " << velocity->get(quarterCell, 0));
            WALBERLA_LOG_RESULT_ON_ROOT("Middle (y=1/2): " << velocity->get(middleCell, 0));
            WALBERLA_LOG_RESULT_ON_ROOT("3/4 (y=3/4): " << velocity->get(threeQuarterCell, 0));
            WALBERLA_LOG_RESULT_ON_ROOT("Top (y=1): " << velocity->get(topCell, 0));
            
            WALBERLA_LOG_RESULT_ON_ROOT("=== TEMPERATURE PROFILE ===");
            WALBERLA_LOG_RESULT_ON_ROOT("Bottom temp: " << temperature->get(bottomCell));
            WALBERLA_LOG_RESULT_ON_ROOT("Middle temp: " << temperature->get(middleCell));
            WALBERLA_LOG_RESULT_ON_ROOT("Top temp: " << temperature->get(topCell));
            
            WALBERLA_LOG_RESULT_ON_ROOT("=== VISCOSITY ===");
            WALBERLA_LOG_RESULT_ON_ROOT("Middle viscosity: " << viscosity->get(middleCell));
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
