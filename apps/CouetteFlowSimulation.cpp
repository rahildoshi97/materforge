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
    
    WALBERLA_LOG_INFO_ON_ROOT("=== FIELD INITIALIZATION DEBUG ===");
    WALBERLA_LOG_INFO_ON_ROOT("Parameters received:");
    WALBERLA_LOG_INFO_ON_ROOT("  - yCells: " << yCells);
    WALBERLA_LOG_INFO_ON_ROOT("  - Hot temperature: " << hotTemperature);
    WALBERLA_LOG_INFO_ON_ROOT("  - Cold temperature: " << coldTemperature);
    WALBERLA_LOG_INFO_ON_ROOT("  - Wall velocity: " << wallVelocity);
    WALBERLA_LOG_INFO_ON_ROOT("  - Reference viscosity: " << refViscosity);
    
    uint_t totalCellsInitialized = 0;
    uint_t blocksProcessed = 0;
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        blocksProcessed++;
        WALBERLA_LOG_INFO_ON_ROOT("Processing block " << blocksProcessed);
        
        // Get field pointers
        VectorField* velocity = block->getData<VectorField>(velocityId);
        ScalarField* temperature = block->getData<ScalarField>(temperatureId);
        ScalarField* density = block->getData<ScalarField>(densityId);
        VectorField* force = block->getData<VectorField>(forceId);
        ScalarField* viscosity = block->getData<ScalarField>(viscosityId);
        
        WALBERLA_LOG_INFO_ON_ROOT("Field pointers obtained successfully");
        WALBERLA_LOG_INFO_ON_ROOT("  - Velocity field: " << (velocity != nullptr ? "OK" : "NULL"));
        WALBERLA_LOG_INFO_ON_ROOT("  - Temperature field: " << (temperature != nullptr ? "OK" : "NULL"));
        WALBERLA_LOG_INFO_ON_ROOT("  - Density field: " << (density != nullptr ? "OK" : "NULL"));
        WALBERLA_LOG_INFO_ON_ROOT("  - Force field: " << (force != nullptr ? "OK" : "NULL"));
        WALBERLA_LOG_INFO_ON_ROOT("  - Viscosity field: " << (viscosity != nullptr ? "OK" : "NULL"));
        
        // Sample some cells for debugging
        uint_t cellsInBlock = 0;
        real_t minTemp = 1e10, maxTemp = -1e10;
        real_t minVel = 1e10, maxVel = -1e10;
        
        for (auto cell = velocity->beginXYZ(); cell != velocity->end(); ++cell) {
            Cell localCell = cell.cell();
            Cell globalCell = localCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, *block);
            
            real_t y_normalized = real_t(globalCell.y()) / real_t(yCells - 1);
            //y_normalized = std::max(0.0, std::min(1.0, y_normalized));
            
            // Linear temperature profile
            real_t temp = coldTemperature + (hotTemperature - coldTemperature) * y_normalized;
            temperature->get(localCell) = temp;
            
            // Initialize linear velocity profile
            real_t velX = wallVelocity * y_normalized;
            velocity->get(localCell, 0) = velX;
            velocity->get(localCell, 1) = 0.0;
            velocity->get(localCell, 2) = 0.0;

            density->get(localCell) = 1.0;
            
            // Initialize force
            force->get(localCell, 0) = 0.0;
            force->get(localCell, 1) = 0.0;
            force->get(localCell, 2) = 0.0;

            viscosity->get(localCell) = refViscosity;
            
            // Track statistics
            minTemp = std::min(minTemp, temp);
            maxTemp = std::max(maxTemp, temp);
            minVel = std::min(minVel, velX);
            maxVel = std::max(maxVel, velX);
            cellsInBlock++;
            totalCellsInitialized++;
            
            // Sample some cells for detailed logging
            if (globalCell.x() == 64 && globalCell.z() == 16 && globalCell.y() % 8 == 0) {
                WALBERLA_LOG_INFO_ON_ROOT("Sample cell (" << globalCell.x() << "," << globalCell.y() << "," << globalCell.z() << "):");
                WALBERLA_LOG_INFO_ON_ROOT("  y_norm=" << y_normalized << ", temp=" << temp << ", vel_x=" << velX);
            }
        }
        
        WALBERLA_LOG_INFO_ON_ROOT("Block " << blocksProcessed << " statistics:");
        WALBERLA_LOG_INFO_ON_ROOT("  - Cells initialized: " << cellsInBlock);
        WALBERLA_LOG_INFO_ON_ROOT("  - Temperature range: " << minTemp << " - " << maxTemp);
        WALBERLA_LOG_INFO_ON_ROOT("  - Velocity range: " << minVel << " - " << maxVel);

        // Initialize PDFs
        WALBERLA_LOG_INFO_ON_ROOT("Initializing PDFs for block " << blocksProcessed);
        try {
            CouetteFlowInit init(densityId, forceId, pdfFieldId, velocityId);
            init(&*block);
            WALBERLA_LOG_INFO_ON_ROOT("PDF initialization successful for block " << blocksProcessed);
        } catch (const std::exception& e) {
            WALBERLA_LOG_INFO_ON_ROOT("ERROR: PDF initialization failed: " << e.what());
            throw;
        }
    }
    
    WALBERLA_LOG_INFO_ON_ROOT("=== FIELD INITIALIZATION SUMMARY ===");
    WALBERLA_LOG_INFO_ON_ROOT("Total blocks processed: " << blocksProcessed);
    WALBERLA_LOG_INFO_ON_ROOT("Total cells initialized: " << totalCellsInitialized);
    WALBERLA_LOG_INFO_ON_ROOT("Expected cells per block: " << (128*32*32));
    WALBERLA_LOG_INFO_ON_ROOT("Field initialization complete");
}

// Debug velocity profile
void debugVelocityProfile(const shared_ptr<StructuredBlockForest>& blocks,
                         BlockDataID velocityId, uint_t timeStep, uint_t xCells, uint_t yCells, uint_t zCells) {
    
    if (timeStep % 500 != 0 && timeStep != 1) return;  // Log every 500 steps + step 1
    
    WALBERLA_LOG_INFO_ON_ROOT("=== VELOCITY DEBUG (Step " << timeStep << ") ===");
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        VectorField* velocity = block->getData<VectorField>(velocityId);
        
        // Sample centerline
        const cell_idx_t centerX = static_cast<cell_idx_t>(xCells/2);
        const cell_idx_t centerZ = static_cast<cell_idx_t>(zCells/2);
        
        WALBERLA_LOG_INFO_ON_ROOT("Centerline velocity profile (x=" << centerX << ", z=" << centerZ << "):");
        
        for (uint_t j = 0; j < yCells; j += yCells/8) {  // Sample 8 points
            Cell sampleCell(centerX, static_cast<cell_idx_t>(j), centerZ);
            real_t vel_x = velocity->get(sampleCell, 0);
            real_t vel_y = velocity->get(sampleCell, 1);
            real_t vel_z = velocity->get(sampleCell, 2);
            real_t y_norm = real_t(j) / real_t(yCells - 1);
            
            WALBERLA_LOG_INFO_ON_ROOT("  y=" << j << " (norm=" << y_norm << "): vel=(" 
                << vel_x << "," << vel_y << "," << vel_z << ")");
        }
        break;  // Only first block
    }
}

// Boundary condition application with debugging
void applySimpleCouetteBoundaries(const shared_ptr<StructuredBlockForest>& blocks,
                                  BlockDataID velocityId, real_t wallVelocity, uint_t timeStep) {
    
    bool shouldLog = (timeStep % 1000 == 0) || (timeStep <= 3);  // Log first 3 steps and every 1000th
    
    if (shouldLog) {
        WALBERLA_LOG_INFO_ON_ROOT("=== BOUNDARY CONDITIONS DEBUG (Step " << timeStep << ") ===");
    }
    
    uint_t totalTopCells = 0, totalBottomCells = 0;
    uint_t blocksWithTopBorder = 0, blocksWithBottomBorder = 0;
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        VectorField* velocity = block->getData<VectorField>(velocityId);
        CellInterval domain = velocity->xyzSize();
        
        if (shouldLog) {
            WALBERLA_LOG_INFO_ON_ROOT("Block domain: " << domain);
        }
        
        // TOP WALL (Moving wall)
        if (blocks->atDomainYMaxBorder(*block)) {
            blocksWithTopBorder++;
            int topY = domain.yMax();
            uint_t topCellsSet = 0;
            
            for (int x = domain.xMin(); x <= domain.xMax(); ++x) {
                for (int z = domain.zMin(); z <= domain.zMax(); ++z) {
                    Cell topCell(x, topY, z);
                    velocity->get(topCell, 0) = wallVelocity;
                    velocity->get(topCell, 1) = 0.0;
                    velocity->get(topCell, 2) = 0.0;
                    topCellsSet++;
                }
            }
            totalTopCells += topCellsSet;
            
            if (shouldLog) {
                WALBERLA_LOG_INFO_ON_ROOT("  Top wall: y=" << topY << ", cells set=" << topCellsSet << ", velocity=" << wallVelocity);
            }
        }
        
        // BOTTOM WALL (Stationary wall)
        if (blocks->atDomainYMinBorder(*block)) {
            blocksWithBottomBorder++;
            int bottomY = domain.yMin();
            uint_t bottomCellsSet = 0;
            
            for (int x = domain.xMin(); x <= domain.xMax(); ++x) {
                for (int z = domain.zMin(); z <= domain.zMax(); ++z) {
                    Cell bottomCell(x, bottomY, z);
                    velocity->get(bottomCell, 0) = 0.0;
                    velocity->get(bottomCell, 1) = 0.0;
                    velocity->get(bottomCell, 2) = 0.0;
                    bottomCellsSet++;
                }
            }
            totalBottomCells += bottomCellsSet;
            
            if (shouldLog) {
                WALBERLA_LOG_INFO_ON_ROOT("  Bottom wall: y=" << bottomY << ", cells set=" << bottomCellsSet << ", velocity=0.0");
            }
        }
    }
    
    if (shouldLog) {
        WALBERLA_LOG_INFO_ON_ROOT("Boundary conditions summary:");
        WALBERLA_LOG_INFO_ON_ROOT("  - Blocks with top border: " << blocksWithTopBorder);
        WALBERLA_LOG_INFO_ON_ROOT("  - Blocks with bottom border: " << blocksWithBottomBorder);
        WALBERLA_LOG_INFO_ON_ROOT("  - Total top cells: " << totalTopCells);
        WALBERLA_LOG_INFO_ON_ROOT("  - Total bottom cells: " << totalBottomCells);
    }
}

// Grid-independent pressure gradient force with debugging
void applyPressureGradientForce(const shared_ptr<StructuredBlockForest>& blocks,
                               BlockDataID forceId, real_t pressureGradient) {
    
    WALBERLA_LOG_INFO_ON_ROOT("=== PRESSURE GRADIENT FORCE DEBUG ===");
    WALBERLA_LOG_INFO_ON_ROOT("Applying pressure gradient: " << pressureGradient);
    
    uint_t totalCells = 0;
    
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
        VectorField* force = block->getData<VectorField>(forceId);
        
        for (auto cell = force->beginXYZ(); cell != force->end(); ++cell) {
            Cell localCell = cell.cell();
            
            // Constant pressure gradient force in flow direction (x)
            force->get(localCell, 0) = pressureGradient;
            force->get(localCell, 1) = 0.0;
            force->get(localCell, 2) = 0.0;
            totalCells++;
        }
    }
    
    WALBERLA_LOG_INFO_ON_ROOT("Pressure gradient applied to " << totalCells << " cells");
    WALBERLA_LOG_INFO_ON_ROOT("Force vector: (" << pressureGradient << ", 0.0, 0.0)");
}

int main(int argc, char** argv) {
    mpi::Environment env(argc, argv);
    
    // Parameters
    const uint_t timesteps = 5000;
    const uint_t vtkWriteFrequency = -50;
    const real_t wallVelocity = 0.1;
    const real_t hotTemperature = 600.0;
    const real_t coldTemperature = 300.0;
    const uint_t xBlocks = 1, yBlocks = 1, zBlocks = 1;
    const uint_t xCells = 128, yCells = 64, zCells = 32;
    const real_t refViscosity = 0.1667;
    // const real_t pressureGradient = 0.0005;

    WALBERLA_LOG_INFO_ON_ROOT("=== THERMAL COUETTE FLOW SIMULATION DEBUG ===");
    WALBERLA_LOG_INFO_ON_ROOT("Simulation parameters:");
    WALBERLA_LOG_INFO_ON_ROOT("  - Grid: " << xCells << " x " << yCells << " x " << zCells);
    WALBERLA_LOG_INFO_ON_ROOT("  - Blocks: " << xBlocks << " x " << yBlocks << " x " << zBlocks);
    WALBERLA_LOG_INFO_ON_ROOT("  - Wall velocity: " << wallVelocity);
    WALBERLA_LOG_INFO_ON_ROOT("  - Temperature range: " << coldTemperature << " - " << hotTemperature);
    WALBERLA_LOG_INFO_ON_ROOT("  - Reference viscosity: " << refViscosity);
    // WALBERLA_LOG_INFO_ON_ROOT("  - Pressure gradient: " << pressureGradient);
    WALBERLA_LOG_INFO_ON_ROOT("  - Timesteps: " << timesteps);
    WALBERLA_LOG_INFO_ON_ROOT("  - VTK frequency: " << vtkWriteFrequency);

    // Create block forest
    WALBERLA_LOG_INFO_ON_ROOT("Creating block forest...");
    auto aabb = math::AABB(0, 0, 0, xBlocks*xCells, yBlocks*yCells, zBlocks*zCells);
    WALBERLA_LOG_INFO_ON_ROOT("Domain AABB: " << aabb);
    
    shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGrid(
        aabb, xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, 
        true,   // oneBlockPerProcess  
        true,   // xPeriodic = true (periodic in x)
        false,  // yPeriodic = false (walls in y)
        true,   // zPeriodic = true (periodic in z)
        false   // keepGlobalBlockInformation = false
    );
    
    WALBERLA_LOG_INFO_ON_ROOT("Block forest created successfully");
    WALBERLA_LOG_INFO_ON_ROOT("  - Total blocks: " << blocks->getNumberOfBlocks());
    WALBERLA_LOG_INFO_ON_ROOT("  - Periodicity: x=true, y=false, z=true");

    // Add fields
    WALBERLA_LOG_INFO_ON_ROOT("Adding fields to storage...");
    BlockDataID densityId = field::addToStorage<ScalarField>(blocks, "density", 1.0, field::fzyx, 1);
    BlockDataID forceId = field::addToStorage<VectorField>(blocks, "force", real_t(0.0), field::fzyx, 1);
    BlockDataID pdfFieldId = field::addToStorage<PdfField>(blocks, "pdfs", real_t(0.0), field::fzyx, 1);
    BlockDataID temperatureId = field::addToStorage<ScalarField>(blocks, "temperature", 300.0, field::fzyx, 1);
    BlockDataID velocityId = field::addToStorage<VectorField>(blocks, "velocity", real_t(0.0), field::fzyx, 1);
    BlockDataID viscosityId = field::addToStorage<ScalarField>(blocks, "viscosity", refViscosity, field::fzyx, 1);
    
    WALBERLA_LOG_INFO_ON_ROOT("All fields added successfully:");
    WALBERLA_LOG_INFO_ON_ROOT("  - Density (scalar): " << densityId);
    WALBERLA_LOG_INFO_ON_ROOT("  - Force (vector): " << forceId);  
    WALBERLA_LOG_INFO_ON_ROOT("  - PDFs (19-component): " << pdfFieldId);
    WALBERLA_LOG_INFO_ON_ROOT("  - Temperature (scalar): " << temperatureId);
    WALBERLA_LOG_INFO_ON_ROOT("  - Velocity (vector): " << velocityId);
    WALBERLA_LOG_INFO_ON_ROOT("  - Viscosity (scalar): " << viscosityId);
    
    // Initialize fields
    WALBERLA_LOG_INFO_ON_ROOT("Initializing fields...");
    initializeFields(blocks, densityId, forceId, pdfFieldId, temperatureId, 
                    velocityId, viscosityId, hotTemperature, coldTemperature, wallVelocity, yCells, refViscosity);
    
    // Apply driving force
    // WALBERLA_LOG_INFO_ON_ROOT("Applying pressure gradient force...");
    // applyPressureGradientForce(blocks, forceId, pressureGradient);

    // Communication scheme
    WALBERLA_LOG_INFO_ON_ROOT("Setting up communication scheme...");
    blockforest::communication::UniformBufferedScheme<stencil::D3Q19> commScheme(blocks);
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(velocityId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(temperatureId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(viscosityId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<ScalarField>>(densityId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<PdfField>>(pdfFieldId));
    commScheme.addPackInfo(make_shared<field::communication::PackInfo<VectorField>>(forceId));
    WALBERLA_LOG_INFO_ON_ROOT("Communication scheme configured");

    // Create sweeps
    WALBERLA_LOG_INFO_ON_ROOT("Creating LBM sweeps...");
    try {
        auto couetteFlowSweep = make_shared<CouetteFlowSweep>(
            // densityId, forceId, pdfFieldId, temperatureId, velocityId, viscosityId);
            densityId, forceId, pdfFieldId, velocityId, viscosityId);
        WALBERLA_LOG_INFO_ON_ROOT("CouetteFlowSweep created successfully");
        WALBERLA_LOG_INFO_ON_ROOT("Constructor parameters order: density, force, pdfs, temperature, velocity, viscosity");
    } catch (const std::exception& e) {
        WALBERLA_LOG_INFO_ON_ROOT("ERROR: Failed to create CouetteFlowSweep: " << e.what());
        return EXIT_FAILURE;
    }

    auto couetteFlowSweep = make_shared<CouetteFlowSweep>(
        // densityId, forceId, pdfFieldId, temperatureId, velocityId, viscosityId);
        densityId, forceId, pdfFieldId, velocityId, viscosityId);

    // Initial velocity debug
    debugVelocityProfile(blocks, velocityId, 0, xCells, yCells, zCells);

    // TIME LOOP with debugging
    WALBERLA_LOG_INFO_ON_ROOT("Setting up time loop...");
    SweepTimeloop timeloop(blocks, timesteps);
    
    uint_t stepCounter = 0;
    
    timeloop.add() 
        << BeforeFunction(commScheme, "Communication")
        << Sweep([couetteFlowSweep](IBlock* block) { 
            (*couetteFlowSweep)(block);
        }, "CouetteFlowSweep")
        << AfterFunction([&]() {
            stepCounter++;
            applySimpleCouetteBoundaries(blocks, velocityId, wallVelocity, stepCounter);
            debugVelocityProfile(blocks, velocityId, stepCounter, xCells, yCells, zCells);
        }, "Boundary Conditions + Debug");
    
    WALBERLA_LOG_INFO_ON_ROOT("Time loop configured");
    
    // VTK output setup
    if (vtkWriteFrequency > 0) {
        WALBERLA_LOG_INFO_ON_ROOT("Setting up VTK output...");
        auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "ThermalCouetteFlow_DEBUG", 
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
        
        WALBERLA_LOG_INFO_ON_ROOT("VTK output configured - files will be written every " << vtkWriteFrequency << " steps");
    }
    
    // Run simulation
    WcTimingPool timeloopTiming;
    WcTimer simTimer;

    WALBERLA_LOG_INFO_ON_ROOT("=== STARTING SIMULATION ===");
    WALBERLA_LOG_INFO_ON_ROOT("Running " << timesteps << " timesteps...");
    simTimer.start();
    timeloop.run(timeloopTiming);
    simTimer.end();
    
    // Results
    real_t simTime = real_c(simTimer.max());
    WALBERLA_MPI_SECTION() { 
        walberla::mpi::reduceInplace(simTime, walberla::mpi::MAX); 
    }
    
    const auto reducedTimeloopTiming = timeloopTiming.getReduced();
    
    WALBERLA_LOG_RESULT_ON_ROOT("=== THERMAL COUETTE RESULTS (DEBUGGING) ===");
    WALBERLA_LOG_RESULT_ON_ROOT("Total time: " << simTime << " seconds");
    WALBERLA_LOG_RESULT_ON_ROOT("Time per timestep: " << (simTime / timesteps * 1000) << " ms");
    WALBERLA_LOG_RESULT_ON_ROOT("Performance timing:\n" << *reducedTimeloopTiming);
    
    // Final comprehensive analysis with debugging
    WALBERLA_ROOT_SECTION() {
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
            ScalarField* temperature = block->getData<ScalarField>(temperatureId);
            ScalarField* viscosity = block->getData<ScalarField>(viscosityId);
            VectorField* velocity = block->getData<VectorField>(velocityId);
            VectorField* force = block->getData<VectorField>(forceId);
            
            Cell quarterCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells/4), static_cast<cell_idx_t>(zCells/2));
            Cell middleCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells/2), static_cast<cell_idx_t>(zCells/2));
            Cell threeQuarterCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(3*yCells/4), static_cast<cell_idx_t>(zCells/2));
            Cell topCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(yCells-1), static_cast<cell_idx_t>(zCells/2));
            Cell bottomCell(static_cast<cell_idx_t>(xCells/2), static_cast<cell_idx_t>(0), static_cast<cell_idx_t>(zCells/2));
            
            WALBERLA_LOG_RESULT_ON_ROOT("=== FINAL VELOCITY PROFILE DEBUG ===");
            WALBERLA_LOG_RESULT_ON_ROOT("Cell locations (center line x=" << xCells/2 << ", z=" << zCells/2 << "):");
            WALBERLA_LOG_RESULT_ON_ROOT("Bottom (y=0): " << velocity->get(bottomCell, 0) << " (expected: 0.000)");
            WALBERLA_LOG_RESULT_ON_ROOT("Quarter (y=" << yCells/4 << "): " << velocity->get(quarterCell, 0) << " (expected: " << wallVelocity*0.25 << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("Middle (y=" << yCells/2 << "): " << velocity->get(middleCell, 0) << " (expected: " << wallVelocity*0.5 << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("3/4 (y=" << 3*yCells/4 << "): " << velocity->get(threeQuarterCell, 0) << " (expected: " << wallVelocity*0.75 << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("Top (y=" << yCells-1 << "): " << velocity->get(topCell, 0) << " (expected: " << wallVelocity << ")");
            
            // Force field check
            WALBERLA_LOG_RESULT_ON_ROOT("=== FORCE FIELD DEBUG ===");
            WALBERLA_LOG_RESULT_ON_ROOT("Force at middle cell: (" 
                << force->get(middleCell, 0) << ", " 
                << force->get(middleCell, 1) << ", " 
                << force->get(middleCell, 2) << ")");
            // WALBERLA_LOG_RESULT_ON_ROOT("Expected force: (" << pressureGradient << ", 0.0, 0.0)");
            
            real_t bottomTemp = temperature->get(bottomCell);
            real_t middleTemp = temperature->get(middleCell);
            real_t topTemp = temperature->get(topCell);
            
            WALBERLA_LOG_RESULT_ON_ROOT("Temperature profile:");
            WALBERLA_LOG_RESULT_ON_ROOT("  Bottom: " << bottomTemp << " K (expected: " << coldTemperature << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("  Middle: " << middleTemp << " K (expected: " << (coldTemperature+hotTemperature)/2.0 << ")");
            WALBERLA_LOG_RESULT_ON_ROOT("  Top: " << topTemp << " K (expected: " << hotTemperature << ")");
            
            // Overall assessment
            real_t maxVelError = std::max({
                std::abs(velocity->get(quarterCell, 0) - wallVelocity*0.25),
                std::abs(velocity->get(middleCell, 0) - wallVelocity*0.5),
                std::abs(velocity->get(threeQuarterCell, 0) - wallVelocity*0.75)
            });
            
            WALBERLA_LOG_RESULT_ON_ROOT("=== SIMULATION ASSESSMENT ===");
            WALBERLA_LOG_RESULT_ON_ROOT("Maximum velocity error: " << maxVelError);
            if (maxVelError < 0.005) {
                WALBERLA_LOG_RESULT_ON_ROOT("✅ VELOCITY PROFILE: EXCELLENT");
            } else if (maxVelError < 0.025) {
                WALBERLA_LOG_RESULT_ON_ROOT("✅ VELOCITY PROFILE: GOOD");
            } else {
                WALBERLA_LOG_RESULT_ON_ROOT("❌ VELOCITY PROFILE: NEEDS IMPROVEMENT");
                WALBERLA_LOG_RESULT_ON_ROOT("   Possible causes: insufficient timesteps, force too weak, or instability");
            }

            break;
        }
    }
    
    WALBERLA_LOG_INFO_ON_ROOT("=== SIMULATION COMPLETED ===");
    return EXIT_SUCCESS;
}

} // namespace couette_flow_thermal
} // namespace walberla

int main(int argc, char** argv) { 
    return walberla::couette_flow_thermal::main(argc, argv); 
}
