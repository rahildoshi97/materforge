# CMakeLists.txt - Unified CPU/GPU Build (Single Executable)

cmake_minimum_required(VERSION 3.24)
project(materforge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# BUILD CONFIGURATION
# ============================================================================

# Single source of truth: TARGET_PLATFORM
# Options: "CPU", "GPU-HIP", "GPU-CUDA"
set(TARGET_PLATFORM "GPU-HIP" CACHE STRING "Target platform: CPU, GPU-HIP, or GPU-CUDA")
set_property(CACHE TARGET_PLATFORM PROPERTY STRINGS "CPU" "GPU-HIP" "GPU-CUDA")

message(STATUS "═══════════════════════════════════════════════════════════")
message(STATUS "TARGET_PLATFORM: ${TARGET_PLATFORM}")
message(STATUS "═══════════════════════════════════════════════════════════")

# ============================================================================
# LANGUAGE SETUP (BEFORE any find_package or add_subdirectory)
# ============================================================================

if(TARGET_PLATFORM STREQUAL "GPU-HIP")
    enable_language(HIP)
    set(HIP_ARCHITECTURES "gfx90a")
    set(CMAKE_HIP_ARCHITECTURES "gfx90a")
    set(AMDGPU_TARGETS "gfx90a")
    message(STATUS "✅ HIP language enabled for gfx90a")
elseif(TARGET_PLATFORM STREQUAL "GPU-CUDA")
    enable_language(CUDA)
    message(STATUS "✅ CUDA language enabled")
endif()

# ============================================================================
# WALBERLA CONFIGURATION (before add_subdirectory)
# ============================================================================

# Unified settings for all builds
set(WALBERLA_BUILD_TESTS OFF CACHE BOOL "")
set(WALBERLA_BUILD_BENCHMARKS OFF CACHE BOOL "")
set(WALBERLA_BUILD_TOOLS OFF CACHE BOOL "")
set(WALBERLA_BUILD_TUTORIALS OFF CACHE BOOL "")
set(WALBERLA_BUILD_SHOWCASES OFF CACHE BOOL "")
set(WALBERLA_BUILD_DOC OFF CACHE BOOL "")
set(WALBERLA_BUILD_EXAMPLES OFF CACHE BOOL "")
set(WALBERLA_BUILD_WITH_PYTHON ON CACHE BOOL "")
set(WALBERLA_ENABLE_SWEEPGEN ON CACHE BOOL "")
set(WALBERLA_SWEEPGEN_MANAGED_VENV ON CACHE BOOL "")

# Configure GPU/CPU based on TARGET_PLATFORM
if(TARGET_PLATFORM STREQUAL "CPU")
    set(WALBERLA_BUILD_WITH_GPU_SUPPORT OFF CACHE BOOL "" FORCE)
    set(WALBERLA_BUILD_WITH_HIP OFF CACHE BOOL "" FORCE)
    set(WALBERLA_BUILD_WITH_CUDA OFF CACHE BOOL "" FORCE)
    message(STATUS "Configured for CPU-only")
    
elseif(TARGET_PLATFORM STREQUAL "GPU-HIP")
    set(WALBERLA_BUILD_WITH_GPU_SUPPORT ON CACHE BOOL "" FORCE)
    set(WALBERLA_BUILD_WITH_HIP ON CACHE BOOL "" FORCE)
    set(WALBERLA_BUILD_WITH_CUDA OFF CACHE BOOL "" FORCE)
    message(STATUS "Configured for HIP GPU")
    
elseif(TARGET_PLATFORM STREQUAL "GPU-CUDA")
    set(WALBERLA_BUILD_WITH_GPU_SUPPORT ON CACHE BOOL "" FORCE)
    set(WALBERLA_BUILD_WITH_HIP OFF CACHE BOOL "" FORCE)
    set(WALBERLA_BUILD_WITH_CUDA ON CACHE BOOL "" FORCE)
    message(STATUS "Configured for CUDA GPU")
endif()

# ============================================================================
# DEPENDENCIES
# ============================================================================

find_package(MPI REQUIRED)
message(STATUS "✅ MPI found")

# Find HIP if needed
if(TARGET_PLATFORM STREQUAL "GPU-HIP")
    find_package(hip REQUIRED)
    message(STATUS "✅ HIP found")
endif()

# Find CUDA if needed
if(TARGET_PLATFORM STREQUAL "GPU-CUDA")
    find_package(CUDA REQUIRED)
    message(STATUS "✅ CUDA found")
endif()

# ============================================================================
# WALBERLA SUBPROJECT
# ============================================================================

add_subdirectory(walberla)
walberla_link_files_to_builddir("*.prm")
sweepgen_venv_require(-e ${CMAKE_SOURCE_DIR}/..)
sweepgen_venv_populate()

# ============================================================================
# UNIFIED EXECUTABLE FUNCTION
# ============================================================================

# Single executable that works on both CPU and GPU
# The sweep script generates code for the target platform automatically
function(add_unified_executable TARGET_NAME SOURCE_FILE SWEEP_SCRIPT)
    add_executable(${TARGET_NAME})
    target_sources(${TARGET_NAME} PRIVATE ${SOURCE_FILE})
    
    # Common libraries for all targets
    target_link_libraries(${TARGET_NAME} PRIVATE
        walberla::blockforest walberla::boundary walberla::core 
        walberla::experimental walberla::gpu walberla::field 
        walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
        walberla::stencil walberla::timeloop walberla::vtk
        walberla::pde
        MPI::MPI_CXX)
    
    # Determine target platform for code generation
    if(TARGET_PLATFORM STREQUAL "CPU")
        set(CODEGEN_TARGET "cpu")
    else()
        set(CODEGEN_TARGET "gpu")
    endif()
    
    # Generate sources with target specification
    if(TARGET_PLATFORM STREQUAL "GPU-HIP")
        target_link_libraries(${TARGET_NAME} PRIVATE hip::device)
        
        walberla_generate_sources(${TARGET_NAME}
            SCRIPTS ${SWEEP_SCRIPT}
            SCRIPT_ARGS --target=${CODEGEN_TARGET}
            FILE_EXTENSIONS .cpp .hpp)
        
        # Set HIP language for generated files
        set_source_files_properties(
            ${CMAKE_CURRENT_BINARY_DIR}/_gen/${TARGET_NAME}/gen/*.cpp
            PROPERTIES LANGUAGE HIP)
        
        target_compile_options(${TARGET_NAME} PRIVATE 
            $<$<COMPILE_LANGUAGE:HIP>:--offload-arch=gfx90a>)
        set_property(TARGET ${TARGET_NAME} PROPERTY HIP_ARCHITECTURES "gfx90a")
        
        # Link GTL library for LUMI
        if(DEFINED ENV{PE_MPICH_GTL_DIR_amd_gfx90a})
            string(REPLACE "-L" "" GTL_DIR "$ENV{PE_MPICH_GTL_DIR_amd_gfx90a}")
            string(REPLACE "-l" "" GTL_LIB "$ENV{PE_MPICH_GTL_LIBS_amd_gfx90a}")
            target_link_directories(${TARGET_NAME} PRIVATE ${GTL_DIR})
            target_link_libraries(${TARGET_NAME} PRIVATE ${GTL_LIB})
            message(STATUS "Linked GTL library for ${TARGET_NAME}")
        endif()
        
    elseif(TARGET_PLATFORM STREQUAL "GPU-CUDA")
        target_link_libraries(${TARGET_NAME} PRIVATE walberla::gpu)
        
        walberla_generate_sources(${TARGET_NAME}
            SCRIPTS ${SWEEP_SCRIPT}
            SCRIPT_ARGS --target=${CODEGEN_TARGET}
            FILE_EXTENSIONS .cu .hpp)
        
        set_source_files_properties(
            ${CMAKE_CURRENT_BINARY_DIR}/_gen/${TARGET_NAME}/gen/*.cu
            PROPERTIES LANGUAGE CUDA)
        
        set_property(TARGET ${TARGET_NAME} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
        
    else()  # CPU
        walberla_generate_sources(${TARGET_NAME}
            SCRIPTS ${SWEEP_SCRIPT}
            SCRIPT_ARGS --target=${CODEGEN_TARGET}
            FILE_EXTENSIONS .cpp .hpp)
    endif()
endfunction()

# ============================================================================
# BUILD UNIFIED COUETTE FLOW SCALING EXECUTABLE
# ============================================================================

# Single executable for CPU/GPU
# Use the unified CouetteFlowScaling.cpp and CouetteFlowSweeps.py
add_unified_executable(CouetteFlowScaling
    CouetteFlowScaling.cpp
    CouetteFlowSweeps.py)

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "═══════════════════════════════════════════════════════════")
if(TARGET_PLATFORM STREQUAL "CPU")
    message(STATUS "✅ CPU executable configured:")
    message(STATUS "  - CouetteFlowScaling (CPU)")
elseif(TARGET_PLATFORM STREQUAL "GPU-HIP")
    message(STATUS "✅ GPU (HIP) executable configured:")
    message(STATUS "  - CouetteFlowScaling (GPU-HIP, unified)")
elseif(TARGET_PLATFORM STREQUAL "GPU-CUDA")
    message(STATUS "✅ GPU (CUDA) executable configured:")
    message(STATUS "  - CouetteFlowScaling (GPU-CUDA, unified)")
endif()
message(STATUS "═══════════════════════════════════════════════════════════")
