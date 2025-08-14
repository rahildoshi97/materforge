cmake_minimum_required(VERSION 3.24)

project(test_materforge LANGUAGES CXX)

SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# waLBerla Options
SET(WALBERLA_BUILD_TESTS OFF CACHE BOOL "")
SET(WALBERLA_BUILD_BENCHMARKS OFF CACHE BOOL "")
SET(WALBERLA_BUILD_TOOLS OFF CACHE BOOL "")
SET(WALBERLA_BUILD_TUTORIALS OFF CACHE BOOL "")
SET(WALBERLA_BUILD_SHOWCASES OFF CACHE BOOL "")
SET(WALBERLA_BUILD_DOC OFF CACHE BOOL "")
SET(WALBERLA_BUILD_WITH_PYTHON OFF CACHE BOOL "")
SET(WALBERLA_BUILD_WITH_CODEGEN OFF CACHE BOOL "")
SET(WALBERLA_BUILD_WITH_OPENMP ON CACHE BOOL "")
SET(WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "")

# Conditional GPU support
option(ENABLE_GPU_SUPPORT "Enable GPU support" ON)

# Check if CUDA was explicitly requested via command line
if(DEFINED WALBERLA_BUILD_WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
    message(STATUS "CUDA explicitly requested via WALBERLA_BUILD_WITH_CUDA")
    set(CUDA_REQUESTED TRUE)
else()
    set(CUDA_REQUESTED FALSE)
endif()

if(ENABLE_GPU_SUPPORT OR CUDA_REQUESTED)
    # Check if CUDA is available
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit QUIET)
        if (CUDAToolkit_FOUND)
            SET(WALBERLA_BUILD_WITH_GPU_SUPPORT ON CACHE BOOL "" FORCE)
            SET(WALBERLA_BUILD_WITH_CUDA ON CACHE BOOL "" FORCE)
            message(STATUS "CUDA found - GPU support enabled")
            set(GPU_AVAILABLE TRUE)
        else()
            message(WARNING "CUDA compiler found but toolkit missing - GPU support disabled")
            SET(WALBERLA_BUILD_WITH_GPU_SUPPORT OFF CACHE BOOL "" FORCE)
            SET(WALBERLA_BUILD_WITH_CUDA OFF CACHE BOOL "" FORCE)
            set(GPU_AVAILABLE FALSE)
        endif()
    else()
        if(CUDA_REQUESTED)
            message(FATAL_ERROR "CUDA was explicitly requested but compiler not found - aborting configuration")
        else()
            message(STATUS "CUDA compiler not found - GPU support disabled")
            SET(WALBERLA_BUILD_WITH_GPU_SUPPORT OFF CACHE BOOL "" FORCE)
            SET(WALBERLA_BUILD_WITH_CUDA OFF CACHE BOOL "" FORCE)
            set(GPU_AVAILABLE FALSE)
        endif()
    endif()
else()
    SET(WALBERLA_BUILD_WITH_GPU_SUPPORT OFF CACHE BOOL "" FORCE)
    SET(WALBERLA_BUILD_WITH_CUDA OFF CACHE BOOL "" FORCE)
    set(GPU_AVAILABLE FALSE)
    message(STATUS "GPU support manually disabled")
endif()

add_subdirectory(walberla)
add_subdirectory(sfg-walberla)

walberla_codegen_venv_require(-e ${CMAKE_SOURCE_DIR}/..)
walberla_codegen_venv_populate()

# CPU-only executables (always built)
add_executable(HeatEquationWithMaterial)
target_sources(HeatEquationWithMaterial PRIVATE CodegenHeatEquationWithMaterial.cpp)
target_link_libraries(HeatEquationWithMaterial PRIVATE
        walberla::blockforest walberla::core walberla::field
        walberla::stencil walberla::timeloop walberla::vtk walberla::pde)

walberla_generate_sources(HeatEquationWithMaterial
        SCRIPTS HeatEquationKernelWithMaterial.py
        FILE_EXTENSIONS .cpp .hpp)

add_executable(HeatEquationWithMaterialCPU)
target_sources(HeatEquationWithMaterialCPU PRIVATE CodegenHeatEquationWithMaterialCPU.cpp)
target_link_libraries(HeatEquationWithMaterialCPU PRIVATE
        walberla::blockforest walberla::core walberla::field
        walberla::stencil walberla::timeloop walberla::vtk walberla::pde)

walberla_generate_sources(HeatEquationWithMaterialCPU
        SCRIPTS HeatEquationKernelWithMaterialCPU.py
        FILE_EXTENSIONS .cpp .hpp)

add_executable(CodegenHeatEquationCPUScaling)
target_sources(CodegenHeatEquationCPUScaling PRIVATE CodegenHeatEquationCPUScaling.cpp)
target_link_libraries(CodegenHeatEquationCPUScaling PRIVATE
        walberla::blockforest walberla::core walberla::field
        walberla::stencil walberla::timeloop walberla::vtk walberla::pde)

walberla_generate_sources(CodegenHeatEquationCPUScaling
        SCRIPTS HeatEquationKernelWithMaterialCPU.py
        FILE_EXTENSIONS .cpp .hpp)

# GPU executables (conditionally built)
if(GPU_AVAILABLE)
    add_executable(HeatEquationWithMaterialGPU)
    target_sources(HeatEquationWithMaterialGPU PRIVATE CodegenHeatEquationWithMaterialGPU.cpp)
    target_link_libraries(HeatEquationWithMaterialGPU PRIVATE
            walberla::blockforest walberla::core walberla::gpu walberla::field
            walberla::stencil walberla::timeloop walberla::vtk walberla::pde)
    # Set CUDA properties
    set_property(TARGET HeatEquationWithMaterialGPU PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    walberla_generate_sources(HeatEquationWithMaterialGPU
            SCRIPTS HeatEquationKernelWithMaterialGPU.py
            FILE_EXTENSIONS .cu .hpp)
    set_source_files_properties(
            ${CMAKE_CURRENT_BINARY_DIR}/_gen/HeatEquationWithMaterialGPU/gen/HeatEquationKernelWithMaterialGPU.cu
            PROPERTIES LANGUAGE CUDA)

    # GPU Scaling executable
    add_executable(CodegenHeatEquationGPUScaling)
    target_sources(CodegenHeatEquationGPUScaling PRIVATE CodegenHeatEquationGPUScaling.cpp)
    target_link_libraries(CodegenHeatEquationGPUScaling PRIVATE
            walberla::blockforest walberla::core walberla::gpu walberla::field
            walberla::stencil walberla::timeloop walberla::vtk walberla::pde)
    # Set CUDA properties for scaling executable
    set_property(TARGET CodegenHeatEquationGPUScaling PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    walberla_generate_sources(CodegenHeatEquationGPUScaling
            SCRIPTS HeatEquationKernelWithMaterialGPU.py
            FILE_EXTENSIONS .cu .hpp)
    set_source_files_properties(
            ${CMAKE_CURRENT_BINARY_DIR}/_gen/CodegenHeatEquationGPUScaling/gen/HeatEquationKernelWithMaterialGPU.cu
            PROPERTIES LANGUAGE CUDA)

    message(STATUS "GPU executables will be built: HeatEquationWithMaterialGPU, CodegenHeatEquationGPUScaling")
else()
    message(STATUS "GPU executables will NOT be built (CUDA not available)")
endif()
