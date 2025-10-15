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
SET(WALBERLA_BUILD_WITH_PYTHON ON CACHE BOOL "")
# SET(WALBERLA_BUILD_WITH_CODEGEN OFF CACHE BOOL "")
# SET(WALBERLA_BUILD_WITH_OPENMP ON CACHE BOOL "")
# SET(WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "")
SET(WALBERLA_ENABLE_SWEEPGEN ON CACHE BOOL "")
SET(WALBERLA_SWEEPGEN_MANAGED_VENV ON CACHE BOOL "")

# Build mode detection
option(BUILD_CPU_ONLY "Build only CPU targets (ignore GPU)" ON) # Default to CPU
option(ENABLE_GPU_SUPPORT "Enable GPU support" OFF) # Default to OFF

# Configure waLBerla based on build mode
if(BUILD_CPU_ONLY)
    message(STATUS "CPU-only build - configuring waLBerla without GPU support")
    SET(WALBERLA_BUILD_WITH_GPU_SUPPORT OFF CACHE BOOL "" FORCE)
    SET(WALBERLA_BUILD_WITH_CUDA OFF CACHE BOOL "" FORCE)
    SET(WALBERLA_BUILD_WITH_HIP OFF CACHE BOOL "" FORCE)
    set(GPU_AVAILABLE FALSE)
else()
    # Clear problematic HIP environment variables BEFORE detection
    unset(ENV{HIP_FLAGS})
    unset(ENV{HIPCXX_FLAGS})
    unset(ENV{CMAKE_HIP_FLAGS})
    
    # Clear CMake HIP flags early
    set(CMAKE_HIP_FLAGS "" CACHE STRING "HIP compiler flags" FORCE)
    set(CMAKE_HIP_FLAGS_DEBUG "" CACHE STRING "HIP debug flags" FORCE) 
    set(CMAKE_HIP_FLAGS_RELEASE "" CACHE STRING "HIP release flags" FORCE)
    
    # GPU build - detect HIP or CUDA and enable language EARLY
    find_package(hip QUIET)
    if(hip_FOUND)
        message(STATUS "HIP found - enabling HIP language and configuring waLBerla for HIP support")
        
        # Enable HIP language BEFORE configuring waLBerla
        enable_language(HIP)
        
        # Set HIP flags after language is enabled
        set(CMAKE_HIP_FLAGS "--offload-arch=gfx90a")
        set(CMAKE_HIP_ARCHITECTURES "gfx90a")
                
        SET(WALBERLA_BUILD_WITH_GPU_SUPPORT ON CACHE BOOL "" FORCE)
        SET(WALBERLA_BUILD_WITH_HIP ON CACHE BOOL "" FORCE)
        SET(WALBERLA_BUILD_WITH_CUDA OFF CACHE BOOL "" FORCE)
        set(GPU_AVAILABLE TRUE)
        set(GPU_TYPE "HIP")
                
    else()
        # Try CUDA
        include(CheckLanguage)
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
            enable_language(CUDA)
            find_package(CUDAToolkit QUIET)
            if(CUDAToolkit_FOUND)
                message(STATUS "CUDA found - configuring waLBerla for CUDA support")
                SET(WALBERLA_BUILD_WITH_GPU_SUPPORT ON CACHE BOOL "" FORCE)
                SET(WALBERLA_BUILD_WITH_CUDA ON CACHE BOOL "" FORCE)
                SET(WALBERLA_BUILD_WITH_HIP OFF CACHE BOOL "" FORCE)
                set(GPU_AVAILABLE TRUE)
                set(GPU_TYPE "CUDA")
            else()
                message(FATAL_ERROR "GPU build requested but no GPU toolkit found")
            endif()
        else()
            message(FATAL_ERROR "GPU build requested but no GPU compiler found")
        endif()
    endif()
endif()

find_package(MPI REQUIRED)
message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")

add_subdirectory(walberla)
walberla_link_files_to_builddir( "*.prm" )
sweepgen_venv_require(-e ${CMAKE_SOURCE_DIR}/..)
sweepgen_venv_populate()

# Build targets based on configuration
if(BUILD_CPU_ONLY)
    # ============== CPU-ONLY EXECUTABLES ==============
    
    # Heat Equation executables
    add_executable(HeatEquationWithMaterial)
    target_sources(HeatEquationWithMaterial PRIVATE CodegenHeatEquationWithMaterial.cpp)
    target_link_libraries(HeatEquationWithMaterial PRIVATE
            walberla::blockforest walberla::core walberla::field
            walberla::stencil walberla::timeloop walberla::vtk walberla::pde
            MPI::MPI_CXX)
    walberla_generate_sources(HeatEquationWithMaterial
            SCRIPTS HeatEquationKernelWithMaterial.py
            FILE_EXTENSIONS .cpp .hpp)

    add_executable(HeatEquationWithMaterialCPU)
    target_sources(HeatEquationWithMaterialCPU PRIVATE CodegenHeatEquationWithMaterialCPU.cpp)
    target_link_libraries(HeatEquationWithMaterialCPU PRIVATE
            walberla::blockforest walberla::core walberla::field
            walberla::stencil walberla::timeloop walberla::vtk walberla::pde
            MPI::MPI_CXX)
    walberla_generate_sources(HeatEquationWithMaterialCPU
            SCRIPTS HeatEquationKernelWithMaterialCPU.py
            FILE_EXTENSIONS .cpp .hpp)

    add_executable(CodegenHeatEquationCPUScaling)
    target_sources(CodegenHeatEquationCPUScaling PRIVATE CodegenHeatEquationCPUScaling.cpp)
    target_link_libraries(CodegenHeatEquationCPUScaling PRIVATE
            walberla::blockforest walberla::core walberla::field
            walberla::stencil walberla::timeloop walberla::vtk walberla::pde
            MPI::MPI_CXX)
    walberla_generate_sources(CodegenHeatEquationCPUScaling
            SCRIPTS HeatEquationKernelWithMaterialCPU.py
            FILE_EXTENSIONS .cpp .hpp)

    # Couette Flow executables
    add_executable(CouetteFlowCodegen)
    target_sources(CouetteFlowCodegen PRIVATE CouetteFlow.cpp)
    target_link_libraries(CouetteFlowCodegen PRIVATE
            walberla::blockforest walberla::boundary walberla::core 
            walberla::experimental walberla::field walberla::lbm 
            walberla::postprocessing walberla::stencil 
            walberla::timeloop walberla::vtk walberla::sqlite
            MPI::MPI_CXX)
    walberla_generate_sources(CouetteFlowCodegen
            SCRIPTS CouetteFlowSweeps.py
            FILE_EXTENSIONS .cpp .hpp)

    add_executable(CouetteFlowCPUScaling)
    target_sources(CouetteFlowCPUScaling PRIVATE CouetteFlowCPUScaling.cpp)
    target_link_libraries(CouetteFlowCPUScaling PRIVATE
            walberla::blockforest walberla::boundary walberla::core 
            walberla::experimental walberla::field walberla::lbm 
            walberla::stencil walberla::timeloop walberla::vtk
            MPI::MPI_CXX)
    walberla_generate_sources(CouetteFlowCPUScaling
            SCRIPTS CouetteFlowSweeps.py
            FILE_EXTENSIONS .cpp .hpp)

    message(STATUS "CPU executables configured:")
    message(STATUS "  - HeatEquationWithMaterial")
    message(STATUS "  - HeatEquationWithMaterialCPU")
    message(STATUS "  - CodegenHeatEquationCPUScaling")
    message(STATUS "  - CouetteFlowCodegen")
    message(STATUS "  - CouetteFlowCPUScaling")

else()
    # ============== GPU EXECUTABLES ==============
    
    if(GPU_TYPE STREQUAL "HIP")
        # Get GTL library information from environment
        set(PE_MPICH_GTL_DIR_amd_gfx90a "$ENV{PE_MPICH_GTL_DIR_amd_gfx90a}")
        set(PE_MPICH_GTL_LIBS_amd_gfx90a "$ENV{PE_MPICH_GTL_LIBS_amd_gfx90a}")
        
        # Strip the -L and -l prefixes for CMake
        if(PE_MPICH_GTL_DIR_amd_gfx90a)
            string(REPLACE "-L" "" GTL_LIB_DIR "${PE_MPICH_GTL_DIR_amd_gfx90a}")
        endif()
        if(PE_MPICH_GTL_LIBS_amd_gfx90a)
            string(REPLACE "-l" "" GTL_LIB_NAME "${PE_MPICH_GTL_LIBS_amd_gfx90a}")
        endif()
        
        # Debug: Print MPI include directories
        message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")
        message(STATUS "GTL_LIB_DIR: ${GTL_LIB_DIR}")
        message(STATUS "GTL_LIB_NAME: ${GTL_LIB_NAME}")
        
        # ========== Heat Equation GPU Executables ==========
        
        add_executable(HeatEquationWithMaterialGPU)
        target_sources(HeatEquationWithMaterialGPU PRIVATE CodegenHeatEquationWithMaterialGPU.cpp)
        target_link_libraries(HeatEquationWithMaterialGPU PRIVATE
                walberla::blockforest walberla::core walberla::gpu walberla::field
                walberla::stencil walberla::timeloop walberla::vtk walberla::pde
                hip::device
                MPI::MPI_CXX)
        
        # Link GTL library correctly
        if(GTL_LIB_DIR AND GTL_LIB_NAME)
            target_link_directories(HeatEquationWithMaterialGPU PRIVATE ${GTL_LIB_DIR})
            target_link_libraries(HeatEquationWithMaterialGPU PRIVATE ${GTL_LIB_NAME})
        endif()

        # Add compile options for HIP
        target_compile_options(HeatEquationWithMaterialGPU PRIVATE 
            $<$<COMPILE_LANGUAGE:HIP>:--offload-arch=gfx90a>)
        
        set_property(TARGET HeatEquationWithMaterialGPU PROPERTY HIP_ARCHITECTURES "gfx90a")
        
        walberla_generate_sources(HeatEquationWithMaterialGPU
                SCRIPTS HeatEquationKernelWithMaterialGPU.py
                FILE_EXTENSIONS .cpp .hpp)
        
        # Set HIP language for generated source
        set_source_files_properties(
                ${CMAKE_CURRENT_BINARY_DIR}/_gen/HeatEquationWithMaterialGPU/gen/HeatEquationKernelWithMaterialGPU.cpp
                PROPERTIES LANGUAGE HIP)

        # CodegenHeatEquationGPUScaling
        add_executable(CodegenHeatEquationGPUScaling)
        target_sources(CodegenHeatEquationGPUScaling PRIVATE CodegenHeatEquationGPUScaling.cpp)
        target_link_libraries(CodegenHeatEquationGPUScaling PRIVATE
                walberla::blockforest walberla::core walberla::gpu walberla::field
                walberla::stencil walberla::timeloop walberla::vtk walberla::pde
                hip::device
                MPI::MPI_CXX)
        
        # Link GTL library correctly
        if(GTL_LIB_DIR AND GTL_LIB_NAME)
            target_link_directories(CodegenHeatEquationGPUScaling PRIVATE ${GTL_LIB_DIR})
            target_link_libraries(CodegenHeatEquationGPUScaling PRIVATE ${GTL_LIB_NAME})
        endif()

        target_compile_options(CodegenHeatEquationGPUScaling PRIVATE 
            $<$<COMPILE_LANGUAGE:HIP>:--offload-arch=gfx90a>)
        
        set_property(TARGET CodegenHeatEquationGPUScaling PROPERTY HIP_ARCHITECTURES "gfx90a")
        
        walberla_generate_sources(CodegenHeatEquationGPUScaling
                SCRIPTS HeatEquationKernelWithMaterialGPU.py
                FILE_EXTENSIONS .cpp .hpp)
                
        set_source_files_properties(
                ${CMAKE_CURRENT_BINARY_DIR}/_gen/CodegenHeatEquationGPUScaling/gen/HeatEquationKernelWithMaterialGPU.cpp
                PROPERTIES LANGUAGE HIP)

        # ========== Couette Flow GPU Executables ==========
        
        add_executable(CouetteFlowGPUScaling)
        target_sources(CouetteFlowGPUScaling PRIVATE CouetteFlowGPUScaling.cpp)
        target_link_libraries(CouetteFlowGPUScaling PRIVATE
                walberla::blockforest walberla::boundary walberla::core 
                walberla::experimental walberla::gpu walberla::field 
                walberla::lbm walberla::stencil walberla::timeloop walberla::vtk
                hip::device
                MPI::MPI_CXX)
        
        if(GTL_LIB_DIR AND GTL_LIB_NAME)
            target_link_directories(CouetteFlowGPUScaling PRIVATE ${GTL_LIB_DIR})
            target_link_libraries(CouetteFlowGPUScaling PRIVATE ${GTL_LIB_NAME})
        endif()

        target_compile_options(CouetteFlowGPUScaling PRIVATE 
            $<$<COMPILE_LANGUAGE:HIP>:--offload-arch=gfx90a>)
        
        set_property(TARGET CouetteFlowGPUScaling PROPERTY HIP_ARCHITECTURES "gfx90a")
        
        walberla_generate_sources(CouetteFlowGPUScaling
                SCRIPTS CouetteFlowSweeps.py
                SCRIPT_ARGS --target=gpu
                FILE_EXTENSIONS .cpp .hpp)
                
        set_source_files_properties(
                ${CMAKE_CURRENT_BINARY_DIR}/_gen/CouetteFlowGPUScaling/gen/CouetteFlowSweeps.cpp
                PROPERTIES LANGUAGE HIP)

        message(STATUS "HIP GPU executables configured:")
        message(STATUS "  - HeatEquationWithMaterialGPU")
        message(STATUS "  - CodegenHeatEquationGPUScaling")
        message(STATUS "  - CouetteFlowGPUScaling")
    
    elseif(GPU_TYPE STREQUAL "CUDA")
        # ========== Heat Equation GPU Executables ==========
        
        add_executable(HeatEquationWithMaterialGPU)
        target_sources(HeatEquationWithMaterialGPU PRIVATE CodegenHeatEquationWithMaterialGPU.cpp)
        target_link_libraries(HeatEquationWithMaterialGPU PRIVATE
                walberla::blockforest walberla::core walberla::gpu walberla::field
                walberla::stencil walberla::timeloop walberla::vtk walberla::pde
                MPI::MPI_CXX)
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
                walberla::stencil walberla::timeloop walberla::vtk walberla::pde
                MPI::MPI_CXX)
        # Set CUDA properties for scaling executable
        set_property(TARGET CodegenHeatEquationGPUScaling PROPERTY CUDA_SEPARABLE_COMPILATION ON)
        walberla_generate_sources(CodegenHeatEquationGPUScaling
                SCRIPTS HeatEquationKernelWithMaterialGPU.py
                FILE_EXTENSIONS .cu .hpp)
        set_source_files_properties(
                ${CMAKE_CURRENT_BINARY_DIR}/_gen/CodegenHeatEquationGPUScaling/gen/HeatEquationKernelWithMaterialGPU.cu
                PROPERTIES LANGUAGE CUDA)

        # ========== Couette Flow GPU Executables ==========
        
        add_executable(CouetteFlowGPUScaling)
        target_sources(CouetteFlowGPUScaling PRIVATE CouetteFlowGPUScaling.cpp)
        target_link_libraries(CouetteFlowGPUScaling PRIVATE
                walberla::blockforest walberla::boundary walberla::core 
                walberla::experimental walberla::gpu walberla::field 
                walberla::lbm walberla::stencil walberla::timeloop walberla::vtk
                MPI::MPI_CXX)
        set_property(TARGET CouetteFlowGPUScaling PROPERTY CUDA_SEPARABLE_COMPILATION ON)
        walberla_generate_sources(CouetteFlowGPUScaling
                SCRIPTS CouetteFlowSweeps.py 
                SCRIPT_ARGS --target=gpu
                FILE_EXTENSIONS .cu .hpp)
        set_source_files_properties(
                ${CMAKE_CURRENT_BINARY_DIR}/_gen/CouetteFlowGPUScaling/gen/CouetteFlowSweeps.cu
                PROPERTIES LANGUAGE CUDA)

        message(STATUS "CUDA GPU executables configured:")
        message(STATUS "  - HeatEquationWithMaterialGPU")
        message(STATUS "  - CodegenHeatEquationGPUScaling")
        message(STATUS "  - CouetteFlowGPUScaling")
    endif()
endif()
