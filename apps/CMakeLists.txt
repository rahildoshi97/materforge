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

# Build mode detection
option(BUILD_CPU_ONLY "Build only CPU targets (ignore GPU)" OFF)
option(ENABLE_GPU_SUPPORT "Enable GPU support" ON)

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

add_subdirectory(walberla)
add_subdirectory(sfg-walberla)
walberla_codegen_venv_require(-e ${CMAKE_SOURCE_DIR}/..)
walberla_codegen_venv_populate()

# Build targets based on configuration
if(BUILD_CPU_ONLY)
    # CPU-only executables
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

    message(STATUS "CPU executables will be built: HeatEquationWithMaterial, HeatEquationWithMaterialCPU, CodegenHeatEquationCPUScaling")

else()
    if(GPU_TYPE STREQUAL "HIP")
        find_package(MPI REQUIRED)
        message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")
        
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
        
        # Create executables with full MPI::MPI_CXX target
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

        message(STATUS "HIP GPU executables will be built: HeatEquationWithMaterialGPU, CodegenHeatEquationGPUScaling")
    
    elseif(GPU_TYPE STREQUAL "CUDA")
        # Find MPI for CUDA builds too
        find_package(MPI REQUIRED)

        # Debug: Print MPI include directories
        message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")

        # Create HeatEquationWithMaterialGPU executable
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

        message(STATUS "CUDA GPU executables will be built: HeatEquationWithMaterialGPU, CodegenHeatEquationGPUScaling")
    endif()
endif()
