cmake_minimum_required( VERSION 3.24 )

project ( test_materforge LANGUAGES CXX)

SET( CMAKE_CXX_STANDARD 20 )
SET( CMAKE_CXX_STANDARD_REQUIRED ON )

#   waLBerla Options
SET ( WALBERLA_BUILD_TESTS OFF CACHE BOOL "" )
SET ( WALBERLA_BUILD_BENCHMARKS OFF CACHE BOOL "" )
SET ( WALBERLA_BUILD_TOOLS OFF CACHE BOOL "" )
SET ( WALBERLA_BUILD_TUTORIALS OFF CACHE BOOL "" )
SET ( WALBERLA_BUILD_SHOWCASES OFF CACHE BOOL "" )
SET ( WALBERLA_BUILD_DOC OFF CACHE BOOL "" )
SET ( WALBERLA_BUILD_WITH_PYTHON OFF CACHE BOOL "" )
# SET(WALBERLA_BUILD_WITH_CODEGEN OFF CACHE BOOL "")
# SET(WALBERLA_BUILD_WITH_OPENMP ON CACHE BOOL "")
# SET(WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "")
SET(WALBERLA_ENABLE_SWEEPGEN ON CACHE BOOL "")
SET(WALBERLA_SWEEPGEN_MANAGED_VENV ON CACHE BOOL "")

find_package(MPI REQUIRED)
message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
message(STATUS "MPI_CXX_COMPILE_FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
message(STATUS "MPI_CXX_LINK_FLAGS: ${MPI_CXX_LINK_FLAGS}")

add_subdirectory(walberla)
sweepgen_venv_require(-e ${CMAKE_SOURCE_DIR}/..)
sweepgen_venv_populate()

add_executable ( HeatEquationWithMaterial )
target_sources( HeatEquationWithMaterial PRIVATE CodegenHeatEquationWithMaterial.cpp )
target_link_libraries( HeatEquationWithMaterial PRIVATE walberla::blockforest walberla::core walberla::field walberla::stencil walberla::timeloop walberla::vtk walberla::pde MPI::MPI_CXX )

walberla_generate_sources( HeatEquationWithMaterial
        SCRIPTS HeatEquationKernelWithMaterial.py
        FILE_EXTENSIONS .cpp .hpp
)
