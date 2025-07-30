cmake_minimum_required( VERSION 3.24 )

project ( test_materforge )

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
SET ( WALBERLA_BUILD_WITH_CODEGEN OFF CACHE BOOL "" )
SET ( WALBERLA_BUILD_WITH_OPENMP ON CACHE BOOL "" )
SET ( WALBERLA_BUILD_WITH_MPI OFF CACHE BOOL "" )


add_subdirectory(walberla)
add_subdirectory(sfg-walberla)

walberla_codegen_venv_require(
        -e ${CMAKE_SOURCE_DIR}/..
)
walberla_codegen_venv_populate()

add_executable ( HeatEquationWithMaterial )
target_sources( HeatEquationWithMaterial PRIVATE CodegenHeatEquationWithMaterial.cpp )
target_link_libraries( HeatEquationWithMaterial PRIVATE walberla::blockforest walberla::core walberla::field walberla::stencil walberla::timeloop walberla::vtk walberla::pde )

walberla_generate_sources( HeatEquationWithMaterial
        SCRIPTS HeatEquationKernelWithMaterial.py
        FILE_EXTENSIONS .cpp .hpp
)
