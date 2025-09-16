cmake_minimum_required(VERSION 3.18.0)
include(ExternalProject)
project(QMCKLBENCH)


#Clone and build trexio, default compilation
ExternalProject_Add(libtrexio
	GIT_REPOSITORY https://github.com/TREX-CoE/trexio.git
	BUILD_IN_SOURCE ON
	PREFIX ${CMAKE_CURRENT_BINARY_DIR}/libtrexio
	CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/libtrexio/src/libtrexio/autogen.sh
	COMMAND ${CMAKE_CURRENT_BINARY_DIR}/libtrexio/src/libtrexio/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/libtrexio/__install
	BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/libtrexio/__install/lib/libtrexio.so)

SET(TREXIO_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/libtrexio/__install/include)
SET(TREXIO_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/libtrexio/__install/lib)

add_library(TREXIO SHARED IMPORTED)
set_target_properties(TREXIO PROPERTIES IMPORTED_LOCATION ${TREXIO_LIB_DIR}/libtrexio.so)

include_directories(${TREXIO_INCLUDE_DIR})


#Clone and build qmckl, forward CMAKE flags and compiler to control its compilation
ExternalProject_Add(libqmckl
	GIT_REPOSITORY https://github.com/TREX-CoE/qmckl.git
	BUILD_IN_SOURCE ON
	PREFIX ${CMAKE_CURRENT_BINARY_DIR}/libqmckl
	CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/libqmckl/src/libqmckl/autogen.sh 
	COMMAND ${CMAKE_CURRENT_BINARY_DIR}/libqmckl/src/libqmckl/configure CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER} FCFLAGS=${CMAKE_Fortran_FLAGS} CFLAGS=${CMAKE_C_FLAGS} --prefix=${CMAKE_CURRENT_BINARY_DIR}/libqmckl/__install --enable-static --enable-hpc --with-openmp --with-trexio=${CMAKE_CURRENT_BINARY_DIR}/libtrexio/__install
	BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/libqmckl/__install/lib/libqmckl.so)

SET(QMCKL_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/libqmckl/__install/include)
SET(QMCKL_LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/libqmckl/__install/lib)

add_dependencies(libqmckl libtrexio)

add_library(QMCKL SHARED IMPORTED)
set_target_properties(QMCKL PROPERTIES IMPORTED_LOCATION ${QMCKL_LIB_DIR}/libqmckl.so)

include_directories(${QMCKL_INCLUDE_DIR})


add_executable(bench_aos src/bench_aos.c)
target_link_libraries(bench_aos QMCKL TREXIO)
add_dependencies(bench_aos libqmckl)

add_executable(bench_mos src/bench_mos.c)
target_link_libraries(bench_mos QMCKL TREXIO)
add_dependencies(bench_mos libqmckl)

add_executable(bench_jastrow src/bench_jastrow.c)
target_link_libraries(bench_jastrow QMCKL TREXIO)
add_dependencies(bench_jastrow libqmckl)

add_executable(bench_prec src/bench_prec.c)
target_link_libraries(bench_prec QMCKL TREXIO)
add_dependencies(bench_prec libqmckl)
