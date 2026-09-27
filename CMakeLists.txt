cmake_minimum_required(VERSION 3.20)
project(magnetosun LANGUAGES Fortran CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# En Debug avec gfortran : contrôle des bornes de tableaux à l'exécution
# (les erreurs d'indice (l, m) sont le bug le plus probable du projet).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -std=f2018 -Wall -Wextra -fcheck=all")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Wall -Wextra -Wpedantic")
endif()

# cfitsio : Homebrew sur macOS ; sur ROMEO, module ou installation dans $HOME
# (ajouter alors son dossier lib/pkgconfig à PKG_CONFIG_PATH).
find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED IMPORTED_TARGET cfitsio)

# LAPACK : Accelerate sur macOS ; OpenBLAS, MKL ou équivalent sur ROMEO.
find_package(LAPACK REQUIRED)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Noyaux numériques Fortran (harmoniques sphériques, moindres carrés, PFSS).
add_library(magnetosun_kernels
  src/fortran/legendre.f90
  src/fortran/sh_fit.f90
  src/fortran/pfss.f90
  src/fortran/capi.f90)
target_include_directories(magnetosun_kernels PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(magnetosun_kernels PUBLIC LAPACK::LAPACK)

add_library(magnetosun src/cpp/synoptic_map.cpp src/cpp/sh_model.cpp src/cpp/insitu.cpp)
target_include_directories(magnetosun PUBLIC src/cpp)
target_link_libraries(magnetosun PUBLIC PkgConfig::CFITSIO magnetosun_kernels)

foreach(program gong_export fit_map pfss_map polarity)
  add_executable(${program} src/cpp/${program}.cpp)
  target_link_libraries(${program} PRIVATE magnetosun)
endforeach()

enable_testing()

add_executable(test_interop tests/interop/test_interop.cpp tests/interop/fill_matrix.f90)
add_test(NAME interop COMMAND test_interop)

foreach(name synoptic_map sh_model insitu)
  add_executable(test_${name} tests/test_${name}.cpp)
  target_link_libraries(test_${name} PRIVATE magnetosun)
  add_test(NAME ${name} COMMAND test_${name})
endforeach()

foreach(name legendre sh_fit pfss)
  add_executable(test_${name} tests/test_${name}.f90)
  target_link_libraries(test_${name} PRIVATE magnetosun_kernels)
  add_test(NAME ${name} COMMAND test_${name})
endforeach()
