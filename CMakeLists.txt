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

add_library(magnetosun src/cpp/synoptic_map.cpp)
target_include_directories(magnetosun PUBLIC src/cpp)
target_link_libraries(magnetosun PUBLIC PkgConfig::CFITSIO)

add_executable(gong_export src/cpp/gong_export.cpp)
target_link_libraries(gong_export PRIVATE magnetosun)

enable_testing()

add_executable(test_interop tests/interop/test_interop.cpp tests/interop/fill_matrix.f90)
add_test(NAME interop COMMAND test_interop)

add_executable(test_synoptic_map tests/test_synoptic_map.cpp)
target_link_libraries(test_synoptic_map PRIVATE magnetosun)
add_test(NAME synoptic_map COMMAND test_synoptic_map)
