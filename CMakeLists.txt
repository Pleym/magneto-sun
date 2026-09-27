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

enable_testing()

add_executable(test_interop tests/interop/test_interop.cpp tests/interop/fill_matrix.f90)
add_test(NAME interop COMMAND test_interop)
