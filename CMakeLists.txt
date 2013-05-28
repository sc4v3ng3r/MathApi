project (testeAPI)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#add_subdirectory(Core)
add_subdirectory(MathBase)
add_subdirectory(LinearAlgebra)
add_subdirectory(NumericalCalculations)

set(SOURCES_C main.c)

add_executable(tst ${SOURCES_C})
target_link_libraries(tst lalib base libnc)