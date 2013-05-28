project(LinearAlgebra)

set(SOURCES_C LinearSystem.c)

add_library(lalib SHARED ${SOURCES_C})
add_dependencies(lalib base)