project(LinearAlgebra)

set(LINEAR_ALGEBRA_SOURCES_C LinearSystem.c)

add_library(la SHARED ${LINEAR_ALGEBRA_SOURCES_C})
add_dependencies(la base tables)