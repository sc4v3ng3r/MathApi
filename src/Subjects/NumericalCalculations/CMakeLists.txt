project (NumericalCalculations)

set(NUMERICAL_SOURCES_C PolynomialResults.c PolynomialRoot.c)

add_library(nc SHARED ${NUMERICAL_SOURCES_C})
add_dependencies(nc base tables)