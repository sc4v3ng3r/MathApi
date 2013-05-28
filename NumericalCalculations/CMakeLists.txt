project (NumericalCalculations)

set(SOURCES_C PolynomialResults.c PolynomialRoot.c)

add_library(libnc ${SOURCES_C})
add_dependencies(libnc base)