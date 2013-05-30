project (MathBase)

set(SOURCES_C Matrix.c OrderedPair.c Polynomial.c ResultsTable.c)
set(MATH -lm)

add_library(base SHARED ${SOURCES_C})
target_link_libraries(base ${MATH})