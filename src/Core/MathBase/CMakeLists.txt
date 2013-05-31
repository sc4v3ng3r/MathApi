project (MathBase)

set(MATH_BASE_SOURCES_C Matrix.c OrderedPair.c Polynomial.c)
set(MATH -lm)

add_library(base SHARED ${MATH_BASE_SOURCES_C})
target_link_libraries(base ${MATH})