project(ResultsTable)

set(RESULTS_TABLES_SOURCES_C ResultsTable.c InteractiveResultsTable.c)

add_library(tables SHARED ${RESULTS_TABLES_SOURCES_C})
add_dependencies(tables base)