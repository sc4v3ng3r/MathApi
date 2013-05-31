project(examples)

set(SOURCE_C main.c)

add_executable(tst ${SOURCE_C})
target_link_libraries(tst base tables nc la)