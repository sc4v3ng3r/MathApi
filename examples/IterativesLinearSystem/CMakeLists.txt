project(interactives)

set(SOURCE_C main.c)

add_executable(interactives ${SOURCE_C})
target_link_libraries(interactives la tables base)