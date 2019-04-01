set(dir ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(optimizers)
add_subdirectory(data_model)
add_subdirectory(parameter)

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
