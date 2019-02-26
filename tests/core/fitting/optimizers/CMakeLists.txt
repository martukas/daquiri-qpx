set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/rosenbrock.h
#  ${dir}/dlib_adapter.cpp
#  ${dir}/gsl_adapter.cpp
#  ${dir}/BFGS.cpp
  ${dir}/optlib_adapter.cpp
  )

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
