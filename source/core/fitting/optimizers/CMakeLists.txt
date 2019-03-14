set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/fit_result.cpp
  ${dir}/BFGS.cpp
#  ${dir}/dlib_adapter.cpp
#  ${dir}/gsl_adapter.cpp
  ${dir}/optlib_adapter.cpp
  )

set(HEADERS
  ${dir}/fittable_function.h
  ${dir}/fit_result.h
  ${dir}/abstract_optimizer.h
  ${dir}/BFGS.h
#  ${dir}/dlib_adapter.h
#  ${dir}/gsl_adapter.h
  ${dir}/optlib_adapter.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

