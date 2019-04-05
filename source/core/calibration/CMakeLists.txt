set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/function_factory.cpp
  ${dir}/calibration.cpp

  ${dir}/abstract_poly.cpp
  ${dir}/polynomial.cpp
  ${dir}/sqrt_poly.cpp
  )

set(HEADERS
  ${dir}/function.h
  ${dir}/function_factory.h
  ${dir}/calibration.h

  ${dir}/abstract_poly.h
  ${dir}/polynomial.h
  ${dir}/sqrt_poly.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
