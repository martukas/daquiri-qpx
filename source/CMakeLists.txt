include(InstallConfig)
include(BuildInfo)

find_package(date REQUIRED)
find_package(fmt REQUIRED)
find_package(h5cpp 0.1.0 REQUIRED)
find_package(Threads REQUIRED)
find_package(GraylogLogger REQUIRED)
add_definitions("-DGRAYLOG")

add_subdirectory(core)
add_subdirectory(consumers)
add_subdirectory(producers)
add_subdirectory(importers)

#ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT COMPONENTS Hist)
if(ROOT_FOUND)
  add_subdirectory(optimizerROOT)
  set(optimizer_LIBRARIES daquiri_optimizerROOT)
endif()

add_subdirectory(optimizerBFGS)
#set(optimizer_LIBRARIES daquiri_optimizerROOT)

add_subdirectory(cmd)
add_subdirectory(gui)
