include(InstallConfig)
include(BuildInfo)

find_package(date REQUIRED)
find_package(fmt REQUIRED)
find_package(h5cpp 0.1.0 REQUIRED)
find_package(Threads REQUIRED)
find_package(GraylogLogger REQUIRED)
add_definitions("-DGRAYLOG")

find_package(GSL REQUIRED)    # See below (2)
find_package(dlib REQUIRED)
set(optimizer_LIBRARIES dlib::dlib GSL::gsl GSL::gslcblas)

add_subdirectory(core)
add_subdirectory(consumers)
add_subdirectory(producers)
add_subdirectory(importers)

add_subdirectory(cmd)
add_subdirectory(gui)
