
CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

INCLUDE(InstallRequiredSystemLibraries)

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
#SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

################################################################################
# ODE Solver
################################################################################

ADD_LIBRARY(OdeSolver INTERFACE)

TARGET_SOURCES(OdeSolver INTERFACE
    Vector.h
    OdeFunction.h
    OdeSolver.h
    Euler.h
    MidPoint.h
    RungeKutta.h
    VelocityVerlet.h
)

TARGET_INCLUDE_DIRECTORIES(OdeSolver INTERFACE ${CMAKE_CURRENT_LIST_DIR})

################################################################################
# Subdirectories
################################################################################

ADD_SUBDIRECTORY(moleculardynamics)
ADD_SUBDIRECTORY(planetdynamics)
