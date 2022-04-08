
################################################################################
# Molecular dynamics
################################################################################

ADD_EXECUTABLE(moleculardynamics main.cpp)

FIND_PACKAGE(Threads)
TARGET_LINK_LIBRARIES(moleculardynamics PRIVATE ${CMAKE_THREAD_LIBS_INIT})

TARGET_LINK_LIBRARIES(moleculardynamics PRIVATE OdeSolver)
