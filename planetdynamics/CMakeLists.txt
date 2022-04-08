
################################################################################
# Planet dynamics
################################################################################

ADD_EXECUTABLE(planetdynamics main.cpp)

FIND_PACKAGE(Threads)
TARGET_LINK_LIBRARIES(planetdynamics PRIVATE ${CMAKE_THREAD_LIBS_INIT})

TARGET_LINK_LIBRARIES(planetdynamics PRIVATE OdeSolver)
