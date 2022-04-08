
################################################################################
# Planet dynamics
################################################################################

ADD_EXECUTABLE(pd main.cpp)

FIND_PACKAGE(Threads)
TARGET_LINK_LIBRARIES(pd PRIVATE ${CMAKE_THREAD_LIBS_INIT})

TARGET_LINK_LIBRARIES(pd PRIVATE OdeSolver)
