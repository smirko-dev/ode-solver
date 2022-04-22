
################################################################################
# Test
################################################################################

ADD_EXECUTABLE(test main.cpp)

FIND_PACKAGE(Threads)
TARGET_LINK_LIBRARIES(test PRIVATE ${CMAKE_THREAD_LIBS_INIT})

TARGET_LINK_LIBRARIES(test PRIVATE OdeSolver)
