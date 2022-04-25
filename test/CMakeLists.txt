
################################################################################
# Test
################################################################################

ADD_EXECUTABLE(runtest main.cpp)

FIND_PACKAGE(Threads)
TARGET_LINK_LIBRARIES(runtest PRIVATE ${CMAKE_THREAD_LIBS_INIT})

TARGET_LINK_LIBRARIES(runtest PRIVATE OdeSolver)

ENABLE_TESTING()

ADD_TEST(NAME runtest COMMAND runtest)
