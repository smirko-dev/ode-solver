
################################################################################
# Lorenz attractor
################################################################################

ADD_EXECUTABLE(la main.cpp)

FIND_PACKAGE(Threads)
TARGET_LINK_LIBRARIES(la PRIVATE ${CMAKE_THREAD_LIBS_INIT})

TARGET_LINK_LIBRARIES(la PRIVATE ode)
