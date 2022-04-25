
################################################################################
# Molecular dynamics
################################################################################

ADD_EXECUTABLE(md main.cpp)

FIND_PACKAGE(Threads)
TARGET_LINK_LIBRARIES(md PRIVATE ${CMAKE_THREAD_LIBS_INIT})

TARGET_LINK_LIBRARIES(md PRIVATE ode)
